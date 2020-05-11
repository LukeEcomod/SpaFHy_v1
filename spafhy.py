# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 16:18:57 2016

@author: slauniai


"""
# import timeit
import os
import numpy as np
from canopygrid import CanopyGrid
from bucketgrid import BucketGrid
from topmodel import Topmodel_Homogenous as Topmodel
from spafhy_io import preprocess_soildata

eps = np.finfo(float).eps  # machine epsilon

""" ************** SpaFHy v1.0 ************************************

Simple spatial hydrology and catchment water balance model.

CONSISTS OF THREE CLASSES, defined in separate modules: 
    CanopyGrid - vegetation and snowpack water storages and flows
    BucketGrid - topsoil bucket model (root zone / topsoil water storage)
    Topmodel - integration to catchment scale using Topmodel -concept
HELPER FUNCTIONS:
    spafhy_parameters - parameter definition file
    spafhy_io - utility functions for data input & output
 
MAIN PROGRAM:   
    spafhy_driver is main program, call it as
    outargs = spathy_driver(spathyparamfile, args)
    
    spathyparamfile - path to parameter file, default is 'spathy_default.ini'
    soil type dependent parameters are in 'soilparam.ini'

NEEDS 2D gis rasters in ascii-grid format

CanopyGrid & BucketGrid can be initialized from gis-data or set to be spatially constant

ToDo:
    CanopyGrid:
        -include topographic shading to radiation received at canopy top
        -radiation-based snowmelt coefficient
        -add simple GPP-model; 2-step Farquhar or LUE-based approach
    BucketGrid:
        -make soil hydrologic properties more realistic e.g. using pedotransfer functions
        -kasvupaikkatyyppi (multi-NFI) --> soil properties
        -add soil frost model, simplest would be Stefan equation with coefficients
         modified based on snow insulation  --> we need snow density algorithm: SWE <-----> depth
    Topmodel:
        -think of definging 'relative m & to grids' (soil-type & elevation-dependent?)
         and calibrate 'catchment averages'
        -topmodel gives 'saturated zone storage deficit in [m]'. This can be
         converted to gwl proxy (?) if: local water retention characteristics are known 
         & hydrostatic equilibrium assumes. Look which water retention models are 
         analytically integrable (Campbell, brooks-corey?)
    
    Graphics and analysis of results:
        -make ready functions

(C) Samuli Launiainen 10/2016-->    

VERSION 07.05.2020 / equations correspond to HESS paper. Some refinements in canopygrid
parameterization (multiple vegetation types) and fix in bucketgrid mass-balance calculations.

"""

def initialize(pgen, pcpy, pbu, ptop, psoil, gisdata, ncf=True,
               cpy_outputs=False, bu_outputs=False, top_outputs=False):
    """ 
    ******************** sets up SpaFHy  **********************
    
    Normal SpaFHy run without parameter optimization
    1) gets parameters as input arguments
    2) reads GIS-data, here predefined format for Seurantaverkko -cathcments
    
    3) creates CanopyGrid (cpy), BucketGrid (bu) and Topmodel (top) -objects 
        within Spathy-object (spa), and temporary outputs
    4) creates netCDF -file for outputs if 'ncf' = True
    
    5) returns following outputs:
        spa - spathy object
        ncf - netcdf file handle (if ncf=True)
        ncf_file - ncf file name (if ncf=True)
        
    IN:
        pgen, pcpy, pbu, psoil - parameter dictionaries
        gisdata - dict of 2d np.arrays containing gis-data with keys:
            cmask - catchment mask; integers within np.Nan outside
            LAI_conif [m2m-2]
            LAI_decid [m2m-2]
            hc, canopy closure [m]
            fc, canopy closure fraction [-]
            soil, soil type integer code 1-5
            flowacc - flow accumulation [units]
            slope - local surface slope [units]
            
            cellsize - gridcell size
            lon0 - x-grid
            lat0 - y-grid
        cpy_outputs, bu_, top_ - True saves cpy, bu and top outputs to lists within each object. 
            Use only for testing, memory issue!
    OUT:
        spa - spathy object
    """

    # start_time = timeit.default_timer()

    # preprocess soildata --> dict used in BucketModel initialization    
    soildata = preprocess_soildata(pbu, psoil, gisdata['soilclass'],
                                   gisdata['cmask'], pgen['spatial_soil'])

    # inputs for CanopyGrid initialization: make sure vegetation data is masked by domain
    cstate = pcpy['state']
        
    for key in cstate.keys():
        if 'LAI_' in key: # reads gisdata LAI-layers to cstate
            cstate[key] = gisdata[key]
            
        cstate[key] *= gisdata['cmask']

    pcpy['state'] = cstate
    del cstate
    
    if ncf:
        flatten = True

    """ greate SpatHy object """
    spa = SpaFHy(pgen, pcpy, ptop, soildata, gisdata, cpy_outputs=cpy_outputs,
                 bu_outputs=bu_outputs, top_outputs=top_outputs, flatten=flatten)

    # create netCDF output file
    dlat, dlon = np.shape(spa.GisData['cmask'])

    resultsfile = os.path.join(spa.pgen['results_folder'], spa.pgen['ncf_file'])

    ncf, ncf_file = initialize_netCDF(ID=spa.id, fname=resultsfile,
                                      lat0=spa.GisData['lat0'],
                                      lon0=spa.GisData['lon0'],
                                      dlat=dlat, dlon=dlon, dtime=None)
            
    print('********* created SpaFHy instance *********')
    #print('Loops total [s]: ', timeit.default_timer() - start_time)
    
    if ncf:
        return spa, ncf, ncf_file
    else:
        return spa

"""
******************************************************************************
            ----------- SpaFHy model class --------
******************************************************************************
"""

class SpaFHy():
    """
    SpaFHy model class
    """
    def __init__(self, pgen, pcpy, ptop, soildata, gisdata,
                 cpy_outputs=False, bu_outputs=False, top_outputs=False,
                 flatten=True):

        self.dt = pgen['dt']  # s
        self.id = pgen['catchment_id']
        self.spinup_end = pgen['spinup_end']
        self.pgen = pgen
        self.step_nr = 0
        self.ncf_file = pgen['ncf_file']

        self.GisData = gisdata
        self.cmask = self.GisData['cmask']
        self.gridshape = np.shape(self.cmask)
        cmask= self.cmask.copy()

        flowacc = gisdata['flowacc'].copy()
        slope = gisdata['slope'].copy()        
        
        """
        flatten=True omits cells outside catchment
        """
        if flatten:
            ix = np.where(np.isfinite(cmask))
            cmask = cmask[ix].copy()
            flowacc = flowacc[ix].copy()
            slope = slope[ix].copy()
            
            for key in pcpy['state']:
                pcpy['state'][key] = pcpy['state'][key][ix].copy()
                        
            for key in soildata:
                soildata[key] = soildata[key][ix].copy()
                
            self.ix = ix  # indices to locate back to 2d grid       
   
        """--- initialize CanopyGrid ---"""
        self.cpy = CanopyGrid(pcpy, pcpy['state'], outputs=cpy_outputs)

        """--- initialize BucketGrid ---"""
        self.bu = BucketGrid(spara=soildata, outputs=bu_outputs)

        """ --- initialize Topmodel --- """
        self.top=Topmodel(ptop, self.GisData['cellsize']**2, cmask,
                          flowacc, slope, outputs=top_outputs)

    def run_timestep(self, forc, ncf=False, flx=False, ave_flx=False):
        """ 
        Runs SpaFHy for one timestep starting from current state
        Args:
            forc - dictionary or pd.DataFrame containing forcing values for the timestep
            ncf - netCDF -file handle, for outputs
            flx - returns flux and state grids to caller as dict
            ave_flx - returns averaged fluxes and states to caller as dict
        Returns:
            optional
        """
        doy = forc['doy']
        ta = forc['T']
        vpd = forc['VPD'] + eps
        rg = forc['Rg']
        par = forc['Par'] + eps
        prec = forc['Prec']
        co2 = forc['CO2']
        u = forc['U'] + eps

        # run Topmodel
        # catchment average ground water recharge [m per unit area]
        RR = self.bu._drainage_to_gw * self.top.CellArea / self.top.CatchmentArea
        qb, _, qr, fsat = self.top.run_timestep(RR)

        # run CanopyGrid
        potinf, trfall, interc, evap, et, transpi, efloor, mbe = \
            self.cpy.run_timestep(doy, self.dt, ta, prec, rg, par, vpd, U=u, CO2=co2,
                                  beta=self.bu.Ree, Rew=self.bu.Rew, P=101300.0)

        # run BucketGrid water balance
        infi, roff, drain, tr, eva, mbes = self.bu.watbal(dt=self.dt, rr=1e-3*potinf, tr=1e-3*transpi,
                                                       evap=1e-3*efloor, retflow=qr)

        # catchment average [m per unit area] saturation excess --> goes to stream
        # as surface runoff
        qs = np.nansum(roff)*self.top.CellArea / self.top.CatchmentArea


        """ outputs """
        # updates state and returns results as dict
        if hasattr(self.top, 'results'):
            self.top.results['Qt'].append(qb + qs + eps)  # total runoff

        if ncf:
            # writes to netCDF -file at every timestep; bit slow - should
            # accumulate into temporary variables and save every 10 days? 
            # for netCDF output, must run in Flatten=True              
            k = self.step_nr
            
            # canopygrid
            ncf['cpy']['LAI'][k,:,:] = self._to_grid(self.cpy.LAI)
            ncf['cpy']['W'][k,:,:] = self._to_grid(self.cpy.W)
            ncf['cpy']['SWE'][k,:,:] = self._to_grid(self.cpy.SWE)
            ncf['cpy']['Trfall'][k,:,:] = self._to_grid(trfall) 
            ncf['cpy']['Potinf'][k,:,:] = self._to_grid(potinf)
            ncf['cpy']['ET'][k,:,:] = self._to_grid(et)
            ncf['cpy']['Transpi'][k,:,:] = self._to_grid(transpi)
            ncf['cpy']['Efloor'][k,:,:] = self._to_grid(efloor)            
            ncf['cpy']['Evap'][k,:,:] = self._to_grid(evap)
            ncf['cpy']['Inter'][k,:,:] = self._to_grid(interc)
            ncf['cpy']['Mbe'][k,:,:] = self._to_grid(mbe)              

            # bucketgrid
            ncf['bu']['Drain'][k,:,:] = self._to_grid(drain)
            ncf['bu']['Infil'][k,:,:] = self._to_grid(infi)
            ncf['bu']['Wliq'][k,:,:] = self._to_grid(self.bu.Wliq)
            ncf['bu']['Wliq_top'][k,:,:] = self._to_grid(self.bu.Wliq_top)            
            ncf['bu']['PondSto'][k,:,:] = self._to_grid(self.bu.PondSto)
            ncf['bu']['Roff'][k,:,:] = self._to_grid(roff)
            ncf['bu']['Retflow'][k,:,:] = self._to_grid(qr)
            ncf['bu']['Mbe'][k,:,:] = self._to_grid(mbes)              

            # topmodel
            ncf['top']['Qb'][k] = qb
            ncf['top']['Qs'][k] = qs
            ncf['top']['Qt'][k] = qb + qs  # total runoff
            ncf['top']['R'][k] = RR
            ncf['top']['fsat'][k] = fsat
            ncf['top']['S'][k] = self.top.S
            ss = self.top.local_s(self.top.S)
            ss[ss < 0] = 0.0
            ncf['top']['Sloc'][k,:,:] = self._to_grid(ss)
            del ss

        # update step number
        self.step_nr += 1
        
        if flx:  # returns flux and state variables as grids
            flx = {'cpy': {'ET': et, 'Transpi': transpi, 'Evap': evap,
                           'Efloor': efloor, 'Inter': interc, 'Trfall': trfall,
                           'Potinf': potinf},
                   'bu': {'Infil': infi, 'Drain': drain},
                   'top': {'Qt': qb + qs, 'Qb': qb, 'Qs': qs, 'fsat': fsat}
                  }
            return flx        

        if ave_flx: # returns average fluxes and state variables
            flx = {
                    'ET': np.nanmean(et),
                    'E': np.nanmean(evap),
                    'Ef': np.nanmean(efloor),
                    'Tr': np.nanmean(transpi),
                    'SWE': np.nanmean(self.cpy.SWE),
                    'Drain': 1e3*RR,
                    'Qt': 1e3 * (qb + qs),
                    'S': self.top.S,
                    'fsat': fsat,
                    'Prec': prec * self.dt,
                    'Rg': rg,
                    'Ta': ta,
                    'VPD': vpd
                   }
            return flx
                         
    def _to_grid(self, x):
        """
        converts variable x back to original grid for NetCDF outputs
        """
        if self.ix:
            a = np.full(self.gridshape, np.NaN)
            a[self.ix] = x
        else: # for non-flattened, return
            a = x
        return a

""" ******* netcdf output file ****** """

def initialize_netCDF(ID, fname, lat0, lon0, dlat, dlon, dtime=None):
    """
    SpatHy netCDF4 format output file initialization
    IN:
        ID -catchment id as str
        fname - filename
        lat0, lon0 - latitude and longitude
        dlat - nr grid cells in lat
        dlon - nr grid cells in lon
        dtime - nr timesteps, dtime=None --> unlimited
    OUT:
        ncf - netCDF file handle. Initializes data
        ff - netCDF filename incl. path
    LAST EDIT 05.10.2018 / Samuli
    """

    from netCDF4 import Dataset #, date2num, num2date
    from datetime import datetime

    print('**** creating SpaFHy netCDF4 file: ' + fname + ' ****')
    
    # create dataset & dimensions
    ncf = Dataset(fname, 'w')
    ncf.description = 'SpatHy results. Catchment : ' + str(ID)
    ncf.history = 'created ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    ncf.source = 'SpaFHy v.1.0'

    ncf.createDimension('dtime', dtime)
    ncf.createDimension('dlon', dlon)
    ncf.createDimension('dlat', dlat)

    # create variables into base and groups 'forc','eval','cpy','bu','top'
    # call as createVariable(varname,type,(dimensions))
    time = ncf.createVariable('time', 'f8', ('dtime',))
    time.units = "days since 0001-01-01 00:00:00.0"
    time.calendar = 'standard'

    lat = ncf.createVariable('lat', 'f4', ('dlat',))
    lat.units = 'ETRS-TM35FIN'
    lon = ncf.createVariable('lon', 'f4', ('dlon',))
    lon.units = 'ETRS-TM35FIN'

    lon[:] = lon0
    lat[:] = lat0
    
    # CanopyGrid outputs
    LAI = ncf.createVariable('/cpy/LAI', 'f4', ('dtime', 'dlat', 'dlon',))
    LAI.units = 'gridcell leaf-area index [m2m-2]'
    W = ncf.createVariable('/cpy/W', 'f4', ('dtime', 'dlat', 'dlon',))
    W.units = 'canopy storage [mm]'
    SWE = ncf.createVariable('/cpy/SWE', 'f4', ('dtime', 'dlat', 'dlon',))
    SWE.units = 'snow water equiv. [mm]'
    Trfall = ncf.createVariable('/cpy/Trfall', 'f4', ('dtime', 'dlat', 'dlon',))
    Trfall.units = 'throughfall [mm]'
    Inter = ncf.createVariable('/cpy/Inter', 'f4', ('dtime', 'dlat', 'dlon',))
    Inter.units = 'interception [mm]'
    Potinf = ncf.createVariable('/cpy/Potinf', 'f4', ('dtime', 'dlat', 'dlon',))
    Potinf.units = 'pot. infiltration [mm]'
    ET = ncf.createVariable('/cpy/ET', 'f4', ('dtime', 'dlat', 'dlon',))
    ET.units = 'dry-canopy et. [mm]'
    Transpi = ncf.createVariable('/cpy/Transpi', 'f4', ('dtime', 'dlat', 'dlon',))
    Transpi.units = 'transpiration [mm]'
    Efloor = ncf.createVariable('/cpy/Efloor', 'f4', ('dtime', 'dlat', 'dlon',))
    Efloor.units = 'forest floor evap. [mm]'
    Evap = ncf.createVariable('/cpy/Evap', 'f4', ('dtime', 'dlat', 'dlon',))
    Evap.units = 'interception evap. [mm]'
    Mbe = ncf.createVariable('/cpy/Mbe', 'f4', ('dtime', 'dlat', 'dlon',))
    Mbe.units = 'mass-balance error [mm]'

    # BucketGrid outputs
    Wliq = ncf.createVariable('/bu/Wliq', 'f4', ('dtime', 'dlat', 'dlon',))
    Wliq.units = 'root zone vol. water cont. [m3m-3]'
    Wliq_top = ncf.createVariable('/bu/Wliq_top', 'f4', ('dtime', 'dlat', 'dlon',))
    Wliq_top.units = 'org. layer vol. water cont. [m3m-3]'
    PondSto = ncf.createVariable('/bu/PondSto', 'f4', ('dtime', 'dlat', 'dlon',))
    PondSto.units = 'pond storage [m]'
    Infil = ncf.createVariable('/bu/Infil', 'f4', ('dtime', 'dlat', 'dlon',))
    Infil.units = 'infiltration [m]'
    Drain = ncf.createVariable('/bu/Drain', 'f4', ('dtime', 'dlat', 'dlon',))
    Drain.units = 'drainage [m]'
    Roff = ncf.createVariable('/bu/Roff', 'f4', ('dtime', 'dlat', 'dlon',))
    Roff.units = 'surface runoff [m]'
    Retflow = ncf.createVariable('/bu/Retflow', 'f4', ('dtime', 'dlat', 'dlon',))
    Retflow.units = 'returnflow [m]'
    Mbe = ncf.createVariable('/bu/Mbe', 'f4', ('dtime', 'dlat', 'dlon',))
    Mbe.units = 'mass-balance error [m]'

    # topmodel outputs
    Qt = ncf.createVariable('/top/Qt', 'f4', ('dtime',))
    Qt.units = 'streamflow[m]'
    Qb = ncf.createVariable('/top/Qb', 'f4', ('dtime',))
    Qb.units = 'baseflow [m]'
    Qr = ncf.createVariable('/top/Qr', 'f4', ('dtime',))
    Qr.units = 'returnflow [m]'
    Qs = ncf.createVariable('/top/Qs', 'f4', ('dtime',))
    Qs.units = 'surface runoff [m]'
    R = ncf.createVariable('/top/R', 'f4', ('dtime',))
    R.units = 'average recharge [m]'
    S = ncf.createVariable('/top/S', 'f4', ('dtime',))
    S.units = 'average sat. deficit [m]'
    
    fsat = ncf.createVariable('/top/fsat', 'f4', ('dtime',))
    fsat.units = 'saturated area fraction [-]'
    Sloc = ncf.createVariable('/top/Sloc', 'f4', ('dtime','dlat','dlon',))
    Sloc.units = 'local sat. deficit [m]'
    
    print('**** netCDF4 file created ****')
    
    return ncf, fname
