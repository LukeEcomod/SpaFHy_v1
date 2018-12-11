# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 16:18:57 2016

@author: slauniai


"""
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import configparser
import timeit

from canopygrid import CanopyGrid
from bucketgrid_2 import BucketGrid
from topmodel import Topmodel_Homogenous as Topmodel

from iotools import create_catchment, read_FMI_weather, read_SVE_runoff

spathy_path = os.path.join('c:', 'c:\Repositories\spathy')
eps = np.finfo(float).eps  # machine epsilon

""" ************** SpatHy ************************************

Simple spatial hydrology and catchment water balance model.

CONSISTS OF THREE CLASSES, imported from spathy.modules.: 
    CanopyGrid - vegetation and snowpack water storages and flows
    BucketGrid - topsoil bucket model (root zone / topsoil water storage)
    Topmodel - integration to catchment scale using Topmodel -concept
HELPER FUNCTIONS:
    in iotools
 
MAIN PROGRAM:   
    spathy_driver is main program, call it as
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
        -make soil hydrologic properties (soiltypes.ini) more realistic
        -kasvupaikkatyyppi (multi-NFI) --> soil properties
        -add soil frost model, simplest would be Stefan equation with coefficients modified based on snow insulation
          --> we need snow density algorithm: SWE <-----> depth
    Topmodel:
        -think of definging 'relative m & to grids' (soil-type & elevation-dependent?) and calibrate 'catchment averages'
        -topmodel gives 'saturated zone storage deficit in [m]'. This can be converted to gwl proxy (?) if: 
        local water retention characteristics are known & hydrostatic equilibrium assumes. 
        Look which water retention model was analytically integrable (Campbell, brooks-corey?)
    
    Graphics and analysis of results:
        -make ready functions


(C) Samuli Launiainen 10/2016-->    

VERSION 08.02.2018

"""


def spathy_driver(setupfile, catchment_id=None, ncf=False, cpy_outputs=False, 
                  bu_outputs=False, top_outputs=False, flatten=False):
    """ 
    ******************** spathy_driver **********************
    
    Normal Spathy run without parameter optimization
    1) reads parameters from 'spathy.ini' (setupfile) 
    2) reads GIS-data, here predefined format for Seurantaverkko -cathcments
    3) reads forcing (defined in 'spathy.ini')
    4) optionally reads runoff measurements (file defined in 'spathy.ini')
    
    5) creates CanopyGrid (cpy), BucketGrid (bu) and Topmodel (top) -objects  within Spathy-object (spa) and temporary outputs
    6) creates netCDF -file for outputs if 'ncf' = True
    
    7) loops through forging timesteps and sequentially solves:
    i) aboveground water budget (cpy)
    ii) top and root zone soil water budget (bu),
    iii) catchment water budget, saturated areas, returnflow and streamflow generation (top)
    
    8) returns following outputs:
        spa - spathy object
        outf - filepath to output netCDF file.
    
    TODO:
        read_Runoff_data(args)
        netCDF writing in 10/30day steps
        implement calibration option - current structure supports that
        netCDF file: open and data access, graphics etc. netCDF can be accessed as dict or DataFrame.
    
    IN:
        setupfile - path to ini-file
        catchment_id - id of catchment, overrides that in ini-file
        ncf - True saves outputs to netCDF-file
        cpy_outputs, bu_, top_ - True saves cpy, bu and top outputs to lists within each object. 
            Use only for testing, memory issue.
        flatten - True flattens 2d arrys to 1d array containing only cells inside catchment. not working with current
                netcdf output file so for calibration only.
    OUT:
        spa - spathy object
        outf - filepath to netCDF-file. if ncf=False, returns None
    """
    setupfile = unicode(os.path.join(spathy_path, 'ini', setupfile))
    start_time = timeit.default_timer()

    """read parameter file into dicts"""
    pgen, pcpy, pbu, ptop = read_setup(setupfile)

    # full path to soil_file
    pgen['soil_file'] = unicode(os.path.join(spathy_path, pgen['soil_file']))

    # if given, override cathcment_id of setupfile
    if catchment_id:
        pgen['catchment_id'] = catchment_id

    gisdata = create_catchment(pgen['catchment_id'], fpath=pgen['gis_folder'],
                               plotgrids=False, plotdistr=False)

    """ greate SpatHy object """
    if ncf:
        flatten=True

    spa = SpatHy(pgen, pcpy, pbu, ptop, gisdata, cpy_outputs=cpy_outputs,
                 bu_outputs=bu_outputs, top_outputs=top_outputs, flatten=flatten)
    Nsteps = spa.Nsteps

    """ create netCDF output file """
    if ncf:
        ncf, outf = initialize_netCDF(spa.id, spa.GisData, spa.FORC,
                                      fpath=spa.pgen['output_folder'],
                                      fname=spa.pgen['ncf_file'])
        print outf
    else:
        outf = None
                                
    """ initialize 3D arrays for temporary outputs """
    r, c = np.shape(spa.cmask)

    #3d array indexing: dim1=time, dim2=rows(lat), dim3=cols(lon). W[1,:,:] --> grid at 1st timestep. 
    #W[:,30,30]=timeserie at cell with index 30,30

    """ ----- MAIN CALCULATION LOOP ----- """
    print('check version')
    print('Reading files and init objects [s]: ', timeit.default_timer() - start_time )
    del start_time
    start_time = timeit.default_timer()

    print '******* Running Spathy ********'
    spa._run(0, Nsteps, calibr=False, ncf=ncf)

    print('Loops total [s]: ', timeit.default_timer() - start_time)
    print '********* done *********'

    return spa, outf


"""
******************************************************************************
            ----------- SpatHy model class --------
******************************************************************************
"""


class SpatHy():
    """
    SpatHy model class
    """
    def __init__(self, pgen, pcpy, pbu, ptop, gisdata, ave_outputs=False,
                 cpy_outputs=False, bu_outputs=False, top_outputs=False, flatten=False):

        self.dt = pgen['dt']  # s
        self.id = pgen['catchment_id']
        self.spinup_end = pgen['spinup_end']
        self.pgen = pgen

        """ read forcing data and catchment runoff file """
        FORC = read_FMI_weather(pgen['catchment_id'],
                                pgen['start_date'],
                                pgen['end_date'],
                                sourcefile=pgen['forcing_file'])
        FORC['Prec'] = FORC['Prec'] / self.dt  # mms-1

        self.FORC = FORC
        self.Nsteps = len(FORC)

        # read runoff measurements
        if pgen['runoff_file'] is not '':
            self.Qmeas = read_SVE_runoff(pgen['catchment_id'],
                                         pgen['start_date'], pgen['end_date'], pgen['runoff_file'])
        else:
            self.Qmeas = None

        self.GisData = gisdata
        self.cmask = self.GisData['cmask']
        self.gridshape = np.shape(self.cmask)
        
        cmask= self.cmask.copy()
        lai_conif = gisdata['LAI_conif'].copy()
        lai_decid = gisdata['LAI_decid'].copy()
        hc = gisdata['hc'].copy()
        cf = gisdata['cf'].copy()
        sdata = gisdata['soil'].copy()
        flowacc = gisdata['flowacc'].copy()
        slope = gisdata['slope'].copy()        
        
        """
        flatten=True omits cells outside catchment
        """
        if flatten:
            ix = np.where(np.isfinite(cmask))
            cmask = cmask[ix].copy()
            lai_conif = lai_conif[ix].copy()
            lai_decid = lai_decid[ix].copy()
            hc = hc[ix].copy()
            cf = cf[ix].copy()
            sdata = sdata[ix].copy()
            flowacc = flowacc[ix].copy()
            slope = slope[ix].copy()
            
            self.ix = ix  # indices
         
        """--- initialize CanopyGrid and BucketGrid ---"""
        if pgen['spatial_cpy'] == 'no':  # spatially constant stand properties
            self.cpy = CanopyGrid(pcpy, cmask=cmask, outputs=cpy_outputs)
        else:
            self.cpy = CanopyGrid(pcpy, lai_conif=lai_conif, lai_decid = lai_decid,
                                  cf=cf, hc=hc, cmask=cmask, outputs=cpy_outputs)

        if pgen['spatial_soil'] == 'no':  # spatially constant soil properties
            self.bu = BucketGrid(pbu, cmask=cmask, outputs=bu_outputs)
        else:
            self.bu = BucketGrid(pbu, soiltypefile=pgen['soil_file'], sdata=sdata, outputs=bu_outputs)

        """ --- initialize homogenous topmodel --- """
        self.top=Topmodel(ptop, self.GisData['cellsize']**2, cmask,
                          flowacc, slope, outputs=top_outputs)
        
        if ave_outputs:
            self.results = {'SWE': np.ones(self.Nsteps)*np.NaN, 'ET': np.ones(self.Nsteps)*np.NaN,
                            'Tr': np.ones(self.Nsteps)*np.NaN, 'Prec': np.ones(self.Nsteps)*np.NaN,
                            'Ef': np.ones(self.Nsteps)*np.NaN, 'Evap': np.ones(self.Nsteps)*np.NaN,
                            'Qt': np.ones(self.Nsteps)*np.NaN, 'S': np.ones(self.Nsteps)*np.NaN, 
                            'fsat': np.ones(self.Nsteps)*np.NaN, 'Wliq': np.ones(self.Nsteps)*np.NaN,
                            'Rg': np.ones(self.Nsteps)*np.NaN, 'Ta': np.ones(self.Nsteps)*np.NaN, 
                            'VPD': np.ones(self.Nsteps)*np.NaN, 'Drain': np.ones(self.Nsteps)*np.NaN,
                            'time': self.FORC.index, 'Qmeas': self.Qmeas
                            }

    def _run(self, fstep, Nsteps, calibr=False, ncf=False, sve=False):
        """ 
        Runs Spathy
        IN:
            fstep - index of starting point [int]
            Nsteps - number of timesteps [int]
            calibr - set True for parameter optimization, returns Qmod
            ncf - netCDF -file handle, for outputs
        OUT:
            res - modeled streamflow at catchment outlet [mm/d]
        """
        dt = self.dt

        # for calibration run, return res
        if calibr:
             res={'Qm':[None]*self.Nsteps}  # 'RR':[None]*self.Nsteps, 'ET':[None]*self.Nsteps, 'Inter':[None]*self.Nsteps, 'Mbet':[None]*self.Nsteps}
             print('M', self.top.M, 'to', self.top.To)
        
        RR = 0.0 # initial value for recharge [m]
        for k in range(fstep, fstep + Nsteps):
            #print 'k=' + str(k)
            # forcing
            doy = self.FORC['doy'].iloc[k]
            ta = self.FORC['T'].iloc[k]
            vpd = self.FORC['VPD'].iloc[k] + eps
            rg = self.FORC['Rg'].iloc[k]
            par = self.FORC['Par'].iloc[k] + eps
            prec = self.FORC['Prec'].iloc[k]
            co2 = self.FORC['CO2'].iloc[k]
            u = 2.0            

            # beta0 = self.bu.WatStoTop / self.bu.MaxStoTop
            # beta0 = self.bu.relcond

            # run Topmodel, take recharge RR from prev. timestep
            # mean baseflow [m], none, returnflow grid [m], sat.area [-]
            qb, _, qr, fsat = self.top.run_timestep(RR)

            # run CanopyGrid
            potinf, trfall, interc, evap, et, transpi, efloor, mbe = \
                self.cpy.run_timestep(doy, dt, ta, prec, rg, par, vpd, U=u, CO2=co2,
                                      beta=self.bu.Ree, Rew=self.bu.Rew, P=101300.0)

            # run BucketGrid water balance

            infi, infi_ex, drain, tr, eva, mbes, rflow_to_rootzone = self.bu.watbal(dt=dt, rr=1e-3*potinf, tr=1e-3*transpi,
                                                           evap=1e-3*efloor, retflow=qr)

            # catchment average [m per unit area] saturation excess --> goes to stream
            # as surface runoff
            qs = np.nansum(infi_ex)*self.top.CellArea / self.top.CatchmentArea

            # catchment average ground water recharge [m per unit area]
            RR = np.nansum(drain)*self.top.CellArea / self.top.CatchmentArea

            """ outputs """
            
            if hasattr(self.top, 'results'):
                self.top.results['Qt'].append(qb + qs + eps)  # total runoff

            if ncf:
                # writes to netCDF -file at every timestep; bit slow - should
                # accumulate into temporary variables and save every 10 days?                

                # canopygrid
                # ncf['cpy']['W'][k,:,:] = self._to_grid(self.cpy.W)
                ncf['cpy']['SWE'][k,:,:] = self._to_grid(self.cpy.SWE)
                # ncf['cpy']['Trfall'][k,:,:] = self._to_grid(trfall) 
                # ncf['cpy']['Potinf'][k,:,:] = self._to_grid(potinf)
                ncf['cpy']['ET'][k,:,:] = self._to_grid(et)
                ncf['cpy']['Transpi'][k,:,:] = self._to_grid(transpi)
                ncf['cpy']['Efloor'][k,:,:] = self._to_grid(efloor)            
                ncf['cpy']['Evap'][k,:,:] = self._to_grid(evap)
                # ncf['cpy']['Inter'][k,:,:] = self._to_grid(interc)
                # ncf['cpy']['Mbe'][k,:,:] = self._to_grid(mbe)              

                # bucketgrid
                ncf['bu']['Drain'][k,:,:] = 1e3 * self._to_grid(drain)
                ncf['bu']['Infil'][k,:,:] = 1e3 * self._to_grid(infi)
                ncf['bu']['Wliq'][k,:,:] = self._to_grid(self.bu.Wliq)
                ncf['bu']['Wliqtop'][k,:,:] = self._to_grid(self.bu.Wliq_top)
                # root zone rescahrge from return flow. 
                ncf['bu']['RetFlow'][k,:,:] = 1e3 * self._to_grid(rflow_to_rootzone)
                # ncf['bu']['Wsto'][k,:,:] = self._to_grid(self.bu.WatSto)
                # ncf['bu']['PondSto'][k,:,:] = self._to_grid(self.bu.PondSto)
                # ncf['bu']['Mbe'][k,:,:] = self._to_grid(mbes)              

                # topmodel
                ncf['top']['Qb'][k] = qb
                ncf['top']['Qs'][k] = qs
                ncf['top']['Qt'][k] = qb + qs  # total runoff
                ncf['top']['R'][k] = RR
                ncf['top']['fsat'][k] = fsat
                ncf['top']['S'][k] = self.top.S

            if calibr:  # calibration run, return only streamflow
                res['Qm'][k] = 1e3*(qb + qs)
            
            if hasattr(self, 'results'):
                self.results['SWE'][k] = np.nanmean(self.cpy.SWE)
                self.results['Wliq'][k] = np.nanmean(self.bu.Wliq)
                self.results['ET'][k] = np.nanmean(transpi + efloor + evap)
                self.results['Tr'][k] = np.nanmean(transpi)
                self.results['Ef'][k] = np.nanmean(efloor)
                self.results['Evap'][k] = np.nanmean(evap)
                self.results['Drain'][k] =  1e3*RR
                self.results['Qt'][k] = 1e3*(qb + qs)
                self.results['S'][k] = self.top.S
                self.results['fsat'][k] = fsat
                self.results['Prec'][k] = prec*self.dt
                self.results['Rg'][k] = rg
                self.results['Ta'][k] = ta
                self.results['VPD'][k] = vpd

        # end of time loop
        if hasattr(self, 'results'):
            res = pd.DataFrame.from_dict(self.results)
            res.index = self.FORC.index
            res = res[res.index > self.spinup_end]
            self.results = res
        if ncf:
            ncf.close()  # close netCDF-file

        if calibr:  # return streamflow in dataframe
            res = pd.DataFrame.from_dict(res)
            res.index = self.FORC.index

            return res

    def _to_grid(self, x):
        """
        converts variable x back to original grid for NetCDF outputs
        """
        a = np.full(self.gridshape, np.NaN)
        a[self.ix] = x
        return a

""" Functions for data input and output """

def read_setup(inifile):
    """
    reads Spathy.ini parameter file into pp dict
    """
    inifile = os.path.join(spathy_path, inifile)
    
    cfg = configparser.ConfigParser()
    cfg.read(inifile)

    pp = {}
    for s in cfg.sections():
        section = s.encode('ascii', 'ignore')
        pp[section] = {}
        for k, v in cfg.items(section):
            key = k.encode('ascii', 'ignore')
            val = v.encode('ascii', 'ignore')
            if section == 'General':  # 'general' section
                pp[section][key] = val
            else:
                pp[section][key] = float(val)
    pp['General']['dt'] = float(pp['General']['dt'])

    pgen = pp['General']
    pcpy = pp['CanopyGrid']
    pbu = pp['BucketGrid']
    ptop = pp['Topmodel']

    return pgen, pcpy, pbu, ptop

#def initialize_netCDF(ID, gis, forc, roff=None, fpath='results', fname=None):
#    """
#    SpatHy netCDF4 format output file initialization
#    IN:
#        ID -catchment id as int or str
#        gis - GisData dict
#        forc - forcing data (pd.dataframe)
#        roff - measured runoff (pd.Series)
#        fpath - path for saving results
#        fname - filename
#    OUT:
#        ncf - netCDF file handle. Initializes data
#        ff - netCDF filename incl. path
#    LAST EDIT 1.11.
#    """
#    from netCDF4 import Dataset, date2num  # , num2date
#    from datetime import datetime
#
#    # dimensions
#    dlat, dlon = np.shape(gis['cmask'])
#    dtime = None
#
#    if fname:
#        ff = os.path.join(spathy_path, fpath, fname)
#        print ff
#    else:
#        ff = os.path.join(spathy_path, fpath, 'Spathy_ch' + str(ID) + '.nc')
#
#    # create dataset & dimensions
#    ncf = Dataset(ff, 'w')
#    ncf.description = 'SpatHy results. Catchment : ' + str(ID)
#    ncf.history = 'created ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S')
#    ncf.source = 'SpatHy -model v.0.99'
#
#    ncf.createDimension('dtime', dtime)
#    ncf.createDimension('dlon', dlon)
#    ncf.createDimension('dlat', dlat)
#
#    # create variables into base and groups 'forc','eval','cpy','bu','top'
#    # call as createVariable(varname,type,(dimensions))
#    time = ncf.createVariable('time', 'f8', ('dtime',))
#    time.units = "days since 0001-01-01 00:00:00.0"
#    time.calendar = 'standard'
#
#    lat = ncf.createVariable('lat', 'f4', ('dlat',))
#    lat.units = 'ETRS-TM35FIN'
#    lon = ncf.createVariable('lon', 'f4', ('dlon',))
#    lon.units = 'ETRS-TM35FIN'
#
#    tvec = [k.to_datetime() for k in forc.index]
#    time[:] = date2num(tvec, units=time.units, calendar=time.calendar)
#    lon[:] = gis['lon0']
#    lat[:] = gis['lat0']
#
#    # evaluation data
#    Roff = ncf.createVariable('/eval/Roff', 'f4', ('dtime',))
#    Roff.units = 'meas. streamflow [mm]'
#    if roff:
#        Roff[:] = roff.values
#
#    # CanopyGrid outputs
#    W = ncf.createVariable('/cpy/W', 'f4', ('dtime', 'dlat', 'dlon',))
#    W.units = 'canopy storage [mm]'
#    SWE = ncf.createVariable('/cpy/SWE', 'f4', ('dtime', 'dlat', 'dlon',))
#    SWE.units = 'snow water equiv. [mm]'
#    Trfall = ncf.createVariable('/cpy/Trfall', 'f4', ('dtime', 'dlat', 'dlon',))
#    Trfall.units = 'throughfall [mm]'
#    Inter = ncf.createVariable('/cpy/Inter', 'f4', ('dtime', 'dlat', 'dlon',))
#    Inter.units = 'interception [mm]'
#    Potinf = ncf.createVariable('/cpy/Potinf', 'f4', ('dtime', 'dlat', 'dlon',))
#    Potinf.units = 'pot. infiltration [mm]'
#    ET = ncf.createVariable('/cpy/ET', 'f4', ('dtime', 'dlat', 'dlon',))
#    ET.units = 'dry-canopy et. [mm]'
#    Transpi = ncf.createVariable('/cpy/Transpi', 'f4', ('dtime', 'dlat', 'dlon',))
#    Transpi.units = 'transpiration [mm]'
#    Efloor = ncf.createVariable('/cpy/Efloor', 'f4', ('dtime', 'dlat', 'dlon',))
#    Efloor.units = 'forest floor evap. [mm]'
#    Evap = ncf.createVariable('/cpy/Evap', 'f4', ('dtime', 'dlat', 'dlon',))
#    Evap.units = 'interception evap. [mm]'
#    Mbe = ncf.createVariable('/cpy/Mbe', 'f4', ('dtime', 'dlat', 'dlon',))
#    Mbe.units = 'mass-balance error [mm]'
#
#    # BucketGrid outputs
#    Wliq = ncf.createVariable('/bu/Wliq', 'f4', ('dtime', 'dlat', 'dlon',))
#    Wliq.units = 'vol. water cont. [m3m-3]'
#    Wsto = ncf.createVariable('/bu/Wsto', 'f4', ('dtime', 'dlat', 'dlon',))
#    Wsto.units = 'water storage [mm]'
#    PondSto = ncf.createVariable('/bu/PondSto', 'f4', ('dtime', 'dlat', 'dlon',))
#    PondSto.units = 'pond storage [mm]'
#    Infil = ncf.createVariable('/bu/Infil', 'f4', ('dtime', 'dlat', 'dlon',))
#    Infil.units = 'infiltration [mm]'
#    Drain = ncf.createVariable('/bu/Drain', 'f4', ('dtime', 'dlat', 'dlon',))
#    Drain.units = 'drainage [mm]'
#    Mbe = ncf.createVariable('/bu/Mbe', 'f4', ('dtime', 'dlat', 'dlon',))
#    Mbe.units = 'mass-balance error [mm]'
#
#    # topmodel outputs
#    Qt = ncf.createVariable('/top/Qt', 'f4', ('dtime',))
#    Qt.units = 'streamflow[m]'
#    Qb = ncf.createVariable('/top/Qb', 'f4', ('dtime',))
#    Qb.units = 'baseflow [m]'
#    Qr = ncf.createVariable('/top/Qr', 'f4', ('dtime',))
#    Qr.units = 'returnflow [m]'
#    Qs = ncf.createVariable('/top/Qs', 'f4', ('dtime',))
#    Qs.units = 'surface runoff [m]'
#    R = ncf.createVariable('/top/R', 'f4', ('dtime',))
#    R.units = 'average recharge [m]'
#    S = ncf.createVariable('/top/S', 'f4', ('dtime',))
#    S.units = 'average sat. deficit [m]'
#    fsat = ncf.createVariable('/top/fsat', 'f4', ('dtime',))
#    fsat.units = 'saturated area fraction [-]'
#
#    return ncf, ff


def initialize_netCDF(ID, gis, forc, roff=None, fpath='results', fname=None):
    """
    SpatHy netCDF4 format output file initialization
    IN:
        ID -catchment id as int or str
        gis - GisData dict
        forc - forcing data (pd.dataframe)
        roff - measured runoff (pd.Series)
        fpath - path for saving results
        fname - filename
    OUT:
        ncf - netCDF file handle. Initializes data
        ff - netCDF filename incl. path
    LAST EDIT 1.11.
    """
    from netCDF4 import Dataset, date2num  # , num2date
    from datetime import datetime

    # dimensions
    dlat, dlon = np.shape(gis['cmask'])
    dtime = None
    print fpath, fname
    if fname:
        ff = os.path.join(spathy_path, fpath, fname)
        print ff
    else:
        ff = os.path.join(spathy_path, fpath, 'Spathy_ch' + str(ID) + '.nc')

    # create dataset & dimensions
    ncf = Dataset(ff, 'w')
    ncf.description = 'SpatHy results. Catchment : ' + str(ID)
    ncf.history = 'created ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    ncf.source = 'SpatHy -model v.0.99'

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

    tvec = [k.to_datetime() for k in forc.index]
    time[:] = date2num(tvec, units=time.units, calendar=time.calendar)
    lon[:] = gis['lon0']
    lat[:] = gis['lat0']

    # evaluation data
    Roff = ncf.createVariable('/eval/Roff', 'f4', ('dtime',))
    Roff.units = 'meas. streamflow [mm]'
    if roff:
        Roff[:] = roff.values

    # CanopyGrid outputs
    # W = ncf.createVariable('/cpy/W', 'f4', ('dtime', 'dlat', 'dlon',))
    # W.units = 'canopy storage [mm]'
    SWE = ncf.createVariable('/cpy/SWE', 'f4', ('dtime', 'dlat', 'dlon',))
    SWE.units = 'snow water equiv. [mm]'
    # Trfall = ncf.createVariable('/cpy/Trfall', 'f4', ('dtime', 'dlat', 'dlon',))
    # Trfall.units = 'throughfall [mm]'
    # Inter = ncf.createVariable('/cpy/Inter', 'f4', ('dtime', 'dlat', 'dlon',))
    # Inter.units = 'interception [mm]'
    # Potinf = ncf.createVariable('/cpy/Potinf', 'f4', ('dtime', 'dlat', 'dlon',))
    # Potinf.units = 'pot. infiltration [mm]'
    ET = ncf.createVariable('/cpy/ET', 'f4', ('dtime', 'dlat', 'dlon',))
    ET.units = 'dry-canopy et. [mm]'
    Transpi = ncf.createVariable('/cpy/Transpi', 'f4', ('dtime', 'dlat', 'dlon',))
    Transpi.units = 'transpiration [mm]'
    Efloor = ncf.createVariable('/cpy/Efloor', 'f4', ('dtime', 'dlat', 'dlon',))
    Efloor.units = 'forest floor evap. [mm]'
    Evap = ncf.createVariable('/cpy/Evap', 'f4', ('dtime', 'dlat', 'dlon',))
    Evap.units = 'interception evap. [mm]'
    # Mbe = ncf.createVariable('/cpy/Mbe', 'f4', ('dtime', 'dlat', 'dlon',))
    # Mbe.units = 'mass-balance error [mm]'

    # BucketGrid outputs
    Wliq = ncf.createVariable('/bu/Wliq', 'f4', ('dtime', 'dlat', 'dlon',))
    Wliq.units = 'vol. water cont. [m3m-3]'
    Wliqtop = ncf.createVariable('/bu/Wliqtop', 'f4', ('dtime', 'dlat', 'dlon',))
    Wliqtop.units = ' top layer vol. water cont. [m3m-3]'
    # Wsto = ncf.createVariable('/bu/Wsto', 'f4', ('dtime', 'dlat', 'dlon',))
    # Wsto.units = 'water storage [mm]'
    # PondSto = ncf.createVariable('/bu/PondSto', 'f4', ('dtime', 'dlat', 'dlon',))
    # PondSto.units = 'pond storage [mm]'
    Infil = ncf.createVariable('/bu/Infil', 'f4', ('dtime', 'dlat', 'dlon',))
    Infil.units = 'infiltration [mm]'
    Drain = ncf.createVariable('/bu/Drain', 'f4', ('dtime', 'dlat', 'dlon',))
    Drain.units = 'drainage [mm]'
    RetFlow = ncf.createVariable('/bu/RetFlow', 'f4', ('dtime', 'dlat', 'dlon',))
    RetFlow.units = 'RetFlow [mm]'
    # Mbe = ncf.createVariable('/bu/Mbe', 'f4', ('dtime', 'dlat', 'dlon',))
    # Mbe.units = 'mass-balance error [mm]'

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

    return ncf, ff