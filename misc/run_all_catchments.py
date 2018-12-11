# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 12:37:40 2018

@author: slauniai
"""
import sys
sys.path.append(r'\repositories\SpaFHy')

import os
import pickle
import numpy as np
from scipy import stats
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import spafhy
from spafhy_parameters_default import soil_properties, parameters
#from spafhy_io import read_FMI_weather, write_AsciiGrid
from spafhy_io import create_catchment, read_FMI_weather, read_SVE_runoff

eps = np.finfo(float).eps


""" catchment data: id, start_date, end_date, spinup_end, top_m """

chm=[['1', '2013-01-01', '2015-12-31', '2013-12-31', 0.025],   # lompolojanganoja 514 ha
     ['2', '2005-01-01', '2008-12-31', '2005-12-31', 0.006],   # liuhapuro 170 ha
     ['3', '2005-01-01', '2015-12-31', '2005-12-31', 0.026],   # porkkavaara 72 ha
     ['10', '2005-01-01', '2013-12-31', '2005-12-31', 0.011],  # kelopuro 74 ha. 2014 gappy, 2015 runoff is low
     ['11', '2014-01-01', '2015-12-31', '2014-12-31', 0.012],  # hauklammenoja 137 ha
     ['13', '2014-01-01', '2015-12-31', '2014-12-31', 0.007],  # rudbacken 436 ha
     ['14', '2005-01-01', '2015-12-31', '2005-12-31', 0.007],  # paunulanpuro 154 ha
     ['16', '2005-01-01', '2015-12-31', '2005-12-31', 0.007],  # huhtisuonoja 500 ha. very flat, large fraction is drained peatlands
     ['17', '2005-01-01', '2015-12-31', '2005-12-31', 0.006],  # kesselinpuro 2100 ha
#     ['18','2011-01-01', '2015-12-31', '2011-12-31'],  # korpijoki, area 12200 ha so not suitable
     ['19', '2005-01-01', '2015-12-31', '2005-12-31', 0.006],  # pahkaoja 2344 ha
     ['20', '2005-01-01', '2015-12-31', '2005-12-31', 0.009],  # vaarajoki 1900 ha
     ['21', '2005-01-01', '2015-12-31', '2005-12-31', 0.01],  # myllypuro 1053 ha
#     ['22', '2005-01-01', '2015-12-31', '2005-12-31', 0.0095],  # vaha-askanjoki 1600 ha
#    [ '23','2011-01-01', '2015-12-31', '2011-12-31'],  # ylijoki 5600 ha, very large and slow
     ['24', '2005-01-01', '2015-12-31', '2005-12-31', 0.0066],  # kotioja 1800 ha
     ['25', '2005-01-01', '2015-12-31', '2005-12-31', 0.0095],  # kohisevanpuro 1070 ha
     ['26', '2005-01-01', '2015-12-31', '2005-12-31', 0.02],  # iittovuoma 1160 ha
     ['27', '2005-01-01', '2015-12-31', '2005-12-31', 0.014],  # laanioja 1362 ha
     ['28', '2013-01-01', '2015-12-31', '2013-12-31', 0.0057],  # kroopinsuo 179 ha
     ['29', '2012-01-01', '2015-12-31', '2012-12-31', 0.0089],  # surnui 71 ha, poor data quality
#    ['30', '2011-01-01', '2015-12-31', '2011-12-31', 0.0064],  # pakopirtti 795 ha, uncertain catchment boundaries
     ['31', '2011-01-01', '2015-12-31', '2011-12-31', 0.0064],  # ojakorpi 33 ha
     ['32', '2011-01-01', '2015-12-31', '2011-12-31', 0.0077],  # rantainrahka 38 ha
     ['33', '2005-01-01', '2015-12-31', '2005-12-31', 0.009],  # kivipuro 54 ha
     ]

# chm = chm[1:3]

# multipliers for ET scenarios
#g1, wmax, wmaxshow, lai
ff = [[1.0, 1.0, 1.0, 1.0], [0.8, 0.8, 0.8, 0.8], [1.2, 1.2, 1.2, 1.2]]  
ss = ['base', 'lowET', 'hiET']


for k in [0, 1, 2]: #, 2]: # loop et-cases
    results = []
    for n in range(0, len(chm)): # loop catchments
        # update parameters
        print('Scenario: ' + str(k) + ' C: ' + chm[n][0])
        
        # default parameters
        pgen, pcpy, pbu, ptop = parameters()
        psoil = soil_properties()
                
        # n -loop parameters
        pgen['catchment_id'] = chm[n][0]
        pgen['start_date'] = chm[n][1]
        pgen['end_date'] = chm[n][2]
        pgen['spinup_end'] = chm[n][3]
        pgen['ncf_file'] = 'Ch' + pgen['catchment_id'] + '-' + ss[k] + '.nc'
        ptop['m'] = chm[n][4]
        
        # k -loop changes
        pcpy['physpara']['g1_conif'] *= ff[k][0]
        pcpy['physpara']['g1_decid'] *= ff[k][0]
        pcpy['interc']['wmax'] *= ff[k][1]
        pcpy['interc']['wmaxsnow'] *= ff[k][2]
        
        # load gis data
        gisdata = create_catchment(pgen['catchment_id'], fpath=pgen['gis_folder'],
                           plotgrids=False, plotdistr=False)
        
        gisdata['LAI_conif'] *= ff[k][3]
        gisdata['LAI_decid'] *= ff[k][3]
        
        # initialize spafhy
        spa = spafhy.initialize(pgen, pcpy, pbu, ptop, psoil, gisdata, cpy_outputs=False, 
                                bu_outputs=False, top_outputs=False, flatten=True)
        # print('LAI', np.nanmean(spa.cpy.LAI))
        # create netCDF output file
        #dlat, dlon = np.shape(spa.GisData['cmask'])
        
        #ncf, ncf_file = spafhy.initialize_netCDF(ID=spa.id, fname=spa.ncf_file, lat0=spa.GisData['lat0'], 
        #                                         lon0=spa.GisData['lon0'], dlat=dlat, dlon=dlon, dtime=None)
        
        # read forcing data
        FORC = read_FMI_weather(pgen['catchment_id'],
                                pgen['start_date'],
                                pgen['end_date'],
                                sourcefile=pgen['forcing_file'])
        FORC['Prec'] = FORC['Prec'] / spa.dt  # mms-1
        FORC['U'] = 2.0 # use constant wind speed ms-1
        
        Nsteps = len(FORC)
        Nspin = np.where(FORC.index == pgen['spinup_end'])[0][0]
        
        # read catchment runoff data
        Qmeas = read_SVE_runoff(pgen['catchment_id'],
                                pgen['start_date'],
                                pgen['end_date'],
                                sourcefile=pgen['runoff_file'])
        
        Qmeas = Qmeas[(Qmeas.index > pgen['spinup_end'])]

        # run spinup
        for j in range(0, Nspin):
            # print('step: ' + str(j))
            forc= FORC[['doy', 'Rg', 'Par', 'T', 'Prec', 'VPD', 'CO2','U']].iloc[j]
            
            spa.run_timestep(forc, ncf=False, ave_flx=False)
            
        
        # run SpaFHy and append results
        N = Nsteps - Nspin -1
        
        res =  {
                'ET': np.zeros(N), 'E': np.zeros(N), 'Ef': np.zeros(N),
                'Tr': np.zeros(N), 'SWE': np.zeros(N), 'Drain': np.zeros(N),
                'Qt': np.zeros(N), 'S': np.zeros(N), 'fsat': np.zeros(N),
                'Prec': np.zeros(N), 'Rg': np.zeros(N), 'Ta': np.zeros(N),
                'VPD': np.zeros(N)
                }
        
        kk = 0
        for j in range(Nspin+1, Nsteps):
            # print('step: ' + str(j))
            forc= FORC[['doy', 'Rg', 'Par', 'T', 'Prec', 'VPD', 'CO2','U']].iloc[j]
            
            flx = spa.run_timestep(forc, ncf=False, ave_flx=True)
            
            for m in res.keys():
                res[m][kk] = flx[m]
            kk += 1

        res['Qmeas'] = Qmeas
        # res['FORC'] = FORC[(FORC.index > pgen['spinup_end'])]
        
        res = pd.DataFrame(data=res, columns=res.keys(), index=Qmeas.index)
        results.append(res)
        
        del pgen, pbu, ptop, spa
    
    # dump into pickle
    ou = os.path.join( 'R-' + ss[k] + '.pkl')
    pickle.dump(results, open(ou, 'wb'))


#%%
    
#def spathy_run_sve(pgen, pcpy, pbu, ptop, ncf=True, flatten=True):
#    """ 
#    Spathy_driver for running sve catchments
#
#    OUT:
#        spa - spathy object
#        outf - filepath to netCDF-file. if ncf=False, returns None
#    """
#
#    gisdata = create_catchment(pgen['catchment_id'], fpath=pgen['gis_folder'],
#                               plotgrids=False, plotdistr=False)
#    gisdata['LAI_conif'] *= pcpy['lai_multip']
#    gisdata['LAI_decid'] *= pcpy['lai_multip']
#    
#    """ greate SpatHy object """
#    spa = SpatHy(pgen, pcpy, pbu, ptop, gisdata, ave_outputs=ave_outputs, flatten=True)
#    Nsteps = spa.Nsteps
#
#    """ create netCDF output file """
#    if ncf:
#        ncf, _= initialize_netCDF(spa.id, spa.GisData, spa.FORC,
#                                      fpath=spa.pgen['output_folder'],
#                                      fname=pgen['outname'])
#                                
#    #3d array indexing: dim1=time, dim2=rows(lat), dim3=cols(lon). W[1,:,:] --> grid at 1st timestep. 
#
#    """ ----- MAIN CALCULATION LOOP ----- """
#
#    print '******* Running Spathy ********'
#    spa._run(0, Nsteps, calibr=False, ncf=ncf)
#    
#    print '********* done *********'
#
#    return spa
