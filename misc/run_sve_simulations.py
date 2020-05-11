# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 12:37:40 2018

@author: slauniai
"""
import os
#import spotpy
import pickle
import numpy as np
from scipy import stats
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from spathy_sve import SpatHy, initialize_netCDF, read_setup
from iotools import create_catchment

eps = np.finfo(float).eps
spathy_path = os.path.join('c:\\repositories\\spathy')
setupfile = os.path.join(spathy_path, 'ini', 'spathy_default.ini')

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
     ['22', '2005-01-01', '2015-12-31', '2005-12-31', 0.0095],  # vaha-askanjoki 1600 ha
#    [ '23','2011-01-01', '2015-12-31', '2011-12-31'],  # ylijoki 5600 ha, very large and slow
     ['24', '2005-01-01', '2015-12-31', '2005-12-31', 0.0066],  # kotioja 1800 ha
     ['25', '2005-01-01', '2015-12-31', '2005-12-31', 0.0095],  # kohisevanpuro 1070 ha
     ['26', '2005-01-01', '2015-12-31', '2005-12-31', 0.02],  # iittovuoma 1160 ha
     ['27', '2005-01-01', '2015-12-31', '2005-12-31', 0.014],  # laanioja 1362 ha
     ['28', '2013-01-01', '2015-12-31', '2013-12-31', 0.0057],  # kroopinsuo 179 ha
     ['29', '2012-01-01', '2015-12-31', '2012-12-31', 0.0089],  # surnui 71 ha, poor data quality
     ['30', '2011-01-01', '2015-12-31', '2011-12-31', 0.0064],  # pakopirtti 795 ha, uncertain catchment boundaries
     ['31', '2011-01-01', '2015-12-31', '2011-12-31', 0.0064],  # ojakorpi 33 ha
     ['32', '2011-01-01', '2015-12-31', '2011-12-31', 0.0077],  # rantainrahka 38 ha
     ['33', '2005-01-01', '2015-12-31', '2005-12-31', 0.009],  # kivipuro 54 ha
     ]

# chm = chm[0:2]

"""read default parameter file into dicts"""

pgen0, pcpy0, pbu0, ptop0 = read_setup(setupfile)

# full path to soil_file
pgen0['soil_file'] = unicode(os.path.join(spathy_path, pgen0['soil_file']))

ff = [[1.0, 1.0, 1.0, 1.0], [0.8, 0.8, 0.8, 0.8], [1.2, 1.2, 1.2, 1.2]]  # gsref, Wmax, Wmaxsnow, LAI multipliers
ss = ['base', 'lowET', 'hiET']

for k in [0, 1, 2]: # loop et-cases
    res = []
    for n in range(0, len(chm)): # loop catchments
        print(chm[n])
        pgen = pgen0.copy()
        pcpy = pcpy0.copy()
        pbu = pbu0.copy()
        ptop = ptop0.copy()
        
        # n -loop parameters
        pgen['catchment_id'] = chm[n][0]
        pgen['start_date'] = chm[n][1]
        pgen['end_date'] = chm[n][2]
        pgen['spinup_end'] = chm[n][3]
        pgen['outname'] = 'Ch' + pgen['catchment_id'] + '-' + ss[k] + '.nc'
        ptop['m'] = chm[n][4]
        
        # k -loop changes
        pcpy['gsref_conif'] *= ff[k][0]
        pcpy['gsref_decid'] *= ff[k][0]
        pcpy['wmax'] *= ff[k][1]
        pcpy['wmaxsnow'] *= ff[k][2]
        pcpy['lai_multip'] =ff[k][3]

        print(k, n, pgen['outname'])

        # run baseline simulation, return nothing
#        if k == 0:
#            spa = spathy_run_sve(cid, pgen, pcpy, pbu, ptop)
#        else:
        spa = spathy_run_sve(pgen, pcpy, pbu, ptop, ncf=False, ave_outputs=True)
        res.append(spa.results)
        del pgen, pbu, ptop
    
    # dump into pickle
    ou = os.path.join(pgen0['output_folder'], 'R-' + ss[k] + '.pkl')
    pickle.dump(res, open(ou, 'wb'))

#%%
    
def spathy_run_sve(pgen, pcpy, pbu, ptop, ncf=True, ave_outputs=True, flatten=True):
    """ 
    Spathy_driver for running sve catchments

    OUT:
        spa - spathy object
        outf - filepath to netCDF-file. if ncf=False, returns None
    """

    gisdata = create_catchment(pgen['catchment_id'], fpath=pgen['gis_folder'],
                               plotgrids=False, plotdistr=False)
    gisdata['LAI_conif'] *= pcpy['lai_multip']
    gisdata['LAI_decid'] *= pcpy['lai_multip']
    
    """ greate SpatHy object """
    spa = SpatHy(pgen, pcpy, pbu, ptop, gisdata, ave_outputs=ave_outputs, flatten=True)
    Nsteps = spa.Nsteps

    """ create netCDF output file """
    if ncf:
        ncf, _= initialize_netCDF(spa.id, spa.GisData, spa.FORC,
                                      fpath=spa.pgen['output_folder'],
                                      fname=pgen['outname'])
                                
    #3d array indexing: dim1=time, dim2=rows(lat), dim3=cols(lon). W[1,:,:] --> grid at 1st timestep. 

    """ ----- MAIN CALCULATION LOOP ----- """

    print '******* Running Spathy ********'
    spa._run(0, Nsteps, calibr=False, ncf=ncf)
    
    print '********* done *********'

    return spa
