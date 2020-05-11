# -*- coding: utf-8 -*-
"""
Created on Tue Oct 02 14:47:14 2018

@author: slauniai

DEMONSTRATES HOW TO RUN SPAFHY FOR A SINGLE CATCHMENT AND SAVE RESULTS TO NCF.
Use data from C3 (Porkkavaara), a small catchment in Eastern Finland

DEMONSTRATES SPAFHY SETUP,  RUNNING, PICKLING AND APPENDING TO NCF.

After model run, use make_figs_C3.py to reproduce Figs. 6-10 in the manuscript

"""
import sys
sys.path.append(r'\repositories\SpaFHy_v1') # path to folder where SpaFHy is

import os
import numpy as np
import pickle
#import matplotlib.pyplot as plt
#from netCDF4 import Dataset #, date2num, num2date
#from datetime import datetime

# import spafhy components
import spafhy
from spafhy_parameters import soil_properties, parameters
from spafhy_io import read_catchment_data, read_FMI_weather, read_SVE_runoff

eps = np.finfo(float).eps

""" set up SpaFHy for the catchment """

# load parameter dictionaries
pgen, pcpy, pbu, ptop = parameters()
psoil = soil_properties()

# read gis data and create necessary inputs for model initialization
gisdata = read_catchment_data(pgen['catchment_id'], fpath=pgen['gis_folder'],
                              plotgrids=False, plotdistr=False)

gisdata['LAI_grass'] = 0.5 * gisdata['LAI_decid']
gisdata['LAI_shrub'] = 0.1 * gisdata['LAI_conif']

# initialize SpaFHy
spa, ncf, ncf_file = spafhy.initialize(pgen, pcpy, pbu, ptop, psoil, gisdata,
                                       ncf=True,
                                       cpy_outputs=False, 
                                       bu_outputs=False,
                                       top_outputs=False)

## create netCDF output file
#dlat, dlon = np.shape(spa.GisData['cmask'])
#
#resultsfile = os.path.join(spa.pgen['results_folder'], spa.pgen['ncf_file'])
#spa.ncf_file = resultsfile
#
#ncf, ncf_file = spafhy.initialize_netCDF(ID=spa.id, fname=resultsfile, lat0=spa.GisData['lat0'], 
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

""" 
run SpaFHy spinup
"""
print('***  Running spinup ***')
for k in range(0, Nspin):
    #print('step: ' + str(k))
    forc= FORC[['doy', 'Rg', 'Par', 'T', 'Prec', 'VPD', 'CO2','U']].iloc[k]
    
    spa.run_timestep(forc, ncf=False)

# dump into pickle, clear local and reload: this is a way to save model object 
# and start later from same state.
#
#print('--- pickling spa, deleting local spa ---')
#pickle.dump(spa, open('spa_state.pk', 'wb'))
#del spa
#ncf.close()
print('--done ---')

"""" 
run spafhy and save results into ncf: note that data in ncf file now starts 
from index Nspin +1; the first ones are Nan's
"""

# load instance from pickle
# 
#print('--- unpickling spa ---')
#spa = pickle.load(open('spa_state.pk', 'rb'))

## open ncf file to append new results
#ncf = Dataset(ncf_file, 'a')

print('***  Running SpaFHy ***')
for k in range(Nspin+1, Nsteps):
    print('step: ' + str(k))
    forc= FORC[['doy', 'Rg', 'Par', 'T', 'Prec', 'VPD', 'CO2','U']].iloc[k]
    
    spa.run_timestep(forc, ncf=ncf)

# close output file
ncf.close()

# pickle spa object for making figures later
with open(os.path.join(spa.pgen['results_folder'], 'C3model.pk'), 'wb') as f:
    run = (spa, Qmeas, FORC)
    pickle.dump(run, f)

print('***  Done: results in file ' + pgen['ncf_file'] + ' ***')
