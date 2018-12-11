# -*- coding: utf-8 -*-
"""
Created on Tue Oct 02 14:47:14 2018

@author: slauniai

DEMONSTRATES SPAFHY SETUP, RUNNING, PICKLING AND APPENDING TO NCF
"""
#import os
import numpy as np
# import pandas as pd
import pickle
import matplotlib.pyplot as plt
from netCDF4 import Dataset #, date2num, num2date
#from datetime import datetime

import spafhy
from spafhy_parameters_default import soil_properties, parameters
#from spafhy_io import read_FMI_weather, write_AsciiGrid
from spafhy_io import create_catchment, read_FMI_weather

eps = np.finfo(float).eps

""" set up SpaFHy for the catchment """

# load parameter dictionaries
pgen, pcpy, pbu, ptop = parameters()
psoil = soil_properties()

# read gis data and create necessary inputs for model initialization
gisdata = create_catchment(pgen['catchment_id'], fpath=pgen['gis_folder'],
                           plotgrids=False, plotdistr=False)
# initialize spafhy
spa = spafhy.initialize(pgen, pcpy, pbu, ptop, psoil, gisdata, cpy_outputs=False, 
                 bu_outputs=False, top_outputs=False, flatten=True)

# create netCDF output file
dlat, dlon = np.shape(spa.GisData['cmask'])

ncf, ncf_file = spafhy.initialize_netCDF(ID=spa.id, fname=spa.ncf_file, lat0=spa.GisData['lat0'], 
                                         lon0=spa.GisData['lon0'], dlat=dlat, dlon=dlon, dtime=None)

# read forcing data
""" read forcing data and catchment runoff file """
FORC = read_FMI_weather(pgen['catchment_id'],
                        pgen['start_date'],
                        pgen['end_date'],
                        sourcefile=pgen['forcing_file'])
FORC['Prec'] = FORC['Prec'] / spa.dt  # mms-1
FORC['U'] = 2.0 # use constant wind speed ms-1
Nsteps = len(FORC)

#%% now test running spafhy in steps, input data to ncf_file and picking/unpickling

# run spafhy for 1st year
for k in range(0, 365):
    print('step: ' + str(k))
    forc= FORC[['doy', 'Rg', 'Par', 'T', 'Prec', 'VPD', 'CO2','U']].iloc[k]
    
    spa.run_timestep(forc, ncf=ncf)

# close output file
ncf.close()

# plot soil water content
plt.figure()
# function spa._to_grid(x) converts x from 1D array back to 2d
plt.imshow(spa._to_grid(spa.bu.Wliq)); plt.colorbar(); plt.title('wliq, step ' +str(spa.step_nr))

#%% dump into pickle, clear local and reload

print('--- pickling spa, deleting local spa ---')
ff = open('spa_state.pk', 'wb')
pickle.dump(spa, ff)
ff.close()
del spa
print('--done ---')

#%% unpicle spa
# HERE ONE COULD CREATE SEVERAL INSTANCES FROM SAME INITIAL STATE - E.G. BY MODIFYING
# CANOPYGRID STATE VARIABLES
print('--- unpickling spa ---')
ff = open('spa_state.pk', 'rb')
spa = pickle.load(ff)
ff.close()
print('--done ---')

# open ncf file for append new results
ncf = Dataset(ncf_file, 'a')

#%%%
# run from current state + 180 days
for k in range(spa.step_nr, spa.step_nr + 365):
    print('step: ' + str(k))
    forc= FORC[['doy', 'Rg', 'Par', 'T', 'Prec', 'VPD', 'CO2','U']].iloc[k]
    
    spa.run_timestep(forc, ncf=ncf)
    
ncf.close()
# plot soil water content
plt.figure()
plt.imshow(spa._to_grid(spa.bu.Wliq)); plt.colorbar(); plt.title('wliq, step ' +str(spa.step_nr))

print('--- done: ncf_file closed. Open for reading as ncf = Dataset(ncf_file, "r")')
