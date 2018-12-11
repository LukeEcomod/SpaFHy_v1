# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 10:13:52 2018

@author: slauniai

DEMO HOW TO RUN POINT-SCALE MODEL FOR A SINGLE OR MULTIPLE SITES.


"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# import model and functions to read data
from spafhy_point import SpaFHy_point
from spafhy_io import read_HydeDaily # , read_FMI_weather 

eps = np.finfo(float).eps  # machine epsilon

""" 
Select demo version: 
1 - run for one site, 
2 = several runs with varying parameters. reproduce Fig. 3
3 = as demo 1 but for 3 different soil types
4 = as demo 3 but in single call
"""
demo = 1

#%%
""" demo 1: """

if demo == 1:
    # set up model for single site; here use FIHy as example
    from spafhy_parameters import parameters_FIHy
    
    # read parameters
    pgen, pcpy, pbu = parameters_FIHy()
    
    # read forcing data to dataframe
    dat, FORC = read_HydeDaily(pgen['forcing_file'])
    FORC['Prec'] = FORC['Prec'] / pgen['dt']  # mms-1
    FORC['T'] = FORC['Ta'].copy()
    
    # np.array cmask is needed to apply model components at a single point
    cmask = np.ones(1)
    
    # create model instance
    model = SpaFHy_point(pgen, pcpy, pbu, FORC, cmask=cmask, 
                         cpy_outputs=True, bu_outputs=True)
    
    
    # simple model run with default parameters:
    Nsteps = len(model.FORC)  # length of forcing
    model._run(0, Nsteps) # calls class function _run to execute model from t=0 to t=Nsteps
    
    # during model run, resutls are stored in two dictionaries: 
    # canopycrid outputs: model.cpy.results
    # bucketgid outputs: model.bu.resutls
    
    # extract them, convert to dataframe and save to csv
    cres = model.cpy.results
    bres = model.bu.results
    del bres['ET'] # also bucket returns ET but in [m] so remove key
    res = {**cres, **bres} # combine into one dict
    res = {key: np.ravel(res[key]) for key in res.keys()} # and ravel
    
    # convert to dataframe and save
    results = pd.DataFrame(data=res, columns=res.keys(), index=model.FORC.index, dtype=float)
    results.to_csv(pgen['output_file'] + '.csv', sep=';')
    del res, cres, bres
    
    #%% now, let's draw timeseries of root zone and organic layer water content and ET
    # components to Fig. 2
    
    plt.figure()
    
    plt.subplot(211)
    plt.plot(results[['Wliq', 'Wliq_top']])
    plt.legend([r'$\theta$', r'$\theta_{org}$'])
    plt.ylabel(r'$\theta$ (m$^{3}$ m$^{-3}$)')
    
    plt.subplot(212)
    plt.plot(results[['Evap', 'Transpi', 'Efloor']])
    plt.legend(['E', 'Tr', 'Ef'])
    plt.ylabel('mm d$^{-1}$')
        

#%%
    
""" demo 2: override parameters and run 3 simulations """
if demo == 2:

    # set up model for single site; here use FIHy as example
    from spafhy_parameters import parameters_FIHy
    
    # read parameters
    pgen, pcpy, pbu = parameters_FIHy()
    
    # read forcing data to dataframe
    dat, FORC = read_HydeDaily(pgen['forcing_file'])
    FORC['Prec'] = FORC['Prec'] / pgen['dt']  # mms-1
    FORC['T'] = FORC['Ta'].copy()
    
    # np.array cmask is needed to apply model components at a single point
    cmask = np.ones(1)
    
    # run model for 3 different parameter combinations: vary
    # g1_conif, g1_decid, wmax, wmaxsnow by  +/- 20%, 20%, 20%, 30%
    
    p =[
        [1.0, 1.0, 1.0, 1.0], # nominal
        [0.8, 0.8, 0.8, 0.7], # low-ET case
        [1.2, 1.2, 1.2, 1.3] # high-ET case
        ]

    # save results to list 
    out = []
    for k in range(3):
        a = p[k]
        # read nominal parameters and modify some of them
        pgen, pcpy, pbu = parameters_FIHy()
    
        pcpy['physpara']['g1_conif'] *= a[0]
        pcpy['physpara']['g1_decid'] *= a[1]
        pcpy['interc']['wmax'] *= a[2]
        pcpy['interc']['wmaxsnow'] *= a[3]
        
        # create model instance and run:
        model = SpaFHy_point(pgen, pcpy, pbu, FORC, cmask=cmask, cpy_outputs=True, bu_outputs=True)
        nsteps=len(FORC)
        model._run(0, nsteps)

        # extract results, convert to dataframe, print to file and append to out
        cres = model.cpy.results
        bres = model.bu.results
        del bres['ET']
        res = {**cres, **bres} # combine into one dict
        res = {key: np.ravel(res[key]) for key in res.keys()} # and convert each variable into 1D array

        results = pd.DataFrame(data=res, columns=res.keys(), index=model.FORC.index)        
        results.to_csv(pgen['output_file'] + '_sim_' + str(k) + '.csv', sep=';')
        
        out.append(results)
        
        del model, res, cres, bres, results, pcpy, pgen, pbu

    #  plot Fig 3 equivalent
    from make_Fig3 import draw_Fig3
    draw_Fig3(dat, out)

#%%

""" demo 3: as demo2 but run for 3 different soil types defined in soil_properties"""

if demo == 3:    
    # soil classes
    soilclass = np.array([1, 2, 3]) # coarse, medium, fine
    cmask = np.ones(1)
    # set up model for single site; here use FIHy as example
    from spafhy_parameters import parameters_FIHy, soil_properties
    from spafhy_io import preprocess_soildata
    
    # read forcing data to dataframe
    dat, FORC = read_HydeDaily(pgen['forcing_file'])
    FORC['Prec'] = FORC['Prec'] / pgen['dt']  # mms-1
    FORC['T'] = FORC['Ta'].copy()
    
    out = []
    # run ofe 3 soil classes
    for k in range(3):

        # read parameters and soil properties
        pgen, pcpy, pbu = parameters_FIHy()
        psoil = soil_properties()

        # get soil properties based on soilclass and update pbu

        pbu = preprocess_soildata(pbu, psoil, soilclass[k], cmask=cmask, spatial=True)
        print(pbu)
        
        # create model instance and run:
        model = SpaFHy_point(pgen, pcpy, pbu, FORC, cmask=cmask, cpy_outputs=True, bu_outputs=True)
        nsteps=len(FORC)
        model._run(0, nsteps)

        # extract results, convert to dataframe, print to file and append to out
        cres = model.cpy.results
        bres = model.bu.results
        del bres['ET']
        res = {**cres, **bres} # combine into one dict
        res = {key: np.ravel(res[key]) for key in res.keys()} # and convert each variable into 1D array

        results = pd.DataFrame(data=res, columns=res.keys(), index=model.FORC.index)        
        results.to_csv(pgen['output_file'] + '_sim_' + str(k) + '.csv', sep=';')
        
        out.append(results)
        
        del model, res, cres, bres, results, pcpy, pgen, pbu
        
    # plot figure of soil water content and Transpiration at each soil class
    plt.figure()
    
    ax1 = plt.subplot(211)
    ax1.plot(out[2]['Wliq'], label='fine text')
    ax1.plot(out[1]['Wliq'], label='medium text')
    ax1.plot(out[0]['Wliq'], label='coarse text')
    ax1.legend()
    ax1.set_ylabel(r'$\theta$ (m$^{3}$ m$^{-3}$)')
    
    ax2 = plt.subplot(212, sharex=ax1)
    ax2.plot(out[2]['Transpi'], label='fine text')
    ax2.plot(out[1]['Transpi'], label='medium text')
    ax2.plot(out[0]['Transpi'], label='coarse text')
    ax2.legend()
    ax2.set_ylabel('mm d$^{-1}$')
        
#%%

""" demo 4: as demo2 but run all 3 different soil types at once"""

if demo == 4:    
    # soil classes
    soilclass = np.array([1, 2, 3]) # coarse, medium, fine
    cmask = np.ones(3)
    # set up model for single site; here use FIHy as example
    from spafhy_parameters import parameters_FIHy, soil_properties
    from spafhy_io import preprocess_soildata

    # read parameters and soil properties
    pgen, pcpy, pbu = parameters_FIHy()
    psoil = soil_properties()
    
    # read forcing data to dataframe
    dat, FORC = read_HydeDaily(pgen['forcing_file'])
    FORC['Prec'] = FORC['Prec'] / pgen['dt']  # mms-1
    FORC['T'] = FORC['Ta'].copy()
    
    # get soil properties based on soilclass and update pbu

    pbu = preprocess_soildata(pbu, psoil, soilclass, cmask=cmask, spatial=True)
    #print(pbu)
    
    # create model instance and run:
    model = SpaFHy_point(pgen, pcpy, pbu, FORC, cmask=cmask, cpy_outputs=True, bu_outputs=True)
    nsteps=len(FORC)
    model._run(0, nsteps)

    # extract results, convert to dataframe, print to file and append to out
    cres = model.cpy.results
    bres = model.bu.results
    del bres['ET']
    res = {**cres, **bres} # combine into one dict
    res = {key: np.array(res[key]) for key in res.keys()}
    
    # now res is dict where each key contains np.array which shape is (nsteps,3)
    # to save each column into separate csv-file and plot figures, we do follwing:
    out = []
    n = 0
    for s in soilclass:
        dummy = {key: res[key][:,n] for key in res.keys()}
        results = pd.DataFrame(data=dummy, columns=res.keys(), index=model.FORC.index)        
        results.to_csv(pgen['output_file'] + '_sim_' + str(s) + '.csv', sep=';')
        n += 1
    
        out.append(results)
    
    del model, res, cres, bres, results, pcpy, pgen, pbu
        
    # plot figure of soil water content and Transpiration at each soil class
    plt.figure()
    
    ax1 = plt.subplot(211)
    ax1.plot(out[2]['Wliq'], label='fine text')
    ax1.plot(out[1]['Wliq'], label='medium text')
    ax1.plot(out[0]['Wliq'], label='coarse text')
    ax1.legend()
    ax1.set_ylabel(r'$\theta$ (m$^{3}$ m$^{-3}$)')
    
    ax2 = plt.subplot(212, sharex=ax1)
    ax2.plot(out[2]['Transpi'], label='fine text')
    ax2.plot(out[1]['Transpi'], label='medium text')
    ax2.plot(out[0]['Transpi'], label='coarse text')
    ax2.legend()
    ax2.set_ylabel('mm d$^{-1}$')
