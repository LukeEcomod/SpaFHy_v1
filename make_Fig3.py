# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 12:26:14 2018

@author: slauniai

Makes Fig 3 of Launiainen et al. 2019 GMD
"""

import os
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns
from scipy import stats

eps = np.finfo(float).eps

#%% plot FIHy figure (Fig. 3)

def draw_Fig3(data, mod):

    # modeled daily values
    ET = mod[0]['ET']
    ET_low = mod[1]['ET']
    ET_high = mod[2]['ET']
    
    Wliq = mod[0]['Wliq']
    Wliq_low = mod[1]['Wliq']
    Wliq_high = mod[2]['Wliq']
    
    SWE = mod[0]['SWE']
    SWE_low = mod[1]['SWE']
    SWE_high = mod[2]['SWE']
    
    # measured values
    et_dry = data['ET']
    Prec = data['Prec']
    Ts = data['Tsh']
    et_dry[Prec > 0.1] = np.NaN   
    tvec = data.index
    doy = data['doy']
    
    # soil moisture
    SWCa = data['SWCa']
    SWCa[Ts <= 0.5] = np.NaN
    SWCa[SWCa > 0.5] = np.NaN

    # SWE
    SWEm = data['SWE']
    SWEm = SWEm.dropna()    
          
    sns.set_style('whitegrid')
    with sns.color_palette('muted'):
        fig = plt.figure()
        
        fig.set_size_inches(6.5, 7.5)
        
        plt.subplot(3,3,(1,2))
                
        plt.plot(tvec, et_dry, 'o', markersize=4, alpha=0.3, label='obs')
        plt.fill_between(tvec, ET_low, ET_high, facecolor='grey', alpha=0.5, label='range')
        plt.plot(tvec, ET, 'k-', alpha=0.6, lw=0.5, label='mod')
        #plt.xlim([pd.datetime(2003, 10, 1), pd.datetime(2011,1,1)])
        #plt.legend(loc=2, fontsize=8)
        plt.setp(plt.gca().get_xticklabels(), fontsize=8)
        plt.setp(plt.gca().get_yticklabels(), fontsize=8)
        plt.ylabel('ET (mm d$^{-1}$)', fontsize=9)
        plt.ylim([-0.05, 5.0])
        plt.xlim([pd.datetime(2002, 1, 1), pd.datetime(2011,1,1)])
        tstr = 'FIHy: %.1f (%.1f) m$^2$m$^{-2}$' % (4.0, 0.5)
        plt.title(tstr, fontsize=9)   
        
        # scatterplot
        plt.subplot(3,3,3)
        xxx = et_dry.copy()
        xxx[doy < 120] = np.NaN
        xxx[doy > 273] = np.NaN
        meas = np.array(xxx.values.tolist())
        ix = np.where(np.isfinite(meas))[0]
        meas = meas[ix].copy()
        mod = ET[ix].copy()
        
        slope, intercept, r_value, p_value, std_err = stats.linregress(meas, mod)
        print('slope', slope, 'interc', intercept)
        # slope, intercept, _, _ = stats.theilslopes(meas, mod, 0.95)
        # force regression through origin
        x = meas[:, np.newaxis]
        slope, res, _, _  = np.linalg.lstsq(x, mod)
        r2 = (1.0 - res / sum((mod - np.mean(mod))**2))[0]
        intercept = 0.0
        rmse = np.sqrt(((mod - meas) ** 2).mean())
        me = np.mean(mod - meas)
        
        xx = np.array([min(meas), max(meas)])
        plt.plot(meas, mod, 'o', markersize=4, alpha=0.3)
    
        plt.plot(xx, slope*xx, 'k-')
        plt.plot([0, 5], [0, 5], 'k--', linewidth=1)
        #plt.text(0.3, 4.5, 'y = %.2fx, R$^2$=%.2f' %(slope, r2), fontsize=8)
        tst = 's=%.2f\nR$^2$=%.2f\nME=%.2f' %(slope, r2, me)
        plt.text(0.3, 3.6, tst, fontsize=8)
        plt.xlim([-0.01, 5]); plt.ylim([-0.01, 5])
        ax = plt.gca()
        ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()
        ax.set_xticks(np.arange(0, 6))
        ax.set_yticks(np.arange(0, 6))
        ax.set_xticklabels(np.arange(0, 6), fontsize=8)
        ax.set_yticklabels(np.arange(0, 6), fontsize=8)
        ax.set_aspect('equal')
        
        plt.ylabel('ET$_{mod}$ (mm d$^{-1}$)', fontsize=9)
        plt.setp(plt.gca().get_yticklabels(), fontsize=8)
        plt.xlabel('ET$_{obs}$ (mm d$^{-1}$)', fontsize=9, labelpad=-3)
        
        # soil moisture
        
        plt.subplot(3,3,(4,5))
        plt.plot(tvec, SWCa, 'o', markersize=4, alpha=0.3,label='obs')
        plt.fill_between(tvec, Wliq_low, Wliq_high, facecolor='grey', alpha=0.6, label='range')
        plt.plot(tvec, Wliq, 'k-',alpha=0.6, lw=0.5, label='mod')        
        #plt.legend(loc=2, fontsize=8)
        #plt.xlim([pd.datetime(2003, 10, 1), pd.datetime(2011,1,1)])
        #plt.legend(loc=2, fontsize=8)
        plt.ylabel('$\\theta$ (m$^3$ m$^{-3}$)', fontsize=9)
        plt.setp(plt.gca().get_xticklabels(), fontsize=8)
        plt.setp(plt.gca().get_yticklabels(), fontsize=8)
        #plt.ylabel('$\\theta$ (m$^3$ m$^{-3}$)', fontsize=8)
        plt.xlim([pd.datetime(2002, 1, 1), pd.datetime(2011,1,1)])
    
    
        # scatterplot
        plt.subplot(3,3,6)
        
        meas = np.array(SWCa.values.tolist())
        ix = np.where(np.isfinite(meas))[0]
        meas = meas[ix].copy()
        mod = Wliq[ix].copy()
        slope, intercept, r_value, p_value, std_err = stats.linregress(meas, mod)
        r2 = r_value**2
        me = np.mean(mod - meas)
        
        #print slope, intercept
        xx = np.array([min(meas), max(meas)])
        plt.plot(meas, mod, 'o', markersize=4, alpha=0.3)
        plt.plot(xx, slope*xx + intercept, 'k-')
        plt.plot([0.05, 0.45], [0.05, 0.45], 'k--', linewidth=1)
        #plt.text( 0.15, 0.08, 'y = %.2f x + %.2f' %(slope, intercept), fontsize=8)
        tst = 's=%.2f\nR$^2$=%.2f\nME=%.2f' %(slope, r2, me)
        plt.text(0.28, 0.08, tst, fontsize=8)
        plt.xlim([0.05, 0.45]); plt.ylim([0.05, 0.45])
        
        ax = plt.gca()
        ax.set_yticks([0.1, 0.2, 0.3, 0.4])
        ax.set_xticks([0.1, 0.2, 0.3, 0.4])
        ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()    
        ax.set_aspect('equal') 
        plt.setp(plt.gca().get_yticklabels(), fontsize=8)
        plt.setp(plt.gca().get_xticklabels(), fontsize=8)
        plt.ylabel('$\\theta_{mod}$ (m$^3$ m$^{-3}$)', fontsize=9)
        plt.xlabel('$\\theta_{obs}$ (m$^3$ m$^{-3}$)', fontsize=9, labelpad=-3)
            
        
        #%plot SWE
        
        plt.subplot(3,3,(7,8))
        plt.plot(SWEm, 'o', markersize=4, alpha=0.3,label='obs')
        plt.fill_between(tvec, SWE_low, SWE_high, facecolor='grey', alpha=0.6, label='range')
        plt.plot(tvec, SWE, 'k-', alpha=0.6, lw=0.7, label='mod')        
        plt.legend(loc=2, fontsize=8)
        plt.ylim([-0.1, 150])
        plt.xlim([pd.datetime(2002, 1, 1), pd.datetime(2011,1,1)])
        plt.ylabel('SWE (mm)', fontsize=9)
        plt.setp(plt.gca().get_yticklabels(), fontsize=8)
        plt.setp(plt.gca().get_xticklabels(), fontsize=8)
        
        x = SWEm.values
        y = SWE
        y = y.loc[SWEm.index].values
        y = y.ravel()
        
        plt.subplot(3,3,9)
        
        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
        r2 = r_value**2
        me = np.mean(y - x)
        
        #print slope, intercept
        xx = np.array([min(x), max(x)])
        plt.plot(x, y, 'o', markersize=4, alpha=0.3)
        plt.plot(xx, slope*xx + intercept, 'k-')
        plt.plot([-0.1, 150], [-0.1, 150], 'k--', linewidth=1)
        #plt.text( 0.15, 0.08, 'y = %.2f x + %.2f' %(slope, intercept), fontsize=8)
        tst = 's=%.2f\nR$^2$=%.2f\nME=%.2f' %(slope, r2, me)
        plt.text(10, 95, tst, fontsize=8)
        plt.xlim([-0.1, 150]); plt.ylim([-0.1, 150])
        ax = plt.gca()
        ax.set_yticks([0, 25, 50, 75, 100, 125])
        ax.set_xticks([0, 25, 50, 75, 100, 125])
        ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()    
        ax.set_aspect('equal') 
        plt.setp(plt.gca().get_yticklabels(), fontsize=8)
        plt.setp(plt.gca().get_xticklabels(), fontsize=8)
        plt.ylabel('SWE$_{mod}$ (mm)', fontsize=9)
        plt.xlabel('SWE$_{obs}$ (mm)', fontsize=9)
            

        tt = 'Fig3_FIHy.png'
        plt.savefig(tt, dpi=600)
        plt.savefig('Fig3_FIHy.pdf')
        
#%%
""" read hyde data, make simulations and plot figure """
from spafhy_point import SpaFHy_point
from spafhy_io import read_HydeDaily
from spafhy_parameters import parameters_FIHy

pgen, pcpy, pbu = parameters_FIHy()

# read forcing data to dataframe
dat, FORC = read_HydeDaily(pgen['forcing_file'])
FORC['Prec'] = FORC['Prec'] / pgen['dt']  # mms-1
FORC['T'] = FORC['Ta'].copy()

# np.array cmask is needed to apply model components at a single point
cmask = np.ones(1)
    
# run model for different parameter combinations, save results into dataframe
# +/- 15%, 15%, 20%, 30%
# amax, g1_conif, wmax, wmaxsnow
p = [[10.0, 2.1, 3.5, 1.3, 4.5],
     [8.5, 1.7, 2.8, 1.05, 3.15],
     [11.5, 2.5, 4.2, 2.6, 5.4]]

out = []
results = []
for k in range(3):
    a = p[k]
    pcpy['amax'] = a[0]
    pcpy['g1_conif'] = a[1]
    pcpy['g1_decid'] = a[2]
    pcpy['wmax'] = a[3]
    pcpy['wmaxsnow'] = a[4]

    model = SpaFHy_point(pgen, pcpy, pbu, FORC, cmask=cmask, cpy_outputs=True, bu_outputs=True)
    nsteps=len(FORC)
    model._run(0, nsteps)
    
    # during model run, results are stored in two dictionaries: 
    # canopycrid outputs: model.cpy.results
    # bucketgid outputs: model.bu.resutls
    
    # extract them, convert to dataframe and save to csv
    cres = model.cpy.results
    bres = model.bu.results
    del bres['ET'] # also bucket returns ET but in [m] so remove key
    res = {**cres, **bres} # combine into one dict
    res = {key: np.ravel(res[key]) for key in res.keys()} # and ravel
    
    # convert to dataframe and save
    resu = pd.DataFrame(data=res, columns=res.keys(), index=model.FORC.index, dtype=float)
    out.append(model)
    results.append(resu)
    
    del model, resu
