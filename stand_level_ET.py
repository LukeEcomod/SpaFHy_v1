# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 14:45:51 2018

@author: slauniai
"""
import numpy as np
import os
# import pickle
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
#import timeit
from scipy import stats

from spafhy_point import SpaFHy_point
from spafhy_parameters import parameters
# from canopygrid import CanopyGrid
# from bucketgrid import BucketGrid

from spafhy_io import read_FMI_weather, read_HydeDaily
  
eps = np.finfo(float).eps  # machine epsilon
spathy_path = os.path.join('c:', r'c:\datat\spathydata')
#ff = os.path.join(r'c:\repositories\spathy\ini', 'spafhy_point_default0.ini')

# dump into pickle
#ou = os.path.join(r'c:\ModelResults\Spathy', 'Site_results_newest.pkl')
#pickle.dump(res, open(ou, 'wb'))
    
def run_forest_sites():
    # runs ECsite_eval for all sites and saves results into pickle
    
    sites=['FIHy', 'FICage4', 'FICage12', 'FISod', 'FIKal', 'FILet', 'SEKno', 'SESky2', 'SENor','FISii', 'FISiiA']
    sites = [sites[1]]
    #print sites
    #sites = ['FIHy']
    res = dict.fromkeys(sites)
    for s in sites:
        out, dat, forc = ECsite_eval(s)
        a = {}
        a['Wliq_mod'] = np.ravel(out[0].bu.results['Wliq'])
        a['Wliq_low'] = np.ravel(out[1].bu.results['Wliq'])
        a['Wliq_high']  = np.ravel(out[2].bu.results['Wliq'])
        a['ET_mod']  = np.ravel(out[0].cpy.results['ET']) 
        a['ET_low']  = np.ravel(out[1].cpy.results['ET'])
        a['ET_high']  = np.ravel(out[2].cpy.results['ET'])
        
        a['Tr_mod']  = np.ravel(out[0].cpy.results['Transpi']) 
        a['Tr_low']  = np.ravel(out[1].cpy.results['Transpi'])
        a['Tr_high']  = np.ravel(out[2].cpy.results['Transpi'])         

        a['Ef_mod'] = np.ravel(out[0].cpy.results['Efloor']) 
        # a['Ef_low'] = np.ravel(out[1].cpy.results['Efloor'])
        # a['Ef_high'] = np.ravel(out[2].cpy.results['Efloor']) 
        
        a['SWE_mod'] =  np.ravel(out[0].cpy.results['SWE'])
        a['SWE_low'] =  np.ravel(out[1].cpy.results['SWE'])
        a['SWE_high'] =  np.ravel(out[2].cpy.results['SWE'])

        a['data'] = dat
        a['forc'] = forc
        
        res[s] = a
        del a, out, forc
        
    return res
        
def ECsite_eval(site):

    p = {'FIHy': {'LAIc': 3.5, 'LAId': 0.5, 'hc': 15.0, 'soil': [0.44, 0.33, 0.13, 2e-6], 'orgd': 0.05, 'fb': True},
         'FICage4': {'LAIc': 0.6, 'LAId': 0.1, 'hc': 0.4, 'soil': [0.44, 0.33, 0.13, 2e-6], 'orgd': 0.05, 'fb': True},
         'FICage12': {'LAIc': 1.4, 'LAId': 0.4, 'hc': 1.7, 'soil': [0.44, 0.30, 0.13, 2e-6], 'orgd': 0.05, 'fb': True},
         'FISod': {'LAIc': 2.1, 'LAId': 0.1, 'hc': 15.0, 'soil': [0.41, 0.21, 0.05, 1e-4], 'orgd': 0.05, 'fb': True},
         'FIKal': {'LAIc': 2.1, 'LAId': 0.1, 'hc': 15.0, 'soil': [0.9, 0.42, 0.11, 5e-10], 'orgd': 0.08, 'fb': True},  # peat
         'FILet': {'LAIc': 4.3, 'LAId': 2.3, 'hc': 15.0, 'soil': [0.9, 0.42, 0.11, 5e-10], 'orgd': 0.08, 'fb': True},  # peat
         'SEKno': {'LAIc': 3.6, 'LAId': 0.2, 'hc': 15.0, 'soil': [0.44, 0.33, 0.13, 2e-6], 'orgd': 0.05, 'fb': True},  # medium textured
         'SESky2': {'LAIc': 5.3, 'LAId': 0.5, 'hc': 15.0, 'soil': [0.44, 0.33, 0.13, 1e-6], 'orgd': 0.05, 'fb': True}, 
         'SENor': {'LAIc': 5.5, 'LAId': 1.3, 'hc': 15.0, 'soil': [0.43, 0.33, 0.02, 1e-6], 'orgd': 0.05, 'fb': True}, 
         'FISii': {'LAIc': 0.01, 'LAId': 0.3, 'hc': 0.3, 'soil': [0.9, 0.42, 0.11, 5e-10], 'orgd': 0.20, 'fb': False},
         'FISiiA': {'LAIc': 0.01, 'LAId': 0.3, 'hc': 0.3, 'soil': [0.9, 0.42, 0.11, 5e-10], 'orgd': 0.20, 'fb': True}, 
         }
    
    pgen, pcpy, pbu, _ = parameters()  # default parameters
    
    pgen['spatial_cpy'] = False
    pgen['spatial_soil'] = False

    pcpy['lai_conif']= p[site]['LAIc']
    pcpy['lai_decid']= p[site]['LAId']
    pcpy['hc']= p[site]['hc']
    L = pcpy['lai_conif'] + pcpy['lai_decid']
    pcpy['cf'] = np.maximum(0.2, 1.5 * L / (1.5 * L + 1.43) - 0.2)

    pbu['poros'] = p[site]['soil'][0]
    pbu['fc'] = p[site]['soil'][1]
    pbu['wp'] = p[site]['soil'][2]
    pbu['ksat'] = p[site]['soil'][3]
    pbu['org_depth'] = p[site]['orgd']
    
    fb = p[site]['fb']
    
    """ read forcing data and evaluation data """
    if site in ['FIHy', 'FICage4', 'FICage12']:
        dat, FORC = read_hydedata(site)
    elif site in ['FISii', 'FISiiA']:
        dat, FORC = read_siikaneva_data()
    else:
        dat, FORC = read_daily_ECdata(site)
    
    FORC['Prec'] = FORC['Prec'] / pgen['dt']  # mms-1
    FORC['T'] = FORC['Ta'].copy()
    cmask = np.ones(1)
    # print(cmask)
    # run model for different parameter combinations, save results into dataframe
    # +/- 20%, 20%, 20%
    p = [[10.0, 2.1, 1.3, 4.5],
         [8.5, 1.7, 1.05, 3.15],
         [11.5, 2.5, 2.6, 5.4]]

    out = []
    for k in range(3):
        a = p[k]
        pcpy['amax'] = a[0]
        pcpy['g1_conif'] = a[1]
        pcpy['g1_decid'] = 1.6 * a[1]
        pcpy['wmax'] = a[2]
        pcpy['wmaxsnow'] = a[3]

        model = SpaFHy_point(pgen, pcpy, pbu, cmask, FORC, cpy_outputs=True, bu_outputs=True)
        nsteps=len(FORC)
        model._run(0, nsteps, soil_feedbacks=fb)
        
        out.append(model)
        del model
        
    # best model
    # Wliq_mod = np.ravel(out[0].bu.results['Wliq'])
    # Wliq_low = np.ravel(out[1].bu.results['Wliq'])
    # Wliq_high = np.ravel(out[2].bu.results['Wliq'])
    ET_mod = np.ravel(out[0].cpy.results['ET']) 
    ET_low = np.ravel(out[1].cpy.results['ET'])
    ET_high = np.ravel(out[2].cpy.results['ET']) 
    
    # E_mod = np.ravel(out[0].cpy.results['Evap']) 
    # E_low = np.ravel(out[1].cpy.results['Evap'])
    # E_high = np.ravel(out[2].cpy.results['Evap'])

    
    tvec = dat.index
    et_dry = dat['ET'].copy()
    et_dry[dat['Prec'] > 0.1] = np.NaN
    
    sns.set_style('whitegrid')
    with sns.color_palette('muted'):
        plt.figure()
        
        plt.subplot(2,3,(1,2))
                
        plt.plot(tvec, et_dry, 'o', markersize=4, alpha=0.3, label='meas')
        plt.fill_between(tvec, ET_low, ET_high, facecolor='grey', alpha=0.6, label='range')
        plt.plot(tvec, ET_mod, 'k-', alpha=0.4, lw=0.5, label='mod')
        #plt.xlim([pd.datetime(2003, 10, 1), pd.datetime(2011,1,1)])
        plt.legend(loc=2, fontsize=8)
        plt.setp(plt.gca().get_xticklabels(), fontsize=8)
        plt.setp(plt.gca().get_yticklabels(), fontsize=8)
        plt.ylabel('ET$_{dry}$ (mm d$^{-1}$)', fontsize=8)
        plt.ylim([-0.05, 5.0])
        #plt.xlim([pd.datetime(2002, 1, 1), pd.datetime(2011,1,1)])
        plt.title(site)
        # sns.despine()        
        
        # scatterplot
        plt.subplot(2,3,3)
        et_dry[dat['doy'] < 0] = np.NaN
        et_dry[dat['doy'] > 366] = np.NaN
        meas = np.array(et_dry.values.tolist())
        ix = np.where(np.isfinite(meas))
        meas = meas[ix].copy()
        mod = ET_mod[ix].copy()
        
        slope, intercept, r_value, p_value, std_err = stats.linregress(meas, mod)
        print('slope', slope, 'interc', intercept)
        # slope, intercept, _, _ = stats.theilslopes(meas, mod, 0.95)
        meas = meas[:, np.newaxis]
        slope, _, _, _ = np.linalg.lstsq(meas, mod)
        intercept = 0.0
        xx = np.array([min(meas), max(meas)])
        plt.plot(meas, mod, 'o', markersize=4, alpha=0.3)

        plt.plot(xx, slope*xx + intercept, 'k-')
        plt.plot([0, 5], [0, 5], 'k--')
        plt.text(0.3, 4.2, 'y = %.2f x + %.2f' %(slope, intercept), fontsize=8)

        plt.xlim([-0.01, 5]); plt.ylim([-0.01, 5])
        plt.setp(plt.gca().get_xticklabels(), fontsize=8)
        ax = plt.gca()
        ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()
        plt.ylabel('ET$_{dry}$ mod (mm d$^{-1}$)', fontsize=8)
        plt.setp(plt.gca().get_yticklabels(), fontsize=8)
        plt.xlabel('ET$_{dry}$ meas (mm d$^{-1}$)', fontsize=8)
        
        os.chdir(os.path.join(r'c:\ModelResults\Spathy'))
        tt = site + '_ET.png'
        plt.savefig(tt)
#        plt.subplot(2,3,(4,5))
#        #plt.plot(tvec, SWCa, 'o', markersize=4, alpha=0.3,label='meas')
#        plt.fill_between(tvec, Wliq_low, Wliq_high, facecolor='grey', alpha=0.6, label='range')
#        plt.plot(tvec, Wliq_mod, 'k-',alpha=0.4, lw=0.5, label='mod')        
#        #plt.xlim([pd.datetime(2003, 10, 1), pd.datetime(2011,1,1)])
#        plt.legend(loc=2, fontsize=8)
#        plt.setp(plt.gca().get_xticklabels(), fontsize=8)
#        plt.setp(plt.gca().get_yticklabels(), fontsize=8)
#        plt.ylabel('$\\theta$ (m$^3$ m$^{-3}$)', fontsize=8)
#        # plt.xlim([pd.datetime(2002, 1, 1), pd.datetime(2011,1,1)])


    return out, dat, FORC


def Hyde_eval(pgen, pcpy, pbu):
    # this deos some figs for hyde
    """ read forcing data and evaluation data """
    fname = r'c:\datat\spathydata\HydeDaily2000-2010.txt'
    dat, FORC = read_HydeDaily(fname)
    FORC['Prec'] = FORC['Prec'] / pgen['dt']  # mms-1
    FORC['T'] = FORC['Ta'].copy()
    cmask = np.ones(1)
    
    # run model for different parameter combinations, save results into dataframe
    # +/- 15%, 15%, 20%, 30%
    # amax, g1_conif, wmax, wmaxsnow
    p = [[10.0, 2.1, 1.3, 4.5],
         [8.5, 1.7, 1.05, 3.15],
         [11.5, 2.5, 2.6, 5.4]]

    out = []
    for k in range(3):
        a = p[k]
        pcpy['amax'] = a[0]
        pcpy['g1_conif'] = a[1]
        pcpy['g1_decid'] = 1.6*a[1]
        pcpy['wmax'] = a[2]
        pcpy['wmaxsnow'] = a[3]
        
        model = SpaFHy_point(pgen, pcpy, pbu, cmask, FORC, cpy_outputs=True, bu_outputs=True)
        nsteps=len(FORC)
        model._run(0, nsteps)
        
        out.append(model)
        del model
        
    # best model
    Wliq_mod = np.ravel(out[0].bu.results['Wliq'])
    Wliq_low = np.ravel(out[1].bu.results['Wliq'])
    Wliq_high = np.ravel(out[2].bu.results['Wliq'])
    ET_mod = np.ravel(out[0].cpy.results['ET']) 
    ET_low = np.ravel(out[1].cpy.results['ET'])
    ET_high = np.ravel(out[2].cpy.results['ET']) 
    
    E_mod = np.ravel(out[0].cpy.results['Evap']) 
    E_low = np.ravel(out[1].cpy.results['Evap'])
    E_high = np.ravel(out[2].cpy.results['Evap'])

#    SWC_mod = np.ravel(out[0].cpy.results['ET']) 
#    SWC_low = np.ravel(out[1].cpy.results['ET'])
#    SWC_high = np.ravel(out[2].cpy.results['ET'])) 
    
    SWCa = dat['SWCa']
    SWCb = dat['SWCb']
    SWCc = dat['SWCc']
    tvec = dat.index
    et_dry = dat['ET']
    et_dry[dat['Prec']>0.1] = np.NaN
    
    sns.set_style('whitegrid')
    with sns.color_palette('muted'):
        plt.figure()
        
        plt.subplot(2,3,(1,2))
                
        plt.plot(tvec, et_dry, 'o', markersize=4, alpha=0.3, label='meas')
        plt.fill_between(tvec, ET_low, ET_high, facecolor='grey', alpha=0.6, label='range')
        plt.plot(tvec, ET_mod, 'k-', alpha=0.4, lw=0.5, label='mod')
        #plt.xlim([pd.datetime(2003, 10, 1), pd.datetime(2011,1,1)])
        plt.legend(loc=2, fontsize=8)
        plt.setp(plt.gca().get_xticklabels(), fontsize=8)
        plt.setp(plt.gca().get_yticklabels(), fontsize=8)
        plt.ylabel('ET$_{dry}$ (mm d$^{-1}$)', fontsize=8)
        plt.ylim([-0.05, 5.0])
        plt.xlim([pd.datetime(2002, 1, 1), pd.datetime(2011,1,1)])

        # sns.despine()        
        plt.subplot(2,3,(4,5))

        plt.plot(tvec, SWCa, 'o', markersize=4, alpha=0.3,label='meas')
        plt.fill_between(tvec, Wliq_low, Wliq_high, facecolor='grey', alpha=0.6, label='range')
        plt.plot(tvec, Wliq_mod, 'k-',alpha=0.4, lw=0.5, label='mod')        
        #plt.xlim([pd.datetime(2003, 10, 1), pd.datetime(2011,1,1)])
        plt.legend(loc=2, fontsize=8)
        plt.setp(plt.gca().get_xticklabels(), fontsize=8)
        plt.setp(plt.gca().get_yticklabels(), fontsize=8)
        plt.ylabel('$\\theta$ (m$^3$ m$^{-3}$)', fontsize=8)
        plt.xlim([pd.datetime(2002, 1, 1), pd.datetime(2011,1,1)])

        # scatterplot
        plt.subplot(2,3,6)

        meas = np.array(SWCa.values.tolist())
        slope, intercept, r_value, p_value, std_err = stats.linregress(meas, Wliq_mod)
        #print slope, intercept
        xx = np.array([min(meas), max(meas)])
        plt.plot(meas, Wliq_mod, 'o', markersize=5, alpha=0.3)
        plt.plot(xx, slope*xx + intercept, 'k-')
        plt.plot([0.05, 0.45], [0.05, 0.45], 'k--')
        plt.text( 0.07, 0.42, 'y = %.2f x + %.2f' %(slope, intercept), fontsize=8)
        plt.xlim([0.05, 0.45]); plt.ylim([0.05, 0.45])
        
        ax = plt.gca()
        ax.set_yticks([0.1, 0.2, 0.3, 0.4])
        ax.set_xticks([0.1, 0.2, 0.3, 0.4])
        ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()
        plt.ylabel('$\\theta$ mod (m$^3$ m$^{-3}$)', fontsize=8)
        
        plt.setp(plt.gca().get_yticklabels(), fontsize=8)
        plt.setp(plt.gca().get_yticklabels(), fontsize=8)
        plt.setp(plt.gca().get_xticklabels(), fontsize=8)
        plt.xlabel('$\\theta$ meas (m$^3$ m$^{-3}$)', fontsize=8)
        
        # scatterplot
        plt.subplot(2,3,3)

        meas = np.array(et_dry.values.tolist())
        ix = np.where(np.isfinite(meas))
        meas=meas[ix].copy()
        mod = ET_mod[ix].copy()
        slope, intercept, r_value, p_value, std_err = stats.linregress(meas, mod)
        xx = np.array([min(meas), max(meas)])
        plt.plot(meas, mod, 'o', markersize=4, alpha=0.3)

        plt.plot(xx, slope*xx + intercept, 'k-')
        plt.plot([0, 5], [0, 5], 'k--')
        plt.text(0.3, 4.2, 'y = %.2f x + %.2f' %(slope, intercept), fontsize=8)

        plt.xlim([-0.01, 5]); plt.ylim([-0.01, 5])
        plt.setp(plt.gca().get_xticklabels(), fontsize=8)
        ax = plt.gca()
        ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()
        plt.ylabel('ET$_{dry}$ mod (mm d$^{-1}$)', fontsize=8)
        plt.setp(plt.gca().get_yticklabels(), fontsize=8)
        plt.xlabel('ET$_{dry}$ meas (mm d$^{-1}$)', fontsize=8)
        
        #plt.savefig('Hyde_validate.pdf')
        #plt.savefig('Hyde_validate.png')
        
        # snowpack and throughfall
        
        plt.figure()
        SWE_mod = np.ravel(out[0].cpy.results['SWE'])        
        SWE_low = np.ravel(out[1].cpy.results['SWE'])
        SWE_hi = np.ravel(out[2].cpy.results['SWE'])
        swe_meas = dat['SWE']

        plt.plot(tvec, swe_meas, 'o', markersize=10, alpha=0.3, label='meas')
        plt.fill_between(tvec, SWE_low, SWE_hi, facecolor='grey', alpha=0.6, label='range')
        plt.plot(tvec, SWE_mod, 'k-', alpha=0.4, lw=0.5, label='mod')
        plt.title('SWE'); plt.ylabel('SWE mm')

    return out, dat, FORC


def read_daily_ECdata(site):
    
    #if site=='FICage4':
        
    if site=='FISod':
        fpath = os.path.join(spathy_path, 'DailyCEIP', 'FMI_Sodankyla')
        yrs = np.arange(2001, 2010)
        fnames = [ 'SodDaily_%4d.dat' %(k) for k in yrs]
        cols=['doy','NEE','GPP','TER','ET','H','NEEflag','ETflag','Hflag','Par','Rnet','Ta','VPD','CO2','Prec','U','Pamb',
              'SWC1','SWC2','Tsoil1', 'Tsoil2', 'Rnetflag', 'Snowdepth']        
    
    if site=='FIKal':
        fpath = os.path.join(spathy_path, 'DailyCEIP', 'FMI_Kalevansuo')
        yrs = [2005, 2006, 2007, 2008]
        fnames = ['KalevansuoDaily_%4d.dat' %(k) for k in yrs]
        cols=['doy','NEE','GPP','TER','ET','H','NEEflag','ETflag','Hflag','Par','Rnet','Ta','VPD','CO2','empty1','Prec','U','Pamb',
              'WTD', 'Snowdepth', 'Rnetflag']        

    if site=='SEKno':
        fpath = os.path.join(spathy_path, 'DailyCEIP', 'Knottasen')
        yrs =[2007]
        fnames = ['KnoDaily_%4d.dat' %(k) for k in yrs]
        cols=['doy','NEE','GPP','TER','ET','H','NEEflag','ETflag','Hflag','Par','Rnet','Ta','VPD','CO2','Prec','U','Pamb',
              'SWC1','SWC2','Tsoil1', 'Tsoil2', 'Rnetflag', 'Snowdepth']   
              
    if site=='SENor':
        fpath = os.path.join(spathy_path, 'DailyCEIP', 'Norunda')
        yrs = [1996, 1997, 1999, 2003]
        fnames = ['NorDaily_%4d.dat' %(k) for k in yrs]
        cols=['doy','NEE','GPP','TER','ET','H','NEEflag','ETflag','Hflag','Par','Rnet','Ta','VPD','CO2','Prec','U','Pamb',
              'SWC1','SWC2','Tsoil1', 'Tsoil2', 'Rnetflag', 'Snowdepth']  

    if site=='SESky2':
        fpath = os.path.join(spathy_path, 'DailyCEIP', 'Skyttorp2')
        yrs = [2005]
        fnames = ['Sky2Daily_2005.dat']
        cols=['doy','NEE','GPP','TER','ET','H','NEEflag','ETflag','Hflag','Par','Rnet','Ta','VPD','CO2','Prec','U','Pamb',
              'SWC1','SWC2','Tsoil1', 'Tsoil2', 'Rnetflag', 'Snowdepth']  
    
    if site=='FILet':
        fpath = os.path.join(spathy_path, 'DailyCEIP', 'FMI_Lettosuo')
        yrs = [2010, 2011, 2012]
        fnames = ['FILet_Daily_%4d.dat' %(k) for k in yrs]
        cols=['year', 'month', 'day', 'doy','NEE','GPP','TER','ET','H','Rnet','Rg', 'Par','Prec_site', 'Prec', 'Ta', 'RH',
              'VPD','CO2','U','Pamb', 'WTD','WTDwest','WTDsouth', 'WTDeast', 'WTDnorth', 'SWC1', 'SWC2', 'empty', 'Ts1', 'Ts2',
              'Ts3', 'Ts4', 'NEEflag', 'ETflag', 'Hflag','Rnetflag']  
    
    if site == 'FICage4':
         cols = ['time','doy','NEE','GPP','TER','ET','H','NEEflag','ETflag','Hflag','Par','Rnet','Ta','VPD','CO2','SWCa','PrecSmear','Prec','U','Pamb']

    dat = pd.DataFrame()
    
    for k in range(len(yrs)):
        fname = os.path.join(fpath, fnames[k])
        tmp = pd.read_csv(fname,sep='\s+',header=None, names=cols)
        tmp['year'] = yrs[k]
        
        dat = dat.append(tmp)
    
    dat['doy'] = dat['doy'].astype(int)
    tvec = pd.to_datetime(dat['year'] * 1000 + dat['doy'], format='%Y%j')
    #dat.drop('year')
    dat.index = tvec
    
    # forcing data
  
    forc = dat[['doy', 'Ta', 'VPD', 'Prec', 'Par', 'U']]
    forc['Par'] = 1./4.6*forc['Par']
    forc['Rg'] = 2.0*forc['Par']
    forc['VPD'][forc['VPD'] <= 0] = eps
    forc = forc.interpolate()  # fills missing values by linear interpolation  
    forc['CO2'] = 380.0
    #relatively extractable water, from soil moisture
    forc['Rew'] = 1.0
    if site=='SEKno':
        forc['Rg'] = 1.4*forc['Rg']
#    fc=0.30
#    wp=0.10
#    Wliq=dat['SWCa']
#    Rew=np.maximum( 0.0, np.minimum( (Wliq-wp)/(fc - wp + eps), 1.0) )
#    forc['Rew']=Rew     
    
    return dat, forc


def read_hydedata(site):
    
    fname = os.path.join(spathy_path, 'HydeDaily2000-2010.txt')
    cols=['time','doy','NEE','GPP','TER','ET','H','NEEflag','ETflag','Hflag','Par','Rnet','Ta','VPD','CO2','PrecSmear','Prec','U','Pamb',
    'SWE0','SWCh','SWCa','SWCb','SWCc', 'Tsh','Tsa','Tsb','Tsc','RnetFlag','Trfall','Snowdepth','Snowdepthstd','SWE','SWEstd','Roff1','Roff2']        
    
    dat=pd.read_csv(fname,sep='\s+',header=None, names=None, parse_dates=[[0,1,2]], keep_date_col=False)
    dat.columns=cols
    dat.index=dat['time']; dat=dat.drop(['time','SWE0'],axis=1)
    
    forc=dat[['doy','Ta','VPD','Prec','Par','U']]; forc['Par']= 1/4.6*forc['Par']; forc['Rg']=2.0*forc['Par']
    forc['VPD'][forc['VPD']<=0]=eps
    
    #relatively extractable water, Hyde A-horizon
    #poros = 0.45    
    fc = 0.30
    wp = 0.10
    Wliq = dat['SWCa']
    Rew = np.maximum( 0.0, np.minimum((Wliq-wp)/(fc - wp + eps), 1.0) )
    forc['Rew'] = Rew
    forc['CO2'] = 380.0
    forc = forc.interpolate()  # fills missing values by linear interpolation

    if site == 'FICage4':
        fname = os.path.join(spathy_path, 'HydeCage4yr-2000.txt')
        cols = ['time','doy','NEE','GPP','TER','ET','H','NEEflag','ETflag','Hflag','Par','Rnet',
              'Ta','VPD','CO2','SWCa','PrecSmear','Prec','U','Pamb']        
    
        dat = pd.read_csv(fname,sep='\s+',header=None, names=None, parse_dates=[[0,1,2]], keep_date_col=False)
        dat.columns=cols
        dat.index = dat['time']
        dat = dat.drop('time',axis=1)
        
        forc = forc.ix[dat.index]
    if site == 'FICage12':
        fname = os.path.join(spathy_path, 'HydeCage12yr-2002.txt')
        cols = ['time','doy','NEE','GPP','TER','ET','H','NEEflag','ETflag','Hflag','Par','Rnet',
              'Ta','VPD','CO2','SWCa','PrecSmear','Prec','U','Pamb']        
    
        dat = pd.read_csv(fname,sep='\s+',header=None, names=None, parse_dates=[[0,1,2]], keep_date_col=False)
        dat.columns=cols
        dat.index = dat['time']
        dat = dat.drop('time',axis=1)
        
        forc = forc.ix[dat.index]

    return dat, forc
    
def read_siikaneva_data():
    # data available from years 2011, 2013 (good), 2014
    yrs = [2011, 2013, 2014]
    fpath = os.path.join(spathy_path, 'DailyCEIP', 'Siikaneva1')
    fname = os.path.join(fpath, 'Sii_%s_daily.dat' % (yrs[0]))
    dat = pd.read_csv(fname, sep=';', header='infer')
    
    start = pd.datetime(yrs[0] ,1, 1)
    end =pd.datetime(yrs[0], 12, 31)
    tm = pd.date_range(start, end, freq='1d')
    dat.index = tm
    
    for k in yrs[1:]:
        fname = os.path.join(fpath, 'Sii_%s_daily.dat' % (k))
        tmp = pd.read_csv(fname,sep=';',header='infer')
        start = pd.datetime(k ,1, 1)
        end = pd.datetime(k, 12, 31)
        tm = pd.date_range(start, end, freq='1d')
        tmp.index = tm
        
        dat = dat.append(tmp)
    

    
#    dat = pd.read_csv(fname, sep=';', header='infer')
#    start = pd.datetime(year ,1, 1)
#    end =pd.datetime(year, 12, 31)
#    tm = pd.date_range(start, end, freq='1d')
#    dat.index = tm
    
    forc = dat[['doy', 'Ta', 'VPD', 'Prec', 'Rg','U']]
    forc['VPD'][forc['VPD'] <= 0] = eps
    forc['Par'] = 0.5 * forc['Rg']
    forc['CO2'] = 380.0
    forc['Rew'] = 1.0
    dat['ET'][dat['q_ET'] > 0.5] = np.NaN
    
    return dat, forc
    # doy;NEE;GPP;TER;ET;H;q_NEE;q_ET;Rglob;Rn;Ta;VPD;Prec;U;Pamb;WTD