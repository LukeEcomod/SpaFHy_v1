# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 12:20:39 2017

@author: L1566
"""

import numpy as np
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import xarray as xr
import matplotlib.lines as mlines
import matplotlib.gridspec as gridspec
from datetime import datetime, timedelta

def plot_streamflow(ID,scen,q_fmi,q_clim,q_clim1,dailymean=True,timeseries=False,fmi_clim=False):
    ### daily average ###
    D_q_fmi = q_fmi.groupby('dtime.dayofyear').mean()
    D_q_clim = q_clim.groupby('dtime.dayofyear').mean()
    D_q_clim1 = q_clim1.groupby('dtime.dayofyear').mean()
    # mean daily value within each month of the year
    Mmean_q_fmi = q_fmi.groupby('dtime.month').mean()
    Mmean_q_clim = q_clim.groupby('dtime.month').mean()
    Mmean_q_clim1 = q_clim1.groupby('dtime.month').mean()
    ### ploting streamflow over 30 years ###
    if dailymean:
        sns.set_style('ticks')
        plt.figure(figsize = (12,8))
        #gs = gridspec.GridSpec(2, 2,width_ratios=[2, 1])
        plt.subplot(211)
        D_q_fmi.plot.line('bo',markersize = 4, alpha = 0.3,label='fmi(1980-2009)')
        D_q_clim.plot.line('r',linewidth = 0.5,label='climate(1980-2009)')
        D_q_clim1.plot.line('k',linewidth = 0.5,label='climate(2070-2099)')
        plt.xlim([0,365])
        plt.xlabel('day of year')
        plt.ylabel('streamflow (mm/d)')
        plt.legend(loc='best')
        ###
        plt.subplot(212)
        Mmean_q_fmi.plot.line('b--',linewidth = 1,label='fmi(1980-2009)')
        Mmean_q_clim.plot.line('r',linewidth = 1,label='climate(1980-2009)')
        Mmean_q_clim1.plot.line('k',linewidth = 1,label='climate(2070-2099)')
        plt.xlabel('month of year')
        plt.ylabel('streamflow (mm/d)')
        plt.legend(loc='best')
        
        sns.despine()
        plt.tight_layout()
        plt.show()
        plt.savefig('C:\Apps\WinPython-64bit-2.7.10.3\SAMULI_ARI\spathy\\results\\climate\\c_'+ID+'_'+scen+'_daily_monthly_mean_flow.pdf',dpi=300)
        
    if fmi_clim:
        sns.set_style('ticks')
        plt.figure(figsize = (12,4))
        q_fmi.plot.line('bo',markersize = 4, alpha = 0.2,label='fmi(1980-2009)')
        q_clim.plot.line('r',linewidth = 0.5,alpha = 0.5,label='climate(1980-2009)')
        plt.xlabel('date')
        plt.ylabel('streamflow')
        plt.legend(loc='best')
        sns.despine()
        plt.show()
        plt.savefig('C:\Apps\WinPython-64bit-2.7.10.3\SAMULI_ARI\spathy\\results\\climate\\fmi_clim_flow_c_'+ID+'_'+scen+'.pdf',dpi=300)

    if timeseries:
        sns.set_style('ticks')
        plt.figure(figsize = (12,4))
        plt.plot(q_clim,'bo',markersize = 4, alpha = 0.2,label='climate(1980-2009)')
        plt.plot(q_clim1,'r',linewidth = 0.5,alpha = 0.5,label='climate(2070-2099)')
        plt.xlabel('number of days')
        plt.ylabel('streamflow')
        plt.legend(loc='best')
        sns.despine()
        plt.show()
        plt.savefig('C:\Apps\WinPython-64bit-2.7.10.3\SAMULI_ARI\spathy\\results\\climate\\c_'+ID+'_'+scen+'_clim_clim1_flow.pdf',dpi=300)
        
def plot_sat_deficit(ID,scen,S_fmi,S_clim,S_clim1): 
    #daily mean
    D_S_fmi = S_fmi.groupby('dtime.dayofyear').mean()
    D_S_clim = S_clim.groupby('dtime.dayofyear').mean()
    D_S_clim1 = S_clim1.groupby('dtime.dayofyear').mean()
    # mean daily value within each month of the year
    Mmean_S_fmi = S_fmi.groupby('dtime.month').mean()
    Mmean_S_clim = S_clim.groupby('dtime.month').mean()
    Mmean_S_clim1 = S_clim1.groupby('dtime.month').mean()
    ### ploting daily sat. deficit over 30 years ###
    sns.set_style('ticks')
    plt.figure(figsize = (12,8))
    plt.subplot(211)
    D_S_fmi.plot.line('bo',markersize = 4, alpha = 0.3,label='fmi(1980-2009)')
    D_S_clim.plot.line('r',linewidth = 0.5,label='climate(1980-2009)')
    D_S_clim1.plot.line('k',linewidth = 0.5,label='climate(2070-2099)')
    plt.xlim([0,365])
    plt.xlabel('day of year')
    plt.ylabel('sat. deficit (m)')
    plt.legend(loc='best')
    ###
    plt.subplot(212)
    Mmean_S_fmi.plot.line('b--',linewidth = 1,label='fmi(1980-2009)')
    Mmean_S_clim.plot.line('r',linewidth = 1,label='climate(1980-2009)')
    Mmean_S_clim1.plot.line('k',linewidth = 1,label='climate(2070-2099)')
    plt.xlabel('month of year')
    plt.ylabel('sat. deficit (m)')
    plt.legend(loc='best')
    sns.despine()
    plt.tight_layout()
    plt.show()
    plt.savefig('C:\Apps\WinPython-64bit-2.7.10.3\SAMULI_ARI\spathy\\results\\climate\\c_'+ID+'_'+scen+'_daily_monthly_mean_sdeficit.pdf',dpi=300)

def mapping_changes(ID,scen,d_fmi,d_clim,d_clim1):
    ## mean ET components over 30 years
    TmET_fmi = (d_fmi.Transpi+d_fmi.Efloor+d_fmi.Evap).mean(dim=['dtime'])
    TmTrans_fmi = d_fmi.Transpi.mean(dim=['dtime'])
    TmEfloor_fmi = d_fmi.Efloor.mean(dim=['dtime'])
    TmEvap_fmi = d_fmi.Evap.mean(dim=['dtime'])
    TmSWE_fmi = d_fmi.SWE.mean(dim=['dtime'])
    ## mean ET components over 30 years
    TmET_clim = (d_clim.Transpi+d_clim.Efloor+d_clim.Evap).mean(dim=['dtime'])
    TmTrans_clim = d_clim.Transpi.mean(dim=['dtime'])
    TmEfloor_clim = d_clim.Efloor.mean(dim=['dtime'])
    TmEvap_clim = d_clim.Evap.mean(dim=['dtime'])
    TmSWE_clim = d_clim.SWE.mean(dim=['dtime'])
    ## mean ET components over 30 years
    TmET_clim1 = (d_clim1.Transpi+d_clim1.Efloor+d_clim1.Evap).mean(dim=['dtime'])
    TmTrans_clim1 = d_clim1.Transpi.mean(dim=['dtime'])
    TmEfloor_clim1 = d_clim1.Efloor.mean(dim=['dtime'])
    TmEvap_clim1 = d_clim1.Evap.mean(dim=['dtime'])
    TmSWE_clim1 = d_clim1.SWE.mean(dim=['dtime'])
    sns.set_style('white')
    plt.figure(figsize = (25,8))
    plt.subplot(251)
    TmET_clim.plot(cmap='jet',vmin=0.45,vmax=0.78)
    plt.xlabel('longitude')
    plt.ylabel('lattitude')
    plt.title('climate(1980-2009)')
    plt.subplot(256)
    TmET_clim1.plot(cmap='jet',vmin=0.45,vmax=0.78)
    plt.xlabel('longitude')
    plt.ylabel('lattitude')
    plt.title('climate(2070-2099)')
    plt.subplot(252)
    TmTrans_clim.plot(cmap='jet',vmin=0.0,vmax=0.72)
    plt.xlabel('longitude')
    plt.ylabel('lattitude')
    plt.title('climate(1980-2009)')
    plt.subplot(257)
    TmTrans_clim1.plot(cmap='jet',vmin=0.0,vmax=0.72)
    plt.xlabel('longitude')
    plt.ylabel('lattitude')
    plt.title('climate(2070-2099)')
    plt.subplot(253)
    TmEfloor_clim.plot(cmap='jet',vmin=0.0,vmax=0.66)
    plt.xlabel('longitude')
    plt.ylabel('lattitude')
    plt.title('climate(1980-2009)')
    plt.subplot(258)
    TmEfloor_clim1.plot(cmap='jet',vmin=0.0,vmax=0.66)
    plt.xlabel('longitude')
    plt.ylabel('lattitude')
    plt.title('climate(2070-2099)')
    plt.subplot(254)
    TmEvap_clim.plot(cmap='jet',vmin=0.0,vmax=0.43)
    plt.xlabel('longitude')
    plt.ylabel('lattitude')
    plt.title('climate(1980-2009)')
    plt.subplot(259)
    TmEvap_clim1.plot(cmap='jet',vmin=0.0,vmax=0.43)
    plt.xlabel('longitude')
    plt.ylabel('lattitude')
    plt.title('climate(2070-2099)')
    plt.subplot(255)
    TmSWE_clim.plot(cmap='jet',vmin=8.0,vmax=36)
    plt.xlabel('longitude')
    plt.ylabel('lattitude')
    plt.title('climate(1980-2009)')
    plt.subplot(2,5,10)
    TmSWE_clim1.plot(cmap='jet',vmin=8.0,vmax=36)
    plt.xlabel('longitude')
    plt.ylabel('lattitude')
    plt.title('climate(2070-2099)')
    plt.tight_layout()
    plt.show()
    
    plt.savefig('C:\Apps\WinPython-64bit-2.7.10.3\SAMULI_ARI\spathy\\results\\climate\\c_'+ID+'_'+scen+'_mapping_changes.pdf',dpi=300)
    
def plot_ET_components(ID,scen,mET_clim,mTrans_clim,mEfloor_clim,mEvap_clim,mSWE_clim,mET_clim1,mTrans_clim1,mEfloor_clim1,mEvap_clim1,mSWE_clim1):
    #1980-2009 daily mean 
    D_mET_clim = mET_clim.groupby('dtime.dayofyear').mean()
    D_mTrans_clim = mTrans_clim.groupby('dtime.dayofyear').mean()
    D_mEfloor_clim = mEfloor_clim.groupby('dtime.dayofyear').mean()
    D_mEvap_clim = mEvap_clim.groupby('dtime.dayofyear').mean()
    D_mSWE_clim = mSWE_clim.groupby('dtime.dayofyear').max()
    #2070-2099
    D_mET_clim1 = mET_clim1.groupby('dtime.dayofyear').mean()
    D_mTrans_clim1 = mTrans_clim1.groupby('dtime.dayofyear').mean()
    D_mEfloor_clim1 = mEfloor_clim1.groupby('dtime.dayofyear').mean()
    D_mEvap_clim1 = mEvap_clim1.groupby('dtime.dayofyear').mean()
    D_mSWE_clim1 = mSWE_clim1.groupby('dtime.dayofyear').max()
    # 1980-2009 mean daily value within each month of the year
    Mmean_mET_clim = mET_clim.groupby('dtime.month').mean()
    Mmean_mTrans_clim = mTrans_clim.groupby('dtime.month').mean()
    Mmean_mEfloor_clim = mEfloor_clim.groupby('dtime.month').mean()
    Mmean_mEvap_clim = mEvap_clim.groupby('dtime.month').mean()
    Mmean_mSWE_clim = mSWE_clim.groupby('dtime.month').mean()
    Mmax_mSWE_clim = mSWE_clim.groupby('dtime.month').max()
    #2070-2099
    Mmean_mET_clim1 = mET_clim1.groupby('dtime.month').mean()
    Mmean_mTrans_clim1 = mTrans_clim1.groupby('dtime.month').mean()
    Mmean_mEfloor_clim1 = mEfloor_clim1.groupby('dtime.month').mean()
    Mmean_mEvap_clim1 = mEvap_clim1.groupby('dtime.month').mean()
    Mmean_mSWE_clim1 = mSWE_clim1.groupby('dtime.month').mean()
    Mmax_mSWE_clim1 = mSWE_clim1.groupby('dtime.month').max()
    ### ploting daily and monthly cat. av. ET, Trans, Efloor, Evap over 30 years ###    
    sns.set_style('ticks')
    plt.figure(figsize = (18,15))
    gs = gridspec.GridSpec(5, 2,width_ratios=[2, 1])
    plt.subplot(gs[0])
    #D_mET_fmi.plot.line('bo',markersize = 4, alpha = 0.3,label='fmi(1980-2009)')
    D_mET_clim.plot.line('r',linewidth = 0.5,label='climate(1980-2009)')
    D_mET_clim1.plot.line('k',linewidth = 0.5,label='climate(2070-2099)')
    plt.xlim([0,365])
    plt.xlabel('day of year')
    plt.ylabel('ET (mm)')
    plt.legend(loc='best')
    ###
    plt.subplot(gs[1])
    #plt.plot(Mmean_mET_fmi.index,Mmean_mET_fmi,'b-',linewidth = 0.5,label='fmi(1980-2009)')
    Mmean_mET_clim.plot.line('r',linewidth = 0.5,label='climate(1980-2009)')
    Mmean_mET_clim1.plot.line('k',linewidth = 0.5,label='climate(2070-2099)')
    plt.xlabel('month')
    plt.ylabel('ET (mm)')
    plt.legend(loc='best')
    
    plt.subplot(gs[2])
    #D_mTrans_fmi.plot.line('bo',markersize = 4, alpha = 0.3,label='fmi(1980-2009)')
    D_mTrans_clim.plot.line('r',linewidth = 0.5,label='climate(1980-2009)')
    D_mTrans_clim1.plot.line('k',linewidth = 0.5,label='climate(2070-2099)')
    plt.xlim([0,365])
    plt.xlabel('day of year')
    plt.ylabel('Transpiration (mm)')
    plt.legend(loc='best')
    ###
    plt.subplot(gs[3])
    #plt.plot(Mmean_mTrans_fmi.index,Mmean_mTrans_fmi,'b-',linewidth = 0.5,label='fmi(1980-2009)')
    Mmean_mTrans_clim.plot.line('r',linewidth = 0.5,label='climate(1980-2009)')
    Mmean_mTrans_clim1.plot.line('k',linewidth = 0.5,label='climate(2070-2099)')
    plt.xlabel('month')
    plt.ylabel('Transpiration (mm)')
    plt.legend(loc='best')
    
    plt.subplot(gs[4])
    #D_mEfloor_fmi.plot.line('bo',markersize = 4, alpha = 0.3,label='fmi(1980-2009)')
    D_mEfloor_clim.plot.line('r',linewidth = 0.5,label='climate(1980-2009)')
    D_mEfloor_clim1.plot.line('k',linewidth = 0.5,label='climate(2070-2099)')
    plt.xlim([0,365])
    plt.xlabel('day of year')
    plt.ylabel('Efloor (mm)')
    plt.legend(loc='best')
    ###
    plt.subplot(gs[5])
    #plt.plot(Mmean_mEfloor_fmi.index,Mmean_mEfloor_fmi,'b-',linewidth = 0.5,label='fmi(1980-2009)')
    Mmean_mEfloor_clim.plot.line('r',linewidth = 0.5,label='climate(1980-2009)')
    Mmean_mEfloor_clim1.plot.line('k',linewidth = 0.5,label='climate(2070-2099)')
    plt.xlabel('month')
    plt.ylabel('Efloor (mm)')
    plt.legend(loc='best')
    
    plt.subplot(gs[6])
    #D_mEvap_fmi.plot.line('bo',markersize = 4, alpha = 0.3,label='fmi(1980-2009)')
    D_mEvap_clim.plot.line('r',linewidth = 0.5,label='climate(1980-2009)')
    D_mEvap_clim1.plot.line('k',linewidth = 0.5,label='climate(2070-2099)')
    plt.xlim([0,365])
    plt.xlabel('day of year')
    plt.ylabel('Evap (mm)')
    plt.legend(loc='best')
    ###
    plt.subplot(gs[7])
    #plt.plot(Mmean_mEvap_fmi.index,Mmean_mEvap_fmi,'b-',linewidth = 0.5,label='fmi(1980-2009)')
    Mmean_mEvap_clim.plot.line('r',linewidth = 0.5,label='climate(1980-2009)')
    Mmean_mEvap_clim1.plot.line('k',linewidth = 0.5,label='climate(2070-2099)')
    plt.xlabel('month')
    plt.ylabel('Evap (mm)')
    plt.legend(loc='best')
    
    plt.subplot(gs[8])
    #D_mSWE_fmi.plot.line('bo',markersize = 4, alpha = 0.3,label='fmi(1980-2009)')
    D_mSWE_clim.plot.line('r',linewidth = 0.5,label='climate(1980-2009)')
    D_mSWE_clim1.plot.line('k',linewidth = 0.5,label='climate(2070-2099)')
    plt.xlim([0,365])
    plt.xlabel('day of year')
    plt.ylabel('SWE (mm)')
    plt.legend(loc='best')
    ###
    plt.subplot(gs[9])
    #plt.plot(Mmean_mSWE_fmi.index,Mmean_mSWE_fmi,'b-',linewidth = 0.5,label='fmi(1980-2009)')
    Mmax_mSWE_clim.plot.line('r',linewidth = 0.5,label='climate(1980-2009)')
    Mmax_mSWE_clim1.plot.line('k',linewidth = 0.5,label='climate(2070-2099)')
    plt.xlabel('month')
    plt.ylabel('SWE (mm)')
    plt.legend(loc='best')
    plt.tight_layout()
    sns.despine()
    plt.show()
    plt.savefig('C:\Apps\WinPython-64bit-2.7.10.3\SAMULI_ARI\spathy\\results\\climate\\c_'+ID+'_'+scen+'_daily_monthly_mean_ETs.pdf',dpi=300)
    
def plt_ETs_fmi_comparison(ID,scen,mET_clim,mTrans_clim,mEfloor_clim,mEvap_clim,mSWE_clim,mET_fmi,mTrans_fmi,mEfloor_fmi,mEvap_fmi,mSWE_fmi):
    #1980-2009 daily mean 
    D_mET_clim = mET_clim.groupby('dtime.dayofyear').mean()
    D_mTrans_clim = mTrans_clim.groupby('dtime.dayofyear').mean()
    D_mEfloor_clim = mEfloor_clim.groupby('dtime.dayofyear').mean()
    D_mEvap_clim = mEvap_clim.groupby('dtime.dayofyear').mean()
    D_mSWE_clim = mSWE_clim.groupby('dtime.dayofyear').mean()
    #fmi data
    D_mET_fmi = mET_fmi.groupby('dtime.dayofyear').mean()
    D_mTrans_fmi = mTrans_fmi.groupby('dtime.dayofyear').mean()
    D_mEfloor_fmi = mEfloor_fmi.groupby('dtime.dayofyear').mean()
    D_mEvap_fmi = mEvap_fmi.groupby('dtime.dayofyear').mean()
    D_mSWE_fmi = mSWE_fmi.groupby('dtime.dayofyear').mean()
    # 1980-2009 mean daily value within each month of the year
    Mmean_mET_clim = mET_clim.groupby('dtime.month').mean()
    Mmean_mTrans_clim = mTrans_clim.groupby('dtime.month').mean()
    Mmean_mEfloor_clim = mEfloor_clim.groupby('dtime.month').mean()
    Mmean_mEvap_clim = mEvap_clim.groupby('dtime.month').mean()
    Mmean_mSWE_clim = mSWE_clim.groupby('dtime.month').mean()
    Mmax_mSWE_clim = mSWE_clim.groupby('dtime.month').max()
    #fmi data
    Mmean_mET_fmi = mET_fmi.groupby('dtime.month').mean()
    Mmean_mTrans_fmi = mTrans_fmi.groupby('dtime.month').mean()
    Mmean_mEfloor_fmi = mEfloor_fmi.groupby('dtime.month').mean()
    Mmean_mEvap_fmi = mEvap_fmi.groupby('dtime.month').mean()
    Mmean_mSWE_fmi = mSWE_fmi.groupby('dtime.month').mean()
    Mmax_mSWE_fmi = mSWE_fmi.groupby('dtime.month').max()
    ### ploting daily and monthly cat. av. ET, Trans, Efloor, Evap over 30 years ###   
    sns.set_style('ticks')
    plt.figure(figsize = (18,15))
    gs = gridspec.GridSpec(5, 2,width_ratios=[2, 1])
    plt.subplot(gs[0])
    D_mET_fmi.plot.line('bo',markersize = 4, alpha = 0.3,label='fmi(1980-2009)')
    D_mET_clim.plot.line('r',linewidth = 0.5,label='climate(1980-2009)')
    plt.xlim([0,365])
    plt.xlabel('day of year')
    plt.ylabel('ET (mm)')
    plt.legend(loc='best')
    ###
    plt.subplot(gs[1])
    Mmean_mET_fmi.plot.line('b-',linewidth = 0.5,label='fmi(1980-2009)')
    Mmean_mET_clim.plot.line('r',linewidth = 0.5,label='climate(1980-2009)')
    plt.xlabel('month')
    plt.ylabel('ET (mm)')
    plt.legend(loc='best')
    
    plt.subplot(gs[2])
    D_mTrans_fmi.plot.line('bo',markersize = 4, alpha = 0.3,label='fmi(1980-2009)')
    D_mTrans_clim.plot.line('r',linewidth = 0.5,label='climate(1980-2009)')
    plt.xlim([0,365])
    plt.xlabel('day of year')
    plt.ylabel('Transpiration (mm)')
    plt.legend(loc='best')
    ###
    plt.subplot(gs[3])
    Mmean_mTrans_fmi.plot.line('b-',linewidth = 0.5,label='fmi(1980-2009)')
    Mmean_mTrans_clim.plot.line('r',linewidth = 0.5,label='climate(1980-2009)')
    plt.xlabel('month')
    plt.ylabel('Transpiration (mm)')
    plt.legend(loc='best')
    
    plt.subplot(gs[4])
    D_mEfloor_fmi.plot.line('bo',markersize = 4, alpha = 0.3,label='fmi(1980-2009)')
    D_mEfloor_clim.plot.line('r',linewidth = 0.5,label='climate(1980-2009)')
    plt.xlim([0,365])
    plt.xlabel('day of year')
    plt.ylabel('Efloor (mm)')
    plt.legend(loc='best')
    ###
    plt.subplot(gs[5])
    Mmean_mEfloor_fmi.plot.line('b-',linewidth = 0.5,label='fmi(1980-2009)')
    Mmean_mEfloor_clim.plot.line('r',linewidth = 0.5,label='climate(1980-2009)')
    plt.xlabel('month')
    plt.ylabel('Efloor (mm)')
    plt.legend(loc='best')
    
    plt.subplot(gs[6])
    D_mEvap_fmi.plot.line('bo',markersize = 4, alpha = 0.3,label='fmi(1980-2009)')
    D_mEvap_clim.plot.line('r',linewidth = 0.5,label='climate(1980-2009)')
    plt.xlim([0,365])
    plt.xlabel('day of year')
    plt.ylabel('Evap (mm)')
    plt.legend(loc='best')
    ###
    plt.subplot(gs[7])
    Mmean_mEvap_fmi.plot.line('b-',linewidth = 0.5,label='fmi(1980-2009)')
    Mmean_mEvap_clim.plot.line('r',linewidth = 0.5,label='climate(1980-2009)')
    plt.xlabel('month')
    plt.ylabel('Evap (mm)')
    plt.legend(loc='best')
    
    plt.subplot(gs[8])
    D_mSWE_fmi.plot.line('bo',markersize = 4, alpha = 0.3,label='fmi(1980-2009)')
    D_mSWE_clim.plot.line('r',linewidth = 0.5,label='climate(1980-2009)')
    plt.xlim([0,365])
    plt.xlabel('day of year')
    plt.ylabel('SWE (mm)')
    plt.legend(loc='best')
    ###
    plt.subplot(gs[9])
    Mmean_mSWE_fmi.plot.line('b-',linewidth = 0.5,label='fmi(1980-2009)')
    Mmean_mSWE_clim.plot.line('r',linewidth = 0.5,label='climate(1980-2009)')
    plt.xlabel('month')
    plt.ylabel('SWE (mm)')
    plt.legend(loc='best')
    plt.tight_layout()
    sns.despine()
    plt.show()   
    plt.savefig('C:\Apps\WinPython-64bit-2.7.10.3\SAMULI_ARI\spathy\\results\\climate\\c_'+ID+'_'+scen+'_daily_monthly_fmi_clim_ETs.pdf',dpi=300)
def Ann_Change(ID,scen,q_clim,mET_clim,mTrans_clim,mEfloor_clim,mEvap_clim,q_clim1,
               mET_clim1,mTrans_clim1,mEfloor_clim1,mEvap_clim1,A_p_clim,A_p_clim1):
    #############################################
    ## variables from climate data simulations (1981-2010)
    #############################################
    A_q_clim = q_clim.groupby('dtime.year').sum()
    A_mET_clim = mET_clim.groupby('dtime.year').sum()
    A_mTrans_clim = mTrans_clim.groupby('dtime.year').sum()
    A_mEfloor_clim = mEfloor_clim.groupby('dtime.year').sum()
    A_mEvap_clim = mEvap_clim.groupby('dtime.year').sum()
    #############################################
    ## variables from climate data simulations (2071-2100)
    #############################################
    A_q_clim1 = q_clim1.groupby('dtime.year').sum()
    A_mET_clim1 = mET_clim1.groupby('dtime.year').sum()
    A_mTrans_clim1 = mTrans_clim1.groupby('dtime.year').sum()
    A_mEfloor_clim1 = mEfloor_clim1.groupby('dtime.year').sum()
    A_mEvap_clim1 = mEvap_clim1.groupby('dtime.year').sum()
    
    ################################################
    ###  data in a Dataframe for boxplot ###########
    ##################################################
    mAQ = np.zeros(len(A_q_clim))
    mAET = np.zeros(len(A_q_clim))
    mATr = np.zeros(len(A_q_clim))
    mAEf = np.zeros(len(A_q_clim))
    mAEv = np.zeros(len(A_q_clim))
    for j in range(len(A_q_clim)):
        mAQ[j]=100.*(A_q_clim1[j]/A_p_clim1[j]-A_q_clim[j]/A_p_clim[j])/(A_q_clim[j]/A_p_clim[j])
        mAET[j]=100.*(A_mET_clim1[j]/A_p_clim1[j]-A_mET_clim[j]/A_p_clim[j])/(A_mET_clim[j]/A_p_clim[j])
        mATr[j]=100.*(A_mTrans_clim1[j]/A_p_clim1[j]-A_mTrans_clim[j]/A_p_clim[j])/(A_mTrans_clim[j]/A_p_clim[j])
        mAEf[j]=100.*(A_mEfloor_clim1[j]/A_p_clim1[j]-A_mEfloor_clim[j]/A_p_clim[j])/(A_mEfloor_clim[j]/A_p_clim[j])
        mAEv[j]=100.*(A_mEvap_clim1[j]/A_p_clim1[j]-A_mEvap_clim[j]/A_p_clim[j])/(A_mEvap_clim[j]/A_p_clim[j])
       
    mAM= pd.DataFrame()
    mAM['Q_P_ANN'] = mAQ
    mAM['ET_P_ANN'] = mAET
    mAM['Tr_P_ANN'] = mATr
    mAM['Ef_P_ANN'] = mAEf
    mAM['Evap_P_ANN'] = mAEv 
    
    ##### plotting boxplot
    sns.set_style('ticks')
    plt.figure(figsize=(6,6))
    sns.boxplot(data=mAM,showmeans=True,showfliers=False,linewidth=0.5,meanprops=dict(markersize=2))
    plt.ylabel('relative change (%)')
    plt.xticks(rotation='vertical')
    plt.axhline(linewidth=0.5)
    sns.despine()
    plt.tight_layout()
    plt.savefig('C:\Apps\WinPython-64bit-2.7.10.3\SAMULI_ARI\spathy\\results\\climate\\c_'+ID+'_'+scen+'_ann_change.pdf',dpi=300)
    
def Seasonal_Change(ID,scen,d_clim,d_clim1):
    #############################################
    ## variables from climate data simulations (1981-2010)
    #############################################
    q_clim = 1e3*d_clim.Qt
    #S_clim = d_clim.S
    mET_clim = (d_clim.Transpi+d_clim.Efloor+d_clim.Evap).mean(dim=['dlat','dlon'])
    mTrans_clim = d_clim.Transpi.mean(dim=['dlat','dlon'])
    mEfloor_clim = d_clim.Efloor.mean(dim=['dlat','dlon'])
    mEvap_clim = d_clim.Evap.mean(dim=['dlat','dlon'])
    mSWE_clim = d_clim.SWE.mean(dim=['dlat','dlon'])
    #############################################
    ## variables from climate data simulations (2071-2100)
    #############################################
    q_clim1 = 1e3*d_clim1.Qt
    #S_clim1 = d_clim1.S
    mET_clim1 = (d_clim1.Transpi+d_clim1.Efloor+d_clim1.Evap).mean(dim=['dlat','dlon'])
    mTrans_clim1 = d_clim1.Transpi.mean(dim=['dlat','dlon'])
    mEfloor_clim1 = d_clim1.Efloor.mean(dim=['dlat','dlon'])
    mEvap_clim1 = d_clim1.Evap.mean(dim=['dlat','dlon'])
    mSWE_clim1 = d_clim1.SWE.mean(dim=['dlat','dlon'])
    ################################################
    ###  seasonal mean daily ###########
    ##################################################
    month_clim = q_clim.groupby('dtime.month').apply(lambda x: x).month
    DJF = (month_clim <= 2) | (month_clim >= 12)
    MAM = (month_clim >= 3) & (month_clim <= 5)
    JJA = (month_clim >= 6) & (month_clim <= 8)
    SON = (month_clim >= 9) & (month_clim <= 11)
    # seasonal daily mean flow
    q_DJF = q_clim.where(DJF).resample('AS-Dec', 'dtime', how='mean')
    q_MAM = q_clim.where(MAM).resample('AS-Mar', 'dtime', how='mean')
    q_JJA = q_clim.where(JJA).resample('AS-Jun', 'dtime', how='mean')
    q_SON = q_clim.where(SON).resample('AS-Sep', 'dtime', how='mean')
    # seasonal daily mean ET
    mET_DJF = mET_clim.where(DJF).resample('AS-Dec', 'dtime', how='mean')
    mET_MAM = mET_clim.where(MAM).resample('AS-Mar', 'dtime', how='mean')
    mET_JJA = mET_clim.where(JJA).resample('AS-Jun', 'dtime', how='mean')
    mET_SON = mET_clim.where(SON).resample('AS-Sep', 'dtime', how='mean')
    # seasonal daily mean Transpi
    mTrans_DJF = mTrans_clim.where(DJF).resample('AS-Dec', 'dtime', how='mean')
    mTrans_MAM = mTrans_clim.where(MAM).resample('AS-Mar', 'dtime', how='mean')
    mTrans_JJA = mTrans_clim.where(JJA).resample('AS-Jun', 'dtime', how='mean')
    mTrans_SON = mTrans_clim.where(SON).resample('AS-Sep', 'dtime', how='mean')
    # seasonal daily mean Efloor
    mEfloor_DJF = mEfloor_clim.where(DJF).resample('AS-Dec', 'dtime', how='mean')
    mEfloor_MAM = mEfloor_clim.where(MAM).resample('AS-Mar', 'dtime', how='mean')
    mEfloor_JJA = mEfloor_clim.where(JJA).resample('AS-Jun', 'dtime', how='mean')
    mEfloor_SON = mEfloor_clim.where(SON).resample('AS-Sep', 'dtime', how='mean')
    # seasonal daily mean Evap
    mEvap_DJF = mEvap_clim.where(DJF).resample('AS-Dec', 'dtime', how='mean')
    mEvap_MAM = mEvap_clim.where(MAM).resample('AS-Mar', 'dtime', how='mean')
    mEvap_JJA = mEvap_clim.where(JJA).resample('AS-Jun', 'dtime', how='mean')
    mEvap_SON = mEvap_clim.where(SON).resample('AS-Sep', 'dtime', how='mean')
    # seasonal daily mean SWE
    mSWE_DJF = mSWE_clim.where(DJF).resample('AS-Dec', 'dtime', how='mean')
    mSWE_MAM = mSWE_clim.where(MAM).resample('AS-Mar', 'dtime', how='mean')
    mSWE_JJA = mSWE_clim.where(JJA).resample('AS-Jun', 'dtime', how='mean')
    mSWE_SON = mSWE_clim.where(SON).resample('AS-Sep', 'dtime', how='mean')
    #####################################    
    # seasonal daily mean flow for 2070-2099
    #####################################
    month_clim = q_clim1.groupby('dtime.month').apply(lambda x: x).month
    DJF = (month_clim <= 2) | (month_clim >= 12)
    MAM = (month_clim >= 3) & (month_clim <= 5)
    JJA = (month_clim >= 6) & (month_clim <= 8)
    SON = (month_clim >= 9) & (month_clim <= 11)
    q1_DJF = q_clim1.where(DJF).resample('AS-Dec', 'dtime', how='mean')
    q1_MAM = q_clim1.where(MAM).resample('AS-Mar', 'dtime', how='mean')
    q1_JJA = q_clim1.where(JJA).resample('AS-Jun', 'dtime', how='mean')
    q1_SON = q_clim1.where(SON).resample('AS-Sep', 'dtime', how='mean')
    # seasonal daily mean ET
    mET1_DJF = mET_clim1.where(DJF).resample('AS-Dec', 'dtime', how='mean')
    mET1_MAM = mET_clim1.where(MAM).resample('AS-Mar', 'dtime', how='mean')
    mET1_JJA = mET_clim1.where(JJA).resample('AS-Jun', 'dtime', how='mean')
    mET1_SON = mET_clim1.where(SON).resample('AS-Sep', 'dtime', how='mean')
    # seasonal daily mean Transpi
    mTrans1_DJF = mTrans_clim1.where(DJF).resample('AS-Dec', 'dtime', how='mean')
    mTrans1_MAM = mTrans_clim1.where(MAM).resample('AS-Mar', 'dtime', how='mean')
    mTrans1_JJA = mTrans_clim1.where(JJA).resample('AS-Jun', 'dtime', how='mean')
    mTrans1_SON = mTrans_clim1.where(SON).resample('AS-Sep', 'dtime', how='mean')
    # seasonal daily mean Efloor
    mEfloor1_DJF = mEfloor_clim1.where(DJF).resample('AS-Dec', 'dtime', how='mean')
    mEfloor1_MAM = mEfloor_clim1.where(MAM).resample('AS-Mar', 'dtime', how='mean')
    mEfloor1_JJA = mEfloor_clim1.where(JJA).resample('AS-Jun', 'dtime', how='mean')
    mEfloor1_SON = mEfloor_clim1.where(SON).resample('AS-Sep', 'dtime', how='mean')
    # seasonal daily mean Evap
    mEvap1_DJF = mEvap_clim1.where(DJF).resample('AS-Dec', 'dtime', how='mean')
    mEvap1_MAM = mEvap_clim1.where(MAM).resample('AS-Mar', 'dtime', how='mean')
    mEvap1_JJA = mEvap_clim1.where(JJA).resample('AS-Jun', 'dtime', how='mean')
    mEvap1_SON = mEvap_clim1.where(SON).resample('AS-Sep', 'dtime', how='mean')
    # seasonal daily mean SWE
    mSWE1_DJF = mSWE_clim1.where(DJF).resample('AS-Dec', 'dtime', how='mean')
    mSWE1_MAM = mSWE_clim1.where(MAM).resample('AS-Mar', 'dtime', how='mean')
    mSWE1_JJA = mSWE_clim1.where(JJA).resample('AS-Jun', 'dtime', how='mean')
    mSWE1_SON = mSWE_clim1.where(SON).resample('AS-Sep', 'dtime', how='mean')
    ################################################
    ###  data in a Dataframe for boxplot ###########
    ##################################################
    ## seasonal mean daily flow
    mQ_DJF = np.zeros(len(q_DJF)-1)
    mQ_MAM = np.zeros(len(q_MAM)-1)
    mQ_JJA = np.zeros(len(q_JJA)-1)
    mQ_SON = np.zeros(len(q_SON)-1)
    for j in range(1,len(q_MAM)):
        if j==1: mQ_DJF[j-1]=np.NaN
        else: mQ_DJF[j-1]=100.*(q1_DJF[j-1]-q_DJF[j-1])/q_DJF[j-1]   
        mQ_MAM[j-1]=100.*(q1_MAM[j]-q_MAM[j])/q_MAM[j]
        mQ_JJA[j-1]=100.*(q1_JJA[j]-q_JJA[j])/q_JJA[j]
        mQ_SON[j-1]=100.*(q1_SON[j]-q_SON[j])/q_SON[j]
    ## seasonal mean daily ET
    mmET_DJF = np.zeros(len(mET_DJF)-1)
    mmET_MAM = np.zeros(len(mET_MAM)-1)
    mmET_JJA = np.zeros(len(mET_JJA)-1)
    mmET_SON = np.zeros(len(mET_SON)-1)
    for j in range(1,len(mET_MAM)):
        if j==1: mmET_DJF[j-1]=np.NaN
        else: mmET_DJF[j-1]=100.*(mET1_DJF[j-1]-mET_DJF[j-1])/mET_DJF[j-1]   
        mmET_MAM[j-1]=100.*(mET1_MAM[j]-mET_MAM[j])/mET_MAM[j]
        mmET_JJA[j-1]=100.*(mET1_JJA[j]-mET_JJA[j])/mET_JJA[j]
        mmET_SON[j-1]=100.*(mET1_SON[j]-mET_SON[j])/mET_SON[j]
    ## seasonal mean daily Trans
    mmTr_DJF = np.zeros(len(mTrans_DJF)-1)
    mmTr_MAM = np.zeros(len(mTrans_MAM)-1)
    mmTr_JJA = np.zeros(len(mTrans_JJA)-1)
    mmTr_SON = np.zeros(len(mTrans_SON)-1)
    for j in range(1,len(mTrans_MAM)):
        if j==1: mmTr_DJF[j-1]=np.NaN
        else: mmTr_DJF[j-1]=100.*(mTrans1_DJF[j-1]-mTrans_DJF[j-1])/mTrans_DJF[j-1]   
        mmTr_MAM[j-1]=100.*(mTrans1_MAM[j]-mTrans_MAM[j])/mTrans_MAM[j]
        mmTr_JJA[j-1]=100.*(mTrans1_JJA[j]-mTrans_JJA[j])/mTrans_JJA[j]
        mmTr_SON[j-1]=100.*(mTrans1_SON[j]-mTrans_SON[j])/mTrans_SON[j] 
    ## seasonal mean daily Efloor
    mmEf_DJF = np.zeros(len(mEfloor_DJF)-1)
    mmEf_MAM = np.zeros(len(mEfloor_MAM)-1)
    mmEf_JJA = np.zeros(len(mEfloor_JJA)-1)
    mmEf_SON = np.zeros(len(mEfloor_SON)-1)
    for j in range(1,len(mEfloor_MAM)):
        if j==1: mmEf_DJF[j-1]=np.NaN
        else: mmEf_DJF[j-1]=100.*(mEfloor1_DJF[j-1]-mEfloor_DJF[j-1])/mEfloor_DJF[j-1]   
        mmEf_MAM[j-1]=100.*(mEfloor1_MAM[j]-mEfloor_MAM[j])/mEfloor_MAM[j]
        mmEf_JJA[j-1]=100.*(mEfloor1_JJA[j]-mEfloor_JJA[j])/mEfloor_JJA[j]
        mmEf_SON[j-1]=100.*(mEfloor1_SON[j]-mEfloor_SON[j])/mEfloor_SON[j]
    mmEf_DJF[mmEf_DJF>1e+03]=np.NaN
    mmEf_MAM[mmEf_MAM>1e+03]=np.NaN
    mmEf_JJA[mmEf_JJA>1e+03]=np.NaN
    mmEf_SON[mmEf_SON>1e+03]=np.NaN
    ## seasonal mean daily Evap
    mmEvap_DJF = np.zeros(len(mEvap_DJF)-1)
    mmEvap_MAM = np.zeros(len(mEvap_MAM)-1)
    mmEvap_JJA = np.zeros(len(mEvap_JJA)-1)
    mmEvap_SON = np.zeros(len(mEvap_SON)-1)
    for j in range(1,len(mEvap_MAM)):
        if j==1: mEvap_DJF[j-1]=np.NaN
        else: mmEvap_DJF[j-1]=100.*(mEvap1_DJF[j-1]-mEvap_DJF[j-1])/mEvap_DJF[j-1]   
        mmEvap_MAM[j-1]=100.*(mEvap1_MAM[j]-mEvap_MAM[j])/mEvap_MAM[j]
        mmEvap_JJA[j-1]=100.*(mEvap1_JJA[j]-mEvap_JJA[j])/mEvap_JJA[j]
        mmEvap_SON[j-1]=100.*(mEvap1_SON[j]-mEvap_SON[j])/mEvap_SON[j]
    ## seasonal mean daily SWE
    mmSWE_DJF = np.zeros(len(mSWE_DJF)-1)
    mmSWE_MAM = np.zeros(len(mSWE_MAM)-1)
    mmSWE_JJA = np.zeros(len(mSWE_JJA)-1)
    mmSWE_SON = np.zeros(len(mSWE_SON)-1)
    for j in range(1,len(mSWE_MAM)):
        if j==1: mSWE_DJF[j-1]=np.NaN
        else: mmSWE_DJF[j-1]=100.*(mSWE1_DJF[j-1]-mSWE_DJF[j-1])/mSWE_DJF[j-1]   
        mmSWE_MAM[j-1]=100.*(mSWE1_MAM[j]-mSWE_MAM[j])/mSWE_MAM[j]
        mmSWE_JJA[j-1]=100.*(mSWE1_JJA[j]-mSWE_JJA[j])/mSWE_JJA[j]
        mmSWE_SON[j-1]=100.*(mSWE1_SON[j]-mSWE_SON[j])/mSWE_SON[j]
       
    mAM= pd.DataFrame()
    mAM['Mean_Q_DJF'] = mQ_DJF
    mAM['Mean_Q_MAM'] = mQ_MAM
    mAM['Mean_Q_JJA'] = mQ_JJA
    mAM['Mean_Q_SON'] = mQ_SON 
    mAM['Mean_ET_DJF'] = mmET_DJF
    mAM['Mean_ET_MAM'] = mmET_MAM
    mAM['Mean_ET_JJA'] = mmET_JJA
    mAM['Mean_ET_SON'] = mmET_SON 
    mAM['Mean_Tr_DJF'] = mmTr_DJF
    mAM['Mean_Tr_MAM'] = mmTr_MAM
    mAM['Mean_Tr_JJA'] = mmTr_JJA
    mAM['Mean_Tr_SON'] = mmTr_SON 
    mAM['Mean_Ef_DJF'] = mmEf_DJF
    mAM['Mean_Ef_MAM'] = mmEf_MAM
    mAM['Mean_Ef_JJA'] = mmEf_JJA
    mAM['Mean_Ef_SON'] = mmEf_SON 
    mAM['Mean_Evap_DJF'] = mmEvap_DJF
    mAM['Mean_Evap_MAM'] = mmEvap_MAM
    mAM['Mean_Evap_JJA'] = mmEvap_JJA
    mAM['Mean_Evap_SON'] = mmEvap_SON 
    mAM['Mean_SWE_DJF'] = mmSWE_DJF
    mAM['Mean_SWE_MAM'] = mmSWE_MAM
    mAM['Mean_SWE_JJA'] = mmSWE_JJA
    mAM['Mean_SWE_SON'] = mmSWE_SON 
    
    ##### plotting boxplot
    sns.set_style('ticks')
    plt.figure(figsize=(8,4))
    sns.boxplot(data=mAM,showmeans=True,showfliers=False,linewidth=0.5,meanprops=dict(markersize=2))
    plt.ylabel('relative change (%)')
    plt.xticks(rotation='vertical')
    plt.axhline(linewidth=0.5)
    sns.despine()
    plt.tight_layout()
    plt.savefig('C:\Apps\WinPython-64bit-2.7.10.3\SAMULI_ARI\spathy\\results\\climate\\c_'+ID+'_'+scen+'_seasonal.pdf',dpi=300)
    
def Ann_Seasonal_Change(ID,scen,d_clim,d_clim1):
    #############################################
    ## variables from climate data simulations (1981-2010)
    #############################################
    q_clim = 1e3*d_clim.Qt
    #S_clim = d_clim.S
    mET_clim = (d_clim.Transpi+d_clim.Efloor+d_clim.Evap).mean(dim=['dlat','dlon'])
    mTrans_clim = d_clim.Transpi.mean(dim=['dlat','dlon'])
    mEfloor_clim = d_clim.Efloor.mean(dim=['dlat','dlon'])
    mEvap_clim = d_clim.Evap.mean(dim=['dlat','dlon'])
    mSWE_clim = d_clim.SWE.mean(dim=['dlat','dlon'])
    # mean daily value within a year
    mA_q_clim = q_clim.groupby('dtime.year').mean()
    #mA_S_clim = S_clim.groupby('dtime.year').mean()
    mA_mET_clim = mET_clim.groupby('dtime.year').mean()
    mA_mTrans_clim = mTrans_clim.groupby('dtime.year').mean()
    mA_mEfloor_clim = mEfloor_clim.groupby('dtime.year').mean()
    mA_mEvap_clim = mEvap_clim.groupby('dtime.year').mean()
    mA_mSWE_clim = mSWE_clim.groupby('dtime.year').mean()
    #############################################
    ## variables from climate data simulations (2071-2100)
    #############################################
    q_clim1 = 1e3*d_clim1.Qt
    #S_clim1 = d_clim1.S
    mET_clim1 = (d_clim1.Transpi+d_clim1.Efloor+d_clim1.Evap).mean(dim=['dlat','dlon'])
    mTrans_clim1 = d_clim1.Transpi.mean(dim=['dlat','dlon'])
    mEfloor_clim1 = d_clim1.Efloor.mean(dim=['dlat','dlon'])
    mEvap_clim1 = d_clim1.Evap.mean(dim=['dlat','dlon'])
    mSWE_clim1 = d_clim1.SWE.mean(dim=['dlat','dlon'])
    ### mean daily value within a year
    mA_q_clim1 = q_clim1.groupby('dtime.year').mean()
    #mA_S_clim1 = S_clim1.groupby('dtime.year').mean()
    mA_mET_clim1 = mET_clim1.groupby('dtime.year').mean()
    mA_mTrans_clim1 = mTrans_clim1.groupby('dtime.year').mean()
    mA_mEfloor_clim1 = mEfloor_clim1.groupby('dtime.year').mean()
    mA_mEvap_clim1 = mEvap_clim1.groupby('dtime.year').mean()
    mA_mSWE_clim1 = mSWE_clim1.groupby('dtime.year').mean()
    ################################################
    ###  seasonal mean daily ###########
    ##################################################
    month_clim = q_clim.groupby('dtime.month').apply(lambda x: x).month
    DJF = (month_clim <= 2) | (month_clim >= 12)
    MAM = (month_clim >= 3) & (month_clim <= 5)
    JJA = (month_clim >= 6) & (month_clim <= 8)
    SON = (month_clim >= 9) & (month_clim <= 11)
    # seasonal daily mean flow
    q_DJF = q_clim.where(DJF).resample('AS-Dec', 'dtime', how='mean')
    q_MAM = q_clim.where(MAM).resample('AS-Mar', 'dtime', how='mean')
    q_JJA = q_clim.where(JJA).resample('AS-Jun', 'dtime', how='mean')
    q_SON = q_clim.where(SON).resample('AS-Sep', 'dtime', how='mean')
    # seasonal daily mean ET
    mET_DJF = mET_clim.where(DJF).resample('AS-Dec', 'dtime', how='mean')
    mET_MAM = mET_clim.where(MAM).resample('AS-Mar', 'dtime', how='mean')
    mET_JJA = mET_clim.where(JJA).resample('AS-Jun', 'dtime', how='mean')
    mET_SON = mET_clim.where(SON).resample('AS-Sep', 'dtime', how='mean')
    # seasonal daily mean Transpi
    mTrans_DJF = mTrans_clim.where(DJF).resample('AS-Dec', 'dtime', how='mean')
    mTrans_MAM = mTrans_clim.where(MAM).resample('AS-Mar', 'dtime', how='mean')
    mTrans_JJA = mTrans_clim.where(JJA).resample('AS-Jun', 'dtime', how='mean')
    mTrans_SON = mTrans_clim.where(SON).resample('AS-Sep', 'dtime', how='mean')
    # seasonal daily mean Efloor
    mEfloor_DJF = mEfloor_clim.where(DJF).resample('AS-Dec', 'dtime', how='mean')
    mEfloor_MAM = mEfloor_clim.where(MAM).resample('AS-Mar', 'dtime', how='mean')
    mEfloor_JJA = mEfloor_clim.where(JJA).resample('AS-Jun', 'dtime', how='mean')
    mEfloor_SON = mEfloor_clim.where(SON).resample('AS-Sep', 'dtime', how='mean')
    # seasonal daily mean Evap
    mEvap_DJF = mEvap_clim.where(DJF).resample('AS-Dec', 'dtime', how='mean')
    mEvap_MAM = mEvap_clim.where(MAM).resample('AS-Mar', 'dtime', how='mean')
    mEvap_JJA = mEvap_clim.where(JJA).resample('AS-Jun', 'dtime', how='mean')
    mEvap_SON = mEvap_clim.where(SON).resample('AS-Sep', 'dtime', how='mean')
    # seasonal daily mean SWE
    mSWE_DJF = mSWE_clim.where(DJF).resample('AS-Dec', 'dtime', how='mean')
    mSWE_MAM = mSWE_clim.where(MAM).resample('AS-Mar', 'dtime', how='mean')
    mSWE_JJA = mSWE_clim.where(JJA).resample('AS-Jun', 'dtime', how='mean')
    mSWE_SON = mSWE_clim.where(SON).resample('AS-Sep', 'dtime', how='mean')
    #####################################    
    # seasonal daily mean flow for 2070-2099
    #####################################
    month_clim = q_clim1.groupby('dtime.month').apply(lambda x: x).month
    DJF = (month_clim <= 2) | (month_clim >= 12)
    MAM = (month_clim >= 3) & (month_clim <= 5)
    JJA = (month_clim >= 6) & (month_clim <= 8)
    SON = (month_clim >= 9) & (month_clim <= 11)
    q1_DJF = q_clim1.where(DJF).resample('AS-Dec', 'dtime', how='mean')
    q1_MAM = q_clim1.where(MAM).resample('AS-Mar', 'dtime', how='mean')
    q1_JJA = q_clim1.where(JJA).resample('AS-Jun', 'dtime', how='mean')
    q1_SON = q_clim1.where(SON).resample('AS-Sep', 'dtime', how='mean')
    # seasonal daily mean ET
    mET1_DJF = mET_clim1.where(DJF).resample('AS-Dec', 'dtime', how='mean')
    mET1_MAM = mET_clim1.where(MAM).resample('AS-Mar', 'dtime', how='mean')
    mET1_JJA = mET_clim1.where(JJA).resample('AS-Jun', 'dtime', how='mean')
    mET1_SON = mET_clim1.where(SON).resample('AS-Sep', 'dtime', how='mean')
    # seasonal daily mean Transpi
    mTrans1_DJF = mTrans_clim1.where(DJF).resample('AS-Dec', 'dtime', how='mean')
    mTrans1_MAM = mTrans_clim1.where(MAM).resample('AS-Mar', 'dtime', how='mean')
    mTrans1_JJA = mTrans_clim1.where(JJA).resample('AS-Jun', 'dtime', how='mean')
    mTrans1_SON = mTrans_clim1.where(SON).resample('AS-Sep', 'dtime', how='mean')
    # seasonal daily mean Efloor
    mEfloor1_DJF = mEfloor_clim1.where(DJF).resample('AS-Dec', 'dtime', how='mean')
    mEfloor1_MAM = mEfloor_clim1.where(MAM).resample('AS-Mar', 'dtime', how='mean')
    mEfloor1_JJA = mEfloor_clim1.where(JJA).resample('AS-Jun', 'dtime', how='mean')
    mEfloor1_SON = mEfloor_clim1.where(SON).resample('AS-Sep', 'dtime', how='mean')
    # seasonal daily mean Evap
    mEvap1_DJF = mEvap_clim1.where(DJF).resample('AS-Dec', 'dtime', how='mean')
    mEvap1_MAM = mEvap_clim1.where(MAM).resample('AS-Mar', 'dtime', how='mean')
    mEvap1_JJA = mEvap_clim1.where(JJA).resample('AS-Jun', 'dtime', how='mean')
    mEvap1_SON = mEvap_clim1.where(SON).resample('AS-Sep', 'dtime', how='mean')
    # seasonal daily mean SWE
    mSWE1_DJF = mSWE_clim1.where(DJF).resample('AS-Dec', 'dtime', how='mean')
    mSWE1_MAM = mSWE_clim1.where(MAM).resample('AS-Mar', 'dtime', how='mean')
    mSWE1_JJA = mSWE_clim1.where(JJA).resample('AS-Jun', 'dtime', how='mean')
    mSWE1_SON = mSWE_clim1.where(SON).resample('AS-Sep', 'dtime', how='mean')
    ################################################
    ###  data in a Dataframe for boxplot ###########
    ##################################################
    mAQ = np.zeros(len(mA_q_clim))
    mAET = np.zeros(len(mA_mET_clim))
    mATr = np.zeros(len(mA_mTrans_clim))
    mAEf = np.zeros(len(mA_mEfloor_clim))
    mAEv = np.zeros(len(mA_mEvap_clim))
    mASWE = np.zeros(len(mA_mSWE_clim))
    for j in range(len(mA_q_clim)):
        mAQ[j]=100.*(mA_q_clim1[j]-mA_q_clim[j])/mA_q_clim[j]
        mAET[j]=100.*(mA_mET_clim1[j]-mA_mET_clim[j])/mA_mET_clim[j]
        mATr[j]=100.*(mA_mTrans_clim1[j]-mA_mTrans_clim[j])/mA_mTrans_clim[j]
        mAEf[j]=100.*(mA_mEfloor_clim1[j]-mA_mEfloor_clim[j])/mA_mEfloor_clim[j]
        mAEv[j]=100.*(mA_mEvap_clim1[j]-mA_mEvap_clim[j])/mA_mEvap_clim[j]
        mASWE[j]=100.*(mA_mSWE_clim1[j]-mA_mSWE_clim[j])/mA_mSWE_clim[j]
    ## seasonal mean daily flow
    mQ_DJF = np.zeros(len(q_DJF)-1)
    mQ_MAM = np.zeros(len(q_MAM)-1)
    mQ_JJA = np.zeros(len(q_JJA)-1)
    mQ_SON = np.zeros(len(q_SON)-1)
    for j in range(1,len(q_MAM)):
        if j==1: mQ_DJF[j-1]=np.NaN
        else: mQ_DJF[j-1]=100.*(q1_DJF[j-1]-q_DJF[j-1])/q_DJF[j-1]   
        mQ_MAM[j-1]=100.*(q1_MAM[j]-q_MAM[j])/q_MAM[j]
        mQ_JJA[j-1]=100.*(q1_JJA[j]-q_JJA[j])/q_JJA[j]
        mQ_SON[j-1]=100.*(q1_SON[j]-q_SON[j])/q_SON[j]
    ## seasonal mean daily ET
    mmET_DJF = np.zeros(len(mET_DJF)-1)
    mmET_MAM = np.zeros(len(mET_MAM)-1)
    mmET_JJA = np.zeros(len(mET_JJA)-1)
    mmET_SON = np.zeros(len(mET_SON)-1)
    for j in range(1,len(mET_MAM)):
        if j==1: mmET_DJF[j-1]=np.NaN
        else: mmET_DJF[j-1]=100.*(mET1_DJF[j-1]-mET_DJF[j-1])/mET_DJF[j-1]   
        mmET_MAM[j-1]=100.*(mET1_MAM[j]-mET_MAM[j])/mET_MAM[j]
        mmET_JJA[j-1]=100.*(mET1_JJA[j]-mET_JJA[j])/mET_JJA[j]
        mmET_SON[j-1]=100.*(mET1_SON[j]-mET_SON[j])/mET_SON[j]
    ## seasonal mean daily Trans
    mmTr_DJF = np.zeros(len(mTrans_DJF)-1)
    mmTr_MAM = np.zeros(len(mTrans_MAM)-1)
    mmTr_JJA = np.zeros(len(mTrans_JJA)-1)
    mmTr_SON = np.zeros(len(mTrans_SON)-1)
    for j in range(1,len(mTrans_MAM)):
        if j==1: mmTr_DJF[j-1]=np.NaN
        else: mmTr_DJF[j-1]=100.*(mTrans1_DJF[j-1]-mTrans_DJF[j-1])/mTrans_DJF[j-1]   
        mmTr_MAM[j-1]=100.*(mTrans1_MAM[j]-mTrans_MAM[j])/mTrans_MAM[j]
        mmTr_JJA[j-1]=100.*(mTrans1_JJA[j]-mTrans_JJA[j])/mTrans_JJA[j]
        mmTr_SON[j-1]=100.*(mTrans1_SON[j]-mTrans_SON[j])/mTrans_SON[j] 
    ## seasonal mean daily Efloor
    mmEf_DJF = np.zeros(len(mEfloor_DJF)-1)
    mmEf_MAM = np.zeros(len(mEfloor_MAM)-1)
    mmEf_JJA = np.zeros(len(mEfloor_JJA)-1)
    mmEf_SON = np.zeros(len(mEfloor_SON)-1)
    for j in range(1,len(mEfloor_MAM)):
        if j==1: mmEf_DJF[j-1]=np.NaN
        else: mmEf_DJF[j-1]=100.*(mEfloor1_DJF[j-1]-mEfloor_DJF[j-1])/mEfloor_DJF[j-1]   
        mmEf_MAM[j-1]=100.*(mEfloor1_MAM[j]-mEfloor_MAM[j])/mEfloor_MAM[j]
        mmEf_JJA[j-1]=100.*(mEfloor1_JJA[j]-mEfloor_JJA[j])/mEfloor_JJA[j]
        mmEf_SON[j-1]=100.*(mEfloor1_SON[j]-mEfloor_SON[j])/mEfloor_SON[j]
    mmEf_DJF[mmEf_DJF>1e+03]=np.NaN
    mmEf_MAM[mmEf_MAM>1e+03]=np.NaN
    mmEf_JJA[mmEf_JJA>1e+03]=np.NaN
    mmEf_SON[mmEf_SON>1e+03]=np.NaN
    ## seasonal mean daily Evap
    mmEvap_DJF = np.zeros(len(mEvap_DJF)-1)
    mmEvap_MAM = np.zeros(len(mEvap_MAM)-1)
    mmEvap_JJA = np.zeros(len(mEvap_JJA)-1)
    mmEvap_SON = np.zeros(len(mEvap_SON)-1)
    for j in range(1,len(mEvap_MAM)):
        if j==1: mEvap_DJF[j-1]=np.NaN
        else: mmEvap_DJF[j-1]=100.*(mEvap1_DJF[j-1]-mEvap_DJF[j-1])/mEvap_DJF[j-1]   
        mmEvap_MAM[j-1]=100.*(mEvap1_MAM[j]-mEvap_MAM[j])/mEvap_MAM[j]
        mmEvap_JJA[j-1]=100.*(mEvap1_JJA[j]-mEvap_JJA[j])/mEvap_JJA[j]
        mmEvap_SON[j-1]=100.*(mEvap1_SON[j]-mEvap_SON[j])/mEvap_SON[j]
    ## seasonal mean daily SWE
    mmSWE_DJF = np.zeros(len(mSWE_DJF)-1)
    mmSWE_MAM = np.zeros(len(mSWE_MAM)-1)
    mmSWE_JJA = np.zeros(len(mSWE_JJA)-1)
    mmSWE_SON = np.zeros(len(mSWE_SON)-1)
    for j in range(1,len(mSWE_MAM)):
        if j==1: mSWE_DJF[j-1]=np.NaN
        else: mmSWE_DJF[j-1]=100.*(mSWE1_DJF[j-1]-mSWE_DJF[j-1])/mSWE_DJF[j-1]   
        mmSWE_MAM[j-1]=100.*(mSWE1_MAM[j]-mSWE_MAM[j])/mSWE_MAM[j]
        mmSWE_JJA[j-1]=100.*(mSWE1_JJA[j]-mSWE_JJA[j])/mSWE_JJA[j]
        mmSWE_SON[j-1]=100.*(mSWE1_SON[j]-mSWE_SON[j])/mSWE_SON[j]
       
    mAM= pd.DataFrame()
    mAM['Mean_Q_ANN'] = mAQ
    mAM['Mean_Q_DJF'] = mQ_DJF
    mAM['Mean_Q_MAM'] = mQ_MAM
    mAM['Mean_Q_JJA'] = mQ_JJA
    mAM['Mean_Q_SON'] = mQ_SON 
    mAM['Mean_ET_ANN'] = mAET
    mAM['Mean_ET_DJF'] = mmET_DJF
    mAM['Mean_ET_MAM'] = mmET_MAM
    mAM['Mean_ET_JJA'] = mmET_JJA
    mAM['Mean_ET_SON'] = mmET_SON 
    mAM['Mean_Tr_ANN'] = mATr
    mAM['Mean_Tr_DJF'] = mmTr_DJF
    mAM['Mean_Tr_MAM'] = mmTr_MAM
    mAM['Mean_Tr_JJA'] = mmTr_JJA
    mAM['Mean_Tr_SON'] = mmTr_SON 
    mAM['Mean_Ef_ANN'] = mAEf
    mAM['Mean_Ef_DJF'] = mmEf_DJF
    mAM['Mean_Ef_MAM'] = mmEf_MAM
    mAM['Mean_Ef_JJA'] = mmEf_JJA
    mAM['Mean_Ef_SON'] = mmEf_SON 
    mAM['Mean_Evap_ANN'] = mAEv
    mAM['Mean_Evap_DJF'] = mmEvap_DJF
    mAM['Mean_Evap_MAM'] = mmEvap_MAM
    mAM['Mean_Evap_JJA'] = mmEvap_JJA
    mAM['Mean_Evap_SON'] = mmEvap_SON 
    mAM['Mean_SWE_ANN'] = mASWE
    mAM['Mean_SWE_DJF'] = mmSWE_DJF
    mAM['Mean_SWE_MAM'] = mmSWE_MAM
    mAM['Mean_SWE_JJA'] = mmSWE_JJA
    mAM['Mean_SWE_SON'] = mmSWE_SON 
    
    ##### plotting boxplot
    sns.set_style('ticks')
    plt.figure(figsize=(8,4))
    sns.boxplot(data=mAM,showmeans=True,showfliers=False,linewidth=0.5,meanprops=dict(markersize=2))
    plt.ylabel('relative change (%)')
    plt.xticks(rotation='vertical')
    plt.axhline(linewidth=0.5)
    sns.despine()
    plt.tight_layout()
    plt.savefig('C:\Apps\WinPython-64bit-2.7.10.3\SAMULI_ARI\spathy\\results\\climate\\c_'+ID+'_'+scen+'_ann_seasonal.pdf',dpi=300)