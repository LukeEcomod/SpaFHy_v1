# -*- coding: utf-8 -*-
"""

MODULE FOR SPATIAL VARIABILITY OF SOIL MOISTURE
Created on Tue May 15 14:49:04 2018

@author: slauniai
"""

import numpy as np
from scipy.stats import spearmanr, rankdata
import os
import matplotlib.pyplot as plt
# from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# import matplotlib.colors as mplcolors
# import matplotlib.cm as mplcm
import pickle
# from netCDF4 import Dataset, date2num
# from spathy_sve import spathy_driver

eps = np.finfo(float).eps

# ou = os.path.join(r'c:\ModelResults\Spathy', 'moisturebudget_data.pkl')
# data = pickle.load(open(ou, 'rb'))
# soil = data['soil']

""" function definitions """


def variance_budget(data, k=30, Lz = 0.4):
    def spatial_deviation(x):
        xprime = np.ones(np.shape(x)) * np.NaN
        xave = np.mean(x, axis=1)
        
        for j in range(0, len(xave)):
            xprime[j,:] = x[j,:] - xave[j]
        return xprime
    
    def spatial_variance(x):
        return np.var(x, axis=1)


    #data = [];  # dimensions: 0 = time, 1=x, 2=y

    ix_x, ix_y = np.where(data['Wliq'][1,:,:] >= 0)
    #ix_x, ix_y = np.where(data['soil'] == 2)
    
    time = data['tvec'] # time
    W = data['Wliq'][:,ix_x, ix_y] # moisture
    I = 1e-3*data['I'] [:,ix_x, ix_y] #infiltration m/d
    T = 1e-3*data['Tr'] [:,ix_x, ix_y] #transpiration m/d
    D = 1e-3*data['D'] [:,ix_x, ix_y] #drainage m/d
    R = 1e-3*data['R'] [:,ix_x, ix_y] #returnflow m/d
    
    # perturbations around spatial average
    w = spatial_deviation(W)
    i = spatial_deviation(I)
    t = spatial_deviation(T)
    d = spatial_deviation(D)
    r = spatial_deviation(R)

    # instantaneous products
    # wi = w*i
    # wt = w*t
    # wd = w*d
    # wr = w*r
    
    # soil moisture variance budget
    m = len(time) / k
    ww = np.zeros(m) + np.NaN
    dww = np.zeros(m) + np.NaN
    wi = np.zeros(m) + np.NaN
    wt = np.zeros(m) + np.NaN
    wd = np.zeros(m) + np.NaN
    wr = np.zeros(m) + np.NaN
    
    n = 0
    for j in range(0, m):
        # ww[j] = np.mean(w[n:n+k,:] * w[n:n+k,:])
        ww[j] = np.var(w[n:n+k,:])
        wi[j] = np.mean(w[n:n+k,:] * i[n:n+k,:])
        wt[j] = np.mean(w[n:n+k,:] * t[n:n+k,:])
        wd[j] = np.mean(w[n:n+k,:] * d[n:n+k,:])
        wr[j] = np.mean(w[n:n+k,:] * r[n:n+k,:])
        
        dww[j] = 2 / Lz * (wi[j] - wt[j] - wd[j] + wr[j])
        n += k
    
    
    # plot figure
    
    plt.figure()
    plt.subplot(211)
    plt.plot(ww, 'ko-'); plt.ylabel('ww & dww')
    plt.subplot(212)
    plt.plot(2/Lz*wi, 'co-', label='wi')
    plt.plot(-2/Lz*wt, 'yo-', label='wt')
    plt.plot(-2/Lz*wd, 'ks--', label='wd')
    plt.plot(+2/Lz*wr, 'bs--', label='wr')
    plt.plot(dww, 'ro-', label='dww')
    plt.plot(np.diff(ww)/k, 'go-')
    plt.legend()
    return dww, ww, wi, wt, wd, wr


def time_stability(w):
    # see Williamns et al. 2009 HESS eq. 1 - 6
    m, i = np.shape(w)  # timesteps, nodes
    w = np.around(w, decimals=2)
    
    delta, delta_ave, delta_std = relative_difference(w)
    rci = rank_change_index(w)
    sigma = np.std(w, axis=1)
    wm = np.mean(w, axis=1)
    return delta, delta_ave, delta_std, rci, wm, sigma, sigma/wm

def relative_difference(x):
    # relative difference and mean relative difference at point i
    delta = np.ones(np.shape(x)) * np.NaN
    xave = np.mean(x, axis=1)
    mm = len(xave)
    for j in range(0, mm):
        delta[j,:] = x[j,:] - xave[j]
    delta_ave = np.mean(delta, axis=0)
    delta_std = sum((delta - delta_ave) / (mm - 1))**0.5
    return delta, delta_ave, delta_std

def spatial_std(x):
    # std at each time point
    return np.std(x, axis=1)
        
def rank_correlation(w, ix0, ix1):
    # rank correlation between time points ix0, ix1
    r, p = spearmanr(w[ix0,:], w[ix1,:])
    return r, p

def rank_change_index(w):
    mm, kk = np.shape(w)
    
    rci = np.zeros(kk)
    Ro = rankdata(w[0,:])
    for j in range(1, mm):
        Rj = rankdata(w[j,:])
        for i in range(0, kk):
            rci[i] += abs(Rj[i] - Ro[i])
        Ro = Rj
    return rci
              
## call functions
#    
## dww, ww, wi, wt, wd, wr = variance_budget(data, k=7, Lz=0.4)
#ix_x, ix_y = np.where(data['Wliq'][1,:,:] >= 0)
##ix_x, ix_y = np.where(data['soil'] == 2)
#
#time = data['tvec'] # time
#W = data['Wliq'][:,ix_x, ix_y] # moisture
#
#delta, delta_ave, delta_std, rci, wm, sigma, cv = time_stability(W)      
#rci = rci / max(rci)  # normalize peak to 1
#plt.figure()
#plt.subplot(221); plt.plot(wm, sigma, 'ko'); plt.xlabel(r'$\langle \theta \rangle$'); plt.ylabel(r'$\sigma$')
#plt.subplot(222); plt.plot(wm, cv, 'ko'); plt.xlabel(r'$\langle \theta \rangle$'); plt.ylabel(r'$CV_{\theta}$')
#plt.subplot(223); plt.plot(delta_ave, rci, 'ko'); plt.xlabel(r'MRD'); plt.ylabel('RCI')
#plt.subplot(224); plt.plot(delta_ave, delta_std, 'ko'); plt.xlabel(r'MRD'); plt.ylabel(r'$\sigma_{MRD}$')
