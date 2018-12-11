# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 18:14:00 2018

@author: slauniai
"""

from iotools import create_catchment

import os
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

eps = np.finfo(float).eps
fpath = r'c:\datat\spathydata\svecatchments'

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


def fit_distr(y, dist_name='norm', plot_figs=False):
    #dist_name = norm, lognorm / fisk
    dist = getattr(stats, dist_name)
    p = dist.fit(y)
    
    if plot_figs:
        plt.figure()
        x = np.linspace(min(y), max(y), 100)
        yhat = dist.pdf(x, p[0], loc=p[1], scale=p[2])
        plt.plot(x, yhat, '-', label=dist_name)
        plt.title('loc=%.2f, scale=%.2f' % (p[1], p[2]))
        plt.xlim([min(y), np.percentile(y, 99.0)])
    
    return p

# cols ['slope',
# 'peatm',
# 'stream',
# 'vol',
# 'soil',
# 'cf',
# 'LAI_decid',
# 'flowacc',
# 'loc',
# 'sitetype',
# 'twi',
# 'maintype',
# 'cmask',
# 'LAI_conif',
# 'bmroot',
# 'ba',
# 'hc',
# 'LAI_spruce',
# 'info',
# 'lat0',
# 'age',
# 'cellsize',
# 'rockm',
# 'LAI_pine',
# 'dem',
# 'lon0']

A = {'id':[], 'lat': [], 'lon': [], 'area': [], 'TWI': [], 'TWI_median': [], 'TWI_b': [],
       'peat': [], 'coarse': [], 'med': [], 'fine': [], 'vol': [], 'LAI': [],
       'LAI_median':[], 'LAI_std': [], 'f_decid': [], 'f_decid_median': [], 'top_m': []
       }
for k in chm:
    gis = create_catchment(k[0],fpath)
    
    A['id'].append(int(k[0]))
    A['top_m'].append(k[-1])
    A['lat'].append(float(gis['loc']['lat']))
    A['lon'].append(float(gis['loc']['lon']))
    A['area'].append(np.nansum(gis['cmask'])*16*16 / (100 *100.0))  # ha
    
    #twi distribution
    
    y = gis['twi']
    y = y[y>0].copy()
    y = y[np.isfinite(y)]
    A['TWI'].append(np.nanmean(y))
    A['TWI_median'].append(np.nanmedian(y))
    p = fit_distr(np.log(y), dist_name='norm')
    A['TWI_b'].append(p[1])  # sigma of normal distr

    A['vol'].append(np.nanmean(gis['vol']))
    laid = gis['LAI_decid']
    laic = gis['LAI_conif']
    lai = laic + laid
    A['LAI'].append(np.nanmean(lai))
    A['LAI_median'].append(np.nanmedian(lai))
    A['LAI_std'].append(np.nanstd(lai))
    A['f_decid'].append(np.nanmean(laid / lai))
    A['f_decid_median'].append(np.nanmedian(laid / lai))
    
    # compute proportion of soil types
    
    cmask = gis['cmask']
    ix = np.where(cmask >0)
    n = float(len(ix[0]))  # gridcells within catchment
    
    soil = gis['soil']
    soil = soil[ix]
    
    x, cnts = np.unique(soil, return_counts=True)
    
    if 4 in x:
        A['peat'].append(float(cnts[x==4] / n))
    else:
        A['peat'].append(0.0)
    if 3 in x:
        A['fine'].append(float(cnts[x==3] / n))
    else:
        A['fine'].append(0.0)
    if 2 in x:
        A['med'].append(float(cnts[x==2] / n))
    else:
        A['med'].append(0.0)
    if 1 in x:
        A['coarse'].append(float(cnts[x==1] / n))
    else:
        A['coarse'].append(0.0) 
    
    #seek open peatlands
#    lai = lai[ix]
#    laid = laid[ix]
#    h = gis['hc']
#    h = h[ix]
#    f = ((soil == 4) & (laid > 0.1) & (h == 0.1))
#    dd = lai[f]
#    print len(f)
#    if len(f) > 0:
#        A['open_peat'].append(float(len(dd) / n))
#    else:
#        A['open_peat'] = 0.0

A = pd.DataFrame(A, index=A['id'])
    
A.to_csv('c:\datat\spathydata\catchment_characteristics.csv', sep=';')

