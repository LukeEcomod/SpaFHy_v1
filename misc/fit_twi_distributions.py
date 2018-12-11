# -*- coding: utf-8 -*-
"""

Fit distributions to TWI and LAI -data at each catchments (LAI todo)
Created on Thu Sep 14 15:57:46 2017

@author: slauniai
"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from scipy import stats
from iotools import create_catchment

eps = np.finfo(float).eps
fpath = r'c:\datat\spathydata\svecatchments'
plt.figure()
k = 1
# read catchment data and compute twi
for catchment in [1, 2, 3, 10, 11, 13, 14, 16, 17, 19, 20, 21, 22, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33]:
    gis = create_catchment(str(catchment), fpath)
#    flowacc = gis['flowacc']
#    slope = np.radians(gis['slope'])  # deg to rad
#    
#    twi = np.log(flowacc / (np.tan(slope) + eps))
    twi = gis['twi']
    twi = twi[twi>0].copy()
    
    # plot distribution and fit log-normal and fisk
    dist_names = ['lognorm', 'fisk']
    param = []
    
    sns.set_style('whitegrid')
    y = twi
    x = np.linspace(min(y), max(y), 100)
    with sns.color_palette('muted'):
        plt.subplot(3,3,k)
        plt.hist(y, 100, alpha=0.6, normed=True)
    
        for dist_name in dist_names:
            dist = getattr(stats, dist_name)
            p = dist.fit(y)
            yhat = dist.pdf(x, p[0], loc=p[1], scale=p[2])
            plt.plot(x, yhat, '-', label=dist_name)
            plt.legend()
            plt.title('C ' + str(catchment))
            plt.xlim([min(y), np.percentile(y, 99.0)])
            plt.xlabel('twi')
            param.append({'id': catchment, 'distr': dist_name, 'param': p})
        k += 1
        if k == 9:
            plt.figure()
            k = 1
            
    #        pdf_fitted = dist.pdf(x, *param[:-2], loc=param[-2], scale=param[-1]) * size
    #        plt.plot(pdf_fitted, label=dist_name)
    #        plt.legend(loc='upper right')
    #        plt.show()