# -*- coding: utf-8 -*-
"""
Created on Mon Dec 03 21:41:07 2018

@author: slauniai
"""

import os
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
#import seaborn as sns
from scipy import stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
from mpl_toolkits.axes_grid1 import make_axes_locatable



from canopygrid import eq_evap
eps = np.finfo(float).eps
#spathy_path = os.path.join('c:', 'c:\datat\spathydata')
#results_path = os.path.join(r'c:\modelresults\spathy')
#fig_path = os.path.join(r'c:\modelresults\spathy\figs')

dt = 86400.0


def calculate_pvalues(df):
    # computes p-values of dataframe
    df = df.dropna()._get_numeric_data()
    dfcols = pd.DataFrame(columns=df.columns)
    pvalues = dfcols.transpose().join(dfcols, how='outer')
    for r in df.columns:
        for c in df.columns:
            pvalues[r][c] = round(stats.pearsonr(df[r], df[c])[1], 4)
    return pvalues

def modified_agreementindex(evaluation, simulation):
    """ modified willmot's agreement index"""
    if len(evaluation) == len(simulation):
        s, e = np.array(simulation), np.array(evaluation)
        # s, e = simulation, evaluation
        
        # compute numerator and denominator
        numerator = sum(abs(e - s))
        denominator = sum(abs(s - np.mean(e)) + abs(e - np.mean(e)))
        obj = 1 - (numerator / (denominator + eps))
        #print obj    
        return obj

ch = [1, 2, 3, 10, 11, 13, 14, 16, 17, 19, 20, 21, 24, 25, 26, 27, 28, 29, 31, 32, 33]
#ch = [1, 3, 10, 11, 13, 14, 16, 17, 19, 20, 21, 24, 25, 26, 27, 28, 31, 33]
yrs = [[range(2014,2016)],
       [[2006, 2008]],
       [range(2006, 2016)],
       [range(2003, 2014)],
       [[2015]],
       [[2015]],
       [range(2007, 2016)],
       [range(2007, 2016)],       
       [range(2006, 2016)],
       [[2008, 2009, 2010, 2011, 2013, 2014, 2015]],
       [range(2006, 2016)],
       [[2006, 2007, 2008, 2009, 2010, 2012, 2013, 2014, 2015]],
       #[range(2006, 2016)],
       [range(2006, 2016)],
       [range(2006, 2016)],
       [range(2006, 2016)],
       [range(2006, 2016)],
       [range(2014, 2016)],
       [[2015]],
       # [range(2012, 2016)],
       [range(2012, 2016)],
       [range(2013, 2016)],
       [[2006, 2007, 2008, 2009, 2012, 2014, 2015]],
       ]
# read chdata from csv file
# chdata = pd.read_csv(r'c:\datat\spathydata\sve_catchment_characteristics.csv', sep=';')

# read data into pickle
d0 = pickle.load(open(r'c:\repositories\Spafhy_v1\results\R-base.pkl', 'rb'))  # baseline
dat0 = []
for k in range(len(d0)):
    dat0.append(pd.DataFrame(data=d0[k], columns=d0[k].keys(), index=d0[k]['Qmeas'].index))

Dj = []
for k in range(0, len(ch)):
    x = dat0[k]
    x = x[x.index.year.isin(yrs[k][0])]
    plt.figure(ch[k])
    plt.plot(x['Qmeas'], 'k-', x['Qt'], 'r-')
    
    y = x['Qmeas'].values
    z = x['Qt'].values
    
    f = np.where(np.isfinite(y))
    Dj.append([ch[k], modified_agreementindex(y[f], z[f])])
    #plt.title(str(Dj[k][1]))

#%% plot all catchments into a subplot
import matplotlib.dates as mdates

fig, ax = plt.subplots(6,3)
fig.set_size_inches(20, 12.5)

yr = np.arange(2013, 2016)
n = 0; m = 0
#ch = [1, 2, 3, 10, 11, 13, 14, 16, 17, 19, 20, 21, 24, 25, 26, 27, 28, 29, 31, 32, 33]
kk = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20]

for k in kk:
    x = dat0[k]
    ix = np.intersect1d(yr, yrs[k][0])
    x = x[x.index.year.isin(ix)]
    y = x['Qmeas'].values
    z = x['Qt'].values
    f = np.where(np.isfinite(y))
    
    objf = modified_agreementindex(y[f], z[f])
    txt = r'C%d: %.2f' %(ch[k], objf)
    ax[n,m].plot(x['Qmeas'], 'k-', x['Qt'].iloc[f], 'r-')
    ax[n,m].text(0.47, 0.85, txt, fontsize=10, transform=ax[n,m].transAxes)
    ax[n,m].set_xlim(['2013-01-01', '2016-01-01'])

    ax[n,m].xaxis.set_major_formatter(mdates.DateFormatter("%y-%m"))
    ax[n,m].xaxis.set_minor_formatter(mdates.DateFormatter("%y-%m"))
    ax[n,m].set_ylabel('Q (mm d$^{-1}$)')

    n +=1
    if n == 6:
        n = 0
        m +=1
        
plt.savefig('All_C_discharge.png', dpi=600)
        