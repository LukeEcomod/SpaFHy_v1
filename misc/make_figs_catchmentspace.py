# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 20:49:51 2018

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

dt = 86400.0

fig_path = r'C:\repositories\SpaFHy\FigsC3'

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
        return obj


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
#     ['22', '2005-01-01', '2015-12-31', '2005-12-31', 0.0095],  # vaha-askanjoki 1600 ha
#    [ '23','2011-01-01', '2015-12-31', '2011-12-31'],  # ylijoki 5600 ha, very large and slow
     ['24', '2005-01-01', '2015-12-31', '2005-12-31', 0.0066],  # kotioja 1800 ha
     ['25', '2005-01-01', '2015-12-31', '2005-12-31', 0.0095],  # kohisevanpuro 1070 ha
     ['26', '2005-01-01', '2015-12-31', '2005-12-31', 0.02],  # iittovuoma 1160 ha
     ['27', '2005-01-01', '2015-12-31', '2005-12-31', 0.014],  # laanioja 1362 ha
     ['28', '2013-01-01', '2015-12-31', '2013-12-31', 0.0057],  # kroopinsuo 179 ha
     ['29', '2012-01-01', '2015-12-31', '2012-12-31', 0.0089],  # surnui 71 ha, poor data quality
#    ['30', '2011-01-01', '2015-12-31', '2011-12-31', 0.0064],  # pakopirtti 795 ha, uncertain catchment boundaries
     ['31', '2011-01-01', '2015-12-31', '2011-12-31', 0.0064],  # ojakorpi 33 ha
     ['32', '2011-01-01', '2015-12-31', '2011-12-31', 0.0077],  # rantainrahka 38 ha
     ['33', '2005-01-01', '2015-12-31', '2005-12-31', 0.009],  # kivipuro 54 ha
     ]

ch = [1, 2, 3, 10, 11, 13, 14, 16, 17, 19, 20, 21, 24, 25, 26, 27, 28, 29,
      31, 32, 33]

# read chdata from csv file
chdata = pd.read_csv('sve_catchment_characteristics.csv', sep=';', encoding='ansi')

# read data into pickle
d0 = pickle.load(open('R-base.pkl', 'rb'))  # baseline
d1 = pickle.load(open('R-lowET.pkl', 'rb'))  # low et scenario
d2 = pickle.load(open('R-hiET.pkl', 'rb'))  # hi et scenario
dat0 = []; dat1 = []; dat2 = []

for k in range(len(d0)):
    dat0.append(pd.DataFrame(data=d0[k], columns=d0[k].keys(), index=d0[k]['Qmeas'].index))
    dat1.append(pd.DataFrame(data=d1[k], columns=d1[k].keys(), index=d1[k]['Qmeas'].index))
    dat2.append(pd.DataFrame(data=d2[k], columns=d2[k].keys(), index=d2[k]['Qmeas'].index))    

del d0, d1, d2

#%%
# compute annual and long-term ET, Q, P and Qm
An = pd.DataFrame(columns=['year', 'ch', 'Prec', 'Qmeas', 'Qt', 'Qt_hi', 'Qt_lo', 'ET', 'ET_hi', 'ET_lo', 'ET0', 'fsnow'], dtype=np.float)
A =  pd.DataFrame(columns=['year', 'ch', 'Prec', 'Qmeas', 'Qt', 'Qt_hi', 'Qt_lo', 'ET', 'ET_hi', 'ET_lo', 'ET0', 'fsnow'], dtype=np.float)

for k in range(0, len(ch)):
    # compute Equilibrium ET and fraction of P falling as snow
    rn = 0.7 * dat0[k]['Rg'].values
    ta = dat0[k]['Ta'].values

    Eo = np.zeros(len(ta))
    for nn in range(0, len(ta)):
        Eo[nn] = dt * eq_evap(rn[nn], ta[nn], units='mm')  # mm/d
    
    
    """ compute fraction of P as snow """
    prec = dat0[k]['Prec'].values  
    
    # ---state of precipitation [as water (fW) or as snow(fS)]
    Tmax = 1.0; Tmin = 0.0
    fW = np.zeros(len(ta))
    fS = np.zeros(len(ta))
    fW[ta >= Tmax] = 1.0
    fS[ta <= Tmin] = 1.0

    ix = np.where((ta > Tmin) & (ta < Tmax))
    fW[ix] = (ta[ix] - Tmin) / (Tmax - Tmin)
    fS[ix] = 1.0 - fW[ix]
    
    fsnow = fS * prec
    del ix, fW, fS
    
    dat0[k]['ET0'] = Eo
    dat0[k]['fsnow'] = fsnow
    y0 = dat0[k][['Prec', 'Qmeas', 'Qt', 'ET', 'Tr', 'Ef', 'E', 'ET0', 'fsnow']].resample('a', how=sum)

    y1 = dat1[k][['Qt', 'ET']].resample('a', how=sum)
    y2 = dat2[k][['Qt', 'ET']].resample('a', how=sum)
    y0[['Qt_hi', 'ET_lo']] = y1[['Qt', 'ET']]
    y0[['Qt_lo', 'ET_hi']] = y2[['Qt', 'ET']]
    y0['ch'] = ch[k]
    y0['year'] = y0.index.year
        
    An = pd.concat([An, y0])
    del y0, y1, y2, Eo, fsnow

An = An.reset_index(drop=True)

#%% filter bad years with poor Qmeas or very suspicious Q/P -ratios
An = An.dropna(axis=0, how='any')
ix = An.Qmeas.values / An.Prec.values
f = (ix > 0.15) & (ix < 0.85) # removes outliers due to anomalously low Qmeas
An = An.iloc[f,:]
del f

f = An[(An.ch == 2) & (An.year == 2007)].index.values
An = An.drop(f)

f = An[(An.ch == 10) & (An.year <= 2012)].index.values
An = An.drop(f)

f = An[(An.ch == 14) & (An.year == 2006)].index.values
An = An.drop(f)

f = An[(An.ch == 16) & (An.year == 2006)].index.values
An = An.drop(f)

f = An[((An.ch == 19) & (An.year == 2007)) | \
       (An.ch == 19) & (An.year == 2012)].index.values
An = An.drop(f)

f = An[(An.ch == 21) & (An.year == 2011)].index.values
An = An.drop(f)

f = An[(An.ch == 29) & (An.year < 2015)].index.values
An = An.drop(f)

f = An[(An.ch == 32) & (An.year == 2012)].index.values
An = An.drop(f)

f = An[(An.ch == 32) & (An.year == 2012)].index.values
An = An.drop(f)

f = An[((An.ch == 33) & (An.year == 2010)) | \
       (An.ch == 33) & (An.year == 2013)].index.values
An = An.drop(f)       

# bad catchments, measurements are very poor quality
#f = An[(An.ch == 22) |  (An.ch == 30)].index.values
#An = An.drop(f)
del f

# remove same from chdata
f = chdata[(chdata.id == 22) |  (chdata.id == 30)].index.values
chdata = chdata.drop(f)
del f
chdata = chdata.rename(columns={'id': 'ch'})
#chdata.index = chdata['ch']

# compute 'new' variables
dA = 0.1  # catchment area uncertainty
dP = 0.05 # precipitation uncertainty
An['Qm_P'] = An.Qmeas / An.Prec
An['Qm_P_lo'] = An['Qm_P'] * (1. - dA) / (1. + dP)
An['Qm_P_hi'] = An['Qm_P'] * (1. + dA) / (1. - dP)

An['Qt_P'] = An.Qt / An.Prec
An['Qt_P_lo'] = An.Qt_lo / An.Prec
An['Qt_P_hi'] = An.Qt_hi / An.Prec

# same for ET
An['ETm_P'] = 1.0 - An.Qm_P
An['ETm_P_lo'] = 1.0 - An.Qm_P_hi
An['ETm_P_hi'] = 1.0 - An.Qm_P_lo

# modeled ET
An['ET_P'] = 1.0 - An.Qt_P
An['ET_P_lo'] = 1.0 - An.Qt_P_hi
An['ET_P_hi'] = 1.0 - An.Qt_P_lo

# ET0
An['ET0_P'] = An.ET0 / An.Prec
An['alpha_m'] = An.ETm_P / An.ET0_P
An['alpha'] = An.ET_P / An.ET0_P

# merge An and chdata
An = pd.merge(An, chdata, on='ch')

# compute catchment averages
A = An.groupby(['ch']).mean()


##%% plot scatterplot of streamflow & ET
#plt.figure()
#x = An.Qmeas.values
#y = An.Qt.values
#s, s0, r, p, se = stats.linregress(x, y)
#x0 = np.array([min(x), max(x)])
#m = s*x0 + s0
#txt = []
#for k in range(0, len(An)):
#    txt.append( '%s-%d' % (An['ch'].iloc[k], An['year'].iloc[k]))
#plt.subplot(121)
#plt.scatter(x, y, c='r', alpha=0.5)
#for t, i, j in zip(txt, x, y):
#    plt.annotate(t, xy = (i, j), xytext = (0, 0), textcoords='offset points', fontsize=6)
#plt.xlim([100, 800])
#plt.ylim([100, 800])
#plt.plot([100, 800], [100, 800], 'k--')
#plt.plot(x0, m, 'k-')
#plt.title('slope=%.2f +/- %.2f, R2=%.2f' % (s, 2*se, r**2))
#plt.subplot(122)
#
#x1 = An.Qmeas / An.Prec
#y1 = An.Qt / An.Prec
#s, s0, r, p, se = stats.linregress(x1, y1)
#x0 = np.array([min(x1), max(x1)])
#m = s*x0 + s0
#plt.scatter(x1, y1, c='g', alpha=0.5)
#for t, i, j in zip(txt, x1, y1):
#    plt.annotate(t, xy = (i, j), xytext = (0, 0), textcoords='offset points', fontsize=6)
#plt.plot(x0, m, 'k-')
#plt.title('slope=%.2f +/- %.2f, R2=%.2f' % (s, 2*se, r**2))
#plt.xlim([0, 1])
#plt.ylim([0, 1])
#plt.plot([0, 1], [0, 1], 'k--')
#

#%% plot ET / P with errorbars

gp = An.groupby('ch')
x = An.ETm_P.values
y = An.ET_P.values

# estimate slope and r2 for line forced through origin
mm = y[:, np.newaxis]
slp, res, _, _ = np.linalg.lstsq(mm, x)
r2 = 1 - res / sum((y - np.mean(y))**2)

# estimate linear regression with uncertainty
# for computing slopes, set origin to mean(x), mean(y)

xx = x - np.mean(x)
yy = y - np.mean(y)

s, s0, r, p, se = stats.linregress(xx, yy)
r2 = r**2
rmse = np.sqrt(((y - x) ** 2).mean())
me = np.mean(y - x)

x0 = np.array([min(x)-0.05, max(x)+0.05])
xx = np.array([min(xx)-0.05, max(xx)+0.05])
m = s*xx + s0 + np.mean(y)
ml = (s + 2*se)*xx + np.mean(y)
mh = (s - 2*se)*xx+ np.mean(y)

tstr = 's = %.2f$\pm$%.2f\nR$^2$ = %.2f\nRMSE = %.2f\nME = %.2f' % (s, 2*se, r2, rmse, me)

n = len(np.unique(An.ch))

fig, ax = plt.subplots()
plt.plot([0, 1], [0, 1], 'k--', alpha=0.3, linewidth=1)
colors = plt.cm.tab20c(np.linspace(0.01, 0.99, n))
ax.set_prop_cycle('color', colors)
ax.set_aspect('equal')
for name, group in gp:
    yerr = abs(np.array([group.ET_P_hi.values - group.ET_P.values, group.ET_P_lo.values - group.ET_P.values]))
    xerr = abs(np.array([group.ETm_P_hi.values - group.ETm_P.values, group.ETm_P_lo.values - group.ETm_P.values]))
    ax.errorbar(1 - group.Qm_P, 1 - group.Qt_P, fmt='o', xerr=xerr, yerr=yerr,  label=name,
                alpha=0.8, linewidth=1)

ax.legend(numpoints=1, loc='upper right', fontsize=8)

plt.plot(x0, m, 'k-', zorder=100)
plt.plot(x0, ml, 'k--', x0, mh, 'k--', linewidth=1, zorder=50)

plt.axis([0, 1.0, 0, 1.0])
plt.ylabel(r'$ \langle \overline{ET}_{mod} \, / \, \overline{P} \rangle $')
plt.xlabel(r'$ \langle \overline{ET}_{wb} \, / \, \overline{P} \rangle $')

plt.text(0.05, 0.8, tstr, fontsize=10)
plt.show()
#plt.savefig(os.path.join(fig_path, 'Fig_5b_cathcmentET_to_P.png'), dpi=600)
#plt.savefig(os.path.join(fig_path, 'Fig_5b_cathcmentET_to_P.pdf'))


#%% statistics for catchment 'means 

x = A.ETm_P.values
y = A.ET_P.values

# estimate slope and r2 for line forced through origin
mm = y[:, np.newaxis]
slp, res, _, _ = np.linalg.lstsq(mm, x)
r2 = 1 - res / sum((y - np.mean(y))**2)

# estimate linear regression with uncertainty
# for computing slopes, set origin to mean(x), mean(y)

xx = x - np.mean(x)
yy = y - np.mean(y)

s, s0, r, p, se = stats.linregress(xx, yy)
r2 = r**2
rmse = np.sqrt(((y - x) ** 2).mean())
me = np.mean(y - x)

print('s',s, '2xse', 2*se, 'r2', r2, 'rmse', rmse, 'me', me)

#%% plot ET / P comparison to ET0/P
#A =  pd.DataFrame(columns=['year', 'ch', 'Prec', 'Qmeas', 'Qt', 'Qt_hi', 'Qt_lo', 'ET', 'ET_hi', 'ET_lo', 'ET0', 'fsnow'])

""" this is for comparing ET in budyko-framework """
dA = 0.1  # catchment area uncertainty
dP = 0.05 # precipitation uncertainty
An['Qm_P_lo'] = An['Qm_P'] * (1. - dA) / (1. + dP)
An['Qm_P_hi'] = An['Qm_P'] * (1. + dA) / (1. - dP)

An['ETm_P'] = 1.0 - An.Qm_P
An['ETm_P_lo'] = 1.0 - An.Qm_P_hi
An['ETm_P_hi'] = 1.0 - An.Qm_P_lo

# modeled ET
An['ET_P'] = 1.0 - An.Qt_P
An['ET_P_lo'] = 1.0 - An.Qt_P_hi
An['ET_P_hi'] = 1.0 - An.Qt_P_lo

# ET0
An['ET0_P'] = An.ET0 / An.Prec

gp = An.groupby('ch')

n = len(np.unique(An.ch))

fig, ax = plt.subplots(2,1)
fig.set_size_inches(6, 11)
#plt.plot([0, 1], [0, 1], 'k--', alpha=0.3, linewidth=1)
colors = plt.cm.tab20c(np.linspace(0.01, 0.99, n))
ax[0].set_prop_cycle('color', colors)
ax[0].set_aspect('equal')
ax[1].set_prop_cycle('color', colors)
ax[1].set_aspect('equal')
for name, group in gp:
    x = group.ET0_P
    
    #plt.subplot(121)
    y1 = group.ETm_P
    yerr1 = abs(np.array([group.ETm_P_lo.values - group.ETm_P.values, group.ETm_P_hi.values - group.ETm_P.values]))
    # xerr = abs(np.array([group.Qm_P_lo.values - group.Qm_P.values, group.Qm_P_hi.values - group.Qm_P.values]))
    ax[0].errorbar(x, y1, fmt='o', yerr=yerr1,  label=name, alpha=0.7)
    ax[0].plot([0, 1], [0, 1], 'k--', alpha=0.5, linewidth=0.5)

    y2 = group.ET_P
    yerr2 = abs(np.array([group.ET_P_lo.values - group.ET_P.values, group.ET_P_hi.values - group.ET_P.values]))
    ax[1].errorbar(x, y2, fmt='o', yerr=yerr1,  label=name, alpha=0.7)    
    ax[1].plot([0, 1], [0, 1], 'k--', alpha=0.5, linewidth=0.5)

ax[0].legend(numpoints=1, loc='upper left', fontsize=8)

#plt.plot(x0, m, 'k-', zorder=100)

ax[0].axis([0, 1.8, 0, 1.0])
ax[0].set_ylabel(r'$\overline{ET}_{obs} \, / \, \overline{P}$ (-)')
ax[0].set_xlabel(r'$\overline{ET_0} \, / \, \overline{P}$ (-)')

ax[1].axis([0, 1.8, 0, 1.0])
ax[1].set_ylabel(r'$\overline{ET}_{mod} \, / \, \overline{P}$ (-)')
ax[1].set_xlabel(r'$\overline{ET_0} \, / \, \overline{P}$ (-)')

#plt.text(0.05, 0.8, tstr, fontsize=10)
plt.show()
#plt.savefig(os.path.join(fig_path, 'Fig_5b_cathcmentET_to_P.png'), dpi=600)
#plt.savefig(os.path.join(fig_path, 'Fig_5b_cathcmentET_to_P.pdf'))

#%% plot ET / P comparison to ET0/P
# As above but with color = ET0 or LAT and size = LAI 

""" this is for comparing ET in budyko-framework """


gp = An.groupby('LAT_deg')


n = len(np.unique(An.LAT_deg))

fig, ax = plt.subplots(2,1)
fig.set_size_inches(6, 11)
#plt.plot([0, 1], [0, 1], 'k--', alpha=0.3, linewidth=1)
colors = plt.cm.RdBu(np.linspace(0.01, 0.99, n))
ax[0].set_prop_cycle('color', colors)
ax[0].set_aspect('equal')
ax[1].set_prop_cycle('color', colors)
ax[1].set_aspect('equal')
for name, group in gp:
    x = group.ET0_P
    
    #plt.subplot(121)
    y1 = group.ETm_P
    yerr1 = abs(np.array([group.ETm_P_lo.values - group.ETm_P.values, group.ETm_P_hi.values - group.ETm_P.values]))
    # xerr = abs(np.array([group.Qm_P_lo.values - group.Qm_P.values, group.Qm_P_hi.values - group.Qm_P.values]))
    ax[0].errorbar(x, y1, fmt='o', yerr=yerr1,  label=name, alpha=0.7)
    ax[0].plot([0, 1], [0, 1], 'k--', alpha=0.5, linewidth=0.5)

    y2 = group.ET_P
    yerr2 = abs(np.array([group.ET_P_lo.values - group.ET_P.values, group.ET_P_hi.values - group.ET_P.values]))
    ax[1].errorbar(x, y2, fmt='o', yerr=yerr1,  label=name, alpha=0.7)    
    ax[1].plot([0, 1], [0, 1], 'k--', alpha=0.5, linewidth=0.5)

ax[0].legend(numpoints=1, loc='upper left', fontsize=8)

#plt.plot(x0, m, 'k-', zorder=100)

ax[0].axis([0, 1.8, 0, 1.0])
ax[0].set_ylabel(r'$\overline{ET}_{obs} \, / \, \overline{P}$ (-)')
ax[0].set_xlabel(r'$\overline{ET_0} \, / \, \overline{P}$ (-)')

ax[1].axis([0, 1.8, 0, 1.0])
ax[1].set_ylabel(r'$\overline{ET}_{mod} \, / \, \overline{P}$ (-)')
ax[1].set_xlabel(r'$\overline{ET_0} \, / \, \overline{P}$ (-)')

#plt.text(0.05, 0.8, tstr, fontsize=10)
plt.show()
#plt.savefig(os.path.join(fig_path, 'Fig_5b_cathcmentET_to_P.png'), dpi=600)
#plt.savefig(os.path.join(fig_path, 'Fig_5b_cathcmentET_to_P.pdf'))

#%% what explains ET/ETo variability?
B = An.sort_values('LAT_deg', axis=0)

relsize = np.array(B.LAI / max(B.LAI))
marksize = 10 + 150*relsize
rr = B.LAT_deg
yerr = abs(np.array([B.ET_P_hi.values - B.ET_P.values, B.ET_P_lo.values - B.ET_P.values]))
xerr = abs(np.array([B.ETm_P_hi.values - B.ETm_P.values, B.ETm_P_lo.values - B.ETm_P.values]))

fig1, ax1 = plt.subplots(1,1)
fig1.set_size_inches(6, 6)

divider = make_axes_locatable(ax1)
cax = divider.append_axes('right', size='5%', pad=0.05)
ax1.set_aspect('equal')
# ax1[1].set_aspect('equal')

# P/E scatterplot
x = B.ETm_P
y = B.ET_P

xx = x - np.mean(x)
yy = y - np.mean(y)

s, s0, r, p, se = stats.linregress(xx, yy)
r2 = r**2
rmse = np.sqrt(((y - x) ** 2).mean())
me = np.mean(y - x)

x0 = np.array([min(x)-0.05, max(x)+0.05])
xx = np.array([min(xx)-0.05, max(xx)+0.05])
m = s*xx + s0 + np.mean(y)
ml = (s + 2*se)*xx + np.mean(y)
mh = (s - 2*se)*xx+ np.mean(y)

tstr = 's = %.2f$\pm$%.2f\nR$^2$ = %.2f\nRMSE = %.2f\nME = %.2f' % (s, 2*se, r2, rmse, me)

#for j in range(len(x)):
#    ax1[0].errorbar(x[j], y[j], yerr[j], xerr[j], marker='None', mec='k', mfc='k', alpha=0.9, zorder=-10)
gp = B.groupby('LAT_deg')

n = len(np.unique(B.LAT_deg))
colors = plt.cm.RdBu(np.linspace(0.01, 0.99, n))
#ax1[0].set_prop_cycle('color', colors)
ax1.set_aspect('equal')
j = 0
for name, group in gp:
    x1 = group.ETm_P
    y1 = group.ET_P

    yerr1 = abs(np.array([group.ET_P_hi.values - group.ET_P.values, group.ET_P_lo.values - group.ET_P.values]))
    xerr1 = abs(np.array([group.ETm_P_hi.values - group.ETm_P.values, group.ETm_P_lo.values - group.ETm_P.values]))
    
#    xerr1 = abs(np.array([group.ETm_P_lo.values - group.ETm_P.values, group.ETm_P_hi.values - group.ETm_P.values]))
#    yerr1 = abs(np.array([group.ET_P_lo.values - group.ET_P.values, group.ET_P_hi.values - group.ET_P.values]))
    ax1.errorbar(x1, y1, yerr=yerr1,  xerr=xerr1, fmt='None', ecolor=colors[j], alpha=0.6, zorder=-10, linewidth=1)
    j += 1

#ax1[0].errorbar(x, y, yerr=yerr, xerr=xerr, fmt='None', ecolor='k', alpha=0.3, zorder=-10)

ax1.plot([0, 1], [0, 1], 'k--', alpha=0.4, linewidth=1)
sc = ax1.scatter(x, y, c=rr, edgecolor = 'k', s=marksize, alpha=0.6, cmap='RdBu')
ax1.text(0.05, 0.75, tstr, fontsize=10)
cb = fig.colorbar(sc, cax=cax, orientation='vertical')
ax1.plot(x0, m, 'k-', zorder=100)
ax1.plot(x0, ml, 'k--', x0, mh, 'k--', linewidth=1, zorder=50)
cb.set_label('LAT (deg)', rotation=90, fontsize=9)

ax1.axis([0, 1.0, 0, 1.0])
ax1.set_xlabel(r'$ \langle \overline{ET}_{wb} \, / \, \overline{P} \rangle $')
ax1.set_ylabel(r'$ \langle \overline{ET}_{mod} \, / \, \overline{P} \rangle $')

plt.show()

plt.savefig('cathcmentET.png', dpi=600)
plt.savefig('cathcmentET.pdf')


#%% E/P vs Eo/P scatterplot
# E/P vs Eo/P scatterplot Modeled values

fig2, ax2 = plt.subplots(1,1)
fig2.set_size_inches(6, 6)

divider = make_axes_locatable(ax2)
cax2 = divider.append_axes('right', size='5%', pad=0.05)
ax2.set_aspect('equal')

x = B.ET0_P
y = B.ETm_P
sc2=ax2.scatter(x, y, c=rr, edgecolor = 'k', s=marksize, alpha=0.6, cmap='RdBu')

cb2 = fig.colorbar(sc2, cax=cax2, orientation='vertical')
cb2.set_label('LAT (deg)', rotation=90, fontsize=9)
ax2.set_xlim([0, 1.3])
ax2.set_ylim([0, 1.0])
ax2.plot([0, 1], [0, 1], 'k--', alpha=0.4, linewidth=1)
ax2.set_ylabel(r'$\overline{ET_{wb}} \, / \, \overline{P}$ (-)')
ax2.set_xlabel(r'$\overline{ET_0} \, / \, \overline{P}$ (-)')


plt.show()

plt.savefig(os.path.join(fig_path, 'Fig6_cathcmentET_Budyko.png'), dpi=600)
plt.savefig(os.path.join(fig_path, 'Fig6_cathcmentET_Budyko.pdf'))

#%% regression model to test ET relations to climatic and catchment variables
del y, x

#y = An.alpha_m
y = An.ET_P * An.Prec.values
y.index = An.ch
x = An[['ET0', 'Prec', 'fsnow', 'LAI', 'f_decid', 'peat']].copy()
x['LAI'] = np.sqrt(x['LAI'])
x.index = An.ch

model = smf.MixedLM(y, x, groups=y.index)
result = model.fit()
print(result.summary())


# simple multiple regression
y = An.alpha_m
x = An[['ET0', 'Prec', 'fsnow', 'LAI', 'f_decid', 'peat']].copy()
x['LAI'] = np.sqrt(x['LAI'])
x = sm.add_constant(x)
est = sm.OLS(y, x).fit()
est.summary()


xx = A[['LAT_deg', 'ET0', 'Prec', 'fsnow', 'LAI', 'f_decid', 'peat', 'Tr', 'Ef', 'E']].copy()
# xx['LAI'] = np.sqrt(xx['LAI'])
print('means', xx.corr())

xxx = An[['LAT_deg', 'ET0', 'Prec', 'fsnow', 'LAI', 'f_decid', 'peat', 'Tr', 'Ef', 'E']].copy()
# xxx['LAI'] = np.sqrt(xxx['LAI'])
print('raw', xxx.corr())

#xm = sm.add_constant(xm)
#estm = sm.OLS(ym, xm).fit()
#estm.summary()
#A =  pd.DataFrame(columns=['year', 'ch', 'Prec', 'Qmeas', 'Qt', 'Qt_hi', 'Qt_lo', 'ET', 'ET_hi', 'ET_lo', 'ET0', 'fsnow'])

#%% multiple regression for model-data residuals

# simple multiple regression
y = An.ET_P - An.ETm_P
x = An[['LAI', 'TWI', 'TWI_b', 'f_decid', 'fine', 'LAT_deg',
       'LON_deg', 'med', 'coarse', 'peat', 'top_m', 'area']].copy()
# x['LAI'] = np.sqrt(x['LAI'])
x = sm.add_constant(x)
est = sm.OLS(y, x).fit()
print(est.summary())

model = smf.MixedLM(y, x, groups=y.index)
result = model.fit()
print(result.summary())

