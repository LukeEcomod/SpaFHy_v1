# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 13:21:26 2017

@author: slauniai

MAKES FIGURES 6-8 IN THE MANUSCRIPT
CHANGE PATHS TO READ .nc FILE AND .pk -FILE.

"""
import sys
sys.path.append(r'\repositories\SpaFHy_v1')

import os
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.colors as mplcolors
import matplotlib.cm as mplcm

import pickle
from netCDF4 import Dataset #, date2num

# import model 
# import spafhy

eps = np.finfo(float).eps

# change working dir
os.chdir(r'c:\repositories\SpaFHy_v1\FigsC3')

# results file
ncf_file = r'c:\repositories\SpaFHy_v1\results\C3.nc'

# pickled model object
pk_file = r'c:\repositories\SpaFHy_v1\results\C3model.pk'

""" load pre-computed results for C3 """
# spa instance
with open(pk_file, 'rb') as ff:
    spa, Qmeas, FORC = pickle.load(ff)

# get time-index when results start
ix = 1 + np.where(FORC.index == spa.pgen['spinup_end'])[0][0]
tvec = FORC.index # date vector


gis = spa.GisData
twi = gis['twi']
LAIc = gis['LAI_conif']
LAId = gis['LAI_decid']
soil = gis['soilclass']

# soil type indexes
peat_ix = np.where(soil == 4)
med_ix = np.where(soil == 2)
coarse_ix = np.where(soil == 1)

# indices for high and low twi
htwi_ix = np.where(twi > 12)
ltwi_ix = np.where(twi < 7)

# open link to results in netcdf:
dat = Dataset(ncf_file, 'r')

# results per sub-model
cres = dat['cpy']   # canopy -submodel
bres = dat['bu']    # bucket -submodel
tres = dat['top']   # topmodel - submodel


#%% FIG. 9 long-term ET/P and components

# NOTE [ix:] ensures we neglect spinup period!

P = np.sum(FORC['Prec'][ix:])*spa.dt  # total precip

# ET components, from spinup end to end of simulation

TR = np.array(cres['Transpi'][ix:, :, :])  # transpi
EF = np.array(cres['Efloor'][ix:, :, :])  # floor evap
IE = np.array(cres['Evap'][ix:, :, :])# interception evap
ET = TR + EF + IE

#--------

fig = plt.figure()
fig.set_size_inches(8.0, 10.0)

# maps to lhs of the figure

plt.subplot(521)
sns.heatmap(LAIc + LAId, cmap='coolwarm',cbar=True, xticklabels=False, yticklabels=False); 
plt.title('LAI (m$^2$m$^{-2}$)')

plt.subplot(523)
sns.heatmap(np.sum(ET, axis=0)/P, cmap='coolwarm',cbar=True, xticklabels=False, yticklabels=False); 
plt.title('$\overline{ET}$/$\overline{P}$ (-)')

plt.subplot(525)
sns.heatmap(np.sum(IE, axis=0)/P, cmap='coolwarm',cbar=True, xticklabels=False, yticklabels=False); 
plt.title('$\overline{E}$/$\overline{P}$ (-)')

plt.subplot(527)
sns.heatmap(np.sum(TR, axis=0)/P, cmap='coolwarm',cbar=True, xticklabels=False, yticklabels=False); 
plt.title('$\overline{T_r}$/$\overline{P}$ (-)')

plt.subplot(529)
sns.heatmap(np.sum(EF, axis=0)/P, cmap='coolwarm',cbar=True, xticklabels=False, yticklabels=False); 
plt.title('$\overline{E_f}$/$\overline{P}$ (-)')

# add LAI -relations to rhs of figure

x = LAId + LAIc

plt.subplot(5,2,4)
y = np.sum(ET, axis=0) / P
rr = soil.copy()
y = np.ravel(y)
rr[rr == 4] = 3 # peat
cm = plt.get_cmap('coolwarm_r', 3)

norm = mplcolors.BoundaryNorm(np.arange(1, 5) - 0.5, 4)
plt.scatter(x, y, c=rr, cmap=cm, norm=norm, s=10, alpha=0.5)
cb1= plt.colorbar(ticks=np.arange(0, 4))
cb1.ax.set_yticklabels(['coarse','med','peat']) 

ax = plt.gca()
ax.grid(linestyle='--', alpha=0.5)
plt.xlim([0, 9])
plt.ylim([0.2, 0.7])
#plt.ylabel('$\overline{ET}$/$\overline{P}$ (-)', labelpad=-1);

# interception
plt.subplot(5,2,6)
rr = LAId / x
y = np.sum(IE, axis=0) / P
sc = plt.scatter(x, y, c=rr, vmin=0, vmax=0.6, s=10, alpha=0.5, cmap='coolwarm')
cb = plt.colorbar(sc, ticks = [0.0, 0.20, 0.40, 0.60])
cb.set_label('LAI$_d$ / LAI', rotation=90, fontsize=9)
plt.xlim([0, 9])
plt.ylim([0, 0.3])
#plt.ylabel('$\overline{E}$/$\overline{P}$ (-)', labelpad=-1);
ax = plt.gca()
ax.grid(linestyle='--', alpha=0.5)


# transpiration / P
plt.subplot(528)
y = np.sum(TR, axis=0) / P
rr = LAId / x # deciduous fraction
sc = plt.scatter(x, y, c=rr, vmin=0, vmax=0.6, s=10, alpha=0.5, cmap='coolwarm')
cb = plt.colorbar(sc, ticks = [0.0, 0.20, 0.40, 0.60])
cb.set_label('LAI$_d$ / LAI', rotation=90, fontsize=9)
plt.xlim([0, 9])
ax = plt.gca()
ax.grid(linestyle='--', alpha=0.5)
plt.xlim([0, 9])
plt.ylim([0, 0.4])
#plt.ylabel('$\overline{T_r}$/$\overline{P}$ (-)', labelpad=-2);

# floor evaporation / P
plt.subplot(5,2,10)
rr = twi
y = np.sum(EF, axis=0) / P

sc = plt.scatter(x, y, c=rr, s=10, alpha=0.5, cmap='coolwarm_r')
cb = plt.colorbar(sc)
cb.set_label('TWI', rotation=90, fontsize=9)
plt.xlim([0, 9])
plt.ylim([0, 0.3])
#plt.ylabel('$\overline{E_f}$/$\overline{P}$ (-)', labelpad=-2)
plt.xlabel('LAI (m$^2$m$^{-2}$)')
ax = plt.gca()
ax.grid(linestyle='--', alpha=0.5)
plt.xlim([0, 9])
plt.ylim([0, 0.3])

plt.savefig('ch3_ETvariability.png', dpi=660)
# plt.savefig('ch3_ETvariability.pdf')

#%% plot ET partitioning as function of LAI

t_et = np.ravel(np.sum(ET, axis=0))
t_tr = np.ravel(np.sum(TR, axis=0))
t_ef = np.ravel(np.sum(EF, axis=0))
t_e = np.ravel(np.sum(IE, axis=0))

fig = plt.figure()
fig.set_size_inches(4.5, 3.2)

x = np.ravel(LAId + LAIc)
cm = plt.get_cmap('Blues')

plt.plot(x, t_e/t_et, 'o', color=cm(0.4), markersize=4, alpha=0.4, label='E')
plt.plot(x, t_tr/t_et, 'o', color=cm(0.7), markersize=4, alpha=0.7, label='T$_r$')
plt.plot(x, t_ef/t_et, 'o', color=cm(1.0), markersize=4, alpha=0.8, label='E$_f$')
plt.legend(fontsize=10, frameon=False)
plt.ylabel('(-)', fontsize=10)
plt.xlabel('LAI (m$^2$m$^{-2}$)', fontsize=10, labelpad=-3)
ax = plt.gca()
ax.grid(linestyle='--', alpha=0.5)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.xlim([0, 9])
plt.ylim([0, 1.0])


plt.savefig(r'ch3_ETpartitioning.png', dpi=660)
# plt.savefig(r'ch3_ETpartitioning.pdf')


#%% FIG. 10: maximum snow water equivalent SWE

# snow water equivalent, seek maximum timing
SWE = np.array(cres['SWE'][ix:, :, :]) # SWE
a = np.nansum(SWE, axis=1)
a = np.nansum(a, axis=1)
swe_max_ix = int(np.where(a == np.nanmax(a))[0][0])
del a

#----
fig = plt.figure()
fig.set_size_inches(8.0, 2.5)

plt.subplot(121)
y = SWE[swe_max_ix,:,:]
sns.heatmap(y, cmap='coolwarm',cbar=True, xticklabels=False, yticklabels=False)
plt.title('max SWE (mm)', fontsize=10)

laiw = LAId * 0.1 + LAIc  # wintertime lai
f = y > 0
yy = y[f]
laiw = laiw[f]
plt.subplot(122)
cm = plt.get_cmap('coolwarm')
plt.plot(laiw, yy/max(yy), 'o', color=cm(0.1), markersize=5, alpha=0.2)
plt.title('relative SWE (-)', fontsize=10)
plt.xlabel('winter LAI  (m$^2$m$^{-2}$)', fontsize=10, labelpad=-30)
plt.xlim([0, 7.5])
plt.ylim([0.75, 1.02])

ax = plt.gca()
ax.yaxis.set_label_position("right")
#ax.yaxis.tick_right()
ax.set_xticks([0,1,2,3,4,5,6,7])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
#ax.set_alpha(0.5)
ax.grid(linestyle='--', alpha=0.5)

plt.savefig(r'ch3_swe_variability.png', dpi=660)
plt.savefig(r'ch3_swe_variability.pdf')


#%% Fig 6: streamflow
f0 = int(np.where(tvec == '2012-01-01')[0])
f1 = int(np.where(tvec == '2014-01-01')[0])

tvec0 = tvec[f0:f1]
Prec = FORC['Prec'].iloc[f0:f1] # precipitation
Wliq = np.array(bres['Wliq'])[f0:f1,:,:] # root zone moisture m3m-3

f, g = np.where(spa.GisData['cmask'] == 1)
swe = np.array(cres['SWE'])
swe = swe[f0:f1, f,g]
a = np.nansum(swe, axis=1)

swe = a / len(f)
max_swe = max(swe)

Qt = 1e3*np.array(tres['Qt'])[f0:f1] # modeled streamflow

# rootzone water storage; seek catchment maximum and minimum timing
a = np.nansum(Wliq, axis=1)
a = np.nansum(a, axis=1)

ix_slow = int(np.where(a == np.nanmin(a))[0])
ix_shi = int(np.where(a == np.nanmax(a))[0])

Qm = Qmeas[(Qmeas.index >='2012-01-01') & (Qmeas.index <'2014-01-01') ]

fig1, ax = plt.subplots()
fig1.set_size_inches(6.5, 4.5)


ax.plot(tvec0, Qm, 'k.', tvec0, Qt, 'r-', linewidth=1)
ax.plot(tvec0[ix_slow], -0.5, marker='o', mec='k', mfc='g', alpha=0.7, ms=8.0)
ax.plot(tvec0[ix_shi], -0.5, marker='o', mec='k', mfc='b', alpha=0.7, ms=8.0)
ax.set_xlim(['2012-01-01', '2014-01-01'])
ax.set_ylabel(r'$\langle Q_f \rangle$ (mm d$^{-1}$)')
ax2=ax.twinx()   
ax2.bar(tvec0, 3600*24.0*Prec, color='k', width=1)
ax2.plot(tvec0, swe / 10.0, 'b-', linewidth=1.5)
# ax2.plot(tvec0, 20*swe, 'r--', linewidth=1)
ax2.set_ylabel(r'P (mm d$^{-1}$) & 0.1 x $\langle SWE \rangle$ (mm)'); plt.ylim([0, 50])
ax2.invert_yaxis()
#plt.xticks(rotation=30)        
for tick in ax.get_xticklabels():
        tick.set_rotation(30)

plt.savefig(r'ch3_discharge_timeseries.png', dpi=660)
#plt.savefig(r'ch3_discharge_timeseries.pdf')

#%% Fig 8: snapshots of soil moisture 

# select hydrologically contrasting years 2012 and 2013 for example
f0 = int(np.where(tvec == '2012-01-01')[0])
f1 = int(np.where(tvec == '2014-01-01')[0])

tvec0 = tvec[f0:f1]
Qt = 1e3*np.array(tres['Qt'])[f0:f1] # modeled streamflow mm/d
Wliq = np.array(bres['Wliq'])[f0:f1,:,:] # root zone moisture m3m-3
S = np.array(tres['S'])[f0:f1]  # saturation deficit 
del f0, f1 


# local saturation deficits
s_hi = 1e3*spa.top.local_s(S[ix_shi])
s_hi[s_hi<0] = 0.0
s_low = 1e3*spa.top.local_s(S[ix_slow])
s_low[s_low<0] = 0.0

# convert 1-D array back to 2-D grid
s_hi = spa._to_grid(s_hi)
s_low = spa._to_grid(s_low)

fig = plt.figure()
fig.set_size_inches(8.0, 4.0)

plt.subplot(221)
sns.heatmap(Wliq[ix_slow,:,:], cmap='RdBu',cbar=True, vmin=0.1, vmax=0.9, 
            xticklabels=False, yticklabels=False);
tt = tvec0[ix_slow]
tt = tt.strftime('%Y-%m-%d')
plt.title('dry: ' + tt, fontsize=10)
plt.xlabel('$\\theta$ (m$^3$m$^{-3}$)', fontsize=10)

plt.subplot(222)
sns.heatmap(Wliq[ix_shi,:,:], cmap='RdBu',cbar=True, vmin=0.1, vmax=0.9, xticklabels=False, yticklabels=False);
tt = tvec0[ix_shi]
tt = tt.strftime('%Y-%m-%d')
plt.title('wet: ' + tt, fontsize=10)
plt.xlabel('$\\theta$ (m$^3$m$^{-3}$)', fontsize=10)

plt.subplot(223)
sns.heatmap(s_low, cmap='RdBu_r',cbar=True, vmin=0.0, vmax=180, xticklabels=False, yticklabels=False);
plt.xlabel('S (mm)', fontsize=10)

plt.subplot(224)
sns.heatmap(s_hi, cmap='RdBu_r',cbar=True, vmin=0.0, vmax=180, xticklabels=False, yticklabels=False);
plt.xlabel('S (mm)', fontsize=10)

plt.savefig('ch3_moisture.png', dpi=660)
#plt.savefig('ch3_moisture.pdf')

#%% Extract data for analysis of soil moisture variability
f0 = int(np.where(tvec == '2012-01-01')[0])
f1 = int(np.where(tvec == '2014-01-01')[0])

tvec0 = tvec[f0:f1]
mo = np.array(tvec0.month)
yr = np.array(tvec0.year)
yr[yr == 2012] = 0
yr[yr == 2013] = 1
doy0 = np.array(tvec0.dayofyear)

f, g = np.where(spa.GisData['cmask'] == 1)

Wliq = np.array(bres['Wliq'])[f0:f1,f,g]
twi0 = np.array(spa.GisData['twi'][f,g])
soil0 = np.array(spa.GisData['soilclass'][f,g])
LAI = LAIc + LAId
lai0 = LAI[f,g]

#%% soil moisture variability
from soil_moisture_budget import time_stability


# relative difference, MRD, std_MRD, rank-change index, mean theta, sigma theta, cv theta
delta, delta_ave, delta_std, rci, wm, sigma, cv = time_stability(Wliq[yr == 0])      
rci = rci / max(rci)  # normalize rci peak to 1
wmd = np.percentile(Wliq[yr == 0], [10.0, 50.0, 90.0], axis=1)
# range of W

delta1, delta_ave1, delta_std1, rci1, wm1, sigma1, cv1 = time_stability(Wliq[yr == 1])      
rci1 = rci1 / max(rci1)  # normalize rci peak to 1
wmd1 = np.percentile(Wliq[yr == 1], [10.0, 50.0, 90.0], axis=1)

fig1, ax = plt.subplots(2,2)
fig1.set_size_inches(8, 8)

# soil moisture timeseries
ax[0,0].fill_between(doy0[yr == 0], wmd[0,:], wmd[2,:], color='k', alpha=0.2)# zorder=-10)
ax[0,0].plot(doy0[yr==0], wm, 'k-', linewidth=1)

ax[0,1].fill_between(doy0[yr == 1], wmd1[0,:], wmd1[2,:], color='k', alpha=0.2)
ax[0,1].plot(doy0[yr==1], wm1, 'k-', linewidth=1)


ax[0,0].set_ylabel(r'$\langle \theta \rangle $', fontsize=12)
ax[0,0].set_title('2012: wet')
ax[0,1].set_yticklabels([])
ax[0,1].set_title('2013: dry')

for tick in ax[0,0].get_xticklabels():
        tick.set_rotation(30)
ax[0,0].set_xticks(np.arange(30, 366, step=45))
ax[0,0].set_ylim([0.1, 0.5]); ax[0,0].set_xlim([0, 366])
ax[0,0].set_xlabel('doy')

for tick in ax[0,1].get_xticklabels():
        tick.set_rotation(30)
ax[0,1].set_xticks(np.arange(30, 366, step=45))
ax[0,1].set_ylim([0.1, 0.5]); ax[0,1].set_xlim([0, 366])
ax[0,1].set_xlabel('doy')

# get color

cmap = mplcm.get_cmap('nipy_spectral')
cc  = cmap(0.15)
ax3=ax[0,0].twinx()
ax3.plot(doy0[yr==0], sigma, color=cc, linestyle='--', linewidth=1)
#ax3.set_ylabel(r'$\sigma_{\theta}$ (-)')# plt.ylim([0, 50])
ax3.set_ylim([0.05, 0.14])
#ax3.yaxis.set_ticks_position('right')
ax3.set_yticklabels([])
ax4=ax[0,1].twinx()
ax4.plot(doy0[yr==1], sigma1, color=cc, linestyle='--', linewidth=1)
ax4.set_ylabel(r'$\sigma_{\theta}$ (-)', fontsize=12) # plt.ylim([0, 50])
ax4.set_ylim([0.05, 0.14])

del cc

ax[1,0].plot(wm, sigma, 'k-', alpha=0.5)
sc = ax[1,0].scatter(wm, sigma, c=doy0[yr == 0], marker='o', edgecolor='k', alpha=0.8, s=40, vmin=1, vmax=365, cmap='nipy_spectral')
ax[1,0].set_xlabel(r'$\langle \theta \rangle$', fontsize=12); 
ax[1,0].set_ylabel(r'$\sigma_{\theta}$', fontsize=12);
#plt.scatter(wm, cv, c=doy0[yr == 0], marker='o', edgecolor='k', alpha=0.5, s=40, vmin=1, vmax=365, cmap='gist_ncar')
#plt.xlabel(r'$\langle \theta \rangle$'); plt.ylabel(r'$CV_{\theta}$'); plt.title('2012: wet')
ax[1,0].text(0.12, 0.13, 'a)')
ax[1,0].arrow(0.32, 0.12, -0.03, -0.01, color='k', alpha=0.5, head_width=0.005)
ax[1,0].set_ylim([0.05, 0.14]);
ax[1,0].set_xlim([0.1, 0.5])

# plot colorbar as inset
cbaxes = inset_axes(ax[1,0], width="10%", height="50%", loc=2) 
cbar=fig.colorbar(sc, cax=cbaxes, orientation='vertical')
cbar.ax.tick_params(labelsize=8)
cbar.set_label('doy', rotation=90, fontsize=9)

cbaxes.yaxis.set_ticks_position('right')
#cb1 = plt.colorbar(sc); #cb1.ax.set_ylabel('doy')


ax[1,1].plot(wm1, sigma1, 'k-', alpha=0.5) 
sc1 = ax[1,1].scatter(wm1, sigma1, c=doy0[yr == 1], marker='o', edgecolor='k', alpha=0.8, s=40, vmin=1, vmax=365, cmap='nipy_spectral'); #map='gist_ncar')
ax[1,1].set_xlabel(r'$\langle \theta \rangle$', fontsize=12);
ax[1,1].yaxis.set_label_position("right")
#plt.scatter(wm1, cv1, c=doy0[yr == 1], marker='o', edgecolor='k', alpha=0.5, s=40, vmin=1, vmax=365, cmap='gist_ncar')
#plt.xlabel(r'$\langle \theta \rangle$'); plt.ylabel(r'$CV_{\theta}$'); plt.title('2013: dry')
#plt.text(0.12, 0.65, 'c)')
#plt.arrow(0.35, 0.35, -0.08, +0.1, color='k', alpha=0.5, head_width=0.01)
#plt.arrow(0.18, 0.35, +0.05, -0.1, color='k', alpha=0.5, head_width=0.01)
#plt.ylim([0.12, 0.7]); plt.xlim([0.1, 0.5])
ax[1,1].text(0.12, 0.13, 'b)')
ax[1,1].arrow(0.30, 0.11, -0.04, +0.0, color='k', alpha=0.5, head_width=0.005)
#plt.arrow(0.18, 0.35, +0.05, -0.1, color='k', alpha=0.5, head_width=0.01)
ax[1,1].set_ylim([0.05, 0.14]); 
ax[1,1].yaxis.set_ticks_position('right')
ax[1,1].set_ylabel(r'$\sigma_{\theta}$', fontsize=12)
ax[1,1].set_xlim([0.1, 0.5])

plt.savefig('ch3_moisture_statistics.png', dpi=660)
plt.savefig('ch3_moisture_statistics.pdf')

##%% soil moisture variability - more figures
#
##from soil_moisture_budget import time_stability
#
#
## relative difference, MRD, std_MRD, rank-change index, mean theta, sigma theta, cv theta
#delta, delta_ave, delta_std, rci, wm, sigma, cv = time_stability(Wliq[yr == 0])      
#rci = rci / max(rci)  # normalize rci peak to 1
#
#delta1, delta_ave1, delta_std1, rci1, wm1, sigma1, cv1 = time_stability(Wliq[yr == 1])      
#rci1 = rci1 / max(rci1)  # normalize rci peak to 1
#
#fig1, ax1 = plt.subplots(2,2)
#fig1.set_size_inches(8, 8)
#
#plt.subplot(221);
#plt.plot(wm, sigma, 'k-', alpha=0.5)
#plt.scatter(wm, sigma, c=doy0[yr == 0], marker='o', edgecolor='k', alpha=0.5, s=40, vmin=1, vmax=365, cmap='gist_ncar')
#plt.xlabel(r'$\langle \theta \rangle$'); plt.ylabel(r'$\sigma_{\theta}$'); plt.title('2012: wet')
##plt.scatter(wm, cv, c=doy0[yr == 0], marker='o', edgecolor='k', alpha=0.5, s=40, vmin=1, vmax=365, cmap='gist_ncar')
##plt.xlabel(r'$\langle \theta \rangle$'); plt.ylabel(r'$CV_{\theta}$'); plt.title('2012: wet')
#plt.text(0.12, 0.13, 'a)')
#plt.arrow(0.32, 0.12, -0.03, -0.01, color='k', alpha=0.5, head_width=0.005)
#plt.ylim([0.05, 0.14]); 
#plt.xlim([0.11, 0.5])
#cb1 = plt.colorbar(); #cb1.ax.set_ylabel('doy')
#
#plt.subplot(222);
#plt.plot(wm1, sigma1, 'k-', alpha=0.5) 
#plt.scatter(wm1, sigma1, c=doy0[yr == 1], marker='o', edgecolor='k', alpha=0.5, s=40, vmin=1, vmax=365, cmap='gist_ncar')
#plt.xlabel(r'$\langle \theta \rangle$'); plt.ylabel(r'$\sigma_{\theta}$'); plt.title('2012: wet')
#
##plt.scatter(wm1, cv1, c=doy0[yr == 1], marker='o', edgecolor='k', alpha=0.5, s=40, vmin=1, vmax=365, cmap='gist_ncar')
##plt.xlabel(r'$\langle \theta \rangle$'); plt.ylabel(r'$CV_{\theta}$'); plt.title('2013: dry')
##plt.text(0.12, 0.65, 'c)')
##plt.arrow(0.35, 0.35, -0.08, +0.1, color='k', alpha=0.5, head_width=0.01)
##plt.arrow(0.18, 0.35, +0.05, -0.1, color='k', alpha=0.5, head_width=0.01)
##plt.ylim([0.12, 0.7]); plt.xlim([0.1, 0.5])
#plt.text(0.12, 0.13, 'b)')
#plt.arrow(0.30, 0.11, -0.04, +0.0, color='k', alpha=0.5, head_width=0.005)
##plt.arrow(0.18, 0.35, +0.05, -0.1, color='k', alpha=0.5, head_width=0.01)
#plt.ylim([0.05, 0.14]); 
#plt.xlim([0.1, 0.5])
#cb2 = plt.colorbar(); cb2.ax.set_ylabel('doy')
##plt.subplot(223);
##plt.plot(wm, sigma, 'k-', alpha=0.5)
##plt.scatter(wm,  sigma, c=doy0[yr == 0], marker='o', edgecolor='k', s=40)
##plt.xlabel(r'$\langle \theta \rangle$'); plt.ylabel(r'$\sigma_{theta}$'); plt.title('2012: wet')
##plt.colorbar()
##
##plt.subplot(224);
##plt.plot(wm1, sigma1, 'k-', alpha=0.5) 
##plt.scatter(wm1, sigma1, c=doy0[yr == 1], marker='o', edgecolor='k', s=40)
##plt.xlabel(r'$\langle \theta \rangle$'); plt.ylabel(r'$\sigma_{theta}$'); plt.title('2013: dry')
##plt.colorbar()
#
##plt.figure()
#
#mark = ['s','o','h']
#n = 0
#for k in np.unique(soil0):
#    f = np.where(soil0 == k)[0]
#
#    plt.subplot(223)
#    plt.scatter(delta_ave[f], rci[f], c=twi0[f], marker=mark[n], edgecolor='k', s=40, alpha=0.5, cmap='coolwarm_r')
#
#    plt.subplot(224)
#    plt.scatter(delta_ave1[f], rci1[f], c=twi0[f], marker=mark[n], edgecolor='k', s=40, alpha=0.5, cmap='coolwarm_r')
#
#    n += 1
#    
#plt.subplot(223)
#plt.xlabel(r'MRD'); plt.ylabel('RCI'); plt.xlim([-0.2, 0.6]); plt.ylim([0, 1])
#plt.text(-0.15, 0.9, 'b)')
#cb3 = plt.colorbar(); #cb3.ax.set_ylabel('TWI')
#
#plt.subplot(224); 
#plt.text(-0.15, 0.9, 'd)')
#plt.xlabel(r'MRD'); plt.ylabel('RCI'); plt.xlim([-0.2, 0.6]); plt.ylim([0, 1])
#cb4 = plt.colorbar(); cb4.ax.set_ylabel('TWI')
#
##plt.savefig(os.path.join(r'c:\ModelResults\Spathy\Figs','ch3_moisture_statistics.png'), dpi=660)
##plt.savefig(os.path.join(r'c:\ModelResults\Spathy\Figs','ch3_moisture_statistics.pdf'))
#
#plt.figure()
#
#mark = ['s','o','h']
#n = 0
#for k in np.unique(soil0):
#    f = np.where(soil0 == k)[0]
#
#    plt.subplot(223)
#    plt.scatter(twi0[f], rci[f], c=lai0[f], marker=mark[n], edgecolor='k', s=40, alpha=0.5, cmap='coolwarm_r')
#
#    plt.subplot(224)
#    plt.scatter(twi0[f], rci1[f], c=lai0[f], marker=mark[n], edgecolor='k', s=40, alpha=0.5, cmap='coolwarm_r')
#
#    n += 1
#
#plt.subplot(223); 
#plt.xlabel(r'TWI'); plt.ylabel('RCI')
#cb4 = plt.colorbar(); cb4.ax.set_ylabel('LAI')
#plt.subplot(224); 
#plt.xlabel(r'TWI'); plt.ylabel('RCI')
#cb4 = plt.colorbar(); cb4.ax.set_ylabel('LAI')

