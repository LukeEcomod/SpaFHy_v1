# -*- coding: utf-8 -*-
"""
Created on Sun Dec  2 18:55:46 2018

@author: slauniai

Creates supplementary figure S2 to Launiainen et al. 2019 GMD 
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle

with open('ch3_peruskartta.pk', 'rb') as f:
    basemap = pickle.load(f) # (raster, cmap)

with open('C3_gis.pk', 'rb') as f:
    gis = pickle.load(f)

print(gis.keys())

LAIc = gis['LAI_conif']
LAId = gis['LAI_decid']
soil = gis['soilclass']
cmask = gis['cmask']
dem = gis['dem']
cf = gis['cf']
hc = gis['hc']
twi = gis['twi']


#%%
fig = plt.figure()
fig.set_size_inches(8.0,8.0)

plt.subplot(4,2,1)

#xcolormap = RdYlGn
#note: bacgroundmap must be cropped to same coordinates than x but resoltion
# can be different
r, c = np.shape(basemap[0])
extent1 = (-0.5, c - 0.5, r - 0.5, -0.5) # extent into where x is re-scaled

#if not fig_nr:
#    h = plt.figure()
#else:
#    h = plt.figure(fig_nr)
#
#if not vmin:
#    vmin = np.nanmin(x)
#if not vmax:
#    vmax = np.nanmax(x)

# peruskartta        
plt.imshow(basemap[0], cmap=basemap[1], alpha=0.8)

x = cmask

vmin = np.nanmin(x)
vmax = np.nanmax(x)
plt.imshow(x, extent=extent1, cmap= plt.get_cmap('coolwarm', 1), alpha=0.1)
cb1= plt.colorbar(ticks=[0])
cb1.ax.set_ylabel('C3')

##
plt.subplot(4,2,3)

# peruskartta        
plt.imshow(basemap[0], cmap=basemap[1], alpha=0.8)

x = LAIc + LAId

vmin = np.nanmin(x)
vmax = np.nanmax(x)
plt.imshow(x, extent=extent1, cmap='coolwarm', alpha=0.6)
cb1= plt.colorbar()
cb1.ax.set_ylabel('m$^{2}$m$^{-2}$')

plt.text(1000,200, 'LAI', fontsize=12)

plt.subplot(4,2,5)

# peruskartta        
plt.imshow(basemap[0], cmap=basemap[1], alpha=0.8)

x =  LAId / (0.01 + LAIc + LAId)

vmin = np.nanmin(x)
vmax = np.nanmax(x)
plt.imshow(x, extent=extent1, cmap='coolwarm', alpha=0.6)
cb1= plt.colorbar()
cb1.ax.set_ylabel('(-)')
plt.text(1000,200, '$f_d$', fontsize=12)

plt.subplot(4,2,7)

# peruskartta        
plt.imshow(basemap[0], cmap=basemap[1], alpha=0.8)

x = hc

vmin = np.nanmin(x)
vmax = np.nanmax(x)
plt.imshow(x, extent=extent1, cmap='coolwarm', alpha=0.6)
cb1= plt.colorbar()
cb1.ax.set_ylabel('(m)')
plt.text(1000,200, '$h_c$', fontsize=12)
## RHS

# dem

plt.subplot(4,2,4)

# peruskartta        
plt.imshow(basemap[0], cmap=basemap[1], alpha=0.8)

x = dem

vmin = np.nanmin(x)
vmax = np.nanmax(x)
plt.imshow(x, extent=extent1, cmap='coolwarm', alpha=0.6)
cb1= plt.colorbar()
cb1.ax.set_ylabel('(m)')
plt.text(1000,200, 'Elev.', fontsize=12)

plt.subplot(4,2,6)

# peruskartta        
plt.imshow(basemap[0], cmap=basemap[1], alpha=0.8)

x = twi

vmin = np.nanmin(x)
vmax = np.nanmax(x)
plt.imshow(x, extent=extent1, cmap='Blues', alpha=0.8)
cb1= plt.colorbar()
cb1.ax.set_ylabel('(-)')
plt.text(1000,200, 'TWI', fontsize=12)

plt.subplot(4,2,8)

# peruskartta        
plt.imshow(basemap[0], cmap=basemap[1], alpha=0.8)

x = soil

vmin = np.nanmin(x)
vmax = np.nanmax(x)
plt.imshow(x, extent=extent1, cmap= plt.get_cmap('coolwarm_r', 3), alpha=0.6)
cb1= plt.colorbar(ticks=[1, 2, 4])
cb1.ax.set_yticklabels(['coarse','med','peat'], rotation=0.0) 
plt.text(900,200, 'Soiltype', fontsize=12)

plt.savefig('C3_gisdata.png', dpi=600)