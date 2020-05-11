# -*- coding: utf-8 -*-
"""
Created on Thu Apr 05 16:12:47 2018

@author: slauniai
"""

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#%matplotlib inline  
import warnings
import matplotlib.cbook
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)

chs = pd.read_csv(r'c:\datat\spathydata\sve_catchment_characteristics.csv', sep=';', encoding = "ISO-8859-1")
f = chs[(chs.id == 22) |  (chs.id == 30)].index.values
chs = chs.drop(f)
chs['xoffs'] = -80000
chs['yoffs'] = 40000
chs['xoffs'].loc[chs.id == 2] = -90000
chs['yoffs'].loc[chs.id == 2] = -20000
chs['xoffs'].loc[chs.id == 31] = 20000
chs['xoffs'].loc[chs.id == 32] = 20000
chs['yoffs'].loc[chs.id == 32] = -40000
chs['yoffs'].loc[chs.id == 13] = -90000
chs['xoffs'].loc[chs.id == 11] = 10000
chs['xoffs'].loc[chs.id == 11] = 10000
chs['xoffs'].loc[chs.id == 3] = 20000
chs['yoffs'].loc[chs.id == 14] = -90000
chs['yoffs'].loc[chs.id == 29] = 20000
del f

ecs = pd.read_csv(r'c:\datat\spathydata\ec_sites_hess.csv', sep=';',encoding = "ISO-8859-1")
ecs = ecs.drop([1, 2, 3, 6])
ec_sitenames = ['FIHy\nFISII\nFICage', 'FISod', 'FIKal\nFILet', 'SEKno', 'SENor\nSESky2']

# vienna cooords
wien_lon = 16.36
wien_lat = 48.21

#%% plot large map
fig = plt.figure(num=None, figsize=(12, 8))
m = Basemap(projection='merc', llcrnrlat=35.0, urcrnrlat=72.0, llcrnrlon=-13, urcrnrlon=37, 
            lat_0=55.0, lon_0=10.0, resolution='i')
m.drawcoastlines(linewidth=0.3)
m.drawcountries(linewidth=0.5, linestyle='solid', color='k' )
m.drawmapboundary(linewidth=1, color='k')
m.drawparallels(np.arange(40, 80, 10), linewidth=0.5, labels=[0,1,0,0]) 
# m.drawrivers(linewidth=0.5, linestyle='solid', color='blue')

# plot EC sites as points
xec, yec = m(ecs.LON_deg.values, ecs.LAT_deg.values)
#x, y = 62.0, 23.0
m.plot(xec, yec, 'bo', markersize=7, alpha=0.7, label='EC sites')
#m.plot(xec[1:], yec[1:], 'bo', markersize=5, alpha=0.7) #, label=tt) #, label=idcode)

#for lbl, xpt, ypt, in zip(ecs.Site, xec, yec):
#    plt.text(xpt + 10000, ypt + 10000, lbl, fontsize=8)
    
# plot catchments as points
x, y = m(chs.LON_deg.values, chs.LAT_deg.values)
#x, y = 62.0, 23.0
labeltexts = chs.id
m.plot(x[0], y[0], 'ro', markersize=5, alpha=0.7, label='Catchments') #, label=tt) #, label=idcode)
m.plot(x[1:], y[1:], 'ro', markersize=5, alpha=0.7) #, label=tt) #, label=idcode)
#for lbl, xpt, ypt, in zip(chs.id.values, x, y):
#    plt.text(xpt + 10000, ypt + 10000, lbl, fontsize=8)

    
x0, y0 = m(wien_lon, wien_lat)
m.plot(x0, y0, 'go', markersize=10)
plt.text(x0 - 100000, y0 - 300000, 'Vienna', fontsize=14)
m.fillcontinents(alpha=0.3)
plt.legend(loc=2)
plt.show()
plt.savefig('Sitemap_large.png', dpi=300)

#%% plot zoomed map
fig, ax = plt.subplots()
fig.set_size_inches(11, 8)
m = Basemap(projection='merc', llcrnrlat=55.0, urcrnrlat=72.0, llcrnrlon=10., urcrnrlon=37, 
            lat_0=65.0, lon_0=25.0, resolution='l')
m.drawcoastlines(linewidth=0.3)
m.drawcountries(linewidth=0.5, linestyle='solid', color='k' )
m.drawmapboundary(linewidth=1, color='k')
m.drawparallels(np.arange(40, 80, 10), linewidth=0.5, labels=[0,1,0,0]) 
# m.drawrivers(linewidth=0.5, linestyle='solid', color='blue')

# plot catchments as points
x, y = m(chs.LON_deg.values, chs.LAT_deg.values)
offsets = np.ones((len(x),2))*10000

m.plot(x, y, 'ro', markersize=5, alpha=0.7, label='Catchments') #, label=tt) #, label=idcode)
#m.plot(x[1:], y[1:], 'ro', markersize=5, alpha=0.7) #, label=tt) #, label=idcode)
for lbl, xpt, ypt, xof, yof in zip(chs.id.values, x, y, chs.xoffs, chs.yoffs):
    plt.text(xpt + xof, ypt + yof, lbl, fontsize=10)


# plot EC sites as points
mrk = ['o', 's','^','v','<','>']

xec, yec = m(ecs.LON_deg.values, ecs.LAT_deg.values)
#x, y = 62.0, 23.0
#colors = plt.cm.viridis(np.linspace(0.01, 0.9, len(xec)))
#ax.set_prop_cycle('color', colors)
nn = 0
for k in range(len(xec)):
    m.plot(xec[k], yec[k], markeredgecolor='k', marker=mrk[nn], markersize=8, linewidth=0.0, alpha=0.5, label=ec_sitenames[k]) #, label=tt) #, label=idcode)
    #nn += 1

#m.fillcontinents(alpha=0.05)

plt.legend(loc=2)
plt.show()
plt.savefig('Sitemap.png', dpi=600)
