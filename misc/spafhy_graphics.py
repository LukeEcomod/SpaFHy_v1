# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 21:07:19 2018

@author: slauniai
"""

import numpy as np
import matplotlib.pyplot as plt

def field_on_map(basemap, x, title=None, bounds=None,
                 xcolormap='RdYlGn', alpha=0.5, vmin=None, vmax=None, fig_nr=None):
    """
    backgroundmap - (pkraster, colormap)
    x - variable n x m np. array
    title - titlestring
    xcolormap - colormap string for x
    """
    #xcolormap = RdYlGn
    #note: bacgroundmap must be cropped to same coordinates than x but resoltion
    # can be different
    r, c = np.shape(basemap[0])
    extent1 = (-0.5, c - 0.5, r - 0.5, -0.5) # extent into where x is re-scaled
    
    if not fig_nr:
        h = plt.figure()
    else:
        h = plt.figure(fig_nr)
    
    if not vmin:
        vmin = np.nanmin(x)
    if not vmax:
        vmax = np.nanmax(x)

    # peruskartta        
    plt.imshow(basemap[0], cmap=basemap[1], alpha=0.8)
    
    plt.imshow(x, extent=extent1, cmap=xcolormap, vmin=vmin, vmax=vmax, alpha=alpha)
    plt.colorbar()
    if bounds is not None:
        plt.imshow(bounds,extent=extent1, cmap='RdYlGn')
    plt.title(title)
    return h
    
def moisture_on_map(basemap, x, title=None,  bounds=None, xcolormap='Blues_r', alpha=0.5, 
                    vmin=None, vmax=None, fig_nr=None):
    """
    FOR PLOTTING MOISTURE MAP.
    backgroundmap - (pkraster, colormap); peruskarttarasteri and its colormap in tuple
    x - clustered saturation deficit, n x m np. array
    title - titlestring
    bounds - raster for plotting calcualtion boundaries
    xcolormap - colormap string used x
    alpha - transparency of x
    vmin, max - for scaling colormap
    fig_nr - figure number
    """
    #xcolormap = RdYlGn
    #note: bacgroundmap must be cropped to same coordinates than x but resoltion
    # can be different
    r, c = np.shape(basemap[0])
    extent1 = (-0.5, c - 0.5, r - 0.5, -0.5) # extent into where x is re-scaled
    
    if not fig_nr:
        h = plt.figure()
    else:
        h = plt.figure(fig_nr)
    
    if not vmin:
        vmin = np.nanmin(x)
    if not vmax:
        vmax = np.nanmax(x)

    # show peruskartta        
    plt.imshow(basemap[0], cmap=basemap[1], alpha=0.8)
    
    ## overlay moisture
    #cmap = plt.cm.get_cmap(xcolormap, 4)
    # 
    plt.imshow(x, extent=extent1, cmap=xcolormap, vmin=vmin, vmax=vmax, alpha=alpha)
    #plt.colorbar()
    # overlay bounds
    if bounds is not None:
        plt.imshow(bounds,extent=extent1, cmap='RdYlGn_r')
    plt.title(title)
    return h