# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 12:29:37 2018

@author: slauniai

This script creates georeferenced peruskarttarasteri clipped to shape of
SpaFHy computation region at Vihti Meolo test site.

"""
import numpy as np
import rasterio
from rasterio.merge import merge
from rasterio.plot import show
from rasterio.mask import mask
import glob
import os
import matplotlib.pyplot as plt
from matplotlib import colors

from shapely.geometry import box
import geopandas as gpd

#%%
def getFeatures(gdf):
    """Function to parse features from GeoDataFrame in such a manner that rasterio wants them"""
    import json
    return [json.loads(gdf.to_json())['features'][0]['geometry']]

#%% read koealamask from asciigrid and save as georeferenced gtif

src_crs=rasterio.crs.CRS.from_epsg(3067) # euref-fin

# ------ read georeferenced catchment mask
ff = r'c:\temp\spafhy\data\vihti\peruskartta\cmask_vihti.dat'
src = rasterio.open(ff)
cmask = src.read(1)
cmask[cmask==src.nodata] = np.NaN
#satdef[np.isnan(satdef)==False] = 1
h, w = np.shape(cmask)
new_data = rasterio.open(r'c:\temp\spafhy\data\vihti\peruskartta\cmask_vihti.tif', 'w', driver='GTiff', 
                         height=h, width=w, count=1, dtype=cmask.dtype, crs=src_crs,
                         transform=src.transform)

new_data.write(cmask,1)
new_data.close()

src.close()

#%% make raster mosaic of 2 peruskarttalehti and save back to tif

dirpath = r'c:\temp\spafhy\data\vihti\peruskartta'
q = os.path.join(dirpath, '*.png')


pngfiles = glob.glob(q)

src_files_to_mosaic = []
for fp in pngfiles:
    src = rasterio.open(fp)
    src_files_to_mosaic.append(src)
    print(src.transform)

mosaic, out_trans = merge(src_files_to_mosaic)

#create colormap
pkcmap = src.colormap(1) # dict, last element in list needs to be removed

tmp = []
for k in pkcmap.keys(): 
#    c.append(list(pkcmap[k][0:3]))
    tmp.append(np.array(pkcmap[k]) / 255.0)

# matplotlib colormap
pk_cmap = colors.ListedColormap(tmp)

# create & update metadata
out_meta = src.meta.copy()

out_meta.update({"driver": "GTiff",
                 "height": mosaic.shape[1],
                 "width": mosaic.shape[2],
                 "transform": out_trans,
                 "crs": src_crs
                 }
                )

# write mosaic to file
out_fp = os.path.join(dirpath,'vihti_whole.tif')
with rasterio.open(out_fp, "w", **out_meta) as dest:
    dest.write(mosaic)

# read back and 
src = rasterio.open(out_fp)
data = src.read(1)

del q, pngfiles, mosaic, out_trans, out_meta, src, src_files_to_mosaic


#%% make raster mosaic of 2 ilmakuva and save back to tif

dirpath = r'c:\temp\spafhy\data\vihti\peruskartta'
q = os.path.join(dirpath, '*.jp2')


pngfiles = glob.glob(q)

src_files_to_mosaic = []
for fp in pngfiles:
    src = rasterio.open(fp)
    src_files_to_mosaic.append(src)
    print(src.transform)
    
mosaic, out_trans = merge(src_files_to_mosaic)

# create & update metadata
out_meta = src.meta.copy()

out_meta.update({"driver": "GTiff",
                 "height": mosaic.shape[1],
                 "width": mosaic.shape[2],
                 "transform": out_trans,
                 "crs": src_crs
                 }
                )

# write mosaic to file
out_fp = os.path.join(dirpath,'vihti_orto_whole.tif')
with rasterio.open(out_fp, "w", **out_meta) as dest:
    dest.write(mosaic)

# read back   
src = rasterio.open(out_fp)
data = src.read(1)

del q, pngfiles, mosaic, out_trans, out_meta, src, src_files_to_mosaic


#%% now read big rasters and clip with smaller

# catchment boundary
cmask = rasterio.open(os.path.join(dirpath,'cmask_vihti.tif'))
#cmask = f.read(1)
# create box shapefile using bbox of cmask: extract box limits and 
# create shapely.geometry.polygon.Polygon
minx, miny, maxx, maxy = cmask.bounds
bbox = box(minx, miny, maxx, maxy)

# insert into geodataframe
geo = gpd.GeoDataFrame({'geometry': bbox}, index=[0], crs=cmask.crs.data)

#re-open gtiff peruskartta
src = rasterio.open(os.path.join(dirpath,'vihti_whole.tif'))
src.meta
#pk = src.read() 

#re-project geodataframe to same coordinate system as raster data
geo = geo.to_crs(crs=src.crs.data)
coords = getFeatures(geo)

#clip raster with polygon and crop
out_img, out_transform = mask(dataset=src, shapes=coords, crop=True)
out_meta = src.meta.copy()
out_meta['height'] = out_img.shape[1]
out_meta['width'] = out_img.shape[2]
out_meta['transform'] = out_transform

with rasterio.open(os.path.join(dirpath, 'peruskartta_vihti.tif'), 'w', **out_meta) as dest:
    dest.write(out_img)

# close orginal raster and open cropped raster
src.close()

#%% repeat with ilmakuva
src = rasterio.open(os.path.join(dirpath,'vihti_orto_whole.tif'))
src.meta
#orto = src.read(1) 

## clip using bbox of cmask: extract box limits and create shapely.geometry.polygon.Polygon
#minx, miny, maxx, maxy = cmask.bounds
#bbox = box(minx, miny, maxx, maxy)

#re-project geodataframe to same coordinate system as raster data
geo = geo.to_crs(crs=src.crs.data)
coords = getFeatures(geo)

#clip raster with polygon and crop
out_img, out_transform = mask(dataset=src, shapes=coords, crop=True)
out_meta = src.meta.copy()
out_meta['height'] = out_img.shape[1]
out_meta['width'] = out_img.shape[2]
out_meta['transform'] = out_transform

with rasterio.open(os.path.join(dirpath, 'orto_vihti.tif'), 'w', **out_meta) as dest:
    dest.write(out_img)

# close orginal raster and open cropped raster
src.close()

#%% create pickle of all necessary input files
cmask = rasterio.open(os.path.join(dirpath,'cmask_vihti.tif'))

data = cmask.read(1)
r, c = np.shape(data)

boundaries = np.ones([r,c])*np.NaN
for k in range(r):
    ix = np.where(data[k,:] > 0)[0]
    #print(type(ix), np.shape(ix))
    #print(min(ix[0]), max(ix[0]))
    if len(ix) > 0:
        boundaries[k, min(ix)] = 1.0
        boundaries[k, max(ix)] = 1.0
    
    #m, n = min(ix, max(ix))


import pickle

sf = os.path.join(dirpath, 'peruskartta_vihti.tif')
pkraster = rasterio.open(sf)

pk = pkraster.read(1)

sf = os.path.join(dirpath, 'orto_vihti.tif')
ortof = rasterio.open(sf)

orto = ortof.read()

# dump into pickle
outfile = open(os.path.join(dirpath, 'vihti_rasters.pk'), 'wb')
out = {'peruskartta': pk, 'orto': orto, 'boundaries': boundaries, 'pkcmap': pk_cmap}
pickle.dump(out, outfile)
outfile.close()