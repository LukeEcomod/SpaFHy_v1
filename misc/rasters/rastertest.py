# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 11:07:03 2018

@author: slauniai
"""
import os
import numpy as np
import pickle
import rasterio as rs
import matplotlib.pyplot as plt
from matplotlib import colors
from rasterio.plot import show
from rasterio.mask import mask

from shapely.geometry import box
import geopandas as gpd
from fiona.crs import from_epsg
import pycrs


#%%
def getFeatures(gdf):
    """Function to parse features from GeoDataFrame in such a manner that rasterio wants them"""
    import json
    return [json.loads(gdf.to_json())['features'][0]['geometry']]

#%%
# define raster coordinates; these are not conded in png-files
src_crs=rs.crs.CRS.from_epsg(3067) # euref-fin

# ------ read georeferenced catchment mask
ff = r'c:\temp\spafhy\rastertests\ch3_mask.tif'
cmask = rs.open(ff)


# ------ peruskartta read raster png (in epsg=3067) and save back to raster file
mapfolder = r'c:\temp\spafhy\data\3\peruskartta'
pkfile = r'Q5321L.png'

ff = os.path.join(mapfolder, pkfile)

# open for reading
raster = rs.open(ff, crs=src_crs)
raster.meta

# read peruskarttarasteri
pk = raster.read()

# read and create colormap for basemap raster
pkcmap = raster.colormap(1) # dict, last element in list needs to be removed

tmp = []
for k in pkcmap.keys(): 
#    c.append(list(pkcmap[k][0:3]))
    tmp.append(np.array(pkcmap[k]) / 255.0)

# matplotlib colormap
pk_cmap = colors.ListedColormap(tmp)


plt.figure()
show(pk, cmap=pk_cmap)

new_data = rs.open(r'c:\temp\spafhy\rastertests\Q5321L_epsg3067.tif', 'w', driver='GTiff', height=12000, width=12000,
                   count=1, dtype=pk.dtype, crs=src_crs, transform=raster.transform )

new_data.write(pk[0],1)
new_data.close()

# --- close files
raster.close()

#re-open gtiff to have crs defined
raster = rs.open(r'c:\temp\spafhy\rastertests\Q5321L_epsg3067.tif', crs=src_crs)
raster.meta
pk = raster.read(1) 
del pk

# clip using bbox of cmask: extract box limits and create shapely.geometry.polygon.Polygon
minx, miny, maxx, maxy = cmask.bounds
bbox = box(minx, miny, maxx, maxy)

# insert into geodataframe
geo = gpd.GeoDataFrame({'geometry': bbox}, index=[0], crs=cmask.crs.data)

#re-project to same coordinate system as raster data
geo = geo.to_crs(crs=raster.crs.data)
coords = getFeatures(geo)

#clip raster with polygon and crop
out_img, out_transform = mask(dataset=raster, shapes=coords, crop=True)
out_meta = raster.meta.copy()
out_meta['height'] = out_img.shape[1]
out_meta['width'] = out_img.shape[2]
out_meta['transform'] = out_transform

with rs.open(r'c:\temp\spafhy\rastertests\peruskartta_ch3.tif', 'w', **out_meta) as dest:
    dest.write(out_img)

# close orginal raster and open cropped raster
raster.close()

sf = r'c:\temp\spafhy\rastertests\peruskartta_ch3.tif'
pkraster = rs.open(sf)

pk = pkraster.read(1)

# dump into pickle
outfile = open(r'c:\temp\spafhy\rastertests\ch3_peruskartta.pk', 'wb')
pickle.dump((pk, tmp), outfile)
outfile.close()

#%% repeat same for spafhy-results
sf = r'c:\temp\spafhy\Slocal365.dat'
ras2 = rs.open(sf)

satdef = ras2.read(1)
satdef[satdef==ras2.nodata] = np.NaN
#satdef[np.isnan(satdef)==False] = 1
h, w = np.shape(satdef)
new_data = rs.open(r'c:\temp\spafhy\rastertests\satdef.tif', 'w', driver='GTiff', height=h, width=w,
                   count=1, dtype=satdef.dtype, crs=src_crs, transform=ras2.transform )

new_data.write(satdef,1)
new_data.close()

ras2.close()

#%% re-open gtiff
ras2 = rs.open(r'c:\temp\spafhy\rastertests\satdef.tif')
ras2.meta
satdef = ras2.read(1)

#%%
pkraster.bounds
r, c = np.shape(pkraster)
extent1 = (-0.5, c - 0.5, r - 0.5, -0.5)
r1, c1 = np.shape(ras2)
extent2 = (-0.5, c - 0.5, r - 0.5, -0.5)

#%% reproject spathy result raster
#
#src = ras2
#src_crs=ras2.crs
#src_transform=ras2.transform
#
#dst_crs='EPSG:3067'
#destin = raster
#profile = src.profile.copy()
#
#transform, width, height = rs.warp.calculate_default_transform(
#    destin.crs, dst_crs, destin.width, destin.height, *destin.bounds)
#
#profile.update({
#    'crs': dst_crs,
#    'transform': transform,
#    'width': width,
#    'height': height
#})
#
#with rs.open(r'c:\temp\spafhy\rastertests\aaa.tif', 'w', **profile) as dst:
#    rs.warp.reproject(
#        source=rs.band(src, 1),
#        destination=rs.band(dst, 1),
#        src_transform=src.transform,
#        src_crs=src.crs,
#        dst_transform=transform,
#        dst_crs=dst_crs,
#        resampling=rs.warp.Resampling.bilinear)
#    
#
#ras3 = rs.open(r'c:\temp\spafhy\rastertests\aaa.tif')
#w,h = np.shape(satdef)
#new_data = rs.open('satdef.tif', 'w', driver='GTiff', height=h, width=w,
#                   count=1, dtype=satdef.dtype, crs=src_crs, transform=raster.transform )
#
#new_data.write(satdef,1)
#new_data.close()
#
#ras2.close()
#
##%% re-open gtiff
#ras2 = rs.open(r'c:\temp\spafhy\rastertests\satdef.tif', crs=src_crs)
#ras2.meta
#data = ras2.read(1)
#ras2.close()
##%%
#src = rs.open(r'c:\temp\spafhy\Data\3\peruskartta\Q5321D.jp2')
#src.meta
#
##%%
#
#new_data = rs.open('testiout.tif', 'w', driver='GTiff', height=12000, width=12000,
#                   count=1, dtype=pk.type, crs=src_crs, transform=raster.transform )