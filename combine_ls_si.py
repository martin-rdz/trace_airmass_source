

import datetime
import argparse
from pathlib import Path

import numpy as np
import rasterio
from rasterio.io import MemoryFile
from rasterio.warp import calculate_default_transform, reproject, Resampling
from affine import Affine
# from pyproj import Proj, transform


parser = argparse.ArgumentParser(usage="usage")
#parser.add_option("interval", help="interval %Y%m%d-%Y%m%d")
parser.add_argument('--date', help='date in the format YYYYMMDD')
args = parser.parse_args()

print('args', args)
day = datetime.datetime.strptime(args.date, "%Y%m%d")

basepath_source = 'data/sea_ice_amsr/' 
f_source = basepath_source + f'asi-AMSR2-s6250-{day:%Y%m%d}-v5.4.tif'
fpath = Path(f_source)
f_out = fpath.with_name(fpath.stem + '_reproj2' + fpath.suffix)
print('output of reprojected ', f_out)
f_ls_shelf = 'data/resampledLCType_custclass_shelfs.tif'
f_final = f'data/ls_si/{day:%Y%m%d}_class_with_si.tif'
print('output final', f_final)


## reprojecting

with rasterio.open(f_source, 'r') as src:
    print(src)
    print(src.meta)
    meta = src.meta
    width = meta['width']
    height = meta['height']
    count = meta['count']
    dtype = meta['dtype']
    self_shape = src.shape
    self_transform = src.transform

    T0 = src.transform
    #p1 = Proj(src.crs)
    print('T0 aka affine transformation ', T0)
    print('src.crs', src.crs)
    print('src.bounds', *src.bounds)

    crs_wgs84 = rasterio.crs.CRS.from_string('EPSG:4326')
    # transform, width, height = calculate_default_transform(
    #    src.crs, crs_wgs84, src.width, src.height, *src.bounds, resolution=(0.1,0.1))
    # print('target transform \n', transform)
    # print('width, height', width, height)
    print('overwrite to the get the full globe')
    transform = Affine(0.1,0,-180, 0,-0.1,90)
    width, height = 3600, 1800
    
    kwargs = src.meta.copy()
    kwargs.update({
        'crs': crs_wgs84,
        'transform': transform,
        'width': width,
        'height': height})
    

    with MemoryFile() as memfile:
        with memfile.open(**kwargs) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=crs_wgs84,
                    resampling=Resampling.nearest,
                    dst_nodata=130
                    )
    
            print('dst.crs', dst.crs)
            meta = dst.meta
            width = meta['width']
            height = meta['height']
            count = meta['count']
            dtype = meta['dtype']
            self_shape = dst.shape
            self_transform = dst.transform

            T0 = dst.transform
            print('T0 aka affine transformation ', T0)
            # allocate memory for image
            im = np.empty([height, width], dtype)
            # read image into memory
            im[:, :] = dst.read(1)

            # for diagnostics save an output file
            profile = dst.profile
            
            with rasterio.open(f_out, 'w', **profile) as dst2:
                dst2.write(im, 1)


## combining
with rasterio.open(f_ls_shelf) as match:
    meta = match.meta
    width = meta['width']
    height = meta['height']
    count = meta['count']
    dtype = meta['dtype']

    ls = np.empty([height, width], dtype)
    ls[:, :] = match.read(1)

    ls_profile = match.profile

with rasterio.open(f_out) as dst:

    meta = dst.meta
    width = meta['width']
    height = meta['height']
    count = meta['count']
    dtype = meta['dtype']
    im = np.empty([height, width], dtype)
    # read image into memory
    im[:, :] = dst.read(1)
    sea_ice_cover = im.copy()

    si_profile = match.profile

# open ocean is also fine
# ls_profile[(sea_ice_cover >= 0) & (sea_ice_cover < 1)] = 1
ls[(sea_ice_cover >= 1) & (sea_ice_cover < 30)] = 7
ls[(sea_ice_cover >= 30) & (sea_ice_cover < 80)] = 8
ls[(sea_ice_cover >= 80) & (sea_ice_cover < 101)] = 9
# the mask class should not be needed
# sea_ice_cover[(sea_ice_cover > 100) & (sea_ice_cover < 127)] = 5

with rasterio.open(f_final, 'w', **si_profile) as dst2:
    dst2.write(ls, 1)
