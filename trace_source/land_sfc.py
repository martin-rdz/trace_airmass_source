#! /usr/bin/env python3
# coding=utf-8
""""""
"""
Author: radenz@tropos.de

"""

import os
import rasterio
import numpy as np
from affine import Affine
from pyproj import Proj, transform
from fastkml import kml
from shapely.geometry.point import Point
import shapely.vectorized
import toml

from numba import jit

@jit(nopython=True)
def nearest(point, array, delta):
    """
    calculate the index directly

    Args:
        point (float): point to search for
        array (np.array): array to search in
        delta: step size in array

    Returns: 
        ``(i, array[i])``
    """
    ## i = bisect.bisect_left(array, point)
    #i = int((point - array[0]) / delta)
    ## print("search nearest ", i, point, " | ", array[max(0,i-3):i+4])

    #nearest = min(array[max(0, i - 10):i + 10], key=lambda t: abs(point - t))

    #i = np.where(array == nearest)[0][0]

    i_calc = min(np.round((point - array[0])/delta), array.shape[0]-1)
    #print(i, i_calc, point, array[i-1:i+2], (point - array[0])/delta)
    #assert i == i_calc
    i = int(i_calc)
    #return (i, nearest)
    return (i, array[i])


def argnearest(value, array, delta):
    """
    """
    i = np.searchsorted(array, value) - 1
    if not i == array.shape[0] - 1:
        if np.abs(array[i] - value) > np.abs(array[i + 1] - value):
            i = i + 1

    return (i, array[i])


#@jit(nopython=True)
def fast_land_sfc(land_sfc_data, lat, lon, lats, longs):
    land_sfc_category = np.zeros((len(lat),))

    land_sfc_category_new = np.zeros((len(lat),))
    ilat = np.round((lat - lats[0])/-0.1).astype(int)
    ilat[ilat >= lats.shape[0]-1] = lats.shape[0]-1
    ilat[ilat < 0] = 0
    ilon = np.round((lon - longs[0])/0.1).astype(int)
    ilon[ilon >= longs.shape[0]-1] = longs.shape[0]-1
    ilon[ilon < 0] = 0
    land_sfc_category = land_sfc_data[ilat, ilon]
    land_sfc_category[lat == -999.] = -1

    return land_sfc_category


class land_sfc():
    """
    load the data from the geotiff ``data/resampledLCType.tif``
    and concatenate to the simplified scheme
    
    ===============   =========================
    MODIS Category     simplified
    ===============   =========================
    0                  0 water
    1,2,3,4,5,6        1 forest
    7,8,9              2 savanna, shrubland
    10, 11, 12, 14     3 grass, cropland
    13                 4 urban
    15                 5 snow
    16                 6 barren
    ===============   =========================

    """
    def __init__(self):
        filename =  os.path.dirname(os.path.abspath(__file__)) +\
                    '/../data/resampledLCType.tif'

        #with rasterio.divers():
        with rasterio.open(filename, 'r') as src:
            meta = src.meta
            width = meta['width']
            height = meta['height']
            count = meta['count']
            dtype = meta['dtype']
            self.shape = src.shape
            self.transform = src.transform

            T0 = src.transform
            p1 = Proj(src.crs)
            print('T0 aka affine transformation ', T0)
            print('src.crs', src.crs)
            print('p1', p1)

            # allocate memory for image
            im = np.empty([height, width], dtype)

            # read image into memory
            print(meta)
            im[:, :] = src.read(1)


        T1 = T0 * Affine.translation(0.5, 0.5)
        # Function to convert pixel row/column index (from 0) to easting/northing at centre
        print('affine transformation', T1)

        etemp, northings = T1 * (np.meshgrid(np.arange(1), np.arange(im.shape[0])))
        eastings, ntemp = T1 * (np.meshgrid(np.arange(im.shape[1]), np.arange(1)))
        #print(eastings, northings)

        # the geotiff provided is already in WGS84
        assert src.crs == 'EPSG:4326'
        # # Project all longitudes, latitudes
        # p2 = Proj(proj='latlong', datum='WGS84')
        # _, lats = transform(p1, p2, etemp, northings)
        # longs, _ = transform(p1, p2, eastings, ntemp)


        # slow version of the coordinate calculation
        # cols, rows = np.meshgrid(np.arange(im.shape[1]), np.arange(im.shape[0]))
        # # Get affine transform for pixel centres
        # T1 = T0 * Affine.translation(0.5, 0.5)
        # # Function to convert pixel row/column index (from 0) to easting/northing at centre
        # rc2en = lambda r, c: (c, r) * T1
        # # All eastings and northings (there is probably a faster way to do this)
        # eastings, northings = np.vectorize(rc2en, otypes=[np.float, np.float])(rows, cols)
        # #print("eastings", eastings.nbytes, eastings)
        # #print("northings", northings.nbytes, northings)
        # # Project all longitudes, latitudes
        # p2 = Proj(proj='latlong', datum='WGS84')
        # longs, lats = transform(p1, p2, eastings, northings)

        if (eastings == eastings[0, :]).all():
            self.longs = eastings[0, :]
        if (northings.T == northings[:, 0]).all():
            self.lats = northings[:, 0]

        self.land_sfc_data = im
        # high green
        self.land_sfc_data[(im >= 1) & (im <= 6)] = 1
        # med green
        self.land_sfc_data[(im >= 7) & (im <= 9)] = 2
        # low green
        self.land_sfc_data[(im >= 10) & (im <= 12)] = 3
        self.land_sfc_data[im == 14] = 3
        # urban
        self.land_sfc_data[im == 13] = 4
        # snow/ice
        self.land_sfc_data[im == 15] = 5
        # desert
        self.land_sfc_data[im == 16] = 6

        self.categories = {0:'water',1:'forest',2:'savanna/shrub',
                           3:'grass/crop',4:'urban',5:'snow',6:'barren'}

    #@jit(nopython=True)
    #@profile
    def get_land_sfc(self, lat, lon):
        """
        get the land use pixel for a given coordinate
        interpolation to the nearest pixel is done

        Args:
            lat (float, array): latitude
            lon (float, array): longitude

        Returns:
            array or int with the land use category
        """

        if isinstance(lat, float):
            lat = [lat]
            lon = [lon]
        if isinstance(lat, np.ma.core.MaskedArray):
            lat = lat.filled(-999.).tolist()
            lon = lon.filled(-999.).tolist()
        elif isinstance(lat, np.ndarray):
            lat = lat.tolist()
            lon = lon.tolist()
        # im[lat, lon]

        land_sfc_data = self.land_sfc_data
        if lat:
            cats = fast_land_sfc(land_sfc_data, lat, lon, self.lats, self.longs)
        else:
            cats = np.array([])
        return cats




    def get_land_sfc_shape(self, shp):
        """
        get the land use inside a given shape

        :param shp: shapely Multipolygon

        .. deprecated:: 0.1
            use the ensemble trajectories instead
        """

        import rasterio.features

        mask = rasterio.features.rasterize(
            shp,
            out_shape=self.shape,
            transform=self.transform)

        masked_land_sfc = np.ma.MaskedArray(self.land_sfc,
                                            mask=np.logical_not(mask),
                                            copy=True)

        #print(masked_land_sfc.compressed())
        return masked_land_sfc


@jit()
def fast_geonames(geo_names_data, polygons, lat, lon):

    geo_id = np.zeros((len(lat),))
    geo_id[:] = -1

    for j, name in geo_names_data.items():
        is_within = shapely.vectorized.contains(polygons[name], lon, lat)
        geo_id[is_within] = j

    return geo_id


class named_geography():
    """
    handle a single trajectory
    
    Args:
        selected_type (str): name in the ``geonames_config.toml``
    """
    def __init__(self, selected_type):
        config_file = 'geonames_config.toml'
        with open(config_file) as config_file:
            self.config = toml.loads(config_file.read())

        filename = os.path.dirname(os.path.abspath(__file__)) +\
                   '/../' + self.config[selected_type]['filename']
        # filename =  os.path.dirname(os.path.abspath(__file__)) +\
        #            '/../prior_examples/geo_names.kml'

        k = kml.KML()
        with open(filename, 'r') as f:
            k.from_string(bytes(bytearray(f.read(), encoding='utf-8')))

        docu = list(k.features())[0]
        polygons = {}
        for p in list(docu.features()):
            print(p.name)
            print(p.geometry)
            polygons[p.name] = p.geometry

        self.polygons = polygons
        # self.geo_names = {0: 'cont_europe', 1: 'sahara', 2: 'arabian_peninsula',
        #                   3: 'far_east_deserts', 4: 'persia', 5: 'india'}
        self.geo_names = {int(k): v for k, v in self.config[selected_type]['geo_names'].items()}
        assert set(polygons.keys()) == set(self.geo_names.values()), \
            ('not all keys in polygon are matched', polygons.keys(), self.geo_names.values())
        print('loaded the named_geography')
        print('available geo names ', self.geo_names)


    #@jit(nopython=True)
    #@profile
    def get_geo_names(self, lat, lon):
        """
        get the names defined in the shapefile for a given coordinate


        Args:
            lat (float, array): latitude
            lon (float, array): longitude

        Returns:
            array or int with the geoname category
        """

        if isinstance(lat, float):
            lat = [lat]
            lon = [lon]
        if isinstance(lat, np.ma.core.MaskedArray):
            lat = lat.filled(-999.).tolist()
            lon = lon.filled(-999.).tolist()
        elif isinstance(lat, np.ndarray):
            lat = lat.tolist()
            lon = lon.tolist()
        # im[lat, lon]

        geo_id = fast_geonames(self.geo_names, self.polygons, lat, lon)

        return geo_id


if __name__ == '__main__':
    ng = named_geography()
    print(ng.get_geo_names(19.67, 22.16))
    print(ng.get_geo_names(25.63, 42.00))
