#!/usr/bin/env python
u"""
polygon.py
Written by Tyler Sutterley (10/2020)
Reads polygons from GeoJSON, kml/kmz or ESRI shapefile files

INPUTS:
    input polygon file

OUTPUT:
    shapely multipolygon object of input file

OPTIONS:
    EPSG: projection identifier for output coordinates
    VARIABLES: reduce to a specific set of identifiers

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    fiona: Python wrapper for vector data access functions from the OGR library
        https://fiona.readthedocs.io/en/latest/manual.html
    geopandas: Python tools for geographic data
        http://geopandas.readthedocs.io/
    shapely: PostGIS-ish operations outside a database context for Python
        http://toblerity.org/shapely/index.html
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/

UPDATE HISTORY:
    Written 10/2020
"""
from __future__ import print_function

import os
import io
import re
import fiona
import pyproj
import zipfile
import osgeo.gdal
import geopandas
import numpy as np
from shapely.geometry import Polygon, MultiPolygon
# enable kml driver for geopandas
fiona.drvsupport.supported_drivers['LIBKML'] = 'rw'

class polygon(object):
    """
    Data class for reading polygon files
    """
    np.seterr(invalid='ignore')
    def __init__(self, epsg=4326):
        self.filename=None
        self.epsg=epsg

    def case_insensitive_filename(self,filename):
        """
        Searches a directory for a filename without case dependence
        """
        self.filename = os.path.expanduser(filename)
        #-- check if file presently exists with input case
        if not os.access(self.filename,os.F_OK):
            #-- search for filename without case dependence
            basename = os.path.basename(filename)
            directory = os.path.dirname(os.path.expanduser(filename))
            f = [f for f in os.listdir(directory) if re.match(basename,f,re.I)]
            if not f:
                raise IOError('{0} not found in file system'.format(filename))
            self.filename = os.path.join(directory,f.pop())
        return self

    def from_geojson(self, filename, variables=None):
        """
        read GeoJSON (.json, .geojson) files
        """
        # set filename
        self.case_insensitive_filename(filename)
        # read the GeoJSON file
        gj = geopandas.read_file(self.filename)

        #-- converting x,y from polygon projection to output EPSG
        crs1 = pyproj.CRS.from_string(gj.crs['init'])
        crs2 = pyproj.CRS.from_string("epsg:{0:d}".format(self.epsg))
        transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)

        # list of polygons
        poly_list = []
        # find features of interest
        geometries = ('LineString','Polygon')
        f = [f for f in gj.iterfeatures() if f['geometry']['type'] in geometries]
        # reduce to variables of interest if specified
        f = [ft for ft in f if ft['id'] in variables] if variables else f

        # for each line string or polygon feature
        for feature in f:
            # extract coordinates for feature
            x,y = np.transpose(feature['geometry']['coordinates'])
            # convert points to EPSG
            xi,yi = transformer.transform(x, y)
            # create shapely polygon
            poly_obj = Polygon(np.c_[xi,yi])
            # cannot have overlapping exterior or interior rings
            if (not poly_obj.is_valid):
                poly_obj = poly_obj.buffer(0)
            poly_list.append(poly_obj)
        # create shapely multipolygon object
        # return the polygon object
        return MultiPolygon(poly_list)

    def from_kml(self, filename, kmz=False, variables=None):
        """
        read keyhole markup language (.kml) files
        """
        # set filename
        self.case_insensitive_filename(filename)
        # if input file is compressed
        if kmz:
            # decompress and parse KMZ file
            z = zipfile.ZipFile(self.filename, 'r')
            kml_file, = [s for s in z.namelist() if re.search(r'\.(kml)$',s)]
            # need to use osgeo virtual file system to add suffix to mmap name
            mmap_name = "/vsimem/{0}".format(kml_file)
            osgeo.gdal.FileFromMemBuffer(mmap_name, z.read(kml_file))
            with fiona.Collection(mmap_name, driver='LIBKML') as f:
                kml = geopandas.GeoDataFrame.from_features(f, crs=f.crs)
        else:
            kml = geopandas.read_file(self.filename)

        #-- converting x,y from polygon projection to output EPSG
        crs1 = pyproj.CRS.from_string(kml.crs['init'])
        crs2 = pyproj.CRS.from_string("epsg:{0:d}".format(self.epsg))
        transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)

        # list of polygons
        poly_list = []

        # find features of interest
        geometries = ('LineString','Polygon')
        f = [f for f in kml.iterfeatures() if f['geometry']['type'] in geometries]
        # reduce to variables of interest if specified
        f = [ft for ft in f if ft['id'] in variables] if variables else f

        # for each line string or polygon feature
        for feature in f:
            # extract coordinates for feature
            coords = np.squeeze(feature['geometry']['coordinates'])
            # convert points to EPSG
            xi,yi = transformer.transform(coords[:,0], coords[:,1])
            # create polygon from coordinate set
            poly_obj = Polygon(np.c_[xi,yi])
            # cannot have overlapping exterior or interior rings
            if (not poly_obj.is_valid):
                poly_obj = poly_obj.buffer(0)
            poly_list.append(poly_obj)
        # create shapely multipolygon object
        # return the polygon object
        return MultiPolygon(poly_list)

    def from_shapefile(self, filename, zip=False, variables=None):
        """
        read ESRI shapefiles
        """
        # set filename
        self.case_insensitive_filename(filename)
        # read input zipfile containing shapefiles
        if zip:
            # read the compressed shapefile and extract entities
            shape = fiona.open('zip://{0}'.format(self.filename))
        else:
            # read the shapefile and extract entities
            shape = fiona.open(self.filename,'r')

        #-- converting x,y from polygon projection to output EPSG
        crs1 = pyproj.CRS.from_string(shape.crs['init'])
        crs2 = pyproj.CRS.from_string("epsg:{0:d}".format(self.epsg))
        transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)

        # find features of interest
        geometries = ('LineString','Polygon')
        f = [f for f in shape.values() if f['geometry']['type'] in geometries]
        # reduce to variables of interest if specified
        f = [ft for ft in f if ft['id'] in variables] if variables else f

        # list of polygons
        poly_list = []
        # for each entity
        for i,ent in enumerate(f):
            # extract coordinates for entity
            for coords in ent['geometry']['coordinates']:
                # convert points to latitude/longitude
                x,y = np.transpose(coords)
                # convert points to EPSG
                xi,yi = transformer.transform(x, y)
                # create shapely polygon
                poly_obj = Polygon(np.c_[xi,yi])
                # cannot have overlapping exterior or interior rings
                if (not poly_obj.is_valid):
                    poly_obj = poly_obj.buffer(0)
                poly_list.append(poly_obj)
        # create shapely multipolygon object
        # return the polygon object
        return MultiPolygon(poly_list)
