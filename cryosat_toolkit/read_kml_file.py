#!/usr/bin/env python
u"""
read_kml_file.py
Written by Tyler Sutterley (09/2019)
Reads polygons from keyhole markup language (.kml or .kmz) files

INPUTS:
    input .kml or .kmz file

OUTPUT:
    shapely multipolygon object of input file

OPTIONS:
    EPSG: projection identifier for output coordinates
    KMZ: input file is compressed
    VARIABLES: reduce to a specific set of identifiers

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        http://www.numpy.org
        http://www.scipy.org/NumPy_for_Matlab_Users
    gdal: Pythonic interface to the Geospatial Data Abstraction Library (GDAL)
        https://pypi.python.org/pypi/GDAL/
    fiona: Python wrapper for vector data access functions from the OGR library
        https://fiona.readthedocs.io/en/latest/manual.html
    geopandas: Python tools for geographic data
        http://geopandas.readthedocs.io/
    shapely: PostGIS-ish operations outside a database context for Python
        http://toblerity.org/shapely/index.html
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/

UPDATE HISTORY:
    Updated 09/2019: made output coordinate system an option (EPSG)
    Updated 07/2019: added option to reduce to specific VARIABLES within file
    Updated 06/2019: convert projection to EPSG:4326 before creating polygons
        only read LineString and Polygon features from the kml/kmz file
    Written 06/2019
"""
from __future__ import print_function

import os
import io
import re
import zipfile
import fiona
import pyproj
import geopandas
import osgeo.gdal
import numpy as np
#-- enable kml driver for geopandas
fiona.drvsupport.supported_drivers['LIBKML'] = 'rw'
from shapely.geometry import Polygon, MultiPolygon

#-- PURPOSE: read keyhole markup language (.kml) files
def read_kml_file(input_file, EPSG=4326, KMZ=False, VARIABLES=None):
    #-- if input file is compressed
    if KMZ:
        #-- decompress and parse KMZ file
        kmz = zipfile.ZipFile(os.path.expanduser(input_file), 'r')
        kml_file, = [s for s in kmz.namelist() if re.search('\.(kml)$',s)]
        #-- need to use osgeo virtual file system to add suffix to mmap name
        mmap_name = "/vsimem/{0}".format(kml_file)
        osgeo.gdal.FileFromMemBuffer(mmap_name, kmz.read(kml_file))
        with fiona.Collection(mmap_name, driver='LIBKML') as f:
            kml = geopandas.GeoDataFrame.from_features(f, crs=f.crs)
    else:
        kml = geopandas.read_file(os.path.expanduser(input_file))

    #-- convert projection to EPSG
    proj1 = pyproj.Proj("+init={0}".format(kml.crs['init']))
    proj2 = pyproj.Proj("+init=EPSG:{0:d}".format(EPSG))

    #-- list of polygons
    poly_list = []

    #-- find features of interest
    geometries = ('LineString','Polygon')
    f = [f for f in kml.iterfeatures() if f['geometry']['type'] in geometries]
    #-- reduce to variables of interest if specified
    f = [ft for ft in f if ft['id'] in VARIABLES] if VARIABLES else f

    #-- for each line string or polygon feature
    for feature in f:
        #-- extract coordinates for feature
        coords = np.squeeze(feature['geometry']['coordinates'])
        #-- convert points to latitude/longitude
        lon,lat = pyproj.transform(proj1, proj2, coords[:,0], coords[:,1])
        #-- create polygon from coordinate set
        poly_obj = Polygon(list(zip(lon, lat)))
        #-- Valid Polygon cannot have overlapping exterior or interior rings
        if (not poly_obj.is_valid):
            poly_obj = poly_obj.buffer(0)
        poly_list.append(poly_obj)
    #-- create shapely multipolygon object
    mpoly_obj = MultiPolygon(poly_list)
    #-- return the polygon object
    return mpoly_obj
