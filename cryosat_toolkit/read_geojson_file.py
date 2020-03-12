#!/usr/bin/env python
u"""
read_geojson_file.py
Written by Tyler Sutterley (09/2019)
Reads polygons from GeoJSON files

INPUTS:
    input GeoJSON file (.json, .geojson)

OUTPUT:
    shapely multipolygon object of input file

OPTIONS:
    EPSG: projection identifier for output coordinates
    VARIABLES: reduce to a specific set of identifiers

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        http://www.numpy.org
        http://www.scipy.org/NumPy_for_Matlab_Users
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
    Updated 06/2019: using geopandas for consistency between read functions
        convert projection to EPSG:4326 before creating polygons
    Written 06/2019
"""
from __future__ import print_function

import os
import pyproj
import geopandas
import numpy as np
from shapely.geometry import Polygon, MultiPolygon

#-- PURPOSE: read GeoJSON (.json, .geojson) files
def read_geojson_file(input_file, EPSG=4326, VARIABLES=None):
    #-- read the GeoJSON file
    gj = geopandas.read_file(input_file)

    #-- convert projection to EPSG
    proj1 = pyproj.Proj("+init={0}".format(gj.crs['init']))
    proj2 = pyproj.Proj("+init=EPSG:{0:d}".format(EPSG))

    #-- list of polygons
    poly_list = []
    #-- find features of interest
    geometries = ('LineString','Polygon')
    f = [f for f in gj.iterfeatures() if f['geometry']['type'] in geometries]
    #-- reduce to variables of interest if specified
    f = [ft for ft in f if ft['id'] in VARIABLES] if VARIABLES else f

    #-- for each line string or polygon feature
    for feature in f:
        #-- extract coordinates for feature
        x,y = np.transpose(feature['geometry']['coordinates'])
        #-- convert points to latitude/longitude
        lon,lat = pyproj.transform(proj1, proj2, x, y)
        #-- create shapely polygon
        poly_obj = Polygon(list(zip(lon, lat)))
        #-- Valid Polygon cannot have overlapping exterior or interior rings
        if (not poly_obj.is_valid):
            poly_obj = poly_obj.buffer(0)
        poly_list.append(poly_obj)
    #-- create shapely multipolygon object
    mpoly_obj = MultiPolygon(poly_list)
    #-- return the polygon object
    return mpoly_obj
