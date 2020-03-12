#!/usr/bin/env python
u"""
read_shapefile.py
Written by Tyler Sutterley (09/2019)
Reads polygons from ESRI shapefiles

INPUTS:
    input shapefile (.shp)

OUTPUT:
    shapely multipolygon object of input file

OPTIONS:
    EPSG: projection identifier for output coordinates
    ZIP: input file is compressed
    VARIABLES: reduce to a specific set of identifiers

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        http://www.numpy.org
        http://www.scipy.org/NumPy_for_Matlab_Users
    fiona: Python wrapper for vector data access functions from the OGR library
        https://fiona.readthedocs.io/en/latest/manual.html
    shapely: PostGIS-ish operations outside a database context for Python
        http://toblerity.org/shapely/index.html
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/

UPDATE HISTORY:
    Updated 09/2019: made output coordinate system an option (EPSG)
    Updated 07/2019: added option to reduce to specific VARIABLES within file
    Updated 06/2019: using fiona for consistency between read functions
        convert projection to EPSG:4326 before creating polygons
    Written 06/2019
"""
from __future__ import print_function

import os
import fiona
import numpy as np
import pyproj
from shapely.geometry import Polygon, MultiPolygon

#-- PURPOSE: read shapefiles
def read_shapefile(input_file, EPSG=4326, ZIP=False, VARIABLES=None):
    #-- read input zipfile containing shapefiles
    if ZIP:
        #-- read the compressed shapefile and extract entities
        shape = fiona.open('zip://{0}'.format(os.path.expanduser(input_file)))
    else:
        #-- read the shapefile and extract entities
        shape = fiona.open(os.path.expanduser(input_file),'r')

    #-- convert projection to EPSG
    proj1 = pyproj.Proj("+init={0}".format(shape.crs['init']))
    proj2 = pyproj.Proj("+init=EPSG:{0:d}".format(EPSG))

    #-- find features of interest
    geometries = ('LineString','Polygon')
    f = [f for f in shape.values() if f['geometry']['type'] in geometries]
    #-- reduce to variables of interest if specified
    f = [ft for ft in f if ft['id'] in VARIABLES] if VARIABLES else f

    #-- list of polygons
    poly_list = []
    #-- for each entity
    for i,ent in enumerate(f):
        #-- extract coordinates for entity
        for coords in ent['geometry']['coordinates']:
            #-- convert points to latitude/longitude
            x,y = np.transpose(coords)
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
