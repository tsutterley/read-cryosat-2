==========
polygon.py
==========

Reads polygons from GeoJSON, kml/kmz or ESRI shapefile files

Calling Sequence
================

Reading an ESRI shapefile containing polygons

.. code-block:: python

    from cryosat_toolkit.polygon import polygon
    mpoly_obj = polygon().from_shapefile(path_to_esri_shapefile)

Reading a kmz file containing polygons

.. code-block:: python

    from cryosat_toolkit.polygon import polygon
    mpoly_obj = polygon().from_kml(path_to_kmz_file,kmz=True)


`Source code`__

.. __: https://github.com/tsutterley/read-cryosat-2/blob/main/cryosat_toolkit/polygon.py


General Attributes and Methods
==============================

.. class:: polygon(object)


    .. attribute:: object.filename

        path of input georeferenced file


    .. attribute:: object.epsg

        spatial projection identifier for output coordinates


    .. method:: object.case_insensitive_filename(filename)

        Searches a directory for a filename without case dependence


    .. method:: object.from_geojson(filename, variables=None)

        Reads polygons from GeoJSON files

        Arguments:

            full path of input GeoJSON file (.json, .geojson)

        Keyword arguments:

            ``variables``: reduce to a specific set of identifiers

        Returns:

            ``mpoly_obj``: shapely multipolygon object


    .. method:: object.from_kml(filename, kmz=False, variables=None)

        Reads polygons from keyhole markup language files

        Arguments:

            full path of input markup file (.kml, .kmz)

        Keyword arguments:

            ``kmz``: input file is compressed

            ``variables``: reduce to a specific set of identifiers

        Returns:

            ``mpoly_obj``: shapely multipolygon object


    .. method:: object.from_shapefile(filename, zip=False, variables=None)

        read ESRI shapefiles

        Arguments:

            full path of input shapefile (*.shp)

        Keyword arguments:

            ``zip`` input file is compressed

            ``variables``: reduce to a specific set of identifiers

        Returns:

            ``mpoly_obj``: shapely multipolygon object
