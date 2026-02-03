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

.. autoclass:: cryosat_toolkit.polygon
   :members:
