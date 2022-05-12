===================
esa_cryosat_sync.py
===================

- Syncs Cryosat-2 Elevation products from the ESA https Science Server
- Will sync all available CryoSat-2 data for a given product and set of years
- Can spatially subset using a bounding box or a georeferenced polygon file (shp, kml, kmz, GeoJSON)

`Source code`__

.. __: https://github.com/tsutterley/read-cryosat-2/blob/main/esa_cryosat_sync.py

.. argparse::
    :filename: ../../scripts/esa_cryosat_sync.py
    :func: arguments
    :prog: esa_cryosat_sync.py
    :nodescription:
    :nodefault:
