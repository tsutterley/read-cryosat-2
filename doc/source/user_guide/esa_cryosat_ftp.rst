==================
esa_cryosat_ftp.py
==================

- Syncs CryoSat-2 Elevation products from the ESA ftp data dissemination server
- Will sync all available CryoSat-2 data for a given product and set of years
- Can spatially subset using a bounding box or a georeferenced polygon file (shp, kml, kmz, GeoJSON)

`Source code`__

.. __: https://github.com/tsutterley/read-cryosat-2/blob/main/esa_cryosat_ftp.py

.. argparse::
    :filename: ../../scripts/esa_cryosat_ftp.py
    :func: arguments
    :prog: esa_cryosat_ftp.py
    :nodescription:
    :nodefault:
