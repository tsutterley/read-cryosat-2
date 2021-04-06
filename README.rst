==============
read-cryosat-2
==============

|Language|
|License|
|Documentation Status|
|Binder|
|Pangeo|

.. |Language| image:: https://img.shields.io/badge/python-v3.7-green.svg
   :target: https://www.python.org/

.. |License| image:: https://img.shields.io/badge/license-MIT-green.svg
   :target: https://github.com/tsutterley/read-cryosat-2/blob/main/LICENSE

.. |Documentation Status| image:: https://readthedocs.org/projects/read-cryosat-2/badge/?version=latest
   :target: https://read-cryosat-2.readthedocs.io/en/latest/?badge=latest

.. |Binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/tsutterley/read-cryosat-2/main

.. |Pangeo| image:: https://binder.pangeo.io/badge.svg
   :target: https://binder.pangeo.io/v2/gh/tsutterley/read-cryosat-2/main

Python tools for reading waveform and geolocated elevation data from the ESA CryoSat-2 mission and writing to HDF5

- `How to access CryoSat-2 data <https://earth.esa.int/web/guest/-/how-to-access-cryosat-data-6842>`_
- `CryoSat-2 Products Overview <https://earth.esa.int/web/guest/-/products-overview-6975>`_
- `CryoSat-2 Products Handbook <https://earth.esa.int/documents/10174/125272/CryoSat_Product_Handbook>`_
- `CryoSat-2 Geographical Mode Masks <https://earth.esa.int/web/guest/-/geographical-mode-mask-7107>`_
- `CryoSat-2 Baseline C Improvements <https://earth.esa.int/documents/10174/1773005/C2-Evolution-BaselineC-Level2-V3>`_
- `CryoSat-2 Baseline D Improvements <https://earth.esa.int/documents/10174/1773005/CryoSat-Baseline-D-Evolutions.pdf>`_

Dependencies
############

- `numpy: Scientific Computing Tools For Python <https://numpy.org>`_
- `scipy: Scientific Tools for Python <https://docs.scipy.org/doc//>`_
- `h5py: Python interface for Hierarchal Data Format 5 (HDF5) <http://h5py.org>`_
- `netCDF4: Python interface to the netCDF C library <https://unidata.github.io/netcdf4-python/netCDF4/index.html>`_
- `gdal: Pythonic interface to the Geospatial Data Abstraction Library (GDAL) <https://pypi.python.org/pypi/GDAL>`_
- `shapely: PostGIS-ish operations outside a database context for Python <http://toblerity.org/shapely/index.html>`_
- `fiona: Python wrapper for vector data access functions from the OGR library <https://fiona.readthedocs.io/en/latest/manual.html>`_
- `geopandas: Python tools for geographic data <http://geopandas.readthedocs.io/>`_
- `pyproj: Python interface to PROJ library <https://pypi.org/project/pyproj/>`_
- `future: Compatibility layer between Python 2 and Python 3 <http://python-future.org/>`_
- `lxml: processing XML and HTML in Python <https://pypi.python.org/pypi/lxml>`_

Download
########

| The program homepage is:
| https://github.com/tsutterley/read-cryosat-2
| A zip archive of the latest version is available directly at:
| https://github.com/tsutterley/read-cryosat-2/archive/main.zip
| Incorporated into the UW-APL ``pointCollection`` utilities at:
| https://github.com/SmithB/pointCollection

Disclaimer
##########

This project contains work and contributions from the `scientific community <./CONTRIBUTORS.rst>`_.
This program is not sponsored or maintained by the Universities Space Research Association (USRA), the European Space Agency (ESA) or NASA.
It is provided here for your convenience but *with no guarantees whatsoever*.
