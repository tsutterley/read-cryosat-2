read-cryosat-2
==============

[![Language](https://img.shields.io/badge/python-v3.7-green.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/tsutterley/read-cryosat-2/blob/main/LICENSE)
[![Documentation Status](https://readthedocs.org/projects/read-cryosat-2/badge/?version=latest)](https://read-cryosat-2.readthedocs.io/en/latest/?badge=latest)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/tsutterley/read-cryosat-2/main)
[![Binder](https://binder.pangeo.io/badge.svg)](https://binder.pangeo.io/v2/gh/tsutterley/read-cryosat-2/main)

Python tools to read waveform and geolocated elevation data from the ESA CryoSat-2 mission and write to HDF5

- [How to access CryoSat-2 data](https://earth.esa.int/web/guest/-/how-to-access-cryosat-data-6842)
- [CryoSat-2 Products Overview](https://earth.esa.int/web/guest/-/products-overview-6975)
- [CryoSat-2 Products Handbook](https://earth.esa.int/documents/10174/125272/CryoSat_Product_Handbook)
- [CryoSat-2 Geographical Mode Masks](https://earth.esa.int/web/guest/-/geographical-mode-mask-7107)
- [CryoSat-2 Baseline C Improvements](https://earth.esa.int/documents/10174/1773005/C2-Evolution-BaselineC-Level2-V3)
- [CryoSat-2 Baseline D Improvements](https://earth.esa.int/documents/10174/1773005/CryoSat-Baseline-D-Evolutions.pdf)

#### Dependencies
- [numpy: Scientific Computing Tools For Python](https://numpy.org)
- [scipy: Scientific Tools for Python](https://docs.scipy.org/doc//)
- [h5py: Python interface for Hierarchal Data Format 5 (HDF5)](http://h5py.org)
- [netCDF4: Python interface to the netCDF C library](https://unidata.github.io/netcdf4-python/netCDF4/index.html)
- [gdal: Pythonic interface to the Geospatial Data Abstraction Library (GDAL)](https://pypi.python.org/pypi/GDAL)
- [shapely: PostGIS-ish operations outside a database context for Python](http://toblerity.org/shapely/index.html)
- [fiona: Python wrapper for vector data access functions from the OGR library](https://fiona.readthedocs.io/en/latest/manual.html)
- [geopandas: Python tools for geographic data](http://geopandas.readthedocs.io/)
- [pyproj: Python interface to PROJ library](https://pypi.org/project/pyproj/)
- [future: Compatibility layer between Python 2 and Python 3](http://python-future.org/)
- [lxml: processing XML and HTML in Python](https://pypi.python.org/pypi/lxml)

#### Download
The program homepage is:  
https://github.com/tsutterley/read-cryosat-2  
A zip archive of the latest version is available directly at:  
https://github.com/tsutterley/read-cryosat-2/archive/main.zip  
Incorporated into the UW-APL pointCollection repository at:  
https://github.com/SmithB/pointCollection  

#### Disclaimer
This program is not sponsored or maintained by the Universities Space Research Association (USRA), the European Space Agency (ESA) or NASA.  It is provided here for your convenience but _with no guarantees whatsoever_.
