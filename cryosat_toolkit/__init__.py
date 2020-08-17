"""
A CryoSat-2 toolkit for Python
========================================
cryosat_toolkit contains Python tools for working with data from the ESA
CryoSat-2 mission

The package works using scientific Python packages (numpy and scipy)
combined with data storage in HDF5

It aims to be a simple and efficient solution for using data from the
CryoSat-2 mission and to support its science applications

Documentation is available at https://read-cryosat-2.readthedocs.io
"""
import cryosat_toolkit.utilities
from cryosat_toolkit.calc_GPS_time import calc_GPS_time
from cryosat_toolkit.read_cryosat_L1b import read_cryosat_L1b
from cryosat_toolkit.read_cryosat_L2 import read_cryosat_L2
from cryosat_toolkit.read_cryosat_L2I import read_cryosat_L2I
from cryosat_toolkit.HDF5_cryosat_L1b import HDF5_cryosat_L1b
from cryosat_toolkit.HDF5_cryosat_L2 import HDF5_cryosat_L2
from cryosat_toolkit.HDF5_cryosat_L2I import HDF5_cryosat_L2I
from cryosat_toolkit.read_shapefile import read_shapefile
from cryosat_toolkit.read_kml_file import read_kml_file
from cryosat_toolkit.read_geojson_file import read_geojson_file
