======================
Setup and Installation
======================

Presently read-cryosat-2 is only available for use as a `GitHub repository <https://github.com/tsutterley/read-cryosat-2>`_.
The contents of the repository can be download as a `zipped file <https://github.com/tsutterley/read-cryosat-2/archive/main.zip>`_  or cloned.
To use this repository, please fork into your own account and then clone onto your system.

.. code-block:: bash

    git clone https://github.com/tsutterley/read-cryosat-2.git


Dependencies
############

read-cryosat-2 is dependent on some open source programs:

- `gdal <https://gdal.org/index.html>`_
- `pyproj <https://download.osgeo.org/proj>`_
- `HDF5 <https://www.hdfgroup.org>`_
- `netCDF <https://www.unidata.ucar.edu/software/netcdf>`_

The version of GDAL used within read-cryosat-2 will match the version of the installed C program.  The path to the C program that will be used with read-cryosat-2 is given by:

.. code-block:: bash

    gdal-config --datadir

The read-cryosat-2 installation uses the ``gdal-config`` routines to set the GDAL package version.

Installation
############

Can then install using ``setuptools``

.. code-block:: bash

    python3 setup.py install

or ``pip``

.. code-block:: bash

    python3 -m pip install --user .

Alternatively can install the cryosat_toolkit utilities from GitHub with ``pip``

.. code-block:: bash

    python3 -m pip install --user git+https://github.com/tsutterley/read-cryosat-2.git

Executable versions of this repository can also be tested using
`Binder <https://mybinder.org/v2/gh/tsutterley/read-cryosat-2/main>`_ and
`Pangeo <https://binder.pangeo.io/v2/gh/tsutterley/read-cryosat-2/main>`_.
