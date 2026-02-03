==================
read_cryosat_L2.py
==================

- Reads Level-2 CryoSat elevation data into a python environment.
- Level-2 data for Baselines A, B and C come as packed `Payload Data System (PDS) files <https://earth.esa.int/documents/10174/125273/CryoSat_L2_Products_Format_Specification>`_.
- Level-2 data for Baseline D comes as `netCDF4 files <https://earth.esa.int/documents/10174/125272/CryoSat-Baseline-D-Product-Handbook>`_.
- Supported CryoSat Modes: LRM, SAR, SARin, FDM, SID, GDR

`Source code`__

.. __: https://github.com/tsutterley/read-cryosat-2/blob/main/cryosat_toolkit/read_cryosat_L2.py

.. autofunction:: cryosat_toolkit.read_cryosat_L2

.. autofunction:: cryosat_toolkit.read_cryosat_L2.cryosat_baseline_AB

.. autofunction:: cryosat_toolkit.read_cryosat_L2.cryosat_baseline_C

.. autofunction:: cryosat_toolkit.read_cryosat_L2.cryosat_baseline_D

.. autofunction:: cryosat_toolkit.read_cryosat_L2.cryosat_scaling_factors

.. autofunction:: cryosat_toolkit.read_cryosat_L2.read_MPH

.. autofunction:: cryosat_toolkit.read_cryosat_L2.read_SPH

.. autofunction:: cryosat_toolkit.read_cryosat_L2.read_DSD

File Convention
###############

The file names use the following convention:
``MM_CCCC_XXXXXXXXXX_yyyymmdd_hhmmss_YYYYMMDD_HHMMSS__bvvv.ttt``

- ``MM`` is the Mission ID (``CS`` for CryoSat-2)
- ``CCCC`` is the file class (see File Class table)
- ``XXXXXXXXXX`` is the file type (see File Type table)
- ``yyyymmdd_hhmmss`` is the start date and time in UTC
- ``YYYYMMDD_HHMMSS`` is the stop date and time in UTC
- ``b`` is the baseline identifier
- ``vvv`` is the version number of the file
- ``ttt`` is the extension
    * ``HDR``: header file
    * ``DBL``: binary data file
    * ``nc``: netCDF data file

File Class | Description
:--------: | -----------
``OFFL`` | Off Line Processing/Systematic
``NRT_`` | Near Real Time
``RPRO`` | Reprocessing
``TEST`` | Testing
``LTA_`` | Long Term Archive

File Type  | Description
:--------: | -----------
``SIR_LRMI2_`` | Intermediate L2 Product from LRM Processing
``SIR_SINI2_`` | Intermediate L2 Product from SIN Processing
``SIR_SIDI2_`` | Intermediate L2 Product from SIN Degraded Process
``SIR_SARI2A`` | Intermediate L2 Product from SAR Step 1 Processing
``SIR_SARI2B`` | Intermediate L2 Product from SAR Step 2 Processing
