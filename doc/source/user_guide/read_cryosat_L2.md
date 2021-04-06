read_cryosat_L2.py
==================

 - Reads Level-2 CryoSat elevation data into a python environment.
 - Level-2 data for Baselines A, B and C come as packed [Payload Data System (PDS) files](https://earth.esa.int/documents/10174/125273/CryoSat_L2_Products_Format_Specification).
 - Level-2 data for Baseline D comes as [netCDF4 files](https://earth.esa.int/documents/10174/125272/CryoSat-Baseline-D-Product-Handbook).
 - Supported CryoSat Modes: LRM, SAR, SARin, FDM, SID, GDR

#### Calling Sequence
```python
import cryosat_toolkit.read_cryosat_L2 as read_cryosat_L2
CS_L2_mds = read_cryosat_L2(full_filename)
```
[Source code](https://github.com/tsutterley/read-cryosat-2/blob/main/cryosat_toolkit/read_cryosat_L2.py)

#### Arguments
 - `full_filename`: full path of CryoSat .DBL or .nc file

#### Keyword arguments
 - `CS_L2_mds`: Python dictionary with groups:
     * 'Data_1Hz': Time and Orbit Parameters
     * 'Corrections': Elevation Corrections and Flags
     * 'Data_20Hz': Geolocation and Elevation Measurements with Quality Parameters
     * 'METADATA': MPH, SPH and DSD Header data
