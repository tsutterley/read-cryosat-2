read_cryosat_L2I.py
===================

 - Reads Level-2 Intermediate CryoSat data into a python environment.
 - Level-2I data for Baselines A, B and C come as packed [Payload Data System (PDS) files](https://earth.esa.int/documents/10174/125273/CryoSat_L2_Products_Format_Specification).
 - Level-2I data for Baseline D comes as [netCDF4 files](https://earth.esa.int/documents/10174/125272/CryoSat-Baseline-D-Product-Handbook).
 - Supported CryoSat Modes: LRM, SAR, SARin, FDM, SID, GDR

#### Calling Sequence
```python
import cryosat_toolkit.read_cryosat_L2I as read_cryosat_L2I
CS_L2I_mds = read_cryosat_L2I(full_filename)
```
[Source code](https://github.com/tsutterley/read-cryosat-2/blob/main/cryosat_toolkit/read_cryosat_L2I.py)

#### Arguments
 - `full_filename`: full path of CryoSat .DBL or .nc file

#### Keyword arguments
 - `CS_L2I_mds`: Python dictionary with groups:
     * 'Location': Time and Orbit Parameters
     * 'Geometry': Elevation Corrections and Flags
     * 'Data': Geolocation and Elevation Measurements with Quality Parameters
     * 'Auxiliary': Auxiliary Data for Elevation Processing
     * 'Instrumental': Intrument Corrections
     * 'METADATA': MPH, SPH and DSD Header data
