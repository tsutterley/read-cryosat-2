HDF5_cryosat_L2.py
==================

- Writes Level-2 CryoSat elevation data into HDF5 files
- Supported CryoSat Modes: LRM, SAR, SARin, FDM, SID, GDR
- Units for each variable should match the original SDS files

#### Calling Sequence
```python
import cryosat_toolkit.read_cryosat_L2 as read_cryosat_L2
import cryosat_toolkit.HDF5_cryosat_L2 as HDF5_cryosat_L2
BASELINE = 'C'
CS_L2_mds = read_cryosat_L2(full_filename)
HDF5_cryosat_L2(CS_L2_mds, BASELINE, FILENAME=full_HDF5_filename)
```
[Source code](https://github.com/tsutterley/read-cryosat-2/blob/main/cryosat_toolkit/HDF5_cryosat_L2.py)

#### Inputs
 1. `CS_L2_mds`: Python dictionary with groups
     * 'Data_1Hz': Time and Orbit Parameters
     * 'Corrections': Elevation Corrections and Flags
     * 'Data_20Hz': Geolocation and Elevation Measurements with Quality Parameters
     * 'METADATA': MPH, SPH and DSD Header data
 2. `BASELINE`: CryoSat-2 baseline (A, B, C, D)

#### Options
 - `FILENAME`: output HDF5 file name
 - `TITLE`: output file description
 - `HEADER`: output CryoSat-2 file headers (MPH, SPH, DSD)
 - `CLOBBER`: overwrite existing HDF5 file
 - `VERBOSE`: print HDF5 structure parameters to screen
