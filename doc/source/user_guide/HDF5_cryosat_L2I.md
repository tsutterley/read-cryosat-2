HDF5_cryosat_L2I.py
===================

- Writes Level-2 Intermediate CryoSat data into HDF5 files
- Supported CryoSat Modes: LRM, SAR, SARin, FDM, SID, GDR
- Units for each variable should match the original SDS files

#### Calling Sequence
```python
import cryosat_toolkit.read_cryosat_L2I as read_cryosat_L2I
import cryosat_toolkit.HDF5_cryosat_L2I as HDF5_cryosat_L2I
MODE = 'SIN'
BASELINE = 'C'
CS_L2I_mds = read_cryosat_L2I(full_filename)
HDF5_cryosat_L2I(CS_L2I_mds, MODE, BASELINE, FILENAME=full_HDF5_filename)
```
[Source code](https://github.com/tsutterley/read-cryosat-2/blob/main/cryosat_toolkit/HDF5_cryosat_L2I.py)

#### Inputs
 1. `CS_L2I_mds`: Python dictionary with groups
     * 'Location': Time and Orbit Parameters
     * 'Geometry': Elevation Corrections and Flags
     * 'Data': Geolocation and Elevation Measurements with Quality Parameters
     * 'Auxiliary': Auxiliary Data for Elevation Processing
     * 'Instrumental': Intrument Corrections
     * 'METADATA': MPH, SPH and DSD Header data
 2. `MODE`: CryoSat-2 Modes  (LRM, SAR, SARin, FDM, SID, GDR)
 3. `BASELINE`: CryoSat-2 baseline (A, B, C, D)

#### Options
 - `FILENAME`: output HDF5 file name
 - `TITLE`: output file description
 - `HEADER`: output CryoSat-2 file headers (MPH, SPH, DSD)
 - `CLOBBER`: overwrite existing HDF5 file
 - `VERBOSE`: print HDF5 structure parameters to screen
