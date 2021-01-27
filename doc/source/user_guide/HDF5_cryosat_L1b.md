HDF5_cryosat_L1b.py
===================

 - Writes Level-1B CryoSat waveform data into HDF5 files
 - Supported CryoSat Modes: LRM, SAR, SARin, FDM, SID, GDR
 - Units for each variable should match the original SDS files

#### Calling Sequence
```python
import cryosat_toolkit.read_cryosat_L1b as read_cryosat_L1b
import cryosat_toolkit.HDF5_cryosat_L1b as HDF5_cryosat_L1b
MODE = 'SIN'
BASELINE = 'C'
CS_L1b_mds = read_cryosat_L1b(full_filename)
HDF5_cryosat_L1b(CS_L1b_mds, MODE, BASELINE, FILENAME=full_HDF5_filename)
```
[Source code](https://github.com/tsutterley/read-cryosat-2/blob/main/cryosat_toolkit/HDF5_cryosat_L1b.py)

#### Inputs
 1. `CS_L1b_mds`: Python dictionary with groups
     * 'Location': Time and Orbit Group
     * 'Data': Measurements Group
     * 'Geometry': External Corrections Group
     * 'Waveform_1Hz': Average Waveforms Group
     * 'Waveform_20Hz': Waveforms Group (with SAR/SARIN Beam Behavior Parameters)
     * 'METADATA': MPH, SPH and DSD Header data
 2. `MODE`: CryoSat-2 Modes  (LRM, SAR, SARin, FDM, SID, GDR)
 3. `BASELINE`: CryoSat-2 baseline (A, B, C, D)

#### Options
 - `FILENAME`: output HDF5 file name
 - `TITLE`: output file description
 - `HEADER`: output CryoSat-2 file headers (MPH, SPH, DSD)
 - `CLOBBER`: overwrite existing HDF5 file
 - `VERBOSE`: print HDF5 structure parameters to screen
