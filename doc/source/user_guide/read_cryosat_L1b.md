read_cryosat_L1b.py
===================

 - Reads Level-1B CryoSat waveform data into a python environment.  
 - Level-1B waveform data for Baselines A, B and C come as packed [Payload Data System (PDS) files](https://earth.esa.int/documents/10174/125273/CryoSat_L1_Products_Format_Specification).  
 - Level-1B data for Baseline D comes as [netCDF4 files](https://earth.esa.int/documents/10174/125272/CryoSat-Baseline-D-Product-Handbook).  
 - Supported CryoSat Modes: LRM, SAR, SARin, FDM, SID, GDR

#### Calling Sequence
```python
import cryosat_toolkit.read_cryosat_L1b as read_cryosat_L1b
CS_L1b_mds = read_cryosat_L1b(full_filename)
```

#### Inputs
 - `full_filename`: full path of CryoSat .DBL or .nc file

#### Options
 - `CS_L1b_mds`: Python dictionary with groups:
     * 'Location': Time and Orbit Group
     * 'Data': Measurements Group
     * 'Geometry': External Corrections Group
     * 'Waveform_1Hz': Average Waveforms Group
     * 'Waveform_20Hz': Waveforms Group (with SAR/SARIN Beam Behavior Parameters)
     * 'METADATA': MPH, SPH and DSD Header data
