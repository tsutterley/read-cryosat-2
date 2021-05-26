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
[Source code](https://github.com/tsutterley/read-cryosat-2/blob/main/cryosat_toolkit/read_cryosat_L1b.py)

#### Arguments
- `full_filename`: full path of CryoSat .DBL or .nc file

#### Keyword arguments
- `CS_L1b_mds`: Python dictionary with groups:
    * `'Location'`: Time and Orbit Group
    * `'Data'`: Measurements Group
    * `'Geometry'`: External Corrections Group
    * `'Waveform_1Hz'`: Average Waveforms Group
    * `'Waveform_20Hz'`: Waveforms Group (with SAR/SARIN Beam Behavior Parameters)
    * `'METADATA'`: Header data
        * `'MPH'`: Main Product Headers
        * `'SPH'`: Specific Product Headers
        * `'DSD'`: Data Set Descriptors

#### Level-1 File Convention
The file names use the following convention:
- `MM_CCCC_XXXXXXXXXX_yyyymmdd_hhmmss_YYYYMMDD_HHMMSS__bvvv.ttt`
- `MM` is the Mission ID (`CS` for CryoSat-2)
- `CCCC` is the file class (see File Class table)
- `XXXXXXXXXX` is the file type (see File Type table)
- `yyyymmdd_hhmmss` is the start date and time in UTC
- `YYYYMMDD_HHMMSS` is the stop date and time in UTC
- `b` is the baseline identifier
- `vvv` is the version number of the file
- `ttt` is the extension
    * `HDR`: header file
    * `DBL`: binary data file
    * `nc`: netCDF data file

File Class | Description
:--------: | -----------
`OFFL` | Off Line Processing/Systematic
`NRT_` | Near Real Time
`RPRO` | Reprocessing
`TEST` | Testing
`TIxx` | Stand alone IPF1 testing
`LTA_` | Long Term Archive

File Type  | Description
:--------: | -----------
`SIR1SAR_FR` | Level 1 FBR SAR Mode (Rx1 Channel)
`SIR2SAR_FR` | Level 1 FBR SAR Mode (Rx2 Channel)
`SIR_SIN_FR` | Level 1 FBR SARin Mode
`SIR_LRM_1B` | Level-1 Product Low Rate Mode
`SIR_FDM_1B` | Level-1 Product Fast Delivery Marine Mode
`SIR_SAR_1B` | Level-1 SAR Mode
`SIR_SIN_1B` | Level-1 SARin Mode
`SIR1LRC11B` | Level-1 CAL1 Low Rate Mode (Rx1 Channel)
`SIR2LRC11B` | Level-1 CAL1 Low Rate Mode (Rx2 Channel)
`SIR1SAC11B` | Level-1 CAL1 SAR Mode (Rx1 Channel)
`SIR2SAC11B` | Level-1 CAL1 SAR Mode (Rx2 Channel)
`SIR_SIC11B` | Level-1 CAL1 SARin Mode
`SIR_SICC1B` | Level-1 CAL1 SARIN Exotic Data
`SIR1SAC21B` | Level-1 CAL2 SAR Mode (Rx1 Channel)
`SIR2SAC21B` | Level-1 CAL2 SAR Mode (Rx2 Channel)
`SIR1SIC21B` | Level-1 CAL2 SARin Mode (Rx1 Channel)
`SIR2SIC21B` | Level-1 CAL2 SARin Mode (Rx1 Channel)
`SIR1LRM_0M` | LRM and TRK Monitoring Data from Rx 1 Channel
`SIR2LRM_0M` | LRM and TRK Monitoring Data from Rx 2 Channel
`SIR1SAR_0M` | SAR Monitoring Data from Rx 1 Channel
`SIR2SAR_0M` | SAR Monitoring Data from Rx 1 Channel
`SIR_SIN_0M` | SARIN Monitoring Data
`SIR_SIC40M` | CAL4 Monitoring Data
