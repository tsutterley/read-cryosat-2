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
    * `'Location'`: Time and Orbit Parameters
    * `'Geometry'`: Elevation Corrections and Flags
    * `'Data'`: Geolocation and Elevation Measurements with Quality Parameters
    * `'Auxiliary'`: Auxiliary Data for Elevation Processing
    * `'Instrumental'`: Intrument Corrections
    * `'METADATA'`: Header data
        * `'MPH'`: Main Product Headers
        * `'SPH'`: Specific Product Headers
        * `'DSD'`: Data Set Descriptors

#### Level-2 Intermediate File Convention
The file names use the following convention:
`MM_CCCC_XXXXXXXXXX_yyyymmdd_hhmmss_YYYYMMDD_HHMMSS__bvvv.ttt`
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
`LTA_` | Long Term Archive

File Type  | Description
:--------: | -----------
`SIR_LRMI2_` | Intermediate L2 Product from LRM Processing
`SIR_SINI2_` | Intermediate L2 Product from SIN Processing
`SIR_SIDI2_` | Intermediate L2 Product from SIN Degraded Process
`SIR_SARI2A` | Intermediate L2 Product from SAR Step 1 Processing
`SIR_SARI2B` | Intermediate L2 Product from SAR Step 2 Processing
