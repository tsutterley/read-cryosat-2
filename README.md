read-cryosat-2
==============

Tools to read waveform and geolocated elevation data from the ESA CryoSat-2 mission and write to HDF5

- [How to access CryoSat-2 data](https://earth.esa.int/web/guest/-/how-to-access-cryosat-data-6842)  
- [CryoSat-2 Products Overview](https://earth.esa.int/web/guest/-/products-overview-6975)  
- [CryoSat-2 Products Handbook](https://earth.esa.int/documents/10174/125272/CryoSat_Product_Handbook)  
- [CryoSat-2 Geographical Mode Masks](https://earth.esa.int/web/guest/-/geographical-mode-mask-7107)  
- [CryoSat-2 Baseline C Improvements](https://earth.esa.int/documents/10174/1773005/C2-Evolution-BaselineC-Level2-V3)  

#### `esa_cryosat_sync.py`
Program to sync Cryosat-2 Elevation products from the ESA data dissemination server.  
Will sync all available CryoSat-2 data for a given product and set of years  
```bash
python esa_cryosat_sync.py --baseline=C --user=<username> --year=2010,2011 SIR_SIN_L2
```
where `<username>` is your ESA data dissemination server username  
your ESA data dissemination password will be entered from the command-line within the program  

##### `read_cryosat_L1b.py`
Program to read Level-1B CryoSat waveform data into a python environment.  
The Level-1B waveform data come as packed [Payload Data System (PDS) files](https://earth.esa.int/documents/10174/125273/CryoSat_L1_Products_Format_Specification).  
Units for each variable will match the original files.  
```python
import cryosat_toolkit.read_cryosat_L1b as read_cryosat_L1b
CS_L1b_mds = read_cryosat_L1b(full_filename)
```

##### `HDF5_cryosat_L1b.py`
Program to write Level-1B CryoSat waveform data into HDF5 files.  
Units for each variable should match the original SDS files.  
```python
import cryosat_toolkit.read_cryosat_L1b as read_cryosat_L1b
import cryosat_toolkit.HDF5_cryosat_L1b as HDF5_cryosat_L1b
MODE = 'SIN'
BASELINE = 'C'
CS_L1b_mds = read_cryosat_L1b(full_filename)
HDF5_cryosat_L1b(CS_L1b_mds, MODE, BASELINE, FILENAME=full_HDF5_filename)
```

##### `read_cryosat_L2.py`
Program to read Level-2 CryoSat elevation data into a python environment.  
The Level-2 data come as packed [Payload Data System (PDS) files](https://earth.esa.int/documents/10174/125273/CryoSat_L2_Products_Format_Specification).  
Units for each variable will match the original files.  
```python
import cryosat_toolkit.read_cryosat_L2 as read_cryosat_L2
CS_L2_mds = read_cryosat_L2(full_filename)
```

##### `HDF5_cryosat_L2.py`
Program to write Level-2 CryoSat elevation data into HDF5 files.  
Units for each variable should match the original SDS files.  
```python
import cryosat_toolkit.read_cryosat_L2 as read_cryosat_L2
import cryosat_toolkit.HDF5_cryosat_L2 as HDF5_cryosat_L2
BASELINE = 'C'
CS_L2_mds = read_cryosat_L2(full_filename)
HDF5_cryosat_L2(CS_L2_mds, BASELINE, FILENAME=full_HDF5_filename)
```

##### `read_cryosat_L2I.py`
Program to read Level-2 Intermediate CryoSat data into a python environment.  
The Level-2I data come as packed Payload Data System (PDS) files.  
Units for each variable will match the original files.  
```python
import cryosat_toolkit.read_cryosat_L2I as read_cryosat_L2I
CS_L2I_mds = read_cryosat_L2I(full_filename)
```

##### `HDF5_cryosat_L2I.py`
Program to write Level-2 Intermediate CryoSat data into HDF5 files.  
Units for each variable should match the original SDS files.  
```python
import cryosat_toolkit.read_cryosat_L2I as read_cryosat_L2I
import cryosat_toolkit.HDF5_cryosat_L2I as HDF5_cryosat_L2I
MODE = 'SIN'
BASELINE = 'C'
CS_L2I_mds = read_cryosat_L2I(full_filename)
read_cryosat_L2I(CS_L2I_mds, MODE, BASELINE, FILENAME=full_HDF5_filename)
```

#### Dependencies
 - [numpy: Scientific Computing Tools For Python](http://www.numpy.org)  
 - [h5py: Python interface for Hierarchal Data Format 5 (HDF5)](http://h5py.org)  

#### Download
The program homepage is:   
https://github.com/tsutterley/read-cryosat-2   
A zip archive of the latest version is available directly at:    
https://github.com/tsutterley/read-cryosat-2/archive/master.zip  

#### Disclaimer  
This program is not sponsored or maintained by the Universities Space Research Association (USRA), the European Space Agency (ESA) or NASA.  It is provided here for your convenience but _with no guarantees whatsoever_.  
