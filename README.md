read-cryosat-2
==============

Reads waveform and geolocated elevation data from the ESA CryoSat-2 mission  

- [How to access CryoSat-2 data](https://earth.esa.int/web/guest/-/how-to-access-cryosat-data-6842)  
- [CryoSat-2 Products Overview](https://earth.esa.int/web/guest/-/products-overview-6975)  
- [CryoSat-2 Products Handbook](https://earth.esa.int/documents/10174/125272/CryoSat_Product_Handbook)  
- [CryoSat-2 Geographical Mode Masks](https://earth.esa.int/web/guest/-/geographical-mode-mask-7107)  
- [CryoSat-2 Baseline C Improvements](https://earth.esa.int/documents/10174/1773005/C2-Evolution-BaselineC-Level2-V3)  

#### `esa_cryosat_sync.py`
Program to sync Cryosat-2 Elevation products from the ESA data dissemination server.  
```bash
python esa_cryosat_sync.py --baseline=C --user=<username> SIR_SIN_L2
```
where `<username>` is your ESA data dissemination server username  

##### `read_cryosat_L1b.py`
Program to read Level-1B CryoSat waveform data into a python environment.  
The Level-1B waveform data come as packed [Payload Data System (PDS) files](https://earth.esa.int/documents/10174/125273/CryoSat_L1_Products_Format_Specification).  
Units for each variable will match the original files.  
```python
import read_cryosat.read_cryosat_L1b as read_cryosat_L1b
CS_L1b_mds = read_cryosat_L1b(full_filename)
```

##### `read_cryosat_L2.py`
Program to read Level-2 CryoSat elevation data into a python environment.  
The Level-2 data come as packed [Payload Data System (PDS) files](https://earth.esa.int/documents/10174/125273/CryoSat_L2_Products_Format_Specification).  
Units for each variable will match the original files.  
```python
import read_cryosat.read_cryosat_L2 as read_cryosat_L2
CS_L2_mds = read_cryosat_L2(full_filename)
```

##### `read_cryosat_L2I.py`
Program to read Level-2 Intermediate CryoSat data into a python environment.  
The Level-2I data come as packed Payload Data System (PDS) files.  
Units for each variable will match the original files.  
```python
import read_cryosat.read_cryosat_L2I as read_cryosat_L2I
CS_L2I_mds = read_cryosat_L2I(full_filename)
```

#### Dependencies
 - [numpy: Scientific Computing Tools For Python](http://www.numpy.org)  

#### Download
The program homepage is:   
https://github.com/tsutterley/read-cryosat-2   
A zip archive of the latest version is available directly at:    
https://github.com/tsutterley/read-cryosat-2/archive/master.zip  

#### Disclaimer  
This program is not sponsored or maintained by the Universities Space Research Association (USRA), the European Space Agency (ESA) or NASA.  It is provided here for your convenience but _with no guarantees whatsoever_.  
