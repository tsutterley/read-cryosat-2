esa_cryosat_ftp.py
==================

Syncs CryoSat-2 Elevation products from the ESA ftp data dissemination server.  
Will sync all available CryoSat-2 data for a given product and set of years  
Can spatially subset using a bounding box or a polygon file (shp, kml, kmz, GeoJSON)  

#### Calling Sequence
```bash
python esa_cryosat_ftp.py --baseline=D --user=<username> --year=2011 SIR_SIN_L2
```
where `<username>` is your ESA data dissemination server username  
your ESA data dissemination password will be entered from the command-line within the program  
[Source code](https://github.com/tsutterley/read-cryosat-2/blob/master/esa_cryosat_ftp.py)  

#### Inputs
 1. CryoSat-2 product to sync with ESA servers

#### Command Line Options
 - `-h`, `--help`: list the command line options
 - `-U X`, `--user=X`: Username for ESA ftp data dissemination
 - `-D X`, `--directory=X`: Working Data Directory
 - `-Y X`, `--year=X`: years to sync separated by commas
 - `-B X`, `--baseline=X`: CryoSat-2 baseline to sync
 - `--bbox=X`: Bounding box (lonmin,latmin,lonmax,latmax)
 - `--polygon=X`: Georeferenced file containing a set of polygons
 - `-M X`, `--mode=X`: Permission mode of directories and files synced
 - `-l`, `--log`: Output log file
 - `-L`, `--list`: Only print files that are to be transferred
 - `--clobber`: Overwrite existing data in transfer
