read_kml_file.py
================

 - Reads polygons from keyhole markup language (.kml or .kmz) files

[Source code](https://github.com/tsutterley/read-cryosat-2/blob/master/cryosat_toolkit/read_kml_file.py)  

#### Inputs
```
 1. input .kml or .kmz file
```

#### Options
 - `EPSG`: projection identifier for output coordinates
 - `KMZ`: input file is compressed
 - `VARIABLES`: reduce to a specific set of identifiers

#### Outputs
 - `poly_obj`: shapely multipolygon object of input file
