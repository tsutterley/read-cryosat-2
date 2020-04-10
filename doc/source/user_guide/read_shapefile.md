read_shapefile.py
=================

 - Reads polygons from ESRI shapefiles

[Source code](https://github.com/tsutterley/read-cryosat-2/blob/master/cryosat_toolkit/read_shapefile.py)  

#### Inputs
```
 1. input shapefile (.shp)
```

#### Options
 - `EPSG`: projection identifier for output coordinates
 - `ZIP`: input file is compressed
 - `VARIABLES`: reduce to a specific set of identifiers

#### Outputs
 - `poly_obj`: shapely multipolygon object of input file
