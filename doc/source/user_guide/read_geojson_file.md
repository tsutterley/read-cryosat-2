read_geojson_file.py
====================

 - Reads polygons from GeoJSON files

[Source code](https://github.com/tsutterley/read-cryosat-2/blob/master/cryosat_toolkit/read_geojson_file.py)  

#### Inputs
```
 1. input GeoJSON file (.json, .geojson)
```

#### Options
 - `EPSG`: projection identifier for output coordinates
 - `VARIABLES`: reduce to a specific set of identifiers

#### Outputs
 - `poly_obj`: shapely multipolygon object of input file
