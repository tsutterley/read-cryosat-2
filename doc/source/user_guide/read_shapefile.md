read_shapefile.py
=================

 - Reads polygons from ESRI shapefiles

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
