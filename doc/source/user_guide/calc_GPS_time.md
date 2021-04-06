calc_GPS_time.py
================

 - Calculates the GPS time for CryoSat-2 data
 - Calculates the number of leap seconds that have passed for each GPS time
 - Can be used to convert from TAI time to UTC

#### Calling Sequence
```python
from cryosat_toolkit.calc_GPS_time import calc_GPS_time, count_leap_seconds
GPS_Time = calc_GPS_time(day, second, micsec)
n_leaps = count_leap_seconds(GPS_Time)
```
[Source code](https://github.com/tsutterley/read-cryosat-2/blob/main/cryosat_toolkit/calc_GPS_time.py)

#### Arguments
 1. `day`: day portion of CryoSat-2 date variable
 2. `second`: seconds portion of CryoSat-2 date variable
 3. `micsec`: microseconds portion of CryoSat-2 date variable

#### Returns
 - GPS time (seconds since 1980-01-06T00:00:00)
