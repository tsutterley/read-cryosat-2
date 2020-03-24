calc_GPS_time.py
================

 - Calculates the GPS time for CryoSat-2 data
 - Calculates the number of leap seconds that have passed for each GPS time
 - Can be used to convert from TAI time to UTC

#### Inputs
 1. `day`: day portion of CryoSat-2 date variable
 2. `second`: seconds portion of CryoSat-2 date variable
 3. `micsec`: microseconds portion of CryoSat-2 date variable

#### Outputs
 - GPS time (seconds since 1980-01-06T00:00:00)
