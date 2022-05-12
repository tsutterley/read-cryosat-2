#!/usr/bin/env python
u"""
calc_GPS_time.py (05/2022)
Calculates the GPS time for CryoSat-2 data

INPUTS:
    day: day portion of date variable
    second: seconds portion of date variable
    micsec: microseconds portion of date variable

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

UPDATE HISTORY:
    Updated 05/2022: added docstrings in numpy documentation format
    Updated 08/2020: use convert delta time function in time module
    Written 10/2017
"""
import numpy as np
import cryosat_toolkit.time

#-- PURPOSE: Calculate the GPS time (seconds since Jan 6, 1980 00:00:00)
def calc_GPS_time(day, second, micsec):
    """
    Calculate the GPS time from CryoSat-2 calendar dates
    (days since 2000-01-01T00:00:00 TAI)

    Parameters
    ----------
    day: float
        Days since 2000-01-01
    second: float
        Seconds of the day
    micsec: float
        Microseconds

    Returns:
    GPS_time: float
        GPS time (seconds since 1980-01-06T00:00:00)
    """
    #-- convert to seconds since 2000-01-01 (TAI) time
    delta_time = day*86400.0 + second.astype('f') + micsec/1e6
    #-- TAI time is ahead of GPS by 19 seconds
    return cryosat_toolkit.time.convert_delta_time(delta_time - 19.0,
        epoch1=(2000,1,1,0,0,0), epoch2=(1980,1,6,0,0,0), scale=1.0)
