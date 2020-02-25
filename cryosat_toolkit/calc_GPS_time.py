#!/usr/bin/env python
u"""
calc_GPS_time.py (10/2017)
Calculate the GPS time for CryoSat-2 data for converting from TAI time to UTC
    by calculating the number of leap seconds that have passed for each GPS time

Based partially on Tiffany Summerscales's PHP conversion algorithm
    https://www.andrews.edu/~tzs/timeconv/timealgorithm.html

INPUTS:
    day: day portion of date variable
    second: seconds portion of date variable
    micsec: microseconds portion of date variable

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (http://www.numpy.org)

UPDATE HISTORY:
    Written 10/2017
"""
import numpy as np

#-- PURPOSE: Calculate the GPS time (seconds since Jan 6, 1980 00:00:00)
def calc_GPS_time(day, second, micsec):
    #-- TAI time is ahead of GPS by 19 seconds
    return (day + 7300.0)*86400.0 + second.astype('f') + micsec/1e6 - 19

#-- PURPOSE: Define GPS leap seconds
def get_leaps():
    leaps = [46828800, 78364801, 109900802, 173059203, 252028804, 315187205,
        346723206, 393984007, 425520008, 457056009, 504489610, 551750411,
        599184012, 820108813, 914803214, 1025136015, 1119744016, 1167264017]
    return leaps

#-- PURPOSE: Count number of leap seconds that have passed for given GPS times
def count_leap_seconds(GPS_Time):
    leaps = get_leaps()
    #-- number of leap seconds prior to GPS_Time
    n_leaps = np.zeros_like(GPS_Time)
    for i,leap in enumerate(leaps):
        count = np.count_nonzero(GPS_Time >= leap)
        if (count > 0):
            i_records,i_blocks = np.nonzero(GPS_Time >= leap)
            n_leaps[i_records,i_blocks] += 1.0
    return n_leaps
