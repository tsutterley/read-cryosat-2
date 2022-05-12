#!/usr/bin/env python
u"""
read_cryosat_L2I.py
Written by Tyler Sutterley (05/2022)

Reads CryoSat Level-2 Intermediate data products from baselines A, B, BC and C
Reads CryoSat Level-2 netCDF4 data products from baseline D and E
Supported CryoSat Modes: LRM, SAR, SARin, FDM, SID, GDR

INPUTS:
    full_filename: full path of CryoSat .DBL or .nc file

OUTPUTS:
    Location: Time and Orbit Parameters
    Geometry: Elevation Corrections and Flags
    Data: Geolocation and Elevation Measurements with Quality Parameters
    Auxiliary: Auxiliary Data for Elevation Processing
    Instrumental: Intrument Corrections
    METADATA: MPH, SPH and DSD Header data

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html

UPDATE HISTORY:
    Updated 05/2022: added docstrings in numpy documentation format
    Updated 05/2021: use raw binary string prefixes (rb) for regular expressions
    Updated 08/2020: flake8 updates for python3
    Updated 02/2020: tilde-expansion of cryosat-2 files before opening
        convert from hard to soft tabulation
    Updated 11/2019: empty placeholder dictionary for baseline D DSD headers
    Updated 09/2019: added netCDF4 read function for baseline D
        will output with same variable names as the binary read functions
    Updated 08/2019: generalize regular expression patterns in read_DSD function
    Updated 10/2018: updated header read functions for python3
    Updated 05/2016: using __future__ print and division functions
    Written 03/2016
"""
from __future__ import print_function
from __future__ import division

import os
import re
import netCDF4
import numpy as np
import scipy.interpolate

#-- PURPOSE: Initiate L2I MDS variables for CryoSat Baselines A and B
def cryosat_baseline_AB(fid,record_size,n_records):
    #-- CryoSat-2 Location Group
    #-- Time and Orbit Parameters plus Measurement Mode
    Location = {}
    #-- Time: day part
    Location['Day'] = np.zeros((n_records),dtype=np.int32)
    #-- Time: second part
    Location['Sec'] = np.zeros((n_records),dtype=np.uint32)
    #-- Time: microsecond part
    Location['Micsec'] = np.zeros((n_records),dtype=np.uint32)
    #-- USO correction factor
    Location['USO_Corr'] = np.zeros((n_records),dtype=np.int32)
    #-- Mode ID
    Location['Mode_ID'] = np.zeros((n_records),dtype=np.uint16)
    #-- Source sequence counter
    Location['SSC'] = np.zeros((n_records),dtype=np.uint16)
    #-- Instrument configuration
    Location['Inst_config'] = np.zeros((n_records),dtype=np.uint32)
    #-- Record Counter
    Location['Rec_Count'] = np.zeros((n_records),dtype=np.uint32)
    #-- Lat: packed units (0.1 micro-degree, 1e-7 degrees)
    Location['Lat'] = np.zeros((n_records),dtype=np.int32)
    #-- Lon: packed units (0.1 micro-degree, 1e-7 degrees)
    Location['Lon'] = np.zeros((n_records),dtype=np.int32)
    #-- Alt: packed units (mm, 1e-3 m)
    #-- Altitude of COG above reference ellipsoid (interpolated value)
    Location['Alt'] = np.zeros((n_records),dtype=np.int32)
    #-- Instantaneous altitude rate derived from orbit: packed units (mm/s, 1e-3 m/s)
    Location['Alt_rate'] = np.zeros((n_records),dtype=np.int32)
    #-- Satellite velocity vector. In ITRF: packed units (mm/s, 1e-3 m/s)
    Location['Sat_velocity'] = np.zeros((n_records,3),dtype=np.int32)
    #-- Real beam direction vector. In CRF: packed units (micro-m/s, 1e-6 m/s)
    Location['Real_beam'] = np.zeros((n_records,3),dtype=np.int32)
    #-- Interferometer baseline vector. In CRF: packed units (micro-m/s, 1e-6 m/s)
    Location['Baseline'] = np.zeros((n_records,3),dtype=np.int32)
    #-- Measurement Confidence Data
    Location['MCD'] = np.zeros((n_records),dtype=np.uint32)

    #-- CryoSat-2 Measurement Group
    #-- Derived from instrument measurement parameters
    Data = {}
    #-- Measured elevation above ellipsoid from retracker: packed units (mm, 1e-3 m)
    Data['Elev'] = np.zeros((n_records),dtype=np.int32)
    #-- Sigma Zero Backscatter for retracker: packed units (1e-2 dB)
    Data['Sig0'] = np.zeros((n_records),dtype=np.int32)
    #-- SWH packed units (mm, 1e-3)
    Data['SWH'] = np.zeros((n_records),dtype=np.int32)
    #-- Peakiness: packed units (1e-2)
    Data['Peakiness'] = np.zeros((n_records),dtype=np.int32)
    #-- Retracked range correction: packed units (mm, 1e-3 m)
    Data['Retrack_range'] = np.zeros((n_records),dtype=np.int32)
    #-- Retracked sigma 0 correction: packed units (1e-2 dB)
    Data['Retrack_sig0'] = np.zeros((n_records),dtype=np.int32)
    #-- Retrackers 3-13 output
    Data['Retrack_3'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_4'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_5'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_6'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_7'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_8'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_9'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_10'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_11'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_12'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_13'] = np.zeros((n_records),dtype=np.int32)
    #-- Power echo shape parameter: packed units (dB/100)
    Data['echo_shape'] = np.zeros((n_records),dtype=np.int32)
    #-- Beam behaviour parameter: unitless code number related to
    #-- surface characteristics
    Data['BB_parameter'] = np.zeros((n_records,50),dtype=np.int16)
    #-- Cross track angle: packed units (micro radians)
    Data['X_Track_Angle'] = np.zeros((n_records),dtype=np.int32)
    #-- Leading edge coherence at retrack point 1/1000
    Data['Coherence'] = np.zeros((n_records),dtype=np.int32)
    #-- Interpolated Ocean Height: packed units (mm above ellipsoid)
    Data['Ocean_ht'] = np.zeros((n_records),dtype=np.int32)
    #-- Freeboard: packed units (mm, 1e-3 m)
    #-- -9999 default value indicates computation has not been performed
    Data['Freeboard'] = np.zeros((n_records),dtype=np.int32)
    #-- Surface Height Anomaly: packed units (mm, 1e-3 m)
    Data['SHA'] = np.zeros((n_records),dtype=np.int32)
    #-- Interpolated Surface Height Anomaly: packed units (mm, 1e-3 m)
    Data['SSHA_interp'] = np.zeros((n_records),dtype=np.int32)
    #-- Error in ocean height interpolation: packed units (mm, 1e-3 m)
    Data['SSHA_interp_RMS'] = np.zeros((n_records),dtype=np.uint16)
    #-- Number of forward records interpolated
    Data['SSHA_interp_count_fwd'] = np.zeros((n_records),dtype=np.uint16)
    #-- Number of backward records interpolated
    Data['SSHA_interp_count_bkwd'] = np.zeros((n_records),dtype=np.uint16)
    #-- Distance in time of most forward record interpolated (milli-seconds)
    Data['SSHA_interp_time_fwd'] = np.zeros((n_records),dtype=np.uint16)
    #-- Distance in time of most backward record interpolated (milli-seconds)
    Data['SSHA_interp_time_bkwd'] = np.zeros((n_records),dtype=np.uint16)
    #-- Interpolation error flag
    Data['SSHA_interp_flag'] = np.zeros((n_records),dtype=np.uint16)
    #-- Measurement mode
    Data['Measurement_Mode'] = np.zeros((n_records),dtype=np.uint32)
    #-- Quality flags
    Data['Quality_flag'] = np.zeros((n_records),dtype=np.uint32)
    #-- Retracker flags
    Data['Retracker_flag'] = np.zeros((n_records),dtype=np.uint32)
    #-- Height calculation details
    #-- Specifies what was applied during the height calculation
    Data['Height_status'] = np.zeros((n_records),dtype=np.uint32)
    #-- SAR freeboard status flag
    Data['Freeboard_status'] = np.zeros((n_records),dtype=np.uint32)
    #-- Number of averaged echoes or beams
    Data['N_avg'] = np.zeros((n_records),dtype=np.uint16)
    #-- Wind Speed packed units (mm/s, 1e-3 m/s)
    Data['Wind_speed'] = np.zeros((n_records),dtype=np.uint16)
    Data['Spares1'] = np.zeros((n_records,3),dtype=np.int32)

    #-- CryoSat-2 Auxiliary Data Group
    Auxiliary = {}
    #-- Ice Concentration packed units (%/1000)
    Auxiliary['Ice_conc'] = np.zeros((n_records),dtype=np.int32)
    #-- Snow Depth packed units (mm, 1e-3 m)
    Auxiliary['Snow_depth'] = np.zeros((n_records),dtype=np.int32)
    #-- Snow Density packed units (kg/m^3)
    Auxiliary['Snow_density'] = np.zeros((n_records),dtype=np.int32)
    #-- Discriminator result
    Auxiliary['Discriminator'] = np.zeros((n_records),dtype=np.int32)
    #-- SARin discriminator parameters 1-10
    Auxiliary['SARIN_disc_1'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['SARIN_disc_2'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['SARIN_disc_3'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['SARIN_disc_4'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['SARIN_disc_5'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['SARIN_disc_6'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['SARIN_disc_7'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['SARIN_disc_8'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['SARIN_disc_9'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['SARIN_disc_10'] = np.zeros((n_records),dtype=np.int32)
    #-- Discriminator flags
    Auxiliary['Discrim_flag'] = np.zeros((n_records),dtype=np.uint32)
    #-- Slope model correction (Attitude of echo in micro-degrees)
    Auxiliary['Attitude'] = np.zeros((n_records),dtype=np.int32)
    #-- Slope model correction (Azimuth of echo in micro-degrees)
    Auxiliary['Azimuth'] = np.zeros((n_records),dtype=np.int32)
    #-- The original latitude of the satellite (micro-degrees)
    Auxiliary['Lat_sat'] = np.zeros((n_records),dtype=np.int32)
    #-- The original longitude of the satellite (micro-degrees)
    Auxiliary['Lon_sat'] = np.zeros((n_records),dtype=np.int32)
    #-- Ambiguity indicator
    Auxiliary['Ambiguity'] = np.zeros((n_records),dtype=np.uint32)
    #-- Mean Sea Surface standard Model: packed units (mm, 1e-3 m)
    Auxiliary['MSS_model'] = np.zeros((n_records),dtype=np.int32)
    #-- Geoid standard Model: packed units (mm, 1e-3 m)
    Auxiliary['Geoid_model'] = np.zeros((n_records),dtype=np.int32)
    #-- ODLE standard Model: packed units (mm, 1e-3 m)
    Auxiliary['ODLE'] = np.zeros((n_records),dtype=np.int32)
    #-- The interpolated elevation value obtained from the DEM (mm)
    Auxiliary['DEM_elev'] = np.zeros((n_records),dtype=np.int32)
    #-- Identification of DEM used in SARin ambiguity test
    Auxiliary['DEM_ID'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['Spares2'] = np.zeros((n_records,4),dtype=np.int32)

    #-- CryoSat-2 External Corrections Group
    Geometry = {}
    #-- Dry Tropospheric Correction packed units (mm, 1e-3 m)
    Geometry['dryTrop'] = np.zeros((n_records),dtype=np.int32)
    #-- Wet Tropospheric Correction packed units (mm, 1e-3 m)
    Geometry['wetTrop'] = np.zeros((n_records),dtype=np.int32)
    #-- Inverse Barometric Correction packed units (mm, 1e-3 m)
    Geometry['InvBar'] = np.zeros((n_records),dtype=np.int32)
    #-- Delta Inverse Barometric Correction packed units (mm, 1e-3 m)
    Geometry['DAC'] = np.zeros((n_records),dtype=np.int32)
    #-- GIM Ionospheric Correction packed units (mm, 1e-3 m)
    Geometry['Iono_GIM'] = np.zeros((n_records),dtype=np.int32)
    #-- Model Ionospheric Correction packed units (mm, 1e-3 m)
    Geometry['Iono_model'] = np.zeros((n_records),dtype=np.int32)
    #-- Ocean tide Correction packed units (mm, 1e-3 m)
    Geometry['ocTideElv'] = np.zeros((n_records),dtype=np.int32)
    #-- Long period equilibrium ocean tide Correction packed units (mm, 1e-3 m)
    Geometry['lpeTideElv'] = np.zeros((n_records),dtype=np.int32)
    #-- Ocean loading tide Correction packed units (mm, 1e-3 m)
    Geometry['olTideElv'] = np.zeros((n_records),dtype=np.int32)
    #-- Solid Earth tide Correction packed units (mm, 1e-3 m)
    Geometry['seTideElv'] = np.zeros((n_records),dtype=np.int32)
    #-- Geocentric Polar tide Correction packed units (mm, 1e-3 m)
    Geometry['gpTideElv'] = np.zeros((n_records),dtype=np.int32)
    #-- Surface Type: Packed in groups of three bits for each of the 20 records
    Geometry['Surf_type'] = np.zeros((n_records),dtype=np.uint32)
    #-- Corrections Status Flag
    Geometry['Corr_status'] = np.zeros((n_records),dtype=np.uint32)
    #-- Correction Error Flag
    Geometry['Corr_error'] = np.zeros((n_records),dtype=np.uint32)
    #-- Sea State Bias Correction packed units (mm, 1e-3 m)
    Geometry['SSB'] = np.zeros((n_records),dtype=np.int32)
    Geometry['Spares3'] = np.zeros((n_records,2),dtype=np.int32)

    #-- CryoSat-2 Internal Corrections Group
    Instrumental = {}
    #-- Doppler range correction: Radial + slope (mm)
    Instrumental['Doppler_range'] = np.zeros((n_records),dtype=np.int32)
    #-- Instrument Range Correction: t-r antenna (mm)
    Instrumental['TR_inst_range'] = np.zeros((n_records),dtype=np.int32)
    #-- Instrument Range Correction: r-only antenna (mm)
    Instrumental['R_inst_range'] = np.zeros((n_records),dtype=np.int32)
    #-- Instrument Sigma 0 Correction: t-r antenna (dB/100)
    Instrumental['TR_inst_gain'] = np.zeros((n_records),dtype=np.int32)
    #-- Instrument Sigma 0 Correction: r-only (dB/100)
    Instrumental['R_inst_gain'] = np.zeros((n_records),dtype=np.int32)
    #-- Internal Phase Correction (milli-radians)
    Instrumental['Internal_phase'] = np.zeros((n_records),dtype=np.int32)
    #-- External Phase Correction (milli-radians)
    Instrumental['External_phase'] = np.zeros((n_records),dtype=np.int32)
    #-- Noise Power measurement
    Instrumental['Noise_power'] = np.zeros((n_records),dtype=np.int32)
    #-- Phase slope correction (microradians)
    Instrumental['Phase_slope'] = np.zeros((n_records),dtype=np.int32)
    Instrumental['Spares4'] = np.zeros((n_records,2),dtype=np.int32)

    #-- for each record in the CryoSat file
    for r in range(n_records):
        #-- get satellite time and orbit parameters for record r
        Location['Day'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Location['Sec'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Location['Micsec'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Location['USO_Corr'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Location['Mode_ID'][r] = np.fromfile(fid,dtype='>u2',count=1)
        Location['SSC'][r] = np.fromfile(fid,dtype='>u2',count=1)
        Location['Inst_config'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Location['Rec_Count'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Location['Lat'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Location['Lon'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Location['Alt'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Location['Alt_rate'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Location['Sat_velocity'][r,:] = np.fromfile(fid,dtype='>i4',count=3)
        Location['Real_beam'][r,:] = np.fromfile(fid,dtype='>i4',count=3)
        Location['Baseline'][r,:] = np.fromfile(fid,dtype='>i4',count=3)
        Location['MCD'][r] = np.fromfile(fid,dtype='>u4',count=1)

        #-- elevation measurements
        Data['Elev'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Sig0'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['SWH'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Peakiness'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_range'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_sig0'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_3'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_4'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_5'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_6'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_7'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_8'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_9'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_10'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_11'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_12'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_13'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['echo_shape'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['BB_parameter'][r,:] = np.fromfile(fid,dtype='>i2',count=50)
        Data['X_Track_Angle'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Coherence'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Ocean_ht'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Freeboard'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['SHA'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['SSHA_interp'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['SSHA_interp_RMS'][r] = np.fromfile(fid,dtype='>u2',count=1)
        Data['SSHA_interp_count_fwd'][r] = np.fromfile(fid,dtype='>u2',count=1)
        Data['SSHA_interp_count_bkwd'][r] = np.fromfile(fid,dtype='>u2',count=1)
        Data['SSHA_interp_time_fwd'][r] = np.fromfile(fid,dtype='>u2',count=1)
        Data['SSHA_interp_time_bkwd'][r] = np.fromfile(fid,dtype='>u2',count=1)
        Data['SSHA_interp_flag'][r] = np.fromfile(fid,dtype='>u2',count=1)
        Data['Measurement_Mode'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Data['Quality_flag'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Data['Retracker_flag'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Data['Height_status'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Data['Freeboard_status'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Data['N_avg'][r] = np.fromfile(fid,dtype='>u2',count=1)
        Data['Wind_speed'][r] = np.fromfile(fid,dtype='>u2',count=1)
        Data['Spares1'][r,:] = np.fromfile(fid,dtype='>i4',count=3)

        #-- Auxiliary Data
        Auxiliary['Ice_conc'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Snow_depth'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Snow_density'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Discriminator'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['SARIN_disc_1'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['SARIN_disc_2'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['SARIN_disc_3'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['SARIN_disc_4'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['SARIN_disc_5'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['SARIN_disc_6'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['SARIN_disc_7'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['SARIN_disc_8'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['SARIN_disc_9'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['SARIN_disc_10'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Discrim_flag'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Auxiliary['Attitude'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Azimuth'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Lat_sat'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Lon_sat'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Ambiguity'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Auxiliary['MSS_model'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Geoid_model'][r] =np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['ODLE'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['DEM_elev'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['DEM_ID'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Spares2'][r,:] = np.fromfile(fid,dtype='>i4',count=4)

        #-- CryoSat-2 External Corrections Group for record r
        Geometry['dryTrop'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['wetTrop'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['InvBar'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['DAC'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['Iono_GIM'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['Iono_model'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['ocTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['lpeTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['olTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['seTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['gpTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['Surf_type'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Geometry['Corr_status'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Geometry['Corr_error'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Geometry['SSB'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['Spares3'][r,:] = np.fromfile(fid,dtype='>i4',count=2)

        #-- CryoSat-2 Internal Corrections Group for record r
        Instrumental['Doppler_range'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Instrumental['TR_inst_range'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Instrumental['R_inst_range'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Instrumental['TR_inst_gain'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Instrumental['R_inst_gain'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Instrumental['Internal_phase'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Instrumental['External_phase'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Instrumental['Noise_power'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Instrumental['Phase_slope'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Instrumental['Spares4'][r,:] = np.fromfile(fid,dtype='>i4',count=2)

    #-- Bind all the bits of the l2i_mds together into a single dictionary
    CS_L2I_mds = {}
    CS_L2I_mds['Location'] = Location
    CS_L2I_mds['Data'] = Data
    CS_L2I_mds['Auxiliary'] = Auxiliary
    CS_L2I_mds['Geometry'] = Geometry
    CS_L2I_mds['Instrumental'] = Instrumental
    #-- return the output dictionary
    return CS_L2I_mds

#-- PURPOSE: Initiate L2I MDS variables for CryoSat Baseline BC
def cryosat_baseline_BC(fid,record_size,n_records):
    #-- CryoSat-2 Location Group
    #-- Time and Orbit Parameters plus Measurement Mode
    Location = {}
    #-- Time: day part
    Location['Day'] = np.zeros((n_records),dtype=np.int32)
    #-- Time: second part
    Location['Sec'] = np.zeros((n_records),dtype=np.uint32)
    #-- Time: microsecond part
    Location['Micsec'] = np.zeros((n_records),dtype=np.uint32)
    #-- USO correction factor
    Location['USO_Corr'] = np.zeros((n_records),dtype=np.int32)
    #-- Mode ID
    Location['Mode_ID'] = np.zeros((n_records),dtype=np.uint16)
    #-- Source sequence counter
    Location['SSC'] = np.zeros((n_records),dtype=np.uint16)
    #-- Instrument configuration
    Location['Inst_config'] = np.zeros((n_records),dtype=np.uint32)
    #-- Record Counter
    Location['Rec_Count'] = np.zeros((n_records),dtype=np.uint32)
    #-- Lat: packed units (0.1 micro-degree, 1e-7 degrees)
    Location['Lat'] = np.zeros((n_records),dtype=np.int32)
    #-- Lon: packed units (0.1 micro-degree, 1e-7 degrees)
    Location['Lon'] = np.zeros((n_records),dtype=np.int32)
    #-- Alt: packed units (mm, 1e-3 m)
    #-- Altitude of COG above reference ellipsoid (interpolated value)
    Location['Alt'] = np.zeros((n_records),dtype=np.int32)
    #-- Instantaneous altitude rate derived from orbit: packed units (mm/s, 1e-3 m/s)
    Location['Alt_rate'] = np.zeros((n_records),dtype=np.int32)
    #-- Satellite velocity vector. In ITRF: packed units (mm/s, 1e-3 m/s)
    Location['Sat_velocity'] = np.zeros((n_records,3),dtype=np.int32)
    #-- Real beam direction vector. In CRF: packed units (micro-m/s, 1e-6 m/s)
    Location['Real_beam'] = np.zeros((n_records,3),dtype=np.int32)
    #-- Interferometer baseline vector. In CRF: packed units (micro-m/s, 1e-6 m/s)
    Location['Baseline'] = np.zeros((n_records,3),dtype=np.int32)
    #-- Star Tracker ID
    Location['ST_ID'] = np.zeros((n_records),dtype=np.int16)
    Location['Spare'] = np.zeros((n_records),dtype=np.int16)
    #-- Roll (Derived from star trackers): packed units (0.1 micro-degree, 1e-7 degrees)
    Location['Roll'] = np.zeros((n_records),dtype=np.int32)
    #-- Pitch (Derived from star trackers): packed units (0.1 micro-degree, 1e-7 degrees)
    Location['Pitch'] = np.zeros((n_records),dtype=np.int32)
    #-- Yaw (Derived from star trackers): packed units (0.1 micro-degree, 1e-7 degrees)
    Location['Yaw'] = np.zeros((n_records),dtype=np.int32)
    #-- Measurement Confidence Data
    Location['MCD'] = np.zeros((n_records),dtype=np.uint32)

    #-- CryoSat-2 Measurement Group
    #-- Derived from instrument measurement parameters
    Data = {}
    #-- Measured elevation above ellipsoid from retracker 1: packed units (mm, 1e-3 m)
    Data['Elev_1'] = np.zeros((n_records),dtype=np.int32)
    #-- Measured elevation above ellipsoid from retracker 2: packed units (mm, 1e-3 m)
    Data['Elev_2'] = np.zeros((n_records),dtype=np.int32)
    #-- Measured elevation above ellipsoid from retracker 3: packed units (mm, 1e-3 m)
    Data['Elev_3'] = np.zeros((n_records),dtype=np.int32)
    #-- Sigma Zero Backscatter for retracker 1: packed units (1e-2 dB)
    Data['Sig0_1'] = np.zeros((n_records),dtype=np.int32)
    #-- Sigma Zero Backscatter for retracker 2: packed units (1e-2 dB)
    Data['Sig0_2'] = np.zeros((n_records),dtype=np.int32)
    #-- Sigma Zero Backscatter for retracker 3: packed units (1e-2 dB)
    Data['Sig0_3'] = np.zeros((n_records),dtype=np.int32)
    #-- SWH packed units (mm, 1e-3)
    Data['SWH'] = np.zeros((n_records),dtype=np.int32)
    #-- Peakiness: packed units (1e-2)
    Data['Peakiness'] = np.zeros((n_records),dtype=np.int32)
    #-- Retracked range correction for retracker 1: packed units (mm, 1e-3 m)
    Data['Range_1'] = np.zeros((n_records),dtype=np.int32)
    #-- Retracked range correction for retracker 2: packed units (mm, 1e-3 m)
    Data['Range_2'] = np.zeros((n_records),dtype=np.int32)
    #-- Retracked range correction for retracker 3: packed units (mm, 1e-3 m)
    Data['Range_3'] = np.zeros((n_records),dtype=np.int32)
    #-- Retracked sigma 0 correction for Retracker 1: packed units (1e-2 dB)
    Data['Retrack_1_sig0'] = np.zeros((n_records),dtype=np.int32)
    #-- Retracked sigma 0 correction for Retracker 2: packed units (1e-2 dB)
    Data['Retrack_2_sig0'] = np.zeros((n_records),dtype=np.int32)
    #-- Retracked sigma 0 correction for Retracker 3: packed units (1e-2 dB)
    Data['Retrack_3_sig0'] = np.zeros((n_records),dtype=np.int32)
    #-- Retracker 1 quality metric
    Data['Quality_1'] = np.zeros((n_records),dtype=np.int32)
    #-- Retracker 2 quality metric
    Data['Quality_2'] = np.zeros((n_records),dtype=np.int32)
    #-- Retracker 3 quality metric
    Data['Quality_3'] = np.zeros((n_records),dtype=np.int32)
    #-- Retrackers 3-23 output
    Data['Retrack_3'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_4'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_5'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_6'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_7'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_8'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_9'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_10'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_11'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_12'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_13'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_14'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_15'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_16'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_17'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_18'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_19'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_20'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_21'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_22'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_23'] = np.zeros((n_records),dtype=np.int32)
    #-- Power echo shape parameter: packed units (dB/100)
    Data['echo_shape'] = np.zeros((n_records),dtype=np.int32)
    #-- Beam behaviour parameter: unitless code number related to
    #-- surface characteristics
    Data['BB_parameter'] = np.zeros((n_records,50),dtype=np.int16)
    #-- Cross track angle: packed units (micro radians)
    Data['X_Track_Angle'] = np.zeros((n_records),dtype=np.int32)
    #-- Cross track angle correction: packed units (micro radians)
    Data['X_Track_Angle_c'] = np.zeros((n_records),dtype=np.int32)
    #-- Leading edge coherence at retrack point 1/1000
    Data['Coherence'] = np.zeros((n_records),dtype=np.int32)
    #-- Interpolated Ocean Height: packed units (mm above ellipsoid)
    Data['Ocean_ht'] = np.zeros((n_records),dtype=np.int32)
    #-- Freeboard: packed units (mm, 1e-3 m)
    #-- -9999 default value indicates computation has not been performed
    Data['Freeboard'] = np.zeros((n_records),dtype=np.int32)
    #-- Surface Height Anomaly: packed units (mm, 1e-3 m)
    Data['SHA'] = np.zeros((n_records),dtype=np.int32)
    #-- Interpolated Surface Height Anomaly: packed units (mm, 1e-3 m)
    Data['SSHA_interp'] = np.zeros((n_records),dtype=np.int32)
    #-- Error in ocean height interpolation: packed units (mm, 1e-3 m)
    Data['SSHA_interp_RMS'] = np.zeros((n_records),dtype=np.uint16)
    #-- Number of forward records interpolated
    Data['SSHA_interp_count_fwd'] = np.zeros((n_records),dtype=np.uint16)
    #-- Number of backward records interpolated
    Data['SSHA_interp_count_bkwd'] = np.zeros((n_records),dtype=np.uint16)
    #-- Distance in time of most forward record interpolated (milli-seconds)
    Data['SSHA_interp_time_fwd'] = np.zeros((n_records),dtype=np.uint16)
    #-- Distance in time of most backward record interpolated (milli-seconds)
    Data['SSHA_interp_time_bkwd'] = np.zeros((n_records),dtype=np.uint16)
    #-- Interpolation error flag
    Data['SSHA_interp_flag'] = np.zeros((n_records),dtype=np.uint16)
    #-- Measurement mode
    Data['Measurement_Mode'] = np.zeros((n_records),dtype=np.uint32)
    #-- Quality flags
    Data['Quality_flag'] = np.zeros((n_records),dtype=np.uint32)
    #-- Retracker flags
    Data['Retracker_flag'] = np.zeros((n_records),dtype=np.uint32)
    #-- Height calculation details
    #-- Specifies what was applied during the height calculation
    Data['Height_status'] = np.zeros((n_records),dtype=np.uint32)
    #-- SAR freeboard status flag
    Data['Freeboard_status'] = np.zeros((n_records),dtype=np.uint32)
    #-- Number of averaged echoes or beams
    Data['N_avg'] = np.zeros((n_records),dtype=np.uint16)
    #-- Wind Speed packed units (mm/s, 1e-3 m/s)
    Data['Wind_speed'] = np.zeros((n_records),dtype=np.uint16)
    Data['Spares1'] = np.zeros((n_records,3),dtype=np.int32)

    #-- CryoSat-2 Auxiliary Data Group
    Auxiliary = {}
    #-- Ice Concentration packed units (%/1000)
    Auxiliary['Ice_conc'] = np.zeros((n_records),dtype=np.int32)
    #-- Snow Depth packed units (mm, 1e-3 m)
    Auxiliary['Snow_depth'] = np.zeros((n_records),dtype=np.int32)
    #-- Snow Density packed units (kg/m^3)
    Auxiliary['Snow_density'] = np.zeros((n_records),dtype=np.int32)
    #-- Discriminator result
    Auxiliary['Discriminator'] = np.zeros((n_records),dtype=np.int32)
    #-- SARin discriminator parameters 1-10
    Auxiliary['SARIN_disc_1'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['SARIN_disc_2'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['SARIN_disc_3'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['SARIN_disc_4'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['SARIN_disc_5'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['SARIN_disc_6'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['SARIN_disc_7'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['SARIN_disc_8'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['SARIN_disc_9'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['SARIN_disc_10'] = np.zeros((n_records),dtype=np.int32)
    #-- Discriminator flags
    Auxiliary['Discrim_flag'] = np.zeros((n_records),dtype=np.uint32)
    #-- Slope model correction (Attitude of echo in micro-degrees)
    Auxiliary['Attitude'] = np.zeros((n_records),dtype=np.int32)
    #-- Slope model correction (Azimuth of echo in micro-degrees)
    Auxiliary['Azimuth'] = np.zeros((n_records),dtype=np.int32)
    #-- Slope doppler correction (mm)
    Auxiliary['Slope_doppler'] = np.zeros((n_records),dtype=np.int32)
    #-- The original latitude of the satellite (micro-degrees)
    Auxiliary['Lat_sat'] = np.zeros((n_records),dtype=np.int32)
    #-- The original longitude of the satellite (micro-degrees)
    Auxiliary['Lon_sat'] = np.zeros((n_records),dtype=np.int32)
    #-- Ambiguity indicator
    Auxiliary['Ambiguity'] = np.zeros((n_records),dtype=np.uint32)
    #-- Mean Sea Surface standard Model: packed units (mm, 1e-3 m)
    Auxiliary['MSS_model'] = np.zeros((n_records),dtype=np.int32)
    #-- Geoid standard Model: packed units (mm, 1e-3 m)
    Auxiliary['Geoid_model'] = np.zeros((n_records),dtype=np.int32)
    #-- ODLE standard Model: packed units (mm, 1e-3 m)
    Auxiliary['ODLE'] = np.zeros((n_records),dtype=np.int32)
    #-- The interpolated elevation value obtained from the DEM (mm)
    Auxiliary['DEM_elev'] = np.zeros((n_records),dtype=np.int32)
    #-- Identification of DEM used in SARin ambiguity test
    Auxiliary['DEM_ID'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['Spares2'] = np.zeros((n_records,4),dtype=np.int32)

    #-- CryoSat-2 External Corrections Group
    Geometry = {}
    #-- Dry Tropospheric Correction packed units (mm, 1e-3 m)
    Geometry['dryTrop'] = np.zeros((n_records),dtype=np.int32)
    #-- Wet Tropospheric Correction packed units (mm, 1e-3 m)
    Geometry['wetTrop'] = np.zeros((n_records),dtype=np.int32)
    #-- Inverse Barometric Correction packed units (mm, 1e-3 m)
    Geometry['InvBar'] = np.zeros((n_records),dtype=np.int32)
    #-- Delta Inverse Barometric Correction packed units (mm, 1e-3 m)
    Geometry['DAC'] = np.zeros((n_records),dtype=np.int32)
    #-- GIM Ionospheric Correction packed units (mm, 1e-3 m)
    Geometry['Iono_GIM'] = np.zeros((n_records),dtype=np.int32)
    #-- Model Ionospheric Correction packed units (mm, 1e-3 m)
    Geometry['Iono_model'] = np.zeros((n_records),dtype=np.int32)
    #-- Ocean tide Correction packed units (mm, 1e-3 m)
    Geometry['ocTideElv'] = np.zeros((n_records),dtype=np.int32)
    #-- Long period equilibrium ocean tide Correction packed units (mm, 1e-3 m)
    Geometry['lpeTideElv'] = np.zeros((n_records),dtype=np.int32)
    #-- Ocean loading tide Correction packed units (mm, 1e-3 m)
    Geometry['olTideElv'] = np.zeros((n_records),dtype=np.int32)
    #-- Solid Earth tide Correction packed units (mm, 1e-3 m)
    Geometry['seTideElv'] = np.zeros((n_records),dtype=np.int32)
    #-- Geocentric Polar tide Correction packed units (mm, 1e-3 m)
    Geometry['gpTideElv'] = np.zeros((n_records),dtype=np.int32)
    #-- Surface Type: Packed in groups of three bits for each of the 20 records
    Geometry['Surf_type'] = np.zeros((n_records),dtype=np.uint32)
    #-- Corrections Status Flag
    Geometry['Corr_status'] = np.zeros((n_records),dtype=np.uint32)
    #-- Correction Error Flag
    Geometry['Corr_error'] = np.zeros((n_records),dtype=np.uint32)
    #-- Sea State Bias Correction packed units (mm, 1e-3 m)
    Geometry['SSB'] = np.zeros((n_records),dtype=np.int32)
    Geometry['Spares3'] = np.zeros((n_records,2),dtype=np.int32)

    #-- CryoSat-2 Internal Corrections Group
    Instrumental = {}
    #-- Doppler range correction: Radial + slope (mm)
    Instrumental['Doppler_range'] = np.zeros((n_records),dtype=np.int32)
    #-- Instrument Range Correction: t-r antenna (mm)
    Instrumental['TR_inst_range'] = np.zeros((n_records),dtype=np.int32)
    #-- Instrument Range Correction: r-only antenna (mm)
    Instrumental['R_inst_range'] = np.zeros((n_records),dtype=np.int32)
    #-- Instrument Sigma 0 Correction: t-r antenna (dB/100)
    Instrumental['TR_inst_gain'] = np.zeros((n_records),dtype=np.int32)
    #-- Instrument Sigma 0 Correction: r-only (dB/100)
    Instrumental['R_inst_gain'] = np.zeros((n_records),dtype=np.int32)
    #-- Internal Phase Correction (milli-radians)
    Instrumental['Internal_phase'] = np.zeros((n_records),dtype=np.int32)
    #-- External Phase Correction (milli-radians)
    Instrumental['External_phase'] = np.zeros((n_records),dtype=np.int32)
    #-- Noise Power measurement
    Instrumental['Noise_power'] = np.zeros((n_records),dtype=np.int32)
    #-- Phase slope correction (microradians)
    Instrumental['Phase_slope'] = np.zeros((n_records),dtype=np.int32)
    Instrumental['Spares4'] = np.zeros((n_records,2),dtype=np.int32)

    #-- for each record in the CryoSat file
    for r in range(n_records):
        #-- CryoSat-2 Location Group for record r
        Location['Day'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Location['Sec'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Location['Micsec'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Location['USO_Corr'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Location['Mode_ID'][r] = np.fromfile(fid,dtype='>u2',count=1)
        Location['SSC'][r] = np.fromfile(fid,dtype='>u2',count=1)
        Location['Inst_config'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Location['Rec_Count'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Location['Lat'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Location['Lon'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Location['Alt'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Location['Alt_rate'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Location['Sat_velocity'][r,:] = np.fromfile(fid,dtype='>i4',count=3)
        Location['Real_beam'][r,:] = np.fromfile(fid,dtype='>i4',count=3)
        Location['Baseline'][r,:] = np.fromfile(fid,dtype='>i4',count=3)
        Location['ST_ID'][r] = np.fromfile(fid,dtype='>i2',count=1)
        Location['Spare'][r] = np.fromfile(fid,dtype='>i2',count=1)
        Location['Roll'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Location['Pitch'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Location['Yaw'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Location['MCD'][r] = np.fromfile(fid,dtype='>u4',count=1)

        #-- CryoSat-2 Measurement Group for record r
        Data['Elev_1'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Elev_2'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Elev_3'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Sig0_1'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Sig0_2'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Sig0_3'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['SWH'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Peakiness'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Range_1'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Range_2'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Range_3'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_1_sig0'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_2_sig0'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_3_sig0'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Quality_1'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Quality_2'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Quality_3'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_3'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_4'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_5'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_6'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_7'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_8'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_9'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_10'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_11'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_12'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_13'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_14'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_15'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_16'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_17'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_18'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_19'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_20'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_21'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_22'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_23'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['echo_shape'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['BB_parameter'][r,:] = np.fromfile(fid,dtype='>i2',count=50)
        Data['X_Track_Angle'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['X_Track_Angle_c'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Coherence'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Ocean_ht'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Freeboard'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['SHA'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['SSHA_interp'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['SSHA_interp_RMS'][r] = np.fromfile(fid,dtype='>u2',count=1)
        Data['SSHA_interp_count_fwd'][r] = np.fromfile(fid,dtype='>u2',count=1)
        Data['SSHA_interp_count_bkwd'][r] = np.fromfile(fid,dtype='>u2',count=1)
        Data['SSHA_interp_time_fwd'][r] = np.fromfile(fid,dtype='>u2',count=1)
        Data['SSHA_interp_time_bkwd'][r] = np.fromfile(fid,dtype='>u2',count=1)
        Data['SSHA_interp_flag'][r] = np.fromfile(fid,dtype='>u2',count=1)
        Data['Measurement_Mode'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Data['Quality_flag'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Data['Retracker_flag'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Data['Height_status'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Data['Freeboard_status'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Data['N_avg'][r] = np.fromfile(fid,dtype='>u2',count=1)
        Data['Wind_speed'][r] = np.fromfile(fid,dtype='>u2',count=1)
        Data['Spares1'][r,:] = np.fromfile(fid,dtype='>i4',count=3)

        #-- CryoSat-2 Auxiliary Data Group for record r
        Auxiliary['Ice_conc'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Snow_depth'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Snow_density'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Discriminator'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['SARIN_disc_1'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['SARIN_disc_2'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['SARIN_disc_3'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['SARIN_disc_4'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['SARIN_disc_5'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['SARIN_disc_6'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['SARIN_disc_7'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['SARIN_disc_8'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['SARIN_disc_9'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['SARIN_disc_10'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Discrim_flag'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Auxiliary['Attitude'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Azimuth'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Slope_doppler'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Lat_sat'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Lon_sat'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Ambiguity'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Auxiliary['MSS_model'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Geoid_model'][r] =np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['ODLE'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['DEM_elev'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['DEM_ID'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Spares2'][r,:] = np.fromfile(fid,dtype='>i4',count=4)

        #-- CryoSat-2 External Corrections Group for record r
        Geometry['dryTrop'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['wetTrop'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['InvBar'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['DAC'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['Iono_GIM'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['Iono_model'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['ocTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['lpeTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['olTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['seTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['gpTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['Surf_type'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Geometry['Corr_status'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Geometry['Corr_error'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Geometry['SSB'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['Spares3'][r,:] = np.fromfile(fid,dtype='>i4',count=2)

        #-- CryoSat-2 Internal Corrections Group for record r
        Instrumental['Doppler_range'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Instrumental['TR_inst_range'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Instrumental['R_inst_range'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Instrumental['TR_inst_gain'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Instrumental['R_inst_gain'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Instrumental['Internal_phase'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Instrumental['External_phase'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Instrumental['Noise_power'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Instrumental['Phase_slope'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Instrumental['Spares4'][r,:] = np.fromfile(fid,dtype='>i4',count=2)

    #-- Bind all the bits of the l2i_mds together into a single dictionary
    CS_L2I_mds = {}
    CS_L2I_mds['Location'] = Location
    CS_L2I_mds['Data'] = Data
    CS_L2I_mds['Auxiliary'] = Auxiliary
    CS_L2I_mds['Geometry'] = Geometry
    CS_L2I_mds['Instrumental'] = Instrumental
    #-- return the output dictionary
    return CS_L2I_mds

#-- PURPOSE: Initiate L2I MDS variables for CryoSat Baseline C
def cryosat_baseline_C(fid,record_size,n_records):
    #-- CryoSat-2 Location Group
    #-- Time and Orbit Parameters plus Measurement Mode
    Location = {}
    #-- Time: day part
    Location['Day'] = np.zeros((n_records),dtype=np.int32)
    #-- Time: second part
    Location['Sec'] = np.zeros((n_records),dtype=np.uint32)
    #-- Time: microsecond part
    Location['Micsec'] = np.zeros((n_records),dtype=np.uint32)
    #-- USO correction factor
    Location['USO_Corr'] = np.zeros((n_records),dtype=np.int32)
    #-- Mode ID
    Location['Mode_ID'] = np.zeros((n_records),dtype=np.uint16)
    #-- Source sequence counter
    Location['SSC'] = np.zeros((n_records),dtype=np.uint16)
    #-- Instrument configuration
    Location['Inst_config'] = np.zeros((n_records),dtype=np.uint32)
    #-- Record Counter
    Location['Rec_Count'] = np.zeros((n_records),dtype=np.uint32)
    #-- Lat: packed units (0.1 micro-degree, 1e-7 degrees)
    Location['Lat'] = np.zeros((n_records),dtype=np.int32)
    #-- Lon: packed units (0.1 micro-degree, 1e-7 degrees)
    Location['Lon'] = np.zeros((n_records),dtype=np.int32)
    #-- Alt: packed units (mm, 1e-3 m)
    #-- Altitude of COG above reference ellipsoid (interpolated value)
    Location['Alt'] = np.zeros((n_records),dtype=np.int32)
    #-- Instantaneous altitude rate derived from orbit: packed units (mm/s, 1e-3 m/s)
    Location['Alt_rate'] = np.zeros((n_records),dtype=np.int32)
    #-- Satellite velocity vector. In ITRF: packed units (mm/s, 1e-3 m/s)
    Location['Sat_velocity'] = np.zeros((n_records,3),dtype=np.int32)
    #-- Real beam direction vector. In CRF: packed units (micro-m/s, 1e-6 m/s)
    Location['Real_beam'] = np.zeros((n_records,3),dtype=np.int32)
    #-- Interferometer baseline vector. In CRF: packed units (micro-m/s, 1e-6 m/s)
    Location['Baseline'] = np.zeros((n_records,3),dtype=np.int32)
    #-- Star Tracker ID
    Location['ST_ID'] = np.zeros((n_records),dtype=np.int16)
    Location['Spare'] = np.zeros((n_records),dtype=np.int16)
    #-- Roll (Derived from star trackers): packed units (0.1 micro-degree, 1e-7 degrees)
    Location['Roll'] = np.zeros((n_records),dtype=np.int32)
    #-- Pitch (Derived from star trackers): packed units (0.1 micro-degree, 1e-7 degrees)
    Location['Pitch'] = np.zeros((n_records),dtype=np.int32)
    #-- Yaw (Derived from star trackers): packed units (0.1 micro-degree, 1e-7 degrees)
    Location['Yaw'] = np.zeros((n_records),dtype=np.int32)
    #-- Measurement Confidence Data
    Location['MCD'] = np.zeros((n_records),dtype=np.uint32)

    #-- CryoSat-2 Measurement Group
    #-- Derived from instrument measurement parameters
    Data = {}
    #-- Measured elevation above ellipsoid from retracker 1: packed units (mm, 1e-3 m)
    Data['Elev_1'] = np.zeros((n_records),dtype=np.int32)
    #-- Measured elevation above ellipsoid from retracker 2: packed units (mm, 1e-3 m)
    Data['Elev_2'] = np.zeros((n_records),dtype=np.int32)
    #-- Measured elevation above ellipsoid from retracker 3: packed units (mm, 1e-3 m)
    Data['Elev_3'] = np.zeros((n_records),dtype=np.int32)
    #-- Sigma Zero Backscatter for retracker 1: packed units (1e-2 dB)
    Data['Sig0_1'] = np.zeros((n_records),dtype=np.int32)
    #-- Sigma Zero Backscatter for retracker 2: packed units (1e-2 dB)
    Data['Sig0_2'] = np.zeros((n_records),dtype=np.int32)
    #-- Sigma Zero Backscatter for retracker 3: packed units (1e-2 dB)
    Data['Sig0_3'] = np.zeros((n_records),dtype=np.int32)
    #-- SWH packed units (mm, 1e-3)
    Data['SWH'] = np.zeros((n_records),dtype=np.int32)
    #-- Peakiness: packed units (1e-2)
    Data['Peakiness'] = np.zeros((n_records),dtype=np.int32)
    #-- Retracked range correction for retracker 1: packed units (mm, 1e-3 m)
    Data['Range_1'] = np.zeros((n_records),dtype=np.int32)
    #-- Retracked range correction for retracker 1: packed units (mm, 1e-3 m)
    Data['Range_2'] = np.zeros((n_records),dtype=np.int32)
    #-- Retracked range correction for retracker 3: packed units (mm, 1e-3 m)
    Data['Range_3'] = np.zeros((n_records),dtype=np.int32)
    Data['Spare2'] = np.zeros((n_records),dtype=np.int32)
    Data['Spare3'] = np.zeros((n_records),dtype=np.int32)
    Data['Spare4'] = np.zeros((n_records),dtype=np.int32)
    #-- Retracker 1 quality metric
    Data['Quality_1'] = np.zeros((n_records),dtype=np.int32)
    #-- Retracker 2 quality metric
    Data['Quality_2'] = np.zeros((n_records),dtype=np.int32)
    #-- Retracker 3 quality metric
    Data['Quality_3'] = np.zeros((n_records),dtype=np.int32)
    #-- Retrackers 3-23 output
    Data['Retrack_3'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_4'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_5'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_6'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_7'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_8'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_9'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_10'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_11'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_12'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_13'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_14'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_15'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_16'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_17'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_18'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_19'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_20'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_21'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_22'] = np.zeros((n_records),dtype=np.int32)
    Data['Retrack_23'] = np.zeros((n_records),dtype=np.int32)
    #-- Power echo shape parameter: packed units (dB/100)
    Data['echo_shape'] = np.zeros((n_records),dtype=np.int32)
    #-- Beam behaviour parameter: unitless code number related to
    #-- surface characteristics
    Data['BB_parameter'] = np.zeros((n_records,50),dtype=np.int16)
    #-- Cross track angle: packed units (micro radians)
    Data['X_Track_Angle'] = np.zeros((n_records),dtype=np.int32)
    #-- Cross track angle correction: packed units (micro radians)
    Data['X_Track_Angle_c'] = np.zeros((n_records),dtype=np.int32)
    #-- Leading edge coherence at retrack point 1/1000
    Data['Coherence'] = np.zeros((n_records),dtype=np.int32)
    #-- Interpolated Ocean Height: packed units (mm above ellipsoid)
    Data['Ocean_ht'] = np.zeros((n_records),dtype=np.int32)
    #-- Freeboard: packed units (mm, 1e-3 m)
    #-- -9999 default value indicates computation has not been performed
    Data['Freeboard'] = np.zeros((n_records),dtype=np.int32)
    #-- Surface Height Anomaly: packed units (mm, 1e-3 m)
    Data['SHA'] = np.zeros((n_records),dtype=np.int32)
    #-- Interpolated Surface Height Anomaly: packed units (mm, 1e-3 m)
    Data['SSHA_interp'] = np.zeros((n_records),dtype=np.int32)
    #-- Error in ocean height interpolation: packed units (mm, 1e-3 m)
    Data['SSHA_interp_RMS'] = np.zeros((n_records),dtype=np.uint16)
    #-- Number of forward records interpolated
    Data['SSHA_interp_count_fwd'] = np.zeros((n_records),dtype=np.uint16)
    #-- Number of backward records interpolated
    Data['SSHA_interp_count_bkwd'] = np.zeros((n_records),dtype=np.uint16)
    #-- Distance in time of most forward record interpolated (milli-seconds)
    Data['SSHA_interp_time_fwd'] = np.zeros((n_records),dtype=np.uint16)
    #-- Distance in time of most backward record interpolated (milli-seconds)
    Data['SSHA_interp_time_bkwd'] = np.zeros((n_records),dtype=np.uint16)
    #-- Interpolation error flag
    Data['SSHA_interp_flag'] = np.zeros((n_records),dtype=np.uint16)
    #-- Measurement mode
    Data['Measurement_Mode'] = np.zeros((n_records),dtype=np.uint32)
    #-- Quality flags
    Data['Quality_flag'] = np.zeros((n_records),dtype=np.uint32)
    #-- Retracker flags
    Data['Retracker_flag'] = np.zeros((n_records),dtype=np.uint32)
    #-- Height calculation details
    #-- Specifies what was applied during the height calculation
    Data['Height_status'] = np.zeros((n_records),dtype=np.uint32)
    #-- SAR freeboard status flag
    Data['Freeboard_status'] = np.zeros((n_records),dtype=np.uint32)
    #-- Number of averaged echoes or beams
    Data['N_avg'] = np.zeros((n_records),dtype=np.uint16)
    #-- Wind Speed packed units (mm/s, 1e-3 m/s)
    Data['Wind_speed'] = np.zeros((n_records),dtype=np.uint16)
    Data['Spares1'] = np.zeros((n_records,3),dtype=np.int32)

    #-- CryoSat-2 Auxiliary Data Group
    Auxiliary = {}
    #-- Ice Concentration packed units (%/1000)
    Auxiliary['Ice_conc'] = np.zeros((n_records),dtype=np.int32)
    #-- Snow Depth packed units (mm, 1e-3 m)
    Auxiliary['Snow_depth'] = np.zeros((n_records),dtype=np.int32)
    #-- Snow Density packed units (kg/m^3)
    Auxiliary['Snow_density'] = np.zeros((n_records),dtype=np.int32)
    #-- Discriminator result
    Auxiliary['Discriminator'] = np.zeros((n_records),dtype=np.int32)
    #-- SARin discriminator parameters 1-10
    Auxiliary['SARIN_disc_1'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['SARIN_disc_2'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['SARIN_disc_3'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['SARIN_disc_4'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['SARIN_disc_5'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['SARIN_disc_6'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['SARIN_disc_7'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['SARIN_disc_8'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['SARIN_disc_9'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['SARIN_disc_10'] = np.zeros((n_records),dtype=np.int32)
    #-- Discriminator flags
    Auxiliary['Discrim_flag'] = np.zeros((n_records),dtype=np.uint32)
    #-- Slope model correction (Attitude of echo in micro-degrees)
    Auxiliary['Attitude'] = np.zeros((n_records),dtype=np.int32)
    #-- Slope model correction (Azimuth of echo in micro-degrees)
    Auxiliary['Azimuth'] = np.zeros((n_records),dtype=np.int32)
    #-- Slope doppler correction (mm)
    Auxiliary['Slope_doppler'] = np.zeros((n_records),dtype=np.int32)
    #-- The original latitude of the satellite (micro-degrees)
    Auxiliary['Lat_sat'] = np.zeros((n_records),dtype=np.int32)
    #-- The original longitude of the satellite (micro-degrees)
    Auxiliary['Lon_sat'] = np.zeros((n_records),dtype=np.int32)
    #-- Ambiguity indicator
    Auxiliary['Ambiguity'] = np.zeros((n_records),dtype=np.uint32)
    #-- Mean Sea Surface Model packed units (mm, 1e-3 m)
    Auxiliary['MSS_model'] = np.zeros((n_records),dtype=np.int32)
    #-- Geoid Model packed units (mm, 1e-3 m)
    Auxiliary['Geoid_model'] = np.zeros((n_records),dtype=np.int32)
    #-- ODLE Model packed units (mm, 1e-3 m)
    Auxiliary['ODLE'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['DEM_elev'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['DEM_ID'] = np.zeros((n_records),dtype=np.int32)
    Auxiliary['Spares2'] = np.zeros((n_records,4),dtype=np.int32)

    #-- CryoSat-2 External Corrections Group
    Geometry = {}
    #-- Dry Tropospheric Correction packed units (mm, 1e-3 m)
    Geometry['dryTrop'] = np.zeros((n_records),dtype=np.int32)
    #-- Wet Tropospheric Correction packed units (mm, 1e-3 m)
    Geometry['wetTrop'] = np.zeros((n_records),dtype=np.int32)
    #-- Inverse Barometric Correction packed units (mm, 1e-3 m)
    Geometry['InvBar'] = np.zeros((n_records),dtype=np.int32)
    #-- Delta Inverse Barometric Correction packed units (mm, 1e-3 m)
    Geometry['DAC'] = np.zeros((n_records),dtype=np.int32)
    #-- GIM Ionospheric Correction packed units (mm, 1e-3 m)
    Geometry['Iono_GIM'] = np.zeros((n_records),dtype=np.int32)
    #-- Model Ionospheric Correction packed units (mm, 1e-3 m)
    Geometry['Iono_model'] = np.zeros((n_records),dtype=np.int32)
    #-- Ocean tide Correction packed units (mm, 1e-3 m)
    Geometry['ocTideElv'] = np.zeros((n_records),dtype=np.int32)
    #-- Long period equilibrium ocean tide Correction packed units (mm, 1e-3 m)
    Geometry['lpeTideElv'] = np.zeros((n_records),dtype=np.int32)
    #-- Ocean loading tide Correction packed units (mm, 1e-3 m)
    Geometry['olTideElv'] = np.zeros((n_records),dtype=np.int32)
    #-- Solid Earth tide Correction packed units (mm, 1e-3 m)
    Geometry['seTideElv'] = np.zeros((n_records),dtype=np.int32)
    #-- Geocentric Polar tide Correction packed units (mm, 1e-3 m)
    Geometry['gpTideElv'] = np.zeros((n_records),dtype=np.int32)
    #-- Surface Type: Packed in groups of three bits for each of the 20 records
    Geometry['Surf_type'] = np.zeros((n_records),dtype=np.uint32)
    #-- Corrections Status Flag
    Geometry['Corr_status'] = np.zeros((n_records),dtype=np.uint32)
    #-- Correction Error Flag
    Geometry['Corr_error'] = np.zeros((n_records),dtype=np.uint32)
    #-- Sea State Bias Correction packed units (mm, 1e-3 m)
    Geometry['SSB'] = np.zeros((n_records),dtype=np.int32)
    Geometry['Spares3'] = np.zeros((n_records,2),dtype=np.int32)

    #-- CryoSat-2 Internal Corrections Group
    Instrumental = {}
    #-- Doppler range correction: Radial + slope (mm)
    Instrumental['Doppler_range'] = np.zeros((n_records),dtype=np.int32)
    #-- Instrument Range Correction: t-r antenna (mm)
    Instrumental['TR_inst_range'] = np.zeros((n_records),dtype=np.int32)
    #-- Instrument Range Correction: r-only antenna (mm)
    Instrumental['R_inst_range'] = np.zeros((n_records),dtype=np.int32)
    #-- Instrument Sigma 0 Correction: t-r antenna (dB/100)
    Instrumental['TR_inst_gain'] = np.zeros((n_records),dtype=np.int32)
    #-- Instrument Sigma 0 Correction: r-only (dB/100)
    Instrumental['R_inst_gain'] = np.zeros((n_records),dtype=np.int32)
    #-- Internal Phase Correction (milli-radians)
    Instrumental['Internal_phase'] = np.zeros((n_records),dtype=np.int32)
    #-- External Phase Correction (milli-radians)
    Instrumental['External_phase'] = np.zeros((n_records),dtype=np.int32)
    #-- Noise Power measurement
    Instrumental['Noise_power'] = np.zeros((n_records),dtype=np.int32)
    #-- Phase slope correction (microradians)
    Instrumental['Phase_slope'] = np.zeros((n_records),dtype=np.int32)
    Instrumental['Spares4'] = np.zeros((n_records,2),dtype=np.int32)

    #-- for each record in the CryoSat file
    for r in range(n_records):
        #-- CryoSat-2 Location Group for record r
        Location['Day'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Location['Sec'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Location['Micsec'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Location['USO_Corr'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Location['Mode_ID'][r] = np.fromfile(fid,dtype='>u2',count=1)
        Location['SSC'][r] = np.fromfile(fid,dtype='>u2',count=1)
        Location['Inst_config'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Location['Rec_Count'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Location['Lat'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Location['Lon'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Location['Alt'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Location['Alt_rate'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Location['Sat_velocity'][r,:] = np.fromfile(fid,dtype='>i4',count=3)
        Location['Real_beam'][r,:] = np.fromfile(fid,dtype='>i4',count=3)
        Location['Baseline'][r,:] = np.fromfile(fid,dtype='>i4',count=3)
        Location['ST_ID'][r] = np.fromfile(fid,dtype='>i2',count=1)
        Location['Spare'][r] = np.fromfile(fid,dtype='>i2',count=1)
        Location['Roll'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Location['Pitch'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Location['Yaw'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Location['MCD'][r] = np.fromfile(fid,dtype='>u4',count=1)

        #-- CryoSat-2 Measurement Group for record r
        Data['Elev_1'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Elev_2'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Elev_3'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Sig0_1'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Sig0_2'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Sig0_3'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['SWH'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Peakiness'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Range_1'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Range_2'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Range_3'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Spare2'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Spare3'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Spare4'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Quality_1'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Quality_2'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Quality_3'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_3'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_4'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_5'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_6'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_7'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_8'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_9'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_10'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_11'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_12'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_13'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_14'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_15'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_16'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_17'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_18'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_19'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_20'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_21'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_22'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Retrack_23'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['echo_shape'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['BB_parameter'][r,:] = np.fromfile(fid,dtype='>i2',count=50)
        Data['X_Track_Angle'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['X_Track_Angle_c'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Coherence'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Ocean_ht'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['Freeboard'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['SHA'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['SSHA_interp'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Data['SSHA_interp_RMS'][r] = np.fromfile(fid,dtype='>u2',count=1)
        Data['SSHA_interp_count_fwd'][r] = np.fromfile(fid,dtype='>u2',count=1)
        Data['SSHA_interp_count_bkwd'][r] = np.fromfile(fid,dtype='>u2',count=1)
        Data['SSHA_interp_time_fwd'][r] = np.fromfile(fid,dtype='>u2',count=1)
        Data['SSHA_interp_time_bkwd'][r] = np.fromfile(fid,dtype='>u2',count=1)
        Data['SSHA_interp_flag'][r] = np.fromfile(fid,dtype='>u2',count=1)
        Data['Measurement_Mode'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Data['Quality_flag'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Data['Retracker_flag'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Data['Height_status'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Data['Freeboard_status'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Data['N_avg'][r] = np.fromfile(fid,dtype='>u2',count=1)
        Data['Wind_speed'][r] = np.fromfile(fid,dtype='>u2',count=1)
        Data['Spares1'][r,:] = np.fromfile(fid,dtype='>i4',count=3)

        #-- CryoSat-2 Auxiliary Data Group for record r
        Auxiliary['Ice_conc'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Snow_depth'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Snow_density'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Discriminator'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['SARIN_disc_1'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['SARIN_disc_2'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['SARIN_disc_3'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['SARIN_disc_4'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['SARIN_disc_5'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['SARIN_disc_6'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['SARIN_disc_7'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['SARIN_disc_8'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['SARIN_disc_9'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['SARIN_disc_10'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Discrim_flag'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Auxiliary['Attitude'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Azimuth'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Slope_doppler'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Lat_sat'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Lon_sat'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Ambiguity'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Auxiliary['MSS_model'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Geoid_model'][r] =np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['ODLE'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['DEM_elev'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['DEM_ID'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Auxiliary['Spares2'][r,:] = np.fromfile(fid,dtype='>i4',count=4)

        #-- CryoSat-2 External Corrections Group for record r
        Geometry['dryTrop'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['wetTrop'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['InvBar'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['DAC'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['Iono_GIM'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['Iono_model'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['ocTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['lpeTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['olTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['seTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['gpTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['Surf_type'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Geometry['Corr_status'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Geometry['Corr_error'][r] = np.fromfile(fid,dtype='>u4',count=1)
        Geometry['SSB'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Geometry['Spares3'][r,:] = np.fromfile(fid,dtype='>i4',count=2)

        #-- CryoSat-2 Internal Corrections Group for record r
        Instrumental['Doppler_range'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Instrumental['TR_inst_range'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Instrumental['R_inst_range'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Instrumental['TR_inst_gain'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Instrumental['R_inst_gain'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Instrumental['Internal_phase'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Instrumental['External_phase'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Instrumental['Noise_power'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Instrumental['Phase_slope'][r] = np.fromfile(fid,dtype='>i4',count=1)
        Instrumental['Spares4'][r,:] = np.fromfile(fid,dtype='>i4',count=2)

    #-- Bind all the bits of the l2i_mds together into a single dictionary
    CS_L2I_mds = {}
    CS_L2I_mds['Location'] = Location
    CS_L2I_mds['Data'] = Data
    CS_L2I_mds['Auxiliary'] = Auxiliary
    CS_L2I_mds['Geometry'] = Geometry
    CS_L2I_mds['Instrumental'] = Instrumental
    #-- return the output dictionary
    return CS_L2I_mds

#-- PURPOSE: Initiate L2I MDS variables for CryoSat Baseline D (netCDF4)
def cryosat_baseline_D(full_filename, UNPACK=False):
    #-- open netCDF4 file for reading
    fid = netCDF4.Dataset(os.path.expanduser(full_filename),'r')
    #-- use original unscaled units unless UNPACK=True
    fid.set_auto_scale(UNPACK)
    #-- get dimensions
    time_20_ku = fid.variables['time_20_ku'][:].copy()
    time_cor_01 = fid.variables['time_cor_01'][:].copy()
    n_records = len(time_20_ku)

    #-- CryoSat-2 Location Group
    #-- Time and Orbit Parameters plus Measurement Mode
    Location = {}

    #-- Data Record Time (MDSR Time Stamp)
    #-- Time (seconds since 2000-01-01)
    Location['Time'] = time_20_ku.copy()
    #-- Time: day part
    Location['Day'] = np.array(time_20_ku/86400.0, dtype=np.int32)
    #-- Time: second part
    Location['Second'] = np.array(time_20_ku -
        Location['Day'][:]*86400.0, dtype=np.uint32)
    #-- Time: microsecond part
    Location['Micsec'] = np.array((time_20_ku -
        Location['Day'][:]*86400.0 -
        Location['Second'][:])*1e6, dtype=np.uint32)
    #-- USO correction factor (2-way range)
    Location['USO_Corr'] = fid.variables['uso_cor_20_ku'][:].copy()
    #-- USO correction factor (1-way range)
    Location['USO_Corr_1way'] = fid.variables['uso_cor_applied_20_ku'][:].copy()
    #-- Mode ID
    Location['Mode_ID'] = fid.variables['flag_instr_mode_op_20_ku'][:].copy()
    #-- Mode Flags
    Location['Mode_flags'] = fid.variables['flag_instr_mode_flags_20_ku'][:].copy()
    #-- Platform attitude control mode
    Location['Att_control'] = fid.variables['flag_instr_mode_att_ctrl_20_ku'][:].copy()
    #-- Instrument configuration
    Location['Inst_config'] = fid.variables['flag_instr_conf_rx_flags_20_ku'][:].copy()
    #-- acquisition band
    Location['Inst_band'] = fid.variables['flag_instr_conf_rx_bwdt_20_ku'][:].copy()
    #-- instrument channel
    Location['Inst_channel'] = fid.variables['flag_instr_conf_rx_in_use_20_ku'][:].copy()
    #-- tracking mode
    Location['Tracking_mode'] = fid.variables['flag_instr_conf_rx_trk_mode_20_ku'][:].copy()
    #-- Source sequence counter
    Location['SSC'] = fid.variables['seq_count_20_ku'][:].copy()
    #-- Record Counter
    Location['Rec_Count'] = fid.variables['rec_count_20_ku'][:].copy()
    #-- Lat: packed units (0.1 micro-degree, 1e-7 degrees)
    Location['Lat'] = fid.variables['lat_20_ku'][:].copy()
    #-- Lon: packed units (0.1 micro-degree, 1e-7 degrees)
    Location['Lon'] = fid.variables['lon_20_ku'][:].copy()
    #-- Alt: packed units (mm, 1e-3 m)
    #-- Altitude of COG above reference ellipsoid (interpolated value)
    Location['Alt'] = fid.variables['alt_20_ku'][:].copy()
    #-- Instantaneous altitude rate derived from orbit: packed units (mm/s, 1e-3 m/s)
    Location['Alt_rate'] = fid.variables['orb_alt_rate_20_ku'][:].copy()
    #-- Satellite velocity vector. In ITRF: packed units (mm/s, 1e-3 m/s)
    Location['Sat_velocity'] = fid.variables['sat_vel_vec_20_ku'][:].copy()
    #-- Real beam direction vector. In CRF: packed units (micro-m/s, 1e-6 m/s)
    Location['Real_beam'] = fid.variables['beam_dir_vec_20_ku'][:].copy()
    #-- Interferometer baseline vector. In CRF: packed units (micro-m/s, 1e-6 m/s)
    Location['Baseline'] = fid.variables['inter_base_vec_20_ku'][:].copy()
    #-- Star Tracker ID
    Location['ST_ID'] = fid.variables['flag_instr_conf_rx_str_in_use_20_ku'][:].copy()
    #-- Roll (Derived from star trackers): packed units (0.1 micro-degree, 1e-7 degrees)
    Location['Roll'] = fid.variables['off_nadir_roll_angle_str_20_ku'][:].copy()
    #-- Pitch (Derived from star trackers): packed units (0.1 micro-degree, 1e-7 degrees)
    Location['Pitch'] = fid.variables['off_nadir_pitch_angle_str_20_ku'][:].copy()
    #-- Yaw (Derived from star trackers): packed units (0.1 micro-degree, 1e-7 degrees)
    Location['Yaw'] = fid.variables['off_nadir_yaw_angle_str_20_ku'][:].copy()
    #-- Measurement Confidence Data Flags
    #-- Generally the MCD flags indicate problems when set
    #-- If MCD is 0 then no problems or non-nominal conditions were detected
    #-- Serious errors are indicated by setting bit 31
    Location['MCD'] = fid.variables['flag_mcd_20_ku'][:].copy()


    #-- CryoSat-2 Measurement Group
    #-- Derived from instrument measurement parameters
    Data = {}
    #-- Measured elevation above ellipsoid from retracker 1: packed units (mm, 1e-3 m)
    Data['Elev_1'] = fid.variables['height_1_20_ku'][:].copy()
    #-- Measured elevation above ellipsoid from retracker 2: packed units (mm, 1e-3 m)
    Data['Elev_2'] = fid.variables['height_2_20_ku'][:].copy()
    #-- Measured elevation above ellipsoid from retracker 3: packed units (mm, 1e-3 m)
    Data['Elev_3'] = fid.variables['height_3_20_ku'][:].copy()
    #-- Sigma Zero Backscatter for retracker 1: packed units (1e-2 dB)
    Data['Sig0_1'] = fid.variables['sig0_1_20_ku'][:].copy()
    #-- Sigma Zero Backscatter for retracker 2: packed units (1e-2 dB)
    Data['Sig0_2'] = fid.variables['sig0_2_20_ku'][:].copy()
    #-- Sigma Zero Backscatter for retracker 3: packed units (1e-2 dB)
    Data['Sig0_3'] = fid.variables['sig0_3_20_ku'][:].copy()
    #-- SWH packed units (mm, 1e-3)
    Data['SWH'] = fid.variables['swh_ocean_20_ku'][:].copy()
    #-- Peakiness: packed units (1e-2)
    Data['Peakiness'] = fid.variables['peakiness_20_ku'][:].copy()
    #-- Retracked range correction for retracker 1: packed units (mm, 1e-3 m)
    Data['Range_1'] = fid.variables['range_1_20_ku'][:].copy()
    #-- Retracked range correction for retracker 1: packed units (mm, 1e-3 m)
    Data['Range_2'] = fid.variables['range_2_20_ku'][:].copy()
    #-- Retracked range correction for retracker 3: packed units (mm, 1e-3 m)
    Data['Range_3'] = fid.variables['range_3_20_ku'][:].copy()
    #-- Retracker 1 quality metric
    Data['Quality_1'] = fid.variables['retracker_1_quality_20_ku'][:].copy()
    #-- Retracker 2 quality metric
    Data['Quality_2'] = fid.variables['retracker_2_quality_20_ku'][:].copy()
    #-- Retracker 3 quality metric
    Data['Quality_3'] = fid.variables['retracker_3_quality_20_ku'][:].copy()
    #-- Retrackers 3-23 output
    Data['Retrack_3'] = fid.variables['retracker_output_3_20_ku'][:].copy()
    Data['Retrack_4'] = fid.variables['retracker_output_4_20_ku'][:].copy()
    Data['Retrack_5'] = fid.variables['retracker_output_5_20_ku'][:].copy()
    Data['Retrack_6'] = fid.variables['retracker_output_6_20_ku'][:].copy()
    Data['Retrack_7'] = fid.variables['retracker_output_7_20_ku'][:].copy()
    Data['Retrack_8'] = fid.variables['retracker_output_8_20_ku'][:].copy()
    Data['Retrack_9'] = fid.variables['retracker_output_9_20_ku'][:].copy()
    Data['Retrack_10'] = fid.variables['retracker_output_10_20_ku'][:].copy()
    Data['Retrack_11'] = fid.variables['retracker_output_11_20_ku'][:].copy()
    Data['Retrack_12'] = fid.variables['retracker_output_12_20_ku'][:].copy()
    Data['Retrack_13'] = fid.variables['retracker_output_13_20_ku'][:].copy()
    Data['Retrack_14'] = fid.variables['retracker_output_14_20_ku'][:].copy()
    Data['Retrack_15'] = fid.variables['retracker_output_15_20_ku'][:].copy()
    Data['Retrack_16'] = fid.variables['retracker_output_16_20_ku'][:].copy()
    Data['Retrack_17'] = fid.variables['retracker_output_17_20_ku'][:].copy()
    Data['Retrack_18'] = fid.variables['retracker_output_18_20_ku'][:].copy()
    Data['Retrack_19'] = fid.variables['retracker_output_19_20_ku'][:].copy()
    Data['Retrack_20'] = fid.variables['retracker_output_20_20_ku'][:].copy()
    Data['Retrack_21'] = fid.variables['retracker_output_21_20_ku'][:].copy()
    Data['Retrack_22'] = fid.variables['retracker_output_22_20_ku'][:].copy()
    Data['Retrack_23'] = fid.variables['retracker_output_23_20_ku'][:].copy()
    #-- Power echo shape parameter: packed units (dB/100)
    Data['echo_shape'] = np.zeros((n_records),dtype=np.int32)
    #-- Beam behaviour parameter: unitless code number related to
    #-- surface characteristics
    Data['Beam'] = {}
    #-- Standard Deviation of Gaussian fit to range integrated stack power.
    Data['Beam']['SD'] =  fid.variables['stack_std_20_ku'][:].copy()
    #-- Stack Center: Mean of Gaussian fit to range integrated stack power.
    Data['Beam']['Center'] = fid.variables['stack_centre_20_ku'][:].copy()
    #-- Stack amplitude parameter scaled in dB/100.
    Data['Beam']['Amplitude'] = fid.variables['stack_scaled_amplitude_20_ku'][:].copy()
    #-- 3rd moment: providing the degree of asymmetry of the range integrated
    #-- stack power distribution.
    Data['Beam']['Skewness'] = fid.variables['stack_skewness_20_ku'][:].copy()
    #-- 4th moment: Measure of peakiness of range integrated stack power distribution.
    Data['Beam']['Kurtosis'] = fid.variables['stack_kurtosis_20_ku'][:].copy()
    #-- Stack peakiness computed from the range integrated power of the single look echoes
    Data['Beam']['Peakiness'] = fid.variables['stack_peakiness_20_ku'][:].copy()
    #-- Stack residuals of Gaussian that fits the range integrated power of the single look echoes
    Data['Beam']['RMS'] = fid.variables['stack_gaussian_fitting_residuals_20_ku'][:].copy()
    #-- Standard deviation as a function of boresight angle (microradians)
    Data['Beam']['SD_boresight_angle'] = fid.variables['stack_std_angle_20_ku'][:].copy()
    #-- Stack Center angle as a function of boresight angle (microradians)
    Data['Beam']['Center_boresight_angle'] = fid.variables['stack_centre_angle_20_ku'][:].copy()
    #-- Stack Center angle as a function of look angle (microradians)
    Data['Beam']['Center_look_angle'] = fid.variables['stack_centre_look_angle_20_ku'][:].copy()
    #-- Number of contributing beams in the stack before weighting
    Data['Beam']['Number'] = fid.variables['stack_number_before_weighting_20_ku'][:].copy()
    #-- Number of contributing beams in the stack after weighting
    Data['Beam']['Weighted_Number'] = fid.variables['stack_number_after_weighting_20_ku'][:].copy()
    #-- Cross track angle: packed units (micro radians)
    Data['X_Track_Angle'] = fid.variables['across_track_angle_20_ku'][:].copy()
    #-- Cross track angle correction: packed units (micro radians)
    Data['X_Track_Angle_c'] = fid.variables['across_track_angle_cor_20_ku'][:].copy()
    #-- Leading edge coherence at retrack point 1/1000
    Data['Coherence'] = fid.variables['coherence_20_ku'][:].copy()
    #-- Interpolated Ocean Height: packed units (mm above ellipsoid)
    Data['Ocean_ht'] = np.zeros((n_records),dtype=np.int32)
    #-- Freeboard: packed units (mm, 1e-3 m)
    #-- -9999 default value indicates computation has not been performed
    Data['Freeboard'] = fid.variables['freeboard_20_ku'][:].copy()
    #-- Sea ice Floe height
    Data['Sea_Ice_Lead'] = fid.variables['height_sea_ice_floe_20_ku'][:].copy()
    #-- Sea ice lead height
    Data['Sea_Ice_Floe'] = fid.variables['height_sea_ice_lead_20_ku'][:].copy()
    #-- Surface Height Anomaly: packed units (mm, 1e-3 m)
    Data['SHA'] = fid.variables['ssha_20_ku'][:].copy()
    #-- Interpolated Surface Height Anomaly: packed units (mm, 1e-3 m)
    Data['SSHA_interp'] = fid.variables['ssha_interp_20_ku'][:].copy()
    #-- Error in ocean height interpolation: packed units (mm, 1e-3 m)
    Data['SSHA_interp_RMS'] = fid.variables['ssha_interp_rms_20_ku'][:].copy()
    #-- Number of forward records interpolated
    Data['SSHA_interp_count_fwd'] = fid.variables['ssha_interp_numval_fwd_20_ku'][:].copy()
    #-- Number of backward records interpolated
    Data['SSHA_interp_count_bkwd'] = fid.variables['ssha_interp_numval_back_20_ku'][:].copy()
    #-- Distance in time of most forward record interpolated (milli-seconds)
    Data['SSHA_interp_time_fwd'] = fid.variables['ssha_interp_time_fwd_20_ku'][:].copy()
    #-- Distance in time of most backward record interpolated (milli-seconds)
    Data['SSHA_interp_time_bkwd'] = fid.variables['ssha_interp_time_back_20_ku'][:].copy()
    #-- Interpolation error flag
    Data['SSHA_interp_flag'] = fid.variables['flag_ssha_interp_20_ku'][:].copy()
    #-- Measurement mode
    Data['Measurement_Mode'] = fid.variables['flag_instr_mode_op_20_ku'][:].copy()
    #-- Quality flags
    Data['Quality_flag'] = fid.variables['flag_quality_20_ku'][:].copy()
    #-- Retracker flags
    Data['Retracker_flag'] = fid.variables['flag_retracker_20_ku'][:].copy()
    #-- Height calculation details
    #-- Specifies what was applied during the height calculation
    Data['Height_status'] = fid.variables['flag_height_20_ku'][:].copy()
    #-- SAR freeboard status flag
    Data['Freeboard_status'] = fid.variables['flag_freeboard_20_ku'][:].copy()
    #-- Number of averaged echoes or beams
    Data['N_avg'] = fid.variables['echo_numval_20_ku'][:].copy()
    #-- Wind Speed packed units (mm/s, 1e-3 m/s)
    Data['Wind_speed'] = fid.variables['wind_speed_alt_20_ku'][:].copy()

    #-- CryoSat-2 Auxiliary Data Group
    Auxiliary = {}
    #-- Ice Concentration packed units (%/1000)
    Auxiliary['Ice_conc'] = fid.variables['sea_ice_concentration_20_ku'][:].copy()
    #-- Snow Depth packed units (mm, 1e-3 m)
    Auxiliary['Snow_depth'] = fid.variables['snow_depth_20_ku'][:].copy()
    #-- Snow Density packed units (kg/m^3)
    Auxiliary['Snow_density'] = fid.variables['snow_density_20_ku'][:].copy()
    #-- Discriminator result
    Auxiliary['Discriminator'] = fid.variables['flag_surf_type_class_20_ku'][:].copy()
    #-- SARin discriminator parameters 1-10
    Auxiliary['SARIN_disc_1'] = fid.variables['sarin_output_1_20_ku'][:].copy()
    Auxiliary['SARIN_disc_2'] = fid.variables['sarin_output_2_20_ku'][:].copy()
    Auxiliary['SARIN_disc_3'] = fid.variables['sarin_output_3_20_ku'][:].copy()
    Auxiliary['SARIN_disc_4'] = fid.variables['sarin_output_4_20_ku'][:].copy()
    Auxiliary['SARIN_disc_5'] = fid.variables['sarin_output_5_20_ku'][:].copy()
    Auxiliary['SARIN_disc_6'] = fid.variables['sarin_output_6_20_ku'][:].copy()
    Auxiliary['SARIN_disc_7'] = fid.variables['sarin_output_7_20_ku'][:].copy()
    Auxiliary['SARIN_disc_8'] = fid.variables['sarin_output_8_20_ku'][:].copy()
    Auxiliary['SARIN_disc_9'] = fid.variables['sarin_output_9_20_ku'][:].copy()
    Auxiliary['SARIN_disc_10'] = fid.variables['sarin_output_10_20_ku'][:].copy()
    #-- Discriminator flags
    Auxiliary['Discrim_flag'] = fid.variables['flag_disc_stat_20_ku'][:].copy()
    #-- Slope model correction (Attitude of echo in micro-degrees)
    Auxiliary['Attitude'] = fid.variables['offset_attitude_20_ku'][:].copy()
    #-- Slope model correction (Azimuth of echo in micro-degrees)
    Auxiliary['Azimuth'] = fid.variables['offset_azimuth_20_ku'][:].copy()
    #-- Slope doppler correction (mm)
    Auxiliary['Slope_doppler'] = fid.variables['slope_dop_cor_20_ku'][:].copy()
    #-- The original latitude of the satellite (micro-degrees)
    Auxiliary['Lat_sat'] = np.zeros((n_records),dtype=np.int32)
    #-- The original longitude of the satellite (micro-degrees)
    Auxiliary['Lon_sat'] = np.zeros((n_records),dtype=np.int32)
    #-- Ambiguity indicator
    Auxiliary['Ambiguity'] = fid.variables['flag_sarin_ambiguity_warning_20_ku'][:].copy()
    #-- Mean Sea Surface Model packed units (mm, 1e-3 m)
    Auxiliary['MSS_model'] = fid.variables['mean_sea_surf_sea_ice_20_ku'][:].copy()
    #-- Geoid Model packed units (mm, 1e-3 m)
    Auxiliary['Geoid_model'] = fid.variables['geoid_20_ku'][:].copy()
    #-- ODLE Model packed units (mm, 1e-3 m)
    Auxiliary['ODLE'] = fid.variables['odle_20_ku'][:].copy()
    Auxiliary['DEM_elev'] = fid.variables['dem_height_20_ku'][:].copy()
    Auxiliary['DEM_ID'] = fid.variables['dem_identifier_20_ku'][:].copy()

    #-- CryoSat-2 External Corrections Group (interpolate 1Hz to 20Hz)
    Geometry = {}
    #-- Dry Tropospheric Correction packed units (mm, 1e-3 m)
    Geometry['dryTrop'] = np.zeros((n_records))
    mod_dry_tropo_cor_01 = fid.variables['mod_dry_tropo_cor_01'][:].copy()
    fint=scipy.interpolate.UnivariateSpline(time_cor_01,mod_dry_tropo_cor_01,k=3,s=0)
    Geometry['dryTrop'][:] = fint(time_20_ku)
    #-- Wet Tropospheric Correction packed units (mm, 1e-3 m)
    Geometry['wetTrop'] = np.zeros((n_records))
    mod_wet_tropo_cor_01 = fid.variables['mod_wet_tropo_cor_01'][:].copy()
    fint=scipy.interpolate.UnivariateSpline(time_cor_01,mod_wet_tropo_cor_01,k=3,s=0)
    Geometry['wetTrop'][:] = fint(time_20_ku)
    #-- Inverse Barometric Correction packed units (mm, 1e-3 m)
    Geometry['InvBar'] = np.zeros((n_records))
    inv_bar_cor_01 = fid.variables['inv_bar_cor_01'][:].copy()
    fint=scipy.interpolate.UnivariateSpline(time_cor_01,inv_bar_cor_01,k=3,s=0)
    Geometry['InvBar'][:] = fint(time_20_ku)
    #-- Delta Inverse Barometric Correction packed units (mm, 1e-3 m)
    Geometry['DAC'] = np.zeros((n_records))
    hf_fluct_total_cor_01 = fid.variables['hf_fluct_total_cor_01'][:].copy()
    fint=scipy.interpolate.UnivariateSpline(time_cor_01,hf_fluct_total_cor_01,k=3,s=0)
    Geometry['DAC'][:] = fint(time_20_ku)
    #-- GIM Ionospheric Correction packed units (mm, 1e-3 m)
    Geometry['Iono_GIM'] = np.zeros((n_records))
    iono_cor_gim_01 = fid.variables['iono_cor_gim_01'][:].copy()
    fint=scipy.interpolate.UnivariateSpline(time_cor_01,iono_cor_gim_01,k=3,s=0)
    Geometry['Iono_GIM'][:] = fint(time_20_ku)
    #-- Model Ionospheric Correction packed units (mm, 1e-3 m)
    Geometry['Iono_model'] = np.zeros((n_records))
    iono_cor_01 = fid.variables['iono_cor_01'][:].copy()
    fint=scipy.interpolate.UnivariateSpline(time_cor_01,iono_cor_01,k=3,s=0)
    Geometry['Iono_model'][:] = fint(time_20_ku)
    #-- Ocean tide Correction packed units (mm, 1e-3 m)
    Geometry['ocTideElv'] = np.zeros((n_records))
    ocean_tide_01 = fid.variables['ocean_tide_01'][:].copy()
    fint=scipy.interpolate.UnivariateSpline(time_cor_01,ocean_tide_01,k=3,s=0)
    Geometry['ocTideElv'][:] = fint(time_20_ku)
    #-- Long period equilibrium ocean tide Correction packed units (mm, 1e-3 m)
    Geometry['lpeTideElv'] = np.zeros((n_records))
    ocean_tide_eq_01 = fid.variables['ocean_tide_eq_01'][:].copy()
    fint=scipy.interpolate.UnivariateSpline(time_cor_01,ocean_tide_eq_01,k=3,s=0)
    Geometry['lpeTideElv'][:] = fint(time_20_ku)
    #-- Ocean loading tide Correction packed units (mm, 1e-3 m)
    Geometry['olTideElv'] = np.zeros((n_records))
    load_tide_01 = fid.variables['load_tide_01'][:].copy()
    fint=scipy.interpolate.UnivariateSpline(time_cor_01,load_tide_01,k=3,s=0)
    Geometry['olTideElv'][:] = fint(time_20_ku)
    #-- Solid Earth tide Correction packed units (mm, 1e-3 m)
    Geometry['seTideElv'] = np.zeros((n_records))
    solid_earth_tide_01 = fid.variables['solid_earth_tide_01'][:].copy()
    fint=scipy.interpolate.UnivariateSpline(time_cor_01,solid_earth_tide_01,k=3,s=0)
    Geometry['seTideElv'][:] = fint(time_20_ku)
    #-- Geocentric Polar tide Correction packed units (mm, 1e-3 m)
    Geometry['gpTideElv'] = np.zeros((n_records))
    pole_tide_01 = fid.variables['pole_tide_01'][:].copy()
    fint=scipy.interpolate.UnivariateSpline(time_cor_01,pole_tide_01,k=3,s=0)
    Geometry['gpTideElv'][:] = fint(time_20_ku)
    #-- Surface Type: Packed in groups of three bits for each of the 20 records
    Geometry['Surf_type'] = fid.variables['surf_type_20_ku'][:].copy()
    #-- Corrections Status Flag
    Geometry['Corr_status'] = fid.variables['flag_cor_status_20_ku'][:].copy()
    #-- Correction Error Flag
    Geometry['Corr_error'] = fid.variables['flag_cor_err_20_ku'][:].copy()
    #-- Sea State Bias Correction packed units (mm, 1e-3 m)
    Geometry['SSB'] = fid.variables['sea_state_bias_20_ku'][:].copy()

    #-- CryoSat-2 Internal Corrections Group
    Instrumental = {}
    #-- Doppler range correction: Radial + slope (mm)
    #-- computed for the component of satellite velocity in the nadir direction
    Instrumental['Doppler_range'] = fid.variables['dop_cor_20_ku'][:].copy()
    #-- Value of Doppler Angle for the first single look echo (1e-7 radians)
    Instrumental['Doppler_angle_start'] = fid.variables['dop_angle_start_20_ku'][:].copy()
    #-- Value of Doppler Angle for the last single look echo (1e-7 radians)
    Instrumental['Doppler_angle_stop'] = fid.variables['dop_angle_stop_20_ku'][:].copy()
    #-- Instrument Range Correction: transmit-receive antenna (mm)
    #-- Calibration correction to range on channel 1 computed from CAL1.
    Instrumental['TR_inst_range'] = fid.variables['instr_cor_range_tx_rx_20_ku'][:].copy()
    #-- Instrument Range Correction: receive-only antenna (mm)
    #-- Calibration correction to range on channel 2 computed from CAL1.
    Instrumental['R_inst_range'] = fid.variables['instr_cor_range_rx_20_ku'][:].copy()
    #-- Instrument Gain Correction: transmit-receive antenna (dB/100)
    #-- Calibration correction to gain on channel 1 computed from CAL1
    Instrumental['TR_inst_gain'] = fid.variables['instr_cor_gain_tx_rx_20_ku'][:].copy()
    #-- Instrument Gain Correction: receive-only (dB/100)
    #-- Calibration correction to gain on channel 2 computed from CAL1
    Instrumental['R_inst_gain'] = fid.variables['instr_cor_gain_rx_20_ku'][:].copy()
    #-- Internal Phase Correction (microradians)
    Instrumental['Internal_phase'] = fid.variables['instr_int_ph_cor_20_ku'][:].copy()
    #-- External Phase Correction (microradians)
    Instrumental['External_phase'] = fid.variables['instr_ext_ph_cor_20_ku'][:].copy()
    #-- Noise Power measurement (dB/100)
    Instrumental['Noise_power'] = fid.variables['noise_power_20_ku'][:].copy()
    #-- Phase slope correction (microradians)
    #-- Computed from the CAL-4 packets during the azimuth impulse response
    #-- amplitude (SARIN only). Set from the latest available CAL-4 packet.
    Instrumental['Phase_slope'] = fid.variables['ph_slope_cor_20_ku'][:].copy()

    #-- Bind all the bits of the l2i_mds together into a single dictionary
    CS_L2I_mds = {}
    CS_L2I_mds['Location'] = Location
    CS_L2I_mds['Data'] = Data
    CS_L2I_mds['Auxiliary'] = Auxiliary
    CS_L2I_mds['Geometry'] = Geometry
    CS_L2I_mds['Instrumental'] = Instrumental

    #-- extract global attributes and assign as MPH and SPH metadata
    CS_L2I_mds['METADATA'] = dict(MPH={},SPH={},DSD={})
    #-- MPH attributes
    CS_L2I_mds['METADATA']['MPH']['PRODUCT'] = fid.product_name
    CS_L2I_mds['METADATA']['MPH']['DOI'] = fid.doi
    CS_L2I_mds['METADATA']['MPH']['PROC_STAGE'] =  fid.processing_stage
    CS_L2I_mds['METADATA']['MPH']['REF_DOC'] =  fid.reference_document
    CS_L2I_mds['METADATA']['MPH']['ACQUISITION_STATION'] = fid.acquisition_station
    CS_L2I_mds['METADATA']['MPH']['PROC_CENTER'] = fid.processing_centre
    CS_L2I_mds['METADATA']['MPH']['PROC_TIME'] = fid.creation_time
    CS_L2I_mds['METADATA']['MPH']['SOFTWARE_VER'] = fid.software_version
    CS_L2I_mds['METADATA']['MPH']['SENSING_START'] = fid.sensing_start
    CS_L2I_mds['METADATA']['MPH']['SENSING_STOP'] = fid.sensing_stop
    CS_L2I_mds['METADATA']['MPH']['PHASE'] = fid.phase
    CS_L2I_mds['METADATA']['MPH']['CYCLE'] = fid.cycle_number
    CS_L2I_mds['METADATA']['MPH']['REL_ORBIT'] = fid.rel_orbit_number
    CS_L2I_mds['METADATA']['MPH']['ABS_ORBIT'] = fid.abs_orbit_number
    CS_L2I_mds['METADATA']['MPH']['STATE_VECTOR_TIME'] = fid.state_vector_time
    CS_L2I_mds['METADATA']['MPH']['DELTA_UT1'] = fid.delta_ut1
    CS_L2I_mds['METADATA']['MPH']['X_POSITION'] = fid.x_position
    CS_L2I_mds['METADATA']['MPH']['Y_POSITION'] = fid.y_position
    CS_L2I_mds['METADATA']['MPH']['Z_POSITION'] = fid.z_position
    CS_L2I_mds['METADATA']['MPH']['X_VELOCITY'] = fid.x_velocity
    CS_L2I_mds['METADATA']['MPH']['Y_VELOCITY'] = fid.y_velocity
    CS_L2I_mds['METADATA']['MPH']['Z_VELOCITY'] = fid.z_velocity
    CS_L2I_mds['METADATA']['MPH']['VECTOR_SOURCE'] = fid.vector_source
    CS_L2I_mds['METADATA']['MPH']['LEAP_UTC'] = fid.leap_utc
    CS_L2I_mds['METADATA']['MPH']['LEAP_SIGN'] = fid.leap_sign
    CS_L2I_mds['METADATA']['MPH']['LEAP_ERR'] = fid.leap_err
    CS_L2I_mds['METADATA']['MPH']['PRODUCT_ERR'] = fid.product_err
    #-- SPH attributes
    CS_L2I_mds['METADATA']['SPH']['START_RECORD_TAI_TIME'] = fid.first_record_time
    CS_L2I_mds['METADATA']['SPH']['STOP_RECORD_TAI_TIME'] = fid.last_record_time
    CS_L2I_mds['METADATA']['SPH']['ABS_ORBIT_START'] = fid.abs_orbit_start
    CS_L2I_mds['METADATA']['SPH']['REL_TIME_ASC_NODE_START'] = fid.rel_time_acs_node_start
    CS_L2I_mds['METADATA']['SPH']['ABS_ORBIT_STOP'] = fid.abs_orbit_stop
    CS_L2I_mds['METADATA']['SPH']['REL_TIME_ASC_NODE_STOP'] = fid.rel_time_acs_node_stop
    CS_L2I_mds['METADATA']['SPH']['EQUATOR_CROSS_TIME_UTC'] = fid.equator_cross_time
    CS_L2I_mds['METADATA']['SPH']['EQUATOR_CROSS_LONG'] = fid.equator_cross_long
    CS_L2I_mds['METADATA']['SPH']['ASCENDING_FLAG'] = fid.ascending_flag
    CS_L2I_mds['METADATA']['SPH']['START_LAT'] = fid.first_record_lat
    CS_L2I_mds['METADATA']['SPH']['START_LONG'] = fid.first_record_lon
    CS_L2I_mds['METADATA']['SPH']['STOP_LAT'] = fid.last_record_lat
    CS_L2I_mds['METADATA']['SPH']['STOP_LONG'] = fid.last_record_lon
    CS_L2I_mds['METADATA']['SPH']['L1_PROC_FLAG'] = fid.l1b_proc_flag
    CS_L2I_mds['METADATA']['SPH']['L1_PROCESSING_QUALITY'] = fid.l1b_processing_quality
    CS_L2I_mds['METADATA']['SPH']['L1_PROC_THRESH'] = fid.l1b_proc_thresh
    CS_L2I_mds['METADATA']['SPH']['INSTR_ID'] = fid.instr_id
    CS_L2I_mds['METADATA']['SPH']['LRM_MODE_PERCENT'] = fid.lrm_mode_percent
    CS_L2I_mds['METADATA']['SPH']['SAR_MODE_PERCENT'] = fid.sar_mode_percent
    CS_L2I_mds['METADATA']['SPH']['SARIN_MODE_PERCENT'] = fid.sarin_mode_percent
    CS_L2I_mds['METADATA']['SPH']['OPEN_OCEAN_PERCENT'] = fid.open_ocean_percent
    CS_L2I_mds['METADATA']['SPH']['CLOSE_SEA_PERCENT'] = fid.close_sea_percent
    CS_L2I_mds['METADATA']['SPH']['CONTINENT_ICE_PERCENT'] = fid.continent_ice_percent
    CS_L2I_mds['METADATA']['SPH']['LAND_PERCENT'] = fid.land_percent
    CS_L2I_mds['METADATA']['SPH']['L2_PROD_STATUS'] = fid.l2_prod_status
    CS_L2I_mds['METADATA']['SPH']['L2_PROC_FLAG'] = fid.l2_proc_flag
    CS_L2I_mds['METADATA']['SPH']['L2_PROCESSING_QUALITY'] = fid.l2_processing_quality
    CS_L2I_mds['METADATA']['SPH']['L2_PROC_THRESH'] = fid.l2_proc_thresh
    CS_L2I_mds['METADATA']['SPH']['SIR_CONFIGURATION'] = fid.sir_configuration
    CS_L2I_mds['METADATA']['SPH']['SIR_OP_MODE'] = fid.sir_op_mode
    CS_L2I_mds['METADATA']['SPH']['ORBIT_FILE'] = fid.xref_orbit
    CS_L2I_mds['METADATA']['SPH']['PROC_CONFIG_PARAMS_FILE'] = fid.xref_pconf
    CS_L2I_mds['METADATA']['SPH']['CONSTANTS_FILE'] = fid.xref_constants
    CS_L2I_mds['METADATA']['SPH']['IPF_RA_DATABASE_FILE'] = fid.xref_siral_characterisation
    CS_L2I_mds['METADATA']['SPH']['DORIS_USO_DRIFT_FILE'] = fid.xref_uso
    CS_L2I_mds['METADATA']['SPH']['STAR_TRACKER_ATTREF_FILE'] = fid.xref_star_tracker_attref
    CS_L2I_mds['METADATA']['SPH']['SIRAL_LEVEL_0_FILE'] = fid.xref_siral_l0
    CS_L2I_mds['METADATA']['SPH']['CALIBRATION_TYPE_1_FILE'] = fid.xref_cal1
    CS_L2I_mds['METADATA']['SPH']['SIR_COMPLEX_CAL1_SARIN'] = fid.xref_cal1_sarin
    CS_L2I_mds['METADATA']['SPH']['SCENARIO_FILE'] = fid.xref_orbit_scenario
    CS_L2I_mds['METADATA']['SPH']['CALIBRATION_TYPE_2_FILE'] = fid.xref_cal2
    CS_L2I_mds['METADATA']['SPH']['SURFACE_PRESSURE_FILE'] = fid.xref_surf_pressure
    CS_L2I_mds['METADATA']['SPH']['MEAN_PRESSURE_FILE'] = fid.xref_mean_pressure
    CS_L2I_mds['METADATA']['SPH']['WET_TROPOSPHERE_FILE'] = fid.xref_wet_trop
    CS_L2I_mds['METADATA']['SPH']['U_WIND_FILE'] = fid.xref_u_wind
    CS_L2I_mds['METADATA']['SPH']['V_WIND_FILE'] = fid.xref_v_wind
    CS_L2I_mds['METADATA']['SPH']['METEO_GRID_DEF_FILE'] = fid.xref_meteo
    CS_L2I_mds['METADATA']['SPH']['S1S2_PRESSURE_00H_MAP'] = fid.xref_s1s2_pressure_00h
    CS_L2I_mds['METADATA']['SPH']['S1S2_PRESSURE_06H_MAP'] = fid.xref_s1s2_pressure_06h
    CS_L2I_mds['METADATA']['SPH']['S1S2_PRESSURE_12H_MAP'] = fid.xref_s1s2_pressure_12h
    CS_L2I_mds['METADATA']['SPH']['S1S2_PRESSURE_18H_MAP'] = fid.xref_s1s2_pressure_18h
    CS_L2I_mds['METADATA']['SPH']['S1_TIDE_AMPLITUDE_MAP'] = fid.xref_s1_tide_amplitude
    CS_L2I_mds['METADATA']['SPH']['S1_TIDE_PHASE_MAP'] = fid.xref_s1_tide_phase
    CS_L2I_mds['METADATA']['SPH']['S2_TIDE_AMPLITUDE_MAP'] = fid.xref_s2_tide_amplitude
    CS_L2I_mds['METADATA']['SPH']['S2_TIDE_PHASE_MAP'] = fid.xref_s2_tide_phase
    CS_L2I_mds['METADATA']['SPH']['GPS_IONO_MAP'] = fid.xref_gim
    CS_L2I_mds['METADATA']['SPH']['MODIFIED_DIP_MAP_FILE'] = fid.xref_dip_map
    CS_L2I_mds['METADATA']['SPH']['IONO_COEFFICENTS_FILE'] = fid.xref_iono_cor
    CS_L2I_mds['METADATA']['SPH']['SAI_FILE'] = fid.xref_sai
    CS_L2I_mds['METADATA']['SPH']['OCEAN_TIDE_FILE'] = fid.xref_ocean_tide
    CS_L2I_mds['METADATA']['SPH']['TIDAL_LOADING_FILE'] = fid.xref_tidal_load
    CS_L2I_mds['METADATA']['SPH']['EARTH_TIDE_FILE'] = fid.xref_earth_tide
    CS_L2I_mds['METADATA']['SPH']['POLE_TIDE_FILE'] = fid.xref_pole_location
    CS_L2I_mds['METADATA']['SPH']['SURFACE_TYPE_FILE'] = fid.xref_surf_type
    CS_L2I_mds['METADATA']['SPH']['AUX_MOG2D'] = fid.xref_mog2d
    CS_L2I_mds['METADATA']['SPH']['SIRAL_LEVEL_1B_FILE'] = fid.xref_siral_l1b
    CS_L2I_mds['METADATA']['SPH']['MEAN_SEA_SURFACE_FILE'] = fid.xref_mss
    CS_L2I_mds['METADATA']['SPH']['GEOID_FILE'] = fid.xref_geoid
    CS_L2I_mds['METADATA']['SPH']['ODLE_FILE'] = fid.xref_odle
    #-- mode dependent attributes
    if ('xref_dem' in fid.ncattrs()):
        CS_L2I_mds['METADATA']['SPH']['DEM_MODEL_FILE'] = fid.xref_dem
    if ('xref_sea_ice' in fid.ncattrs()):
        CS_L2I_mds['METADATA']['SPH']['SEA_ICE_FILE'] = fid.xref_sea_ice
    if ('xref_snow_depth' in fid.ncattrs()):
        CS_L2I_mds['METADATA']['SPH']['SNOW_DEPTH_FILE'] = fid.xref_snow_depth

    #-- return the output dictionary
    return CS_L2I_mds

#-- PURPOSE: Get scaling factors for converting unpacked units in binary files
def cryosat_scaling_factors():
    #-- dictionary of scale factors for CryoSat-2 variables
    CS_l2i_scale = {}

    #-- CryoSat-2 Location Group
    #-- Time and Orbit Parameters plus Measurement Mode
    CS_l2i_scale['Location'] = {}
    #-- Time: day part
    CS_l2i_scale['Location']['Day'] = 1.0
    #-- Time: second part
    CS_l2i_scale['Location']['Sec'] = 1.0
    #-- Time: microsecond part
    CS_l2i_scale['Location']['Micsec'] = 1.0
    #-- USO correction factor
    CS_l2i_scale['Location']['USO_Corr'] = 1.0
    #-- Mode ID
    CS_l2i_scale['Location']['Mode_ID'] = 1
    #-- Source sequence counter
    CS_l2i_scale['Location']['SSC'] = 1
    #-- Instrument configuration
    CS_l2i_scale['Location']['Inst_config'] = 1
    #-- Record Counter
    CS_l2i_scale['Location']['Rec_Count'] = 1
    #-- Lat: packed units (0.1 micro-degree, 1e-7 degrees)
    CS_l2i_scale['Location']['Lat'] = 1e-7
    #-- Lon: packed units (0.1 micro-degree, 1e-7 degrees)
    CS_l2i_scale['Location']['Lon'] = 1e-7
    #-- Alt: packed units (mm, 1e-3 m)
    #-- Altitude of COG above reference ellipsoid (interpolated value)
    CS_l2i_scale['Location']['Alt'] = 1e-3
    #-- Instantaneous altitude rate derived from orbit: packed units (mm/s, 1e-3 m/s)
    CS_l2i_scale['Location']['Alt_rate'] = 1e-3
    #-- Satellite velocity vector. In ITRF: packed units (mm/s, 1e-3 m/s)
    CS_l2i_scale['Location']['Sat_velocity'] = 1e-3
    #-- Real beam direction vector. In CRF: packed units (micro-m/s, 1e-6 m/s)
    CS_l2i_scale['Location']['Real_beam'] = 1e-6
    #-- Interferometer baseline vector. In CRF: packed units (micro-m/s, 1e-6 m/s)
    CS_l2i_scale['Location']['Baseline'] = 1e-6
    #-- Star Tracker ID
    CS_l2i_scale['Location']['ST_ID'] = 1
    CS_l2i_scale['Location']['Spare'] = 1
    #-- Roll (Derived from star trackers): packed units (0.1 micro-degree, 1e-7 degrees)
    CS_l2i_scale['Location']['Roll'] = 1e-7
    #-- Pitch (Derived from star trackers): packed units (0.1 micro-degree, 1e-7 degrees)
    CS_l2i_scale['Location']['Pitch'] = 1e-7
    #-- Yaw (Derived from star trackers): packed units (0.1 micro-degree, 1e-7 degrees)
    CS_l2i_scale['Location']['Yaw'] = 1e-7
    #-- Measurement Confidence Data
    CS_l2i_scale['Location']['MCD'] = 1

    #-- CryoSat-2 Measurement Group
    #-- Derived from instrument measurement parameters
    CS_l2i_scale['Data'] = {}
    #-- Measured elevation above ellipsoid from retracker 1: packed units (mm, 1e-3 m)
    CS_l2i_scale['Data']['Elev_1'] = 1e-3
    #-- Measured elevation above ellipsoid from retracker 2: packed units (mm, 1e-3 m)
    CS_l2i_scale['Data']['Elev_2'] = 1e-3
    #-- Measured elevation above ellipsoid from retracker 3: packed units (mm, 1e-3 m)
    CS_l2i_scale['Data']['Elev_3'] = 1e-3
    #-- Sigma Zero Backscatter for retracker 1: packed units (1e-2 dB)
    CS_l2i_scale['Data']['Sig0_1'] = 1e-2
    #-- Sigma Zero Backscatter for retracker 2: packed units (1e-2 dB)
    CS_l2i_scale['Data']['Sig0_2'] = 1e-2
    #-- Sigma Zero Backscatter for retracker 3: packed units (1e-2 dB)
    CS_l2i_scale['Data']['Sig0_3'] = 1e-2
    #-- SWH packed units (mm, 1e-3)
    CS_l2i_scale['Data']['SWH'] = 1e-3
    #-- Peakiness: packed units (1e-2)
    CS_l2i_scale['Data']['Peakiness'] = 1e-2
    #-- Retracked range correction for retracker 1: packed units (mm, 1e-3 m)
    CS_l2i_scale['Data']['Range_1'] = 1e-3
    #-- Retracked range correction for retracker 2: packed units (mm, 1e-3 m)
    CS_l2i_scale['Data']['Range_2'] = 1e-3
    #-- Retracked range correction for retracker 3: packed units (mm, 1e-3 m)
    CS_l2i_scale['Data']['Range_3'] = 1e-3
    #-- Retracked sigma 0 correction for Retracker 1: packed units (1e-2 dB)
    CS_l2i_scale['Data']['Retrack_1_sig0'] = 1e-2
    #-- Retracked sigma 0 correction for Retracker 2: packed units (1e-2 dB)
    CS_l2i_scale['Data']['Retrack_2_sig0'] = 1e-2
    #-- Retracked sigma 0 correction for Retracker 3: packed units (1e-2 dB)
    CS_l2i_scale['Data']['Retrack_3_sig0'] = 1e-2
    #-- Retracker 1 quality metric
    CS_l2i_scale['Data']['Quality_1'] = 1
    #-- Retracker 2 quality metric
    CS_l2i_scale['Data']['Quality_2'] = 1
    #-- Retracker 3 quality metric
    CS_l2i_scale['Data']['Quality_3'] = 1
    #-- Retrackers 3-23 output
    CS_l2i_scale['Data']['Retrack_3'] = 1
    CS_l2i_scale['Data']['Retrack_4'] = 1
    CS_l2i_scale['Data']['Retrack_5'] = 1
    CS_l2i_scale['Data']['Retrack_6'] = 1
    CS_l2i_scale['Data']['Retrack_7'] = 1
    CS_l2i_scale['Data']['Retrack_8'] = 1
    CS_l2i_scale['Data']['Retrack_9'] = 1
    CS_l2i_scale['Data']['Retrack_10'] = 1
    CS_l2i_scale['Data']['Retrack_11'] = 1
    CS_l2i_scale['Data']['Retrack_12'] = 1
    CS_l2i_scale['Data']['Retrack_13'] = 1
    CS_l2i_scale['Data']['Retrack_14'] = 1
    CS_l2i_scale['Data']['Retrack_15'] = 1
    CS_l2i_scale['Data']['Retrack_16'] = 1
    CS_l2i_scale['Data']['Retrack_17'] = 1
    CS_l2i_scale['Data']['Retrack_18'] = 1
    CS_l2i_scale['Data']['Retrack_19'] = 1
    CS_l2i_scale['Data']['Retrack_20'] = 1
    CS_l2i_scale['Data']['Retrack_21'] = 1
    CS_l2i_scale['Data']['Retrack_22'] = 1
    CS_l2i_scale['Data']['Retrack_23'] = 1
    #-- Power echo shape parameter: packed units (dB/100)
    CS_l2i_scale['Data']['echo_shape'] = 1e-2
    #-- Beam behaviour parameter: unitless code number related to
    #-- surface characteristics
    CS_l2i_scale['Data']['BB_parameter'] = 1
    #-- Cross track angle: packed units (micro radians)
    CS_l2i_scale['Data']['X_Track_Angle'] = 1e-6
    #-- Cross track angle correction: packed units (micro radians)
    CS_l2i_scale['Data']['X_Track_Angle_c'] = 1e-6
    #-- Leading edge coherence at retrack point 1/1000
    CS_l2i_scale['Data']['Coherence'] = 1e-3
    #-- Interpolated Ocean Height: packed units (mm above ellipsoid)
    CS_l2i_scale['Data']['Ocean_ht'] = 1e-3
    #-- Freeboard: packed units (mm, 1e-3 m)
    #-- -9999 default value indicates computation has not been performed
    CS_l2i_scale['Data']['Freeboard'] = 1e-3
    #-- Surface Height Anomaly: packed units (mm, 1e-3 m)
    CS_l2i_scale['Data']['SHA'] = 1e-3
    #-- Interpolated Surface Height Anomaly: packed units (mm, 1e-3 m)
    CS_l2i_scale['Data']['SSHA_interp'] = 1e-3
    #-- Error in ocean height interpolation: packed units (mm, 1e-3 m)
    CS_l2i_scale['Data']['SSHA_interp_RMS'] = 1e-3
    #-- Number of forward records interpolated
    CS_l2i_scale['Data']['SSHA_interp_count_fwd'] = 1
    #-- Number of backward records interpolated
    CS_l2i_scale['Data']['SSHA_interp_count_bkwd'] = 1
    #-- Distance in time of most forward record interpolated (milli-seconds)
    CS_l2i_scale['Data']['SSHA_interp_time_fwd'] = 1e-3
    #-- Distance in time of most backward record interpolated (milli-seconds)
    CS_l2i_scale['Data']['SSHA_interp_time_bkwd'] = 1e-3
    #-- Interpolation error flag
    CS_l2i_scale['Data']['SSHA_interp_flag'] = 1
    #-- Measurement mode
    CS_l2i_scale['Data']['Measurement_Mode'] = 1
    #-- Quality flags
    CS_l2i_scale['Data']['Quality_flag'] = 1
    #-- Retracker flags
    CS_l2i_scale['Data']['Retracker_flag'] = 1
    #-- Height calculation details
    #-- Specifies what was applied during the height calculation
    CS_l2i_scale['Data']['Height_status'] = 1
    #-- SAR freeboard status flag
    CS_l2i_scale['Data']['Freeboard_status'] = 1
    #-- Number of averaged echoes or beams
    CS_l2i_scale['Data']['N_avg'] = 1
    #-- Wind Speed packed units (mm/s, 1e-3 m/s)
    CS_l2i_scale['Data']['Wind_speed'] = 1e-3
    CS_l2i_scale['Data']['Spares1'] = 1

    #-- CryoSat-2 Auxiliary Data Group
    CS_l2i_scale['Auxiliary'] = {}
    #-- Ice Concentration packed units (%/1000)
    CS_l2i_scale['Auxiliary']['Ice_conc'] = 1e-3
    #-- Snow Depth packed units (mm, 1e-3 m)
    CS_l2i_scale['Auxiliary']['Snow_depth'] = 1e-3
    #-- Snow Density packed units (kg/m^3)
    CS_l2i_scale['Auxiliary']['Snow_density'] = 1.0
    #-- Discriminator result
    CS_l2i_scale['Auxiliary']['Discriminator'] = 1
    #-- SARin discriminator parameters 1-10
    CS_l2i_scale['Auxiliary']['SARIN_disc_1'] = 1
    CS_l2i_scale['Auxiliary']['SARIN_disc_2'] = 1
    CS_l2i_scale['Auxiliary']['SARIN_disc_3'] = 1
    CS_l2i_scale['Auxiliary']['SARIN_disc_4'] = 1
    CS_l2i_scale['Auxiliary']['SARIN_disc_5'] = 1
    CS_l2i_scale['Auxiliary']['SARIN_disc_6'] = 1
    CS_l2i_scale['Auxiliary']['SARIN_disc_7'] = 1
    CS_l2i_scale['Auxiliary']['SARIN_disc_8'] = 1
    CS_l2i_scale['Auxiliary']['SARIN_disc_9'] = 1
    CS_l2i_scale['Auxiliary']['SARIN_disc_10'] = 1
    #-- Discriminator flags
    CS_l2i_scale['Auxiliary']['Discrim_flag'] = 1
    #-- Slope model correction (Attitude of echo in micro-degrees)
    CS_l2i_scale['Auxiliary']['Attitude'] = 1e-6
    #-- Slope model correction (Azimuth of echo in micro-degrees)
    CS_l2i_scale['Auxiliary']['Azimuth'] = 1e-6
    #-- Slope doppler correction (mm)
    CS_l2i_scale['Auxiliary']['Slope_doppler'] = 1e-3
    #-- The original latitude of the satellite (micro-degrees)
    CS_l2i_scale['Auxiliary']['Lat_sat'] = 1e-6
    #-- The original longitude of the satellite (micro-degrees)
    CS_l2i_scale['Auxiliary']['Lon_sat'] = 1e-6
    #-- Ambiguity indicator
    CS_l2i_scale['Auxiliary']['Ambiguity'] = 1
    #-- Mean Sea Surface standard Model: packed units (mm, 1e-3 m)
    CS_l2i_scale['Auxiliary']['MSS_model'] = 1e-3
    #-- Geoid standard Model: packed units (mm, 1e-3 m)
    CS_l2i_scale['Auxiliary']['Geoid_model'] = 1e-3
    #-- ODLE standard Model: packed units (mm, 1e-3 m)
    CS_l2i_scale['Auxiliary']['ODLE'] = 1e-3
    #-- The interpolated elevation value obtained from the DEM (mm)
    CS_l2i_scale['Auxiliary']['DEM_elev'] = 1e-3
    #-- Identification of DEM used in SARin ambiguity test
    CS_l2i_scale['Auxiliary']['DEM_ID'] = 1
    CS_l2i_scale['Auxiliary']['Spares2'] = 1

    #-- CryoSat-2 External Corrections Group
    CS_l2i_scale['Geometry'] = {}
    #-- Dry Tropospheric Correction packed units (mm, 1e-3 m)
    CS_l2i_scale['Geometry']['dryTrop'] = 1e-3
    #-- Wet Tropospheric Correction packed units (mm, 1e-3 m)
    CS_l2i_scale['Geometry']['wetTrop'] = 1e-3
    #-- Inverse Barometric Correction packed units (mm, 1e-3 m)
    CS_l2i_scale['Geometry']['InvBar'] = 1e-3
    #-- Delta Inverse Barometric Correction packed units (mm, 1e-3 m)
    CS_l2i_scale['Geometry']['DAC'] = 1e-3
    #-- GIM Ionospheric Correction packed units (mm, 1e-3 m)
    CS_l2i_scale['Geometry']['Iono_GIM'] = 1e-3
    #-- Model Ionospheric Correction packed units (mm, 1e-3 m)
    CS_l2i_scale['Geometry']['Iono_model'] = 1e-3
    #-- Ocean tide Correction packed units (mm, 1e-3 m)
    CS_l2i_scale['Geometry']['ocTideElv'] = 1e-3
    #-- Long period equilibrium ocean tide Correction packed units (mm, 1e-3 m)
    CS_l2i_scale['Geometry']['lpeTideElv'] = 1e-3
    #-- Ocean loading tide Correction packed units (mm, 1e-3 m)
    CS_l2i_scale['Geometry']['olTideElv'] = 1e-3
    #-- Solid Earth tide Correction packed units (mm, 1e-3 m)
    CS_l2i_scale['Geometry']['seTideElv'] = 1e-3
    #-- Geocentric Polar tide Correction packed units (mm, 1e-3 m)
    CS_l2i_scale['Geometry']['gpTideElv'] = 1e-3
    #-- Surface Type: Packed in groups of three bits for each of the 20 records
    CS_l2i_scale['Geometry']['Surf_type'] = 1
    #-- Corrections Status Flag
    CS_l2i_scale['Geometry']['Corr_status'] = 1
    #-- Correction Error Flag
    CS_l2i_scale['Geometry']['Corr_error'] = 1
    #-- Sea State Bias Correction packed units (mm, 1e-3 m)
    CS_l2i_scale['Geometry']['SSB'] = 1e-3
    CS_l2i_scale['Geometry']['Spares3'] = 1

    #-- CryoSat-2 Internal Corrections Group
    CS_l2i_scale['Instrumental'] = {}
    #-- Doppler range correction: Radial + slope (mm)
    CS_l2i_scale['Instrumental']['Doppler_range'] = 1e-3
    #-- Instrument Range Correction: t-r antenna (mm)
    CS_l2i_scale['Instrumental']['TR_inst_range'] = 1e-3
    #-- Instrument Range Correction: r-only antenna (mm)
    CS_l2i_scale['Instrumental']['R_inst_range'] = 1e-3
    #-- Instrument Sigma 0 Correction: t-r antenna (dB/100)
    CS_l2i_scale['Instrumental']['TR_inst_gain'] = 1e-2
    #-- Instrument Sigma 0 Correction: r-only (dB/100)
    CS_l2i_scale['Instrumental']['R_inst_gain'] = 1e-2
    #-- Internal Phase Correction (milli-radians)
    CS_l2i_scale['Instrumental']['Internal_phase'] = 1e-3
    #-- External Phase Correction (milli-radians)
    CS_l2i_scale['Instrumental']['External_phase'] = 1e-3
    #-- Noise Power measurement
    CS_l2i_scale['Instrumental']['Noise_power'] = 1
    #-- Phase slope correction (microradians)
    CS_l2i_scale['Instrumental']['Phase_slope'] = 1e-6
    CS_l2i_scale['Instrumental']['Spares4'] = 1

    #-- return the scaling factors
    return CS_l2i_scale

#-- PURPOSE: Read ASCII Main Product Header (MPH) block from an ESA PDS file
def read_MPH(full_filename):
    #-- read input data file
    with open(os.path.expanduser(full_filename), 'rb') as fid:
        file_contents = fid.read().splitlines()

    #-- Define constant values associated with PDS file formats
    #-- number of text lines in standard MPH
    n_MPH_lines = 41
    #-- check that first line of header matches PRODUCT
    if not bool(re.match(br'PRODUCT\=\"(.*)(?=\")',file_contents[0])):
        raise IOError('File does not start with a valid PDS MPH')
    #-- read MPH header text
    s_MPH_fields = {}
    for i in range(n_MPH_lines):
        #-- use regular expression operators to read headers
        if bool(re.match(br'(.*?)\=\"(.*)(?=\")',file_contents[i])):
            #-- data fields within quotes
            field,value=re.findall(br'(.*?)\=\"(.*)(?=\")',file_contents[i]).pop()
            s_MPH_fields[field.decode('utf-8')] = value.decode('utf-8').rstrip()
        elif bool(re.match(br'(.*?)\=(.*)',file_contents[i])):
            #-- data fields without quotes
            field,value=re.findall(br'(.*?)\=(.*)',file_contents[i]).pop()
            s_MPH_fields[field.decode('utf-8')] = value.decode('utf-8').rstrip()

    #-- Return block name array to calling function
    return s_MPH_fields

#-- PURPOSE: Read ASCII Specific Product Header (SPH) block from a PDS file
def read_SPH(full_filename,j_sph_size):
    #-- read input data file
    with open(os.path.expanduser(full_filename), 'rb') as fid:
        file_contents = fid.read().splitlines()

    #-- Define constant values associated with PDS file formats
    #-- number of text lines in standard MPH
    n_MPH_lines = 41
    #-- compile regular expression operator for reading headers
    rx = re.compile(br'(.*?)\=\"?(.*)',re.VERBOSE)
    #-- check first line of header matches SPH_DESCRIPTOR
    if not bool(re.match(br'SPH\_DESCRIPTOR\=',file_contents[n_MPH_lines+1])):
        raise IOError('File does not have a valid PDS DSD')
    #-- read SPH header text (no binary control characters)
    s_SPH_lines = [li for li in file_contents[n_MPH_lines+1:] if rx.match(li)
        and not re.search(br'[^\x20-\x7e]+',li)]

    #-- extract SPH header text
    s_SPH_fields = {}
    c = 0
    while (c < len(s_SPH_lines)):
        #-- check if line is within DS_NAME portion of SPH header
        if bool(re.match(br'DS_NAME',s_SPH_lines[c])):
            #-- add dictionary for DS_NAME
            field,value=re.findall(br'(.*?)\=\"(.*)(?=\")',s_SPH_lines[c]).pop()
            key = value.decode('utf-8').rstrip()
            s_SPH_fields[key] = {}
            for line in s_SPH_lines[c+1:c+7]:
                if bool(re.match(br'(.*?)\=\"(.*)(?=\")',line)):
                    #-- data fields within quotes
                    dsfield,dsvalue=re.findall(br'(.*?)\=\"(.*)(?=\")',line).pop()
                    s_SPH_fields[key][dsfield.decode('utf-8')] = dsvalue.decode('utf-8').rstrip()
                elif bool(re.match(br'(.*?)\=(.*)',line)):
                    #-- data fields without quotes
                    dsfield,dsvalue=re.findall(br'(.*?)\=(.*)',line).pop()
                    s_SPH_fields[key][dsfield.decode('utf-8')] = dsvalue.decode('utf-8').rstrip()
            #-- add 6 to counter to go to next entry
            c += 6
        #-- use regular expression operators to read headers
        elif bool(re.match(br'(.*?)\=\"(.*)(?=\")',s_SPH_lines[c])):
            #-- data fields within quotes
            field,value=re.findall(br'(.*?)\=\"(.*)(?=\")',s_SPH_lines[c]).pop()
            s_SPH_fields[field.decode('utf-8')] = value.decode('utf-8').rstrip()
        elif bool(re.match(br'(.*?)\=(.*)',s_SPH_lines[c])):
            #-- data fields without quotes
            field,value=re.findall(br'(.*?)\=(.*)',s_SPH_lines[c]).pop()
            s_SPH_fields[field.decode('utf-8')] = value.decode('utf-8').rstrip()
        #-- add 1 to counter to go to next line
        c += 1

    #-- Return block name array to calling function
    return s_SPH_fields

#-- PURPOSE: Read ASCII Data Set Descriptors (DSD) block from a PDS file
def read_DSD(full_filename):
    #-- read input data file
    with open(os.path.expanduser(full_filename), 'rb') as fid:
        file_contents = fid.read().splitlines()

    #-- Define constant values associated with PDS file formats
    #-- number of text lines in standard MPH
    n_MPH_lines = 41
    #-- number of text lines in a DSD header
    n_DSD_lines = 8

    #-- Level-2 CryoSat DS_NAMES within files
    regex_patterns = []
    regex_patterns.append(br'DS_NAME\="SIR_LRM_L2(_I)?[\s+]*"')
    regex_patterns.append(br'DS_NAME\="SIR_LRMIL2[\s+]*"')
    regex_patterns.append(br'DS_NAME\="SIR_SAR_L2(A|B)?(_I)?[\s+]*"')
    regex_patterns.append(br'DS_NAME\="SIR_SARIL2(A|B)?[\s+]*"')
    regex_patterns.append(br'DS_NAME\="SIR_FDM_L2[\s+]*"')
    regex_patterns.append(br'DS_NAME\="SIR_SIN_L2(_I)?[\s+]*"')
    regex_patterns.append(br'DS_NAME\="SIR_SINIL2[\s+]*"')
    regex_patterns.append(br'DS_NAME\="SIR_SID_L2(_I)?[\s+]*"')
    regex_patterns.append(br'DS_NAME\="SIR_SIDIL2[\s+]*"')
    regex_patterns.append(br'DS_NAME\="SIR_GDR_2(A|B|_)?[\s+]*"')
    #-- find the DSD starting line within the SPH header
    c = 0
    Flag = False
    while ((Flag is False) and (c < len(regex_patterns))):
        #-- find indice within
        indice = [i for i,line in enumerate(file_contents[n_MPH_lines+1:]) if
            re.search(regex_patterns[c],line)]
        if indice:
            Flag = True
        else:
            c+=1
    #-- check that valid indice was found within header
    if not indice:
        raise IOError('Can not find correct DSD field')

    #-- extract s_DSD_fields info
    DSD_START = n_MPH_lines + indice[0] + 1
    s_DSD_fields = {}
    for i in range(DSD_START,DSD_START+n_DSD_lines):
        #-- use regular expression operators to read headers
        if bool(re.match(br'(.*?)\=\"(.*)(?=\")',file_contents[i])):
            #-- data fields within quotes
            field,value=re.findall(br'(.*?)\=\"(.*)(?=\")',file_contents[i]).pop()
            s_DSD_fields[field.decode('utf-8')] = value.decode('utf-8').rstrip()
        elif bool(re.match(br'(.*?)\=(.*)',file_contents[i])):
            #-- data fields without quotes
            field,value=re.findall(br'(.*?)\=(.*)',file_contents[i]).pop()
            s_DSD_fields[field.decode('utf-8')] = value.decode('utf-8').rstrip()

    #-- Return block name array to calling function
    return s_DSD_fields

#-- PURPOSE: read CryoSat Level-2 Intermediate data
def read_cryosat_L2I(full_filename, VERBOSE=False):
    #-- file basename and file extension of input file
    fileBasename,fileExtension=os.path.splitext(os.path.basename(full_filename))

    #-- CryoSat file class
    #-- OFFL (Off Line Processing/Systematic)
    #-- NRT_ (Near Real Time)
    #-- RPRO (ReProcessing)
    #-- TEST (Testing)
    #-- LTA_ (Long Term Archive)
    regex_class = r'OFFL|NRT_|RPRO|TEST|LTA_'
    #-- CryoSat mission products
    #-- SIR_LRM_2 L2 Product from Low Resolution Mode Processing
    #-- SIR_FDM_2 L2 Product from Fast Delivery Marine Mode Processing
    #-- SIR_SIN_2 L2 Product from SAR Interferometric Processing
    #-- SIR_SID_2 L2 Product from SIN Degraded Processing
    #-- SIR_SAR_2 L2 Product from SAR Processing
    #-- SIR_GDR_2 L2 Consolidated Product
    #-- SIR_LRMI2 In-depth L2 Product from LRM Processing
    #-- SIR_SINI2 In-depth L2 Product from SIN Processing
    #-- SIR_SIDI2 In-depth L2 Product from SIN Degraded Process.
    #-- SIR_SARI2 In-depth L2 Product from SAR Processing
    regex_products = (r'SIR_LRM_2|SIR_FDM_2|SIR_SIN_2|SIR_SID_2|'
        r'SIR_SAR_2|SIR_GDR_2|SIR_LRMI2|SIR_SINI2|SIR_SIDI2|SIR_SARI2')
    #-- CRYOSAT LEVEL-2 PRODUCTS NAMING RULES
    #-- Mission Identifier
    #-- File Class
    #-- File Product
    #-- Validity Start Date and Time
    #-- Validity Stop Date and Time
    #-- Baseline Identifier
    #-- Version Number
    regex_pattern = r'(.*?)_({0})_({1})__(\d+T?\d+)_(\d+T?\d+)_(.*?)(\d+)'
    rx = re.compile(regex_pattern.format(regex_class,regex_products),re.VERBOSE)
    #-- extract file information from filename
    MI,CLASS,PRODUCT,START,STOP,BASELINE,VERSION=rx.findall(fileBasename).pop()

    #-- check if input file is original binary *.DBL or new netCDF4 *.nc format
    if (fileExtension == '.nc'):
        print(fileBasename) if VERBOSE else None
        CS_L2I_mds = cryosat_baseline_D(full_filename, UNPACK=False)
    elif (fileExtension == '.DBL'):
        #-- Record sizes
        CS_L2I_MDS_REC_SIZE = 556
        CS_L2I_BC_MDS_REC_SIZE = 664
        CS_L2I_C_MDS_REC_SIZE = 664
        #-- check baseline from file to set i_record_size and allocation function
        if (BASELINE == 'C'):
            i_record_size = CS_L2I_C_MDS_REC_SIZE
            read_cryosat_variables = cryosat_baseline_C
        elif (BASELINE == 'BC'):
            i_record_size = CS_L2I_BC_MDS_REC_SIZE
            read_cryosat_variables = cryosat_baseline_BC
        else:
            i_record_size = CS_L2I_MDS_REC_SIZE
            read_cryosat_variables = cryosat_baseline_AB

        #-- read the input file to get file information
        fid = os.open(os.path.expanduser(full_filename),os.O_RDONLY)
        file_info = os.fstat(fid)
        os.close(fid)

        #-- num DSRs from SPH
        j_num_DSR = np.int32(file_info.st_size//i_record_size)
        #-- print file information
        if VERBOSE:
            print(fileBasename)
            print('{0:d} {1:d} {2:d}'.format(j_num_DSR,file_info.st_size,i_record_size))
            #-- Check if MPH/SPH/DSD headers
            if (j_num_DSR*i_record_size == file_info.st_size):
                print('No Header on file')
                print('The number of DSRs is: {0:d}'.format(j_num_DSR))
            else:
                print('Header on file')

        #-- Check if MPH/SPH/DSD headers
        if (j_num_DSR*i_record_size != file_info.st_size):
            #-- If there are MPH/SPH/DSD headers
            s_MPH_fields = read_MPH(full_filename)
            j_sph_size = np.int32(re.findall('[-+]?\d+',s_MPH_fields['SPH_SIZE']).pop())
            s_SPH_fields = read_SPH(full_filename,j_sph_size)
            #-- extract information from DSD fields
            s_DSD_fields = read_DSD(full_filename)
            #-- extract DS_OFFSET
            j_DS_start = np.int32(re.findall('[-+]?\d+',s_DSD_fields['DS_OFFSET']).pop())
            #-- extract number of DSR in the file
            j_num_DSR = np.int32(re.findall('[-+]?\d+',s_DSD_fields['NUM_DSR']).pop())
            #-- check the record size
            j_DSR_size = np.int32(re.findall('[-+]?\d+',s_DSD_fields['DSR_SIZE']).pop())
            #--  minimum size is start of the read plus number of records to read
            j_check_size = j_DS_start +(j_DSR_size*j_num_DSR)
            if VERBOSE:
                print('The offset of the DSD is: {0:d} bytes'.format(j_DS_start))
                print('The number of DSRs is {0:d}'.format(j_num_DSR))
                print('The size of the DSR is {0:d}'.format(j_DSR_size))
            #-- check if invalid file size
            if (j_check_size > file_info.st_size):
                raise IOError('File size error')
            #-- extract binary data from input CryoSat data file (skip headers)
            fid = open(os.path.expanduser(full_filename), 'rb')
            cryosat_header = fid.read(j_DS_start)
            #-- iterate through CryoSat file and fill output variables
            CS_L2I_mds = read_cryosat_variables(fid,i_record_size,j_num_DSR)
            #-- add headers to output dictionary as METADATA
            CS_L2I_mds['METADATA'] = {}
            CS_L2I_mds['METADATA']['MPH'] = s_MPH_fields
            CS_L2I_mds['METADATA']['SPH'] = s_SPH_fields
            CS_L2I_mds['METADATA']['DSD'] = s_DSD_fields
            #-- close the input CryoSat binary file
            fid.close()
        else:
            #-- If there are not MPH/SPH/DSD headers
            #-- extract binary data from input CryoSat data file
            fid = open(os.path.expanduser(full_filename), 'rb')
            #-- iterate through CryoSat file and fill output variables
            CS_L2I_mds = read_cryosat_variables(fid,i_record_size,j_num_DSR)
            #-- close the input CryoSat binary file
            fid.close()

    #-- return the data and headers
    return CS_L2I_mds
