#!/usr/bin/env python
u"""
read_cryosat_L2.py
Written by Tyler Sutterley (09/2019)

Reads CryoSat Level-2 data products from baselines A, B and C
Reads CryoSat Level-2 netCDF4 data products from baseline D
Supported CryoSat Modes: LRM, SAR, SARin, FDM, SID, GDR

INPUTS:
	full_filename: full path of CryoSat .DBL or .nc file

OUTPUTS:
	Data_1Hz: Time and Orbit Parameters
	Corrections: Elevation Corrections and Flags
	Data_20Hz: Geolocation and Elevation Measurements with Quality Parameters
	METADATA: MPH, SPH and DSD Header data

PYTHON DEPENDENCIES:
	numpy: Scientific Computing Tools For Python
		http://www.numpy.org
		http://www.scipy.org/NumPy_for_Matlab_Users
	netCDF4: Python interface to the netCDF C library
	 	https://unidata.github.io/netcdf4-python/netCDF4/index.html

UPDATE HISTORY:
Updated 09/2019: added netCDF4 read function for baseline D
	will output with same variable names and units as the binary read functions
	output 20Hz data as masked arrays for all baselines
Updated 08/2019: generalize regular expression patterns in read_DSD function
Updated 10/2018: updated header read functions for python3
Updated 11/2016: added Abs_Orbit and Ascending_Flg to Data_1Hz outputs
	Abs_Orbit should be same as in read_cryosat_ground_tracks.py
	Ascending_Flg can use in surface regression fits following McMillan (2014)
Updated 05/2016: using __future__ print and division functions
Written 03/2016
"""
from __future__ import print_function
from __future__ import division

import os
import re
import netCDF4
import numpy as np

#-- PURPOSE: Initiate L2 MDS variables for CryoSat Baselines A and B
def cryosat_baseline_AB(fid,record_size,n_records):
	#-- CryoSat-2 1 Hz data fields (Location Group)
	#-- Time and Orbit Parameters plus Measurement Mode
	Data_1Hz = {}
	#-- Time: day part
	Data_1Hz['Day'] = np.zeros((n_records),dtype=np.int32)
	#-- Time: second part
	Data_1Hz['Second'] = np.zeros((n_records),dtype=np.int32)
	#-- Time: microsecond part
	Data_1Hz['Micsec'] = np.zeros((n_records),dtype=np.int32)
	#-- SIRAL mode
	Data_1Hz['Siral_mode'] = np.zeros((n_records),dtype=np.uint64)
	#-- Lat_1Hz: packed units (0.1 micro-degree, 1e-7 degrees)
	Data_1Hz['Lat_1Hz'] = np.zeros((n_records),dtype=np.int32)
	#-- Lon_1Hz: packed units (0.1 micro-degree, 1e-7 degrees)
	Data_1Hz['Lon_1Hz'] = np.zeros((n_records),dtype=np.int32)
	#-- Alt_1Hz: packed units (mm, 1e-3 m)
	#-- Altitude of COG above reference ellipsoid (interpolated value)
	Data_1Hz['Alt_1Hz'] = np.zeros((n_records),dtype=np.int32)
	#-- Mispointing: packed units (millidegrees, 1e-3 degrees)
	Data_1Hz['Mispointing'] = np.zeros((n_records),dtype=np.int16)
	#-- Number of valid records in the block of twenty that contain data
	#-- Last few records of the last block of a dataset may be blank blocks
	#-- inserted to bring the file up to a multiple of twenty.
	Data_1Hz['N_valid'] = np.zeros((n_records),dtype=np.int16)

	#-- CryoSat-2 geophysical corrections (External Corrections Group)
	Corrections = {}
	#-- Dry Tropospheric Correction packed units (mm, 1e-3 m)
	Corrections['dryTrop'] = np.zeros((n_records),dtype=np.int16)
	#-- Wet Tropospheric Correction packed units (mm, 1e-3 m)
	Corrections['wetTrop'] = np.zeros((n_records),dtype=np.int16)
	#-- Inverse Barometric Correction packed units (mm, 1e-3 m)
	Corrections['InvBar'] = np.zeros((n_records),dtype=np.int16)
	#-- Dynamic Atmosphere Correction packed units (mm, 1e-3 m)
	Corrections['DAC'] = np.zeros((n_records),dtype=np.int16)
	#-- Ionospheric Correction packed units (mm, 1e-3 m)
	Corrections['Iono'] = np.zeros((n_records),dtype=np.int16)
	#-- Sea State Bias Correction packed units (mm, 1e-3 m)
	Corrections['SSB'] = np.zeros((n_records),dtype=np.int16)
	#-- Ocean tide Correction packed units (mm, 1e-3 m)
	Corrections['ocTideElv'] = np.zeros((n_records),dtype=np.int16)
	#-- Long period equilibrium ocean tide Correction packed units (mm, 1e-3 m)
	Corrections['lpeTideElv'] = np.zeros((n_records),dtype=np.int16)
	#-- Ocean loading tide Correction packed units (mm, 1e-3 m)
	Corrections['olTideElv'] = np.zeros((n_records),dtype=np.int16)
	#-- Solid Earth tide Correction packed units (mm, 1e-3 m)
	Corrections['seTideElv'] = np.zeros((n_records),dtype=np.int16)
	#-- Geocentric Polar tide Correction packed units (mm, 1e-3 m)
	Corrections['gpTideElv'] = np.zeros((n_records),dtype=np.int16)
	Corrections['Spare1'] = np.zeros((n_records),dtype=np.int16)
	#-- Surface Type: Packed in groups of three bits for each of the 20 records
	Corrections['Surf_type'] = np.zeros((n_records),dtype=np.uint64)
	#-- Mean Sea Surface or Geoid packed units (mm, 1e-3 m)
	Corrections['MSS_Geoid'] = np.zeros((n_records),dtype=np.int32)
	#-- Ocean Depth/Land Elevation Model (ODLE) packed units (mm, 1e-3 m)
	Corrections['ODLE'] = np.zeros((n_records),dtype=np.int32)
	#-- Ice Concentration packed units (%/100)
	Corrections['Ice_conc'] = np.zeros((n_records),dtype=np.int16)
	#-- Snow Depth packed units (mm, 1e-3 m)
	Corrections['Snow_depth'] = np.zeros((n_records),dtype=np.int16)
	#-- Snow Density packed units (kg/m^3)
	Corrections['Snow_density'] = np.zeros((n_records),dtype=np.int16)
	Corrections['Spare2'] = np.zeros((n_records),dtype=np.int16)
	#-- Corrections Status Flag
	Corrections['C_status'] = np.zeros((n_records),dtype=np.uint32)
	#-- Significant Wave Height (SWH) packed units (mm, 1e-3)
	Corrections['SWH'] = np.zeros((n_records),dtype=np.int16)
	#-- Wind Speed packed units (mm/s, 1e-3 m/s)
	Corrections['Wind_speed'] = np.zeros((n_records),dtype=np.uint16)
	Corrections['Spare3'] = np.zeros((n_records),dtype=np.int16)
	Corrections['Spare4'] = np.zeros((n_records),dtype=np.int16)
	Corrections['Spare5'] = np.zeros((n_records),dtype=np.int16)
	Corrections['Spare6'] = np.zeros((n_records),dtype=np.int16)

	#-- CryoSat-2 20 Hz data fields (Measurement Group)
	#-- Derived from instrument measurement parameters
	n_blocks = 20
	Data_20Hz = {}
	#-- Delta between the timestamps for 20Hz record and the 1Hz record
	#-- D_time_mics packed units (microseconds)
	Data_20Hz['D_time_mics'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
	Data_20Hz['D_time_mics'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- Lat: packed units (0.1 micro-degree, 1e-7 degrees)
	Data_20Hz['Lat'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
	Data_20Hz['Lat'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- Lon: packed units (0.1 micro-degree, 1e-7 degrees)
	Data_20Hz['Lon'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
	Data_20Hz['Lon'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- Measured elevation above ellipsoid from retracker: packed units (mm, 1e-3 m)
	Data_20Hz['Elev'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
	Data_20Hz['Elev'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- Interpolated Sea Surface Height Anomaly: packed units (mm, 1e-3 m)
	Data_20Hz['SSHA_interp'] = np.ma.zeros((n_records,n_blocks),dtype=np.int16)
	Data_20Hz['SSHA_interp'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- Interpolated Sea Surface Height measurement count
	Data_20Hz['SSHA_num'] = np.ma.zeros((n_records,n_blocks),dtype=np.int16)
	Data_20Hz['SSHA_num'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- Interpolation quality estimate RSS: packed units (mm, 1e-3 m)
	Data_20Hz['SSHA_qual'] = np.ma.zeros((n_records,n_blocks),dtype=np.int16)
	Data_20Hz['SSHA_qual'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- Sigma Zero Backscatter for retracker: packed units (1e-2 dB)
	Data_20Hz['Sig0'] = np.ma.zeros((n_records,n_blocks),dtype=np.int16)
	Data_20Hz['Sig0'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- Peakiness: packed units (1e-2)
	Data_20Hz['Peakiness'] = np.ma.zeros((n_records,n_blocks),dtype=np.uint16)
	Data_20Hz['Peakiness'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- Freeboard: packed units (mm, 1e-3 m)
	#-- -9999 default value indicates computation has not been performed
	Data_20Hz['Freeboard'] = np.ma.zeros((n_records,n_blocks),dtype=np.int16)
	Data_20Hz['Freeboard'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- Number of averaged echoes or beams
	Data_20Hz['N_avg'] = np.ma.zeros((n_records,n_blocks),dtype=np.int16)
	Data_20Hz['N_avg'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	Data_20Hz['Spare1'] = np.ma.zeros((n_records,n_blocks),dtype=np.int16)
	Data_20Hz['Spare1'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- Quality flags
	Data_20Hz['Quality_Flg'] = np.ma.zeros((n_records,n_blocks),dtype=np.uint32)
	Data_20Hz['Quality_Flg'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	Data_20Hz['Spare2'] = np.ma.zeros((n_records,n_blocks),dtype=np.int16)
	Data_20Hz['Spare2'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	Data_20Hz['Spare3'] = np.ma.zeros((n_records,n_blocks),dtype=np.int16)
	Data_20Hz['Spare3'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	Data_20Hz['Spare4'] = np.ma.zeros((n_records,n_blocks),dtype=np.int16)
	Data_20Hz['Spare4'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	Data_20Hz['Spare5'] = np.ma.zeros((n_records,n_blocks),dtype=np.int16)
	Data_20Hz['Spare5'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- for each record in the CryoSat file
	for r in range(n_records):
		#-- CryoSat-2 Location Group for record r
		Data_1Hz['Day'][r] = np.fromfile(fid,dtype='>i4',count=1)
		Data_1Hz['Second'][r] = np.fromfile(fid,dtype='>i4',count=1)
		Data_1Hz['Micsec'][r] = np.fromfile(fid,dtype='>i4',count=1)
		Data_1Hz['Siral_mode'][r] = np.fromfile(fid,dtype='>u8',count=1)
		Data_1Hz['Lat_1Hz'][r] = np.fromfile(fid,dtype='>i4',count=1)
		Data_1Hz['Lon_1Hz'][r] = np.fromfile(fid,dtype='>i4',count=1)
		Data_1Hz['Alt_1Hz'][r] = np.fromfile(fid,dtype='>i4',count=1)
		Data_1Hz['Mispointing'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Data_1Hz['N_valid'][r] = np.fromfile(fid,dtype='>i2',count=1)
		#-- CryoSat-2 External Corrections Group for record r
		Corrections['dryTrop'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['wetTrop'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['InvBar'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['DAC'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['Iono'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['SSB'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['ocTideElv'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['lpeTideElv'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['olTideElv'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['seTideElv'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['gpTideElv'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['Spare1'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['Surf_type'][r] = np.fromfile(fid,dtype='>u8',count=1)
		Corrections['MSS_Geoid'][r] = np.fromfile(fid,dtype='>i4',count=1)
		Corrections['ODLE'][r] = np.fromfile(fid,dtype='>i4',count=1)
		Corrections['Ice_conc'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['Snow_depth'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['Snow_density'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['Spare2'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['C_status'][r] = np.fromfile(fid,dtype='>u4',count=1)
		Corrections['SWH'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['Wind_speed'][r] = np.fromfile(fid,dtype='>u2',count=1)
		Corrections['Spare3'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['Spare4'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['Spare5'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['Spare6'][r] = np.fromfile(fid,dtype='>i2',count=1)
		#-- CryoSat-2 Measurements Group for record r and block b
		for b in range(n_blocks):
			Data_20Hz['D_time_mics'].data[r,b] = np.fromfile(fid,dtype='>i4',count=1)
			Data_20Hz['Lat'].data[r,b] = np.fromfile(fid,dtype='>i4',count=1)
			Data_20Hz['Lon'].data[r,b] = np.fromfile(fid,dtype='>i4',count=1)
			Data_20Hz['Elev'].data[r,b] = np.fromfile(fid,dtype='>i4',count=1)
			Data_20Hz['SSHA_interp'].data[r,b] = np.fromfile(fid,dtype='>i2',count=1)
			Data_20Hz['SSHA_num'].data[r,b] = np.fromfile(fid,dtype='>i2',count=1)
			Data_20Hz['SSHA_qual'].data[r,b] = np.fromfile(fid,dtype='>i2',count=1)
			Data_20Hz['Sig0'].data[r,b] = np.fromfile(fid,dtype='>i2',count=1)
			Data_20Hz['Peakiness'].data[r,b] = np.fromfile(fid,dtype='>u2',count=1)
			Data_20Hz['Freeboard'].data[r,b] = np.fromfile(fid,dtype='>i2',count=1)
			Data_20Hz['N_avg'].data[r,b] = np.fromfile(fid,dtype='>i2',count=1)
			Data_20Hz['Spare1'].data[r,b] = np.fromfile(fid,dtype='>i2',count=1)
			Data_20Hz['Quality_Flg'].data[r,b] = np.fromfile(fid,dtype='>u4',count=1)
			Data_20Hz['Spare2'].data[r,b] = np.fromfile(fid,dtype='>i2',count=1)
			Data_20Hz['Spare3'].data[r,b] = np.fromfile(fid,dtype='>i2',count=1)
			Data_20Hz['Spare4'].data[r,b] = np.fromfile(fid,dtype='>i2',count=1)
			Data_20Hz['Spare5'].data[r,b] = np.fromfile(fid,dtype='>i2',count=1)
		#-- Set CryoSat-2 Measurements Group Masks for record r
		Data_20Hz['D_time_mics'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['Lat'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['Lon'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['Elev'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['SSHA_interp'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['SSHA_num'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['SSHA_qual'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['Sig0'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['Peakiness'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['Freeboard'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['N_avg'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['Spare1'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['Quality_Flg'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['Spare2'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['Spare3'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['Spare4'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['Spare5'].mask[r,:Data_1Hz['N_valid'][r]] = False

	#-- Bind all the bits of the l2_mds together into a single dictionary
	CS_l2_mds = {}
	CS_l2_mds['Data_1Hz'] = Data_1Hz
	CS_l2_mds['Corrections'] = Corrections
	CS_l2_mds['Data_20Hz'] = Data_20Hz
	#-- return the output dictionary
	return CS_l2_mds

#-- PURPOSE: Initiate L2 MDS variables for CryoSat Baseline C
def cryosat_baseline_C(fid,record_size,n_records):
	#-- CryoSat-2 1 Hz data fields (Location Group)
	#-- Time and Orbit Parameters plus Measurement Mode
	Data_1Hz = {}
	#-- Time: day part
	Data_1Hz['Day'] = np.zeros((n_records),dtype=np.int32)
	#-- Time: second part
	Data_1Hz['Second'] = np.zeros((n_records),dtype=np.int32)
	#-- Time: microsecond part
	Data_1Hz['Micsec'] = np.zeros((n_records),dtype=np.int32)
	#-- SIRAL mode
	Data_1Hz['Siral_mode'] = np.zeros((n_records),dtype=np.uint64)
	#-- Lat_1Hz: packed units (0.1 micro-degree, 1e-7 degrees)
	Data_1Hz['Lat_1Hz'] = np.zeros((n_records),dtype=np.int32)
	#-- Lon_1Hz: packed units (0.1 micro-degree, 1e-7 degrees)
	Data_1Hz['Lon_1Hz'] = np.zeros((n_records),dtype=np.int32)
	#-- Alt_1Hz: packed units (mm, 1e-3 m)
	#-- Altitude of COG above reference ellipsoid (interpolated value)
	Data_1Hz['Alt_1Hz'] = np.zeros((n_records),dtype=np.int32)
	#-- Roll: packed units (0.1 micro-degree, 1e-7 degrees)
	Data_1Hz['Roll'] = np.zeros((n_records),dtype=np.int32)
	#-- Pitch: packed units (0.1 micro-degree, 1e-7 degrees)
	Data_1Hz['Pitch'] = np.zeros((n_records),dtype=np.int32)
	#-- Yaw: packed units (0.1 micro-degree, 1e-7 degrees)
	Data_1Hz['Yaw'] = np.zeros((n_records),dtype=np.int32)
	Data_1Hz['Spare'] = np.zeros((n_records),dtype=np.int16)
	#-- Number of valid records in the block of twenty that contain data
	#-- Last few records of the last block of a dataset may be blank blocks
	#-- inserted to bring the file up to a multiple of twenty.
	Data_1Hz['N_valid'] = np.zeros((n_records),dtype=np.int16)

	#-- CryoSat-2 geophysical corrections (External Corrections Group)
	Corrections = {}
	#-- Dry Tropospheric Correction packed units (mm, 1e-3 m)
	Corrections['dryTrop'] = np.zeros((n_records),dtype=np.int16)
	#-- Wet Tropospheric Correction packed units (mm, 1e-3 m)
	Corrections['wetTrop'] = np.zeros((n_records),dtype=np.int16)
	#-- Inverse Barometric Correction packed units (mm, 1e-3 m)
	Corrections['InvBar'] = np.zeros((n_records),dtype=np.int16)
	#-- Dynamic Atmosphere Correction packed units (mm, 1e-3 m)
	Corrections['DAC'] = np.zeros((n_records),dtype=np.int16)
	#-- Ionospheric Correction packed units (mm, 1e-3 m)
	Corrections['Iono'] = np.zeros((n_records),dtype=np.int16)
	#-- Sea State Bias Correction packed units (mm, 1e-3 m)
	Corrections['SSB'] = np.zeros((n_records),dtype=np.int16)
	#-- Ocean tide Correction packed units (mm, 1e-3 m)
	Corrections['ocTideElv'] = np.zeros((n_records),dtype=np.int16)
	#-- Long period equilibrium ocean tide Correction packed units (mm, 1e-3 m)
	Corrections['lpeTideElv'] = np.zeros((n_records),dtype=np.int16)
	#-- Ocean loading tide Correction packed units (mm, 1e-3 m)
	Corrections['olTideElv'] = np.zeros((n_records),dtype=np.int16)
	#-- Solid Earth tide Correction packed units (mm, 1e-3 m)
	Corrections['seTideElv'] = np.zeros((n_records),dtype=np.int16)
	#-- Geocentric Polar tide Correction packed units (mm, 1e-3 m)
	Corrections['gpTideElv'] = np.zeros((n_records),dtype=np.int16)
	Corrections['Spare1'] = np.zeros((n_records),dtype=np.int16)
	#-- Surface Type: Packed in groups of three bits for each of the 20 records
	Corrections['Surf_type'] = np.zeros((n_records),dtype=np.uint64)
	#-- Mean Sea Surface or Geoid packed units (mm, 1e-3 m)
	Corrections['MSS_Geoid'] = np.zeros((n_records),dtype=np.int32)
	#-- Ocean Depth/Land Elevation Model (ODLE) packed units (mm, 1e-3 m)
	Corrections['ODLE'] = np.zeros((n_records),dtype=np.int32)
	#-- Ice Concentration packed units (%/100)
	Corrections['Ice_conc'] = np.zeros((n_records),dtype=np.int16)
	#-- Snow Depth packed units (mm, 1e-3 m)
	Corrections['Snow_depth'] = np.zeros((n_records),dtype=np.int16)
	#-- Snow Density packed units (kg/m^3)
	Corrections['Snow_density'] = np.zeros((n_records),dtype=np.int16)
	Corrections['Spare2'] = np.zeros((n_records),dtype=np.int16)
	#-- Corrections Status Flag
	Corrections['C_status'] = np.zeros((n_records),dtype=np.uint32)
	#-- Significant Wave Height (SWH) packed units (mm, 1e-3)
	Corrections['SWH'] = np.zeros((n_records),dtype=np.int16)
	#-- Wind Speed packed units (mm/s, 1e-3 m/s)
	Corrections['Wind_speed'] = np.zeros((n_records),dtype=np.uint16)
	Corrections['Spare3'] = np.zeros((n_records),dtype=np.int16)
	Corrections['Spare4'] = np.zeros((n_records),dtype=np.int16)
	Corrections['Spare5'] = np.zeros((n_records),dtype=np.int16)
	Corrections['Spare6'] = np.zeros((n_records),dtype=np.int16)

	#-- CryoSat-2 20 Hz data fields (Measurement Group)
	#-- Derived from instrument measurement parameters
	n_blocks = 20
	Data_20Hz = {}
	#-- Delta between the timestamps for 20Hz record and the 1Hz record
	#-- D_time_mics packed units (microseconds)
	Data_20Hz['D_time_mics'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
	Data_20Hz['D_time_mics'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- Lat: packed units (0.1 micro-degree, 1e-7 degrees)
	Data_20Hz['Lat'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
	Data_20Hz['Lat'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- Lon: packed units (0.1 micro-degree, 1e-7 degrees)
	Data_20Hz['Lon'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
	Data_20Hz['Lon'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- Measured elevation above ellipsoid from retracker 1: packed units (mm, 1e-3 m)
	Data_20Hz['Elev_1'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
	Data_20Hz['Elev_1'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- Measured elevation above ellipsoid from retracker 2: packed units (mm, 1e-3 m)
	Data_20Hz['Elev_2'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
	Data_20Hz['Elev_2'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- Measured elevation above ellipsoid from retracker 3: packed units (mm, 1e-3 m)
	Data_20Hz['Elev_3'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
	Data_20Hz['Elev_3'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- Sigma Zero Backscatter for retracker 1: packed units (1e-2 dB)
	Data_20Hz['Sig0_1'] = np.ma.zeros((n_records,n_blocks),dtype=np.int16)
	Data_20Hz['Sig0_1'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- Sigma Zero Backscatter for retracker 2: packed units (1e-2 dB)
	Data_20Hz['Sig0_2'] = np.ma.zeros((n_records,n_blocks),dtype=np.int16)
	Data_20Hz['Sig0_2'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- Sigma Zero Backscatter for retracker 3: packed units (1e-2 dB)
	Data_20Hz['Sig0_3'] = np.ma.zeros((n_records,n_blocks),dtype=np.int16)
	Data_20Hz['Sig0_3'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- Freeboard: packed units (mm, 1e-3 m)
	#-- -9999 default value indicates computation has not been performed
	Data_20Hz['Freeboard'] = np.ma.zeros((n_records,n_blocks),dtype=np.int16)
	Data_20Hz['Freeboard'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- Interpolated Sea Surface Height Anomaly: packed units (mm, 1e-3 m)
	Data_20Hz['SSHA_interp'] = np.ma.zeros((n_records,n_blocks),dtype=np.int16)
	Data_20Hz['SSHA_interp'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- Interpolated Sea Surface Height measurement count
	Data_20Hz['SSHA_num'] = np.ma.zeros((n_records,n_blocks),dtype=np.int16)
	Data_20Hz['SSHA_num'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- Interpolation quality estimate RSS: packed units (mm, 1e-3 m)
	Data_20Hz['SSHA_qual'] = np.ma.zeros((n_records,n_blocks),dtype=np.int16)
	Data_20Hz['SSHA_qual'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- Peakiness: packed units (1e-2)
	Data_20Hz['Peakiness'] = np.ma.zeros((n_records,n_blocks),dtype=np.uint16)
	Data_20Hz['Peakiness'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- Number of averaged echoes or beams
	Data_20Hz['N_avg'] = np.ma.zeros((n_records,n_blocks),dtype=np.int16)
	Data_20Hz['N_avg'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	Data_20Hz['Spare1'] = np.ma.zeros((n_records,n_blocks),dtype=np.int16)
	Data_20Hz['Spare1'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- Quality flags
	Data_20Hz['Quality_Flg'] = np.ma.zeros((n_records,n_blocks),dtype=np.uint32)
	Data_20Hz['Quality_Flg'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- Corrections Application Flag
	Data_20Hz['Corrections_Flg'] = np.ma.zeros((n_records,n_blocks),dtype=np.uint32)
	Data_20Hz['Corrections_Flg'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- Quality metric for retracker 1
	Data_20Hz['Quality_1'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
	Data_20Hz['Quality_1'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- Quality metric for retracker 2
	Data_20Hz['Quality_2'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
	Data_20Hz['Quality_2'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- Quality metric for retracker 3
	Data_20Hz['Quality_3'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
	Data_20Hz['Quality_3'].mask = np.ones((n_records,n_blocks),dtype=np.bool)
	#-- for each record in the CryoSat file
	for r in range(n_records):
		#-- CryoSat-2 Location Group for record r
		Data_1Hz['Day'][r] = np.fromfile(fid,dtype='>i4',count=1)
		Data_1Hz['Second'][r] = np.fromfile(fid,dtype='>i4',count=1)
		Data_1Hz['Micsec'][r] = np.fromfile(fid,dtype='>i4',count=1)
		Data_1Hz['Siral_mode'][r] = np.fromfile(fid,dtype='>u8',count=1)
		Data_1Hz['Lat_1Hz'][r] = np.fromfile(fid,dtype='>i4',count=1)
		Data_1Hz['Lon_1Hz'][r] = np.fromfile(fid,dtype='>i4',count=1)
		Data_1Hz['Alt_1Hz'][r] = np.fromfile(fid,dtype='>i4',count=1)
		Data_1Hz['Roll'][r] = np.fromfile(fid,dtype='>i4',count=1)
		Data_1Hz['Pitch'][r] = np.fromfile(fid,dtype='>i4',count=1)
		Data_1Hz['Yaw'][r] = np.fromfile(fid,dtype='>i4',count=1)
		Data_1Hz['Spare'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Data_1Hz['N_valid'][r] = np.fromfile(fid,dtype='>i2',count=1)
		#-- CryoSat-2 External Corrections Group for record r
		Corrections['dryTrop'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['wetTrop'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['InvBar'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['DAC'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['Iono'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['SSB'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['ocTideElv'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['lpeTideElv'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['olTideElv'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['seTideElv'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['gpTideElv'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['Spare1'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['Surf_type'][r] = np.fromfile(fid,dtype='>u8',count=1)
		Corrections['MSS_Geoid'][r] = np.fromfile(fid,dtype='>i4',count=1)
		Corrections['ODLE'][r] = np.fromfile(fid,dtype='>i4',count=1)
		Corrections['Ice_conc'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['Snow_depth'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['Snow_density'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['Spare2'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['C_status'][r] = np.fromfile(fid,dtype='>u4',count=1)
		Corrections['SWH'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['Wind_speed'][r] = np.fromfile(fid,dtype='>u2',count=1)
		Corrections['Spare3'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['Spare4'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['Spare5'][r] = np.fromfile(fid,dtype='>i2',count=1)
		Corrections['Spare6'][r] = np.fromfile(fid,dtype='>i2',count=1)
		#-- CryoSat-2 Measurements Group for record r and block b
		for b in range(n_blocks):
			Data_20Hz['D_time_mics'].data[r,b] = np.fromfile(fid,dtype='>i4',count=1)
			Data_20Hz['Lat'].data[r,b] = np.fromfile(fid,dtype='>i4',count=1)
			Data_20Hz['Lon'].data[r,b] = np.fromfile(fid,dtype='>i4',count=1)
			Data_20Hz['Elev_1'].data[r,b] = np.fromfile(fid,dtype='>i4',count=1)
			Data_20Hz['Elev_2'].data[r,b] = np.fromfile(fid,dtype='>i4',count=1)
			Data_20Hz['Elev_3'].data[r,b] = np.fromfile(fid,dtype='>i4',count=1)
			Data_20Hz['Sig0_1'].data[r,b] = np.fromfile(fid,dtype='>i2',count=1)
			Data_20Hz['Sig0_2'].data[r,b] = np.fromfile(fid,dtype='>i2',count=1)
			Data_20Hz['Sig0_3'].data[r,b] = np.fromfile(fid,dtype='>i2',count=1)
			Data_20Hz['Freeboard'].data[r,b] = np.fromfile(fid,dtype='>i2',count=1)
			Data_20Hz['SSHA_interp'].data[r,b] = np.fromfile(fid,dtype='>i2',count=1)
			Data_20Hz['SSHA_num'].data[r,b] = np.fromfile(fid,dtype='>i2',count=1)
			Data_20Hz['SSHA_qual'].data[r,b] = np.fromfile(fid,dtype='>i2',count=1)
			Data_20Hz['Peakiness'].data[r,b] = np.fromfile(fid,dtype='>u2',count=1)
			Data_20Hz['N_avg'].data[r,b] = np.fromfile(fid,dtype='>i2',count=1)
			Data_20Hz['Spare1'].data[r,b] = np.fromfile(fid,dtype='>i2',count=1)
			Data_20Hz['Quality_Flg'].data[r,b] = np.fromfile(fid,dtype='>u4',count=1)
			Data_20Hz['Corrections_Flg'].data[r,b] = np.fromfile(fid,dtype='>u4',count=1)
			Data_20Hz['Quality_1'].data[r,b] = np.fromfile(fid,dtype='>i4',count=1)
			Data_20Hz['Quality_2'].data[r,b] = np.fromfile(fid,dtype='>i4',count=1)
			Data_20Hz['Quality_3'].data[r,b] = np.fromfile(fid,dtype='>i4',count=1)
		#-- Set CryoSat-2 Measurements Group Masks for record r
		Data_20Hz['D_time_mics'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['Lat'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['Lon'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['Elev_1'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['Elev_2'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['Elev_3'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['Sig0_1'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['Sig0_2'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['Sig0_3'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['Freeboard'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['SSHA_interp'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['SSHA_num'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['SSHA_qual'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['Peakiness'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['N_avg'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['Spare1'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['Quality_Flg'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['Corrections_Flg'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['Quality_1'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['Quality_2'].mask[r,:Data_1Hz['N_valid'][r]] = False
		Data_20Hz['Quality_3'].mask[r,:Data_1Hz['N_valid'][r]] = False

	#-- Bind all the bits of the l2_mds together into a single dictionary
	CS_l2_mds = {}
	CS_l2_mds['Data_1Hz'] = Data_1Hz
	CS_l2_mds['Corrections'] = Corrections
	CS_l2_mds['Data_20Hz'] = Data_20Hz
	#-- return the output dictionary
	return CS_l2_mds

#-- PURPOSE: Initiate L2 MDS variables for CryoSat Baseline D (netCDF4)
def cryosat_baseline_D(full_filename, UNPACK=False):
	#-- open netCDF4 file for reading
	fid = netCDF4.Dataset(os.path.expanduser(full_filename),'r')
	#-- use original unscaled units unless UNPACK=True
	fid.set_auto_scale(UNPACK)
	#-- get dimensions
	n_records, = fid.variables['time_cor_01'].shape
	time_cor_01 = fid.variables['time_cor_01'][:].copy()
	#-- CryoSat-2 1 Hz data fields (Location Group)
	#-- Time and Orbit Parameters plus Measurement Mode
	Data_1Hz = {}
	#-- Time (seconds since 2000-01-01)
	Data_1Hz['Time'] = time_cor_01.copy()
	#-- Time: day part
	Data_1Hz['Day'] = np.array(time_cor_01/86400.0,dtype=np.int32)
	#-- Time: second part
	Data_1Hz['Second'] = np.array(time_cor_01-Data_1Hz['Day'][:]*86400.0,dtype=np.int32)
	#-- Time: microsecond part
	Data_1Hz['Micsec'] = np.array((time_cor_01-Data_1Hz['Day'][:]*86400.0-
		Data_1Hz['Second'][:])*1e6,dtype=np.int32)
	#-- Lat_1Hz: packed units (0.1 micro-degree, 1e-7 degrees)
	Data_1Hz['Lat_1Hz'] = fid.variables['lat_01'][:].copy()
	#-- Lon_1Hz: packed units (0.1 micro-degree, 1e-7 degrees)
	Data_1Hz['Lon_1Hz'] = fid.variables['lon_01'][:].copy()
	#-- Alt_1Hz: packed units (mm, 1e-3 m)
	#-- Altitude of COG above reference ellipsoid (interpolated value)
	Data_1Hz['Alt_1Hz'] = fid.variables['alt_01'][:].copy()
	#-- Roll: packed units (0.1 micro-degree, 1e-7 degrees)
	Data_1Hz['Roll'] = fid.variables['off_nadir_roll_angle_str_01'][:].copy()
	#-- Pitch: packed units (0.1 micro-degree, 1e-7 degrees)
	Data_1Hz['Pitch'] = fid.variables['off_nadir_pitch_angle_str_01'][:].copy()
	#-- Yaw: packed units (0.1 micro-degree, 1e-7 degrees)
	Data_1Hz['Yaw'] = fid.variables['off_nadir_yaw_angle_str_01'][:].copy()
	#-- Number of valid records in the block of twenty that contain data
	#-- Last few records of the last block of a dataset may be blank blocks
	#-- inserted to bring the file up to a multiple of twenty.
	Data_1Hz['N_valid'] = fid.variables['num_valid_01'][:].copy()
	#-- add absolute orbit number to 1Hz data
	Data_1Hz['Abs_Orbit'] = np.zeros((n_records),dtype=np.uint32)
	Data_1Hz['Abs_Orbit'][:] = np.uint32(fid.abs_orbit_number)
	#-- add ascending/descending flag to 1Hz data (A=ascending,D=descending)
	Data_1Hz['Ascending_Flg'] = np.zeros((n_records),dtype=np.bool)
	Data_1Hz['Ascending_Flg'][:] = (fid.ascending_flag == 'A')

	#-- CryoSat-2 geophysical corrections (External Corrections Group)
	Corrections = {}
	#-- Dry Tropospheric Correction packed units (mm, 1e-3 m)
	Corrections['dryTrop'] = fid.variables['mod_dry_tropo_cor_01'][:].copy()
	#-- Wet Tropospheric Correction packed units (mm, 1e-3 m)
	Corrections['wetTrop'] = fid.variables['mod_wet_tropo_cor_01'][:].copy()
	#-- Inverse Barometric Correction packed units (mm, 1e-3 m)
	Corrections['InvBar'] = fid.variables['inv_bar_cor_01'][:].copy()
	#-- Dynamic Atmosphere Correction packed units (mm, 1e-3 m)
	Corrections['DAC'] = fid.variables['hf_fluct_total_cor_01'][:].copy()
	#-- Ionospheric Correction packed units (mm, 1e-3 m)
	Corrections['Iono'] = fid.variables['iono_cor_01'][:].copy()
	Corrections['Iono_GIM'] = fid.variables['iono_cor_gim_01'][:].copy()
	#-- Sea State Bias Correction packed units (mm, 1e-3 m)
	Corrections['SSB'] = fid.variables['sea_state_bias_01_ku'][:].copy()
	#-- Ocean tide Correction packed units (mm, 1e-3 m)
	Corrections['ocTideElv'] = fid.variables['ocean_tide_01'][:].copy()
	#-- Long period equilibrium ocean tide Correction packed units (mm, 1e-3 m)
	Corrections['lpeTideElv'] = fid.variables['ocean_tide_eq_01'][:].copy()
	#-- Ocean loading tide Correction packed units (mm, 1e-3 m)
	Corrections['olTideElv'] = fid.variables['load_tide_01'][:].copy()
	#-- Solid Earth tide Correction packed units (mm, 1e-3 m)
	Corrections['seTideElv'] = fid.variables['solid_earth_tide_01'][:].copy()
	#-- Geocentric Polar tide Correction packed units (mm, 1e-3 m)
	Corrections['gpTideElv'] = fid.variables['pole_tide_01'][:].copy()
	#-- Mean Sea Surface and Geoid packed units (mm, 1e-3 m)
	Corrections['Geoid'] = fid.variables['geoid_01'][:].copy()
	Corrections['MSS'] = fid.variables['mean_sea_surf_sea_ice_01'][:].copy()
	#-- Ocean Depth/Land Elevation Model (ODLE) packed units (mm, 1e-3 m)
	Corrections['ODLE'] = fid.variables['odle_01'][:].copy()
	#-- Ice Concentration packed units (%/100)
	Corrections['Ice_conc'] = fid.variables['sea_ice_concentration_01'][:].copy()
	#-- Snow Depth packed units (mm, 1e-3 m)
	Corrections['Snow_depth'] = fid.variables['snow_depth_01'][:].copy()
	#-- Snow Density packed units (kg/m^3)
	Corrections['Snow_density'] = fid.variables['snow_density_01'][:].copy()
	#-- Corrections Status Flag
	Corrections['C_status'] = fid.variables['flag_cor_err_01'][:].copy()
	#-- Significant Wave Height (SWH) packed units (mm, 1e-3)
	Corrections['SWH'] = fid.variables['swh_ocean_01_ku'][:].copy()
	#-- Wind Speed packed units (mm/s, 1e-3 m/s)
	Corrections['Wind_speed'] = fid.variables['wind_speed_alt_01_ku'][:].copy()

	#-- CryoSat-2 20 Hz data fields (Measurement Group)
	#-- Derived from instrument measurement parameters
	n_blocks = 20
	Data_20Hz = {}
	#-- Time (seconds since 2000-01-01)
	Data_20Hz['Time'] = np.ma.zeros((n_records,n_blocks))
	Data_20Hz['Time'].mask = np.ma.ones((n_records,n_blocks),dtype=np.bool)
	time_20_ku = fid.variables['time_20_ku'][:].copy()
	#-- Delta between the timestamps for 20Hz record and the 1Hz record
	#-- D_time_mics packed units (microseconds)
	Data_20Hz['D_time_mics'] = np.ma.zeros((n_records,n_blocks))
	Data_20Hz['D_time_mics'].mask = np.ma.ones((n_records,n_blocks),dtype=np.bool)
	#-- Lat: packed units
	Data_20Hz['Lat'] = np.ma.zeros((n_records,n_blocks))
	Data_20Hz['Lat'].mask = np.ma.ones((n_records,n_blocks),dtype=np.bool)
	lat_poca_20_ku = fid.variables['lat_poca_20_ku'][:].copy()
	#-- Lon: packed units
	Data_20Hz['Lon'] = np.ma.zeros((n_records,n_blocks))
	Data_20Hz['Lon'].mask = np.ma.ones((n_records,n_blocks),dtype=np.bool)
	lon_poca_20_ku = fid.variables['lon_poca_20_ku'][:].copy()
	#-- Measured elevation above ellipsoid from retracker 1
	Data_20Hz['Elev_1'] = np.ma.zeros((n_records,n_blocks))
	Data_20Hz['Elev_1'].mask = np.ma.ones((n_records,n_blocks),dtype=np.bool)
	height_1_20_ku = fid.variables['height_1_20_ku'][:].copy()
	#-- Measured elevation above ellipsoid from retracker 2
	Data_20Hz['Elev_2'] = np.ma.zeros((n_records,n_blocks))
	Data_20Hz['Elev_2'].mask = np.ma.ones((n_records,n_blocks),dtype=np.bool)
	height_2_20_ku = fid.variables['height_2_20_ku'][:].copy()
	#-- Measured elevation above ellipsoid from retracker 3
	Data_20Hz['Elev_3'] = np.ma.zeros((n_records,n_blocks))
	Data_20Hz['Elev_3'].mask = np.ma.ones((n_records,n_blocks),dtype=np.bool)
	height_3_20_ku = fid.variables['height_3_20_ku'][:].copy()
	#-- Sigma Zero Backscatter for retracker 1
	Data_20Hz['Sig0_1'] = np.ma.zeros((n_records,n_blocks))
	Data_20Hz['Sig0_1'].mask = np.ma.ones((n_records,n_blocks),dtype=np.bool)
	sig0_1_20_ku = fid.variables['sig0_1_20_ku'][:].copy()
	#-- Sigma Zero Backscatter for retracker 2
	Data_20Hz['Sig0_2'] = np.ma.zeros((n_records,n_blocks))
	Data_20Hz['Sig0_2'].mask = np.ma.ones((n_records,n_blocks),dtype=np.bool)
	sig0_2_20_ku = fid.variables['sig0_2_20_ku'][:].copy()
	#-- Sigma Zero Backscatter for retracker 3
	Data_20Hz['Sig0_3'] = np.ma.zeros((n_records,n_blocks))
	Data_20Hz['Sig0_3'].mask = np.ma.ones((n_records,n_blocks),dtype=np.bool)
	sig0_3_20_ku = fid.variables['sig0_3_20_ku'][:].copy()
	#-- Measured range from the satellite CoM to the surface from retracker 1
	Data_20Hz['Range_1'] = np.ma.zeros((n_records,n_blocks))
	Data_20Hz['Range_1'].mask = np.ma.ones((n_records,n_blocks),dtype=np.bool)
	range_1_20_ku = fid.variables['range_1_20_ku'][:].copy()
	#-- Measured range from the satellite CoM to the surface from retracker 2
	Data_20Hz['Range_2'] = np.ma.zeros((n_records,n_blocks))
	Data_20Hz['Range_2'].mask = np.ma.ones((n_records,n_blocks),dtype=np.bool)
	range_2_20_ku = fid.variables['range_2_20_ku'][:].copy()
	#-- Measured range from the satellite CoM to the surface from retracker 3
	Data_20Hz['Range_3'] = np.ma.zeros((n_records,n_blocks))
	Data_20Hz['Range_3'].mask = np.ma.ones((n_records,n_blocks),dtype=np.bool)
	range_3_20_ku = fid.variables['range_3_20_ku'][:].copy()
	#-- Freeboard
	Data_20Hz['Freeboard'] = np.ma.zeros((n_records,n_blocks))
	Data_20Hz['Freeboard'].mask = np.ma.ones((n_records,n_blocks),dtype=np.bool)
	freeboard_20_ku = fid.variables['freeboard_20_ku'][:].copy()
	#-- Sea ice Floe height
	Data_20Hz['Sea_Ice_Lead'] = np.ma.zeros((n_records,n_blocks))
	Data_20Hz['Sea_Ice_Lead'].mask = np.ma.ones((n_records,n_blocks),dtype=np.bool)
	height_sea_ice_floe_20_ku = fid.variables['height_sea_ice_floe_20_ku'][:].copy()
	#-- Sea ice lead height
	Data_20Hz['Sea_Ice_Floe'] = np.ma.zeros((n_records,n_blocks))
	Data_20Hz['Sea_Ice_Floe'].mask = np.ma.ones((n_records,n_blocks),dtype=np.bool)
	height_sea_ice_lead_20_ku = fid.variables['height_sea_ice_lead_20_ku'][:].copy()
	#-- Interpolated Sea Surface Height Anomaly
	Data_20Hz['SSHA_interp'] = np.ma.zeros((n_records,n_blocks))
	Data_20Hz['SSHA_interp'].mask = np.ma.ones((n_records,n_blocks),dtype=np.bool)
	ssha_interp_20_ku = fid.variables['ssha_interp_20_ku'][:].copy()
	#-- Interpolated Sea Surface Height measurement count
	Data_20Hz['SSHA_num'] = np.ma.zeros((n_records,n_blocks))
	Data_20Hz['SSHA_num'].mask = np.ma.ones((n_records,n_blocks),dtype=np.bool)
	ssha_interp_numval_20_ku = fid.variables['ssha_interp_numval_20_ku'][:].copy()
	#-- Interpolation quality estimate RSS
	Data_20Hz['SSHA_qual'] = np.ma.zeros((n_records,n_blocks))
	Data_20Hz['SSHA_qual'].mask = np.ma.ones((n_records,n_blocks),dtype=np.bool)
	ssha_interp_rms_20_ku = fid.variables['ssha_interp_rms_20_ku'][:].copy()
	#-- Peakiness
	Data_20Hz['Peakiness'] = np.ma.zeros((n_records,n_blocks))
	Data_20Hz['Peakiness'].mask = np.ma.ones((n_records,n_blocks),dtype=np.bool)
	peakiness_20_ku = fid.variables['peakiness_20_ku'][:].copy()
	#-- Number of averaged echoes or beams
	Data_20Hz['N_avg'] = np.ma.zeros((n_records,n_blocks))
	Data_20Hz['N_avg'].mask = np.ma.ones((n_records,n_blocks),dtype=np.bool)
	echo_avg_numval_20_ku = fid.variables['echo_avg_numval_20_ku'][:].copy()
	#-- Quality flags
	Data_20Hz['Quality_Flg'] = np.ma.zeros((n_records,n_blocks))
	Data_20Hz['Quality_Flg'].mask = np.ma.ones((n_records,n_blocks),dtype=np.bool)
	flag_prod_status_20_ku = fid.variables['flag_prod_status_20_ku'][:].copy()
	#-- Corrections Application Flag
	Data_20Hz['Corrections_Flg'] = np.ma.zeros((n_records,n_blocks))
	Data_20Hz['Corrections_Flg'].mask = np.ma.ones((n_records,n_blocks),dtype=np.bool)
	flag_cor_applied_20_ku = fid.variables['flag_cor_applied_20_ku'][:].copy()
	#-- Measurement mode
	Data_20Hz['Measurement_Mode'] = np.ma.zeros((n_records,n_blocks))
	Data_20Hz['Measurement_Mode'].mask = np.ma.ones((n_records,n_blocks),dtype=np.bool)
	flag_instr_mode_op_20_ku = fid.variables['flag_instr_mode_op_20_ku'][:].copy()
	#-- Surface Type
	Data_20Hz['Surf_type'] = np.ma.zeros((n_records,n_blocks))
	Data_20Hz['Surf_type'].mask = np.ma.ones((n_records,n_blocks),dtype=np.bool)
	surf_type_20_ku = fid.variables['surf_type_20_ku'][:].copy()
	#-- Quality metric for retracker 1
	Data_20Hz['Quality_1'] = np.ma.zeros((n_records,n_blocks))
	Data_20Hz['Quality_1'].mask = np.ma.ones((n_records,n_blocks),dtype=np.bool)
	retracker_1_quality_20_ku = fid.variables['retracker_1_quality_20_ku'][:].copy()
	#-- Quality metric for retracker 2
	Data_20Hz['Quality_2'] = np.ma.zeros((n_records,n_blocks))
	Data_20Hz['Quality_2'].mask = np.ma.ones((n_records,n_blocks),dtype=np.bool)
	retracker_2_quality_20_ku = fid.variables['retracker_2_quality_20_ku'][:].copy()
	#-- Quality metric for retracker 3
	Data_20Hz['Quality_3'] = np.ma.zeros((n_records,n_blocks))
	Data_20Hz['Quality_3'].mask = np.ma.ones((n_records,n_blocks),dtype=np.bool)
	retracker_3_quality_20_ku = fid.variables['retracker_3_quality_20_ku'][:].copy()
	#-- for each record in the CryoSat file
	for r in range(n_records):
		#-- index for record r
		idx = fid.variables['ind_first_meas_20hz_01'][r].copy()
		#-- number of valid blocks in record r
		cnt = np.copy(fid.variables['num_valid_01'][r])
		#-- CryoSat-2 Measurements Group for record r
		Data_20Hz['Time'].data[r,:cnt] = time_20_ku[idx:idx+cnt].copy()
		Data_20Hz['Time'].mask[r,:cnt] = False
		Data_20Hz['D_time_mics'].data[r,:cnt] = 1e6*(time_20_ku[idx:idx+cnt] - time_cor_01[r])
		Data_20Hz['D_time_mics'].mask[r,:cnt] = False
		Data_20Hz['Lat'].data[r,:cnt] = lat_poca_20_ku[idx:idx+cnt].copy()
		Data_20Hz['Lat'].mask[r,:cnt] = False
		Data_20Hz['Lon'].data[r,:cnt] = lon_poca_20_ku[idx:idx+cnt].copy()
		Data_20Hz['Lon'].mask[r,:cnt] = False
		Data_20Hz['Elev_1'].data[r,:cnt] = height_1_20_ku[idx:idx+cnt].copy()
		Data_20Hz['Elev_1'].mask[r,:cnt] = False
		Data_20Hz['Elev_2'].data[r,:cnt] = height_2_20_ku[idx:idx+cnt].copy()
		Data_20Hz['Elev_2'].mask[r,:cnt] = False
		Data_20Hz['Elev_3'].data[r,:cnt] = height_3_20_ku[idx:idx+cnt].copy()
		Data_20Hz['Elev_3'].mask[r,:cnt] = False
		Data_20Hz['Sig0_1'].data[r,:cnt] = sig0_1_20_ku[idx:idx+cnt].copy()
		Data_20Hz['Sig0_1'].mask[r,:cnt] = False
		Data_20Hz['Sig0_2'].data[r,:cnt] = sig0_2_20_ku[idx:idx+cnt].copy()
		Data_20Hz['Sig0_2'].mask[r,:cnt] = False
		Data_20Hz['Sig0_3'].data[r,:cnt] = sig0_3_20_ku[idx:idx+cnt].copy()
		Data_20Hz['Sig0_3'].mask[r,:cnt] = False
		Data_20Hz['Range_1'].data[r,:cnt] = range_1_20_ku[idx:idx+cnt].copy()
		Data_20Hz['Range_1'].mask[r,:cnt] = False
		Data_20Hz['Range_2'].data[r,:cnt] = range_2_20_ku[idx:idx+cnt].copy()
		Data_20Hz['Range_2'].mask[r,:cnt] = False
		Data_20Hz['Range_3'].data[r,:cnt] = range_3_20_ku[idx:idx+cnt].copy()
		Data_20Hz['Range_3'].mask[r,:cnt] = False
		Data_20Hz['Freeboard'].data[r,:cnt] = freeboard_20_ku[idx:idx+cnt].copy()
		Data_20Hz['Freeboard'].mask[r,:cnt] = False
		Data_20Hz['Sea_Ice_Floe'].data[r,:cnt] = height_sea_ice_floe_20_ku[idx:idx+cnt].copy()
		Data_20Hz['Sea_Ice_Floe'].mask[r,:cnt] = False
		Data_20Hz['Sea_Ice_Lead'].data[r,:cnt] = height_sea_ice_lead_20_ku[idx:idx+cnt].copy()
		Data_20Hz['Sea_Ice_Lead'].mask[r,:cnt] = False
		Data_20Hz['SSHA_interp'].data[r,:cnt] = ssha_interp_20_ku[idx:idx+cnt].copy()
		Data_20Hz['SSHA_interp'].mask[r,:cnt] = False
		Data_20Hz['SSHA_num'].data[r,:cnt] = ssha_interp_numval_20_ku[idx:idx+cnt].copy()
		Data_20Hz['SSHA_num'].mask[r,:cnt] = False
		Data_20Hz['SSHA_qual'].data[r,:cnt] = ssha_interp_rms_20_ku[idx:idx+cnt].copy()
		Data_20Hz['SSHA_qual'].mask[r,:cnt] = False
		Data_20Hz['Peakiness'].data[r,:cnt] = peakiness_20_ku[idx:idx+cnt].copy()
		Data_20Hz['Peakiness'].mask[r,:cnt] = False
		Data_20Hz['N_avg'].data[r,:cnt] = echo_avg_numval_20_ku[idx:idx+cnt].copy()
		Data_20Hz['N_avg'].mask[r,:cnt] = False
		Data_20Hz['Quality_Flg'].data[r,:cnt] = flag_prod_status_20_ku[idx:idx+cnt].copy()
		Data_20Hz['Quality_Flg'].mask[r,:cnt] = False
		Data_20Hz['Corrections_Flg'].data[r,:cnt] = flag_cor_applied_20_ku[idx:idx+cnt].copy()
		Data_20Hz['Corrections_Flg'].mask[r,:cnt] = False
		Data_20Hz['Measurement_Mode'].data[r,:cnt] = flag_instr_mode_op_20_ku[idx:idx+cnt].copy()
		Data_20Hz['Measurement_Mode'].mask[r,:cnt] = False
		Data_20Hz['Surf_type'].data[r,:cnt] = surf_type_20_ku[idx:idx+cnt].copy()
		Data_20Hz['Surf_type'].mask[r,:cnt] = False
		Data_20Hz['Quality_1'].data[r,:cnt] = retracker_1_quality_20_ku[idx:idx+cnt].copy()
		Data_20Hz['Quality_1'].mask[r,:cnt] = False
		Data_20Hz['Quality_2'].data[r,:cnt] = retracker_2_quality_20_ku[idx:idx+cnt].copy()
		Data_20Hz['Quality_2'].mask[r,:cnt] = False
		Data_20Hz['Quality_3'].data[r,:cnt] = retracker_3_quality_20_ku[idx:idx+cnt].copy()
		Data_20Hz['Quality_3'].mask[r,:cnt] = False

	#-- Bind all the variables of the l2_mds together into a single dictionary
	CS_l2_mds = {}
	CS_l2_mds['Data_1Hz'] = Data_1Hz
	CS_l2_mds['Corrections'] = Corrections
	CS_l2_mds['Data_20Hz'] = Data_20Hz

	#-- extract global attributes and assign as MPH and SPH metadata
	CS_l2_mds['METADATA'] = dict(MPH={},SPH={})
	#-- MPH attributes
	CS_l2_mds['METADATA']['MPH']['PRODUCT'] = fid.product_name
	CS_l2_mds['METADATA']['MPH']['DOI'] = fid.doi
	CS_l2_mds['METADATA']['MPH']['PROC_STAGE'] =  fid.processing_stage
	CS_l2_mds['METADATA']['MPH']['REF_DOC'] =  fid.reference_document
	CS_l2_mds['METADATA']['MPH']['ACQUISITION_STATION'] = fid.acquisition_station
	CS_l2_mds['METADATA']['MPH']['PROC_CENTER'] = fid.processing_centre
	CS_l2_mds['METADATA']['MPH']['PROC_TIME'] = fid.creation_time
	CS_l2_mds['METADATA']['MPH']['SOFTWARE_VER'] = fid.software_version
	CS_l2_mds['METADATA']['MPH']['SENSING_START'] = fid.sensing_start
	CS_l2_mds['METADATA']['MPH']['SENSING_STOP'] = fid.sensing_stop
	CS_l2_mds['METADATA']['MPH']['PHASE'] = fid.phase
	CS_l2_mds['METADATA']['MPH']['CYCLE'] = fid.cycle_number
	CS_l2_mds['METADATA']['MPH']['REL_ORBIT'] = fid.rel_orbit_number
	CS_l2_mds['METADATA']['MPH']['ABS_ORBIT'] = fid.abs_orbit_number
	CS_l2_mds['METADATA']['MPH']['STATE_VECTOR_TIME'] = fid.state_vector_time
	CS_l2_mds['METADATA']['MPH']['DELTA_UT1'] = fid.delta_ut1
	CS_l2_mds['METADATA']['MPH']['X_POSITION'] = fid.x_position
	CS_l2_mds['METADATA']['MPH']['Y_POSITION'] = fid.y_position
	CS_l2_mds['METADATA']['MPH']['Z_POSITION'] = fid.z_position
	CS_l2_mds['METADATA']['MPH']['X_VELOCITY'] = fid.x_velocity
	CS_l2_mds['METADATA']['MPH']['Y_VELOCITY'] = fid.y_velocity
	CS_l2_mds['METADATA']['MPH']['Z_VELOCITY'] = fid.z_velocity
	CS_l2_mds['METADATA']['MPH']['VECTOR_SOURCE'] = fid.vector_source
	CS_l2_mds['METADATA']['MPH']['LEAP_UTC'] = fid.leap_utc
	CS_l2_mds['METADATA']['MPH']['LEAP_SIGN'] = fid.leap_sign
	CS_l2_mds['METADATA']['MPH']['LEAP_ERR'] = fid.leap_err
	CS_l2_mds['METADATA']['MPH']['PRODUCT_ERR'] = fid.product_err
	#-- SPH attributes
	CS_l2_mds['METADATA']['SPH']['START_RECORD_TAI_TIME'] = fid.first_record_time
	CS_l2_mds['METADATA']['SPH']['STOP_RECORD_TAI_TIME'] = fid.last_record_time
	CS_l2_mds['METADATA']['SPH']['ABS_ORBIT_START'] = fid.abs_orbit_start
	CS_l2_mds['METADATA']['SPH']['REL_TIME_ASC_NODE_START'] = fid.rel_time_acs_node_start
	CS_l2_mds['METADATA']['SPH']['ABS_ORBIT_STOP'] = fid.abs_orbit_stop
	CS_l2_mds['METADATA']['SPH']['REL_TIME_ASC_NODE_STOP'] = fid.rel_time_acs_node_stop
	CS_l2_mds['METADATA']['SPH']['EQUATOR_CROSS_TIME_UTC'] = fid.equator_cross_time
	CS_l2_mds['METADATA']['SPH']['EQUATOR_CROSS_LONG'] = fid.equator_cross_long
	CS_l2_mds['METADATA']['SPH']['ASCENDING_FLAG'] = fid.ascending_flag
	CS_l2_mds['METADATA']['SPH']['START_LAT'] = fid.first_record_lat
	CS_l2_mds['METADATA']['SPH']['START_LONG'] = fid.first_record_lon
	CS_l2_mds['METADATA']['SPH']['STOP_LAT'] = fid.last_record_lat
	CS_l2_mds['METADATA']['SPH']['STOP_LONG'] = fid.last_record_lon
	CS_l2_mds['METADATA']['SPH']['L1_PROC_FLAG'] = fid.l1b_proc_flag
	CS_l2_mds['METADATA']['SPH']['L1_PROCESSING_QUALITY'] = fid.l1b_processing_quality
	CS_l2_mds['METADATA']['SPH']['L1_PROC_THRESH'] = fid.l1b_proc_thresh
	CS_l2_mds['METADATA']['SPH']['INSTR_ID'] = fid.instr_id
	CS_l2_mds['METADATA']['SPH']['LRM_MODE_PERCENT'] = fid.lrm_mode_percent
	CS_l2_mds['METADATA']['SPH']['SAR_MODE_PERCENT'] = fid.sar_mode_percent
	CS_l2_mds['METADATA']['SPH']['SARIN_MODE_PERCENT'] = fid.sarin_mode_percent
	CS_l2_mds['METADATA']['SPH']['OPEN_OCEAN_PERCENT'] = fid.open_ocean_percent
	CS_l2_mds['METADATA']['SPH']['CLOSE_SEA_PERCENT'] = fid.close_sea_percent
	CS_l2_mds['METADATA']['SPH']['CONTINENT_ICE_PERCENT'] = fid.continent_ice_percent
	CS_l2_mds['METADATA']['SPH']['LAND_PERCENT'] = fid.land_percent
	CS_l2_mds['METADATA']['SPH']['L2_PROD_STATUS'] = fid.l2_prod_status
	CS_l2_mds['METADATA']['SPH']['L2_PROC_FLAG'] = fid.l2_proc_flag
	CS_l2_mds['METADATA']['SPH']['L2_PROCESSING_QUALITY'] = fid.l2_processing_quality
	CS_l2_mds['METADATA']['SPH']['L2_PROC_THRESH'] = fid.l2_proc_thresh
	CS_l2_mds['METADATA']['SPH']['SIR_CONFIGURATION'] = fid.sir_configuration
	CS_l2_mds['METADATA']['SPH']['SIR_OP_MODE'] = fid.sir_op_mode
	CS_l2_mds['METADATA']['SPH']['ORBIT_FILE'] = fid.xref_orbit
	CS_l2_mds['METADATA']['SPH']['PROC_CONFIG_PARAMS_FILE'] = fid.xref_pconf
	CS_l2_mds['METADATA']['SPH']['CONSTANTS_FILE'] = fid.xref_constants
	CS_l2_mds['METADATA']['SPH']['IPF_RA_DATABASE_FILE'] = fid.xref_siral_characterisation
	CS_l2_mds['METADATA']['SPH']['DORIS_USO_DRIFT_FILE'] = fid.xref_uso
	CS_l2_mds['METADATA']['SPH']['STAR_TRACKER_ATTREF_FILE'] = fid.xref_star_tracker_attref
	CS_l2_mds['METADATA']['SPH']['SIRAL_LEVEL_0_FILE'] = fid.xref_siral_l0
	CS_l2_mds['METADATA']['SPH']['CALIBRATION_TYPE_1_FILE'] = fid.xref_cal1
	CS_l2_mds['METADATA']['SPH']['SIR_COMPLEX_CAL1_SARIN'] = fid.xref_cal1_sarin
	CS_l2_mds['METADATA']['SPH']['SCENARIO_FILE'] = fid.xref_orbit_scenario
	CS_l2_mds['METADATA']['SPH']['CALIBRATION_TYPE_2_FILE'] = fid.xref_cal2
	CS_l2_mds['METADATA']['SPH']['SURFACE_PRESSURE_FILE'] = fid.xref_surf_pressure
	CS_l2_mds['METADATA']['SPH']['MEAN_PRESSURE_FILE'] = fid.xref_mean_pressure
	CS_l2_mds['METADATA']['SPH']['WET_TROPOSPHERE_FILE'] = fid.xref_wet_trop
	CS_l2_mds['METADATA']['SPH']['U_WIND_FILE'] = fid.xref_u_wind
	CS_l2_mds['METADATA']['SPH']['V_WIND_FILE'] = fid.xref_v_wind
	CS_l2_mds['METADATA']['SPH']['METEO_GRID_DEF_FILE'] = fid.xref_meteo
	CS_l2_mds['METADATA']['SPH']['S1S2_PRESSURE_00H_MAP'] = fid.xref_s1s2_pressure_00h
	CS_l2_mds['METADATA']['SPH']['S1S2_PRESSURE_06H_MAP'] = fid.xref_s1s2_pressure_06h
	CS_l2_mds['METADATA']['SPH']['S1S2_PRESSURE_12H_MAP'] = fid.xref_s1s2_pressure_12h
	CS_l2_mds['METADATA']['SPH']['S1S2_PRESSURE_18H_MAP'] = fid.xref_s1s2_pressure_18h
	CS_l2_mds['METADATA']['SPH']['S1_TIDE_AMPLITUDE_MAP'] = fid.xref_s1_tide_amplitude
	CS_l2_mds['METADATA']['SPH']['S1_TIDE_PHASE_MAP'] = fid.xref_s1_tide_phase
	CS_l2_mds['METADATA']['SPH']['S2_TIDE_AMPLITUDE_MAP'] = fid.xref_s2_tide_amplitude
	CS_l2_mds['METADATA']['SPH']['S2_TIDE_PHASE_MAP'] = fid.xref_s2_tide_phase
	CS_l2_mds['METADATA']['SPH']['GPS_IONO_MAP'] = fid.xref_gim
	CS_l2_mds['METADATA']['SPH']['MODIFIED_DIP_MAP_FILE'] = fid.xref_dip_map
	CS_l2_mds['METADATA']['SPH']['IONO_COEFFICENTS_FILE'] = fid.xref_iono_cor
	CS_l2_mds['METADATA']['SPH']['SAI_FILE'] = fid.xref_sai
	CS_l2_mds['METADATA']['SPH']['OCEAN_TIDE_FILE'] = fid.xref_ocean_tide
	CS_l2_mds['METADATA']['SPH']['TIDAL_LOADING_FILE'] = fid.xref_tidal_load
	CS_l2_mds['METADATA']['SPH']['EARTH_TIDE_FILE'] = fid.xref_earth_tide
	CS_l2_mds['METADATA']['SPH']['POLE_TIDE_FILE'] = fid.xref_pole_location
	CS_l2_mds['METADATA']['SPH']['SURFACE_TYPE_FILE'] = fid.xref_surf_type
	CS_l2_mds['METADATA']['SPH']['AUX_MOG2D'] = fid.xref_mog2d
	CS_l2_mds['METADATA']['SPH']['SIRAL_LEVEL_1B_FILE'] = fid.xref_siral_l1b
	CS_l2_mds['METADATA']['SPH']['MEAN_SEA_SURFACE_FILE'] = fid.xref_mss
	CS_l2_mds['METADATA']['SPH']['GEOID_FILE'] = fid.xref_geoid
	CS_l2_mds['METADATA']['SPH']['ODLE_FILE'] = fid.xref_odle
	#-- mode dependent attributes
	if ('xref_dem' in fid.ncattrs()):
		CS_l2_mds['METADATA']['SPH']['DEM_MODEL_FILE'] = fid.xref_dem
	if ('xref_sea_ice' in fid.ncattrs()):
		CS_l2_mds['METADATA']['SPH']['SEA_ICE_FILE'] = fid.xref_sea_ice
	if ('xref_snow_depth' in fid.ncattrs()):
		CS_l2_mds['METADATA']['SPH']['SNOW_DEPTH_FILE'] = fid.xref_snow_depth

	#-- close the netCDF4 file
	fid.close()
	#-- return the output dictionary
	return CS_l2_mds

#-- PURPOSE: Read ASCII Main Product Header (MPH) block from an ESA PDS file
def read_MPH(full_filename):
	#-- read input data file
	with open(full_filename, 'rb') as fid:
		file_contents = fid.read().splitlines()

	#-- Define constant values associated with PDS file formats
	#-- number of text lines in standard MPH
	n_MPH_lines	= 41
	#-- check that first line of header matches PRODUCT
	if not bool(re.match(b'PRODUCT\=\"(.*)(?=\")',file_contents[0])):
		raise IOError('File does not start with a valid PDS MPH')
	#-- read MPH header text
	s_MPH_fields = {}
	for i in range(n_MPH_lines):
		#-- use regular expression operators to read headers
		if bool(re.match(b'(.*?)\=\"(.*)(?=\")',file_contents[i])):
			#-- data fields within quotes
			field,value=re.findall(b'(.*?)\=\"(.*)(?=\")',file_contents[i]).pop()
			s_MPH_fields[field.decode('utf-8')] = value.decode('utf-8').rstrip()
		elif bool(re.match(b'(.*?)\=(.*)',file_contents[i])):
			#-- data fields without quotes
			field,value=re.findall(b'(.*?)\=(.*)',file_contents[i]).pop()
			s_MPH_fields[field.decode('utf-8')] = value.decode('utf-8').rstrip()

	#-- Return block name array to calling function
	return s_MPH_fields

#-- PURPOSE: Read ASCII Specific Product Header (SPH) block from a PDS file
def read_SPH(full_filename,j_sph_size):
	#-- read input data file
	with open(full_filename, 'rb') as fid:
		file_contents = fid.read().splitlines()

	#-- Define constant values associated with PDS file formats
	#-- number of text lines in standard MPH
	n_MPH_lines	= 41
	#-- compile regular expression operator for reading headers
	rx = re.compile(b'(.*?)\=\"?(.*)',re.VERBOSE)
	#-- check first line of header matches SPH_DESCRIPTOR
	if not bool(re.match(b'SPH\_DESCRIPTOR\=',file_contents[n_MPH_lines+1])):
		raise IOError('File does not have a valid PDS DSD')
	#-- read SPH header text (no binary control characters)
	s_SPH_lines = [li for li in file_contents[n_MPH_lines+1:] if rx.match(li)
		and not re.search(b'[^\x20-\x7e]+',li)]

	#-- extract SPH header text
	s_SPH_fields = {}
	c = 0
	while (c < len(s_SPH_lines)):
		#-- check if line is within DS_NAME portion of SPH header
		if bool(re.match(b'DS_NAME',s_SPH_lines[c])):
			#-- add dictionary for DS_NAME
			field,value=re.findall(b'(.*?)\=\"(.*)(?=\")',s_SPH_lines[c]).pop()
			key = value.decode('utf-8').rstrip()
			s_SPH_fields[key] = {}
			for line in s_SPH_lines[c+1:c+7]:
				if bool(re.match(b'(.*?)\=\"(.*)(?=\")',line)):
					#-- data fields within quotes
					dsfield,dsvalue=re.findall(b'(.*?)\=\"(.*)(?=\")',line).pop()
					s_SPH_fields[key][dsfield.decode('utf-8')] = dsvalue.decode('utf-8').rstrip()
				elif bool(re.match(b'(.*?)\=(.*)',line)):
					#-- data fields without quotes
					dsfield,dsvalue=re.findall(b'(.*?)\=(.*)',line).pop()
					s_SPH_fields[key][dsfield.decode('utf-8')] = dsvalue.decode('utf-8').rstrip()
			#-- add 6 to counter to go to next entry
			c += 6
		#-- use regular expression operators to read headers
		elif bool(re.match(b'(.*?)\=\"(.*)(?=\")',s_SPH_lines[c])):
			#-- data fields within quotes
			field,value=re.findall(b'(.*?)\=\"(.*)(?=\")',s_SPH_lines[c]).pop()
			s_SPH_fields[field.decode('utf-8')] = value.decode('utf-8').rstrip()
		elif bool(re.match(b'(.*?)\=(.*)',s_SPH_lines[c])):
			#-- data fields without quotes
			field,value=re.findall(b'(.*?)\=(.*)',s_SPH_lines[c]).pop()
			s_SPH_fields[field.decode('utf-8')] = value.decode('utf-8').rstrip()
		#-- add 1 to counter to go to next line
		c += 1

	#-- Return block name array to calling function
	return s_SPH_fields

#-- PURPOSE: Read ASCII Data Set Descriptors (DSD) block from a PDS file
def read_DSD(full_filename):
	#-- read input data file
	with open(full_filename, 'rb') as fid:
		file_contents = fid.read().splitlines()

	#-- Define constant values associated with PDS file formats
	#-- number of text lines in standard MPH
	n_MPH_lines	= 41
	#-- number of text lines in a DSD header
	n_DSD_lines = 8

	#-- Level-2 CryoSat DS_NAMES within files
	regex_patterns = []
	regex_patterns.append(b'DS_NAME\="SIR_LRM_L2(_I)?[\s+]*"')
	regex_patterns.append(b'DS_NAME\="SIR_LRMIL2[\s+]*"')
	regex_patterns.append(b'DS_NAME\="SIR_SAR_L2(A|B)?(_I)?[\s+]*"')
	regex_patterns.append(b'DS_NAME\="SIR_SARIL2(A|B)?[\s+]*"')
	regex_patterns.append(b'DS_NAME\="SIR_FDM_L2[\s+]*"')
	regex_patterns.append(b'DS_NAME\="SIR_SIN_L2(_I)?[\s+]*"')
	regex_patterns.append(b'DS_NAME\="SIR_SINIL2[\s+]*"')
	regex_patterns.append(b'DS_NAME\="SIR_SID_L2(_I)?[\s+]*"')
	regex_patterns.append(b'DS_NAME\="SIR_SIDIL2[\s+]*"')
	regex_patterns.append(b'DS_NAME\="SIR_GDR_2(A|B|_)?[\s+]*"')
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
		if bool(re.match(b'(.*?)\=\"(.*)(?=\")',file_contents[i])):
			#-- data fields within quotes
			field,value=re.findall(b'(.*?)\=\"(.*)(?=\")',file_contents[i]).pop()
			s_DSD_fields[field.decode('utf-8')] = value.decode('utf-8').rstrip()
		elif bool(re.match(b'(.*?)\=(.*)',file_contents[i])):
			#-- data fields without quotes
			field,value=re.findall(b'(.*?)\=(.*)',file_contents[i]).pop()
			s_DSD_fields[field.decode('utf-8')] = value.decode('utf-8').rstrip()

	#-- Return block name array to calling function
	return s_DSD_fields

#-- PURPOSE: read CryoSat Level-2 data
def read_cryosat_L2(full_filename, VERBOSE=False):
	#-- file basename and file extension of input file
	fileBasename,fileExtension=os.path.splitext(os.path.basename(full_filename))

	#-- CryoSat file class
	#-- OFFL (Off Line Processing/Systematic)
	#-- NRT_ (Near Real Time)
	#-- RPRO (ReProcessing)
	#-- TEST (Testing)
	#-- LTA_ (Long Term Archive)
	regex_class = 'OFFL|NRT_|RPRO|TEST|LTA_'
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
	regex_products = ('SIR_LRM_2|SIR_FDM_2|SIR_SIN_2|SIR_SID_2|'
	'SIR_SAR_2|SIR_GDR_2|SIR_LRMI2|SIR_SINI2|SIR_SIDI2|SIR_SARI2')
	#-- CRYOSAT LEVEL-2 PRODUCTS NAMING RULES
	#-- Mission Identifier
	#-- File Class
	#-- File Product
	#-- Validity Start Date and Time
	#-- Validity Stop Date and Time
	#-- Baseline Identifier
	#-- Version Number
	regex_pattern = '(.*?)_({0})_({1})__(\d+T?\d+)_(\d+T?\d+)_(.*?)(\d+)'
	rx = re.compile(regex_pattern.format(regex_class,regex_products),re.VERBOSE)
	#-- extract file information from filename
	MI,CLASS,PRODUCT,START,STOP,BASELINE,VERSION=rx.findall(fileBasename).pop()

	#-- check if input file is original binary *.DBL or new netCDF4 *.nc format
	if (fileExtension == '.nc'):
		print(fileBasename) if VERBOSE else None
		CS_L2_mds = cryosat_baseline_D(full_filename, UNPACK=False)
	elif (fileExtension == '.DBL'):
		#-- Record sizes
		CS_L2_MDS_REC_SIZE = 980
		CS_L2_C_MDS_REC_SIZE = 1392
		#-- check baseline from file to set i_record_size and allocation function
		if (BASELINE == 'C'):
			i_record_size = CS_L2_C_MDS_REC_SIZE
			read_cryosat_variables = cryosat_baseline_C
		else:
			i_record_size = CS_L2_MDS_REC_SIZE
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
			fid = open(full_filename, 'rb')
			cryosat_header = fid.read(j_DS_start)
			#-- iterate through CryoSat file and fill output variables
			CS_L2_mds = read_cryosat_variables(fid,i_record_size,j_num_DSR)
			#-- add headers to output dictionary as METADATA
			CS_L2_mds['METADATA'] = {}
			CS_L2_mds['METADATA']['MPH'] = s_MPH_fields
			CS_L2_mds['METADATA']['SPH'] = s_SPH_fields
			CS_L2_mds['METADATA']['DSD'] = s_DSD_fields
			#-- add absolute orbit number to 1Hz data
			CS_L2_mds['Data_1Hz']['Abs_Orbit']=np.zeros((j_num_DSR),dtype=np.uint32)
			CS_L2_mds['Data_1Hz']['Abs_Orbit'][:]=np.uint32(s_MPH_fields['ABS_ORBIT'])
			#-- add ascending/descending flag to 1Hz data (A=ascending,D=descending)
			CS_L2_mds['Data_1Hz']['Ascending_Flg']=np.zeros((j_num_DSR),dtype=np.bool)
			if (s_SPH_fields['ASCENDING_FLAG'] == 'A'):
				CS_L2_mds['Data_1Hz']['Ascending_Flg'][:] = True
			#-- close the input CryoSat binary file
			fid.close()
		else:
			#-- If there are not MPH/SPH/DSD headers
			#-- extract binary data from input CryoSat data file
			fid = open(full_filename, 'rb')
			#-- iterate through CryoSat file and fill output variables
			CS_L2_mds = read_cryosat_variables(fid,i_record_size,j_num_DSR)
			#-- close the input CryoSat binary file
			fid.close()

	#-- return the data and headers
	return CS_L2_mds
