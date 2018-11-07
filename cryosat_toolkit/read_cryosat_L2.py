#!/usr/bin/env python
u"""
read_cryosat_L2.py
Written by Tyler Sutterley (10/2018)

Reads CryoSat Level-2 data products from baselines A, B and C
Supported CryoSat Modes: LRM, SAR, SARin, FDM, SID, GDR

INPUTS:
	full_filename: full path of CryoSat .DBL file

OUTPUTS:
	Data_1Hz: Time and Orbit Parameters
	Corrections: Elevation Corrections and Flags
	Data_20Hz: Geolocation and Elevation Measurements with Quality Parameters
	METADATA: MPH, SPH and DSD Header data

UPDATE HISTORY:
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
import numpy as np

#-- PURPOSE: Initiate L2 MDS variables for CryoSat Baselines A and B
def cryosat_baseline_AB(fid,record_size,n_records):
	#-- CryoSat-2 1 Hz data fields (Location Group)
	#-- Time and Orbit Parameters plus Measurement Mode
	L2_1Hz_parameters = {}
	#-- Time: day part
	L2_1Hz_parameters['Day'] = np.zeros((n_records),dtype=np.int32)
	#-- Time: second part
	L2_1Hz_parameters['Second'] = np.zeros((n_records),dtype=np.int32)
	#-- Time: microsecond part
	L2_1Hz_parameters['Micsec'] = np.zeros((n_records),dtype=np.int32)
	#-- SIRAL mode
	L2_1Hz_parameters['Siral_mode'] = np.zeros((n_records),dtype=np.uint64)
	#-- Lat_1Hz: packed units (0.1 micro-degree, 1e-7 degrees)
	L2_1Hz_parameters['Lat_1Hz'] = np.zeros((n_records),dtype=np.int32)
	#-- Lon_1Hz: packed units (0.1 micro-degree, 1e-7 degrees)
	L2_1Hz_parameters['Lon_1Hz'] = np.zeros((n_records),dtype=np.int32)
	#-- Alt_1Hz: packed units (mm, 1e-3 m)
	#-- Altitude of COG above reference ellipsoid (interpolated value)
	L2_1Hz_parameters['Alt_1Hz'] = np.zeros((n_records),dtype=np.int32)
	#-- Mispointing: packed units (millidegrees, 1e-3 degrees)
	L2_1Hz_parameters['Mispointing'] = np.zeros((n_records),dtype=np.int16)
	#-- Number of valid records in the block of twenty that contain data
	#-- Last few records of the last block of a dataset may be blank blocks
	#-- inserted to bring the file up to a multiple of twenty.
	L2_1Hz_parameters['N_valid'] = np.zeros((n_records),dtype=np.int16)

	#-- CryoSat-2 geophysical corrections (External Corrections Group)
	L2_final_corrections = {}
	#-- Dry Tropospheric Correction packed units (mm, 1e-3 m)
	L2_final_corrections['dryTrop'] = np.zeros((n_records),dtype=np.int16)
	#-- Wet Tropospheric Correction packed units (mm, 1e-3 m)
	L2_final_corrections['wetTrop'] = np.zeros((n_records),dtype=np.int16)
	#-- Inverse Barometric Correction packed units (mm, 1e-3 m)
	L2_final_corrections['InvBar'] = np.zeros((n_records),dtype=np.int16)
	#-- Dynamic Atmosphere Correction packed units (mm, 1e-3 m)
	L2_final_corrections['DynAtm'] = np.zeros((n_records),dtype=np.int16)
	#-- Ionospheric Correction packed units (mm, 1e-3 m)
	L2_final_corrections['Iono'] = np.zeros((n_records),dtype=np.int16)
	#-- Sea State Bias Correction packed units (mm, 1e-3 m)
	L2_final_corrections['SSB'] = np.zeros((n_records),dtype=np.int16)
	#-- Ocean tide Correction packed units (mm, 1e-3 m)
	L2_final_corrections['ocTideElv'] = np.zeros((n_records),dtype=np.int16)
	#-- Long period equilibrium ocean tide Correction packed units (mm, 1e-3 m)
	L2_final_corrections['lpeTideElv'] = np.zeros((n_records),dtype=np.int16)
	#-- Ocean loading tide Correction packed units (mm, 1e-3 m)
	L2_final_corrections['olTideElv'] = np.zeros((n_records),dtype=np.int16)
	#-- Solid Earth tide Correction packed units (mm, 1e-3 m)
	L2_final_corrections['seTideElv'] = np.zeros((n_records),dtype=np.int16)
	#-- Geocentric Polar tide Correction packed units (mm, 1e-3 m)
	L2_final_corrections['gpTideElv'] = np.zeros((n_records),dtype=np.int16)
	L2_final_corrections['Spare1'] = np.zeros((n_records),dtype=np.int16)
	#-- Surface Type: Packed in groups of three bits for each of the 20 records
	L2_final_corrections['Surf_type'] = np.zeros((n_records),dtype=np.uint64)
	#-- Mean Sea Surface or Geoid packed units (mm, 1e-3 m)
	L2_final_corrections['MSS_Geoid'] = np.zeros((n_records),dtype=np.int32)
	#-- Ocean Depth/Land Elevation Model (ODLE) packed units (mm, 1e-3 m)
	L2_final_corrections['ODLE'] = np.zeros((n_records),dtype=np.int32)
	#-- Ice Concentration packed units (%/100)
	L2_final_corrections['Ice_conc'] = np.zeros((n_records),dtype=np.int16)
	#-- Snow Depth packed units (mm, 1e-3 m)
	L2_final_corrections['Snow_depth'] = np.zeros((n_records),dtype=np.int16)
	#-- Snow Density packed units (kg/m^3)
	L2_final_corrections['Snow_density'] = np.zeros((n_records),dtype=np.int16)
	L2_final_corrections['Spare2'] = np.zeros((n_records),dtype=np.int16)
	#-- Corrections Status Flag
	L2_final_corrections['C_status'] = np.zeros((n_records),dtype=np.uint32)
	#-- Significant Wave Height (SWH) packed units (mm, 1e-3)
	L2_final_corrections['SWH'] = np.zeros((n_records),dtype=np.int16)
	#-- Wind Speed packed units (mm/s, 1e-3 m/s)
	L2_final_corrections['Wind_speed'] = np.zeros((n_records),dtype=np.uint16)
	L2_final_corrections['Spare3'] = np.zeros((n_records),dtype=np.int16)
	L2_final_corrections['Spare4'] = np.zeros((n_records),dtype=np.int16)
	L2_final_corrections['Spare5'] = np.zeros((n_records),dtype=np.int16)
	L2_final_corrections['Spare6'] = np.zeros((n_records),dtype=np.int16)

	#-- CryoSat-2 20 Hz data fields (Measurement Group)
	#-- Derived from instrument measurement parameters
	n_blocks = 20
	L2_final_measurements = {}
	#-- Delta between the timestamps for 20Hz record and the 1Hz record
	#-- D_time_mics packed units (microseconds)
	L2_final_measurements['D_time_mics'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Lat: packed units (0.1 micro-degree, 1e-7 degrees)
	L2_final_measurements['Lat'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Lon: packed units (0.1 micro-degree, 1e-7 degrees)
	L2_final_measurements['Lon'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Measured elevation above ellipsoid from retracker: packed units (mm, 1e-3 m)
	L2_final_measurements['Elev'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Interpolated Sea Surface Height Anomaly: packed units (mm, 1e-3 m)
	L2_final_measurements['SSHA_interp'] = np.zeros((n_records,n_blocks),dtype=np.int16)
	#-- Interpolated Sea Surface Height measurement count
	L2_final_measurements['SSHA_num'] = np.zeros((n_records,n_blocks),dtype=np.int16)
	#-- Interpolation quality estimate RSS: packed units (mm, 1e-3 m)
	L2_final_measurements['SSHA_qual'] = np.zeros((n_records,n_blocks),dtype=np.int16)
	#-- Sigma Zero Backscatter for retracker: packed units (1e-2 dB)
	L2_final_measurements['Sig0'] = np.zeros((n_records,n_blocks),dtype=np.int16)
	#-- Peakiness: packed units (1e-2)
	L2_final_measurements['Peakiness'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
	#-- Freeboard: packed units (mm, 1e-3 m)
	#-- -9999 default value indicates computation has not been performed
	L2_final_measurements['Freeboard'] = np.zeros((n_records,n_blocks),dtype=np.int16)
	#-- Number of averaged echoes or beams
	L2_final_measurements['N_avg'] = np.zeros((n_records,n_blocks),dtype=np.int16)
	L2_final_measurements['Spare1'] = np.zeros((n_records,n_blocks),dtype=np.int16)
	#-- Quality flags
	L2_final_measurements['Quality_Flg'] = np.zeros((n_records,n_blocks),dtype=np.uint32)
	L2_final_measurements['Spare2'] = np.zeros((n_records,n_blocks),dtype=np.int16)
	L2_final_measurements['Spare3'] = np.zeros((n_records,n_blocks),dtype=np.int16)
	L2_final_measurements['Spare4'] = np.zeros((n_records,n_blocks),dtype=np.int16)
	L2_final_measurements['Spare5'] = np.zeros((n_records,n_blocks),dtype=np.int16)
	#-- for each record in the CryoSat file
	for r in range(n_records):
		#-- CryoSat-2 Location Group for record r
		L2_1Hz_parameters['Day'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2_1Hz_parameters['Second'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2_1Hz_parameters['Micsec'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2_1Hz_parameters['Siral_mode'][r] = np.fromfile(fid,dtype='>u8',count=1)
		L2_1Hz_parameters['Lat_1Hz'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2_1Hz_parameters['Lon_1Hz'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2_1Hz_parameters['Alt_1Hz'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2_1Hz_parameters['Mispointing'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_1Hz_parameters['N_valid'][r] = np.fromfile(fid,dtype='>i2',count=1)
		#-- CryoSat-2 External Corrections Group for record r
		L2_final_corrections['dryTrop'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_final_corrections['wetTrop'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_final_corrections['InvBar'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_final_corrections['DynAtm'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_final_corrections['Iono'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_final_corrections['SSB'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_final_corrections['ocTideElv'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_final_corrections['lpeTideElv'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_final_corrections['olTideElv'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_final_corrections['seTideElv'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_final_corrections['gpTideElv'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_final_corrections['Spare1'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_final_corrections['Surf_type'][r] = np.fromfile(fid,dtype='>u8',count=1)
		L2_final_corrections['MSS_Geoid'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2_final_corrections['ODLE'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2_final_corrections['Ice_conc'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_final_corrections['Snow_depth'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_final_corrections['Snow_density'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_final_corrections['Spare2'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_final_corrections['C_status'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2_final_corrections['SWH'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_final_corrections['Wind_speed'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2_final_corrections['Spare3'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_final_corrections['Spare4'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_final_corrections['Spare5'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_final_corrections['Spare6'][r] = np.fromfile(fid,dtype='>i2',count=1)
		#-- CryoSat-2 Measurements Group for record r and block b
		for b in range(n_blocks):
			L2_final_measurements['D_time_mics'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L2_final_measurements['Lat'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L2_final_measurements['Lon'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L2_final_measurements['Elev'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L2_final_measurements['SSHA_interp'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
			L2_final_measurements['SSHA_num'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
			L2_final_measurements['SSHA_qual'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
			L2_final_measurements['Sig0'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
			L2_final_measurements['Peakiness'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
			L2_final_measurements['Freeboard'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
			L2_final_measurements['N_avg'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
			L2_final_measurements['Spare1'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
			L2_final_measurements['Quality_Flg'][r,b] = np.fromfile(fid,dtype='>u4',count=1)
			L2_final_measurements['Spare2'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
			L2_final_measurements['Spare3'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
			L2_final_measurements['Spare4'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
			L2_final_measurements['Spare5'][r,b] = np.fromfile(fid,dtype='>i2',count=1)

	#-- Bind all the bits of the l2_mds together into a single dictionary
	CS_l2_mds = {}
	CS_l2_mds['Data_1Hz'] = L2_1Hz_parameters
	CS_l2_mds['Corrections'] = L2_final_corrections
	CS_l2_mds['Data_20Hz'] = L2_final_measurements
	#-- return the output dictionary
	return CS_l2_mds

#-- PURPOSE: Initiate L2 MDS variables for CryoSat Baseline C
def cryosat_baseline_C(fid,record_size,n_records):
	#-- CryoSat-2 1 Hz data fields (Location Group)
	#-- Time and Orbit Parameters plus Measurement Mode
	L2_c_1Hz_parameters = {}
	#-- Time: day part
	L2_c_1Hz_parameters['Day'] = np.zeros((n_records),dtype=np.int32)
	#-- Time: second part
	L2_c_1Hz_parameters['Second'] = np.zeros((n_records),dtype=np.int32)
	#-- Time: microsecond part
	L2_c_1Hz_parameters['Micsec'] = np.zeros((n_records),dtype=np.int32)
	#-- SIRAL mode
	L2_c_1Hz_parameters['Siral_mode'] = np.zeros((n_records),dtype=np.uint64)
	#-- Lat_1Hz: packed units (0.1 micro-degree, 1e-7 degrees)
	L2_c_1Hz_parameters['Lat_1Hz'] = np.zeros((n_records),dtype=np.int32)
	#-- Lon_1Hz: packed units (0.1 micro-degree, 1e-7 degrees)
	L2_c_1Hz_parameters['Lon_1Hz'] = np.zeros((n_records),dtype=np.int32)
	#-- Alt_1Hz: packed units (mm, 1e-3 m)
	#-- Altitude of COG above reference ellipsoid (interpolated value)
	L2_c_1Hz_parameters['Alt_1Hz'] = np.zeros((n_records),dtype=np.int32)
	#-- Roll: packed units (0.1 micro-degree, 1e-7 degrees)
	L2_c_1Hz_parameters['Roll'] = np.zeros((n_records),dtype=np.int32)
	#-- Pitch: packed units (0.1 micro-degree, 1e-7 degrees)
	L2_c_1Hz_parameters['Pitch'] = np.zeros((n_records),dtype=np.int32)
	#-- Yaw: packed units (0.1 micro-degree, 1e-7 degrees)
	L2_c_1Hz_parameters['Yaw'] = np.zeros((n_records),dtype=np.int32)
	L2_c_1Hz_parameters['Spare'] = np.zeros((n_records),dtype=np.int16)
	#-- Number of valid records in the block of twenty that contain data
	#-- Last few records of the last block of a dataset may be blank blocks
	#-- inserted to bring the file up to a multiple of twenty.
	L2_c_1Hz_parameters['N_valid'] = np.zeros((n_records),dtype=np.int16)

	#-- CryoSat-2 geophysical corrections (External Corrections Group)
	L2_c_final_corrections = {}
	#-- Dry Tropospheric Correction packed units (mm, 1e-3 m)
	L2_c_final_corrections['dryTrop'] = np.zeros((n_records),dtype=np.int16)
	#-- Wet Tropospheric Correction packed units (mm, 1e-3 m)
	L2_c_final_corrections['wetTrop'] = np.zeros((n_records),dtype=np.int16)
	#-- Inverse Barometric Correction packed units (mm, 1e-3 m)
	L2_c_final_corrections['InvBar'] = np.zeros((n_records),dtype=np.int16)
	#-- Dynamic Atmosphere Correction packed units (mm, 1e-3 m)
	L2_c_final_corrections['DynAtm'] = np.zeros((n_records),dtype=np.int16)
	#-- Ionospheric Correction packed units (mm, 1e-3 m)
	L2_c_final_corrections['Iono'] = np.zeros((n_records),dtype=np.int16)
	#-- Sea State Bias Correction packed units (mm, 1e-3 m)
	L2_c_final_corrections['SSB'] = np.zeros((n_records),dtype=np.int16)
	#-- Ocean tide Correction packed units (mm, 1e-3 m)
	L2_c_final_corrections['ocTideElv'] = np.zeros((n_records),dtype=np.int16)
	#-- Long period equilibrium ocean tide Correction packed units (mm, 1e-3 m)
	L2_c_final_corrections['lpeTideElv'] = np.zeros((n_records),dtype=np.int16)
	#-- Ocean loading tide Correction packed units (mm, 1e-3 m)
	L2_c_final_corrections['olTideElv'] = np.zeros((n_records),dtype=np.int16)
	#-- Solid Earth tide Correction packed units (mm, 1e-3 m)
	L2_c_final_corrections['seTideElv'] = np.zeros((n_records),dtype=np.int16)
	#-- Geocentric Polar tide Correction packed units (mm, 1e-3 m)
	L2_c_final_corrections['gpTideElv'] = np.zeros((n_records),dtype=np.int16)
	L2_c_final_corrections['Spare1'] = np.zeros((n_records),dtype=np.int16)
	#-- Surface Type: Packed in groups of three bits for each of the 20 records
	L2_c_final_corrections['Surf_type'] = np.zeros((n_records),dtype=np.uint64)
	#-- Mean Sea Surface or Geoid packed units (mm, 1e-3 m)
	L2_c_final_corrections['MSS_Geoid'] = np.zeros((n_records),dtype=np.int32)
	#-- Ocean Depth/Land Elevation Model (ODLE) packed units (mm, 1e-3 m)
	L2_c_final_corrections['ODLE'] = np.zeros((n_records),dtype=np.int32)
	#-- Ice Concentration packed units (%/100)
	L2_c_final_corrections['Ice_conc'] = np.zeros((n_records),dtype=np.int16)
	#-- Snow Depth packed units (mm, 1e-3 m)
	L2_c_final_corrections['Snow_depth'] = np.zeros((n_records),dtype=np.int16)
	#-- Snow Density packed units (kg/m^3)
	L2_c_final_corrections['Snow_density'] = np.zeros((n_records),dtype=np.int16)
	L2_c_final_corrections['Spare2'] = np.zeros((n_records),dtype=np.int16)
	#-- Corrections Status Flag
	L2_c_final_corrections['C_status'] = np.zeros((n_records),dtype=np.uint32)
	#-- Significant Wave Height (SWH) packed units (mm, 1e-3)
	L2_c_final_corrections['SWH'] = np.zeros((n_records),dtype=np.int16)
	#-- Wind Speed packed units (mm/s, 1e-3 m/s)
	L2_c_final_corrections['Wind_speed'] = np.zeros((n_records),dtype=np.uint16)
	L2_c_final_corrections['Spare3'] = np.zeros((n_records),dtype=np.int16)
	L2_c_final_corrections['Spare4'] = np.zeros((n_records),dtype=np.int16)
	L2_c_final_corrections['Spare5'] = np.zeros((n_records),dtype=np.int16)
	L2_c_final_corrections['Spare6'] = np.zeros((n_records),dtype=np.int16)

	#-- CryoSat-2 20 Hz data fields (Measurement Group)
	#-- Derived from instrument measurement parameters
	n_blocks = 20
	L2_c_final_measurements = {}
	#-- Delta between the timestamps for 20Hz record and the 1Hz record
	#-- D_time_mics packed units (microseconds)
	L2_c_final_measurements['D_time_mics'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Lat: packed units (0.1 micro-degree, 1e-7 degrees)
	L2_c_final_measurements['Lat'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Lon: packed units (0.1 micro-degree, 1e-7 degrees)
	L2_c_final_measurements['Lon'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Measured elevation above ellipsoid from retracker 1: packed units (mm, 1e-3 m)
	L2_c_final_measurements['Elev_1'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Measured elevation above ellipsoid from retracker 2: packed units (mm, 1e-3 m)
	L2_c_final_measurements['Elev_2'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Measured elevation above ellipsoid from retracker 3: packed units (mm, 1e-3 m)
	L2_c_final_measurements['Elev_3'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Sigma Zero Backscatter for retracker 1: packed units (1e-2 dB)
	L2_c_final_measurements['Sig0_1'] = np.zeros((n_records,n_blocks),dtype=np.int16)
	#-- Sigma Zero Backscatter for retracker 2: packed units (1e-2 dB)
	L2_c_final_measurements['Sig0_2'] = np.zeros((n_records,n_blocks),dtype=np.int16)
	#-- Sigma Zero Backscatter for retracker 3: packed units (1e-2 dB)
	L2_c_final_measurements['Sig0_3'] = np.zeros((n_records,n_blocks),dtype=np.int16)
	#-- Freeboard: packed units (mm, 1e-3 m)
	#-- -9999 default value indicates computation has not been performed
	L2_c_final_measurements['Freeboard'] = np.zeros((n_records,n_blocks),dtype=np.int16)
	#-- Interpolated Sea Surface Height Anomaly: packed units (mm, 1e-3 m)
	L2_c_final_measurements['SSHA_interp'] = np.zeros((n_records,n_blocks),dtype=np.int16)
	#-- Interpolated Sea Surface Height measurement count
	L2_c_final_measurements['SSHA_num'] = np.zeros((n_records,n_blocks),dtype=np.int16)
	#-- Interpolation quality estimate RSS: packed units (mm, 1e-3 m)
	L2_c_final_measurements['SSHA_qual'] = np.zeros((n_records,n_blocks),dtype=np.int16)
	#-- Peakiness: packed units (1e-2)
	L2_c_final_measurements['Peakiness'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
	#-- Number of averaged echoes or beams
	L2_c_final_measurements['N_avg'] = np.zeros((n_records,n_blocks),dtype=np.int16)
	L2_c_final_measurements['Spare1'] = np.zeros((n_records,n_blocks),dtype=np.int16)
	#-- Quality flags
	L2_c_final_measurements['Quality_Flg'] = np.zeros((n_records,n_blocks),dtype=np.uint32)
	#-- Corrections Application Flag
	L2_c_final_measurements['Corrections_Flg'] = np.zeros((n_records,n_blocks),dtype=np.uint32)
	#-- Quality metric for retracker 1
	L2_c_final_measurements['Quality_1'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Quality metric for retracker 2
	L2_c_final_measurements['Quality_2'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Quality metric for retracker 3
	L2_c_final_measurements['Quality_3'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- for each record in the CryoSat file
	for r in range(n_records):
		#-- CryoSat-2 Location Group for record r
		L2_c_1Hz_parameters['Day'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2_c_1Hz_parameters['Second'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2_c_1Hz_parameters['Micsec'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2_c_1Hz_parameters['Siral_mode'][r] = np.fromfile(fid,dtype='>u8',count=1)
		L2_c_1Hz_parameters['Lat_1Hz'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2_c_1Hz_parameters['Lon_1Hz'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2_c_1Hz_parameters['Alt_1Hz'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2_c_1Hz_parameters['Roll'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2_c_1Hz_parameters['Pitch'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2_c_1Hz_parameters['Yaw'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2_c_1Hz_parameters['Spare'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_c_1Hz_parameters['N_valid'][r] = np.fromfile(fid,dtype='>i2',count=1)
		#-- CryoSat-2 External Corrections Group for record r
		L2_c_final_corrections['dryTrop'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_c_final_corrections['wetTrop'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_c_final_corrections['InvBar'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_c_final_corrections['DynAtm'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_c_final_corrections['Iono'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_c_final_corrections['SSB'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_c_final_corrections['ocTideElv'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_c_final_corrections['lpeTideElv'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_c_final_corrections['olTideElv'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_c_final_corrections['seTideElv'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_c_final_corrections['gpTideElv'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_c_final_corrections['Spare1'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_c_final_corrections['Surf_type'][r] = np.fromfile(fid,dtype='>u8',count=1)
		L2_c_final_corrections['MSS_Geoid'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2_c_final_corrections['ODLE'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2_c_final_corrections['Ice_conc'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_c_final_corrections['Snow_depth'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_c_final_corrections['Snow_density'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_c_final_corrections['Spare2'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_c_final_corrections['C_status'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2_c_final_corrections['SWH'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_c_final_corrections['Wind_speed'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2_c_final_corrections['Spare3'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_c_final_corrections['Spare4'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_c_final_corrections['Spare5'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2_c_final_corrections['Spare6'][r] = np.fromfile(fid,dtype='>i2',count=1)
		#-- CryoSat-2 Measurements Group for record r and block b
		for b in range(n_blocks):
			L2_c_final_measurements['D_time_mics'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L2_c_final_measurements['Lat'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L2_c_final_measurements['Lon'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L2_c_final_measurements['Elev_1'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L2_c_final_measurements['Elev_2'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L2_c_final_measurements['Elev_3'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L2_c_final_measurements['Sig0_1'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
			L2_c_final_measurements['Sig0_2'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
			L2_c_final_measurements['Sig0_3'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
			L2_c_final_measurements['Freeboard'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
			L2_c_final_measurements['SSHA_interp'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
			L2_c_final_measurements['SSHA_num'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
			L2_c_final_measurements['SSHA_qual'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
			L2_c_final_measurements['Peakiness'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
			L2_c_final_measurements['N_avg'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
			L2_c_final_measurements['Spare1'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
			L2_c_final_measurements['Quality_Flg'][r,b] = np.fromfile(fid,dtype='>u4',count=1)
			L2_c_final_measurements['Corrections_Flg'][r,b] = np.fromfile(fid,dtype='>u4',count=1)
			L2_c_final_measurements['Quality_1'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L2_c_final_measurements['Quality_2'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L2_c_final_measurements['Quality_3'][r,b] = np.fromfile(fid,dtype='>i4',count=1)

	#-- Bind all the bits of the l2_mds together into a single dictionary
	CS_l2_c_mds = {}
	CS_l2_c_mds['Data_1Hz'] = L2_c_1Hz_parameters
	CS_l2_c_mds['Corrections'] = L2_c_final_corrections
	CS_l2_c_mds['Data_20Hz'] = L2_c_final_measurements
	#-- return the output dictionary
	return CS_l2_c_mds

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
	regex_patterns.append(b'DS_NAME\="SIR_LRM_L2[\s+]*"')
	regex_patterns.append(b'DS_NAME\="SIR_SAR_L2B[\s+]*"')
	regex_patterns.append(b'DS_NAME\="SIR_SAR_L2[\s+]*"')
	regex_patterns.append(b'DS_NAME\="SIR_FDM_L2[\s+]*"')
	regex_patterns.append(b'DS_NAME\="SIR_SARIL2B[\s+]*"')
	regex_patterns.append(b'DS_NAME\="SIR_SARIL2[\s+]*"')
	regex_patterns.append(b'DS_NAME\="SIR_SAR_L2B_I[\s+]*"')
	regex_patterns.append(b'DS_NAME\="SIR_SAR_L2A[\s+]*"')
	regex_patterns.append(b'DS_NAME\="SIR_SIN_L2[\s+]*"')
	regex_patterns.append(b'DS_NAME\="SIR_SID_L2[\s+]*"')
	regex_patterns.append(b'DS_NAME\="SIR_LRMIL2[\s+]*"')
	regex_patterns.append(b'DS_NAME\="SIR_LRM_L2_I[\s+]*"')
	regex_patterns.append(b'DS_NAME\="SIR_SARIL2A[\s+]*"')
	regex_patterns.append(b'DS_NAME\="SIR_SAR_L2A_I[\s+]*"')
	regex_patterns.append(b'DS_NAME\="SIR_SAR_L2_I[\s+]*"')
	regex_patterns.append(b'DS_NAME\="SIR_SINIL2[\s+]*"')
	regex_patterns.append(b'DS_NAME\="SIR_SIN_L2_I[\s+]*"')
	regex_patterns.append(b'DS_NAME\="SIR_SIDIL2[\s+]*"')
	regex_patterns.append(b'DS_NAME\="SIR_SID_L2_I[\s+]*"')
	regex_patterns.append(b'DS_NAME\="SIR_GDR_2A[\s+]*"')
	regex_patterns.append(b'DS_NAME\="SIR_GDR_2B[\s+]*"')
	regex_patterns.append(b'DS_NAME\="SIR_GDR_2[\s+]*"')
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
	regex_pattern = '(.*?)_({0})_({1})__(\d+T?\d+)_(\d+T?\d+)_(.*?)(\d+)'.format(
		regex_class, regex_products)
	rx = re.compile(regex_pattern, re.VERBOSE)
	#-- extract file information from filename
	MI,CLASS,PRODUCT,START,STOP,BASELINE,VERSION=rx.findall(fileBasename).pop()
	#-- Extract Date information
	start_yr,start_mon,start_day=np.array([START[:4],START[4:6],START[6:8]],dtype=np.uint16)
	start_hh,start_mm,start_ss=np.array([START[-6:-4],START[-4:-2],START[-2:]],dtype=np.uint8)
	stop_yr,stop_mon,stop_day=np.array([STOP[:4],STOP[4:6],STOP[6:8]],dtype=np.uint16)
	stop_hh,stop_mm,stop_ss=np.array([STOP[-6:-4],STOP[-4:-2],STOP[-2:]],dtype=np.uint8)

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
