#!/usr/bin/env python
u"""
read_cryosat_L2I.py
Written by Tyler Sutterley (10/2018)

Reads CryoSat Level-2 Intermediate data products from baselines A, B, BC and C
Supported CryoSat Modes: LRM, SAR, SARin, FDM, SID, GDR

INPUTS:
	full_filename: full path of CryoSat .DBL file

OUTPUTS:
	Location: Time and Orbit Parameters
	Geometry: Elevation Corrections and Flags
	Data: Geolocation and Elevation Measurements with Quality Parameters
	Auxiliary: Auxiliary Data for Elevation Processing
	Instrumental: Intrument Corrections
	METADATA: MPH, SPH and DSD Header data

UPDATE HISTORY:
Updated 10/2018: updated header read functions for python3
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
	#-- CryoSat-2 Location Group
	#-- Time and Orbit Parameters plus Measurement Mode
	L2I_location_parameters = {}
	#-- Time: day part
	L2I_location_parameters['Day'] = np.zeros((n_records),dtype=np.int32)
	#-- Time: second part
	L2I_location_parameters['Sec'] = np.zeros((n_records),dtype=np.uint32)
	#-- Time: microsecond part
	L2I_location_parameters['Micsec'] = np.zeros((n_records),dtype=np.uint32)
	#-- USO correction factor
	L2I_location_parameters['USO_Corr'] = np.zeros((n_records),dtype=np.int32)
	#-- Mode ID
	L2I_location_parameters['Mode_ID'] = np.zeros((n_records),dtype=np.uint16)
	#-- Source sequence counter
	L2I_location_parameters['SSC'] = np.zeros((n_records),dtype=np.uint16)
	#-- Instrument configuration
	L2I_location_parameters['Inst_config'] = np.zeros((n_records),dtype=np.uint32)
	#-- Record Counter
	L2I_location_parameters['Rec_Count'] = np.zeros((n_records),dtype=np.uint32)
	#-- Lat: packed units (0.1 micro-degree, 1e-7 degrees)
	L2I_location_parameters['Lat'] = np.zeros((n_records),dtype=np.int32)
	#-- Lon: packed units (0.1 micro-degree, 1e-7 degrees)
	L2I_location_parameters['Lon'] = np.zeros((n_records),dtype=np.int32)
	#-- Alt: packed units (mm, 1e-3 m)
	#-- Altitude of COG above reference ellipsoid (interpolated value)
	L2I_location_parameters['Alt'] = np.zeros((n_records),dtype=np.int32)
	#-- Instantaneous altitude rate derived from orbit: packed units (mm/s, 1e-3 m/s)
	L2I_location_parameters['Alt_rate'] = np.zeros((n_records),dtype=np.int32)
	#-- Satellite velocity vector. In ITRF: packed units (mm/s, 1e-3 m/s)
	L2I_location_parameters['Sat_velocity'] = np.zeros((n_records,3),dtype=np.int32)
	#-- Real beam direction vector. In CRF: packed units (micro-m/s, 1e-6 m/s)
	L2I_location_parameters['Real_beam'] = np.zeros((n_records,3),dtype=np.int32)
	#-- Interferometer baseline vector. In CRF: packed units (micro-m/s, 1e-6 m/s)
	L2I_location_parameters['Baseline'] = np.zeros((n_records,3),dtype=np.int32)
	#-- Measurement Confidence Data
	L2I_location_parameters['L2_MCD'] = np.zeros((n_records),dtype=np.uint32)

	#-- CryoSat-2 Measurement Group
	#-- Derived from instrument measurement parameters
	L2I_measurements = {}
	#-- Measured elevation above ellipsoid from retracker: packed units (mm, 1e-3 m)
	L2I_measurements['Elev'] = np.zeros((n_records),dtype=np.int32)
	#-- Sigma Zero Backscatter for retracker: packed units (1e-2 dB)
	L2I_measurements['Sig0'] = np.zeros((n_records),dtype=np.int32)
	#-- SWH packed units (mm, 1e-3)
	L2I_measurements['SWH'] = np.zeros((n_records),dtype=np.int32)
	#-- Peakiness: packed units (1e-2)
	L2I_measurements['Peakiness'] = np.zeros((n_records),dtype=np.int32)
	#-- Retracked range correction: packed units (mm, 1e-3 m)
	L2I_measurements['Retrack_range'] = np.zeros((n_records),dtype=np.int32)
	#-- Retracked sigma 0 correction: packed units (1e-2 dB)
	L2I_measurements['Retrack_sig0'] = np.zeros((n_records),dtype=np.int32)
	#-- Retrackers 3-13 output
	L2I_measurements['Retrack_3'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_4'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_5'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_6'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_7'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_8'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_9'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_10'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_11'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_12'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_13'] = np.zeros((n_records),dtype=np.int32)
	#-- Power echo shape parameter: packed units (dB/100)
	L2I_measurements['echo_shape'] = np.zeros((n_records),dtype=np.int32)
	#-- Beam behaviour parameter: unitless code number related to
	#-- surface characteristics
	L2I_measurements['BB_parameter'] = np.zeros((n_records,50),dtype=np.int16)
	#-- Cross track angle: packed units (micro radians)
	L2I_measurements['X_Track_Angle'] = np.zeros((n_records),dtype=np.int32)
	#-- Leading edge coherence at retrack point 1/1000
	L2I_measurements['Coherence'] = np.zeros((n_records),dtype=np.int32)
	#-- Interpolated Ocean Height: packed units (mm above ellipsoid)
	L2I_measurements['Ocean_ht'] = np.zeros((n_records),dtype=np.int32)
	#-- Freeboard: packed units (mm, 1e-3 m)
	#-- -9999 default value indicates computation has not been performed
	L2I_measurements['Freeboard'] = np.zeros((n_records),dtype=np.int32)
	#-- Surface Height Anomaly: packed units (mm, 1e-3 m)
	L2I_measurements['SHA'] = np.zeros((n_records),dtype=np.int32)
	#-- Interpolated Surface Height Anomaly: packed units (mm, 1e-3 m)
	L2I_measurements['SHA_interp'] = np.zeros((n_records),dtype=np.int32)
	#-- Error in ocean height interpolation: packed units (mm, 1e-3 m)
	L2I_measurements['Interp_Err'] = np.zeros((n_records),dtype=np.uint16)
	#-- Number of forward records interpolated
	L2I_measurements['Interp_Count_Fwd'] = np.zeros((n_records),dtype=np.uint16)
	#-- Number of backward records interpolated
	L2I_measurements['Interp_Count_Bkwd'] = np.zeros((n_records),dtype=np.uint16)
	#-- Distance in time of most forward record interpolated (milli-seconds)
	L2I_measurements['Interp_Time_Fwd'] = np.zeros((n_records),dtype=np.uint16)
	#-- Distance in time of most backward record interpolated (milli-seconds)
	L2I_measurements['Interp_Time_Bkwd'] = np.zeros((n_records),dtype=np.uint16)
	#-- Interpolation error flag
	L2I_measurements['Interp_Error_Flg'] = np.zeros((n_records),dtype=np.uint16)
	#-- Measurement mode
	L2I_measurements['Meas_Mode'] = np.zeros((n_records),dtype=np.uint32)
	#-- Quality flags
	L2I_measurements['Quality_Flg'] = np.zeros((n_records),dtype=np.uint32)
	#-- Retracker flags
	L2I_measurements['Retracker_Flg'] = np.zeros((n_records),dtype=np.uint32)
	#-- Height calculation details
	#-- Specifies what was applied during the height calculation
	L2I_measurements['Height_status'] = np.zeros((n_records),dtype=np.uint32)
	#-- SAR freeboard status flag
	L2I_measurements['Freeboard_status'] = np.zeros((n_records),dtype=np.uint32)
	#-- Number of averaged echoes or beams
	L2I_measurements['N_avg'] = np.zeros((n_records),dtype=np.uint16)
	#-- Wind Speed packed units (mm/s, 1e-3 m/s)
	L2I_measurements['Wind_speed'] = np.zeros((n_records),dtype=np.uint16)
	L2I_measurements['Spares1'] = np.zeros((n_records,3),dtype=np.int32)

	#-- CryoSat-2 Auxiliary Data Group
	L2I_aux_data = {}
	#-- Ice Concentration packed units (%/1000)
	L2I_aux_data['Ice_conc'] = np.zeros((n_records),dtype=np.int32)
	#-- Snow Depth packed units (mm, 1e-3 m)
	L2I_aux_data['Snow_depth'] = np.zeros((n_records),dtype=np.int32)
	#-- Snow Density packed units (kg/m^3)
	L2I_aux_data['Snow_density'] = np.zeros((n_records),dtype=np.int32)
	#-- Discriminator result
	L2I_aux_data['Discriminator'] = np.zeros((n_records),dtype=np.int32)
	#-- SARin discriminator parameters 1-10
	L2I_aux_data['SARIN_disc_1'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['SARIN_disc_2'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['SARIN_disc_3'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['SARIN_disc_4'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['SARIN_disc_5'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['SARIN_disc_6'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['SARIN_disc_7'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['SARIN_disc_8'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['SARIN_disc_9'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['SARIN_disc_10'] = np.zeros((n_records),dtype=np.int32)
	#-- Discriminator flags
	L2I_aux_data['Discrim_Flg'] = np.zeros((n_records),dtype=np.uint32)
	#-- Slope model correction (Attitude of echo in micro-degrees)
	L2I_aux_data['Attitude'] = np.zeros((n_records),dtype=np.int32)
	#-- Slope model correction (Azimuth of echo in micro-degrees)
	L2I_aux_data['Azimuth'] = np.zeros((n_records),dtype=np.int32)
	#-- The original latitude of the satellite (micro-degrees)
	L2I_aux_data['Lat_sat'] = np.zeros((n_records),dtype=np.int32)
	#-- The original longitude of the satellite (micro-degrees)
	L2I_aux_data['Lon_sat'] = np.zeros((n_records),dtype=np.int32)
	#-- Ambiguity indicator
	L2I_aux_data['Ambiguity'] = np.zeros((n_records),dtype=np.uint32)
	#-- Mean Sea Surface standard Model: packed units (mm, 1e-3 m)
	L2I_aux_data['MSS_model'] = np.zeros((n_records),dtype=np.int32)
	#-- Geoid standard Model: packed units (mm, 1e-3 m)
	L2I_aux_data['Geoid_model'] = np.zeros((n_records),dtype=np.int32)
	#-- ODLE standard Model: packed units (mm, 1e-3 m)
	L2I_aux_data['ODLE_model'] = np.zeros((n_records),dtype=np.int32)
	#-- The interpolated elevation value obtained from the DEM (mm)
	L2I_aux_data['DEM_elev'] = np.zeros((n_records),dtype=np.int32)
	#-- Identification of DEM used in SARin ambiguity test
	L2I_aux_data['DEM_ID'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['Spares2'] = np.zeros((n_records,4),dtype=np.int32)

	#-- CryoSat-2 External Corrections Group
	L2I_geo_corrections = {}
	#-- Dry Tropospheric Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['dryTrop'] = np.zeros((n_records),dtype=np.int32)
	#-- Wet Tropospheric Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['wetTrop'] = np.zeros((n_records),dtype=np.int32)
	#-- Inverse Barometric Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['InvBar'] = np.zeros((n_records),dtype=np.int32)
	#-- Delta Inverse Barometric Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['DAC'] = np.zeros((n_records),dtype=np.int32)
	#-- GIM Ionospheric Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['Iono_GIM'] = np.zeros((n_records),dtype=np.int32)
	#-- Model Ionospheric Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['Iono_model'] = np.zeros((n_records),dtype=np.int32)
	#-- Ocean tide Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['ocTideElv'] = np.zeros((n_records),dtype=np.int32)
	#-- Long period equilibrium ocean tide Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['lpeTideElv'] = np.zeros((n_records),dtype=np.int32)
	#-- Ocean loading tide Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['olTideElv'] = np.zeros((n_records),dtype=np.int32)
	#-- Solid Earth tide Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['seTideElv'] = np.zeros((n_records),dtype=np.int32)
	#-- Geocentric Polar tide Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['gpTideElv'] = np.zeros((n_records),dtype=np.int32)
	#-- Surface Type: Packed in groups of three bits for each of the 20 records
	L2I_geo_corrections['Surf_type'] = np.zeros((n_records),dtype=np.uint32)
	#-- Corrections Status Flag
	L2I_geo_corrections['Corr_status'] = np.zeros((n_records),dtype=np.uint32)
	#-- Correction Error Flag
	L2I_geo_corrections['Corr_error'] = np.zeros((n_records),dtype=np.uint32)
	#-- Sea State Bias Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['SSB'] = np.zeros((n_records),dtype=np.int32)
	L2I_geo_corrections['Spares3'] = np.zeros((n_records,2),dtype=np.int32)

	#-- CryoSat-2 Internal Corrections Group
	L2I_inst_corrections = {}
	#-- Doppler range correction: Radial + slope (mm)
	L2I_inst_corrections['Doppler_range'] = np.zeros((n_records),dtype=np.int32)
	#-- Instrument Range Correction: t-r antenna (mm)
	L2I_inst_corrections['TR_inst_range'] = np.zeros((n_records),dtype=np.int32)
	#-- Instrument Range Correction: r-only antenna (mm)
	L2I_inst_corrections['R_inst_range'] = np.zeros((n_records),dtype=np.int32)
	#-- Instrument Sigma 0 Correction: t-r antenna (dB/100)
	L2I_inst_corrections['TR_inst_gain'] = np.zeros((n_records),dtype=np.int32)
	#-- Instrument Sigma 0 Correction: r-only (dB/100)
	L2I_inst_corrections['R_inst_gain'] = np.zeros((n_records),dtype=np.int32)
	#-- Internal Phase Correction (milli-radians)
	L2I_inst_corrections['Internal_phase'] = np.zeros((n_records),dtype=np.int32)
	#-- External Phase Correction (milli-radians)
	L2I_inst_corrections['External_phase'] = np.zeros((n_records),dtype=np.int32)
	#-- Noise Power measurement
	L2I_inst_corrections['Noise_power'] = np.zeros((n_records),dtype=np.int32)
	#-- Phase slope correction (microradians)
	L2I_inst_corrections['Phase_slope'] = np.zeros((n_records),dtype=np.int32)
	L2I_inst_corrections['Spares4'] = np.zeros((n_records,2),dtype=np.int32)

	#-- for each record in the CryoSat file
	for r in range(n_records):
		#-- get satellite time and orbit parameters for record r
		L2I_location_parameters['Day'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_location_parameters['Sec'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_location_parameters['Micsec'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_location_parameters['USO_Corr'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_location_parameters['Mode_ID'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2I_location_parameters['SSC'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2I_location_parameters['Inst_config'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_location_parameters['Rec_Count'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_location_parameters['Lat'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_location_parameters['Lon'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_location_parameters['Alt'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_location_parameters['Alt_rate'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_location_parameters['Sat_velocity'][r,:] = np.fromfile(fid,dtype='>i4',count=3)
		L2I_location_parameters['Real_beam'][r,:] = np.fromfile(fid,dtype='>i4',count=3)
		L2I_location_parameters['Baseline'][r,:] = np.fromfile(fid,dtype='>i4',count=3)
		L2I_location_parameters['L2_MCD'][r] = np.fromfile(fid,dtype='>u4',count=1)

		#-- elevation measurements
		L2I_measurements['Elev'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Sig0'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['SWH'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Peakiness'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_range'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_sig0'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_3'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_4'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_5'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_6'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_7'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_8'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_9'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_10'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_11'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_12'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_13'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['echo_shape'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['BB_parameter'][r,:] = np.fromfile(fid,dtype='>i2',count=50)
		L2I_measurements['X_Track_Angle'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Coherence'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Ocean_ht'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Freeboard'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['SHA'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['SHA_interp'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Interp_Err'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2I_measurements['Interp_Count_Fwd'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2I_measurements['Interp_Count_Bkwd'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2I_measurements['Interp_Time_Fwd'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2I_measurements['Interp_Time_Bkwd'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2I_measurements['Interp_Error_Flg'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2I_measurements['Meas_Mode'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_measurements['Quality_Flg'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_measurements['Retracker_Flg'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_measurements['Height_status'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_measurements['Freeboard_status'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_measurements['N_avg'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2I_measurements['Wind_speed'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2I_measurements['Spares1'][r,:] = np.fromfile(fid,dtype='>i4',count=3)

		#-- Auxiliary Data
		L2I_aux_data['Ice_conc'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Snow_depth'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Snow_density'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Discriminator'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['SARIN_disc_1'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['SARIN_disc_2'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['SARIN_disc_3'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['SARIN_disc_4'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['SARIN_disc_5'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['SARIN_disc_6'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['SARIN_disc_7'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['SARIN_disc_8'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['SARIN_disc_9'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['SARIN_disc_10'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Discrim_Flg'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_aux_data['Attitude'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Azimuth'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Lat_sat'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Lon_sat'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Ambiguity'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_aux_data['MSS_model'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Geoid_model'][r] =np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['ODLE_model'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['DEM_elev'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['DEM_ID'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Spares2'][r,:] = np.fromfile(fid,dtype='>i4',count=4)

		#-- CryoSat-2 External Corrections Group for record r
		L2I_geo_corrections['dryTrop'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['wetTrop'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['InvBar'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['DAC'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['Iono_GIM'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['Iono_model'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['ocTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['lpeTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['olTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['seTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['gpTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['Surf_type'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_geo_corrections['Corr_status'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_geo_corrections['Corr_error'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_geo_corrections['SSB'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['Spares3'][r,:] = np.fromfile(fid,dtype='>i4',count=2)

		#-- CryoSat-2 Internal Corrections Group for record r
		L2I_inst_corrections['Doppler_range'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_inst_corrections['TR_inst_range'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_inst_corrections['R_inst_range'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_inst_corrections['TR_inst_gain'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_inst_corrections['R_inst_gain'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_inst_corrections['Internal_phase'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_inst_corrections['External_phase'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_inst_corrections['Noise_power'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_inst_corrections['Phase_slope'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_inst_corrections['Spares4'][r,:] = np.fromfile(fid,dtype='>i4',count=2)

	#-- Bind all the bits of the l2i_mds together into a single dictionary
	CS_L2I_mds = {}
	CS_L2I_mds['Location'] = L2I_location_parameters
	CS_L2I_mds['Data'] = L2I_measurements
	CS_L2I_mds['Auxiliary'] = L2I_aux_data
	CS_L2I_mds['Geometry'] = L2I_geo_corrections
	CS_L2I_mds['Instrumental'] = L2I_inst_corrections
	#-- return the output dictionary
	return CS_L2I_mds

#-- PURPOSE: Initiate L2 MDS variables for CryoSat Baseline BC
def cryosat_baseline_BC(fid,record_size,n_records):
	#-- CryoSat-2 Location Group
	#-- Time and Orbit Parameters plus Measurement Mode
	L2I_location_parameters = {}
	#-- Time: day part
	L2I_location_parameters['Day'] = np.zeros((n_records),dtype=np.int32)
	#-- Time: second part
	L2I_location_parameters['Sec'] = np.zeros((n_records),dtype=np.uint32)
	#-- Time: microsecond part
	L2I_location_parameters['Micsec'] = np.zeros((n_records),dtype=np.uint32)
	#-- USO correction factor
	L2I_location_parameters['USO_Corr'] = np.zeros((n_records),dtype=np.int32)
	#-- Mode ID
	L2I_location_parameters['Mode_ID'] = np.zeros((n_records),dtype=np.uint16)
	#-- Source sequence counter
	L2I_location_parameters['SSC'] = np.zeros((n_records),dtype=np.uint16)
	#-- Instrument configuration
	L2I_location_parameters['Inst_config'] = np.zeros((n_records),dtype=np.uint32)
	#-- Record Counter
	L2I_location_parameters['Rec_Count'] = np.zeros((n_records),dtype=np.uint32)
	#-- Lat: packed units (0.1 micro-degree, 1e-7 degrees)
	L2I_location_parameters['Lat'] = np.zeros((n_records),dtype=np.int32)
	#-- Lon: packed units (0.1 micro-degree, 1e-7 degrees)
	L2I_location_parameters['Lon'] = np.zeros((n_records),dtype=np.int32)
	#-- Alt: packed units (mm, 1e-3 m)
	#-- Altitude of COG above reference ellipsoid (interpolated value)
	L2I_location_parameters['Alt'] = np.zeros((n_records),dtype=np.int32)
	#-- Instantaneous altitude rate derived from orbit: packed units (mm/s, 1e-3 m/s)
	L2I_location_parameters['Alt_rate'] = np.zeros((n_records),dtype=np.int32)
	#-- Satellite velocity vector. In ITRF: packed units (mm/s, 1e-3 m/s)
	L2I_location_parameters['Sat_velocity'] = np.zeros((n_records,3),dtype=np.int32)
	#-- Real beam direction vector. In CRF: packed units (micro-m/s, 1e-6 m/s)
	L2I_location_parameters['Real_beam'] = np.zeros((n_records,3),dtype=np.int32)
	#-- Interferometer baseline vector. In CRF: packed units (micro-m/s, 1e-6 m/s)
	L2I_location_parameters['Baseline'] = np.zeros((n_records,3),dtype=np.int32)
	#-- Star Tracker ID
	L2I_location_parameters['ST_ID'] = np.zeros((n_records),dtype=np.int16)
	L2I_location_parameters['Spare'] = np.zeros((n_records),dtype=np.int16)
	#-- Roll (Derived from star trackers): packed units (0.1 micro-degree, 1e-7 degrees)
	L2I_location_parameters['Roll'] = np.zeros((n_records),dtype=np.int32)
	#-- Pitch (Derived from star trackers): packed units (0.1 micro-degree, 1e-7 degrees)
	L2I_location_parameters['Pitch'] = np.zeros((n_records),dtype=np.int32)
	#-- Yaw (Derived from star trackers): packed units (0.1 micro-degree, 1e-7 degrees)
	L2I_location_parameters['Yaw'] = np.zeros((n_records),dtype=np.int32)
	#-- Measurement Confidence Data
	L2I_location_parameters['L2_MCD'] = np.zeros((n_records),dtype=np.uint32)

	#-- CryoSat-2 Measurement Group
	#-- Derived from instrument measurement parameters
	L2I_measurements = {}
	#-- Measured elevation above ellipsoid from retracker 1: packed units (mm, 1e-3 m)
	L2I_measurements['Elev_1'] = np.zeros((n_records),dtype=np.int32)
	#-- Measured elevation above ellipsoid from retracker 2: packed units (mm, 1e-3 m)
	L2I_measurements['Elev_2'] = np.zeros((n_records),dtype=np.int32)
	#-- Measured elevation above ellipsoid from retracker 3: packed units (mm, 1e-3 m)
	L2I_measurements['Elev_3'] = np.zeros((n_records),dtype=np.int32)
	#-- Sigma Zero Backscatter for retracker 1: packed units (1e-2 dB)
	L2I_measurements['Sig0_1'] = np.zeros((n_records),dtype=np.int32)
	#-- Sigma Zero Backscatter for retracker 2: packed units (1e-2 dB)
	L2I_measurements['Sig0_2'] = np.zeros((n_records),dtype=np.int32)
	#-- Sigma Zero Backscatter for retracker 3: packed units (1e-2 dB)
	L2I_measurements['Sig0_3'] = np.zeros((n_records),dtype=np.int32)
	#-- SWH packed units (mm, 1e-3)
	L2I_measurements['SWH'] = np.zeros((n_records),dtype=np.int32)
	#-- Peakiness: packed units (1e-2)
	L2I_measurements['Peakiness'] = np.zeros((n_records),dtype=np.int32)
	#-- Retracked range correction for retracker 1: packed units (mm, 1e-3 m)
	L2I_measurements['Retrack_1_range'] = np.zeros((n_records),dtype=np.int32)
	#-- Retracked range correction for retracker 2: packed units (mm, 1e-3 m)
	L2I_measurements['Retrack_2_range'] = np.zeros((n_records),dtype=np.int32)
	#-- Retracked range correction for retracker 3: packed units (mm, 1e-3 m)
	L2I_measurements['Retrack_3_range'] = np.zeros((n_records),dtype=np.int32)
	#-- Retracked sigma 0 correction for Retracker 1: packed units (1e-2 dB)
	L2I_measurements['Retrack_1_sig0'] = np.zeros((n_records),dtype=np.int32)
	#-- Retracked sigma 0 correction for Retracker 2: packed units (1e-2 dB)
	L2I_measurements['Retrack_2_sig0'] = np.zeros((n_records),dtype=np.int32)
	#-- Retracked sigma 0 correction for Retracker 3: packed units (1e-2 dB)
	L2I_measurements['Retrack_3_sig0'] = np.zeros((n_records),dtype=np.int32)
	#-- Retracker 1 quality metric
	L2I_measurements['Quality_1'] = np.zeros((n_records),dtype=np.int32)
	#-- Retracker 2 quality metric
	L2I_measurements['Quality_2'] = np.zeros((n_records),dtype=np.int32)
	#-- Retracker 3 quality metric
	L2I_measurements['Quality_3'] = np.zeros((n_records),dtype=np.int32)
	#-- Retrackers 3-23 output
	L2I_measurements['Retrack_3'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_4'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_5'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_6'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_7'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_8'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_9'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_10'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_11'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_12'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_13'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_14'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_15'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_16'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_17'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_18'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_19'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_20'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_21'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_22'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_23'] = np.zeros((n_records),dtype=np.int32)
	#-- Power echo shape parameter: packed units (dB/100)
	L2I_measurements['echo_shape'] = np.zeros((n_records),dtype=np.int32)
	#-- Beam behaviour parameter: unitless code number related to
	#-- surface characteristics
	L2I_measurements['BB_parameter'] = np.zeros((n_records,50),dtype=np.int16)
	#-- Cross track angle: packed units (micro radians)
	L2I_measurements['X_Track_Angle'] = np.zeros((n_records),dtype=np.int32)
	#-- Cross track angle correction: packed units (micro radians)
	L2I_measurements['X_Track_Angle_c'] = np.zeros((n_records),dtype=np.int32)
	#-- Leading edge coherence at retrack point 1/1000
	L2I_measurements['Coherence'] = np.zeros((n_records),dtype=np.int32)
	#-- Interpolated Ocean Height: packed units (mm above ellipsoid)
	L2I_measurements['Ocean_ht'] = np.zeros((n_records),dtype=np.int32)
	#-- Freeboard: packed units (mm, 1e-3 m)
	#-- -9999 default value indicates computation has not been performed
	L2I_measurements['Freeboard'] = np.zeros((n_records),dtype=np.int32)
	#-- Surface Height Anomaly: packed units (mm, 1e-3 m)
	L2I_measurements['SHA'] = np.zeros((n_records),dtype=np.int32)
	#-- Interpolated Surface Height Anomaly: packed units (mm, 1e-3 m)
	L2I_measurements['SHA_interp'] = np.zeros((n_records),dtype=np.int32)
	#-- Error in ocean height interpolation: packed units (mm, 1e-3 m)
	L2I_measurements['Interp_Err'] = np.zeros((n_records),dtype=np.uint16)
	#-- Number of forward records interpolated
	L2I_measurements['Interp_Count_Fwd'] = np.zeros((n_records),dtype=np.uint16)
	#-- Number of backward records interpolated
	L2I_measurements['Interp_Count_Bkwd'] = np.zeros((n_records),dtype=np.uint16)
	#-- Distance in time of most forward record interpolated (milli-seconds)
	L2I_measurements['Interp_Time_Fwd'] = np.zeros((n_records),dtype=np.uint16)
	#-- Distance in time of most backward record interpolated (milli-seconds)
	L2I_measurements['Interp_Time_Bkwd'] = np.zeros((n_records),dtype=np.uint16)
	#-- Interpolation error flag
	L2I_measurements['Interp_Error_Flg'] = np.zeros((n_records),dtype=np.uint16)
	#-- Measurement mode
	L2I_measurements['Meas_Mode'] = np.zeros((n_records),dtype=np.uint32)
	#-- Quality flags
	L2I_measurements['Quality_Flg'] = np.zeros((n_records),dtype=np.uint32)
	#-- Retracker flags
	L2I_measurements['Retracker_Flg'] = np.zeros((n_records),dtype=np.uint32)
	#-- Height calculation details
	#-- Specifies what was applied during the height calculation
	L2I_measurements['Height_status'] = np.zeros((n_records),dtype=np.uint32)
	#-- SAR freeboard status flag
	L2I_measurements['Freeboard_status'] = np.zeros((n_records),dtype=np.uint32)
	#-- Number of averaged echoes or beams
	L2I_measurements['N_avg'] = np.zeros((n_records),dtype=np.uint16)
	#-- Wind Speed packed units (mm/s, 1e-3 m/s)
	L2I_measurements['Wind_speed'] = np.zeros((n_records),dtype=np.uint16)
	L2I_measurements['Spares1'] = np.zeros((n_records,3),dtype=np.int32)

	#-- CryoSat-2 Auxiliary Data Group
	L2I_aux_data = {}
	#-- Ice Concentration packed units (%/1000)
	L2I_aux_data['Ice_conc'] = np.zeros((n_records),dtype=np.int32)
	#-- Snow Depth packed units (mm, 1e-3 m)
	L2I_aux_data['Snow_depth'] = np.zeros((n_records),dtype=np.int32)
	#-- Snow Density packed units (kg/m^3)
	L2I_aux_data['Snow_density'] = np.zeros((n_records),dtype=np.int32)
	#-- Discriminator result
	L2I_aux_data['Discriminator'] = np.zeros((n_records),dtype=np.int32)
	#-- SARin discriminator parameters 1-10
	L2I_aux_data['SARIN_disc_1'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['SARIN_disc_2'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['SARIN_disc_3'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['SARIN_disc_4'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['SARIN_disc_5'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['SARIN_disc_6'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['SARIN_disc_7'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['SARIN_disc_8'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['SARIN_disc_9'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['SARIN_disc_10'] = np.zeros((n_records),dtype=np.int32)
	#-- Discriminator flags
	L2I_aux_data['Discrim_Flg'] = np.zeros((n_records),dtype=np.uint32)
	#-- Slope model correction (Attitude of echo in micro-degrees)
	L2I_aux_data['Attitude'] = np.zeros((n_records),dtype=np.int32)
	#-- Slope model correction (Azimuth of echo in micro-degrees)
	L2I_aux_data['Azimuth'] = np.zeros((n_records),dtype=np.int32)
	#-- Slope doppler correction (mm)
	L2I_aux_data['Slope_doppler'] = np.zeros((n_records),dtype=np.int32)
	#-- The original latitude of the satellite (micro-degrees)
	L2I_aux_data['Lat_sat'] = np.zeros((n_records),dtype=np.int32)
	#-- The original longitude of the satellite (micro-degrees)
	L2I_aux_data['Lon_sat'] = np.zeros((n_records),dtype=np.int32)
	#-- Ambiguity indicator
	L2I_aux_data['Ambiguity'] = np.zeros((n_records),dtype=np.uint32)
	#-- Mean Sea Surface standard Model: packed units (mm, 1e-3 m)
	L2I_aux_data['MSS_model'] = np.zeros((n_records),dtype=np.int32)
	#-- Geoid standard Model: packed units (mm, 1e-3 m)
	L2I_aux_data['Geoid_model'] = np.zeros((n_records),dtype=np.int32)
	#-- ODLE standard Model: packed units (mm, 1e-3 m)
	L2I_aux_data['ODLE_model'] = np.zeros((n_records),dtype=np.int32)
	#-- The interpolated elevation value obtained from the DEM (mm)
	L2I_aux_data['DEM_elev'] = np.zeros((n_records),dtype=np.int32)
	#-- Identification of DEM used in SARin ambiguity test
	L2I_aux_data['DEM_ID'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['Spares2'] = np.zeros((n_records,4),dtype=np.int32)

	#-- CryoSat-2 External Corrections Group
	L2I_geo_corrections = {}
	#-- Dry Tropospheric Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['dryTrop'] = np.zeros((n_records),dtype=np.int32)
	#-- Wet Tropospheric Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['wetTrop'] = np.zeros((n_records),dtype=np.int32)
	#-- Inverse Barometric Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['InvBar'] = np.zeros((n_records),dtype=np.int32)
	#-- Delta Inverse Barometric Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['DAC'] = np.zeros((n_records),dtype=np.int32)
	#-- GIM Ionospheric Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['Iono_GIM'] = np.zeros((n_records),dtype=np.int32)
	#-- Model Ionospheric Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['Iono_model'] = np.zeros((n_records),dtype=np.int32)
	#-- Ocean tide Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['ocTideElv'] = np.zeros((n_records),dtype=np.int32)
	#-- Long period equilibrium ocean tide Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['lpeTideElv'] = np.zeros((n_records),dtype=np.int32)
	#-- Ocean loading tide Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['olTideElv'] = np.zeros((n_records),dtype=np.int32)
	#-- Solid Earth tide Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['seTideElv'] = np.zeros((n_records),dtype=np.int32)
	#-- Geocentric Polar tide Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['gpTideElv'] = np.zeros((n_records),dtype=np.int32)
	#-- Surface Type: Packed in groups of three bits for each of the 20 records
	L2I_geo_corrections['Surf_type'] = np.zeros((n_records),dtype=np.uint32)
	#-- Corrections Status Flag
	L2I_geo_corrections['Corr_status'] = np.zeros((n_records),dtype=np.uint32)
	#-- Correction Error Flag
	L2I_geo_corrections['Corr_error'] = np.zeros((n_records),dtype=np.uint32)
	#-- Sea State Bias Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['SSB'] = np.zeros((n_records),dtype=np.int32)
	L2I_geo_corrections['Spares3'] = np.zeros((n_records,2),dtype=np.int32)

	#-- CryoSat-2 Internal Corrections Group
	L2I_inst_corrections = {}
	#-- Doppler range correction: Radial + slope (mm)
	L2I_inst_corrections['Doppler_range'] = np.zeros((n_records),dtype=np.int32)
	#-- Instrument Range Correction: t-r antenna (mm)
	L2I_inst_corrections['TR_inst_range'] = np.zeros((n_records),dtype=np.int32)
	#-- Instrument Range Correction: r-only antenna (mm)
	L2I_inst_corrections['R_inst_range'] = np.zeros((n_records),dtype=np.int32)
	#-- Instrument Sigma 0 Correction: t-r antenna (dB/100)
	L2I_inst_corrections['TR_inst_gain'] = np.zeros((n_records),dtype=np.int32)
	#-- Instrument Sigma 0 Correction: r-only (dB/100)
	L2I_inst_corrections['R_inst_gain'] = np.zeros((n_records),dtype=np.int32)
	#-- Internal Phase Correction (milli-radians)
	L2I_inst_corrections['Internal_phase'] = np.zeros((n_records),dtype=np.int32)
	#-- External Phase Correction (milli-radians)
	L2I_inst_corrections['External_phase'] = np.zeros((n_records),dtype=np.int32)
	#-- Noise Power measurement
	L2I_inst_corrections['Noise_power'] = np.zeros((n_records),dtype=np.int32)
	#-- Phase slope correction (microradians)
	L2I_inst_corrections['Phase_slope'] = np.zeros((n_records),dtype=np.int32)
	L2I_inst_corrections['Spares4'] = np.zeros((n_records,2),dtype=np.int32)

	#-- for each record in the CryoSat file
	for r in range(n_records):
		#-- CryoSat-2 Location Group for record r
		L2I_location_parameters['Day'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_location_parameters['Sec'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_location_parameters['Micsec'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_location_parameters['USO_Corr'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_location_parameters['Mode_ID'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2I_location_parameters['SSC'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2I_location_parameters['Inst_config'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_location_parameters['Rec_Count'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_location_parameters['Lat'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_location_parameters['Lon'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_location_parameters['Alt'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_location_parameters['Alt_rate'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_location_parameters['Sat_velocity'][r,:] = np.fromfile(fid,dtype='>i4',count=3)
		L2I_location_parameters['Real_beam'][r,:] = np.fromfile(fid,dtype='>i4',count=3)
		L2I_location_parameters['Baseline'][r,:] = np.fromfile(fid,dtype='>i4',count=3)
		L2I_location_parameters['ST_ID'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2I_location_parameters['Spare'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2I_location_parameters['Roll'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_location_parameters['Pitch'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_location_parameters['Yaw'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_location_parameters['L2_MCD'][r] = np.fromfile(fid,dtype='>u4',count=1)

		#-- CryoSat-2 Measurement Group for record r
		L2I_measurements['Elev_1'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Elev_2'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Elev_3'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Sig0_1'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Sig0_2'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Sig0_3'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['SWH'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Peakiness'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_1_range'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_2_range'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_3_range'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_1_sig0'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_2_sig0'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_3_sig0'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Quality_1'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Quality_2'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Quality_3'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_3'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_4'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_5'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_6'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_7'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_8'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_9'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_10'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_11'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_12'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_13'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_14'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_15'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_16'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_17'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_18'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_19'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_20'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_21'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_22'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_23'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['echo_shape'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['BB_parameter'][r,:] = np.fromfile(fid,dtype='>i2',count=50)
		L2I_measurements['X_Track_Angle'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['X_Track_Angle_c'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Coherence'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Ocean_ht'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Freeboard'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['SHA'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['SHA_interp'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Interp_Err'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2I_measurements['Interp_Count_Fwd'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2I_measurements['Interp_Count_Bkwd'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2I_measurements['Interp_Time_Fwd'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2I_measurements['Interp_Time_Bkwd'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2I_measurements['Interp_Error_Flg'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2I_measurements['Meas_Mode'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_measurements['Quality_Flg'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_measurements['Retracker_Flg'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_measurements['Height_status'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_measurements['Freeboard_status'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_measurements['N_avg'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2I_measurements['Wind_speed'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2I_measurements['Spares1'][r,:] = np.fromfile(fid,dtype='>i4',count=3)

		#-- CryoSat-2 Auxiliary Data Group for record r
		L2I_aux_data['Ice_conc'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Snow_depth'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Snow_density'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Discriminator'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['SARIN_disc_1'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['SARIN_disc_2'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['SARIN_disc_3'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['SARIN_disc_4'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['SARIN_disc_5'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['SARIN_disc_6'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['SARIN_disc_7'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['SARIN_disc_8'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['SARIN_disc_9'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['SARIN_disc_10'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Discrim_Flg'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_aux_data['Attitude'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Azimuth'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Slope_doppler'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Lat_sat'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Lon_sat'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Ambiguity'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_aux_data['MSS_model'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Geoid_model'][r] =np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['ODLE_model'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['DEM_elev'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['DEM_ID'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Spares2'][r,:] = np.fromfile(fid,dtype='>i4',count=4)

		#-- CryoSat-2 External Corrections Group for record r
		L2I_geo_corrections['dryTrop'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['wetTrop'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['InvBar'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['DAC'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['Iono_GIM'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['Iono_model'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['ocTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['lpeTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['olTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['seTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['gpTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['Surf_type'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_geo_corrections['Corr_status'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_geo_corrections['Corr_error'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_geo_corrections['SSB'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['Spares3'][r,:] = np.fromfile(fid,dtype='>i4',count=2)

		#-- CryoSat-2 Internal Corrections Group for record r
		L2I_inst_corrections['Doppler_range'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_inst_corrections['TR_inst_range'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_inst_corrections['R_inst_range'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_inst_corrections['TR_inst_gain'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_inst_corrections['R_inst_gain'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_inst_corrections['Internal_phase'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_inst_corrections['External_phase'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_inst_corrections['Noise_power'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_inst_corrections['Phase_slope'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_inst_corrections['Spares4'][r,:] = np.fromfile(fid,dtype='>i4',count=2)

	#-- Bind all the bits of the l2i_mds together into a single dictionary
	CS_L2I_mds = {}
	CS_L2I_mds['Location'] = L2I_location_parameters
	CS_L2I_mds['Data'] = L2I_measurements
	CS_L2I_mds['Auxiliary'] = L2I_aux_data
	CS_L2I_mds['Geometry'] = L2I_geo_corrections
	CS_L2I_mds['Instrumental'] = L2I_inst_corrections
	#-- return the output dictionary
	return CS_L2I_mds

#-- PURPOSE: Initiate L2 MDS variables for CryoSat Baseline C
def cryosat_baseline_C(fid,record_size,n_records):
	#-- CryoSat-2 Location Group
	#-- Time and Orbit Parameters plus Measurement Mode
	L2I_location_parameters = {}
	#-- Time: day part
	L2I_location_parameters['Day'] = np.zeros((n_records),dtype=np.int32)
	#-- Time: second part
	L2I_location_parameters['Sec'] = np.zeros((n_records),dtype=np.uint32)
	#-- Time: microsecond part
	L2I_location_parameters['Micsec'] = np.zeros((n_records),dtype=np.uint32)
	#-- USO correction factor
	L2I_location_parameters['USO_Corr'] = np.zeros((n_records),dtype=np.int32)
	#-- Mode ID
	L2I_location_parameters['Mode_ID'] = np.zeros((n_records),dtype=np.uint16)
	#-- Source sequence counter
	L2I_location_parameters['SSC'] = np.zeros((n_records),dtype=np.uint16)
	#-- Instrument configuration
	L2I_location_parameters['Inst_config'] = np.zeros((n_records),dtype=np.uint32)
	#-- Record Counter
	L2I_location_parameters['Rec_Count'] = np.zeros((n_records),dtype=np.uint32)
	#-- Lat: packed units (0.1 micro-degree, 1e-7 degrees)
	L2I_location_parameters['Lat'] = np.zeros((n_records),dtype=np.int32)
	#-- Lon: packed units (0.1 micro-degree, 1e-7 degrees)
	L2I_location_parameters['Lon'] = np.zeros((n_records),dtype=np.int32)
	#-- Alt: packed units (mm, 1e-3 m)
	#-- Altitude of COG above reference ellipsoid (interpolated value)
	L2I_location_parameters['Alt'] = np.zeros((n_records),dtype=np.int32)
	#-- Instantaneous altitude rate derived from orbit: packed units (mm/s, 1e-3 m/s)
	L2I_location_parameters['Alt_rate'] = np.zeros((n_records),dtype=np.int32)
	#-- Satellite velocity vector. In ITRF: packed units (mm/s, 1e-3 m/s)
	L2I_location_parameters['Sat_velocity'] = np.zeros((n_records,3),dtype=np.int32)
	#-- Real beam direction vector. In CRF: packed units (micro-m/s, 1e-6 m/s)
	L2I_location_parameters['Real_beam'] = np.zeros((n_records,3),dtype=np.int32)
	#-- Interferometer baseline vector. In CRF: packed units (micro-m/s, 1e-6 m/s)
	L2I_location_parameters['Baseline'] = np.zeros((n_records,3),dtype=np.int32)
	#-- Star Tracker ID
	L2I_location_parameters['ST_ID'] = np.zeros((n_records),dtype=np.int16)
	L2I_location_parameters['Spare'] = np.zeros((n_records),dtype=np.int16)
	#-- Roll (Derived from star trackers): packed units (0.1 micro-degree, 1e-7 degrees)
	L2I_location_parameters['Roll'] = np.zeros((n_records),dtype=np.int32)
	#-- Pitch (Derived from star trackers): packed units (0.1 micro-degree, 1e-7 degrees)
	L2I_location_parameters['Pitch'] = np.zeros((n_records),dtype=np.int32)
	#-- Yaw (Derived from star trackers): packed units (0.1 micro-degree, 1e-7 degrees)
	L2I_location_parameters['Yaw'] = np.zeros((n_records),dtype=np.int32)
	#-- Measurement Confidence Data
	L2I_location_parameters['L2_MCD'] = np.zeros((n_records),dtype=np.uint32)

	#-- CryoSat-2 Measurement Group
	#-- Derived from instrument measurement parameters
	L2I_measurements = {}
	#-- Measured elevation above ellipsoid from retracker 1: packed units (mm, 1e-3 m)
	L2I_measurements['Elev_1'] = np.zeros((n_records),dtype=np.int32)
	#-- Measured elevation above ellipsoid from retracker 2: packed units (mm, 1e-3 m)
	L2I_measurements['Elev_2'] = np.zeros((n_records),dtype=np.int32)
	#-- Measured elevation above ellipsoid from retracker 3: packed units (mm, 1e-3 m)
	L2I_measurements['Elev_3'] = np.zeros((n_records),dtype=np.int32)
	#-- Sigma Zero Backscatter for retracker 1: packed units (1e-2 dB)
	L2I_measurements['Sig0_1'] = np.zeros((n_records),dtype=np.int32)
	#-- Sigma Zero Backscatter for retracker 2: packed units (1e-2 dB)
	L2I_measurements['Sig0_2'] = np.zeros((n_records),dtype=np.int32)
	#-- Sigma Zero Backscatter for retracker 3: packed units (1e-2 dB)
	L2I_measurements['Sig0_3'] = np.zeros((n_records),dtype=np.int32)
	#-- SWH packed units (mm, 1e-3)
	L2I_measurements['SWH'] = np.zeros((n_records),dtype=np.int32)
	#-- Peakiness: packed units (1e-2)
	L2I_measurements['Peakiness'] = np.zeros((n_records),dtype=np.int32)
	#-- Retracked range correction for retracker 1: packed units (mm, 1e-3 m)
	L2I_measurements['Retrack_1_range'] = np.zeros((n_records),dtype=np.int32)
	#-- Retracked range correction for retracker 1: packed units (mm, 1e-3 m)
	L2I_measurements['Retrack_2_range'] = np.zeros((n_records),dtype=np.int32)
	#-- Retracked range correction for retracker 3: packed units (mm, 1e-3 m)
	L2I_measurements['Retrack_3_range'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Spare2'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Spare3'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Spare4'] = np.zeros((n_records),dtype=np.int32)
	#-- Retracker 1 quality metric
	L2I_measurements['Quality_1'] = np.zeros((n_records),dtype=np.int32)
	#-- Retracker 2 quality metric
	L2I_measurements['Quality_2'] = np.zeros((n_records),dtype=np.int32)
	#-- Retracker 3 quality metric
	L2I_measurements['Quality_3'] = np.zeros((n_records),dtype=np.int32)
	#-- Retrackers 3-23 output
	L2I_measurements['Retrack_3'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_4'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_5'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_6'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_7'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_8'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_9'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_10'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_11'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_12'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_13'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_14'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_15'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_16'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_17'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_18'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_19'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_20'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_21'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_22'] = np.zeros((n_records),dtype=np.int32)
	L2I_measurements['Retrack_23'] = np.zeros((n_records),dtype=np.int32)
	#-- Power echo shape parameter: packed units (dB/100)
	L2I_measurements['echo_shape'] = np.zeros((n_records),dtype=np.int32)
	#-- Beam behaviour parameter: unitless code number related to
	#-- surface characteristics
	L2I_measurements['BB_parameter'] = np.zeros((n_records,50),dtype=np.int16)
	#-- Cross track angle: packed units (micro radians)
	L2I_measurements['X_Track_Angle'] = np.zeros((n_records),dtype=np.int32)
	#-- Cross track angle correction: packed units (micro radians)
	L2I_measurements['X_Track_Angle_c'] = np.zeros((n_records),dtype=np.int32)
	#-- Leading edge coherence at retrack point 1/1000
	L2I_measurements['Coherence'] = np.zeros((n_records),dtype=np.int32)
	#-- Interpolated Ocean Height: packed units (mm above ellipsoid)
	L2I_measurements['Ocean_ht'] = np.zeros((n_records),dtype=np.int32)
	#-- Freeboard: packed units (mm, 1e-3 m)
	#-- -9999 default value indicates computation has not been performed
	L2I_measurements['Freeboard'] = np.zeros((n_records),dtype=np.int32)
	#-- Surface Height Anomaly: packed units (mm, 1e-3 m)
	L2I_measurements['SHA'] = np.zeros((n_records),dtype=np.int32)
	#-- Interpolated Surface Height Anomaly: packed units (mm, 1e-3 m)
	L2I_measurements['SHA_interp'] = np.zeros((n_records),dtype=np.int32)
	#-- Error in ocean height interpolation: packed units (mm, 1e-3 m)
	L2I_measurements['Interp_Err'] = np.zeros((n_records),dtype=np.uint16)
	#-- Number of forward records interpolated
	L2I_measurements['Interp_Count_Fwd'] = np.zeros((n_records),dtype=np.uint16)
	#-- Number of backward records interpolated
	L2I_measurements['Interp_Count_Bkwd'] = np.zeros((n_records),dtype=np.uint16)
	#-- Distance in time of most forward record interpolated (milli-seconds)
	L2I_measurements['Interp_Time_Fwd'] = np.zeros((n_records),dtype=np.uint16)
	#-- Distance in time of most backward record interpolated (milli-seconds)
	L2I_measurements['Interp_Time_Bkwd'] = np.zeros((n_records),dtype=np.uint16)
	#-- Interpolation error flag
	L2I_measurements['Interp_Error_Flg'] = np.zeros((n_records),dtype=np.uint16)
	#-- Measurement mode
	L2I_measurements['Meas_Mode'] = np.zeros((n_records),dtype=np.uint32)
	#-- Quality flags
	L2I_measurements['Quality_Flg'] = np.zeros((n_records),dtype=np.uint32)
	#-- Retracker flags
	L2I_measurements['Retracker_Flg'] = np.zeros((n_records),dtype=np.uint32)
	#-- Height calculation details
	#-- Specifies what was applied during the height calculation
	L2I_measurements['Height_status'] = np.zeros((n_records),dtype=np.uint32)
	#-- SAR freeboard status flag
	L2I_measurements['Freeboard_status'] = np.zeros((n_records),dtype=np.uint32)
	#-- Number of averaged echoes or beams
	L2I_measurements['N_avg'] = np.zeros((n_records),dtype=np.uint16)
	#-- Wind Speed packed units (mm/s, 1e-3 m/s)
	L2I_measurements['Wind_speed'] = np.zeros((n_records),dtype=np.uint16)
	L2I_measurements['Spares1'] = np.zeros((n_records,3),dtype=np.int32)

	#-- CryoSat-2 Auxiliary Data Group
	L2I_aux_data = {}
	#-- Ice Concentration packed units (%/1000)
	L2I_aux_data['Ice_conc'] = np.zeros((n_records),dtype=np.int32)
	#-- Snow Depth packed units (mm, 1e-3 m)
	L2I_aux_data['Snow_depth'] = np.zeros((n_records),dtype=np.int32)
	#-- Snow Density packed units (kg/m^3)
	L2I_aux_data['Snow_density'] = np.zeros((n_records),dtype=np.int32)
	#-- Discriminator result
	L2I_aux_data['Discriminator'] = np.zeros((n_records),dtype=np.int32)
	#-- SARin discriminator parameters 1-10
	L2I_aux_data['SARIN_disc_1'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['SARIN_disc_2'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['SARIN_disc_3'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['SARIN_disc_4'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['SARIN_disc_5'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['SARIN_disc_6'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['SARIN_disc_7'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['SARIN_disc_8'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['SARIN_disc_9'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['SARIN_disc_10'] = np.zeros((n_records),dtype=np.int32)
	#-- Discriminator flags
	L2I_aux_data['Discrim_Flg'] = np.zeros((n_records),dtype=np.uint32)
	#-- Slope model correction (Attitude of echo in micro-degrees)
	L2I_aux_data['Attitude'] = np.zeros((n_records),dtype=np.int32)
	#-- Slope model correction (Azimuth of echo in micro-degrees)
	L2I_aux_data['Azimuth'] = np.zeros((n_records),dtype=np.int32)
	#-- Slope doppler correction (mm)
	L2I_aux_data['Slope_doppler'] = np.zeros((n_records),dtype=np.int32)
	#-- The original latitude of the satellite (micro-degrees)
	L2I_aux_data['Lat_sat'] = np.zeros((n_records),dtype=np.int32)
	#-- The original longitude of the satellite (micro-degrees)
	L2I_aux_data['Lon_sat'] = np.zeros((n_records),dtype=np.int32)
	#-- Ambiguity indicator
	L2I_aux_data['Ambiguity'] = np.zeros((n_records),dtype=np.uint32)
	#-- Mean Sea Surface Model packed units (mm, 1e-3 m)
	L2I_aux_data['MSS_model'] = np.zeros((n_records),dtype=np.int32)
	#-- Geoid Model packed units (mm, 1e-3 m)
	L2I_aux_data['Geoid_model'] = np.zeros((n_records),dtype=np.int32)
	#-- ODLE Model packed units (mm, 1e-3 m)
	L2I_aux_data['ODLE_model'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['DEM_elev'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['DEM_ID'] = np.zeros((n_records),dtype=np.int32)
	L2I_aux_data['Spares2'] = np.zeros((n_records,4),dtype=np.int32)

	#-- CryoSat-2 External Corrections Group
	L2I_geo_corrections = {}
	#-- Dry Tropospheric Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['dryTrop'] = np.zeros((n_records),dtype=np.int32)
	#-- Wet Tropospheric Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['wetTrop'] = np.zeros((n_records),dtype=np.int32)
	#-- Inverse Barometric Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['InvBar'] = np.zeros((n_records),dtype=np.int32)
	#-- Delta Inverse Barometric Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['DAC'] = np.zeros((n_records),dtype=np.int32)
	#-- GIM Ionospheric Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['Iono_GIM'] = np.zeros((n_records),dtype=np.int32)
	#-- Model Ionospheric Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['Iono_model'] = np.zeros((n_records),dtype=np.int32)
	#-- Ocean tide Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['ocTideElv'] = np.zeros((n_records),dtype=np.int32)
	#-- Long period equilibrium ocean tide Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['lpeTideElv'] = np.zeros((n_records),dtype=np.int32)
	#-- Ocean loading tide Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['olTideElv'] = np.zeros((n_records),dtype=np.int32)
	#-- Solid Earth tide Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['seTideElv'] = np.zeros((n_records),dtype=np.int32)
	#-- Geocentric Polar tide Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['gpTideElv'] = np.zeros((n_records),dtype=np.int32)
	#-- Surface Type: Packed in groups of three bits for each of the 20 records
	L2I_geo_corrections['Surf_type'] = np.zeros((n_records),dtype=np.uint32)
	#-- Corrections Status Flag
	L2I_geo_corrections['Corr_status'] = np.zeros((n_records),dtype=np.uint32)
	#-- Correction Error Flag
	L2I_geo_corrections['Corr_error'] = np.zeros((n_records),dtype=np.uint32)
	#-- Sea State Bias Correction packed units (mm, 1e-3 m)
	L2I_geo_corrections['SSB'] = np.zeros((n_records),dtype=np.int32)
	L2I_geo_corrections['Spares3'] = np.zeros((n_records,2),dtype=np.int32)

	#-- CryoSat-2 Internal Corrections Group
	L2I_inst_corrections = {}
	#-- Doppler range correction: Radial + slope (mm)
	L2I_inst_corrections['Doppler_range'] = np.zeros((n_records),dtype=np.int32)
	#-- Instrument Range Correction: t-r antenna (mm)
	L2I_inst_corrections['TR_inst_range'] = np.zeros((n_records),dtype=np.int32)
	#-- Instrument Range Correction: r-only antenna (mm)
	L2I_inst_corrections['R_inst_range'] = np.zeros((n_records),dtype=np.int32)
	#-- Instrument Sigma 0 Correction: t-r antenna (dB/100)
	L2I_inst_corrections['TR_inst_gain'] = np.zeros((n_records),dtype=np.int32)
	#-- Instrument Sigma 0 Correction: r-only (dB/100)
	L2I_inst_corrections['R_inst_gain'] = np.zeros((n_records),dtype=np.int32)
	#-- Internal Phase Correction (milli-radians)
	L2I_inst_corrections['Internal_phase'] = np.zeros((n_records),dtype=np.int32)
	#-- External Phase Correction (milli-radians)
	L2I_inst_corrections['External_phase'] = np.zeros((n_records),dtype=np.int32)
	#-- Noise Power measurement
	L2I_inst_corrections['Noise_power'] = np.zeros((n_records),dtype=np.int32)
	#-- Phase slope correction (microradians)
	L2I_inst_corrections['Phase_slope'] = np.zeros((n_records),dtype=np.int32)
	L2I_inst_corrections['Spares4'] = np.zeros((n_records,2),dtype=np.int32)

	#-- for each record in the CryoSat file
	for r in range(n_records):
		#-- CryoSat-2 Location Group for record r
		L2I_location_parameters['Day'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_location_parameters['Sec'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_location_parameters['Micsec'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_location_parameters['USO_Corr'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_location_parameters['Mode_ID'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2I_location_parameters['SSC'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2I_location_parameters['Inst_config'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_location_parameters['Rec_Count'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_location_parameters['Lat'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_location_parameters['Lon'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_location_parameters['Alt'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_location_parameters['Alt_rate'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_location_parameters['Sat_velocity'][r,:] = np.fromfile(fid,dtype='>i4',count=3)
		L2I_location_parameters['Real_beam'][r,:] = np.fromfile(fid,dtype='>i4',count=3)
		L2I_location_parameters['Baseline'][r,:] = np.fromfile(fid,dtype='>i4',count=3)
		L2I_location_parameters['ST_ID'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2I_location_parameters['Spare'][r] = np.fromfile(fid,dtype='>i2',count=1)
		L2I_location_parameters['Roll'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_location_parameters['Pitch'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_location_parameters['Yaw'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_location_parameters['L2_MCD'][r] = np.fromfile(fid,dtype='>u4',count=1)

		#-- CryoSat-2 Measurement Group for record r
		L2I_measurements['Elev_1'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Elev_2'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Elev_3'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Sig0_1'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Sig0_2'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Sig0_3'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['SWH'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Peakiness'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_1_range'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_2_range'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_3_range'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Spare2'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Spare3'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Spare4'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Quality_1'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Quality_2'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Quality_3'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_3'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_4'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_5'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_6'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_7'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_8'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_9'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_10'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_11'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_12'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_13'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_14'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_15'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_16'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_17'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_18'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_19'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_20'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_21'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_22'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Retrack_23'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['echo_shape'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['BB_parameter'][r,:] = np.fromfile(fid,dtype='>i2',count=50)
		L2I_measurements['X_Track_Angle'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['X_Track_Angle_c'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Coherence'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Ocean_ht'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Freeboard'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['SHA'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['SHA_interp'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_measurements['Interp_Err'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2I_measurements['Interp_Count_Fwd'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2I_measurements['Interp_Count_Bkwd'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2I_measurements['Interp_Time_Fwd'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2I_measurements['Interp_Time_Bkwd'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2I_measurements['Interp_Error_Flg'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2I_measurements['Meas_Mode'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_measurements['Quality_Flg'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_measurements['Retracker_Flg'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_measurements['Height_status'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_measurements['Freeboard_status'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_measurements['N_avg'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2I_measurements['Wind_speed'][r] = np.fromfile(fid,dtype='>u2',count=1)
		L2I_measurements['Spares1'][r,:] = np.fromfile(fid,dtype='>i4',count=3)

		#-- CryoSat-2 Auxiliary Data Group for record r
		L2I_aux_data['Ice_conc'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Snow_depth'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Snow_density'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Discriminator'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['SARIN_disc_1'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['SARIN_disc_2'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['SARIN_disc_3'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['SARIN_disc_4'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['SARIN_disc_5'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['SARIN_disc_6'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['SARIN_disc_7'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['SARIN_disc_8'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['SARIN_disc_9'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['SARIN_disc_10'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Discrim_Flg'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_aux_data['Attitude'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Azimuth'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Slope_doppler'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Lat_sat'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Lon_sat'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Ambiguity'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_aux_data['MSS_model'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Geoid_model'][r] =np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['ODLE_model'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['DEM_elev'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['DEM_ID'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_aux_data['Spares2'][r,:] = np.fromfile(fid,dtype='>i4',count=4)

		#-- CryoSat-2 External Corrections Group for record r
		L2I_geo_corrections['dryTrop'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['wetTrop'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['InvBar'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['DAC'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['Iono_GIM'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['Iono_model'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['ocTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['lpeTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['olTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['seTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['gpTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['Surf_type'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_geo_corrections['Corr_status'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_geo_corrections['Corr_error'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L2I_geo_corrections['SSB'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_geo_corrections['Spares3'][r,:] = np.fromfile(fid,dtype='>i4',count=2)

		#-- CryoSat-2 Internal Corrections Group for record r
		L2I_inst_corrections['Doppler_range'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_inst_corrections['TR_inst_range'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_inst_corrections['R_inst_range'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_inst_corrections['TR_inst_gain'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_inst_corrections['R_inst_gain'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_inst_corrections['Internal_phase'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_inst_corrections['External_phase'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_inst_corrections['Noise_power'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_inst_corrections['Phase_slope'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L2I_inst_corrections['Spares4'][r,:] = np.fromfile(fid,dtype='>i4',count=2)

	#-- Bind all the bits of the l2i_mds together into a single dictionary
	CS_L2I_mds = {}
	CS_L2I_mds['Location'] = L2I_location_parameters
	CS_L2I_mds['Data'] = L2I_measurements
	CS_L2I_mds['Auxiliary'] = L2I_aux_data
	CS_L2I_mds['Geometry'] = L2I_geo_corrections
	CS_L2I_mds['Instrumental'] = L2I_inst_corrections
	#-- return the output dictionary
	return CS_L2I_mds

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
		fid = open(full_filename, 'rb')
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
		fid = open(full_filename, 'rb')
		#-- iterate through CryoSat file and fill output variables
		CS_L2I_mds = read_cryosat_variables(fid,i_record_size,j_num_DSR)
		#-- close the input CryoSat binary file
		fid.close()

	#-- return the data and headers
	return CS_L2I_mds
