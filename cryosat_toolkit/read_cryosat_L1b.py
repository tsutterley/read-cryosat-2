#!/usr/bin/env python
u"""
read_cryosat_L1b.py
Written by Tyler Sutterley (03/2016)

Reads CryoSat Level-1b data products from baselines A, B and C
Supported CryoSat Modes: LRM, SAR, SARin, FDM, SID, GDR

INPUTS:
	full_filename: full path of CryoSat .DBL file

OUTPUTS:
	Location: Time and Orbit Group
	Data: Measurements Group
	Geometry: External Corrections Group
	Waveform_1Hz: Average Waveforms Group
	Waveform_20Hz: Waveforms Group (with SAR/SARIN Beam Behavior Parameters)
	METADATA: MPH, SPH and DSD Header data

UPDATE HISTORY:
Updated 05/2016: using __future__ print and division functions
Written 03/2016
"""
from __future__ import print_function
from __future__ import division

import os
import re
import numpy as np

#-- PURPOSE: Initiate L1b MDS variables for CryoSat Baselines A and B
def cryosat_baseline_AB(fid, n_records, MODE):
	n_SARIN_RW = 512
	n_SAR_RW = 128
	n_LRM_RW = 128
	n_blocks = 20
	n_BeamBehaviourParams = 50

	#-- CryoSat-2 Time and Orbit Group
	L1b_location_parameters = {}
	#-- Time: day part
	L1b_location_parameters['Day'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Time: second part
	L1b_location_parameters['Second'] = np.zeros((n_records,n_blocks),dtype=np.uint32)
	#-- Time: microsecond part
	L1b_location_parameters['Micsec'] = np.zeros((n_records,n_blocks),dtype=np.uint32)
	#-- USO correction factor
	L1b_location_parameters['USO_Corr'] = np.zeros((n_records,n_blocks),dtype=np.uint32)
	#-- Mode ID
	L1b_location_parameters['Mode_ID'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
	#-- Source sequence counter
	L1b_location_parameters['SSC'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
	#-- Instrument configuration
	L1b_location_parameters['Inst_config'] = np.zeros((n_records,n_blocks),dtype=np.uint32)
	#-- Record Counter
	L1b_location_parameters['Rec_Count'] = np.zeros((n_records,n_blocks),dtype=np.uint32)
	#-- Lat: packed units (0.1 micro-degree, 1e-7 degrees)
	L1b_location_parameters['Lat'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Lon: packed units (0.1 micro-degree, 1e-7 degrees)
	L1b_location_parameters['Lon'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Alt: packed units (mm, 1e-3 m)
	#-- Altitude of COG above reference ellipsoid (interpolated value)
	L1b_location_parameters['Alt'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Instantaneous altitude rate derived from orbit: packed units (mm/s, 1e-3 m/s)
	L1b_location_parameters['Alt_rate'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Satellite velocity vector. In ITRF: packed units (mm/s, 1e-3 m/s)
	#-- ITRF= International Terrestrial Reference Frame
	L1b_location_parameters['Sat_velocity'] = np.zeros((n_records,n_blocks,3),dtype=np.int32)
	#-- Real beam direction vector. In CRF: packed units (micro-m, 1e-6 m)
	#-- CRF= CryoSat Reference Frame.
	L1b_location_parameters['Real_beam'] = np.zeros((n_records,n_blocks,3),dtype=np.int32)
	#-- Interferometric baseline vector. In CRF: packed units (micro-m, 1e-6 m)
	L1b_location_parameters['Baseline'] = np.zeros((n_records,n_blocks,3),dtype=np.int32)
	#-- Measurement Confidence Data Flags
	#-- Generally the MCD flags indicate problems when set
	#-- If MCD is 0 then no problems or non-nominal conditions were detected
	#-- Serious errors are indicated by setting bit 31
	L1b_location_parameters['MCD'] = np.zeros((n_records,n_blocks),dtype=np.uint32)

	#-- CryoSat-2 Measurement Group
	#-- Derived from instrument measurement parameters
	L1b_measurements = {}
	#-- Window Delay reference (two-way) corrected for instrument delays
	L1b_measurements['TD'] = np.zeros((n_records,n_blocks),dtype=np.int64)
	#-- H0 Initial Height Word from telemetry
	L1b_measurements['H_0'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- COR2 Height Rate: on-board tracker height rate over the radar cycle
	L1b_measurements['COR2'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Coarse Range Word (LAI) derived from telemetry
	L1b_measurements['LAI'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Fine Range Word (FAI) derived from telemetry
	L1b_measurements['FAI'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Automatic Gain Control Channel 1: AGC gain applied on Rx channel 1.
	#-- Gain calibration corrections are applied (Sum of AGC stages 1 and 2
	#-- plus the corresponding corrections) (dB/100)
	L1b_measurements['AGC_CH1'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Automatic Gain Control Channel 2: AGC gain applied on Rx channel 2.
	#-- Gain calibration corrections are applied (dB/100)
	L1b_measurements['AGC_CH2'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Total Fixed Gain On Channel 1: gain applied by the RF unit. (dB/100)
	L1b_measurements['TR_gain_CH1'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Total Fixed Gain On Channel 2: gain applied by the RF unit. (dB/100)
	L1b_measurements['TR_gain_CH2'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Transmit Power in microWatts
	L1b_measurements['TX_Power'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Doppler range correction: Radial component (mm)
	#-- computed for the component of satellite velocity in the nadir direction
	L1b_measurements['Doppler_range'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Instrument Range Correction: transmit-receive antenna (mm)
	#-- Calibration correction to range on channel 1 computed from CAL1.
	L1b_measurements['TR_inst_range'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Instrument Range Correction: receive-only antenna (mm)
	#-- Calibration correction to range on channel 2 computed from CAL1.
	L1b_measurements['R_inst_range'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Instrument Gain Correction: transmit-receive antenna (dB/100)
	#-- Calibration correction to gain on channel 1 computed from CAL1
	L1b_measurements['TR_inst_gain'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Instrument Gain Correction: receive-only (dB/100)
	#-- Calibration correction to gain on channel 2 computed from CAL1
	L1b_measurements['R_inst_gain'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Internal Phase Correction (microradians)
	L1b_measurements['Internal_phase'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- External Phase Correction (microradians)
	L1b_measurements['External_phase'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Noise Power measurement (dB/100): converted from telemetry units to be
	#-- the noise floor of FBR measurement echoes.
	#-- Set to -9999.99 when the telemetry contains zero.
	L1b_measurements['Noise_power'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Phase slope correction (microradians)
	#-- Computed from the CAL-4 packets during the azimuth impulse response
	#-- amplitude (SARIN only). Set from the latest available CAL-4 packet.
	L1b_measurements['Phase_slope'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	L1b_measurements['Spares1'] = np.zeros((n_records,n_blocks,4),dtype=np.int8)

	#-- CryoSat-2 External Corrections Group
	L1b_geo_corrections = {}
	#-- Dry Tropospheric Correction packed units (mm, 1e-3 m)
	L1b_geo_corrections['dryTrop'] = np.zeros((n_records),dtype=np.int32)
	#-- Wet Tropospheric Correction packed units (mm, 1e-3 m)
	L1b_geo_corrections['wetTrop'] = np.zeros((n_records),dtype=np.int32)
	#-- Inverse Barometric Correction packed units (mm, 1e-3 m)
	L1b_geo_corrections['InvBar'] = np.zeros((n_records),dtype=np.int32)
	#-- Delta Inverse Barometric Correction packed units (mm, 1e-3 m)
	L1b_geo_corrections['DAC'] = np.zeros((n_records),dtype=np.int32)
	#-- GIM Ionospheric Correction packed units (mm, 1e-3 m)
	L1b_geo_corrections['Iono_GIM'] = np.zeros((n_records),dtype=np.int32)
	#-- Model Ionospheric Correction packed units (mm, 1e-3 m)
	L1b_geo_corrections['Iono_model'] = np.zeros((n_records),dtype=np.int32)
	#-- Ocean tide Correction packed units (mm, 1e-3 m)
	L1b_geo_corrections['ocTideElv'] = np.zeros((n_records),dtype=np.int32)
	#-- Long period equilibrium ocean tide Correction packed units (mm, 1e-3 m)
	L1b_geo_corrections['lpeTideElv'] = np.zeros((n_records),dtype=np.int32)
	#-- Ocean loading tide Correction packed units (mm, 1e-3 m)
	L1b_geo_corrections['olTideElv'] = np.zeros((n_records),dtype=np.int32)
	#-- Solid Earth tide Correction packed units (mm, 1e-3 m)
	L1b_geo_corrections['seTideElv'] = np.zeros((n_records),dtype=np.int32)
	#-- Geocentric Polar tide Correction packed units (mm, 1e-3 m)
	L1b_geo_corrections['gpTideElv'] = np.zeros((n_records),dtype=np.int32)
	#-- Surface Type: enumerated key to classify surface at nadir
	#-- 0 = Open Ocean
	#-- 1 = Closed Sea
	#-- 2 = Continental Ice
	#-- 3 = Land
	L1b_geo_corrections['Surf_type'] = np.zeros((n_records),dtype=np.uint32)
	L1b_geo_corrections['Spare1'] = np.zeros((n_records,4),dtype=np.int8)
	#-- Corrections Status Flag
	L1b_geo_corrections['Corr_status'] = np.zeros((n_records),dtype=np.uint32)
	#-- Correction Error Flag
	L1b_geo_corrections['Corr_error'] = np.zeros((n_records),dtype=np.uint32)
	L1b_geo_corrections['Spare2'] = np.zeros((n_records,4),dtype=np.int8)

	#-- CryoSat-2 Average Waveforms Groups
	#-- Low-Resolution Mode
	L1b_1Hz_LRM_waveform = {}
	#-- Data Record Time (MDSR Time Stamp)
	L1b_1Hz_LRM_waveform['Day_1Hz'] = np.zeros((n_records),dtype=np.int32)
	L1b_1Hz_LRM_waveform['Sec_1Hz'] = np.zeros((n_records),dtype=np.uint32)
	L1b_1Hz_LRM_waveform['Micsec_1Hz'] = np.zeros((n_records),dtype=np.uint32)
	#-- Lat: packed units (0.1 micro-degree, 1e-7 degrees)
	L1b_1Hz_LRM_waveform['Lat_1Hz'] = np.zeros((n_records),dtype=np.int32)
	#-- Lon: packed units (0.1 micro-degree, 1e-7 degrees)
	L1b_1Hz_LRM_waveform['Lon_1Hz'] = np.zeros((n_records),dtype=np.int32)
	#-- Alt: packed units (mm, 1e-3 m)
	#-- Altitude of COG above reference ellipsoid (interpolated value)
	L1b_1Hz_LRM_waveform['Alt_1Hz'] = np.zeros((n_records),dtype=np.int32)
	#-- Window Delay (two-way) corrected for instrument delays
	L1b_1Hz_LRM_waveform['TD_1Hz'] = np.zeros((n_records),dtype=np.int64)
	#-- 1 Hz Averaged Power Echo Waveform
	L1b_1Hz_LRM_waveform['Waveform'] = np.zeros((n_records,n_LRM_RW),dtype=np.uint16)
	#-- Echo Scale Factor (to scale echo to watts)
	L1b_1Hz_LRM_waveform['Linear_Wfm_Multiplier'] = np.zeros((n_records),dtype=np.int32)
	#-- Echo Scale Power (a power of 2)
	L1b_1Hz_LRM_waveform['Power2_Wfm_Multiplier'] = np.zeros((n_records),dtype=np.int32)
	#-- Number of echoes averaged
	L1b_1Hz_LRM_waveform['N_avg_echoes'] = np.zeros((n_records),dtype=np.uint16)
	L1b_1Hz_LRM_waveform['Flags'] = np.zeros((n_records),dtype=np.uint16)

	#-- SAR Mode
	L1b_1Hz_SAR_waveform = {}
	#-- Data Record Time (MDSR Time Stamp)
	L1b_1Hz_SAR_waveform['Day_1Hz'] = np.zeros((n_records),dtype=np.int32)
	L1b_1Hz_SAR_waveform['Sec_1Hz'] = np.zeros((n_records),dtype=np.uint32)
	L1b_1Hz_SAR_waveform['Micsec_1Hz'] = np.zeros((n_records),dtype=np.uint32)
	#-- Lat: packed units (0.1 micro-degree, 1e-7 degrees)
	L1b_1Hz_SAR_waveform['Lat_1Hz'] = np.zeros((n_records),dtype=np.int32)
	#-- Lon: packed units (0.1 micro-degree, 1e-7 degrees)
	L1b_1Hz_SAR_waveform['Lon_1Hz'] = np.zeros((n_records),dtype=np.int32)
	#-- Alt: packed units (mm, 1e-3 m)
	#-- Altitude of COG above reference ellipsoid (interpolated value)
	L1b_1Hz_SAR_waveform['Alt_1Hz'] = np.zeros((n_records),dtype=np.int32)
	#-- Window Delay (two-way) corrected for instrument delays
	L1b_1Hz_SAR_waveform['TD_1Hz'] = np.zeros((n_records),dtype=np.int64)
	#-- 1 Hz Averaged Power Echo Waveform
	L1b_1Hz_SAR_waveform['Waveform'] = np.zeros((n_records,n_SAR_RW),dtype=np.uint16)
	#-- Echo Scale Factor (to scale echo to watts)
	L1b_1Hz_SAR_waveform['Linear_Wfm_Multiplier'] = np.zeros((n_records),dtype=np.int32)
	#-- Echo Scale Power (a power of 2)
	L1b_1Hz_SAR_waveform['Power2_Wfm_Multiplier'] = np.zeros((n_records),dtype=np.int32)
	#-- Number of echoes averaged
	L1b_1Hz_SAR_waveform['N_avg_echoes'] = np.zeros((n_records),dtype=np.uint16)
	L1b_1Hz_SAR_waveform['Flags'] = np.zeros((n_records),dtype=np.uint16)

	#-- SARIN Mode
	#-- Same as the LRM/SAR groups but the waveform array is 512 bins instead of
	#-- 128 and the number of echoes averaged is different.
	L1b_1Hz_SARIN_waveform = {}
	#-- Data Record Time (MDSR Time Stamp)
	L1b_1Hz_SARIN_waveform['Day'] = np.zeros((n_records),dtype=np.int32)
	L1b_1Hz_SARIN_waveform['Sec'] = np.zeros((n_records),dtype=np.uint32)
	L1b_1Hz_SARIN_waveform['Micsec'] = np.zeros((n_records),dtype=np.uint32)
	#-- Lat: packed units (0.1 micro-degree, 1e-7 degrees)
	L1b_1Hz_SARIN_waveform['Lat'] = np.zeros((n_records),dtype=np.int32)
	#-- Lon: packed units (0.1 micro-degree, 1e-7 degrees)
	L1b_1Hz_SARIN_waveform['Lon'] = np.zeros((n_records),dtype=np.int32)
	#-- Alt: packed units (mm, 1e-3 m)
	#-- Altitude of COG above reference ellipsoid (interpolated value)
	L1b_1Hz_SARIN_waveform['Alt'] = np.zeros((n_records),dtype=np.int32)
	#-- Window Delay (two-way) corrected for instrument delays
	L1b_1Hz_SARIN_waveform['TD'] = np.zeros((n_records),dtype=np.int64)
	#-- 1 Hz Averaged Power Echo Waveform
	L1b_1Hz_SARIN_waveform['Waveform'] = np.zeros((n_records,n_SARIN_RW),dtype=np.uint16)
	#-- Echo Scale Factor (to scale echo to watts)
	L1b_1Hz_SARIN_waveform['Linear_Wfm_Multiplier'] = np.zeros((n_records),dtype=np.int32)
	#-- Echo Scale Power (a power of 2)
	L1b_1Hz_SARIN_waveform['Power2_Wfm_Multiplier'] = np.zeros((n_records),dtype=np.int32)
	#-- Number of echoes averaged
	L1b_1Hz_SARIN_waveform['N_avg_echoes'] = np.zeros((n_records),dtype=np.uint16)
	L1b_1Hz_SARIN_waveform['Flags'] = np.zeros((n_records),dtype=np.uint16)

	#-- CryoSat-2 Waveforms Groups
	#-- Beam Behavior Parameters
	L1b_Beam_Behavior = {}
	#-- Standard Deviation of Gaussian fit to range integrated stack power.
	L1b_Beam_Behavior['SD'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
	#-- Stack Center: Mean of Gaussian fit to range integrated stack power.
	L1b_Beam_Behavior['Center'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
	#-- Stack amplitude parameter scaled in dB/100.
	L1b_Beam_Behavior['Amplitude'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
	#-- 3rd moment: providing the degree of asymmetry of the range integrated
	#-- stack power distribution.
	L1b_Beam_Behavior['Skewness'] = np.zeros((n_records,n_blocks),dtype=np.int16)
	#-- 4th moment: Measure of peakiness of range integrated stack power distribution.
	L1b_Beam_Behavior['Kurtosis'] = np.zeros((n_records,n_blocks),dtype=np.int16)
	L1b_Beam_Behavior['Spare'] = np.zeros((n_records,n_blocks,n_BeamBehaviourParams-5),dtype=np.int16)

	#-- Low-Resolution Mode
	L1b_LRM_waveform = {}
	#-- Averaged Power Echo Waveform [128]
	L1b_LRM_waveform['Waveform'] = np.zeros((n_records,n_blocks,n_LRM_RW),dtype=np.uint16)
	#-- Echo Scale Factor (to scale echo to watts)
	L1b_LRM_waveform['Linear_Wfm_Multiplier'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Echo Scale Power (a power of 2 to scale echo to Watts)
	L1b_LRM_waveform['Power2_Wfm_Multiplier'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Number of echoes averaged
	L1b_LRM_waveform['N_avg_echoes'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
	L1b_LRM_waveform['Flags'] = np.zeros((n_records,n_blocks),dtype=np.uint16)

	#-- SAR Mode
	L1b_SAR_waveform = {}
	#-- Averaged Power Echo Waveform [128]
	L1b_SAR_waveform['Waveform'] = np.zeros((n_records,n_blocks,n_SAR_RW),dtype=np.uint16)
	#-- Echo Scale Factor (to scale echo to watts)
	L1b_SAR_waveform['Linear_Wfm_Multiplier'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Echo Scale Power (a power of 2 to scale echo to Watts)
	L1b_SAR_waveform['Power2_Wfm_Multiplier'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Number of echoes averaged
	L1b_SAR_waveform['N_avg_echoes'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
	L1b_SAR_waveform['Flags'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
	#-- Beam behaviour parameters
	L1b_SAR_waveform['Beam'] = L1b_Beam_Behavior

	#-- SARIN Mode
	L1b_SARIN_waveform = {}
	#-- Averaged Power Echo Waveform [512]
	L1b_SARIN_waveform['Waveform'] = np.zeros((n_records,n_blocks,n_SARIN_RW),dtype=np.uint16)
	#-- Echo Scale Factor (to scale echo to watts)
	L1b_SARIN_waveform['Linear_Wfm_Multiplier'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Echo Scale Power (a power of 2 to scale echo to Watts)
	L1b_SARIN_waveform['Power2_Wfm_Multiplier'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Number of echoes averaged
	L1b_SARIN_waveform['N_avg_echoes'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
	L1b_SARIN_waveform['Flags'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
	#-- Beam behaviour parameters
	L1b_SARIN_waveform['Beam'] = L1b_Beam_Behavior
	#-- Coherence [512]: packed units (1/1000)
	L1b_SARIN_waveform['Coherence'] = np.zeros((n_records,n_blocks,n_SARIN_RW),dtype=np.int16)
	#-- Phase Difference [512]: packed units (microradians)
	L1b_SARIN_waveform['Phase_diff'] = np.zeros((n_records,n_blocks,n_SARIN_RW),dtype=np.int32)

	#-- for each record in the CryoSat file
	for r in range(n_records):
		#-- CryoSat-2 Time and Orbit Group
		for b in range(n_blocks):
			L1b_location_parameters['Day'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_location_parameters['Second'][r,b] = np.fromfile(fid,dtype='>u4',count=1)
			L1b_location_parameters['Micsec'][r,b] = np.fromfile(fid,dtype='>u4',count=1)
			L1b_location_parameters['USO_Corr'][r,b] = np.fromfile(fid,dtype='>u4',count=1)
			L1b_location_parameters['Mode_ID'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
			L1b_location_parameters['SSC'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
			L1b_location_parameters['Inst_config'][r,b] = np.fromfile(fid,dtype='>u4',count=1)
			L1b_location_parameters['Rec_Count'][r,b] = np.fromfile(fid,dtype='>u4',count=1)
			L1b_location_parameters['Lat'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_location_parameters['Lon'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_location_parameters['Alt'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_location_parameters['Alt_rate'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_location_parameters['Sat_velocity'][r,b,:] = np.fromfile(fid,dtype='>i4',count=3)
			L1b_location_parameters['Real_beam'][r,b,:] = np.fromfile(fid,dtype='>i4',count=3)
			L1b_location_parameters['Baseline'][r,b,:] = np.fromfile(fid,dtype='>i4',count=3)
			L1b_location_parameters['MCD'][r,b] = np.fromfile(fid,dtype='>u4',count=1)

		#-- CryoSat-2 Measurement Group
		#-- Derived from instrument measurement parameters
		for b in range(n_blocks):
			L1b_measurements['TD'][r,b] = np.fromfile(fid,dtype='>i8',count=1)
			L1b_measurements['H_0'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_measurements['COR2'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_measurements['LAI'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_measurements['FAI'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_measurements['AGC_CH1'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_measurements['AGC_CH2'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_measurements['TR_gain_CH1'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_measurements['TR_gain_CH2'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_measurements['TX_Power'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_measurements['Doppler_range'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_measurements['TR_inst_range'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_measurements['R_inst_range'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_measurements['TR_inst_gain'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_measurements['R_inst_gain'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_measurements['Internal_phase'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_measurements['External_phase'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_measurements['Noise_power'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_measurements['Phase_slope'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_measurements['Spares1'][r,b,:] = np.fromfile(fid,dtype='>i1',count=4)

		#-- CryoSat-2 External Corrections Group
		L1b_geo_corrections['dryTrop'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L1b_geo_corrections['wetTrop'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L1b_geo_corrections['InvBar'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L1b_geo_corrections['DAC'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L1b_geo_corrections['Iono_GIM'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L1b_geo_corrections['Iono_model'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L1b_geo_corrections['ocTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L1b_geo_corrections['lpeTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L1b_geo_corrections['olTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L1b_geo_corrections['seTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L1b_geo_corrections['gpTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L1b_geo_corrections['Surf_type'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L1b_geo_corrections['Spare1'][r,:] = np.fromfile(fid,dtype='>i1',count=4)
		L1b_geo_corrections['Corr_status'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L1b_geo_corrections['Corr_error'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L1b_geo_corrections['Spare2'][r,:] = np.fromfile(fid,dtype='>i1',count=4)

		#-- CryoSat-2 Average Waveforms Groups
		if (MODE == 'LRM'):
			#-- Low-Resolution Mode
			L1b_1Hz_LRM_waveform['Day_1Hz'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_1Hz_LRM_waveform['Sec_1Hz'][r] = np.fromfile(fid,dtype='>u4',count=1)
			L1b_1Hz_LRM_waveform['Micsec_1Hz'][r] = np.fromfile(fid,dtype='>u4',count=1)
			L1b_1Hz_LRM_waveform['Lat_1Hz'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_1Hz_LRM_waveform['Lon_1Hz'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_1Hz_LRM_waveform['Alt_1Hz'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_1Hz_LRM_waveform['TD_1Hz'][r] = np.fromfile(fid,dtype='>i8',count=1)
			L1b_1Hz_LRM_waveform['Waveform'][r,:] = np.fromfile(fid,dtype='>u2',count=n_LRM_RW)
			L1b_1Hz_LRM_waveform['Linear_Wfm_Multiplier'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_1Hz_LRM_waveform['Power2_Wfm_Multiplier'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_1Hz_LRM_waveform['N_avg_echoes'][r] = np.fromfile(fid,dtype='>u2',count=1)
			L1b_1Hz_LRM_waveform['Flags'][r] = np.fromfile(fid,dtype='>u2',count=1)
		elif (MODE == 'SAR'):
			#-- SAR Mode
			L1b_1Hz_SAR_waveform['Day_1Hz'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_1Hz_SAR_waveform['Sec_1Hz'][r] = np.fromfile(fid,dtype='>u4',count=1)
			L1b_1Hz_SAR_waveform['Micsec_1Hz'][r] = np.fromfile(fid,dtype='>u4',count=1)
			L1b_1Hz_SAR_waveform['Lat_1Hz'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_1Hz_SAR_waveform['Lon_1Hz'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_1Hz_SAR_waveform['Alt_1Hz'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_1Hz_SAR_waveform['TD_1Hz'][r] = np.fromfile(fid,dtype='>i8',count=1)
			L1b_1Hz_SAR_waveform['Waveform'][r,:] = np.fromfile(fid,dtype='>u2',count=n_SAR_RW)
			L1b_1Hz_SAR_waveform['Linear_Wfm_Multiplier'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_1Hz_SAR_waveform['Power2_Wfm_Multiplier'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_1Hz_SAR_waveform['N_avg_echoes'][r] = np.fromfile(fid,dtype='>u2',count=1)
			L1b_1Hz_SAR_waveform['Flags'][r] = np.fromfile(fid,dtype='>u2',count=1)
		elif (MODE == 'SIN'):
			#-- SARIN Mode
			L1b_1Hz_SARIN_waveform['Day'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_1Hz_SARIN_waveform['Sec'][r] = np.fromfile(fid,dtype='>u4',count=1)
			L1b_1Hz_SARIN_waveform['Micsec'][r] = np.fromfile(fid,dtype='>u4',count=1)
			L1b_1Hz_SARIN_waveform['Lat'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_1Hz_SARIN_waveform['Lon'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_1Hz_SARIN_waveform['Alt'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_1Hz_SARIN_waveform['TD'][r] = np.fromfile(fid,dtype='>i8',count=1)
			L1b_1Hz_SARIN_waveform['Waveform'][r,:] = np.fromfile(fid,dtype='>u2',count=n_SARIN_RW)
			L1b_1Hz_SARIN_waveform['Linear_Wfm_Multiplier'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_1Hz_SARIN_waveform['Power2_Wfm_Multiplier'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_1Hz_SARIN_waveform['N_avg_echoes'][r] = np.fromfile(fid,dtype='>u2',count=1)
			L1b_1Hz_SARIN_waveform['Flags'][r] = np.fromfile(fid,dtype='>u2',count=1)

		#-- CryoSat-2 Waveforms Groups
		if (MODE == 'LRM'):
			#-- Low-Resolution Mode
			for b in range(n_blocks):
				L1b_LRM_waveform['Waveform'][r,b,:] = np.fromfile(fid,dtype='>u2',count=n_LRM_RW)
				L1b_LRM_waveform['Linear_Wfm_Multiplier'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
				L1b_LRM_waveform['Power2_Wfm_Multiplier'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
				L1b_LRM_waveform['N_avg_echoes'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
				L1b_LRM_waveform['Flags'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
		elif (MODE == 'SAR'):
			#-- SAR Mode
			for b in range(n_blocks):
				L1b_SAR_waveform['Waveform'][r,b,:] = np.fromfile(fid,dtype='>u2',count=n_SAR_RW)
				L1b_SAR_waveform['Linear_Wfm_Multiplier'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
				L1b_SAR_waveform['Power2_Wfm_Multiplier'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
				L1b_SAR_waveform['N_avg_echoes'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
				L1b_SAR_waveform['Flags'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
				L1b_SAR_waveform['Beam']['SD'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
				L1b_SAR_waveform['Beam']['Center'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
				L1b_SAR_waveform['Beam']['Amplitude'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
				L1b_SAR_waveform['Beam']['Skewness'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
				L1b_SAR_waveform['Beam']['Kurtosis'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
				L1b_SAR_waveform['Beam']['Spare'][r,b,:] = np.fromfile(fid,dtype='>i2',count=(n_BeamBehaviourParams-5))
		elif (MODE == 'SIN'):
			#-- SARIN Mode
			for b in range(n_blocks):
				L1b_SARIN_waveform['Waveform'][r,b,:] = np.fromfile(fid,dtype='>u2',count=n_SARIN_RW)
				L1b_SARIN_waveform['Linear_Wfm_Multiplier'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
				L1b_SARIN_waveform['Power2_Wfm_Multiplier'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
				L1b_SARIN_waveform['N_avg_echoes'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
				L1b_SARIN_waveform['Flags'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
				L1b_SARIN_waveform['Beam']['SD'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
				L1b_SARIN_waveform['Beam']['Center'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
				L1b_SARIN_waveform['Beam']['Amplitude'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
				L1b_SARIN_waveform['Beam']['Skewness'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
				L1b_SARIN_waveform['Beam']['Kurtosis'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
				L1b_SARIN_waveform['Beam']['Spare'][r,b,:] = np.fromfile(fid,dtype='>i2',count=(n_BeamBehaviourParams-5))
				L1b_SARIN_waveform['Coherence'][r,b,:] = np.fromfile(fid,dtype='>i2',count=n_SARIN_RW)
				L1b_SARIN_waveform['Phase_diff'][r,b,:] = np.fromfile(fid,dtype='>i4',count=n_SARIN_RW)

	#-- Bind all the bits of the l1b_mds together into a single dictionary
	CS_l1b_mds = {}
	CS_l1b_mds['Location'] = L1b_location_parameters
	CS_l1b_mds['Data'] = L1b_measurements
	CS_l1b_mds['Geometry'] = L1b_geo_corrections
	if (MODE == 'LRM'):
		CS_l1b_mds['Waveform_1Hz'] = L1b_1Hz_LRM_waveform
		CS_l1b_mds['Waveform_20Hz'] = L1b_LRM_waveform
	elif (MODE == 'SAR'):
		CS_l1b_mds['Waveform_1Hz'] = L1b_1Hz_SAR_waveform
		CS_l1b_mds['Waveform_20Hz'] = L1b_SAR_waveform
	elif (MODE == 'SIN'):
		CS_l1b_mds['Waveform_1Hz'] = L1b_1Hz_SARIN_waveform
		CS_l1b_mds['Waveform_20Hz'] = L1b_SARIN_waveform

	#-- return the output dictionary
	return CS_l1b_mds

#-- PURPOSE: Initiate L1b MDS variables for CryoSat Baseline C
def cryosat_baseline_C(fid, n_records, MODE):
	n_SARIN_BC_RW = 1024
	n_SARIN_RW = 512
	n_SAR_BC_RW = 256
	n_SAR_RW = 128
	n_LRM_RW = 128
	n_blocks = 20
	n_BeamBehaviourParams = 50

	#-- CryoSat-2 Time and Orbit Group
	L1b_c_location_parameters = {}
	#-- Time: day part
	L1b_c_location_parameters['Day'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Time: second part
	L1b_c_location_parameters['Second'] = np.zeros((n_records,n_blocks),dtype=np.uint32)
	#-- Time: microsecond part
	L1b_c_location_parameters['Micsec'] = np.zeros((n_records,n_blocks),dtype=np.uint32)
	#-- USO correction factor
	L1b_c_location_parameters['USO_Corr'] = np.zeros((n_records,n_blocks),dtype=np.uint32)
	#-- Mode ID
	L1b_c_location_parameters['Mode_ID'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
	#-- Source sequence counter
	L1b_c_location_parameters['SSC'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
	#-- Instrument configuration
	L1b_c_location_parameters['Inst_config'] = np.zeros((n_records,n_blocks),dtype=np.uint32)
	#-- Record Counter
	L1b_c_location_parameters['Rec_Count'] = np.zeros((n_records,n_blocks),dtype=np.uint32)
	#-- Lat: packed units (0.1 micro-degree, 1e-7 degrees)
	L1b_c_location_parameters['Lat'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Lon: packed units (0.1 micro-degree, 1e-7 degrees)
	L1b_c_location_parameters['Lon'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Alt: packed units (mm, 1e-3 m)
	#-- Altitude of COG above reference ellipsoid (interpolated value)
	L1b_c_location_parameters['Alt'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Instantaneous altitude rate derived from orbit: packed units (mm/s, 1e-3 m/s)
	L1b_c_location_parameters['Alt_rate'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Satellite velocity vector. In ITRF: packed units (mm/s, 1e-3 m/s)
	#-- ITRF= International Terrestrial Reference Frame
	L1b_c_location_parameters['Sat_velocity'] = np.zeros((n_records,n_blocks,3),dtype=np.int32)
	#-- Real beam direction vector. In CRF: packed units (micro-m/s, 1e-6 m/s)
	#-- CRF= CryoSat Reference Frame.
	L1b_c_location_parameters['Real_beam'] = np.zeros((n_records,n_blocks,3),dtype=np.int32)
	#-- Interferometric baseline vector. In CRF: packed units (micro-m/s, 1e-6 m/s)
	L1b_c_location_parameters['Baseline'] = np.zeros((n_records,n_blocks,3),dtype=np.int32)
	#-- Star Tracker ID
	L1b_c_location_parameters['ST_ID'] = np.zeros((n_records,n_blocks),dtype=np.int16)
	#-- Antenna Bench Roll Angle (Derived from star trackers)
	#-- packed units (0.1 micro-degree, 1e-7 degrees)
	L1b_c_location_parameters['Roll'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Antenna Bench Pitch Angle (Derived from star trackers)
	#-- packed units (0.1 micro-degree, 1e-7 degrees)
	L1b_c_location_parameters['Pitch'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Antenna Bench Yaw Angle (Derived from star trackers)
	#-- packed units (0.1 micro-degree, 1e-7 degrees)
	L1b_c_location_parameters['Yaw'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Measurement Confidence Data Flags
	#-- Generally the MCD flags indicate problems when set
	#-- If MCD is 0 then no problems or non-nominal conditions were detected
	#-- Serious errors are indicated by setting bit 31
	L1b_c_location_parameters['MCD'] = np.zeros((n_records,n_blocks),dtype=np.uint32)
	L1b_c_location_parameters['Spares'] = np.zeros((n_records,n_blocks,2),dtype=np.int16)

	#-- CryoSat-2 Measurement Group
	#-- Derived from instrument measurement parameters
	L1b_c_measurements = {}
	#-- Window Delay reference (two-way) corrected for instrument delays
	L1b_c_measurements['TD'] = np.zeros((n_records,n_blocks),dtype=np.int64)
	#-- H0 Initial Height Word from telemetry
	L1b_c_measurements['H_0'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- COR2 Height Rate: on-board tracker height rate over the radar cycle
	L1b_c_measurements['COR2'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Coarse Range Word (LAI) derived from telemetry
	L1b_c_measurements['LAI'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Fine Range Word (FAI) derived from telemetry
	L1b_c_measurements['FAI'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Automatic Gain Control Channel 1: AGC gain applied on Rx channel 1.
	#-- Gain calibration corrections are applied (Sum of AGC stages 1 and 2
	#-- plus the corresponding corrections) (dB/100)
	L1b_c_measurements['AGC_CH1'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Automatic Gain Control Channel 2: AGC gain applied on Rx channel 2.
	#-- Gain calibration corrections are applied (dB/100)
	L1b_c_measurements['AGC_CH2'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Total Fixed Gain On Channel 1: gain applied by the RF unit. (dB/100)
	L1b_c_measurements['TR_gain_CH1'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Total Fixed Gain On Channel 2: gain applied by the RF unit. (dB/100)
	L1b_c_measurements['TR_gain_CH2'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Transmit Power in microWatts
	L1b_c_measurements['TX_Power'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Doppler range correction: Radial component (mm)
	#-- computed for the component of satellite velocity in the nadir direction
	L1b_c_measurements['Doppler_range'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Instrument Range Correction: transmit-receive antenna (mm)
	#-- Calibration correction to range on channel 1 computed from CAL1.
	L1b_c_measurements['TR_inst_range'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Instrument Range Correction: receive-only antenna (mm)
	#-- Calibration correction to range on channel 2 computed from CAL1.
	L1b_c_measurements['R_inst_range'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Instrument Gain Correction: transmit-receive antenna (dB/100)
	#-- Calibration correction to gain on channel 1 computed from CAL1
	L1b_c_measurements['TR_inst_gain'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Instrument Gain Correction: receive-only (dB/100)
	#-- Calibration correction to gain on channel 2 computed from CAL1
	L1b_c_measurements['R_inst_gain'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Internal Phase Correction (microradians)
	L1b_c_measurements['Internal_phase'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- External Phase Correction (microradians)
	L1b_c_measurements['External_phase'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Noise Power measurement (dB/100)
	L1b_c_measurements['Noise_power'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Phase slope correction (microradians)
	#-- Computed from the CAL-4 packets during the azimuth impulse response
	#-- amplitude (SARIN only). Set from the latest available CAL-4 packet.
	L1b_c_measurements['Phase_slope'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	L1b_c_measurements['Spares1'] = np.zeros((n_records,n_blocks,4),dtype=np.int8)

	#-- CryoSat-2 External Corrections Group
	L1b_c_geo_corrections = {}
	#-- Dry Tropospheric Correction packed units (mm, 1e-3 m)
	L1b_c_geo_corrections['dryTrop'] = np.zeros((n_records),dtype=np.int32)
	#-- Wet Tropospheric Correction packed units (mm, 1e-3 m)
	L1b_c_geo_corrections['wetTrop'] = np.zeros((n_records),dtype=np.int32)
	#-- Inverse Barometric Correction packed units (mm, 1e-3 m)
	L1b_c_geo_corrections['InvBar'] = np.zeros((n_records),dtype=np.int32)
	#-- Delta Inverse Barometric Correction packed units (mm, 1e-3 m)
	L1b_c_geo_corrections['DAC'] = np.zeros((n_records),dtype=np.int32)
	#-- GIM Ionospheric Correction packed units (mm, 1e-3 m)
	L1b_c_geo_corrections['Iono_GIM'] = np.zeros((n_records),dtype=np.int32)
	#-- Model Ionospheric Correction packed units (mm, 1e-3 m)
	L1b_c_geo_corrections['Iono_model'] = np.zeros((n_records),dtype=np.int32)
	#-- Ocean tide Correction packed units (mm, 1e-3 m)
	L1b_c_geo_corrections['ocTideElv'] = np.zeros((n_records),dtype=np.int32)
	#-- Long period equilibrium ocean tide Correction packed units (mm, 1e-3 m)
	L1b_c_geo_corrections['lpeTideElv'] = np.zeros((n_records),dtype=np.int32)
	#-- Ocean loading tide Correction packed units (mm, 1e-3 m)
	L1b_c_geo_corrections['olTideElv'] = np.zeros((n_records),dtype=np.int32)
	#-- Solid Earth tide Correction packed units (mm, 1e-3 m)
	L1b_c_geo_corrections['seTideElv'] = np.zeros((n_records),dtype=np.int32)
	#-- Geocentric Polar tide Correction packed units (mm, 1e-3 m)
	L1b_c_geo_corrections['gpTideElv'] = np.zeros((n_records),dtype=np.int32)
	#-- Surface Type: enumerated key to classify surface at nadir
	#-- 0 = Open Ocean
	#-- 1 = Closed Sea
	#-- 2 = Continental Ice
	#-- 3 = Land
	L1b_c_geo_corrections['Surf_type'] = np.zeros((n_records),dtype=np.uint32)
	L1b_c_geo_corrections['Spare1'] = np.zeros((n_records,4),dtype=np.int8)
	#-- Corrections Status Flag
	L1b_c_geo_corrections['Corr_status'] = np.zeros((n_records),dtype=np.uint32)
	#-- Correction Error Flag
	L1b_c_geo_corrections['Corr_error'] = np.zeros((n_records),dtype=np.uint32)
	L1b_c_geo_corrections['Spare2'] = np.zeros((n_records,4),dtype=np.int8)

	#-- CryoSat-2 Average Waveforms Groups
	#-- Low-Resolution Mode
	L1b_c_1Hz_LRM_waveform = {}
	#-- Data Record Time (MDSR Time Stamp)
	L1b_c_1Hz_LRM_waveform['Day_1Hz'] = np.zeros((n_records),dtype=np.int32)
	L1b_c_1Hz_LRM_waveform['Sec_1Hz'] = np.zeros((n_records),dtype=np.uint32)
	L1b_c_1Hz_LRM_waveform['Micsec_1Hz'] = np.zeros((n_records),dtype=np.uint32)
	#-- Lat: packed units (0.1 micro-degree, 1e-7 degrees)
	L1b_c_1Hz_LRM_waveform['Lat_1Hz'] = np.zeros((n_records),dtype=np.int32)
	#-- Lon: packed units (0.1 micro-degree, 1e-7 degrees)
	L1b_c_1Hz_LRM_waveform['Lon_1Hz'] = np.zeros((n_records),dtype=np.int32)
	#-- Alt: packed units (mm, 1e-3 m)
	#-- Altitude of COG above reference ellipsoid (interpolated value)
	L1b_c_1Hz_LRM_waveform['Alt_1Hz'] = np.zeros((n_records),dtype=np.int32)
	#-- Window Delay (two-way) corrected for instrument delays
	L1b_c_1Hz_LRM_waveform['TD_1Hz'] = np.zeros((n_records),dtype=np.int64)
	#-- 1 Hz Averaged Power Echo Waveform
	L1b_c_1Hz_LRM_waveform['Waveform'] = np.zeros((n_records,n_LRM_RW),dtype=np.uint16)
	#-- Echo Scale Factor (to scale echo to watts)
	L1b_c_1Hz_LRM_waveform['Linear_Wfm_Multiplier'] = np.zeros((n_records),dtype=np.int32)
	#-- Echo Scale Power (a power of 2 to scale echo to Watts)
	L1b_c_1Hz_LRM_waveform['Power2_Wfm_Multiplier'] = np.zeros((n_records),dtype=np.int32)
	#-- Number of echoes averaged
	L1b_c_1Hz_LRM_waveform['N_avg_echoes'] = np.zeros((n_records),dtype=np.uint16)
	L1b_c_1Hz_LRM_waveform['Flags'] = np.zeros((n_records),dtype=np.uint16)

	#-- SAR Mode
	L1b_c_1Hz_SAR_waveform = {}
	#-- Data Record Time (MDSR Time Stamp)
	L1b_c_1Hz_SAR_waveform['Day_1Hz'] = np.zeros((n_records),dtype=np.int32)
	L1b_c_1Hz_SAR_waveform['Sec_1Hz'] = np.zeros((n_records),dtype=np.uint32)
	L1b_c_1Hz_SAR_waveform['Micsec_1Hz'] = np.zeros((n_records),dtype=np.uint32)
	#-- Lat: packed units (0.1 micro-degree, 1e-7 degrees)
	L1b_c_1Hz_SAR_waveform['Lat_1Hz'] = np.zeros((n_records),dtype=np.int32)
	#-- Lon: packed units (0.1 micro-degree, 1e-7 degrees)
	L1b_c_1Hz_SAR_waveform['Lon_1Hz'] = np.zeros((n_records),dtype=np.int32)
	#-- Alt: packed units (mm, 1e-3 m)
	#-- Altitude of COG above reference ellipsoid (interpolated value)
	L1b_c_1Hz_SAR_waveform['Alt_1Hz'] = np.zeros((n_records),dtype=np.int32)
	#-- Window Delay (two-way) corrected for instrument delays
	L1b_c_1Hz_SAR_waveform['TD_1Hz'] = np.zeros((n_records),dtype=np.int64)
	#-- 1 Hz Averaged Power Echo Waveform
	L1b_c_1Hz_SAR_waveform['Waveform'] = np.zeros((n_records,n_SAR_RW),dtype=np.uint16)
	#-- Echo Scale Factor (to scale echo to watts)
	L1b_c_1Hz_SAR_waveform['Linear_Wfm_Multiplier'] = np.zeros((n_records),dtype=np.int32)
	#-- Echo Scale Power (a power of 2 to scale echo to Watts)
	L1b_c_1Hz_SAR_waveform['Power2_Wfm_Multiplier'] = np.zeros((n_records),dtype=np.int32)
	#-- Number of echoes averaged
	L1b_c_1Hz_SAR_waveform['N_avg_echoes'] = np.zeros((n_records),dtype=np.uint16)
	L1b_c_1Hz_SAR_waveform['Flags'] = np.zeros((n_records),dtype=np.uint16)

	#-- SARIN Mode
	#-- Same as the LRM/SAR groups but the waveform array is 512 bins instead of
	#-- 128 and the number of echoes averaged is different.
	L1b_c_1Hz_SARIN_waveform = {}
	#-- Data Record Time (MDSR Time Stamp)
	L1b_c_1Hz_SARIN_waveform['Day'] = np.zeros((n_records),dtype=np.int32)
	L1b_c_1Hz_SARIN_waveform['Sec'] = np.zeros((n_records),dtype=np.uint32)
	L1b_c_1Hz_SARIN_waveform['Micsec'] = np.zeros((n_records),dtype=np.uint32)
	#-- Lat: packed units (0.1 micro-degree, 1e-7 degrees)
	L1b_c_1Hz_SARIN_waveform['Lat'] = np.zeros((n_records),dtype=np.int32)
	#-- Lon: packed units (0.1 micro-degree, 1e-7 degrees)
	L1b_c_1Hz_SARIN_waveform['Lon'] = np.zeros((n_records),dtype=np.int32)
	#-- Alt: packed units (mm, 1e-3 m)
	#-- Altitude of COG above reference ellipsoid (interpolated value)
	L1b_c_1Hz_SARIN_waveform['Alt'] = np.zeros((n_records),dtype=np.int32)
	#-- Window Delay (two-way) corrected for instrument delays
	L1b_c_1Hz_SARIN_waveform['TD'] = np.zeros((n_records),dtype=np.int64)
	#-- 1 Hz Averaged Power Echo Waveform
	L1b_c_1Hz_SARIN_waveform['Waveform'] = np.zeros((n_records,n_SARIN_RW),dtype=np.uint16)
	#-- Echo Scale Factor (to scale echo to watts)
	L1b_c_1Hz_SARIN_waveform['Linear_Wfm_Multiplier'] = np.zeros((n_records),dtype=np.int32)
	#-- Echo Scale Power (a power of 2 to scale echo to Watts)
	L1b_c_1Hz_SARIN_waveform['Power2_Wfm_Multiplier'] = np.zeros((n_records),dtype=np.int32)
	#-- Number of echoes averaged
	L1b_c_1Hz_SARIN_waveform['N_avg_echoes'] = np.zeros((n_records),dtype=np.uint16)
	L1b_c_1Hz_SARIN_waveform['Flags'] = np.zeros((n_records),dtype=np.uint16)

	#-- CryoSat-2 Waveforms Groups
	#-- Beam Behavior Parameters
	L1b_c_Beam_Behavior = {}
	#-- Standard Deviation of Gaussian fit to range integrated stack power.
	L1b_c_Beam_Behavior['SD'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
	#-- Stack Center: Mean of Gaussian fit to range integrated stack power.
	L1b_c_Beam_Behavior['Center'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
	#-- Stack amplitude parameter scaled in dB/100.
	L1b_c_Beam_Behavior['Amplitude'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
	#-- 3rd moment: providing the degree of asymmetry of the range integrated
	#-- stack power distribution.
	L1b_c_Beam_Behavior['Skewness'] = np.zeros((n_records,n_blocks),dtype=np.int16)
	#-- 4th moment: Measure of peakiness of range integrated stack power distribution.
	L1b_c_Beam_Behavior['Kurtosis'] = np.zeros((n_records,n_blocks),dtype=np.int16)
	#-- Standard deviation as a function of boresight angle (microradians)
	L1b_c_Beam_Behavior['SD_boresight_angle'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
	#-- Stack Center angle as a function of boresight angle (microradians)
	L1b_c_Beam_Behavior['Center_boresight_angle'] = np.zeros((n_records,n_blocks),dtype=np.int16)
	L1b_c_Beam_Behavior['Spare'] = np.zeros((n_records,n_blocks,n_BeamBehaviourParams-7),dtype=np.int16)

	#-- Low-Resolution Mode
	L1b_c_LRM_waveform = {}
	#-- Averaged Power Echo Waveform [128]
	L1b_c_LRM_waveform['Waveform'] = np.zeros((n_records,n_blocks,n_LRM_RW),dtype=np.uint16)
	#-- Echo Scale Factor (to scale echo to watts)
	L1b_c_LRM_waveform['Linear_Wfm_Multiplier'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Echo Scale Power (a power of 2)
	L1b_c_LRM_waveform['Power2_Wfm_Multiplier'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Number of echoes averaged
	L1b_c_LRM_waveform['N_avg_echoes'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
	L1b_c_LRM_waveform['Flags'] = np.zeros((n_records,n_blocks),dtype=np.uint16)

	#-- SAR Mode
	L1b_c_SAR_waveform = {}
	#-- Averaged Power Echo Waveform [256]
	L1b_c_SAR_waveform['Waveform'] = np.zeros((n_records,n_blocks,n_SAR_BC_RW),dtype=np.uint16)
	#-- Echo Scale Factor (to scale echo to watts)
	L1b_c_SAR_waveform['Linear_Wfm_Multiplier'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Echo Scale Power (a power of 2)
	L1b_c_SAR_waveform['Power2_Wfm_Multiplier'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Number of echoes averaged
	L1b_c_SAR_waveform['N_avg_echoes'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
	L1b_c_SAR_waveform['Flags'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
	#-- Beam behaviour parameters
	L1b_c_SAR_waveform['Beam'] = L1b_c_Beam_Behavior

	#-- SARIN Mode
	L1b_c_SARIN_waveform = {}
	#-- Averaged Power Echo Waveform [1024]
	L1b_c_SARIN_waveform['Waveform'] = np.zeros((n_records,n_blocks,n_SARIN_BC_RW),dtype=np.uint16)
	#-- Echo Scale Factor (to scale echo to watts)
	L1b_c_SARIN_waveform['Linear_Wfm_Multiplier'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Echo Scale Power (a power of 2)
	L1b_c_SARIN_waveform['Power2_Wfm_Multiplier'] = np.zeros((n_records,n_blocks),dtype=np.int32)
	#-- Number of echoes averaged
	L1b_c_SARIN_waveform['N_avg_echoes'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
	L1b_c_SARIN_waveform['Flags'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
	#-- Beam behaviour parameters
	L1b_c_SARIN_waveform['Beam'] = L1b_c_Beam_Behavior
	#-- Coherence [1024]: packed units (1/1000)
	L1b_c_SARIN_waveform['Coherence'] = np.zeros((n_records,n_blocks,n_SARIN_BC_RW),dtype=np.int16)
	#-- Phase Difference [1024]: packed units (microradians)
	L1b_c_SARIN_waveform['Phase_diff'] = np.zeros((n_records,n_blocks,n_SARIN_BC_RW),dtype=np.int32)

	#-- for each record in the CryoSat file
	for r in range(n_records):
		#-- CryoSat-2 Time and Orbit Group
		for b in range(n_blocks):
			L1b_c_location_parameters['Day'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_location_parameters['Second'][r,b] = np.fromfile(fid,dtype='>u4',count=1)
			L1b_c_location_parameters['Micsec'][r,b] = np.fromfile(fid,dtype='>u4',count=1)
			L1b_c_location_parameters['USO_Corr'][r,b] = np.fromfile(fid,dtype='>u4',count=1)
			L1b_c_location_parameters['Mode_ID'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
			L1b_c_location_parameters['SSC'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
			L1b_c_location_parameters['Inst_config'][r,b] = np.fromfile(fid,dtype='>u4',count=1)
			L1b_c_location_parameters['Rec_Count'][r,b] = np.fromfile(fid,dtype='>u4',count=1)
			L1b_c_location_parameters['Lat'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_location_parameters['Lon'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_location_parameters['Alt'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_location_parameters['Alt_rate'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_location_parameters['Sat_velocity'][r,b,:] = np.fromfile(fid,dtype='>i4',count=3)
			L1b_c_location_parameters['Real_beam'][r,b,:] = np.fromfile(fid,dtype='>i4',count=3)
			L1b_c_location_parameters['Baseline'][r,b,:] = np.fromfile(fid,dtype='>i4',count=3)
			L1b_c_location_parameters['ST_ID'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
			L1b_c_location_parameters['Roll'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_location_parameters['Pitch'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_location_parameters['Yaw'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_location_parameters['MCD'][r,b] = np.fromfile(fid,dtype='>u4',count=1)
			L1b_c_location_parameters['Spares'][r,b,:] = np.fromfile(fid,dtype='>i2',count=2)

		#-- CryoSat-2 Measurement Group
		#-- Derived from instrument measurement parameters
		for b in range(n_blocks):
			L1b_c_measurements['TD'][r,b] = np.fromfile(fid,dtype='>i8',count=1)
			L1b_c_measurements['H_0'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_measurements['COR2'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_measurements['LAI'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_measurements['FAI'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_measurements['AGC_CH1'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_measurements['AGC_CH2'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_measurements['TR_gain_CH1'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_measurements['TR_gain_CH2'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_measurements['TX_Power'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_measurements['Doppler_range'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_measurements['TR_inst_range'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_measurements['R_inst_range'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_measurements['TR_inst_gain'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_measurements['R_inst_gain'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_measurements['Internal_phase'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_measurements['External_phase'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_measurements['Noise_power'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_measurements['Phase_slope'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_measurements['Spares1'][r,b,:] = np.fromfile(fid,dtype='>i1',count=4)

		#-- CryoSat-2 External Corrections Group
		L1b_c_geo_corrections['dryTrop'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L1b_c_geo_corrections['wetTrop'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L1b_c_geo_corrections['InvBar'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L1b_c_geo_corrections['DAC'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L1b_c_geo_corrections['Iono_GIM'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L1b_c_geo_corrections['Iono_model'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L1b_c_geo_corrections['ocTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L1b_c_geo_corrections['lpeTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L1b_c_geo_corrections['olTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L1b_c_geo_corrections['seTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L1b_c_geo_corrections['gpTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
		L1b_c_geo_corrections['Surf_type'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L1b_c_geo_corrections['Spare1'][r,:] = np.fromfile(fid,dtype='>i1',count=4)
		L1b_c_geo_corrections['Corr_status'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L1b_c_geo_corrections['Corr_error'][r] = np.fromfile(fid,dtype='>u4',count=1)
		L1b_c_geo_corrections['Spare2'][r,:] = np.fromfile(fid,dtype='>i1',count=4)

		#-- CryoSat-2 Average Waveforms Groups
		if (MODE == 'LRM'):
			#-- Low-Resolution Mode
			L1b_c_1Hz_LRM_waveform['Day_1Hz'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_1Hz_LRM_waveform['Sec_1Hz'][r] = np.fromfile(fid,dtype='>u4',count=1)
			L1b_c_1Hz_LRM_waveform['Micsec_1Hz'][r] = np.fromfile(fid,dtype='>u4',count=1)
			L1b_c_1Hz_LRM_waveform['Lat_1Hz'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_1Hz_LRM_waveform['Lon_1Hz'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_1Hz_LRM_waveform['Alt_1Hz'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_1Hz_LRM_waveform['TD_1Hz'][r] = np.fromfile(fid,dtype='>i8',count=1)
			L1b_c_1Hz_LRM_waveform['Waveform'][r,:] = np.fromfile(fid,dtype='>u2',count=n_LRM_RW)
			L1b_c_1Hz_LRM_waveform['Linear_Wfm_Multiplier'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_1Hz_LRM_waveform['Power2_Wfm_Multiplier'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_1Hz_LRM_waveform['N_avg_echoes'][r] = np.fromfile(fid,dtype='>u2',count=1)
			L1b_c_1Hz_LRM_waveform['Flags'][r] = np.fromfile(fid,dtype='>u2',count=1)
		elif (MODE == 'SAR'):
			#-- SAR Mode
			L1b_c_1Hz_SAR_waveform['Day_1Hz'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_1Hz_SAR_waveform['Sec_1Hz'][r] = np.fromfile(fid,dtype='>u4',count=1)
			L1b_c_1Hz_SAR_waveform['Micsec_1Hz'][r] = np.fromfile(fid,dtype='>u4',count=1)
			L1b_c_1Hz_SAR_waveform['Lat_1Hz'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_1Hz_SAR_waveform['Lon_1Hz'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_1Hz_SAR_waveform['Alt_1Hz'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_1Hz_SAR_waveform['TD_1Hz'][r] = np.fromfile(fid,dtype='>i8',count=1)
			L1b_c_1Hz_SAR_waveform['Waveform'][r,:] = np.fromfile(fid,dtype='>u2',count=n_SAR_RW)
			L1b_c_1Hz_SAR_waveform['Linear_Wfm_Multiplier'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_1Hz_SAR_waveform['Power2_Wfm_Multiplier'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_1Hz_SAR_waveform['N_avg_echoes'][r] = np.fromfile(fid,dtype='>u2',count=1)
			L1b_c_1Hz_SAR_waveform['Flags'][r] = np.fromfile(fid,dtype='>u2',count=1)
		elif (MODE == 'SIN'):
			#-- SARIN Mode
			L1b_c_1Hz_SARIN_waveform['Day'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_1Hz_SARIN_waveform['Sec'][r] = np.fromfile(fid,dtype='>u4',count=1)
			L1b_c_1Hz_SARIN_waveform['Micsec'][r] = np.fromfile(fid,dtype='>u4',count=1)
			L1b_c_1Hz_SARIN_waveform['Lat'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_1Hz_SARIN_waveform['Lon'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_1Hz_SARIN_waveform['Alt'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_1Hz_SARIN_waveform['TD'][r] = np.fromfile(fid,dtype='>i8',count=1)
			L1b_c_1Hz_SARIN_waveform['Waveform'][r,:] = np.fromfile(fid,dtype='>u2',count=n_SARIN_RW)
			L1b_c_1Hz_SARIN_waveform['Linear_Wfm_Multiplier'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_1Hz_SARIN_waveform['Power2_Wfm_Multiplier'][r] = np.fromfile(fid,dtype='>i4',count=1)
			L1b_c_1Hz_SARIN_waveform['N_avg_echoes'][r] = np.fromfile(fid,dtype='>u2',count=1)
			L1b_c_1Hz_SARIN_waveform['Flags'][r] = np.fromfile(fid,dtype='>u2',count=1)

		#-- CryoSat-2 Waveforms Groups
		if (MODE == 'LRM'):
			#-- Low-Resolution Mode
			for b in range(n_blocks):
				L1b_c_LRM_waveform['Waveform'][r,b,:] = np.fromfile(fid,dtype='>u2',count=n_LRM_RW)
				L1b_c_LRM_waveform['Linear_Wfm_Multiplier'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
				L1b_c_LRM_waveform['Power2_Wfm_Multiplier'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
				L1b_c_LRM_waveform['N_avg_echoes'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
				L1b_c_LRM_waveform['Flags'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
		elif (MODE == 'SAR'):
			#-- SAR Mode
			for b in range(n_blocks):
				L1b_c_SAR_waveform['Waveform'][r,b,:] = np.fromfile(fid,dtype='>u2',count=n_SAR_BC_RW)
				L1b_c_SAR_waveform['Linear_Wfm_Multiplier'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
				L1b_c_SAR_waveform['Power2_Wfm_Multiplier'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
				L1b_c_SAR_waveform['N_avg_echoes'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
				L1b_c_SAR_waveform['Flags'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
				L1b_c_SAR_waveform['Beam']['SD'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
				L1b_c_SAR_waveform['Beam']['Center'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
				L1b_c_SAR_waveform['Beam']['Amplitude'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
				L1b_c_SAR_waveform['Beam']['Skewness'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
				L1b_c_SAR_waveform['Beam']['Kurtosis'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
				L1b_c_SAR_waveform['Beam']['SD_boresight_angle'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
				L1b_c_SAR_waveform['Beam']['Center_boresight_angle'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
				L1b_c_SAR_waveform['Beam']['Spare'][r,b,:] = np.fromfile(fid,dtype='>i2',count=(n_BeamBehaviourParams-7))
		elif (MODE == 'SIN'):
			#-- SARIN Mode
			for b in range(n_blocks):
				L1b_c_SARIN_waveform['Waveform'][r,b,:] = np.fromfile(fid,dtype='>u2',count=n_SARIN_BC_RW)
				L1b_c_SARIN_waveform['Linear_Wfm_Multiplier'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
				L1b_c_SARIN_waveform['Power2_Wfm_Multiplier'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
				L1b_c_SARIN_waveform['N_avg_echoes'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
				L1b_c_SARIN_waveform['Flags'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
				L1b_c_SARIN_waveform['Beam']['SD'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
				L1b_c_SARIN_waveform['Beam']['Center'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
				L1b_c_SARIN_waveform['Beam']['Amplitude'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
				L1b_c_SARIN_waveform['Beam']['Skewness'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
				L1b_c_SARIN_waveform['Beam']['Kurtosis'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
				L1b_c_SARIN_waveform['Beam']['SD_boresight_angle'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
				L1b_c_SARIN_waveform['Beam']['Center_boresight_angle'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
				L1b_c_SARIN_waveform['Beam']['Spare'][r,b,:] = np.fromfile(fid,dtype='>i2',count=(n_BeamBehaviourParams-7))
				L1b_c_SARIN_waveform['Coherence'][r,b,:] = np.fromfile(fid,dtype='>i2',count=n_SARIN_BC_RW)
				L1b_c_SARIN_waveform['Phase_diff'][r,b,:] = np.fromfile(fid,dtype='>i4',count=n_SARIN_BC_RW)

	#-- Bind all the bits of the l1b_mds together into a single dictionary
	CS_l1b_c_mds = {}
	CS_l1b_c_mds['Location'] = L1b_c_location_parameters
	CS_l1b_c_mds['Data'] = L1b_c_measurements
	CS_l1b_c_mds['Geometry'] = L1b_c_geo_corrections
	if (MODE == 'LRM'):
		CS_l1b_c_mds['Waveform_1Hz'] = L1b_c_1Hz_LRM_waveform
		CS_l1b_c_mds['Waveform_20Hz'] = L1b_c_LRM_waveform
	elif (MODE == 'SAR'):
		CS_l1b_c_mds['Waveform_1Hz'] = L1b_c_1Hz_SAR_waveform
		CS_l1b_c_mds['Waveform_20Hz'] = L1b_c_SAR_waveform
	elif (MODE == 'SIN'):
		CS_l1b_c_mds['Waveform_1Hz'] = L1b_c_1Hz_SARIN_waveform
		CS_l1b_c_mds['Waveform_20Hz'] = L1b_c_SARIN_waveform

	#-- return the output dictionary
	return CS_l1b_c_mds

#-- PURPOSE: Read ASCII Main Product Header (MPH) block from an ESA PDS file
def read_MPH(full_filename):
	#-- read input data file
	with open(full_filename, 'rb') as fid:
		file_contents = fid.read().splitlines()

	#-- Define constant values associated with PDS file formats
	#-- number of text lines in standard MPH
	n_MPH_lines	= 41
	#-- check that first line of header matches PRODUCT
	if not bool(re.match('PRODUCT\=\"(.*)(?=\")',file_contents[0])):
		raise IOError('File does not start with a valid PDS MPH')
	#-- read MPH header text
	s_MPH_fields = {}
	for i in range(n_MPH_lines):
		#-- use regular expression operators to read headers
		if bool(re.match('(.*?)\=\"(.*)(?=\")',file_contents[i])):
			#-- data fields within quotes
			field,value=re.findall('(.*?)\=\"(.*)(?=\")',file_contents[i]).pop()
			s_MPH_fields[field] = value.rstrip()
		elif bool(re.match('(.*?)\=(.*)',file_contents[i])):
			#-- data fields without quotes
			field,value=re.findall('(.*?)\=(.*)',file_contents[i]).pop()
			s_MPH_fields[field] = value.rstrip()

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
	rx = re.compile('(.*?)\=\"?(.*)',re.VERBOSE)
	#-- check first line of header matches SPH_DESCRIPTOR
	if not bool(re.match('SPH\_DESCRIPTOR\=',file_contents[n_MPH_lines+1])):
		raise IOError('File does not have a valid PDS DSD')
	#-- read SPH header text (no binary control characters)
	s_SPH_lines = [li for li in file_contents[n_MPH_lines+1:] if rx.match(li)
		and not re.search(r'[^\x20-\x7e]+',li)]

	#-- extract SPH header text
	s_SPH_fields = {}
	c = 0
	while (c < len(s_SPH_lines)):
		#-- check if line is within DS_NAME portion of SPH header
		if bool(re.match('DS_NAME',s_SPH_lines[c])):
			#-- add dictionary for DS_NAME
			field,value=re.findall('(.*?)\=\"(.*)(?=\")',s_SPH_lines[c]).pop()
			s_SPH_fields[value.rstrip()] = {}
			for line in s_SPH_lines[c+1:c+7]:
				if bool(re.match('(.*?)\=\"(.*)(?=\")',line)):
					#-- data fields within quotes
					dsfield,dsvalue=re.findall('(.*?)\=\"(.*)(?=\")',line).pop()
					s_SPH_fields[value.rstrip()][dsfield] = dsvalue.rstrip()
				elif bool(re.match('(.*?)\=(.*)',line)):
					#-- data fields without quotes
					dsfield,dsvalue=re.findall('(.*?)\=(.*)',line).pop()
					s_SPH_fields[value.rstrip()][dsfield] = dsvalue.rstrip()
			#-- add 6 to counter to go to next entry
			c += 6
		#-- use regular expression operators to read headers
		elif bool(re.match('(.*?)\=\"(.*)(?=\")',s_SPH_lines[c])):
			#-- data fields within quotes
			field,value=re.findall('(.*?)\=\"(.*)(?=\")',s_SPH_lines[c]).pop()
			s_SPH_fields[field] = value.rstrip()
		elif bool(re.match('(.*?)\=(.*)',s_SPH_lines[c])):
			#-- data fields without quotes
			field,value=re.findall('(.*?)\=(.*)',s_SPH_lines[c]).pop()
			s_SPH_fields[field] = value.rstrip()
		#-- add 1 to counter to go to next line
		c += 1

	#-- Return block name array to calling function
	return s_SPH_fields

#-- PURPOSE: Read ASCII Data Set Descriptors (DSD) block from a PDS file
def read_DSD(full_filename, DS_TYPE=None):
	#-- read input data file
	with open(full_filename, 'rb') as fid:
		file_contents = fid.read().splitlines()

	#-- Define constant values associated with PDS file formats
	#-- number of text lines in standard MPH
	n_MPH_lines	= 41
	#-- number of text lines in a DSD header
	n_DSD_lines = 8

	#-- Level-1b CryoSat DS_NAMES within files
	regex_patterns = []
	if (DS_TYPE == 'CS_L1B'):
		regex_patterns.append('DS_NAME\="SIR_L1B_LRM[\s+]*"')
		regex_patterns.append('DS_NAME\="SIR_L1B_SAR[\s+]*"')
		regex_patterns.append('DS_NAME\="SIR_L1B_SARIN[\s+]*"')
	elif (DS_TYPE == 'SIR_L1B_FDM'):
		regex_patterns.append('DS_NAME\="SIR_L1B_FDM[\s+]*"')
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
		if bool(re.match('(.*?)\=\"(.*)(?=\")',file_contents[i])):
			#-- data fields within quotes
			field,value=re.findall('(.*?)\=\"(.*)(?=\")',file_contents[i]).pop()
			s_DSD_fields[field] = value.rstrip()
		elif bool(re.match('(.*?)\=(.*)',file_contents[i])):
			#-- data fields without quotes
			field,value=re.findall('(.*?)\=(.*)',file_contents[i]).pop()
			s_DSD_fields[field] = value.rstrip()

	#-- Return block name array to calling function
	return s_DSD_fields

#-- PURPOSE: read CryoSat Level-1b data
def read_cryosat_L1b(full_filename, VERBOSE=False):
	#-- file basename and file extension of input file
	fileBasename,fileExtension=os.path.splitext(os.path.basename(full_filename))

	#-- CryoSat file class
	#-- OFFL (Off Line Processing/Systematic)
	#-- NRT_ (Near Real Time)
	#-- RPRO (ReProcessing)
	#-- TEST (Testing)
	#-- TIxx (Stand alone IPF1 testing)
	#-- LTA_ (Long Term Archive)
	regex_class = 'OFFL|NRT_|RPRO|TEST|TIxx|LTA_'
	#-- CryoSat mission products
	#-- SIR1SAR_FR: Level 1 FBR SAR Mode (Rx1 Channel)
	#-- SIR2SAR_FR: Level 1 FBR SAR Mode (Rx2 Channel)
	#-- SIR_SIN_FR: Level 1 FBR SARin Mode
	#-- SIR_LRM_1B: Level-1 Product Low Rate Mode
	#-- SIR_FDM_1B: Level-1 Product Fast Delivery Marine Mode
	#-- SIR_SAR_1B: Level-1 SAR Mode
	#-- SIR_SIN_1B: Level-1 SARin Mode
	#-- SIR1LRC11B: Level-1 CAL1 Low Rate Mode (Rx1 Channel)
	#-- SIR2LRC11B: Level-1 CAL1 Low Rate Mode (Rx2 Channel)
	#-- SIR1SAC11B: Level-1 CAL1 SAR Mode (Rx1 Channel)
	#-- SIR2SAC11B: Level-1 CAL1 SAR Mode (Rx2 Channel)
	#-- SIR_SIC11B: Level-1 CAL1 SARin Mode
	#-- SIR_SICC1B: Level-1 CAL1 SARIN Exotic Data
	#-- SIR1SAC21B: Level-1 CAL2 SAR Mode (Rx1 Channel)
	#-- SIR2SAC21B: Level-1 CAL2 SAR Mode (Rx2 Channel)
	#-- SIR1SIC21B: Level-1 CAL2 SARin Mode (Rx1 Channel)
	#-- SIR2SIC21B: Level-1 CAL2 SARin Mode (Rx1 Channel)
	#-- SIR1LRM_0M: LRM and TRK Monitoring Data from Rx 1 Channel
	#-- SIR2LRM_0M: LRM and TRK Monitoring Data from Rx 2 Channel
	#-- SIR1SAR_0M: SAR Monitoring Data from Rx 1 Channel
	#-- SIR2SAR_0M: SAR Monitoring Data from Rx 1 Channel
	#-- SIR_SIN_0M: SARIN Monitoring Data
	#-- SIR_SIC40M: CAL4 Monitoring Data
	regex_products = ('SIR1SAR_FR|SIR2SAR_FR|SIR_SIN_FR|SIR_LRM_1B|SIR_FDM_1B|'
		'SIR_SAR_1B|SIR_SIN_1B|SIR1LRC11B|SIR2LRC11B|SIR1SAC11B|SIR2SAC11B|'
		'SIR_SIC11B|SIR_SICC1B|SIR1SAC21B|SIR2SAC21B|SIR1SIC21B|SIR2SIC21B|'
		'SIR1LRM_0M|SIR2LRM_0M|SIR1SAR_0M|SIR2SAR_0M|SIR_SIN_0M|SIR_SIC40M')
	#-- CRYOSAT LEVEL-1b PRODUCTS NAMING RULES
	#-- Mission Identifier
	#-- File Class
	#-- File Product
	#-- Validity Start Date and Time
	#-- Validity Stop Date and Time
	#-- Baseline Identifier
	#-- Version Number
	regex_pattern = '(.*?)_({0})_({1})_(\d+T?\d+)_(\d+T?\d+)_(.*?)(\d+)'.format(
		regex_class, regex_products)
	rx = re.compile(regex_pattern, re.VERBOSE)
	#-- extract file information from filename
	MI,CLASS,PRODUCT,START,STOP,BASELINE,VERSION=rx.findall(fileBasename).pop()
	#-- Extract Date information
	start_yr,start_mon,start_day=np.array([START[:4],START[4:6],START[6:8]],dtype=np.uint16)
	start_hh,start_mm,start_ss=np.array([START[-6:-4],START[-4:-2],START[-2:]],dtype=np.uint8)
	stop_yr,stop_mon,stop_day=np.array([STOP[:4],STOP[4:6],STOP[6:8]],dtype=np.uint16)
	stop_hh,stop_mm,stop_ss=np.array([STOP[-6:-4],STOP[-4:-2],STOP[-2:]],dtype=np.uint8)

	#-- CryoSat-2 Mode record sizes
	i_size_timestamp = 12
	n_SARIN_BC_RW = 1024
	n_SARIN_RW = 512
	n_SAR_BC_RW = 256
	n_SAR_RW = 125
	n_LRM_RW = 128
	n_blocks = 20
	n_BeamBehaviourParams = 50
	#-- check baseline from file to set i_record_size and allocation function
	if (BASELINE == 'C'):
		#-- calculate total record sizes of each dataset group
		i_size_timegroup = i_size_timestamp + 4 + 2*2 + 6*4 + 3*3*4 + 3*2 + 4*4
		i_size_measuregroup = 8 + 4*17 + 8
		i_size_external_corr = 4*13 + 12
		i_size_1Hz_LRM = i_size_timestamp + 3*4 + 8 + n_LRM_RW*2 + 2*4 + 2*2
		i_size_1Hz_SAR = i_size_timestamp + 4*3 + 8 + n_SAR_RW*2 + 4 + 4 + 2 + 2
		i_size_1Hz_SARIN = i_size_timestamp + 4*3 + 8 + n_SARIN_RW*2 + 4 + 4 + 2 + 2
		i_size_LRM_waveform = n_LRM_RW*2 + 4 + 4 + 2 + 2
		i_size_SAR_waveform = n_SAR_BC_RW*2 + 4 + 4 + 2 + 2 + n_BeamBehaviourParams*2
		i_size_SARIN_waveform = n_SARIN_BC_RW*2 + 4 + 4 + 2 + 2 + n_SARIN_BC_RW*2 + \
			n_SARIN_BC_RW*4 + n_BeamBehaviourParams*2
		#-- Low-Resolution Mode Record Size
		i_record_size_LRM_L1b = n_blocks * (i_size_timegroup + \
			i_size_measuregroup + i_size_LRM_waveform) + i_size_external_corr + \
			i_size_1Hz_LRM
		#-- SAR Mode Record Size
		i_record_size_SAR_L1b = n_blocks * (i_size_timegroup + \
			i_size_measuregroup + i_size_SAR_waveform) + i_size_external_corr + \
			i_size_1Hz_SAR
		#-- SARIN Mode Record Size
		i_record_size_SARIN_L1b = n_blocks * (i_size_timegroup + \
			i_size_measuregroup + i_size_SARIN_waveform) + i_size_external_corr + \
			i_size_1Hz_SARIN
		#-- set read function for Baseline C
		read_cryosat_variables = cryosat_baseline_C
	else:
		#-- calculate total record sizes of each dataset group
		i_size_timegroup = i_size_timestamp + 4 + 2*2+ 6*4 + 3*3*4 + 4
		i_size_measuregroup = 8 + 4*17 + 8
		i_size_external_corr = 4*13 + 12
		i_size_1Hz_LRM = i_size_timestamp + 3*4 + 8 + n_LRM_RW*2 + 2*4 + 2*2
		i_size_1Hz_SAR = i_size_timestamp + 4*3 + 8 + n_SAR_RW*2 + 4 + 4 + 2 + 2
		i_size_1Hz_SARIN = i_size_timestamp + 4*3 + 8 + n_SARIN_RW*2 + 4 + 4 + 2 + 2
		i_size_LRM_waveform = n_LRM_RW*2 + 4 + 4 + 2 + 2
		i_size_SAR_waveform = n_SAR_RW*2 + 4 + 4 + 2 + 2 + n_BeamBehaviourParams*2
		i_size_SARIN_waveform = n_SARIN_RW*2 + 4 + 4 + 2 + 2 + n_SARIN_RW*2 + \
			n_SARIN_RW*4 + n_BeamBehaviourParams*2
		#-- Low-Resolution Mode Record Size
		i_record_size_LRM_L1b = n_blocks * (i_size_timegroup + \
			i_size_measuregroup + i_size_LRM_waveform) + i_size_external_corr + \
			i_size_1Hz_LRM
		#-- SAR Mode Record Size
		i_record_size_SAR_L1b = n_blocks * (i_size_timegroup + \
			i_size_measuregroup + i_size_SAR_waveform) + i_size_external_corr + \
			i_size_1Hz_SAR
		#-- SARIN Mode Record Size
		i_record_size_SARIN_L1b = n_blocks * (i_size_timegroup + \
			i_size_measuregroup + i_size_SARIN_waveform) + i_size_external_corr + \
			i_size_1Hz_SARIN
		#-- set read function for Baselines A and B
		read_cryosat_variables = cryosat_baseline_AB

	#-- get dataset MODE from PRODUCT portion of file name
	#-- set record sizes and DS_TYPE for read_DSD function
	MODE = re.findall('(LRM|FDM|SAR|SIN)', PRODUCT).pop()
	if (MODE == 'LRM'):
		i_record_size = i_record_size_LRM_L1b
		DS_TYPE = 'CS_L1B'
	elif (MODE == 'FDM'):
		i_record_size = i_record_size_FDM_L1b
		DS_TYPE = 'SIR_L1B_FDM'
	elif (MODE == 'SAR'):
		i_record_size = i_record_size_SAR_L1b
		DS_TYPE = 'CS_L1B'
	elif (MODE == 'SIN'):
		i_record_size = i_record_size_SARIN_L1b
		DS_TYPE = 'CS_L1B'

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
		s_SPH_fields = read_SPH(full_filename, j_sph_size)
		#-- extract information from DSD fields
		s_DSD_fields = read_DSD(full_filename, DS_TYPE=DS_TYPE)
		#-- extract DS_OFFSET
		j_DS_start = np.int32(re.findall('[-+]?\d+',s_DSD_fields['DS_OFFSET']).pop())
		#-- extract number of DSR in the file
		j_num_DSR = np.int32(re.findall('[-+]?\d+',s_DSD_fields['NUM_DSR']).pop())
		#-- check the record size
		j_DSR_size = np.int32(re.findall('[-+]?\d+',s_DSD_fields['DSR_SIZE']).pop())
		#--  minimum size is start of the read plus number of records to read
		j_check_size = j_DS_start + (j_DSR_size*j_num_DSR)
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
		CS_L1b_mds = read_cryosat_variables(fid, j_num_DSR, MODE)
		#-- add headers to output dictionary as METADATA
		CS_L1b_mds['METADATA'] = {}
		CS_L1b_mds['METADATA']['MPH'] = s_MPH_fields
		CS_L1b_mds['METADATA']['SPH'] = s_SPH_fields
		CS_L1b_mds['METADATA']['DSD'] = s_DSD_fields
		#-- close the input CryoSat binary file
		fid.close()
	else:
		#-- If there are not MPH/SPH/DSD headers
		#-- extract binary data from input CryoSat data file
		fid = open(full_filename, 'rb')
		#-- iterate through CryoSat file and fill output variables
		CS_L1b_mds = read_cryosat_variables(fid, j_num_DSR, MODE)
		#-- close the input CryoSat binary file
		fid.close()

	#-- return the data and headers
	return CS_L1b_mds
