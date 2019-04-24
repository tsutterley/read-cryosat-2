#!/usr/bin/env python
u"""
HDF5_cryosat_L2.py (04/2019)
Reads and Writes HDF5 files for CryoSat-2 Level-2 data products
Supported CryoSat Modes: LRM, SAR, SARin, FDM, SID, GDR

OUTPUTS a formatted HDF5 file with:
	Data_1Hz: Time and Orbit Parameters
	Corrections: Elevation Corrections and Flags
	Data_20Hz: Geolocation and Elevation Measurements with Quality Parameters
	METADATA: MPH, SPH and DSD Header data

OPTIONS:
	BASELINE (HDF5_cryosat_L2): CryoSat-2 baseline (A, B, C)
	FILENAME (HDF5_cryosat_L2): output HDF5 file name
	TITLE (HDF5_cryosat_L2): output file description
	HEADER (HDF5_cryosat_L2): output CryoSat-2 file headers (MPH, SPH, DSD)
		1: for single CryoSat-2 files
		2: for merged CryoSat-2 files from convert_cryosat_L2.py
	CLOBBER (HDF5_cryosat_L2): overwrite existing HDF5 file
	VERBOSE: print HDF5 structure parameters to screen
	ATTRIBUTES (read_HDF5_cryosat_L2): input variable attributes from HDF5 file

PYTHON DEPENDENCIES:
	numpy: Scientific Computing Tools For Python (http://www.numpy.org)
	h5py: Python interface for Hierarchal Data Format 5 (HDF5)
		(http://h5py.org)

UPDATE HISTORY:
	Updated 04/2019: print HDF5 keys from list for python3 compatibility
	Updated 06/2018: use items instead of iteritems for python3 compatibility
	Updated 11/2016: added Abs_Orbit and Ascending_Flg to Data_1Hz outputs
	Updated 05/2016: using __future__ print function
	Written 03/2016
"""
from __future__ import print_function

import os
import re
import h5py

#-- PURPOSE: write CryoSat-2 HDF5 files
def HDF5_cryosat_L2(CS_l2_mds, BASELINE, FILENAME='', TITLE='', HEADER=0,
	CLOBBER='Y', VERBOSE='N'):
	#-- setting HDF5 clobber attribute
	if CLOBBER in ('Y','y'):
		clobber = 'w'
	else:
		clobber = 'w-'

	#-- getting HDF5 dataset attributes for each variable
	CS_l2_attrib = cryosat_L2_attributes(BASELINE)

	#-- open output HDF5 file
	fileID = h5py.File(os.path.expanduser(FILENAME), clobber)
	#-- create sub-groups within HDF5 file
	fileID.create_group('Data_1Hz')
	fileID.create_group('Data_20Hz')
	fileID.create_group('Corrections')

	#-- Dimensions of parameters
	n_records,n_blocks = CS_l2_mds['Data_20Hz']['D_time_mics'].shape

	#-- find keys to output (do not output empty Spares variables)
	Data_1Hz_keys = [key for key in CS_l2_mds['Data_1Hz'].keys() if not
		re.search('Spare',key)]
	Data_20Hz_keys = [key for key in CS_l2_mds['Data_20Hz'].keys() if not
		re.search('Spare',key)]
	correction_keys = [key for key in CS_l2_mds['Corrections'].keys() if not
		re.search('Spare',key)]

	#-- create HDF5 records
	h5 = {}
	h5['Data_1Hz'] = {}
	h5['Data_20Hz'] = {}
	h5['Corrections'] = {}
	#-- CryoSat-2 1 Hz data fields (Location Group)
	#-- Time and Orbit Parameters plus Measurement Mode
	for key in Data_1Hz_keys:
		val = CS_l2_mds['Data_1Hz'][key]
		#-- Defining the HDF5 dataset variables
		h5['Data_1Hz'][key] = fileID.create_dataset('Data_1Hz/{0}'.format(key),
			(n_records,), data=val, dtype=val.dtype, compression='gzip')
		#-- add HDF5 variable attributes
		for att_name,att_val in CS_l2_attrib['Data_1Hz'][key].items():
			h5['Data_1Hz'][key].attrs[att_name] = att_val
		#-- attach dimensions
		h5['Data_1Hz'][key].dims[0].label='CS_L2_MDS_REC_SIZE'
	#-- CryoSat-2 20 Hz data fields (Measurement Group)
	#-- Derived from instrument measurement parameters
	for key in Data_20Hz_keys:
		val = CS_l2_mds['Data_20Hz'][key]
		#-- Defining the HDF5 dataset variables
		h5['Data_20Hz'][key] = fileID.create_dataset('Data_20Hz/{0}'.format(key),
			(n_records,n_blocks,), data=val, dtype=val.dtype, compression='gzip')
		#-- add HDF5 variable attributes
		for att_name,att_val in CS_l2_attrib['Data_20Hz'][key].items():
			h5['Data_20Hz'][key].attrs[att_name] = att_val
		#-- attach dimensions
		h5['Data_20Hz'][key].dims[0].label='CS_L2_MDS_REC_SIZE'
		h5['Data_20Hz'][key].dims[1].label='CS_L2_MDS_BLOCK_SIZE'
	#-- CryoSat-2 geophysical corrections (External Corrections Group)
	for key in correction_keys:
		val = CS_l2_mds['Corrections'][key]
		#-- Defining the HDF5 dataset variables
		h5['Corrections'][key] = fileID.create_dataset('Corrections/{0}'.format(key),
			(n_records,), data=val, dtype=val.dtype, compression='gzip')
		#-- add HDF5 variable attributes
		for att_name,att_val in CS_l2_attrib['Corrections'][key].items():
			h5['Corrections'][key].attrs[att_name] = att_val
		#-- attach dimensions
		h5['Corrections'][key].dims[0].label='CS_L2_MDS_REC_SIZE'

	#-- output MPH/SPH/DSD headers as group attributes
	if (HEADER == 1):
		#-- HEADER 1 is for single CryoSat-2 files
		fileID.create_group('METADATA')
		fileID['METADATA'].create_group('MPH')
		fileID['METADATA'].create_group('SPH')
		fileID['METADATA'].create_group('DSD')
		#-- Main Product Header (MPH) are all strings
		for att_name,att_val in CS_l2_mds['METADATA']['MPH'].items():
			fileID['METADATA']['MPH'].attrs[att_name] = att_val
		#-- Specific Product Header (SPH) are both strings and dictionaries
		for att_name,att_val in CS_l2_mds['METADATA']['SPH'].items():
			if isinstance(att_val,dict):
				#-- if att_val is dictionary
				fileID['METADATA']['SPH'].create_group(att_name)
				for ds_name,ds_val in att_val.items():
					fileID['METADATA']['SPH'][att_name].attrs[ds_name] = ds_val
			elif isinstance(att_val,basestring) and att_name:
				#-- if att_val is string
				fileID['METADATA']['SPH'].attrs[att_name] = att_val
		#-- Data Set Descriptors (DSD) are all strings
		for att_name,att_val in CS_l2_mds['METADATA']['DSD'].items():
			fileID['METADATA']['DSD'].attrs[att_name] = att_val
	elif (HEADER == 2):
		#-- HEADER 2 is for merged CryoSat-2 files from convert_cryosat_L2.py
		fileID.create_group('METADATA')
		fileID['METADATA'].create_group('MPH')
		fileID['METADATA'].create_group('SPH')
		fileID['METADATA'].create_group('DSD')
		#-- Main Product Header (MPH) are all strings
		for fi in CS_l2_mds['METADATA']['MPH'].keys():
			fileID['METADATA']['MPH'].create_group(fi)
			for att_name,att_val in CS_l2_mds['METADATA']['MPH'][fi].items():
				fileID['METADATA']['MPH'][fi].attrs[att_name] = att_val
		#-- Specific Product Header (SPH) are both strings and dictionaries
		for fi in CS_l2_mds['METADATA']['SPH'].keys():
			fileID['METADATA']['SPH'].create_group(fi)
			for att_name,att_val in CS_l2_mds['METADATA']['SPH'][fi].items():
				if isinstance(att_val,dict):
					#-- if att_val is dictionary
					fileID['METADATA']['SPH'][fi].create_group(att_name)
					for dsn,dsv in att_val.items():
						fileID['METADATA']['SPH'][fi][att_name].attrs[dsn] = dsv
				elif isinstance(att_val,basestring) and att_name:
					#-- if att_val is string
					fileID['METADATA']['SPH'][fi].attrs[att_name] = att_val
		#-- Data Set Descriptors (DSD) are all strings
		for fi in CS_l2_mds['METADATA']['DSD'].keys():
			fileID['METADATA']['DSD'].create_group(fi)
			for att_name,att_val in CS_l2_mds['METADATA']['DSD'][fi].items():
				fileID['METADATA']['DSD'][fi].attrs[att_name] = att_val

	#-- output file title
	fileID.attrs['description'] = TITLE

	#-- Output HDF5 structure information
	if VERBOSE in ('Y','y'):
		print(FILENAME)
		print(list(fileID.keys()))

	#-- Closing the HDF5 file
	fileID.close()

#-- PURPOSE: read CryoSat-2 HDF5 files
def read_HDF5_cryosat_L2(FILENAME, ATTRIBUTES='Y', VERBOSE='N'):
	#-- Open the HDF5 file for reading
	fileID = h5py.File(os.path.expanduser(FILENAME), 'r')

	#-- Output HDF5 file information
	if VERBOSE in ('Y','y'):
		print(fileID.filename)
		print(list(fileID.keys()))

	#-- allocate python dictionaries for output CS_l2_mds variables
	CS_l2_mds = {}
	CS_l2_mds['Data_1Hz'] = {}
	CS_l2_mds['Data_20Hz'] = {}
	CS_l2_mds['Corrections'] = {}
	#-- get each HDF5 variable
	#-- CryoSat-2 Location Group
	for key in fileID['Data_1Hz'].keys():
		CS_l2_mds['Data_1Hz'][key] = fileID['Data_1Hz'][key][:]
	#-- CryoSat-2 Measurement Group
	for key in fileID['Data_20Hz'].keys():
		CS_l2_mds['Data_20Hz'][key] = fileID['Data_20Hz'][key][:,:]
	#-- CryoSat-2 External Corrections Group
	for key in fileID['Corrections'].keys():
		CS_l2_mds['Corrections'][key] = fileID['Corrections'][key][:]

	#-- Getting attributes of included variables
	if ATTRIBUTES in ('Y','y'):
		#-- allocate python dictionaries for output CS_l2_mds attributes
		CS_l2_mds['Attributes'] = {}
		CS_l2_mds['Attributes']['Data_1Hz'] = {}
		CS_l2_mds['Attributes']['Data_20Hz'] = {}
		CS_l2_mds['Attributes']['Corrections'] = {}
		#-- CryoSat-2 Location Group
		for key in fileID['Data_1Hz'].keys():
			CS_l2_mds['Attributes']['Data_1Hz'][key] = {}
			for att_name,att_val in fileID['Data_1Hz'][key].attrs.items():
				CS_l2_mds['Attributes']['Data_1Hz'][key][att_name] = att_val
		#-- CryoSat-2 Measurement Group
		for key in fileID['Data_20Hz'].keys():
			CS_l2_mds['Attributes']['Data_20Hz'][key] = {}
			for att_name,att_val in fileID['Data_20Hz'][key].attrs.items():
				CS_l2_mds['Attributes']['Data_20Hz'][key][att_name] = att_val
		#-- CryoSat-2 External Corrections Group
		for key in fileID['Corrections'].keys():
			CS_l2_mds['Attributes']['Corrections'][key] = {}
			for att_name,att_val in fileID['Corrections'][key].attrs.items():
				CS_l2_mds['Attributes']['Corrections'][key][att_name] = att_val
		#-- Global attribute description
		CS_l2_mds['Attributes']['title'] = fileID.attrs['description']

	#-- Closing the HDF5 file
	fileID.close()
	return CS_l2_mds

#-- PURPOSE: get the number of records and number of blocks in an HDF5 file
def HDF5_cryosat_L2_shape(FILENAME):
	#-- Open the HDF5 file for reading
	fileID = h5py.File(os.path.expanduser(FILENAME), 'r')
	n_records,n_blocks = fileID['Data_20Hz']['D_time_mics'].shape
	#-- Closing the HDF5 file
	fileID.close()
	return (n_records,n_blocks)

#-- PURPOSE: get attribute names for baseline
def cryosat_L2_attributes(BASELINE):
	#-- CryoSat-2 1 Hz data fields (Location Group)
	#-- Time and Orbit Parameters plus Measurement Mode
	L2_1Hz_attributes = {}
	#-- Time: day part
	L2_1Hz_attributes['Day'] = {}
	L2_1Hz_attributes['Day']['long_name'] = 'MDSR time stamp days'
	L2_1Hz_attributes['Day']['units'] = 'days since 2000-01-01 00:00:00 TAI'
	L2_1Hz_attributes['Day']['hertz'] = 1
	#-- Time: second part
	L2_1Hz_attributes['Second'] = {}
	L2_1Hz_attributes['Second']['long_name'] = 'MDSR time stamp seconds'
	L2_1Hz_attributes['Second']['units'] = 'seconds'
	L2_1Hz_attributes['Second']['hertz'] = 1
	#-- Time: microsecond part
	L2_1Hz_attributes['Micsec'] = {}
	L2_1Hz_attributes['Micsec']['long_name'] = 'MDSR time stamp microseconds'
	L2_1Hz_attributes['Micsec']['units'] = 'microseconds'
	L2_1Hz_attributes['Micsec']['hertz'] = 1
	#-- SIRAL mode
	L2_1Hz_attributes['Siral_mode'] = {}
	L2_1Hz_attributes['Siral_mode']['long_name'] = 'SIRAL mode'
	L2_1Hz_attributes['Siral_mode']['hertz'] = 1
	#-- Lat_1Hz: packed units (0.1 micro-degree, 1e-7 degrees)
	L2_1Hz_attributes['Lat_1Hz'] = {}
	L2_1Hz_attributes['Lat_1Hz']['long_name'] = 'Latitude of measurement'
	L2_1Hz_attributes['Lat_1Hz']['description'] = ('Corresponding to the nadir '
		'position at the time of the 1Hz time stamp')
	L2_1Hz_attributes['Lat_1Hz']['units'] = '0.1 micro-degree'
	L2_1Hz_attributes['Lat_1Hz']['valid_min'] = -9e8
	L2_1Hz_attributes['Lat_1Hz']['valid_max'] = 9e8
	L2_1Hz_attributes['Lat_1Hz']['hertz'] = 1
	#-- Lon_1Hz: packed units (0.1 micro-degree, 1e-7 degrees)
	L2_1Hz_attributes['Lon_1Hz'] = {}
	L2_1Hz_attributes['Lon_1Hz']['long_name'] = 'Longitude of measurement'
	L2_1Hz_attributes['Lon_1Hz']['description'] = ('Corresponding to the nadir '
		'position at the time of the 1Hz time stamp')
	L2_1Hz_attributes['Lon_1Hz']['units'] = '0.1 micro-degree'
	L2_1Hz_attributes['Lon_1Hz']['valid_min'] = -18e8
	L2_1Hz_attributes['Lon_1Hz']['valid_max'] = 18e8
	L2_1Hz_attributes['Lon_1Hz']['hertz'] = 1
	#-- Alt_1Hz: packed units (mm, 1e-3 m)
	#-- Altitude of COG above reference ellipsoid (interpolated value)
	L2_1Hz_attributes['Alt_1Hz'] = {}
	L2_1Hz_attributes['Alt_1Hz']['long_name'] = 'Altitude'
	L2_1Hz_attributes['Alt_1Hz']['description'] = ('Altitude of Satellite COG '
		'above reference ellipsoid at nadir (interpolated value)')
	L2_1Hz_attributes['Alt_1Hz']['units'] = 'millimeters'
	L2_1Hz_attributes['Alt_1Hz']['hertz'] = 1

	#-- Spacecraft mispointing (Roll, Pitch, and Yaw for Baseline-C)
	if (BASELINE == 'C'):
		#-- Roll: packed units (0.1 micro-degree, 1e-7 degrees)
		L2_1Hz_attributes['Roll'] = {}
		L2_1Hz_attributes['Roll']['long_name'] = ('Spacecraft roll angle '
			'derived from star trackers')
		L2_1Hz_attributes['Roll']['units'] = '0.1 micro-degree'
		L2_1Hz_attributes['Roll']['hertz'] = 1
		#-- Pitch: packed units (0.1 micro-degree, 1e-7 degrees)
		L2_1Hz_attributes['Pitch'] = {}
		L2_1Hz_attributes['Pitch']['long_name'] = ('Spacecraft pitch angle '
			'derived from star trackers')
		L2_1Hz_attributes['Pitch']['units'] = '0.1 micro-degree'
		L2_1Hz_attributes['Pitch']['hertz'] = 1
		#-- Yaw: packed units (0.1 micro-degree, 1e-7 degrees)
		L2_1Hz_attributes['Yaw'] = {}
		L2_1Hz_attributes['Yaw']['long_name'] = ('Spacecraft yaw angle '
			'derived from star trackers')
		L2_1Hz_attributes['Yaw']['units'] = '0.1 micro-degree'
		L2_1Hz_attributes['Yaw']['hertz'] = 1
	else:
		#-- Mispointing: packed units (millidegrees, 1e-3 degrees)
		L2_1Hz_attributes['Mispointing'] = {}
		L2_1Hz_attributes['Mispointing']['long_name'] =('Spacecraft mispointing'
			'derived from star trackers')
		L2_1Hz_attributes['Mispointing']['units'] = 'millidegrees'
		L2_1Hz_attributes['Mispointing']['hertz'] = 1

	#-- Number of valid records in the block of twenty that contain data
	#-- Last few records of the last block of a datset may be blank blocks
	#-- inserted to bring the file up to a multiple of twenty.
	L2_1Hz_attributes['N_valid'] = {}
	L2_1Hz_attributes['N_valid']['long_name'] = 'Number of valid measurements'
	L2_1Hz_attributes['N_valid']['description'] = ('The number of records in '
		'the block of twenty that actually contain data. The last few records '
		'of the last block of a datset may be blank blocks inserted to bring '
		'the file up to a multiple of twenty.')
	L2_1Hz_attributes['N_valid']['hertz'] = 1

	#-- Absolute Orbit Number from MPH Header
	L2_1Hz_attributes['Abs_Orbit'] = {}
	L2_1Hz_attributes['Abs_Orbit']['long_name'] = 'Absolute Orbit Number'
	L2_1Hz_attributes['Abs_Orbit']['units'] = 'Count'
	L2_1Hz_attributes['Abs_Orbit']['hertz'] = 1
	#-- Ascending Flag from SPH Header
	L2_1Hz_attributes['Ascending_Flg'] = {}
	L2_1Hz_attributes['Ascending_Flg']['long_name'] = 'Ascending Flag'
	L2_1Hz_attributes['Ascending_Flg']['description'] = ('If True: satellite '
		'is in ascending orbit. If False: satellite is in descending orbit.')
	L2_1Hz_attributes['Ascending_Flg']['hertz'] = 1

	#-- CryoSat-2 geophysical corrections (External Corrections Group)
	L2_corr_attributes = {}
	#-- Dry Tropospheric Correction packed units (mm, 1e-3 m)
	L2_corr_attributes['dryTrop'] = {}
	L2_corr_attributes['dryTrop']['long_name'] = 'Dry Tropospheric Correction'
	L2_corr_attributes['dryTrop']['description'] = ('Altimeter range correction'
		' due to the dry-gas component of the troposphere')
	L2_corr_attributes['dryTrop']['units'] = 'millimeters'
	L2_corr_attributes['dryTrop']['hertz'] = 1
	#-- Wet Tropospheric Correction packed units (mm, 1e-3 m)
	L2_corr_attributes['wetTrop'] = {}
	L2_corr_attributes['wetTrop']['long_name'] = 'Wet Tropospheric Correction'
	L2_corr_attributes['wetTrop']['description'] = ('Altimeter range correction'
		' due to the water component of the troposphere')
	L2_corr_attributes['wetTrop']['units'] = 'millimeters'
	L2_corr_attributes['wetTrop']['hertz'] = 1
	#-- Inverse Barometric Correction packed units (mm, 1e-3 m)
	L2_corr_attributes['InvBar'] = {}
	L2_corr_attributes['InvBar']['long_name'] = 'Inverse Barometric Correction'
	L2_corr_attributes['InvBar']['description'] = ('Altimeter range correction '
		'for the depression of the ocean surface caused by the local barometric '
		'pressure')
	L2_corr_attributes['InvBar']['units'] = 'millimeters'
	L2_corr_attributes['InvBar']['hertz'] = 1
	#-- Dynamic Atmosphere Correction packed units (mm, 1e-3 m)
	L2_corr_attributes['DynAtm'] = {}
	L2_corr_attributes['DynAtm']['long_name'] = 'Dynamic Atmosphere Correction'
	L2_corr_attributes['DynAtm']['description'] = ('Altimeter range correction '
		'for the dynamic component of the wind effect on the ocean '
		'(applicable only on ocean surface without sea-ice)')
	L2_corr_attributes['DynAtm']['units'] = 'millimeters'
	L2_corr_attributes['DynAtm']['hertz'] = 1
	#-- Ionospheric Correction packed units (mm, 1e-3 m)
	L2_corr_attributes['Iono'] = {}
	L2_corr_attributes['Iono']['long_name'] = 'Ionospheric Correction'
	L2_corr_attributes['Iono']['description'] = ('Altimeter range correction '
		'for the delay of the radar pulse caused by free electrons in the '
		'ionosphere. (Computed from a simple model in NRT data or from GPS '
		'satellite derived (GIM) map in normal processing.)')
	L2_corr_attributes['Iono']['units'] = 'millimeters'
	L2_corr_attributes['Iono']['hertz'] = 1
	#-- Sea State Bias Correction packed units (mm, 1e-3 m)
	L2_corr_attributes['SSB'] = {}
	L2_corr_attributes['SSB']['long_name'] = 'Sea State Bias Correction'
	L2_corr_attributes['SSB']['description'] = ('Empirical altimeter range '
		'correction proportional to the significant wave height which '
		'compensates for the asymmetric shape of ocean waves. '
		'(computed by the geophysical CFI library)')
	L2_corr_attributes['SSB']['units'] = 'millimeters'
	L2_corr_attributes['SSB']['hertz'] = 1
	#-- Ocean tide Correction packed units (mm, 1e-3 m)
	L2_corr_attributes['ocTideElv'] = {}
	L2_corr_attributes['ocTideElv']['long_name'] = 'Ocean Tide'
	L2_corr_attributes['ocTideElv']['description'] = ('Removes the effect of '
		'local tide and adjusts the measurement to the mean sea surface')
	L2_corr_attributes['ocTideElv']['units'] = 'millimeters'
	L2_corr_attributes['ocTideElv']['hertz'] = 1
	#-- Long period equilibrium ocean tide Correction packed units (mm, 1e-3 m)
	L2_corr_attributes['lpeTideElv'] = {}
	L2_corr_attributes['lpeTideElv']['long_name'] = ('Long-Period Equilibrium '
		'Ocean Tide')
	L2_corr_attributes['lpeTideElv']['description'] = ('Removes the effect of '
		'local tide and adjusts the measurement to the mean sea surface')
	L2_corr_attributes['lpeTideElv']['units'] = 'millimeters'
	L2_corr_attributes['lpeTideElv']['hertz'] = 1
	#-- Ocean loading tide Correction packed units (mm, 1e-3 m)
	L2_corr_attributes['olTideElv'] = {}
	L2_corr_attributes['olTideElv']['long_name'] = 'Ocean Loading Tide'
	L2_corr_attributes['olTideElv']['description'] = ('Removes the effect of '
		'local tidal distortion of the Earth crust')
	L2_corr_attributes['olTideElv']['units'] = 'millimeters'
	L2_corr_attributes['olTideElv']['hertz'] = 1
	#-- Solid Earth tide Correction packed units (mm, 1e-3 m)
	L2_corr_attributes['seTideElv'] = {}
	L2_corr_attributes['seTideElv']['long_name'] = 'Solid Earth Tide'
	L2_corr_attributes['seTideElv']['description'] = ('Removes the effect of '
		'local tidal distortion in the Earth crust')
	L2_corr_attributes['seTideElv']['units'] = 'millimeters'
	L2_corr_attributes['seTideElv']['hertz'] = 1
	#-- Geocentric Polar tide Correction packed units (mm, 1e-3 m)
	L2_corr_attributes['gpTideElv'] = {}
	L2_corr_attributes['gpTideElv']['long_name'] = 'Geocentric Polar Tide'
	L2_corr_attributes['gpTideElv']['description'] = ('Removes a long-period '
		'distortion of the Earth crust caused by variations in the centrifugal '
		'force as the Earth rotational axis moves its geographic location')
	L2_corr_attributes['gpTideElv']['units'] = 'millimeters'
	L2_corr_attributes['gpTideElv']['hertz'] = 1
	#-- Surface Type: Packed in groups of three bits for each of the 20 records
	L2_corr_attributes['Surf_type'] = {}
	L2_corr_attributes['Surf_type']['long_name'] = 'Surface Type Flag'
	L2_corr_attributes['Surf_type']['description'] = ('Enumerated key to '
		'classify surface at nadir provided by a model: (0=Open Ocean, '
		'1=Closed Sea, 2=Continental Ice, 3=Land, 4-7=currently unused)')
	L2_corr_attributes['Surf_type']['hertz'] = 1
	#-- Mean Sea Surface or Geoid packed units (mm, 1e-3 m)
	L2_corr_attributes['MSS_Geoid'] = {}
	L2_corr_attributes['MSS_Geoid']['long_name'] = 'Mean Sea Surface or Geoid'
	L2_corr_attributes['MSS_Geoid']['description'] = ('Over Ocean (Surface'
		'Types 0 and 1): surface height from the Mean Sea Surface. Over Land '
		'Surface Types 2 and 3): geoid height from the CFI library.')
	L2_corr_attributes['MSS_Geoid']['units'] = 'millimeters'
	L2_corr_attributes['MSS_Geoid']['hertz'] = 1
	#-- Ocean Depth/Land Elevation Model (ODLE) packed units (mm, 1e-3 m)
	L2_corr_attributes['ODLE'] = {}
	L2_corr_attributes['ODLE']['long_name'] = 'ODLE from model'
	L2_corr_attributes['ODLE']['description'] = ('Ocean Depth / Land Elevation'
		'model supplied in the Geophysical Corrections CFI library')
	L2_corr_attributes['ODLE']['units'] = 'millimeters'
	L2_corr_attributes['ODLE']['hertz'] = 1
	#-- Ice Concentration packed units (%/100)
	L2_corr_attributes['Ice_conc'] = {}
	L2_corr_attributes['Ice_conc']['long_name'] = 'Ice Concentration'
	L2_corr_attributes['Ice_conc']['description'] = ('Obtained from a dynamic '
		'auxiliary file if data is available for the current time period. '
		'If data is not available, a climatology model is used.')
	L2_corr_attributes['Ice_conc']['units'] = '1/100 of a percent'
	L2_corr_attributes['Ice_conc']['hertz'] = 1
	#-- Snow Depth packed units (mm, 1e-3 m)
	L2_corr_attributes['Snow_depth'] = {}
	L2_corr_attributes['Snow_depth']['long_name'] = 'Snow Depth'
	L2_corr_attributes['Snow_depth']['description'] = ('Obtained from'
		'climatology model data. Can be used to adjust the freeboard estimate '
		'to account for snow-loading.')
	L2_corr_attributes['Snow_depth']['units'] = 'millimeters'
	L2_corr_attributes['Snow_depth']['hertz'] = 1
	#-- Snow Density packed units (kg/m^3)
	L2_corr_attributes['Snow_density'] = {}
	L2_corr_attributes['Snow_density']['long_name'] = 'Snow Density'
	L2_corr_attributes['Snow_density']['description'] = ('Obtained from'
		'climatology model data. Can be used to adjust the freeboard estimate '
		'to account for snow-loading.')
	L2_corr_attributes['Snow_density']['units'] = 'kg/m^3'
	L2_corr_attributes['Snow_density']['hertz'] = 1
	#-- Corrections Status Flag
	L2_corr_attributes['C_status'] = {}
	L2_corr_attributes['C_status']['long_name'] = 'Corrections Status Flag'
	L2_corr_attributes['C_status']['description'] = ('Used to show validity of '
		'the 1Hz corrections.  See table 2.3.3.1-4 of the "L2 Products Format '
		'Specification" document')
	L2_corr_attributes['C_status']['hertz'] = 1
	#-- Significant Wave Height (SWH) packed units (mm, 1e-3)
	L2_corr_attributes['SWH'] = {}
	L2_corr_attributes['SWH']['long_name'] = 'Significant Wave Height (SWH)'
	L2_corr_attributes['SWH']['units'] = 'millimeters'
	L2_corr_attributes['SWH']['hertz'] = 1
	#-- Wind Speed packed units (mm/s, 1e-3 m/s)
	L2_corr_attributes['Wind_speed'] = {}
	L2_corr_attributes['Wind_speed']['long_name'] = 'Altimetric wind speed'
	L2_corr_attributes['Wind_speed']['description'] = ('Calculated using a '
		'model by the geocorrections CFI')
	L2_corr_attributes['Wind_speed']['units'] = 'mm/second'
	L2_corr_attributes['Wind_speed']['hertz'] = 1

	#-- CryoSat-2 20 Hz data fields (Measurement Group)
	#-- Derived from instrument measurement parameters
	L2_20Hz_attributes = {}
	#-- Delta between the timestamps for 20Hz record and the 1Hz record
	#-- D_time_mics packed units (microseconds)
	L2_20Hz_attributes['D_time_mics'] = {}
	L2_20Hz_attributes['D_time_mics']['long_name'] = 'Delta Time'
	L2_20Hz_attributes['D_time_mics']['description'] = ('Adds to the 1 Hz time '
		'stamp to give the correct time for each of the 20Hz measurement '
		'(Set to zero if block is partially empty)')
	L2_20Hz_attributes['D_time_mics']['units'] = 'microseconds'
	L2_20Hz_attributes['D_time_mics']['hertz'] = 20
	#-- Lat: packed units (0.1 micro-degree, 1e-7 degrees)
	L2_20Hz_attributes['Lat'] = {}
	L2_20Hz_attributes['Lat']['long_name'] = 'Latitude of measurement'
	L2_20Hz_attributes['Lat']['description'] = ('Measurement Latitude of the '
		'echoing point position.  Includes the x-track offset for SIN '
		'measurements and the slope-corrected position for LRM measurements.')
	L2_20Hz_attributes['Lat']['units'] = '0.1 micro-degree'
	L2_20Hz_attributes['Lat']['valid_min'] = -9e8
	L2_20Hz_attributes['Lat']['valid_max'] = 9e8
	L2_20Hz_attributes['Lat']['hertz'] = 20
	#-- Lon: packed units (0.1 micro-degree, 1e-7 degrees)
	L2_20Hz_attributes['Lon'] = {}
	L2_20Hz_attributes['Lon']['long_name'] = 'Longitude of measurement'
	L2_20Hz_attributes['Lon']['description'] = ('Measurement Longitude of the '
		'echoing point position.  Includes the x-track offset for SIN '
		'measurements and the slope-corrected position for LRM measurements.')
	L2_20Hz_attributes['Lon']['units'] = '0.1 micro-degree'
	L2_20Hz_attributes['Lon']['valid_min'] = -18e8
	L2_20Hz_attributes['Lon']['valid_max'] = 18e8
	L2_20Hz_attributes['Lon']['hertz'] = 20

	#-- Measured elevation and Backscatter (for 3 retrackers with Baseline-C)
	if (BASELINE == 'C'):
		#-- Measured elevation above ellipsoid from retracker 1
		#-- packed units (mm, 1e-3 m)
		L2_20Hz_attributes['Elev_1'] = {}
		L2_20Hz_attributes['Elev_1']['long_name'] = 'Surface Elevation'
		L2_20Hz_attributes['Elev_1']['description'] = ('Height of surface at '
			'measurement point with respect to the reference ellipsoid (WGS84). '
			'It is calculated with the retracked and geocorrected range')
		L2_20Hz_attributes['Elev_1']['units'] = 'millimeters'
		L2_20Hz_attributes['Elev_1']['retracker'] = 1
		L2_20Hz_attributes['Elev_1']['hertz'] = 20
		#-- Measured elevation above ellipsoid from retracker 2
		#-- packed units (mm, 1e-3 m)
		L2_20Hz_attributes['Elev_2'] = {}
		L2_20Hz_attributes['Elev_2']['long_name'] = 'Surface Elevation'
		L2_20Hz_attributes['Elev_2']['description'] = ('Height of surface at '
			'measurement point with respect to the reference ellipsoid (WGS84). '
			'It is calculated with the retracked and geocorrected range')
		L2_20Hz_attributes['Elev_2']['units'] = 'millimeters'
		L2_20Hz_attributes['Elev_2']['retracker'] = 2
		L2_20Hz_attributes['Elev_2']['hertz'] = 20
		#-- Measured elevation above ellipsoid from retracker 3
		#-- packed units (mm, 1e-3 m)
		L2_20Hz_attributes['Elev_3'] = {}
		L2_20Hz_attributes['Elev_3']['long_name'] = 'Surface Elevation'
		L2_20Hz_attributes['Elev_3']['description'] = ('Height of surface at '
			'measurement point with respect to the reference ellipsoid (WGS84). '
			'It is calculated with the retracked and geocorrected range')
		L2_20Hz_attributes['Elev_3']['units'] = 'millimeters'
		L2_20Hz_attributes['Elev_3']['retracker'] = 3
		L2_20Hz_attributes['Elev_3']['hertz'] = 20
		#-- Sigma Zero Backscatter for retracker 1: packed units (1e-2 dB)
		L2_20Hz_attributes['Sig0_1'] = {}
		L2_20Hz_attributes['Sig0_1']['long_name'] = 'Sigma Zero Backscatter'
		L2_20Hz_attributes['Sig0_1']['description'] = ('Fully corrected '
			'including instrument gain corrections and retracker correction')
		L2_20Hz_attributes['Sig0_1']['units'] = 'dB/100'
		L2_20Hz_attributes['Sig0_1']['retracker'] = 1
		L2_20Hz_attributes['Sig0_1']['hertz'] = 20
		#-- Sigma Zero Backscatter for retracker 2: packed units (1e-2 dB)
		L2_20Hz_attributes['Sig0_2'] = {}
		L2_20Hz_attributes['Sig0_2']['long_name'] = 'Sigma Zero Backscatter'
		L2_20Hz_attributes['Sig0_2']['description'] = ('Fully corrected '
			'including instrument gain corrections and retracker correction')
		L2_20Hz_attributes['Sig0_2']['units'] = 'dB/100'
		L2_20Hz_attributes['Sig0_2']['retracker'] = 2
		L2_20Hz_attributes['Sig0_2']['hertz'] = 20
		#-- Sigma Zero Backscatter for retracker 3: packed units (1e-2 dB)
		L2_20Hz_attributes['Sig0_3'] = {}
		L2_20Hz_attributes['Sig0_3']['long_name'] = 'Sigma Zero Backscatter'
		L2_20Hz_attributes['Sig0_3']['description'] = ('Fully corrected '
			'including instrument gain corrections and retracker correction')
		L2_20Hz_attributes['Sig0_3']['units'] = 'dB/100'
		L2_20Hz_attributes['Sig0_3']['retracker'] = 3
		L2_20Hz_attributes['Sig0_3']['hertz'] = 20
	else:
		#-- Measured elevation above ellipsoid from retracker
		#-- packed units (mm, 1e-3 m)
		L2_20Hz_attributes['Elev'] = {}
		L2_20Hz_attributes['Elev']['long_name'] = 'Surface Elevation'
		L2_20Hz_attributes['Elev']['description'] = ('Height of surface at '
			'measurement point with respect to the reference ellipsoid (WGS84). '
			'It is calculated with the retracked and geocorrected range')
		L2_20Hz_attributes['Elev']['units'] = 'millimeters'
		L2_20Hz_attributes['Elev']['hertz'] = 20
		#-- Sigma Zero Backscatter for retracker: packed units (1e-2 dB)
		L2_20Hz_attributes['Sig0'] = {}
		L2_20Hz_attributes['Sig0']['long_name'] = 'Sigma Zero Backscatter'
		L2_20Hz_attributes['Sig0']['description'] = ('Fully corrected '
			'including instrument gain corrections and retracker correction')
		L2_20Hz_attributes['Sig0']['units'] = 'dB/100'
		L2_20Hz_attributes['Sig0']['hertz'] = 20

	#-- Freeboard: packed units (mm, 1e-3 m)
	#-- -9999 default value indicates computation has not been performed
	L2_20Hz_attributes['Freeboard'] = {}
	L2_20Hz_attributes['Freeboard']['long_name'] = 'Freeboard'
	L2_20Hz_attributes['Freeboard']['description'] = ('SAR mode computed '
		'freeboard of the Sea Ice. Freeboard can be a small negative value '
		'when there is sufficient snow loading on thin ice. Set to 0 in SIN '
		'and LRM modes.')
	L2_20Hz_attributes['Freeboard']['units'] = 'millimeters'
	L2_20Hz_attributes['Freeboard']['_FillValue'] = -9999
	L2_20Hz_attributes['Freeboard']['hertz'] = 20
	#-- Interpolated Sea Surface Height Anomaly: packed units (mm, 1e-3 m)
	L2_20Hz_attributes['SSHA_interp'] = {}
	L2_20Hz_attributes['SSHA_interp']['long_name'] = ('Interpolated Sea '
		'Surface Height Anomaly')
	L2_20Hz_attributes['SSHA_interp']['description'] = ('Ocean height anomaly '
		'defined by comparing the interpolated ocean height from the SAR '
		'processing with the MSS from the model')
	L2_20Hz_attributes['SSHA_interp']['units'] = 'millimeters'
	L2_20Hz_attributes['SSHA_interp']['hertz'] = 20
	#-- Interpolated Sea Surface Height measurement count
	L2_20Hz_attributes['SSHA_num'] = {}
	L2_20Hz_attributes['SSHA_num']['long_name'] = ('Interpolated Sea '
		'Surface Height measurement count')
	L2_20Hz_attributes['SSHA_num']['description'] = ('Number of records used '
		'to create the fit in the SSHA calculation')
	L2_20Hz_attributes['SSHA_num']['hertz'] = 20
	#-- Interpolation quality estimate RSS: packed units (mm, 1e-3 m)
	L2_20Hz_attributes['SSHA_qual'] = {}
	L2_20Hz_attributes['SSHA_qual']['long_name'] = ('Interpolated Sea '
		'Surface Height quality estimate')
	L2_20Hz_attributes['SSHA_qual']['description'] = ('Root mean square (RMS)'
		'of the residuals of the SSHA fit')
	L2_20Hz_attributes['SSHA_qual']['units'] = 'millimeters'
	L2_20Hz_attributes['SSHA_qual']['hertz'] = 20
	#-- Peakiness: packed units (1e-2)
	L2_20Hz_attributes['Peakiness'] = {}
	L2_20Hz_attributes['Peakiness']['long_name'] = 'Peakiness'
	L2_20Hz_attributes['Peakiness']['description'] = ('Peakiness of the echo '
		'in the L1b product. Requires different interpretation for SAR and SIN '
		'echoes which do not have the usual pulse-limited echo shape.')
	L2_20Hz_attributes['Peakiness']['units'] = '1/100'
	L2_20Hz_attributes['Peakiness']['hertz'] = 20
	#-- Number of averaged echoes or beams
	L2_20Hz_attributes['N_avg'] = {}
	L2_20Hz_attributes['N_avg']['long_name'] = 'Number of Echoes/Beams averaged'
	L2_20Hz_attributes['N_avg']['description'] = ('In LRM mode: number of '
		'echoes which have been averaged to make one measurement (normally). '
		'In SAR and SIN modes: number of Doppler beams which have been stacked '
		'to derive each measurement. Near the begining and end of each section '
		'of SAR or SIN mode operation, this number reduces below the nominal '
		'value and there is a corresponding decrease in the signal to noise '
		'ratio of the waveform.')
	L2_20Hz_attributes['N_avg']['hertz'] = 20
	#-- Quality flags
	L2_20Hz_attributes['Quality_Flg'] = {}
	L2_20Hz_attributes['Quality_Flg']['long_name'] = 'Measurement Quality Flags'
	L2_20Hz_attributes['Quality_Flg']['description'] = ('Used to show validity '
		'of the 20Hz measurements.  See table 2.3.3.1-5 of the "L2 Products '
		'Format Specification" document')
	L2_20Hz_attributes['Quality_Flg']['hertz'] = 20

	#-- Corrections Flag and Quality metrics for 3 retrackers with Baseline-C
	if (BASELINE == 'C'):
		#-- Corrections Application Flag
		L2_20Hz_attributes['Corrections_Flg'] = {}
		L2_20Hz_attributes['Corrections_Flg']['long_name'] = ('Corrections '
			'Application Flag')
		L2_20Hz_attributes['Corrections_Flg']['hertz'] = 20
		#-- Quality metric for retracker 1
		L2_20Hz_attributes['Quality_1'] = {}
		L2_20Hz_attributes['Quality_1']['long_name'] = ('Quality metric for '
			'retracker')
		L2_20Hz_attributes['Quality_1']['retracker'] = 1
		L2_20Hz_attributes['Quality_1']['hertz'] = 20
		#-- Quality metric for retracker 2
		L2_20Hz_attributes['Quality_2'] = {}
		L2_20Hz_attributes['Quality_2']['long_name'] = ('Quality metric for '
			'retracker')
		L2_20Hz_attributes['Quality_2']['retracker'] = 2
		L2_20Hz_attributes['Quality_2']['hertz'] = 20
		#-- Quality metric for retracker 3
		L2_20Hz_attributes['Quality_3'] = {}
		L2_20Hz_attributes['Quality_3']['long_name'] = ('Quality metric for '
			'retracker')
		L2_20Hz_attributes['Quality_3']['retracker'] = 3
		L2_20Hz_attributes['Quality_3']['hertz'] = 20

	#-- Bind all the l2 attributes together into single dictionary
	CS_l2_attrib = {}
	CS_l2_attrib['Data_1Hz'] = L2_1Hz_attributes
	CS_l2_attrib['Corrections'] = L2_corr_attributes
	CS_l2_attrib['Data_20Hz'] = L2_20Hz_attributes
	#-- return the output dictionary
	return CS_l2_attrib
