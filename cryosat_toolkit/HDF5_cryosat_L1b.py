#!/usr/bin/env python
u"""
HDF5_cryosat_L1b.py (08/2020)
Reads and Writes HDF5 files for CryoSat-2 Level-1b data products
Supported CryoSat Modes: LRM, SAR, SARin, FDM, SID, GDR

OUTPUTS a formatted HDF5 file with:
    Location: Time and Orbit Group
    Data: Measurements Group
    Geometry: External Corrections Group
    Waveform_1Hz: Average Waveforms Group
    Waveform_20Hz: Waveforms Group (with SAR/SARIN Beam Behavior Parameters)
    METADATA: MPH, SPH and DSD Header data

OPTIONS:
    BASELINE (HDF5_cryosat_L1b): CryoSat-2 baseline (A, B, C)
    FILENAME (HDF5_cryosat_L1b): output HDF5 file name
    TITLE (HDF5_cryosat_L1b): output file description
    HEADER (HDF5_cryosat_L1b): output CryoSat-2 file headers (MPH, SPH, DSD)
        1: for single CryoSat-2 files
        2: for merged CryoSat-2 files from convert_cryosat_L1b.py
    CLOBBER (HDF5_cryosat_L1b): overwrite existing HDF5 file
    VERBOSE: print HDF5 structure parameters to screen
    ATTRIBUTES (read_HDF5_cryosat_L1b): input variable attributes from HDF5 file

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        (https://www.h5py.org/)

UPDATE HISTORY:
    Updated 08/2020: flake8 updates for python3
    Updated 02/2020: convert from hard to soft tabulation
    Updated 10/2019: changing Y/N flags to True/False
    Updated 09/2019: updates for Baseline D
    Updated 04/2019: print HDF5 keys from list for python3 compatibility
    Updated 06/2018: use items instead of iteritems for python3 compatibility
    Updated 05/2016: using __future__ print function
    Updated 04/2016: fixed read attributes for Beam Behavior Parameters
    Written 03/2016
"""
from __future__ import print_function

import os
import re
import h5py

#-- PURPOSE: write CryoSat-2 HDF5 files
def HDF5_cryosat_L1b(CS_l1b_mds, MODE, BASELINE, FILENAME='', TITLE='',
    HEADER=0, CLOBBER=True, VERBOSE=False):
    #-- setting HDF5 clobber attribute
    if CLOBBER:
        clobber = 'w'
    else:
        clobber = 'w-'

    #-- getting HDF5 dataset attributes for each variable
    CS_l1b_attrib = cryosat_L1b_attributes(MODE, BASELINE)

    #-- open output HDF5 file
    fileID = h5py.File(os.path.expanduser(FILENAME), clobber)
    #-- create sub-groups within HDF5 file
    fileID.create_group('Location')
    fileID.create_group('Data')
    fileID.create_group('Geometry')
    fileID.create_group('Waveform_1Hz')
    fileID.create_group('Waveform_20Hz')
    #-- for SAR and SARIN modes: add subgroup for Beam Behavior Parameters
    #-- within Waveform_20Hz group
    if MODE in ('SAR','SIN'):
        fileID['Waveform_20Hz'].create_group('Beam')

    #-- Dimensions of parameters
    n_records,n_blocks = CS_l1b_mds['Location']['Day'].shape
    n_1Hz_wfm = CS_l1b_mds['Waveform_1Hz']['Waveform'].shape[1]
    n_20Hz_wfm = CS_l1b_mds['Waveform_20Hz']['Waveform'].shape[2]

    #-- find keys to output (do not output empty Spares variables)
    Location_keys = [key for key in CS_l1b_mds['Location'].keys() if not
        re.search('Spare',key)]
    Data_keys = [key for key in CS_l1b_mds['Data'].keys() if not
        re.search('Spare',key)]
    Geometry_keys = [key for key in CS_l1b_mds['Geometry'].keys() if not
        re.search('Spare',key)]
    Wfm_1Hz_keys = [key for key in CS_l1b_mds['Waveform_1Hz'].keys() if not
        re.search('Spare',key)]
    Wfm_20Hz_keys = [key for key in CS_l1b_mds['Waveform_20Hz'].keys() if not
        re.search('Spare',key)]
    if MODE in ('SAR','SIN'):
        beam_keys = [key for key in CS_l1b_mds['Waveform_20Hz']['Beam'].keys()
            if not re.search('Spare',key)]

    #-- create HDF5 records
    h5 = {}
    h5['Location'] = {}
    h5['Data'] = {}
    h5['Geometry'] = {}
    h5['Waveform_1Hz'] = {}
    h5['Waveform_20Hz'] = {}

    #-- CryoSat-2 Time and Orbit Group
    for key in Location_keys:
        val = CS_l1b_mds['Location'][key]
        if key in ('Sat_velocity','Real_beam','Baseline'):
            #-- Defining the HDF5 dataset variables
            h5['Location'][key] = fileID.create_dataset('Location/{0}'.format(key),
                (n_records,n_blocks,3,), data=val, dtype=val.dtype,
                compression='gzip')
            #-- attach dimensions
            h5['Location'][key].dims[0].label='CS_L1b_MDS_REC_SIZE'
            h5['Location'][key].dims[1].label='CS_L1b_MDS_BLOCK_SIZE'
            h5['Location'][key].dims[2].label='CS_L1b_MDS_VECTOR_SIZE'
        else:
            #-- Defining the HDF5 dataset variables
            h5['Location'][key] = fileID.create_dataset('Location/{0}'.format(key),
                (n_records,n_blocks,), data=val, dtype=val.dtype,
                compression='gzip')
            #-- attach dimensions
            h5['Location'][key].dims[0].label='CS_L1b_MDS_REC_SIZE'
            h5['Location'][key].dims[1].label='CS_L1b_MDS_BLOCK_SIZE'
        #-- add HDF5 variable attributes
        for att_name,att_val in CS_l1b_attrib['Location'][key].items():
            h5['Location'][key].attrs[att_name] = att_val

    #-- CryoSat-2 Measurement Group
    #-- Derived from instrument measurement parameters
    for key in Data_keys:
        val = CS_l1b_mds['Data'][key]
        #-- Defining the HDF5 dataset variables
        h5['Data'][key] = fileID.create_dataset('Data/{0}'.format(key),
            (n_records,n_blocks,), data=val, dtype=val.dtype, compression='gzip')
        #-- attach dimensions
        h5['Data'][key].dims[0].label='CS_L1b_MDS_REC_SIZE'
        h5['Data'][key].dims[1].label='CS_L1b_MDS_BLOCK_SIZE'
        #-- add HDF5 variable attributes
        for att_name,att_val in CS_l1b_attrib['Data'][key].items():
            h5['Data'][key].attrs[att_name] = att_val

    #-- CryoSat-2 External Corrections Group
    for key in Geometry_keys:
        val = CS_l1b_mds['Geometry'][key]
        #-- Defining the HDF5 dataset variables
        h5['Geometry'][key] = fileID.create_dataset('Geometry/{0}'.format(key),
            (n_records,), data=val, dtype=val.dtype, compression='gzip')
        #-- attach dimensions
        h5['Geometry'][key].dims[0].label='CS_L1b_MDS_REC_SIZE'
        #-- add HDF5 variable attributes
        for att_name,att_val in CS_l1b_attrib['Geometry'][key].items():
            h5['Geometry'][key].attrs[att_name] = att_val

    #-- CryoSat-2 Average Waveforms Group
    for key in Wfm_1Hz_keys:
        val = CS_l1b_mds['Waveform_1Hz'][key]
        if key in ('Waveform'):
            #-- Defining the HDF5 dataset variables
            h5['Waveform_1Hz'][key] = fileID.create_dataset('Waveform_1Hz/{0}'.format(key),
                (n_records,n_1Hz_wfm), data=val, dtype=val.dtype,
                compression='gzip')
            #-- attach dimensions
            h5['Waveform_1Hz'][key].dims[0].label='CS_L1b_MDS_REC_SIZE'
            h5['Waveform_1Hz'][key].dims[1].label='CS_L1b_MDS_1Hz_WAVEFORM_SIZE'
        else:
            #-- Defining the HDF5 dataset variables
            h5['Waveform_1Hz'][key] = fileID.create_dataset('Waveform_1Hz/{0}'.format(key),
                (n_records,), data=val, dtype=val.dtype, compression='gzip')
            #-- attach dimensions
            h5['Waveform_1Hz'][key].dims[0].label='CS_L1b_MDS_REC_SIZE'
        #-- add HDF5 variable attributes
        for att_name,att_val in CS_l1b_attrib['Waveform_1Hz'][key].items():
            h5['Waveform_1Hz'][key].attrs[att_name] = att_val

    #-- CryoSat-2 Waveforms Group with Beam Behavior Parameters
    for key in Wfm_20Hz_keys:
        val = CS_l1b_mds['Waveform_20Hz'][key]
        if key in ('Waveform','Coherence','Phase_diff'):
            #-- Defining the HDF5 dataset variables
            h5['Waveform_20Hz'][key] = fileID.create_dataset('Waveform_20Hz/{0}'.format(key),
                (n_records,n_blocks,n_20Hz_wfm,), data=val, dtype=val.dtype,
                compression='gzip')
            #-- attach dimensions
            h5['Waveform_20Hz'][key].dims[0].label='CS_L1b_MDS_REC_SIZE'
            h5['Waveform_20Hz'][key].dims[1].label='CS_L1b_MDS_BLOCK_SIZE'
            h5['Waveform_20Hz'][key].dims[2].label='CS_L1b_MDS_20Hz_WAVEFORM_SIZE'
            #-- add HDF5 variable attributes
            for att_name,att_val in CS_l1b_attrib['Waveform_20Hz'][key].items():
                h5['Waveform_20Hz'][key].attrs[att_name] = att_val
        elif key in ('Beam'):
            h5['Waveform_20Hz'][key] = {}
            for ds_name in beam_keys:
                ds_val = val[ds_name]
                #-- Defining the HDF5 dataset variables
                h5['Waveform_20Hz'][key][ds_name] = fileID.create_dataset(
                    'Waveform_20Hz/{0}/{1}'.format(key,ds_name), (n_records,n_blocks,),
                    data=ds_val, dtype=ds_val.dtype, compression='gzip')
                #-- attach dimensions
                h5['Waveform_20Hz'][key][ds_name].dims[0].label='CS_L1b_MDS_REC_SIZE'
                h5['Waveform_20Hz'][key][ds_name].dims[1].label='CS_L1b_MDS_BLOCK_SIZE'
                #-- add HDF5 variable attributes
                for att_name,att_val in CS_l1b_attrib['Waveform_20Hz'][key][ds_name].items():
                    h5['Waveform_20Hz'][key][ds_name].attrs[att_name] = att_val
        else:
            #-- Defining the HDF5 dataset variables
            h5['Waveform_20Hz'][key] = fileID.create_dataset('Waveform_20Hz/{0}'.format(key),
                (n_records,n_blocks,), data=val, dtype=val.dtype, compression='gzip')
            #-- attach dimensions
            h5['Waveform_20Hz'][key].dims[0].label='CS_L1b_MDS_REC_SIZE'
            h5['Waveform_20Hz'][key].dims[1].label='CS_L1b_MDS_BLOCK_SIZE'
            #-- add HDF5 variable attributes
            for att_name,att_val in CS_l1b_attrib['Waveform_20Hz'][key].items():
                h5['Waveform_20Hz'][key].attrs[att_name] = att_val

    #-- output MPH/SPH/DSD headers as group attributes
    if (HEADER == 1):
        #-- HEADER 1 is for single CryoSat-2 files
        fileID.create_group('METADATA')
        fileID['METADATA'].create_group('MPH')
        fileID['METADATA'].create_group('SPH')
        fileID['METADATA'].create_group('DSD')
        #-- Main Product Header (MPH) are all strings
        for att_name,att_val in CS_l1b_mds['METADATA']['MPH'].items():
            fileID['METADATA']['MPH'].attrs[att_name] = att_val
        #-- Specific Product Header (SPH) are both strings and dictionaries
        for att_name,att_val in CS_l1b_mds['METADATA']['SPH'].items():
            if isinstance(att_val,dict):
                #-- if att_val is dictionary
                fileID['METADATA']['SPH'].create_group(att_name)
                for ds_name,ds_val in att_val.items():
                    fileID['METADATA']['SPH'][att_name].attrs[ds_name] = ds_val
            elif isinstance(att_val,str) and att_name:
                #-- if att_val is string
                fileID['METADATA']['SPH'].attrs[att_name] = att_val
        #-- Data Set Descriptors (DSD) are all strings
        for att_name,att_val in CS_l1b_mds['METADATA']['DSD'].items():
            fileID['METADATA']['DSD'].attrs[att_name] = att_val
    elif (HEADER == 2):
        #-- HEADER 2 is for merged CryoSat-2 files from convert_cryosat_L1b.py
        fileID.create_group('METADATA')
        fileID['METADATA'].create_group('MPH')
        fileID['METADATA'].create_group('SPH')
        fileID['METADATA'].create_group('DSD')
        #-- Main Product Header (MPH) are all strings
        for fi in CS_l1b_mds['METADATA']['MPH'].keys():
            fileID['METADATA']['MPH'].create_group(fi)
            for att_name,att_val in CS_l1b_mds['METADATA']['MPH'][fi].items():
                fileID['METADATA']['MPH'][fi].attrs[att_name] = att_val
        #-- Specific Product Header (SPH) are both strings and dictionaries
        for fi in CS_l1b_mds['METADATA']['SPH'].keys():
            fileID['METADATA']['SPH'].create_group(fi)
            for att_name,att_val in CS_l1b_mds['METADATA']['SPH'][fi].items():
                if isinstance(att_val,dict):
                    #-- if att_val is dictionary
                    fileID['METADATA']['SPH'][fi].create_group(att_name)
                    for dsn,dsv in att_val.items():
                        fileID['METADATA']['SPH'][fi][att_name].attrs[dsn] = dsv
                elif isinstance(att_val,str) and att_name:
                    #-- if att_val is string
                    fileID['METADATA']['SPH'][fi].attrs[att_name] = att_val
        #-- Data Set Descriptors (DSD) are all strings
        for fi in CS_l1b_mds['METADATA']['DSD'].keys():
            fileID['METADATA']['DSD'].create_group(fi)
            for att_name,att_val in CS_l1b_mds['METADATA']['DSD'][fi].items():
                fileID['METADATA']['DSD'][fi].attrs[att_name] = att_val

    #-- output file title
    fileID.attrs['description'] = TITLE

    #-- Output HDF5 structure information
    if VERBOSE:
        print(FILENAME)
        print(list(fileID.keys()))

    #-- Closing the HDF5 file
    fileID.close()

#-- PURPOSE: read CryoSat-2 HDF5 files
def read_HDF5_cryosat_L1b(FILENAME, ATTRIBUTES=True, VERBOSE=False):
    #-- Open the HDF5 file for reading
    fileID = h5py.File(os.path.expanduser(FILENAME), 'r')

    #-- Output HDF5 file information
    if VERBOSE:
        print(fileID.filename)
        print(list(fileID.keys()))

    #-- allocate python dictionaries for output CS_l1b_mds variables
    CS_l1b_mds = {}
    CS_l1b_mds['Location'] = {}
    CS_l1b_mds['Data'] = {}
    CS_l1b_mds['Geometry'] = {}
    CS_l1b_mds['Waveform_1Hz'] = {}
    CS_l1b_mds['Waveform_20Hz'] = {}
    #-- get each HDF5 variable
    #-- CryoSat-2 Location Group
    for key in fileID['Location'].keys():
        if key in ('Sat_velocity','Real_beam','Baseline'):
            CS_l1b_mds['Location'][key] = fileID['Location'][key][:,:,:]
        else:
            CS_l1b_mds['Location'][key] = fileID['Location'][key][:,:]
    #-- CryoSat-2 Measurement Group
    for key in fileID['Data'].keys():
        CS_l1b_mds['Data'][key] = fileID['Data'][key][:,:]
    #-- CryoSat-2 External Corrections Group
    for key in fileID['Geometry'].keys():
        CS_l1b_mds['Geometry'][key] = fileID['Geometry'][key][:]
    #-- CryoSat-2 Average Waveform Group
    for key in fileID['Waveform_1Hz'].keys():
        if key in ('Waveform'):
            CS_l1b_mds['Waveform_1Hz'][key] = fileID['Waveform_1Hz'][key][:,:]
        else:
            CS_l1b_mds['Waveform_1Hz'][key] = fileID['Waveform_1Hz'][key][:]
    #-- CryoSat-2 Waveform Group
    for key in fileID['Waveform_20Hz'].keys():
        if key in ('Waveform','Coherence','Phase_diff'):
            CS_l1b_mds['Waveform_20Hz'][key] = fileID['Waveform_20Hz'][key][:,:,:]
        elif key in ('Beam'):
            CS_l1b_mds['Waveform_20Hz'][key] = {}
            for ds_name,ds_val in fileID['Waveform_20Hz'][key].items():
                CS_l1b_mds['Waveform_20Hz'][key][ds_name] = ds_val[:,:]
        else:
            CS_l1b_mds['Waveform_20Hz'][key] = fileID['Waveform_20Hz'][key][:,:]

    #-- Getting attributes of included variables
    if ATTRIBUTES:
        #-- allocate python dictionaries for output CS_l1b_mds attributes
        CS_l1b_mds['Attributes'] = {}
        CS_l1b_mds['Attributes']['Location'] = {}
        CS_l1b_mds['Attributes']['Data'] = {}
        CS_l1b_mds['Attributes']['Geometry'] = {}
        CS_l1b_mds['Attributes']['Waveform_1Hz'] = {}
        CS_l1b_mds['Attributes']['Waveform_20Hz'] = {}
        #-- CryoSat-2 Location Group
        for key in fileID['Location'].keys():
            CS_l1b_mds['Attributes']['Location'][key] = {}
            for att_name,att_val in fileID['Location'][key].attrs.items():
                CS_l1b_mds['Attributes']['Location'][key][att_name] = att_val
        #-- CryoSat-2 Measurement Group
        for key in fileID['Data'].keys():
            CS_l1b_mds['Attributes']['Data'][key] = {}
            for att_name,att_val in fileID['Data'][key].attrs.items():
                CS_l1b_mds['Attributes']['Data'][key][att_name] = att_val
        #-- CryoSat-2 External Corrections Group
        for key in fileID['Geometry'].keys():
            CS_l1b_mds['Attributes']['Geometry'][key] = {}
            for att_name,att_val in fileID['Geometry'][key].attrs.items():
                CS_l1b_mds['Attributes']['Geometry'][key][att_name] = att_val
        #-- CryoSat-2 Average Waveform Group
        for key in fileID['Waveform_1Hz'].keys():
            CS_l1b_mds['Attributes']['Waveform_1Hz'][key] = {}
            for att_name,att_val in fileID['Waveform_1Hz'][key].attrs.items():
                CS_l1b_mds['Attributes']['Waveform_1Hz'][key][att_name] = att_val
        #-- CryoSat-2 Waveform Group
        for key in fileID['Waveform_20Hz'].keys():
            if key in ('Beam'):
                CS_l1b_mds['Attributes']['Waveform_20Hz'][key] = {}
                for dsn in fileID['Waveform_20Hz'][key].keys():
                    CS_l1b_mds['Attributes']['Waveform_20Hz'][key][dsn] = {}
                    for atn,atv in fileID['Waveform_20Hz'][key][dsn].attrs.items():
                        CS_l1b_mds['Attributes']['Waveform_20Hz'][key][dsn][atn] = atv
            else:
                CS_l1b_mds['Attributes']['Waveform_20Hz'][key] = {}
                for atn,atv in fileID['Waveform_20Hz'][key].attrs.items():
                    CS_l1b_mds['Attributes']['Waveform_20Hz'][key][atn] = atv
        #-- Global attribute description
        CS_l1b_mds['Attributes']['title'] = fileID.attrs['description']

    #-- Closing the HDF5 file
    fileID.close()
    return CS_l1b_mds

#-- PURPOSE: get the number of records and number of blocks in an HDF5 file
def HDF5_cryosat_L1b_shape(FILENAME):
    #-- Open the HDF5 file for reading
    with h5py.File(os.path.expanduser(FILENAME), 'r') as fid:
        n_records,n_blocks = fid['Location']['Day'].shape
        n_1Hz_wfm = fid['Waveform_1Hz']['Waveform'].shape[1]
        n_20Hz_wfm = fid['Waveform_20Hz']['Waveform'].shape[2]
    return (n_records,n_blocks,n_1Hz_wfm,n_20Hz_wfm)

#-- PURPOSE: get attribute names for baseline
def cryosat_L1b_attributes(MODE, BASELINE):
    #-- CryoSat-2 Time and Orbit Group
    L1b_location_attributes = {}
    #-- Time: day part
    L1b_location_attributes['Day'] = {}
    L1b_location_attributes['Day']['long_name'] = 'MDSR time stamp days'
    L1b_location_attributes['Day']['units'] = 'days since 2000-01-01 00:00:00 TAI'
    L1b_location_attributes['Day']['hertz'] = 20
    #-- Time: second part
    L1b_location_attributes['Second'] = {}
    L1b_location_attributes['Second']['long_name'] = 'MDSR time stamp seconds'
    L1b_location_attributes['Second']['units'] = 'seconds'
    L1b_location_attributes['Second']['hertz'] = 20
    #-- Time: microsecond part
    L1b_location_attributes['Micsec'] = {}
    L1b_location_attributes['Micsec']['long_name'] = 'MDSR time stamp microseconds'
    L1b_location_attributes['Micsec']['units'] = 'microseconds'
    L1b_location_attributes['Micsec']['hertz'] = 20
    #-- USO correction factor
    L1b_location_attributes['USO_Corr'] = {}
    L1b_location_attributes['USO_Corr']['long_name'] = ('DORIS Ultra Stable '
        'Oscillator drift correction factor')
    L1b_location_attributes['USO_Corr']['description'] = ('USO_Corr_Factor = '
        'USO_freq_nominal / (USO_freq_nominal + model_freq_deviation). '
        'USO_freq_nominal is the nominal frequency provided in the IPF database.'
        'model_freq_deviation is the modelled frequency deviation provided by '
        'the DORIS USO drift file')
    L1b_location_attributes['USO_Corr']['units'] = '1e-15'
    L1b_location_attributes['USO_Corr']['hertz'] = 20
    #-- Mode ID
    L1b_location_attributes['Mode_ID'] = {}
    L1b_location_attributes['Mode_ID']['long_name'] = 'Mode ID'
    L1b_location_attributes['Mode_ID']['description'] = ('Identifies the SIRAL '
        'instrument measurement mode. See table 2.3.3-2 of the "L1b Products '
        'Format Specification" document')
    L1b_location_attributes['Mode_ID']['units'] = '1e-15'
    L1b_location_attributes['Mode_ID']['flag_meanings'] = 'lrm sar sarin'
    L1b_location_attributes['Mode_ID']['hertz'] = 20
    #-- Mode Flags
    L1b_location_attributes['Mode_flags'] = {}
    L1b_location_attributes['Mode_flags']['long_name'] = 'Mode flags'
    L1b_location_attributes['Mode_flags']['description'] = ('Flags related to '
        'sub-modes of SARIn mode from instrument configuration bits in L0. '
        'Identifies the sarin degraded case and the CAL4 flag')
    L1b_location_attributes['Mode_flags']['flag_meanings'] = ('sarin_degraded_case '
        'cal4_packet_detection')
    L1b_location_attributes['Mode_flags']['hertz'] = 20
    #-- Platform attitude control mode
    L1b_location_attributes['Att_control'] = {}
    L1b_location_attributes['Att_control']['long_name'] = 'Platform Attitude Control'
    L1b_location_attributes['Att_control']['description'] = ('Platform attitude '
        'control mode from instrument configuration bits in L0.')
    L1b_location_attributes['Att_control']['flag_meanings'] = ('unknown '
        'local_normal_pointing yaw_steering')
    L1b_location_attributes['Att_control']['hertz'] = 20
    #-- Source sequence counter
    L1b_location_attributes['SSC'] = {}
    L1b_location_attributes['SSC']['long_name'] = 'Source sequence counter'
    L1b_location_attributes['SSC']['description'] = ('Read from the L0 echo '
        'telemetry packet (of the master channel in the case of SARin).'
        'This is a 16384 cyclic modulo counter, starting from 0, incrementing '
        'by 1. A separate counter is maintained for each instrument mode')
    L1b_location_attributes['SSC']['hertz'] = 20
    #-- Instrument configuration
    L1b_location_attributes['Inst_config'] = {}
    L1b_location_attributes['Inst_config']['long_name'] = ('Instrument '
        'Configuration flag')
    L1b_location_attributes['Inst_config']['description'] = ('This is derived '
        'from flags in the L0 packets for tracking and the echo. See table '
        '2.3.3-3 of the "L1b Products Format Specification" document.')
    L1b_location_attributes['Inst_config']['flag_meanings'] = ('siral_redundant '
        'external_cal open_loop loss_of_echo real_time_error echo_saturation '
        'rx_band_attenuated cycle_report_error')
    L1b_location_attributes['Inst_config']['hertz'] = 20
    #-- acquisition band
    L1b_location_attributes['Inst_band'] = {}
    L1b_location_attributes['Inst_band']['long_name'] = 'Acquisition Band'
    L1b_location_attributes['Inst_band']['description'] = ('This flag '
        'contains the acquisition band of the SIRAL instrument.')
    L1b_location_attributes['Inst_band']['flag_meanings'] = ('unknown '
        '320_mhz 40_mhz')
    L1b_location_attributes['Inst_band']['hertz'] = 20
    #-- instrument channel
    L1b_location_attributes['Inst_channel'] = {}
    L1b_location_attributes['Inst_channel']['long_name'] = 'Rx Channel'
    L1b_location_attributes['Inst_channel']['description'] = ('This flag '
        'contains the SIRAL instrument channel in use.')
    L1b_location_attributes['Inst_channel']['flag_meanings'] = ('unknown '
        'rx1 rx2 both')
    L1b_location_attributes['Inst_channel']['hertz'] = 20
    #-- tracking mode
    L1b_location_attributes['Tracking_mode'] = {}
    L1b_location_attributes['Tracking_mode']['long_name'] = 'Tracking Mode'
    L1b_location_attributes['Tracking_mode']['description'] = ('This flag '
        'indicates the tracking mode of the SIRAL instrument.')
    L1b_location_attributes['Tracking_mode']['flag_meanings'] = ('unknown '
        'lrm sar sarin')
    L1b_location_attributes['Tracking_mode']['hertz'] = 20
    #-- Record Counter
    L1b_location_attributes['Rec_Count'] = {}
    L1b_location_attributes['Rec_Count']['long_name'] = ('Instrument '
        'Configuration flag')
    L1b_location_attributes['Rec_Count']['description'] = ('Progressive counter '
        'incremented by 1 for each data block. Hence the first full MDS record '
        'contains the numbers 1-20, the second 21-40, etc.')
    L1b_location_attributes['Rec_Count']['hertz'] = 20
    #-- Lat: packed units (0.1 micro-degree, 1e-7 degrees)
    L1b_location_attributes['Lat'] = {}
    L1b_location_attributes['Lat']['long_name'] = 'Latitude of measurement'
    L1b_location_attributes['Lat']['description'] = ('Corresponding to the '
        'position at the MDSR Time Stamp')
    L1b_location_attributes['Lat']['units'] = '0.1 micro-degree'
    L1b_location_attributes['Lat']['valid_min'] = -9e8
    L1b_location_attributes['Lat']['valid_max'] = 9e8
    L1b_location_attributes['Lat']['hertz'] = 20
    #-- Lon: packed units (0.1 micro-degree, 1e-7 degrees)
    L1b_location_attributes['Lon'] = {}
    L1b_location_attributes['Lon']['long_name'] = 'Longitude of measurement'
    L1b_location_attributes['Lon']['description'] = ('Corresponding to the '
        'position at the MDSR Time Stamp')
    L1b_location_attributes['Lon']['units'] = '0.1 micro-degree'
    L1b_location_attributes['Lon']['valid_min'] = -18e8
    L1b_location_attributes['Lon']['valid_max'] = 18e8
    L1b_location_attributes['Lon']['hertz'] = 20
    #-- Alt: packed units (mm, 1e-3 m)
    #-- Altitude of COG above reference ellipsoid
    L1b_location_attributes['Alt'] = {}
    L1b_location_attributes['Alt']['long_name'] = 'Altitude'
    L1b_location_attributes['Alt']['description'] = ('Altitude of Satellite '
        'COG above reference ellipsoid corresponding to the MDSR Time Stamp')
    L1b_location_attributes['Alt']['units'] = 'millimeters'
    L1b_location_attributes['Alt']['hertz'] = 20
    #-- Instantaneous altitude rate derived from orbit: packed units (mm/s, 1e-3 m/s)
    L1b_location_attributes['Alt_rate'] = {}
    L1b_location_attributes['Alt_rate']['long_name'] = 'Altitude Rate'
    L1b_location_attributes['Alt_rate']['description'] = ('Instantaneous '
        'altitude rate derived from orbit corresponding to the MDSR Time Stamp')
    L1b_location_attributes['Alt_rate']['units'] = 'millimeters/second'
    L1b_location_attributes['Alt_rate']['hertz'] = 20
    #-- Satellite velocity vector. In ITRF: packed units (mm/s, 1e-3 m/s)
    #-- ITRF= International Terrestrial Reference Frame
    L1b_location_attributes['Sat_velocity'] = {}
    L1b_location_attributes['Sat_velocity']['long_name'] = ('Satellite velocity '
        'vector')
    L1b_location_attributes['Sat_velocity']['description'] = ('In the '
        'International Terrestrial Reference Frame (ITRF) in the International '
        'Earth Fixed System. From Orbit CFI call. This is not a unit vector as '
        'the velocity magnitude is also required.')
    L1b_location_attributes['Sat_velocity']['units'] = 'millimeters/second'
    L1b_location_attributes['Sat_velocity']['hertz'] = 20
    #-- Real beam direction vector. In CRF: packed units (micro-m, 1e-6 m)
    #-- CRF= CryoSat Reference Frame.
    L1b_location_attributes['Real_beam'] = {}
    L1b_location_attributes['Real_beam']['long_name'] = ('Real beam direction '
        'vector')
    L1b_location_attributes['Real_beam']['description'] = ('In the '
        'CryoSat Reference Frame (CRF). This is a unit vector.')
    L1b_location_attributes['Real_beam']['units'] = 'micrometers'
    L1b_location_attributes['Real_beam']['hertz'] = 20
    #-- Interferometric baseline vector. In CRF: packed units (micro-m, 1e-6 m)
    L1b_location_attributes['Baseline'] = {}
    L1b_location_attributes['Baseline']['long_name'] = ('Interferometric '
        'baseline vector')
    L1b_location_attributes['Baseline']['description'] = ('In the '
        'CryoSat Reference Frame (CRF). This is a unit vector.')
    L1b_location_attributes['Baseline']['units'] = 'micrometers'
    L1b_location_attributes['Baseline']['hertz'] = 20
    #-- Star Tracker ID and Spacecraft mispointing for Baseline-C and D
    if BASELINE in ('C','D'):
        #-- Star Tracker ID
        L1b_location_attributes['ST_ID'] = {}
        L1b_location_attributes['ST_ID']['long_name'] = 'Star Tracker ID'
        L1b_location_attributes['ST_ID']['hertz'] = 20
        #-- Antenna Bench Roll Angle (Derived from star trackers)
        #-- packed units (0.1 micro-degree, 1e-7 degrees)
        L1b_location_attributes['Roll'] = {}
        L1b_location_attributes['Roll']['long_name'] = ('Antenna Bench Roll '
            'Angle derived from star trackers')
        L1b_location_attributes['Roll']['units'] = '0.1 micro-degree'
        L1b_location_attributes['Roll']['hertz'] = 20
        #-- Antenna Bench Pitch Angle (Derived from star trackers)
        #-- packed units (0.1 micro-degree, 1e-7 degrees)
        L1b_location_attributes['Pitch'] = {}
        L1b_location_attributes['Pitch']['long_name'] = ('Antenna Bench Pitch '
            'Angle derived from star trackers')
        L1b_location_attributes['Pitch']['units'] = '0.1 micro-degree'
        L1b_location_attributes['Pitch']['hertz'] = 20
        #-- Antenna Bench Yaw Angle (Derived from star trackers)
        #-- packed units (0.1 micro-degree, 1e-7 degrees)
        L1b_location_attributes['Yaw'] = {}
        L1b_location_attributes['Yaw']['long_name'] = ('Antenna Bench Yaw '
            'Angle derived from star trackers')
        L1b_location_attributes['Yaw']['units'] = '0.1 micro-degree'
        L1b_location_attributes['Yaw']['hertz'] = 20
    #-- Measurement Confidence Data Flags
    L1b_location_attributes['MCD'] = {}
    L1b_location_attributes['MCD']['long_name'] = ('Measurement Confidence '
        'Data Flags')
    L1b_location_attributes['MCD']['description'] = ('Generally the MCD flags '
        'indicate problems when set. If MCD is 0 then no problems or non-nominal '
        'conditions were detected. Serious errors are indicated by setting bit 31')
    L1b_location_attributes['MCD']['flag_meanings'] = ('block_degraded '
        'blank_block datation_degraded orbit_prop_error orbit_file_change '
        'orbit_gap echo_saturated other_echo_error sarin_rx1_error '
        'sarin_rx2_error window_delay_error agc_error cal1_missing '
        'cal1_default doris_uso_missing ccal1_default trk_echo_error '
        'echo_rx1_error echo_rx2_error npm_error cal1_pwr_corr_type '
        'phase_pert_cor_missing cal2_missing cal2_default power_scale_error '
        'attitude_cor_missing phase_pert_cor_default')
    L1b_location_attributes['MCD']['hertz'] = 20

    #-- CryoSat-2 Measurement Group
    L1b_measurement_attributes = {}
    #-- Window Delay reference (two-way) corrected for instrument delays
    L1b_measurement_attributes['TD'] = {}
    L1b_measurement_attributes['TD']['long_name'] = 'Window Delay'
    L1b_measurement_attributes['TD']['description'] = ('Window delay from the '
        'telemetry converted to physical units. This is a 2-way measurement: '
        'the time taken for the radar pulse to travel to the surface and back.'
        'COM offset is applied and Calibration correction from CAL1 is applied')
    L1b_measurement_attributes['TD']['units'] = 'picoseconds'
    L1b_measurement_attributes['TD']['hertz'] = 20
    #-- H0 Initial Height Word from telemetry
    L1b_measurement_attributes['H_0'] = {}
    L1b_measurement_attributes['H_0']['long_name'] = 'H0 Initial Height Word'
    L1b_measurement_attributes['H_0']['units'] = '48.8 ps'
    L1b_measurement_attributes['H_0']['hertz'] = 20
    #-- COR2 Height Rate: on-board tracker height rate over the radar cycle
    L1b_measurement_attributes['COR2'] = {}
    L1b_measurement_attributes['COR2']['long_name'] = 'COR2 Height Rate'
    L1b_measurement_attributes['COR2']['description'] = ('On-board tracker '
        'height rate over the radar cycle')
    L1b_measurement_attributes['COR2']['units'] = '3.05 ps/rc'
    L1b_measurement_attributes['COR2']['hertz'] = 20
    #-- Coarse Range Word (LAI) derived from telemetry
    L1b_measurement_attributes['LAI'] = {}
    L1b_measurement_attributes['LAI']['long_name'] = 'Coarse Range Word'
    L1b_measurement_attributes['LAI']['units'] = '12.5 ns'
    L1b_measurement_attributes['LAI']['hertz'] = 20
    #-- Fine Range Word (FAI) derived from telemetry
    L1b_measurement_attributes['FAI'] = {}
    L1b_measurement_attributes['FAI']['long_name'] = 'Fine Range Word'
    L1b_measurement_attributes['FAI']['units'] = '12.5/256 ns'
    L1b_measurement_attributes['FAI']['hertz'] = 20
    #-- Automatic Gain Control Channel 1: AGC gain applied on Rx channel 1.
    #-- Gain calibration corrections are applied (Sum of AGC stages 1 and 2
    #-- plus the corresponding corrections) (dB/100)
    L1b_measurement_attributes['AGC_CH1'] = {}
    L1b_measurement_attributes['AGC_CH1']['long_name'] = ('Automatic Gain '
        'Control Channel 1')
    L1b_measurement_attributes['AGC_CH1']['description'] = ('AGC gain applied '
        'on Rx channel 1. Gain calibration corrections are applied (Sum of AGC '
        'stages 1 and 2 plus the corresponding corrections)')
    L1b_measurement_attributes['AGC_CH1']['units'] = 'dB/100'
    L1b_measurement_attributes['AGC_CH1']['hertz'] = 20
    #-- Automatic Gain Control Channel 2: AGC gain applied on Rx channel 2.
    #-- Gain calibration corrections are applied (dB/100)
    L1b_measurement_attributes['AGC_CH2'] = {}
    L1b_measurement_attributes['AGC_CH2']['long_name'] = ('Automatic Gain '
        'Control Channel 2')
    L1b_measurement_attributes['AGC_CH2']['description'] = ('AGC gain applied '
        'on Rx channel 2.')
    L1b_measurement_attributes['AGC_CH2']['units'] = 'dB/100'
    L1b_measurement_attributes['AGC_CH2']['hertz'] = 20
    #-- Total Fixed Gain On Channel 1: gain applied by the RF unit. (dB/100)
    L1b_measurement_attributes['TR_gain_CH1'] = {}
    L1b_measurement_attributes['TR_gain_CH1']['long_name'] = ('Total Fixed Gain '
        'On Channel 1')
    L1b_measurement_attributes['TR_gain_CH1']['description'] = ('Gain applied '
        'by the RF unit.')
    L1b_measurement_attributes['TR_gain_CH1']['units'] = 'dB/100'
    L1b_measurement_attributes['TR_gain_CH1']['hertz'] = 20
    #-- Total Fixed Gain On Channel 2: gain applied by the RF unit. (dB/100)
    L1b_measurement_attributes['TR_gain_CH2'] = {}
    L1b_measurement_attributes['TR_gain_CH2']['long_name'] = ('Total Fixed Gain '
        'On Channel 2')
    L1b_measurement_attributes['TR_gain_CH2']['description'] = ('Gain applied '
        'by the RF unit.')
    L1b_measurement_attributes['TR_gain_CH2']['units'] = 'dB/100'
    L1b_measurement_attributes['TR_gain_CH2']['hertz'] = 20
    #-- Transmit Power in microWatts
    L1b_measurement_attributes['TX_Power'] = {}
    L1b_measurement_attributes['TX_Power']['long_name'] = 'Transmit Power'
    L1b_measurement_attributes['TX_Power']['units'] = 'microWatts'
    L1b_measurement_attributes['TX_Power']['hertz'] = 20
    #-- Doppler range correction: Radial component (mm)
    #-- computed for the component of satellite velocity in the nadir direction
    L1b_measurement_attributes['Doppler_range'] = {}
    L1b_measurement_attributes['Doppler_range']['long_name'] = ('Doppler range '
        'correction')
    L1b_measurement_attributes['Doppler_range']['description'] = ('Radial '
        'component computed for the component of satellite velocity in the '
        'nadir direction.')
    L1b_measurement_attributes['Doppler_range']['units'] = 'mm'
    L1b_measurement_attributes['Doppler_range']['hertz'] = 20
    #-- Value of Doppler Angle for the first single look echo (1e-7 radians)
    L1b_measurement_attributes['Doppler_angle_start'] = {}
    L1b_measurement_attributes['Doppler_angle_start']['long_name'] = ('Doppler '
        'angle start')
    L1b_measurement_attributes['Doppler_angle_start']['description'] = ('Value '
        'of Doppler Angle for the first single look echo in the stack. It is '
        'the angle between: (a) direction perpendicular to the velocity '
        'vector, (b) direction from satellite to surface location. The Doppler '
        'angle depends on velocity vector and on geometry.')
    L1b_measurement_attributes['Doppler_angle_start']['units'] = '1e-7 radians'
    L1b_measurement_attributes['Doppler_angle_start']['hertz'] = 20
    #-- Value of Doppler Angle for the last single look echo (1e-7 radians)
    L1b_measurement_attributes['Doppler_angle_stop'] = {}
    L1b_measurement_attributes['Doppler_angle_stop']['long_name'] = ('Doppler '
        'angle stop')
    L1b_measurement_attributes['Doppler_angle_stop']['description'] = ('Value '
        'of Doppler Angle for the last single look echo in the stack. It is '
        'the angle between: (a) direction perpendicular to the velocity '
        'vector, (b) direction from satellite to surface location. The Doppler '
        'angle depends on velocity vector and on geometry.')
    L1b_measurement_attributes['Doppler_angle_stop']['units'] = '1e-7 radians'
    L1b_measurement_attributes['Doppler_angle_stop']['hertz'] = 20
    #-- Instrument Range Correction: transmit-receive antenna (mm)
    #-- Calibration correction to range on channel 1 computed from CAL1.
    L1b_measurement_attributes['TR_inst_range'] = {}
    L1b_measurement_attributes['TR_inst_range']['long_name'] = ('Instrument '
        'Range Correction: transmit-receive antenna')
    L1b_measurement_attributes['TR_inst_range']['description'] = ('Calibration '
        'correction to range on channel 1 computed from CAL1.')
    L1b_measurement_attributes['TR_inst_range']['units'] = 'mm'
    L1b_measurement_attributes['TR_inst_range']['hertz'] = 20
    #-- Instrument Range Correction: receive-only antenna (mm)
    #-- Calibration correction to range on channel 2 computed from CAL1.
    L1b_measurement_attributes['R_inst_range'] = {}
    L1b_measurement_attributes['R_inst_range']['long_name'] = ('Instrument '
        'Range Correction: receive-only antenna')
    L1b_measurement_attributes['R_inst_range']['description'] = ('Calibration '
        'correction to range on channel 2 computed from CAL1.')
    L1b_measurement_attributes['R_inst_range']['units'] = 'mm'
    L1b_measurement_attributes['R_inst_range']['hertz'] = 20
    #-- Instrument Gain Correction: transmit-receive antenna (dB/100)
    #-- Calibration correction to gain on channel 1 computed from CAL1
    L1b_measurement_attributes['TR_inst_gain'] = {}
    L1b_measurement_attributes['TR_inst_gain']['long_name'] = ('Instrument '
        'Gain Correction: transmit-receive antenna')
    L1b_measurement_attributes['TR_inst_gain']['description'] = ('Calibration '
        'correction to gain on channel 1 computed from CAL1.')
    L1b_measurement_attributes['TR_inst_gain']['units'] = 'dB/100'
    L1b_measurement_attributes['TR_inst_gain']['hertz'] = 20
    #-- Instrument Gain Correction: receive-only (dB/100)
    #-- Calibration correction to gain on channel 2 computed from CAL1
    L1b_measurement_attributes['R_inst_gain'] = {}
    L1b_measurement_attributes['R_inst_gain']['long_name'] = ('Instrument '
        'Gain Correction: receive-only antenna')
    L1b_measurement_attributes['R_inst_gain']['description'] = ('Calibration '
        'correction to gain on channel 2 computed from CAL1.')
    L1b_measurement_attributes['R_inst_gain']['units'] = 'dB/100'
    L1b_measurement_attributes['R_inst_gain']['hertz'] = 20
    #-- Internal Phase Correction (microradians)
    L1b_measurement_attributes['Internal_phase'] = {}
    L1b_measurement_attributes['Internal_phase']['long_name'] = ('Internal '
        'Phase Correction')
    L1b_measurement_attributes['Internal_phase']['description'] = ('Set to '
        'zero due to no availability of correction until specialized FBR-L1B '
        'processing.')
    L1b_measurement_attributes['Internal_phase']['units'] = 'microradians'
    L1b_measurement_attributes['Internal_phase']['hertz'] = 20
    #-- External Phase Correction (microradians)
    L1b_measurement_attributes['External_phase'] = {}
    L1b_measurement_attributes['External_phase']['long_name'] = ('External '
        'Phase Correction')
    L1b_measurement_attributes['External_phase']['description'] = ('Taken from '
        'the IPFDB file (SARIN only) to be added to the internal phase '
        'correction term. The external phase correction is the temperature-'
        'averaged component of external inter-channel phase difference derived '
        'from phase difference sensitive antenna subsystem, waveguides and '
        'instrument waveguide switches. The external phase correction does not '
        'contain internal instrument temperature dependent effects of '
        'calibration coupler and duplexer which are dealt with by the CAL-4 '
        'signal. These CAL-4 data are processed to compute the internal phase '
        'correction parameter.')
    L1b_measurement_attributes['External_phase']['units'] = 'microradians'
    L1b_measurement_attributes['External_phase']['hertz'] = 20
    #-- Noise Power measurement (dB/100): converted from telemetry units to be
    #-- the noise floor of FBR measurement echoes.
    #-- Set to -9999.99 when the telemetry contains zero.
    L1b_measurement_attributes['Noise_power'] = {}
    L1b_measurement_attributes['Noise_power']['long_name'] = 'Noise power'
    L1b_measurement_attributes['Noise_power']['description'] = ('Noise power '
        'measurement converted from telemetry units to be the noise floor of '
        'FBR measurement echoes. This field is set to the default value equal '
        'to -9999.99 when the telemetry contains zero.')
    L1b_measurement_attributes['Noise_power']['units'] = 'dB/100'
    L1b_measurement_attributes['Noise_power']['_FillValue'] = -9999.99
    L1b_measurement_attributes['Noise_power']['hertz'] = 20
    #-- Phase slope correction (microradians)
    #-- Computed from the CAL-4 packets during the azimuth impulse response
    #-- amplitude (SARIN only). Set from the latest available CAL-4 packet.
    L1b_measurement_attributes['Phase_slope'] = {}
    L1b_measurement_attributes['Phase_slope']['long_name'] = ('Phase Slope '
        'Correction')
    L1b_measurement_attributes['Phase_slope']['description'] = ('Differential '
        'group delay phase difference slope correction (across the whole '
        'bandwidth): fixed and variable group delays introduce a phase '
        'difference slope across the instrument bandwidth. Fixed elements of '
        'differential group delay have been determined during ground testing '
        'and characterisation and cover the elements of antenna, calibration '
        'coupler, Louis waveguide. These fixed elements can be retrieved from '
        'the IPFDB. Variable elements cover differences between the CAL-1 and '
        'CAL-4 paths and can be computed by processing the CAL-1 and CAL-4 data. '
        'SIR_SAR_1B and SIR_LRM_1B products contain this parameter but is set '
        'to zero. Since the correction can only be made at the rate of the '
        'CAL-4 which is 1 Hz and the measurement group high rate blocks are '
        'provided at 20 Hz the product provides the closest in time to FBR '
        'value of slope correction.')
    L1b_measurement_attributes['Phase_slope']['units'] = 'microradians'
    L1b_measurement_attributes['Phase_slope']['hertz'] = 20

    #-- CryoSat-2 External Corrections Group
    L1b_corr_attributes = {}
    #-- Dry Tropospheric Correction packed units (mm, 1e-3 m)
    L1b_corr_attributes['dryTrop'] = {}
    L1b_corr_attributes['dryTrop']['long_name'] = 'Dry Tropospheric Correction'
    L1b_corr_attributes['dryTrop']['description'] = ('Altimeter range correction'
        ' due to the dry-gas component of the Earths atmosphere')
    L1b_corr_attributes['dryTrop']['units'] = 'millimeters'
    L1b_corr_attributes['dryTrop']['hertz'] = 1
    #-- Wet Tropospheric Correction packed units (mm, 1e-3 m)
    L1b_corr_attributes['wetTrop'] = {}
    L1b_corr_attributes['wetTrop']['long_name'] = 'Wet Tropospheric Correction'
    L1b_corr_attributes['wetTrop']['description'] = ('Altimeter range correction'
        ' due to the water component of the Earths atmosphere')
    L1b_corr_attributes['wetTrop']['units'] = 'millimeters'
    L1b_corr_attributes['wetTrop']['hertz'] = 1
    #-- Inverse Barometric Correction packed units (mm, 1e-3 m)
    L1b_corr_attributes['InvBar'] = {}
    L1b_corr_attributes['InvBar']['long_name'] = 'Inverse Barometric Correction'
    L1b_corr_attributes['InvBar']['description'] = ('Altimeter range correction '
        'for the depression of the ocean surface caused by the local barometric '
        'pressure')
    L1b_corr_attributes['InvBar']['units'] = 'millimeters'
    L1b_corr_attributes['InvBar']['hertz'] = 1
    #-- Dynamic Atmosphere Correction packed units (mm, 1e-3 m)
    L1b_corr_attributes['DAC'] = {}
    L1b_corr_attributes['DAC']['long_name'] = 'Dynamic Atmosphere Correction'
    L1b_corr_attributes['DAC']['description'] = ('Altimeter range correction '
        'for both the Inverse Barometric effect and the high-frequency dynamic '
        'component of the wind effect on the ocean. Only one of inverse '
        'barometric correction and DAC have to be used as they are alternatives')
    L1b_corr_attributes['DAC']['units'] = 'millimeters'
    L1b_corr_attributes['DAC']['hertz'] = 1
    #-- GIM Ionospheric Correction packed units (mm, 1e-3 m)
    L1b_corr_attributes['Iono_GIM'] = {}
    L1b_corr_attributes['Iono_GIM']['long_name'] = 'Ionospheric Correction'
    L1b_corr_attributes['Iono_GIM']['description'] = ('Altimeter range correction '
        'for the delay of the radar pulse caused by free electrons in the '
        'ionosphere. Computed a GPS satellite-derived (GIM) map')
    L1b_corr_attributes['Iono_GIM']['units'] = 'millimeters'
    L1b_corr_attributes['Iono_GIM']['hertz'] = 1
    #-- Model Ionospheric Correction packed units (mm, 1e-3 m)
    L1b_corr_attributes['Iono_model'] = {}
    L1b_corr_attributes['Iono_model']['long_name'] = 'Ionospheric Correction'
    L1b_corr_attributes['Iono_model']['description'] = ('Altimeter range '
        'correction for the delay of the radar pulse caused by free electrons '
        'in the ionosphere. Computed from a simple ionospheric model.')
    L1b_corr_attributes['Iono_model']['units'] = 'millimeters'
    L1b_corr_attributes['Iono_model']['hertz'] = 1
    #-- Ocean tide Correction packed units (mm, 1e-3 m)
    L1b_corr_attributes['ocTideElv'] = {}
    L1b_corr_attributes['ocTideElv']['long_name'] = 'Elastic Ocean Tide'
    L1b_corr_attributes['ocTideElv']['description'] = ('Removes the effect of '
        'local tide and adjusts the measurement to the mean sea surface')
    L1b_corr_attributes['ocTideElv']['units'] = 'millimeters'
    L1b_corr_attributes['ocTideElv']['_FillValue'] = 32767
    L1b_corr_attributes['ocTideElv']['hertz'] = 1
    #-- Long period equilibrium ocean tide Correction packed units (mm, 1e-3 m)
    L1b_corr_attributes['lpeTideElv'] = {}
    L1b_corr_attributes['lpeTideElv']['long_name'] = ('Long-Period Equilibrium '
        'Ocean Tide')
    L1b_corr_attributes['lpeTideElv']['description'] = ('Removes the effect of '
        'the oceanic response to the single tidal forcing.')
    L1b_corr_attributes['lpeTideElv']['units'] = 'millimeters'
    L1b_corr_attributes['lpeTideElv']['_FillValue'] = 32767
    L1b_corr_attributes['lpeTideElv']['hertz'] = 1
    #-- Ocean loading tide Correction packed units (mm, 1e-3 m)
    L1b_corr_attributes['olTideElv'] = {}
    L1b_corr_attributes['olTideElv']['long_name'] = 'Ocean Loading Tide'
    L1b_corr_attributes['olTideElv']['description'] = ('Removes the effect of '
        'local tidal distortion of the Earth crust')
    L1b_corr_attributes['olTideElv']['units'] = 'millimeters'
    L1b_corr_attributes['olTideElv']['hertz'] = 1
    #-- Solid Earth tide Correction packed units (mm, 1e-3 m)
    L1b_corr_attributes['seTideElv'] = {}
    L1b_corr_attributes['seTideElv']['long_name'] = 'Solid Earth Tide'
    L1b_corr_attributes['seTideElv']['description'] = ('Removes the effect of '
        'local tidal distortion in the Earth crust')
    L1b_corr_attributes['seTideElv']['units'] = 'millimeters'
    L1b_corr_attributes['seTideElv']['hertz'] = 1
    #-- Geocentric Polar tide Correction packed units (mm, 1e-3 m)
    L1b_corr_attributes['gpTideElv'] = {}
    L1b_corr_attributes['gpTideElv']['long_name'] = 'Geocentric Polar Tide'
    L1b_corr_attributes['gpTideElv']['description'] = ('Removes a long-period '
        'distortion of the Earth crust caused by variations in the centrifugal '
        'force as the Earth rotational axis moves its geographic location')
    L1b_corr_attributes['gpTideElv']['units'] = 'millimeters'
    L1b_corr_attributes['gpTideElv']['hertz'] = 1
    #-- Surface Type: enumerated key to classify surface at nadir
    #-- 0 = Open Ocean
    #-- 1 = Closed Sea
    #-- 2 = Continental Ice
    #-- 3 = Land
    L1b_corr_attributes['Surf_type'] = {}
    L1b_corr_attributes['Surf_type']['long_name'] = 'Surface Type Flag'
    L1b_corr_attributes['Surf_type']['description'] = ('Enumerated key to '
        'classify surface at nadir provided by a model: (0=Open Ocean, '
        '1=Closed Sea, 2=Continental Ice, 3=Land, 4-7=currently unused)')
    L1b_corr_attributes['Surf_type']['hertz'] = 1
    #-- Corrections Status Flag
    L1b_corr_attributes['Corr_status'] = {}
    L1b_corr_attributes['Corr_status']['long_name'] = 'Corrections Status Flag'
    L1b_corr_attributes['Corr_status']['description'] = ('Shows correction '
        'algorithms called in processing. See table 2.3.3-5 of the "L1b '
        'Products Format Specification" document')
    L1b_corr_attributes['Corr_status']['hertz'] = 1
    #-- Correction Error Flag
    L1b_corr_attributes['Corr_error'] = {}
    L1b_corr_attributes['Corr_error']['long_name'] = 'Correction Error Flag'
    L1b_corr_attributes['Corr_error']['description'] = ('Shows if a correction '
        'algorithm returned an error when called. See table 2.3.3-6 of the '
        '"L1b Products Format Specification" document')
    L1b_corr_attributes['Corr_error']['hertz'] = 1

    #-- CryoSat-2 Average Waveforms Groups
    #-- Low-Resolution Mode (LRM)
    L1b_1Hz_LRM_wfm_attributes = {}
    #-- Data Record Time (MDSR Time Stamp)
    L1b_1Hz_LRM_wfm_attributes['Day'] = {}
    L1b_1Hz_LRM_wfm_attributes['Day']['long_name'] = 'MDSR time stamp days'
    L1b_1Hz_LRM_wfm_attributes['Day']['units'] = 'days since 2000-01-01 00:00:00 TAI'
    L1b_1Hz_LRM_wfm_attributes['Day']['description'] = ('Corresponding to '
        'the middle of group of pulses')
    L1b_1Hz_LRM_wfm_attributes['Day']['hertz'] = 1
    #-- Time: second part
    L1b_1Hz_LRM_wfm_attributes['Second'] = {}
    L1b_1Hz_LRM_wfm_attributes['Second']['long_name'] = 'MDSR time stamp seconds'
    L1b_1Hz_LRM_wfm_attributes['Second']['units'] = 'seconds'
    L1b_1Hz_LRM_wfm_attributes['Second']['description'] = ('Corresponding to '
        'the middle of group of pulses')
    L1b_1Hz_LRM_wfm_attributes['Second']['hertz'] = 1
    #-- Time: microsecond part
    L1b_1Hz_LRM_wfm_attributes['Micsec'] = {}
    L1b_1Hz_LRM_wfm_attributes['Micsec']['long_name'] = 'MDSR time stamp microseconds'
    L1b_1Hz_LRM_wfm_attributes['Micsec']['units'] = 'microseconds'
    L1b_1Hz_LRM_wfm_attributes['Micsec']['description'] = ('Corresponding '
        'to the middle of group of pulses')
    L1b_1Hz_LRM_wfm_attributes['Micsec']['hertz'] = 1
    #-- Lat: packed units (0.1 micro-degree, 1e-7 degrees)
    L1b_1Hz_LRM_wfm_attributes['Lat'] = {}
    L1b_1Hz_LRM_wfm_attributes['Lat']['long_name'] = 'Latitude of measurement'
    L1b_1Hz_LRM_wfm_attributes['Lat']['description'] = ('Corresponding to the '
        'position at the MDSR Time Stamp')
    L1b_1Hz_LRM_wfm_attributes['Lat']['units'] = '0.1 micro-degree'
    L1b_1Hz_LRM_wfm_attributes['Lat']['valid_min'] = -9e8
    L1b_1Hz_LRM_wfm_attributes['Lat']['valid_max'] = 9e8
    L1b_1Hz_LRM_wfm_attributes['Lat']['hertz'] = 1
    #-- Lon: packed units (0.1 micro-degree, 1e-7 degrees)
    L1b_1Hz_LRM_wfm_attributes['Lon'] = {}
    L1b_1Hz_LRM_wfm_attributes['Lon']['long_name'] = 'Longitude of measurement'
    L1b_1Hz_LRM_wfm_attributes['Lon']['description'] = ('Corresponding to the '
        'position at the MDSR Time Stamp')
    L1b_1Hz_LRM_wfm_attributes['Lon']['units'] = '0.1 micro-degree'
    L1b_1Hz_LRM_wfm_attributes['Lon']['valid_min'] = -18e8
    L1b_1Hz_LRM_wfm_attributes['Lon']['valid_max'] = 18e8
    L1b_1Hz_LRM_wfm_attributes['Lon']['hertz'] = 1
    #-- Alt: packed units (mm, 1e-3 m)
    #-- Altitude of COG above reference ellipsoid
    L1b_1Hz_LRM_wfm_attributes['Alt'] = {}
    L1b_1Hz_LRM_wfm_attributes['Alt']['long_name'] = 'Altitude'
    L1b_1Hz_LRM_wfm_attributes['Alt']['description'] = ('Altitude of Satellite '
        'COG above reference ellipsoid corresponding to the MDSR Time Stamp')
    L1b_1Hz_LRM_wfm_attributes['Alt']['units'] = 'millimeters'
    L1b_1Hz_LRM_wfm_attributes['Alt']['hertz'] = 1
    #-- Window Delay (two-way) corrected for instrument delays
    L1b_1Hz_LRM_wfm_attributes['TD'] = {}
    L1b_1Hz_LRM_wfm_attributes['TD']['long_name'] = 'Altitude'
    L1b_1Hz_LRM_wfm_attributes['TD']['description'] = ('Window Delay '
        '(two-way) from the telemetry corrected for instrument delays')
    L1b_1Hz_LRM_wfm_attributes['TD']['units'] = 'picoseconds'
    L1b_1Hz_LRM_wfm_attributes['TD']['hertz'] = 1
    #-- 1 Hz Averaged Power Echo Waveform
    L1b_1Hz_LRM_wfm_attributes['Waveform'] = {}
    L1b_1Hz_LRM_wfm_attributes['Waveform']['long_name'] = 'Averaged Power Echo'
    L1b_1Hz_LRM_wfm_attributes['Waveform']['description'] = ('Array of 128 bins. '
        'Averaged from all individual L0 echoes in approx 1 second (20 for LRM). '
        'Converted to Watts by using the scaling parameters. Power in Watts = '
        'counts*(A*1e-9)*2^B')
    L1b_1Hz_LRM_wfm_attributes['Waveform']['valid_min'] = 0
    L1b_1Hz_LRM_wfm_attributes['Waveform']['valid_max'] = 65535
    L1b_1Hz_LRM_wfm_attributes['Waveform']['hertz'] = 1
    #-- Echo Scale Factor (to scale echo to watts)
    L1b_1Hz_LRM_wfm_attributes['Linear_Wfm_Multiplier'] = {}
    L1b_1Hz_LRM_wfm_attributes['Linear_Wfm_Multiplier']['long_name'] = ('Echo '
        'Scale factor "A"')
    L1b_1Hz_LRM_wfm_attributes['Linear_Wfm_Multiplier']['description'] = ('"A" '
        'Parameter to scale echo to Watts. Power in Watts = counts*(A*1e-9)*2^B')
    L1b_1Hz_LRM_wfm_attributes['Linear_Wfm_Multiplier']['hertz'] = 1
    #-- Echo Scale Power (a power of 2 to scale echo to Watts)
    L1b_1Hz_LRM_wfm_attributes['Power2_Wfm_Multiplier'] = {}
    L1b_1Hz_LRM_wfm_attributes['Power2_Wfm_Multiplier']['long_name'] = ('Echo '
        'Scale power "B"')
    L1b_1Hz_LRM_wfm_attributes['Power2_Wfm_Multiplier']['description'] = ('"B" '
        'Parameter to scale echo to Watts. Power in Watts = counts*(A*1e-9)*2^B')
    L1b_1Hz_LRM_wfm_attributes['Power2_Wfm_Multiplier']['hertz'] = 1
    #-- Number of echoes averaged
    L1b_1Hz_LRM_wfm_attributes['N_avg_echoes'] = {}
    L1b_1Hz_LRM_wfm_attributes['N_avg_echoes']['long_name'] = ('Number of'
        'echoes averaged')
    L1b_1Hz_LRM_wfm_attributes['N_avg_echoes']['description'] = ('Normally '
        '1820 for LRM (= 91 averaged on-board *20))')
    L1b_1Hz_LRM_wfm_attributes['N_avg_echoes']['hertz'] = 1
    #-- Flags
    L1b_1Hz_LRM_wfm_attributes['Flags'] = {}
    L1b_1Hz_LRM_wfm_attributes['Flags']['long_name'] = 'Flags'
    L1b_1Hz_LRM_wfm_attributes['Flags']['description'] = ('For errors or '
        'information about echoes. See table 2.3.4-3a of the "L1b Products '
        'Format Specification" document')
    L1b_1Hz_LRM_wfm_attributes['Flags']['flag_meanings'] = \
        '1_hz_echo_error_not_computed mispointing_bad_angles'
    L1b_1Hz_LRM_wfm_attributes['Flags']['hertz'] = 1

    #-- SAR Mode
    L1b_1Hz_SAR_wfm_attributes = {}
    #-- Data Record Time (MDSR Time Stamp)
    L1b_1Hz_SAR_wfm_attributes['Day'] = {}
    L1b_1Hz_SAR_wfm_attributes['Day']['long_name'] = 'MDSR time stamp days'
    L1b_1Hz_SAR_wfm_attributes['Day']['units'] = 'days since 2000-01-01 00:00:00 TAI'
    L1b_1Hz_SAR_wfm_attributes['Day']['description'] = ('Corresponding to '
        'ground bounce time of the individual pulse')
    L1b_1Hz_SAR_wfm_attributes['Day']['hertz'] = 1
    #-- Time: second part
    L1b_1Hz_SAR_wfm_attributes['Second'] = {}
    L1b_1Hz_SAR_wfm_attributes['Second']['long_name'] = 'MDSR time stamp seconds'
    L1b_1Hz_SAR_wfm_attributes['Second']['units'] = 'seconds'
    L1b_1Hz_SAR_wfm_attributes['Second']['description'] = ('Corresponding to '
        'ground bounce time of the individual pulse')
    L1b_1Hz_SAR_wfm_attributes['Second']['hertz'] = 1
    #-- Time: microsecond part
    L1b_1Hz_SAR_wfm_attributes['Micsec'] = {}
    L1b_1Hz_SAR_wfm_attributes['Micsec']['long_name'] = 'MDSR time stamp microseconds'
    L1b_1Hz_SAR_wfm_attributes['Micsec']['units'] = 'microseconds'
    L1b_1Hz_SAR_wfm_attributes['Micsec']['description'] = ('Corresponding '
        'ground bounce time of the individual pulse')
    L1b_1Hz_SAR_wfm_attributes['Micsec']['hertz'] = 1
    #-- Lat: packed units (0.1 micro-degree, 1e-7 degrees)
    L1b_1Hz_SAR_wfm_attributes['Lat'] = {}
    L1b_1Hz_SAR_wfm_attributes['Lat']['long_name'] = 'Latitude of measurement'
    L1b_1Hz_SAR_wfm_attributes['Lat']['description'] = ('Corresponding to the '
        'position at the MDSR Time Stamp')
    L1b_1Hz_SAR_wfm_attributes['Lat']['units'] = '0.1 micro-degree'
    L1b_1Hz_SAR_wfm_attributes['Lat']['valid_min'] = -9e8
    L1b_1Hz_SAR_wfm_attributes['Lat']['valid_max'] = 9e8
    L1b_1Hz_SAR_wfm_attributes['Lat']['hertz'] = 1
    #-- Lon: packed units (0.1 micro-degree, 1e-7 degrees)
    L1b_1Hz_SAR_wfm_attributes['Lon'] = {}
    L1b_1Hz_SAR_wfm_attributes['Lon']['long_name'] = 'Longitude of measurement'
    L1b_1Hz_SAR_wfm_attributes['Lon']['description'] = ('Corresponding to the '
        'position at the MDSR Time Stamp')
    L1b_1Hz_SAR_wfm_attributes['Lon']['units'] = '0.1 micro-degree'
    L1b_1Hz_SAR_wfm_attributes['Lon']['valid_min'] = -18e8
    L1b_1Hz_SAR_wfm_attributes['Lon']['valid_max'] = 18e8
    L1b_1Hz_SAR_wfm_attributes['Lon']['hertz'] = 1
    #-- Alt: packed units (mm, 1e-3 m)
    #-- Altitude of COG above reference ellipsoid
    L1b_1Hz_SAR_wfm_attributes['Alt'] = {}
    L1b_1Hz_SAR_wfm_attributes['Alt']['long_name'] = 'Altitude'
    L1b_1Hz_SAR_wfm_attributes['Alt']['description'] = ('Altitude of Satellite '
        'COG above reference ellipsoid corresponding to the MDSR Time Stamp')
    L1b_1Hz_SAR_wfm_attributes['Alt']['units'] = 'millimeters'
    L1b_1Hz_SAR_wfm_attributes['Alt']['hertz'] = 1
    #-- Window Delay (two-way) corrected for instrument delays
    L1b_1Hz_SAR_wfm_attributes['TD'] = {}
    L1b_1Hz_SAR_wfm_attributes['TD']['long_name'] = 'Altitude'
    L1b_1Hz_SAR_wfm_attributes['TD']['description'] = ('Window Delay '
        '(two-way) from the telemetry corrected for instrument delays')
    L1b_1Hz_SAR_wfm_attributes['TD']['units'] = 'picoseconds'
    L1b_1Hz_SAR_wfm_attributes['TD']['hertz'] = 1
    #-- 1 Hz Averaged Power Echo Waveform
    L1b_1Hz_SAR_wfm_attributes['Waveform'] = {}
    L1b_1Hz_SAR_wfm_attributes['Waveform']['long_name'] = 'Averaged Power Echo'
    L1b_1Hz_SAR_wfm_attributes['Waveform']['description'] = ('Array of 128 bins. '
        'Averaged from all individual L0 echoes in approx 1 second (5120 for SAR). '
        'Converted to Watts by using the scaling parameters. Power in Watts = '
        'counts*(A*1e-9)*2^B. The last 1Hz average waveform of the product '
        'is meaningless in most of the cases because there are not enough FBR '
        'samples to be used in the averaging operation. When this happens the '
        'waveform is flagged as invalid')
    L1b_1Hz_SAR_wfm_attributes['Waveform']['valid_min'] = 0
    L1b_1Hz_SAR_wfm_attributes['Waveform']['valid_max'] = 65535
    L1b_1Hz_SAR_wfm_attributes['Waveform']['hertz'] = 1
    #-- Echo Scale Factor (to scale echo to watts)
    L1b_1Hz_SAR_wfm_attributes['Linear_Wfm_Multiplier'] = {}
    L1b_1Hz_SAR_wfm_attributes['Linear_Wfm_Multiplier']['long_name'] = ('Echo '
        'Scale factor "A"')
    L1b_1Hz_SAR_wfm_attributes['Linear_Wfm_Multiplier']['description'] = ('"A" '
        'Parameter to scale echo to Watts. Power in Watts = counts*(A*1e-9)*2^B')
    L1b_1Hz_SAR_wfm_attributes['Linear_Wfm_Multiplier']['hertz'] = 1
    #-- Echo Scale Power (a power of 2 to scale echo to Watts)
    L1b_1Hz_SAR_wfm_attributes['Power2_Wfm_Multiplier'] = {}
    L1b_1Hz_SAR_wfm_attributes['Power2_Wfm_Multiplier']['long_name'] = ('Echo '
        'Scale power "B"')
    L1b_1Hz_SAR_wfm_attributes['Power2_Wfm_Multiplier']['description'] = ('"B" '
        'Parameter to scale echo to Watts. Power in Watts = counts*(A*1e-9)*2^B')
    L1b_1Hz_SAR_wfm_attributes['Power2_Wfm_Multiplier']['hertz'] = 1
    #-- Number of echoes averaged
    L1b_1Hz_SAR_wfm_attributes['N_avg_echoes'] = {}
    L1b_1Hz_SAR_wfm_attributes['N_avg_echoes']['long_name'] = ('Number of'
        'echoes averaged')
    L1b_1Hz_SAR_wfm_attributes['N_avg_echoes']['description'] = ('Normally 5120'
        ' for SAR. May be lower if individual echoes are missing or rejected')
    L1b_1Hz_SAR_wfm_attributes['N_avg_echoes']['hertz'] = 1
    #-- Flags
    L1b_1Hz_SAR_wfm_attributes['Flags'] = {}
    L1b_1Hz_SAR_wfm_attributes['Flags']['long_name'] = 'Flags'
    L1b_1Hz_SAR_wfm_attributes['Flags']['description'] = ('For errors or '
        'information about echoes. See table 2.3.4-3b of the "L1b Products '
        'Format Specification" document')
    L1b_1Hz_SAR_wfm_attributes['Flags']['flag_meanings'] = \
        '1_hz_echo_error_not_computed mispointing_bad_angles'
    L1b_1Hz_SAR_wfm_attributes['Flags']['hertz'] = 1

    #-- SARIN Mode
    #-- Same as the LRM/SAR groups but the waveform array is 512 bins instead of
    #-- 128 and the number of echoes averaged is different.
    L1b_1Hz_SARIN_wfm_attributes = {}
    #-- Data Record Time (MDSR Time Stamp)
    L1b_1Hz_SARIN_wfm_attributes['Day'] = {}
    L1b_1Hz_SARIN_wfm_attributes['Day']['long_name'] = 'MDSR time stamp days'
    L1b_1Hz_SARIN_wfm_attributes['Day']['units'] = 'days since 2000-01-01 00:00:00 TAI'
    L1b_1Hz_SARIN_wfm_attributes['Day']['description'] = ('Corresponding to '
        'ground bounce time of the individual pulse')
    L1b_1Hz_SARIN_wfm_attributes['Day']['hertz'] = 1
    #-- Time: second part
    L1b_1Hz_SARIN_wfm_attributes['Second'] = {}
    L1b_1Hz_SARIN_wfm_attributes['Second']['long_name'] = 'MDSR time stamp seconds'
    L1b_1Hz_SARIN_wfm_attributes['Second']['units'] = 'seconds'
    L1b_1Hz_SARIN_wfm_attributes['Second']['description'] = ('Corresponding to '
        'ground bounce time of the individual pulse')
    L1b_1Hz_SARIN_wfm_attributes['Second']['hertz'] = 1
    #-- Time: microsecond part
    L1b_1Hz_SARIN_wfm_attributes['Micsec'] = {}
    L1b_1Hz_SARIN_wfm_attributes['Micsec']['long_name'] = 'MDSR time stamp microseconds'
    L1b_1Hz_SARIN_wfm_attributes['Micsec']['units'] = 'microseconds'
    L1b_1Hz_SARIN_wfm_attributes['Micsec']['description'] = ('Corresponding '
        'ground bounce time of the individual pulse')
    L1b_1Hz_SARIN_wfm_attributes['Micsec']['hertz'] = 1
    #-- Lat: packed units (0.1 micro-degree, 1e-7 degrees)
    L1b_1Hz_SARIN_wfm_attributes['Lat'] = {}
    L1b_1Hz_SARIN_wfm_attributes['Lat']['long_name'] = 'Latitude of measurement'
    L1b_1Hz_SARIN_wfm_attributes['Lat']['description'] = ('Corresponding to the '
        'position at the MDSR Time Stamp')
    L1b_1Hz_SARIN_wfm_attributes['Lat']['units'] = '0.1 micro-degree'
    L1b_1Hz_SARIN_wfm_attributes['Lat']['valid_min'] = -9e8
    L1b_1Hz_SARIN_wfm_attributes['Lat']['valid_max'] = 9e8
    L1b_1Hz_SARIN_wfm_attributes['Lat']['hertz'] = 1
    #-- Lon: packed units (0.1 micro-degree, 1e-7 degrees)
    L1b_1Hz_SARIN_wfm_attributes['Lon'] = {}
    L1b_1Hz_SARIN_wfm_attributes['Lon']['long_name'] = 'Longitude of measurement'
    L1b_1Hz_SARIN_wfm_attributes['Lon']['description'] = ('Corresponding to the '
        'position at the MDSR Time Stamp')
    L1b_1Hz_SARIN_wfm_attributes['Lon']['units'] = '0.1 micro-degree'
    L1b_1Hz_SARIN_wfm_attributes['Lon']['valid_min'] = -18e8
    L1b_1Hz_SARIN_wfm_attributes['Lon']['valid_max'] = 18e8
    L1b_1Hz_SARIN_wfm_attributes['Lon']['hertz'] = 1
    #-- Alt: packed units (mm, 1e-3 m)
    #-- Altitude of COG above reference ellipsoid
    L1b_1Hz_SARIN_wfm_attributes['Alt'] = {}
    L1b_1Hz_SARIN_wfm_attributes['Alt']['long_name'] = 'Altitude'
    L1b_1Hz_SARIN_wfm_attributes['Alt']['description'] = ('Altitude of Satellite '
        'COG above reference ellipsoid corresponding to the MDSR Time Stamp')
    L1b_1Hz_SARIN_wfm_attributes['Alt']['units'] = 'millimeters'
    L1b_1Hz_SARIN_wfm_attributes['Alt']['hertz'] = 1
    #-- Window Delay (two-way) corrected for instrument delays
    L1b_1Hz_SARIN_wfm_attributes['TD'] = {}
    L1b_1Hz_SARIN_wfm_attributes['TD']['long_name'] = 'Altitude'
    L1b_1Hz_SARIN_wfm_attributes['TD']['description'] = ('Window Delay '
        '(two-way) from the telemetry corrected for instrument delays')
    L1b_1Hz_SARIN_wfm_attributes['TD']['units'] = 'picoseconds'
    L1b_1Hz_SARIN_wfm_attributes['TD']['hertz'] = 1
    #-- 1 Hz Averaged Power Echo Waveform
    L1b_1Hz_SARIN_wfm_attributes['Waveform'] = {}
    L1b_1Hz_SARIN_wfm_attributes['Waveform']['long_name'] = 'Averaged Power Echo'
    L1b_1Hz_SARIN_wfm_attributes['Waveform']['description'] = ('Array of 512 bins. '
        'Averaged from all individual 1280 L0 echoes in SARin mode. '
        'Converted to Watts by using the scaling parameters. Power in Watts = '
        'counts*(A*1e-9)*2^B. The last 1Hz average waveform of the product '
        'is meaningless in most of the cases because there are not enough FBR '
        'samples to be used in the averaging operation. When this happens the '
        'waveform is flagged as invalid')
    L1b_1Hz_SARIN_wfm_attributes['Waveform']['valid_min'] = 0
    L1b_1Hz_SARIN_wfm_attributes['Waveform']['valid_max'] = 65535
    L1b_1Hz_SARIN_wfm_attributes['Waveform']['hertz'] = 1
    #-- Echo Scale Factor (to scale echo to watts)
    L1b_1Hz_SARIN_wfm_attributes['Linear_Wfm_Multiplier'] = {}
    L1b_1Hz_SARIN_wfm_attributes['Linear_Wfm_Multiplier']['long_name'] = ('Echo '
        'Scale factor "A"')
    L1b_1Hz_SARIN_wfm_attributes['Linear_Wfm_Multiplier']['description'] = ('"A" '
        'Parameter to scale echo to Watts. Power in Watts = counts*(A*1e-9)*2^B')
    L1b_1Hz_SARIN_wfm_attributes['Linear_Wfm_Multiplier']['hertz'] = 1
    #-- Echo Scale Power (a power of 2 to scale echo to Watts)
    L1b_1Hz_SARIN_wfm_attributes['Power2_Wfm_Multiplier'] = {}
    L1b_1Hz_SARIN_wfm_attributes['Power2_Wfm_Multiplier']['long_name'] = ('Echo '
        'Scale power "B"')
    L1b_1Hz_SARIN_wfm_attributes['Power2_Wfm_Multiplier']['description'] = ('"B" '
        'Parameter to scale echo to Watts. Power in Watts = counts*(A*1e-9)*2^B')
    L1b_1Hz_SARIN_wfm_attributes['Power2_Wfm_Multiplier']['hertz'] = 1
    #-- Number of echoes averaged
    L1b_1Hz_SARIN_wfm_attributes['N_avg_echoes'] = {}
    L1b_1Hz_SARIN_wfm_attributes['N_avg_echoes']['long_name'] = ('Number of'
        'echoes averaged')
    L1b_1Hz_SARIN_wfm_attributes['N_avg_echoes']['description'] = ('Normally 1280'
        ' for SARIN. May be lower if individual echoes are missing or rejected')
    L1b_1Hz_SARIN_wfm_attributes['N_avg_echoes']['hertz'] = 1
    #-- Flags
    L1b_1Hz_SARIN_wfm_attributes['Flags'] = {}
    L1b_1Hz_SARIN_wfm_attributes['Flags']['long_name'] = 'Flags'
    L1b_1Hz_SARIN_wfm_attributes['Flags']['description'] = ('For errors or '
        'information about echoes. See table 2.3.4-3b of the "L1b Products '
        'Format Specification" document')
    L1b_1Hz_SARIN_wfm_attributes['Flags']['flag_meanings'] = \
        '1_hz_echo_error_not_computed mispointing_bad_angles'
    L1b_1Hz_SARIN_wfm_attributes['Flags']['hertz'] = 1

    #-- CryoSat-2 Waveforms Groups
    #-- Beam Behavior Parameters
    L1b_BB_attributes = {}
    #-- Standard Deviation of Gaussian fit to range integrated stack power.
    L1b_BB_attributes['SD'] = {}
    L1b_BB_attributes['SD']['long_name'] = 'Standard Deviation'
    L1b_BB_attributes['SD']['description'] = ('Standard Deviation of Gaussian '
        'fit to range integrated stack power')
    L1b_BB_attributes['SD']['units'] = '1/100'
    L1b_BB_attributes['SD']['hertz'] = 20
    #-- Stack Center: Mean of Gaussian fit to range integrated stack power.
    L1b_BB_attributes['Center'] = {}
    L1b_BB_attributes['Center']['long_name'] = 'Stack Center'
    L1b_BB_attributes['Center']['description'] = ('Mean of Gaussian fit to '
        'range integrated stack power')
    L1b_BB_attributes['Center']['units'] = '1/100'
    L1b_BB_attributes['Center']['hertz'] = 20
    #-- Stack amplitude parameter scaled in dB/100.
    L1b_BB_attributes['Amplitude'] = {}
    L1b_BB_attributes['Amplitude']['long_name'] = 'Stack amplitude'
    L1b_BB_attributes['Amplitude']['description'] = 'Amplitude Parameter'
    L1b_BB_attributes['Amplitude']['units'] = 'dB/100'
    L1b_BB_attributes['Amplitude']['hertz'] = 20
    #-- 3rd moment: providing the degree of asymmetry of the range integrated
    #-- stack power distribution.
    L1b_BB_attributes['Skewness'] = {}
    L1b_BB_attributes['Skewness']['long_name'] = 'Stack Skewness'
    L1b_BB_attributes['Skewness']['description'] = ('3rd moment: providing the '
        'degree of asymmetry of the range integrated stack power distribution')
    L1b_BB_attributes['Skewness']['units'] = '1/100'
    L1b_BB_attributes['Skewness']['_FillValue'] = -99900
    L1b_BB_attributes['Skewness']['hertz'] = 20
    #-- 4th moment: Measure of peakiness of range integrated stack power distribution.
    L1b_BB_attributes['Kurtosis'] = {}
    L1b_BB_attributes['Kurtosis']['long_name'] = 'Stack Kurtosis'
    L1b_BB_attributes['Kurtosis']['description'] = ('4th moment: measure of '
        'peakiness of range integrated stack power distribution')
    L1b_BB_attributes['Kurtosis']['units'] = '1/100'
    L1b_BB_attributes['Kurtosis']['_FillValue'] = -99900
    L1b_BB_attributes['Kurtosis']['hertz'] = 20
    #-- Stack peakiness computed from the range integrated power of the single look echoes
    L1b_BB_attributes['Peakiness'] = {}
    L1b_BB_attributes['Peakiness']['long_name'] = 'Stack Peakiness'
    L1b_BB_attributes['Peakiness']['description'] = ('Stack peakiness computed '
        'from the range integrated power of the single look echoes within a '
        'stack. Stack peakiness is defined as the inverse of the average of '
        'the range integrated power normalized for the power at zero look angle')
    L1b_BB_attributes['Peakiness']['units'] = '1/100'
    L1b_BB_attributes['Peakiness']['_FillValue'] = -99900
    L1b_BB_attributes['Peakiness']['hertz'] = 20
    #-- Stack residuals of Gaussian that fits the range integrated power of the single look echoes
    L1b_BB_attributes['RMS'] = {}
    L1b_BB_attributes['RMS']['long_name'] = 'Gaussian Power Fit Residuals'
    L1b_BB_attributes['RMS']['description'] = ('Residuals of Gaussian that '
        'fits the range integrated power of the single look echoes within a '
        'stack. It is the root mean squared error between the Gaussian fitting '
        'and the range integrated power of the single look echoes within a stack')
    L1b_BB_attributes['RMS']['units'] = 'dbW'
    L1b_BB_attributes['RMS']['_FillValue'] = -99900
    L1b_BB_attributes['RMS']['hertz'] = 20
    #-- Standard deviation as a function of boresight angle (microradians)
    L1b_BB_attributes['SD_boresight_angle'] = {}
    L1b_BB_attributes['SD_boresight_angle']['long_name'] = 'Standard Deviation'
    L1b_BB_attributes['SD_boresight_angle']['description'] = ('Standard '
        'deviation as a function of boresight angle')
    L1b_BB_attributes['SD_boresight_angle']['units'] = 'microradians'
    L1b_BB_attributes['SD_boresight_angle']['valid_min'] = 0
    L1b_BB_attributes['SD_boresight_angle']['valid_max'] = 65525
    L1b_BB_attributes['SD_boresight_angle']['hertz'] = 20
    #-- Stack Center angle as a function of boresight angle (microradians)
    L1b_BB_attributes['Center_boresight_angle'] = {}
    L1b_BB_attributes['Center_boresight_angle']['long_name'] = 'Stack Centre Angle'
    L1b_BB_attributes['Center_boresight_angle']['description'] = ('Stack Centre '
        'angle as a function of boresight angle')
    L1b_BB_attributes['Center_boresight_angle']['units'] = 'microradians'
    L1b_BB_attributes['Center_boresight_angle']['valid_min'] = -32768
    L1b_BB_attributes['Center_boresight_angle']['valid_max'] = 32768
    L1b_BB_attributes['Center_boresight_angle']['hertz'] = 20
    #-- Stack Center angle as a function of look angle (microradians)
    L1b_BB_attributes['Center_look_angle'] = {}
    L1b_BB_attributes['Center_look_angle']['long_name'] = 'Stack Centre Angle'
    L1b_BB_attributes['Center_look_angle']['description'] = ('Stack Centre '
        'angle as a function of look angle')
    L1b_BB_attributes['Center_look_angle']['units'] = 'microradians'
    L1b_BB_attributes['Center_look_angle']['valid_min'] = -32768
    L1b_BB_attributes['Center_look_angle']['valid_max'] = 32768
    L1b_BB_attributes['Center_look_angle']['hertz'] = 20
    #-- Number of contributing beams in the stack before weighting
    L1b_BB_attributes['Number'] = {}
    L1b_BB_attributes['Number']['long_name'] = ('Number of contributing '
        'beams before weighting')
    L1b_BB_attributes['Number']['description'] = ('Number of contributing '
        'beams in the stack before weighting: number of single look echoes '
        'in the stack before the Surface Sample Stack weighting is applied')
    L1b_BB_attributes['Number']['units'] = 'count'
    L1b_BB_attributes['Number']['hertz'] = 20
    #-- Number of contributing beams in the stack after weighting
    L1b_BB_attributes['Weighted_Number'] = {}
    L1b_BB_attributes['Weighted_Number']['long_name'] = ('Number of contributing '
        'beams after weighting')
    L1b_BB_attributes['Weighted_Number']['description'] = ('Number of contributing '
        'beams in the stack after weighting: number of single look echoes '
        'in the stack after the Surface Sample Stack weighting is applied')
    L1b_BB_attributes['Weighted_Number']['units'] = 'count'
    L1b_BB_attributes['Weighted_Number']['hertz'] = 20

    #-- Low-Resolution Mode
    L1b_20Hz_LRM_wfm_attributes = {}
    #-- Averaged Power Echo Waveform
    L1b_20Hz_LRM_wfm_attributes['Waveform'] = {}
    L1b_20Hz_LRM_wfm_attributes['Waveform']['long_name'] = 'Averaged Power Echo'
    L1b_20Hz_LRM_wfm_attributes['Waveform']['description'] = ('Array of 128 bins. '
        'Averaged (on-board) from 91 individual in pulse limited mode (LRM). '
        'Converted to Watts by using the scaling parameters. Power in Watts = '
        'counts*(A*1e-9)*2^B')
    L1b_20Hz_LRM_wfm_attributes['Waveform']['valid_min'] = 0
    L1b_20Hz_LRM_wfm_attributes['Waveform']['valid_max'] = 65535
    L1b_20Hz_LRM_wfm_attributes['Waveform']['hertz'] = 20
    #-- Echo Scale Factor (to scale echo to watts)
    L1b_20Hz_LRM_wfm_attributes['Linear_Wfm_Multiplier'] = {}
    L1b_20Hz_LRM_wfm_attributes['Linear_Wfm_Multiplier']['long_name'] = ('Echo '
        'Scale factor "A"')
    L1b_20Hz_LRM_wfm_attributes['Linear_Wfm_Multiplier']['description'] = ('"A" '
        'Parameter to scale echo to Watts. Power in Watts = counts*(A*1e-9)*2^B')
    L1b_20Hz_LRM_wfm_attributes['Linear_Wfm_Multiplier']['hertz'] = 20
    #-- Echo Scale Power (a power of 2 to scale echo to Watts)
    L1b_20Hz_LRM_wfm_attributes['Power2_Wfm_Multiplier'] = {}
    L1b_20Hz_LRM_wfm_attributes['Power2_Wfm_Multiplier']['long_name'] = ('Echo '
        'Scale power "B"')
    L1b_20Hz_LRM_wfm_attributes['Power2_Wfm_Multiplier']['description'] = ('"B" '
        'Parameter to scale echo to Watts. Power in Watts = counts*(A*1e-9)*2^B')
    L1b_20Hz_LRM_wfm_attributes['Power2_Wfm_Multiplier']['hertz'] = 20
    #-- Number of echoes averaged
    L1b_20Hz_LRM_wfm_attributes['N_avg_echoes'] = {}
    L1b_20Hz_LRM_wfm_attributes['N_avg_echoes']['long_name'] = ('Number of'
        'echoes averaged')
    L1b_20Hz_LRM_wfm_attributes['N_avg_echoes']['description'] = ('Normally '
        '91 for LRM')
    L1b_20Hz_LRM_wfm_attributes['N_avg_echoes']['hertz'] = 20
    #-- Flags
    L1b_20Hz_LRM_wfm_attributes['Flags'] = {}
    L1b_20Hz_LRM_wfm_attributes['Flags']['long_name'] = 'Flags'
    L1b_20Hz_LRM_wfm_attributes['Flags']['description'] = ('TRK cycle report '
        '(as extracted from the L0). See table 2.3.4-4a of the "L1b Products '
        'Format Specification" document')
    L1b_20Hz_LRM_wfm_attributes['Flags']['hertz'] = 20

    #-- SAR Mode
    L1b_20Hz_SAR_wfm_attributes = {}
    #-- Averaged Power Echo Waveform
    L1b_20Hz_SAR_wfm_attributes['Waveform'] = {}
    L1b_20Hz_SAR_wfm_attributes['Waveform']['long_name'] = 'Averaged Power Echo'
    L1b_20Hz_SAR_wfm_attributes['Waveform']['description'] = ('Array of 256 bins. '
        'Averaged from a set of Doppler beam echoes formed at a common surface '
        'location. Converted to Watts by using the scaling parameters. '
        'Power in Watts = counts*(A*1e-9)*2^B.')
    L1b_20Hz_SAR_wfm_attributes['Waveform']['valid_min'] = 0
    L1b_20Hz_SAR_wfm_attributes['Waveform']['valid_max'] = 65535
    L1b_20Hz_SAR_wfm_attributes['Waveform']['hertz'] = 20
    #-- Echo Scale Factor (to scale echo to watts)
    L1b_20Hz_SAR_wfm_attributes['Linear_Wfm_Multiplier'] = {}
    L1b_20Hz_SAR_wfm_attributes['Linear_Wfm_Multiplier']['long_name'] = ('Echo '
        'Scale factor "A"')
    L1b_20Hz_SAR_wfm_attributes['Linear_Wfm_Multiplier']['description'] = ('"A" '
        'Parameter to scale echo to Watts. Power in Watts = counts*(A*1e-9)*2^B')
    L1b_20Hz_SAR_wfm_attributes['Linear_Wfm_Multiplier']['hertz'] = 20
    #-- Echo Scale Power (a power of 2 to scale echo to Watts)
    L1b_20Hz_SAR_wfm_attributes['Power2_Wfm_Multiplier'] = {}
    L1b_20Hz_SAR_wfm_attributes['Power2_Wfm_Multiplier']['long_name'] = ('Echo '
        'Scale power "B"')
    L1b_20Hz_SAR_wfm_attributes['Power2_Wfm_Multiplier']['description'] = ('"B" '
        'Parameter to scale echo to Watts. Power in Watts = counts*(A*1e-9)*2^B')
    L1b_20Hz_SAR_wfm_attributes['Power2_Wfm_Multiplier']['hertz'] = 20
    #-- Number of echoes averaged
    L1b_20Hz_SAR_wfm_attributes['N_avg_echoes'] = {}
    L1b_20Hz_SAR_wfm_attributes['N_avg_echoes']['long_name'] = ('Number of'
        'echoes averaged')
    L1b_20Hz_SAR_wfm_attributes['N_avg_echoes']['description'] = ('Normally 280'
        ' for SAR')
    L1b_20Hz_SAR_wfm_attributes['N_avg_echoes']['hertz'] = 20
    #-- Flags
    L1b_20Hz_SAR_wfm_attributes['Flags'] = {}
    L1b_20Hz_SAR_wfm_attributes['Flags']['long_name'] = 'Flags'
    L1b_20Hz_SAR_wfm_attributes['Flags']['description'] = ('For errors or '
        'information about echoes. See table 2.3.4-4b of the "L1b Products '
        'Format Specification" document')
    L1b_20Hz_SAR_wfm_attributes['Flags']['hertz'] = 20
    #-- Beam Behavior Parameters
    L1b_20Hz_SAR_wfm_attributes['Beam'] = L1b_BB_attributes

    #-- SARIN Mode
    L1b_20Hz_SARIN_wfm_attributes = {}
    #-- Averaged Power Echo Waveform
    L1b_20Hz_SARIN_wfm_attributes['Waveform'] = {}
    L1b_20Hz_SARIN_wfm_attributes['Waveform']['long_name'] = 'Averaged Power Echo'
    L1b_20Hz_SARIN_wfm_attributes['Waveform']['description'] = ('Array of 1024 bins. '
        'Averaged from 2 sets of Doppler beam echoes formed (on 2 receive '
        'channels) at a common surface location. Converted to Watts by using '
        'the scaling parameters. Power in Watts = counts*(A*1e-9)*2^B.')
    L1b_20Hz_SARIN_wfm_attributes['Waveform']['valid_min'] = 0
    L1b_20Hz_SARIN_wfm_attributes['Waveform']['valid_max'] = 65535
    L1b_20Hz_SARIN_wfm_attributes['Waveform']['hertz'] = 20
    #-- Echo Scale Factor (to scale echo to watts)
    L1b_20Hz_SARIN_wfm_attributes['Linear_Wfm_Multiplier'] = {}
    L1b_20Hz_SARIN_wfm_attributes['Linear_Wfm_Multiplier']['long_name'] = ('Echo '
        'Scale factor "A"')
    L1b_20Hz_SARIN_wfm_attributes['Linear_Wfm_Multiplier']['description'] = ('"A" '
        'Parameter to scale echo to Watts. Power in Watts = counts*(A*1e-9)*2^B')
    L1b_20Hz_SARIN_wfm_attributes['Linear_Wfm_Multiplier']['hertz'] = 20
    #-- Echo Scale Power (a power of 2 to scale echo to Watts)
    L1b_20Hz_SARIN_wfm_attributes['Power2_Wfm_Multiplier'] = {}
    L1b_20Hz_SARIN_wfm_attributes['Power2_Wfm_Multiplier']['long_name'] = ('Echo '
        'Scale power "B"')
    L1b_20Hz_SARIN_wfm_attributes['Power2_Wfm_Multiplier']['description'] = ('"B" '
        'Parameter to scale echo to Watts. Power in Watts = counts*(A*1e-9)*2^B')
    L1b_20Hz_SARIN_wfm_attributes['Power2_Wfm_Multiplier']['hertz'] = 20
    #-- Number of echoes averaged
    L1b_20Hz_SARIN_wfm_attributes['N_avg_echoes'] = {}
    L1b_20Hz_SARIN_wfm_attributes['N_avg_echoes']['long_name'] = ('Number of'
        'echoes averaged')
    L1b_20Hz_SARIN_wfm_attributes['N_avg_echoes']['description'] = ('Normally 70'
        ' for SARIN')
    L1b_20Hz_SARIN_wfm_attributes['N_avg_echoes']['hertz'] = 20
    #-- Flags
    L1b_20Hz_SARIN_wfm_attributes['Flags'] = {}
    L1b_20Hz_SARIN_wfm_attributes['Flags']['long_name'] = 'Flags'
    L1b_20Hz_SARIN_wfm_attributes['Flags']['description'] = ('For errors or '
        'information about echoes. See table 2.3.4-4b of the "L1b Products '
        'Format Specification" document')
    L1b_20Hz_SARIN_wfm_attributes['Flags']['hertz'] = 20
    #-- Coherence [1024]: packed units (1/1000)
    L1b_20Hz_SARIN_wfm_attributes['Coherence'] = {}
    L1b_20Hz_SARIN_wfm_attributes['Coherence']['long_name'] = 'Coherence'
    L1b_20Hz_SARIN_wfm_attributes['Coherence']['description'] = ('Array of '
        '1024 bins. Computed from the complex echoes on the 2 Rx channels')
    L1b_20Hz_SARIN_wfm_attributes['Coherence']['units'] = '1/1000'
    L1b_20Hz_SARIN_wfm_attributes['Coherence']['hertz'] = 20
    #-- Phase Difference [1024]: packed units (microradians)
    L1b_20Hz_SARIN_wfm_attributes['Phase_diff'] = {}
    L1b_20Hz_SARIN_wfm_attributes['Phase_diff']['long_name'] = 'Phase Difference'
    L1b_20Hz_SARIN_wfm_attributes['Phase_diff']['description'] = ('Array of '
        '1024 bins. Computed from the complex echoes on the 2 Rx channels')
    L1b_20Hz_SARIN_wfm_attributes['Phase_diff']['units'] = 'microradians'
    L1b_20Hz_SARIN_wfm_attributes['Phase_diff']['valid_min'] = -3141593
    L1b_20Hz_SARIN_wfm_attributes['Phase_diff']['valid_max'] = 3141593
    L1b_20Hz_SARIN_wfm_attributes['Phase_diff']['hertz'] = 20
    #-- Beam Behavior Parameters
    L1b_20Hz_SARIN_wfm_attributes['Beam'] = L1b_BB_attributes

    #-- Bind all the l1b attributes together into single dictionary
    CS_l1b_attrib = {}
    CS_l1b_attrib['Location'] = L1b_location_attributes
    CS_l1b_attrib['Data'] = L1b_measurement_attributes
    CS_l1b_attrib['Geometry'] = L1b_corr_attributes
    if (MODE == 'LRM'):
        CS_l1b_attrib['Waveform_1Hz'] = L1b_1Hz_LRM_wfm_attributes
        CS_l1b_attrib['Waveform_20Hz'] = L1b_20Hz_LRM_wfm_attributes
    elif (MODE == 'SAR'):
        CS_l1b_attrib['Waveform_1Hz'] = L1b_1Hz_SAR_wfm_attributes
        CS_l1b_attrib['Waveform_20Hz'] = L1b_20Hz_SAR_wfm_attributes
    elif (MODE == 'SIN'):
        CS_l1b_attrib['Waveform_1Hz'] = L1b_1Hz_SARIN_wfm_attributes
        CS_l1b_attrib['Waveform_20Hz'] = L1b_20Hz_SARIN_wfm_attributes

    #-- return the output dictionary
    return CS_l1b_attrib
