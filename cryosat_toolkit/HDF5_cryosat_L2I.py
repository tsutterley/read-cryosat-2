#!/usr/bin/env python
u"""
HDF5_cryosat_L2I.py (05/2022)
Reads and Writes HDF5 files for CryoSat-2 Level-2 intermediary data products
Supported CryoSat Modes: LRM, SAR, SARin, FDM, SID, GDR

OUTPUTS a formatted HDF5 file with:
    Location: Time and Orbit Parameters
    Geometry: Elevation Corrections and Flags
    Data: Geolocation and Elevation Measurements with Quality Parameters
    Auxiliary: Auxiliary Data for Elevation Processing
    Instrumental: Intrument Corrections
    METADATA: MPH, SPH and DSD Header data

OPTIONS:
    BASELINE (HDF5_cryosat_L2I): CryoSat-2 baseline (A, B, C)
    FILENAME (HDF5_cryosat_L2I): output HDF5 file name
    TITLE (HDF5_cryosat_L2I): output file description
    HEADER (HDF5_cryosat_L2I): output CryoSat-2 file headers (MPH, SPH, DSD)
        1: for single CryoSat-2 files
        2: for merged CryoSat-2 files from convert_cryosat_L2I.py
    CLOBBER (HDF5_cryosat_L2I): overwrite existing HDF5 file
    VERBOSE: print HDF5 structure parameters to screen
    ATTRIBUTES (read_HDF5_cryosat_L2I): input variable attributes from HDF5 file

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        (https://www.h5py.org/)

UPDATE HISTORY:
    Updated 05/2022: added docstrings in numpy documentation format
    Updated 08/2020: flake8 updates for python3
    Updated 02/2020: convert from hard to soft tabulation
    Updated 10/2019: changing Y/N flags to True/False
    Updated 04/2019: print HDF5 keys from list for python3 compatibility
    Updated 06/2018: use items instead of iteritems for python3 compatibility
    Updated 05/2016: using __future__ print function
    Written 04/2016
"""
from __future__ import print_function

import os
import re
import h5py

#-- PURPOSE: write CryoSat-2 HDF5 files
def HDF5_cryosat_L2I(CS_l2I_mds, MODE, BASELINE, FILENAME='', TITLE='',
    HEADER=0, CLOBBER=True, VERBOSE=False):
    #-- setting HDF5 clobber attribute
    if CLOBBER:
        clobber = 'w'
    else:
        clobber = 'w-'

    # #-- getting HDF5 dataset attributes for each variable
    # CS_l2I_attrib = cryosat_L2I_attributes(MODE, BASELINE)

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
    n_records,n_blocks = CS_l2I_mds['Location']['Day'].shape
    n_1Hz_wfm = CS_l2I_mds['Waveform_1Hz']['Waveform'].shape[1]
    n_20Hz_wfm = CS_l2I_mds['Waveform_20Hz']['Waveform'].shape[2]

    #-- find keys to output (do not output empty Spares variables)
    Location_keys = [key for key in CS_l2I_mds['Location'].keys() if not
        re.search('Spare',key)]
    Geometry_keys = [key for key in CS_l2I_mds['Geometry'].keys() if not
        re.search('Spare',key)]
    Data_keys = [key for key in CS_l2I_mds['Data'].keys() if not
        re.search('Spare',key)]
    Auxiliary_keys = [key for key in CS_l2I_mds['Auxiliary'].keys() if not
        re.search('Spare',key)]
    Instrumental_keys = [key for key in CS_l2I_mds['Instrumental'].keys() if not
        re.search('Spare',key)]

    #-- create HDF5 records
    h5 = {}
    h5['Location'] = {}
    h5['Data'] = {}
    h5['Auxiliary'] = {}
    h5['Geometry'] = {}
    h5['Instrumental'] = {}

    #-- CryoSat-2 Time and Orbit Group
    for key in Location_keys:
        val = CS_l2I_mds['Location'][key]
        if key in ('Sat_velocity','Real_beam','Baseline'):
            #-- Defining the HDF5 dataset variables
            h5['Location'][key] = fileID.create_dataset('Location/{0}'.format(key),
                (n_records,3,), data=val, dtype=val.dtype,
                compression='gzip')
            #-- attach dimensions
            h5['Location'][key].dims[0].label='CS_L2I_MDS_REC_SIZE'
            h5['Location'][key].dims[1].label='CS_L2I_MDS_VECTOR_SIZE'
        else:
            #-- Defining the HDF5 dataset variables
            h5['Location'][key] = fileID.create_dataset('Location/{0}'.format(key),
                (n_records,), data=val, dtype=val.dtype,
                compression='gzip')
            #-- attach dimensions
            h5['Location'][key].dims[0].label='CS_L2I_MDS_REC_SIZE'
        # #-- add HDF5 variable attributes
        # for att_name,att_val in CS_l2I_attrib['Location'][key].items():
        #     h5['Location'][key].attrs[att_name] = att_val

    #-- CryoSat-2 Measurement Group
    for key in Geometry_keys:
        val = CS_l2I_mds['Data'][key]
        if key in ('BB_parameter'):
            #-- Defining the HDF5 dataset Beam Behavior variables
            h5['Data'][key] = fileID.create_dataset('Data/{0}'.format(key),
                (n_records,50,), data=val, dtype=val.dtype,
                compression='gzip')
            #-- attach dimensions
            h5['Data'][key].dims[0].label='CS_L2I_MDS_REC_SIZE'
            h5['Data'][key].dims[1].label='CS_L2I_MDS_PARAMETER_SIZE'
        else:
            #-- Defining the HDF5 dataset variables
            h5['Data'][key] = fileID.create_dataset('Data/{0}'.format(key),
                (n_records,), data=val, dtype=val.dtype,
                compression='gzip')
            #-- attach dimensions
            h5['Data'][key].dims[0].label='CS_L2I_MDS_REC_SIZE'
        # #-- add HDF5 variable attributes
        # for att_name,att_val in CS_l2I_attrib['Data'][key].items():
        #     h5['Data'][key].attrs[att_name] = att_val

    #-- CryoSat-2 Auxiliary Data Group
    for key in Auxiliary_keys:
        val = CS_l2I_mds['Auxiliary'][key]
        #-- Defining the HDF5 dataset Beam Behavior variables
        h5['Auxiliary'][key] = fileID.create_dataset('Auxiliary/{0}'.format(key),
            (n_records,), data=val, dtype=val.dtype, compression='gzip')
        #-- attach dimensions
        h5['Auxiliary'][key].dims[0].label='CS_L2I_MDS_REC_SIZE'

    #-- CryoSat-2 External Corrections Group
    for key in Geometry_keys:
        val = CS_l2I_mds['Geometry'][key]
        #-- Defining the HDF5 dataset Beam Behavior variables
        h5['Geometry'][key] = fileID.create_dataset('Geometry/{0}'.format(key),
            (n_records,), data=val, dtype=val.dtype, compression='gzip')
        #-- attach dimensions
        h5['Geometry'][key].dims[0].label='CS_L2I_MDS_REC_SIZE'

    #-- CryoSat-2 Internal Corrections Group
    for key in Instrumental_keys:
        val = CS_l2I_mds['Instrumental'][key]
        #-- Defining the HDF5 dataset Beam Behavior variables
        h5['Instrumental'][key] = fileID.create_dataset('Instrumental/{0}'.format(key),
            (n_records,), data=val, dtype=val.dtype, compression='gzip')
        #-- attach dimensions
        h5['Instrumental'][key].dims[0].label='CS_L2I_MDS_REC_SIZE'

    #-- output MPH/SPH/DSD headers as group attributes
    if (HEADER == 1):
        #-- HEADER 1 is for single CryoSat-2 files
        fileID.create_group('METADATA')
        fileID['METADATA'].create_group('MPH')
        fileID['METADATA'].create_group('SPH')
        fileID['METADATA'].create_group('DSD')
        #-- Main Product Header (MPH) are all strings
        for att_name,att_val in CS_l2I_mds['METADATA']['MPH'].items():
            fileID['METADATA']['MPH'].attrs[att_name] = att_val
        #-- Specific Product Header (SPH) are both strings and dictionaries
        for att_name,att_val in CS_l2I_mds['METADATA']['SPH'].items():
            if isinstance(att_val,dict):
                #-- if att_val is dictionary
                fileID['METADATA']['SPH'].create_group(att_name)
                for ds_name,ds_val in att_val.items():
                    fileID['METADATA']['SPH'][att_name].attrs[ds_name] = ds_val
            elif isinstance(att_val,str) and att_name:
                #-- if att_val is string
                fileID['METADATA']['SPH'].attrs[att_name] = att_val
        #-- Data Set Descriptors (DSD) are all strings
        for att_name,att_val in CS_l2I_mds['METADATA']['DSD'].items():
            fileID['METADATA']['DSD'].attrs[att_name] = att_val
    elif (HEADER == 2):
        #-- HEADER 2 is for merged CryoSat-2 files from convert_cryosat_L2I.py
        fileID.create_group('METADATA')
        fileID['METADATA'].create_group('MPH')
        fileID['METADATA'].create_group('SPH')
        fileID['METADATA'].create_group('DSD')
        #-- Main Product Header (MPH) are all strings
        for fi in CS_l2I_mds['METADATA']['MPH'].keys():
            fileID['METADATA']['MPH'].create_group(fi)
            for att_name,att_val in CS_l2I_mds['METADATA']['MPH'][fi].items():
                fileID['METADATA']['MPH'][fi].attrs[att_name] = att_val
        #-- Specific Product Header (SPH) are both strings and dictionaries
        for fi in CS_l2I_mds['METADATA']['SPH'].keys():
            fileID['METADATA']['SPH'].create_group(fi)
            for att_name,att_val in CS_l2I_mds['METADATA']['SPH'][fi].items():
                if isinstance(att_val,dict):
                    #-- if att_val is dictionary
                    fileID['METADATA']['SPH'][fi].create_group(att_name)
                    for dsn,dsv in att_val.items():
                        fileID['METADATA']['SPH'][fi][att_name].attrs[dsn] = dsv
                elif isinstance(att_val,str) and att_name:
                    #-- if att_val is string
                    fileID['METADATA']['SPH'][fi].attrs[att_name] = att_val
        #-- Data Set Descriptors (DSD) are all strings
        for fi in CS_l2I_mds['METADATA']['DSD'].keys():
            fileID['METADATA']['DSD'].create_group(fi)
            for att_name,att_val in CS_l2I_mds['METADATA']['DSD'][fi].items():
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
def read_HDF5_cryosat_L2I(FILENAME, ATTRIBUTES=True, VERBOSE=False):
    #-- Open the HDF5 file for reading
    fileID = h5py.File(os.path.expanduser(FILENAME), 'r')

    #-- Output HDF5 file information
    if VERBOSE:
        print(fileID.filename)
        print(list(fileID.keys()))

    #-- allocate python dictionaries for output CS_l2I_mds variables
    CS_l2I_mds = {}
    CS_l2I_mds['Location'] = {}
    CS_l2I_mds['Data'] = {}
    CS_l2I_mds['Auxiliary'] = {}
    CS_l2I_mds['Geometry'] = {}
    CS_l2I_mds['Instrumental'] = {}

    #-- get each HDF5 variable
    #-- CryoSat-2 Location Group
    for key in fileID['Location'].keys():
        if key in ('Sat_velocity','Real_beam','Baseline'):
            CS_l2I_mds['Location'][key] = fileID['Location'][key][:,:]
        else:
            CS_l2I_mds['Location'][key] = fileID['Location'][key][:]
    #-- CryoSat-2 Measurement Group
    for key in fileID['Data'].keys():
        if key in ('BB_parameter'):
            CS_l2I_mds['Data'][key] = fileID['Data'][key][:,:]
        else:
            CS_l2I_mds['Location'][key] = fileID['Location'][key][:]
    #-- CryoSat-2 Auxiliary Data Group
    for key in fileID['Auxiliary'].keys():
        CS_l2I_mds['Auxiliary'][key] = fileID['Auxiliary'][key][:]
    #-- CryoSat-2 External Corrections Group
    for key in fileID['Geometry'].keys():
        CS_l2I_mds['Geometry'][key] = fileID['Geometry'][key][:]
    #-- CryoSat-2 Internal Corrections Group
    for key in fileID['Instrumental'].keys():
        CS_l2I_mds['Instrumental'][key] = fileID['Instrumental'][key][:]

    #-- Getting attributes of included variables
    if ATTRIBUTES:
        #-- allocate python dictionaries for output CS_l2I_mds attributes
        CS_l2I_mds['Attributes'] = {}
        CS_l2I_mds['Attributes']['Location'] = {}
        CS_l2I_mds['Attributes']['Data'] = {}
        CS_l2I_mds['Attributes']['Auxiliary'] = {}
        CS_l2I_mds['Attributes']['Geometry'] = {}
        CS_l2I_mds['Attributes']['Instrumental'] = {}
        #-- CryoSat-2 Location Group
        for key in fileID['Location'].keys():
            CS_l2I_mds['Attributes']['Location'][key] = {}
            for att_name,att_val in fileID['Location'][key].attrs.items():
                CS_l2I_mds['Attributes']['Location'][key][att_name] = att_val
        #-- CryoSat-2 Measurement Group
        for key in fileID['Data'].keys():
            CS_l2I_mds['Attributes']['Data'][key] = {}
            for att_name,att_val in fileID['Data'][key].attrs.items():
                CS_l2I_mds['Attributes']['Data'][key][att_name] = att_val
        #-- CryoSat-2 Auxiliary Data Group
        for key in fileID['Auxiliary'].keys():
            CS_l2I_mds['Attributes']['Auxiliary'][key] = {}
            for att_name,att_val in fileID['Auxiliary'][key].attrs.items():
                CS_l2I_mds['Attributes']['Auxiliary'][key][att_name] = att_val
        #-- CryoSat-2 External Corrections Group
        for key in fileID['Geometry'].keys():
            CS_l2I_mds['Attributes']['Geometry'][key] = {}
            for att_name,att_val in fileID['Geometry'][key].attrs.items():
                CS_l2I_mds['Attributes']['Geometry'][key][att_name] = att_val
        #-- CryoSat-2 Internal Corrections Group
        for key in fileID['Instrumental'].keys():
            CS_l2I_mds['Attributes']['Instrumental'][key] = {}
            for att_name,att_val in fileID['Instrumental'][key].attrs.items():
                CS_l2I_mds['Attributes']['Instrumental'][key][att_name] = att_val
        #-- Global attribute description
        CS_l2I_mds['Attributes']['title'] = fileID.attrs['description']

    #-- Closing the HDF5 file
    fileID.close()
    return CS_l2I_mds

#-- PURPOSE: get the number of records in an HDF5 file
def HDF5_cryosat_L2I_shape(FILENAME):
    #-- Open the HDF5 file for reading
    with h5py.File(os.path.expanduser(FILENAME), 'r') as fid:
        n_records,n_blocks = fid['Location']['Day'].shape
    return n_records
