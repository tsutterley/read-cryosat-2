#!/usr/bin/env python
u"""
test_download_and_read.py (08/2020)
Tests the read program to verify that coefficients are being extracted
"""
import sys
import warnings
import pytest
import posixpath
import cryosat_toolkit.utilities
from cryosat_toolkit.read_cryosat_L1b import read_cryosat_L1b
from cryosat_toolkit.read_cryosat_L2 import read_cryosat_L2
if sys.version_info[0] == 2:
    from urllib import urlencode
else:
    from urllib.parse import urlencode

# filepaths for Baseline-E L1/L2 LRM/SAR/SIN files
filepath = []
# LRM L1
filepath.append(['SIR_LRM_L1','2022','04',
    'CS_OFFL_SIR_LRM_1B_20220416T000235_20220416T000655_E001.nc'])
# LRM L2
filepath.append(['SIR_LRM_L2','2022','04',
    'CS_OFFL_SIR_LRM_2__20220416T000235_20220416T000655_E001.nc'])
# SAR L1
filepath.append(['SIR_SAR_L1','2022','04',
    'CS_OFFL_SIR_SAR_1B_20220416T000940_20220416T001027_E001.nc'])
# SAR L2
filepath.append(['SIR_SAR_L2','2022','04',
    'CS_OFFL_SIR_SAR_2__20220416T000940_20220416T001027_E001.nc'])
# SARin L1
filepath.append(['SIR_SIN_L1','2022','04',
    'CS_OFFL_SIR_SIN_1B_20220416T000655_20220416T000940_E001.nc'])
# SARin L2
filepath.append(['SIR_SIN_L2','2022','04',
    'CS_OFFL_SIR_SINI2__20220416T000655_20220416T000940_E001.nc'])

# parameterize files to download from https
@pytest.mark.parametrize("filepath", filepath)
# PURPOSE: Download files from ESA https and check read program for file
@pytest.mark.skip(reason='HTTPS Science Server Currently Offline')
def test_https_download_and_read(filepath):
    # encode file path
    url=urlencode({'file':posixpath.join('Cry0Sat2_data',*filepath)})
    HOST=['https://science-pds.cryosat.esa.int','?do=download&{0}'.format(url)]
    # download file from https
    cryosat_toolkit.utilities.from_http(HOST,local=filepath[-1],verbose=True)
    # read file with level-1b or level-2 programs
    if filepath[0] in ('SIR_LRM_L1','SIR_SAR_L1','SIR_SIN_L1'):
        CS_L1b_mds = read_cryosat_L1b(filepath[-1], VERBOSE=True)
    elif filepath[0] in ('SIR_LRM_L2','SIR_SAR_L2','SIR_SIN_L2'):
        CS_L2_mds = read_cryosat_L2(filepath[-1], VERBOSE=True)

# filepaths for Baseline-C L1/L2 LRM/SAR/SIN files
filepath = []
# LRM L1
filepath.append(['Ice_Baseline_C','SIR_LRM_L1','2010','07',
    'CS_LTA__SIR_LRM_1B_20100716T001203_20100716T001716_C001.DBL'])
# LRM L2
filepath.append(['Ice_Baseline_C','SIR_LRM_L2','2010','07',
    'CS_LTA__SIR_LRM_2__20100716T001203_20100716T001716_C001.DBL'])
# SAR L1
filepath.append(['Ice_Baseline_C','SIR_SAR_L1','2010','07',
    'CS_LTA__SIR_SAR_1B_20100716T000457_20100716T000636_C001.DBL'])
# SAR L2
filepath.append(['Ice_Baseline_C','SIR_SAR_L2','2010','07',
    'CS_LTA__SIR_SAR_2__20100716T000457_20100716T000636_C001.DBL'])
# SARin L1
filepath.append(['Ice_Baseline_C','SIR_SIN_L1','2010','07',
    'CS_LTA__SIR_SIN_1B_20100716T000848_20100716T001203_C001.DBL'])
# SARin L2
filepath.append(['Ice_Baseline_C','SIR_SIN_L2','2010','07',
    'CS_LTA__SIR_SIN_2__20100716T000848_20100716T001203_C001.DBL'])
# parameterize files to download from ftp
@pytest.mark.parametrize("filepath", filepath)
# PURPOSE: Download files from ESA ftp and check read program for file
@pytest.mark.skip(reason='Deprecated Version')
def test_ftp_download_and_read(username,password,filepath):
    # download file from ftp
    cryosat_toolkit.utilities.from_ftp(['science-pds.cryosat.esa.int',*filepath],
        username=username,password=password,local=filepath[-1],verbose=True)
    # read file with level-1b or level-2 programs
    if filepath[0] in ('SIR_LRM_L1','SIR_SAR_L1','SIR_SIN_L1'):
        CS_L1b_mds = read_cryosat_L1b(filepath[-1], VERBOSE=True)
    elif filepath[0] in ('SIR_LRM_L2','SIR_SAR_L2','SIR_SIN_L2'):
        CS_L2_mds = read_cryosat_L2(filepath[-1], VERBOSE=True)
