#!/usr/bin/env python
u"""
esa_cryosat_sync.py
Written by Tyler Sutterley (07/2020)

This program syncs Cryosat Elevation products
From the ESA CryoSat-2 Science Server:
https://earth.esa.int/web/guest/-/how-to-access-cryosat-data-6842
https://earth.esa.int/web/guest/-/products-overview-6975

INPUTS:
    CryoSat-2 product to sync with ESA servers
        SIR_LRM_L2: CryoSat-2 Low-Resolution Mode
        SIR_SAR_L2: CryoSat-2 SAR Mode
        SIR_SIN_L2: CryoSat-2 SARin Mode

CALLING SEQUENCE:
    python esa_cryosat_sync.py --baseline=D SIR_SIN_L2

COMMAND LINE OPTIONS:
    --help: list the command line options
    -Y X, --year=X: years to sync separated by commas
    -B X, --baseline=X: CryoSat-2 baseline to sync
    --directory: working data directory (default: current working directory)
    --bbox=X: Bounding box (lonmin,latmin,lonmax,latmax)
    --polygon=X: Georeferenced file containing a set of polygons
    -M X, --mode=X: Local permissions mode of the directories and files synced
    --log: output log of files downloaded
    --list: print files to be transferred, but do not execute transfer
    --clobber: Overwrite existing data in transfer

PYTHON DEPENDENCIES:
    lxml: Pythonic XML and HTML processing library using libxml2/libxslt
        https://lxml.de/
        https://github.com/lxml/lxml
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    gdal: Pythonic interface to the Geospatial Data Abstraction Library (GDAL)
        https://pypi.python.org/pypi/GDAL/
    fiona: Python wrapper for vector data access functions from the OGR library
        https://fiona.readthedocs.io/en/latest/manual.html
    geopandas: Python tools for geographic data
        http://geopandas.readthedocs.io/
    shapely: PostGIS-ish operations outside a database context for Python
        http://toblerity.org/shapely/index.html
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/
    future: Compatibility layer between Python 2 and Python 3
        (http://python-future.org/)

UPDATE HISTORY:
    Updated 07/2020: create compiled lists of filenames and last modified times
    Updated 06/2020: use json to decode API output
    Updated 03/2020: add spatial subsetting to reduce files to sync
    Updated 02/2020: convert from hard to soft tabulation
    Updated 09/2019: increase urlopen timeout. place regex compilations together
    Updated 08/2019: include baseline in regular expression patterns
        updated for https Cryosat-2 Science Server
        ftp program renamed esa_cryosat_ftp.py
    Updated 06/2018: using python3 compatible octal and input
    Updated 05/2018 for public release.
    Updated 05/2017: exception if ESA Cryosat-2 credentials weren't entered
        using os.makedirs to recursively create directories
    Updated 04/2017: minor changes to check_connection function to use ftplib
        updated regular expression for months to include year of interest
    Updated 02/2017: switching username and password to login command
    Written 11/2016
"""
from __future__ import print_function

import sys
import re
import os
import ssl
import json
import getopt
import shutil
import base64
import builtins
import posixpath
import lxml.etree
import calendar, time
import shapely.geometry
from cryosat_toolkit.read_shapefile import read_shapefile
from cryosat_toolkit.read_kml_file import read_kml_file
from cryosat_toolkit.read_geojson_file import read_geojson_file
if sys.version_info[0] == 2:
    from cookielib import CookieJar
    from urllib import urlencode
    import urllib2
else:
    from http.cookiejar import CookieJar
    from urllib.parse import urlencode
    import urllib.request as urllib2

#-- PURPOSE: check internet connection
def check_connection():
    #-- attempt to connect to https ESA CryoSat-2 Science Server
    try:
        HOST = 'https://science-pds.cryosat.esa.int'
        urllib2.urlopen(HOST,timeout=20,context=ssl.SSLContext())
    except urllib2.URLError:
        raise RuntimeError('Check internet connection')
    else:
        return True

#-- PURPOSE: compile regular expression operator to find CryoSat-2 files
def compile_regex_pattern(PRODUCT, BASELINE, START='\d+T?\d+', STOP='\d+T?\d+',
    SUFFIX='(DBL|HDR|nc)'):
    #-- CryoSat file class
    #-- OFFL (Off Line Processing/Systematic)
    #-- NRT_ (Near Real Time)
    #-- RPRO (ReProcessing)
    #-- TEST (Testing)
    #-- LTA_ (Long Term Archive)
    regex_class = 'OFFL|NRT_|RPRO|TEST|LTA_'
    #-- CryoSat mission products
    #-- SIR_LRM_1B L1B Product from Low Resolution Mode Processing
    #-- SIR_FDM_1B L1B Product from Fast Delivery Marine Mode Processing
    #-- SIR_SIN_1B L1B Product from SAR Interferometric Processing
    #-- SIR_SID_1B L1B Product from SIN Degraded Processing
    #-- SIR_SAR_1B L1B Product from SAR Processing
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
    regex_products = {}
    regex_products['SIR_LRM_L1'] = 'SIR_LRM_1B'
    regex_products['SIR_FDM_L1'] = 'SIR_FDM_1B'
    regex_products['SIR_SIN_L1'] = 'SIR_SIN_1B'
    regex_products['SIR_SID_L1'] = 'SIR_SID_1B'
    regex_products['SIR_SAR_L1'] = 'SIR_SAR_1B'
    regex_products['SIR_LRM_L2'] = 'SIR_LRM_2_'
    regex_products['SIR_FDM_L2'] = 'SIR_FDM_2_'
    regex_products['SIR_SIN_L2'] = 'SIR_SIN_2_'
    regex_products['SIR_SID_L2'] = 'SIR_SID_2_'
    regex_products['SIR_SAR_L2'] = 'SIR_SAR_2_'
    regex_products['SIR_GDR_L2'] = 'SIR_GDR_2_'
    regex_products['SIR_LRM_L2I'] = 'SIR_LRMI2_'
    regex_products['SIR_SIN_L2I'] = 'SIR_SINI2_'
    regex_products['SIR_SID_L2I'] = 'SIR_SIDI2_'
    regex_products['SIR_SAR_L2I'] = 'SIR_SARI2_'
    #-- Cryosat baseline Identifier
    regex_baseline = '({0})'.format(BASELINE) if BASELINE else '(.*?)'
    #-- CRYOSAT LEVEL-2 PRODUCTS NAMING RULES
    #-- Mission Identifier
    #-- File Class
    #-- File Product
    #-- Validity Start Date and Time
    #-- Validity Stop Date and Time
    #-- Baseline Identifier
    #-- Version Number
    regex_pattern = '({0})_({1})_({2})_({3})_({4})_{5}(\d+).{6}$'.format('CS',
        regex_class,regex_products[PRODUCT],START,STOP,regex_baseline,SUFFIX)
    return re.compile(regex_pattern, re.VERBOSE)

#-- PURPOSE: sync local Cryosat-2 files with ESA server
def esa_cryosat_sync(PRODUCT, YEARS, BASELINE=None, DIRECTORY=None, BBOX=None,
    POLYGON=None, LOG=False, LIST=False, MODE=None, CLOBBER=False):

    #-- create log file with list of synchronized files (or print to terminal)
    if LOG:
        #-- check if log directory exists and recursively create if not
        os.makedirs(DIRECTORY,MODE) if not os.path.exists(DIRECTORY) else None
        #-- format: ESA_CS_SIR_SIN_L2_sync_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = 'ESA_CS_{0}_sync_{1}.log'.format(PRODUCT,today)
        fid1 = open(os.path.join(DIRECTORY,LOGFILE),'w')
        print('ESA CryoSat-2 Sync Log ({0})'.format(today), file=fid1)
        print('PRODUCT={0}'.format(PRODUCT), file=fid1)
    else:
        #-- standard output (terminal output)
        fid1 = sys.stdout

    #-- CryoSat-2 Science Server url
    #-- using the JSON api protocols to retrieve files
    #-- static site is no longer available
    HOST = posixpath.join('https://science-pds.cryosat.esa.int')
    #-- compile xml parsers for lxml
    XMLparser = lxml.etree.XMLParser()
    #-- Create cookie jar for storing cookies
    cookie_jar = CookieJar()
    #-- create "opener" (OpenerDirector instance)
    opener = urllib2.build_opener(
        urllib2.HTTPSHandler(context=ssl.SSLContext()),
        urllib2.HTTPCookieProcessor(cookie_jar))
    #-- Now all calls to urllib2.urlopen use our opener.
    urllib2.install_opener(opener)
    #-- All calls to urllib2.urlopen will now use handler
    #-- Make sure not to include the protocol in with the URL, or
    #-- HTTPPasswordMgrWithDefaultRealm will be confused.

    #-- compile regular expression operator for years to sync
    regex_years = '|'.join('{0:d}'.format(y) for y in YEARS)
    R1 = re.compile('({0})'.format(regex_years), re.VERBOSE)
    #-- regular expression pattern for months of the year
    regex_months = '|'.join('{0:02d}'.format(m) for m in range(1,13))
    R2 = re.compile('({0})'.format(regex_months), re.VERBOSE)

    #-- compile the regular expression operator to find CryoSat-2 files
    #-- spatially subset data using bounding box or polygon file
    if BBOX:
        #-- if using a bounding box to spatially subset data
        #-- only find header files to extract latitude and longitude coordinates
        R3 = compile_regex_pattern(PRODUCT, BASELINE, SUFFIX='(HDR)')
        #-- min_lon,min_lat,max_lon,max_lat
        lon = [BBOX[0],BBOX[2],BBOX[2],BBOX[0],BBOX[0]]
        lat = [BBOX[1],BBOX[1],BBOX[3],BBOX[3],BBOX[1]]
        #-- create shapely polygon
        poly_obj = shapely.geometry.Polygon(list(zip(lon, lat)))
        #-- Valid Polygon cannot have overlapping exterior or interior rings
        if (not poly_obj.is_valid):
            poly_obj = poly_obj.buffer(0)
    elif POLYGON:
        #-- if using a polygon file to spatially subset data
        #-- only find header files to extract latitude and longitude coordinates
        R3 = compile_regex_pattern(PRODUCT, BASELINE, SUFFIX='(HDR)')
        #-- read shapefile, kml/kmz file or GeoJSON file
        fileBasename,fileExtension = os.path.splitext(POLYGON)
        #-- extract file name and subsetter indices lists
        match_object = re.match('(.*?)(\[(.*?)\])?$',POLYGON)
        FILE = os.path.expanduser(match_object.group(1))
        #-- read specific variables of interest
        v = match_object.group(3).split(',') if match_object.group(2) else None
        #-- get MultiPolygon object from input spatial file
        if fileExtension in ('.shp','.zip'):
            #-- if reading a shapefile or a zipped directory with a shapefile
            ZIP = (fileExtension == '.zip')
            m = read_shapefile(os.path.expanduser(FILE), VARIABLES=v, ZIP=ZIP)
        elif fileExtension in ('.kml','.kmz'):
            #-- if reading a keyhole markup language (can be compressed)
            KMZ = (fileExtension == '.kmz')
            m = read_kml_file(os.path.expanduser(FILE), VARIABLES=v, KMZ=KMZ)
        elif fileExtension in ('.json','.geojson'):
            #-- if reading a GeoJSON file
            m = read_geojson_file(os.path.expanduser(FILE), VARIABLES=v)
        else:
            raise IOError('Unlisted polygon type ({0})'.format(fileExtension))
        #-- calculate the convex hull of the MultiPolygon object for subsetting
        poly_obj = m.convex_hull
        #-- Valid Polygon cannot have overlapping exterior or interior rings
        if (not poly_obj.is_valid):
            poly_obj = poly_obj.buffer(0)
    else:
        R3 = compile_regex_pattern(PRODUCT, BASELINE)

    #-- open connection with Cryosat-2 science server at remote directory
    #-- [sic Cry0Sat2_data]
    parameters = {'file':posixpath.join('Cry0Sat2_data',PRODUCT)}
    url = posixpath.join(HOST,'?do=list&{0}'.format(urlencode(parameters)))
    request = urllib2.Request(url=url)
    response = urllib2.urlopen(request,timeout=60)
    table = json.loads(response.read().decode())
    #-- find remote yearly directories for PRODUCT within YEARS
    YRS = [t['name'] for t in table['results'] if R1.match(t['name'])]
    for Y in YRS:
        #-- open connection with Cryosat-2 science server at remote directory
        #-- [sic Cry0Sat2_data]
        parameters = {'file':posixpath.join('Cry0Sat2_data',PRODUCT,Y)}
        url = posixpath.join(HOST,'?do=list&{0}'.format(urlencode(parameters)))
        request = urllib2.Request(url=url)
        response = urllib2.urlopen(request,timeout=360)
        table = json.loads(response.read().decode())
        #-- find remote monthly directories for PRODUCT within year
        MNS = [t['name'] for t in table['results'] if R2.match(t['name'])]
        for M in MNS:
            #-- local directory for data product of year and month
            local_dir = os.path.join(DIRECTORY,PRODUCT,Y,M)
            #-- check if local directory exists and recursively create if not
            os.makedirs(local_dir,MODE) if not os.path.exists(local_dir) else None

            #-- create combined list of filenames and last modified times
            #-- find remote files for PRODUCT within year and month
            colnames = []
            collastmod = []
            #-- position, maximum number of files to list, flag to check if done
            pos,maxfiles,prevmax = (0,500,500)
            #-- iterate to get a compiled list of files
            #-- will iterate until there are no more files to add to the lists
            while (maxfiles == prevmax):
                #-- set previous flag to maximum 
                prevmax = maxfiles
                #-- open connection with Cryosat-2 science server at remote directory
                #-- to list maxfiles number of files at position [sic Cry0Sat2_data]
                parameters = {'maxfiles':prevmax,'pos':pos,
                    'file':posixpath.join('Cry0Sat2_data',PRODUCT,Y,M)}
                url=posixpath.join(HOST,'?do=list&{0}'.format(urlencode(parameters)))
                request = urllib2.Request(url=url)
                response = urllib2.urlopen(request,timeout=360)
                table = json.loads(response.read().decode())
                #-- extend lists with new files
                colnames.extend([t['name'] for t in table['results']])
                collastmod.extend([t['mtime'] for t in table['results']])
                #-- update maximum number of files
                maxfiles = len(table['results'])
                #-- update position
                pos += maxfiles

            #-- if spatially subsetting
            if BBOX or POLYGON:
                #-- find names of valid header files
                header_files = [f for i,f in enumerate(colnames) if R3.match(f)]
                for hf in sorted(header_files):
                    #-- remote and local versions of the file
                    #-- [sic Cry0Sat2_data]
                    parameters = {'file':posixpath.join('Cry0Sat2_data',
                        PRODUCT,Y,M,hf)}
                    remote_file = posixpath.join(HOST,
                        '?do=download&{0}'.format(urlencode(parameters)))
                    #-- extract information from filename
                    MI,CLASS,PRD,START,STOP,BSLN,VERS,SFX = R3.findall(hf).pop()
                    #-- read XML header file and check if intersecting
                    if parse_xml_file(remote_file, poly_obj, XMLparser):
                        #-- compile regular expression operator for times
                        R4 = compile_regex_pattern(PRODUCT, BASELINE,
                            START=START, STOP=STOP)
                        #-- will sync both the data and header files
                        subset=[i for i,f in enumerate(colnames) if R4.match(f)]
                        for i in subset:
                            #-- remote and local versions of the file
                            #-- [sic Cry0Sat2_data]
                            parameters = {'file':posixpath.join('Cry0Sat2_data',
                                PRODUCT,Y,M,colnames[i])}
                            remote_file = posixpath.join(HOST,
                                '?do=download&{0}'.format(urlencode(parameters)))
                            local_file = os.path.join(local_dir,colnames[i])
                            #-- get last modified date in unix time
                            remote_mtime = collastmod[i]
                            http_pull_file(fid1, remote_file, remote_mtime,
                                local_file, LIST, CLOBBER, MODE)
            else:
                #-- find lines of valid files
                valid_lines = [i for i,f in enumerate(colnames) if R3.match(f)]
                #-- for each data and header file
                for i in valid_lines:
                    #-- remote and local versions of the file
                    #-- [sic Cry0Sat2_data]
                    parameters = {'file':posixpath.join('Cry0Sat2_data',
                        PRODUCT,Y,M,colnames[i])}
                    remote_file = posixpath.join(HOST,
                        '?do=download&{0}'.format(urlencode(parameters)))
                    local_file = os.path.join(local_dir,colnames[i])
                    #-- get last modified date in unix time
                    remote_mtime = collastmod[i]
                    #-- check that file is not in file system unless overwriting
                    http_pull_file(fid1, remote_file, remote_mtime, local_file,
                        LIST, CLOBBER, MODE)

    #-- close log file and set permissions level to MODE
    if LOG:
        fid1.close()
        os.chmod(os.path.join(DIRECTORY,LOGFILE), MODE)

#-- PURPOSE: pull and parse xml header file to check if intersecting subsetter
def parse_xml_file(remote_file, poly_obj, XMLparser):
    request = urllib2.Request(remote_file)
    response = urllib2.urlopen(request,timeout=20)
    tree = lxml.etree.parse(response,XMLparser)
    #-- extract starting latitude/longitude and ending latitude/longitude
    Start_Lat, = tree.xpath('//Product_Location/Start_Lat/text()')
    Start_Long, = tree.xpath('//Product_Location/Start_Long/text()')
    Stop_Lat, = tree.xpath('//Product_Location/Stop_Lat/text()')
    Stop_Long, = tree.xpath('//Product_Location/Stop_Long/text()')
    #-- convert latitude and longitude to degrees from microdegrees
    lon = [float(Start_Long)/1e6,float(Stop_Long)/1e6]
    lat = [float(Start_Lat)/1e6,float(Stop_Lat)/1e6]
    #-- create shapely line string object and buffer by 1 degree
    line_obj = shapely.geometry.LineString(list(zip(lon, lat))).buffer(1)
    #-- check if intersecting and return flag
    return poly_obj.intersects(line_obj)

#-- PURPOSE: pull file from a remote host checking if file exists locally
#-- and if the remote file is newer than the local file
def http_pull_file(fid,remote_file,remote_mtime,local_file,LIST,CLOBBER,MODE):
    #-- if file exists in file system: check if remote file is newer
    TEST = False
    OVERWRITE = ' (clobber)'
    #-- check if local version of file exists
    if os.access(local_file, os.F_OK):
        #-- check last modification time of local file
        local_mtime = os.stat(local_file).st_mtime
        #-- if remote file is newer: overwrite the local file
        if (remote_mtime > local_mtime):
            TEST = True
            OVERWRITE = ' (overwrite)'
    else:
        TEST = True
        OVERWRITE = ' (new)'
    #-- if file does not exist locally, is to be overwritten, or CLOBBER is set
    if TEST or CLOBBER:
        #-- Printing files transferred
        print('{0} --> '.format(remote_file), file=fid)
        print('\t{0}{1}\n'.format(local_file,OVERWRITE), file=fid)
        #-- if executing copy command (not only printing the files)
        if not LIST:
            #-- Create and submit request. There are a wide range of exceptions
            #-- that can be thrown here, including HTTPError and URLError.
            req = urllib2.Request(remote_file)
            resp = urllib2.urlopen(req,timeout=360)
            #-- chunked transfer encoding size
            CHUNK = 16 * 1024
            #-- copy contents to local file using chunked transfer encoding
            #-- transfer should work properly with ascii and binary data formats
            with open(local_file, 'wb') as f:
                shutil.copyfileobj(resp, f, CHUNK)
            #-- keep remote modification time of file and local access time
            os.utime(local_file, (os.stat(local_file).st_atime, remote_mtime))
            os.chmod(local_file, MODE)

#-- PURPOSE: help module to describe the optional input parameters
def usage():
    print('\nHelp: {}'.format(os.path.basename(sys.argv[0])))
    print(' -Y X, --year=X\t\tYears to sync separated by commas')
    print(' -B X, --baseline=X\tCryoSat Baseline to run')
    print(' --directory=X\t\tWorking data directory')
    print(' --bbox=X\t\tBounding box (lonmin,latmin,lonmax,latmax)')
    print(' --polygon=X\t\tGeoreferenced file containing a set of polygons')
    print(' -M X, --mode=X\t\tPermission mode of directories and files synced')
    print(' -L, --list\t\tOnly print files that are to be transferred')
    print(' -C, --clobber\t\tOverwrite existing data in transfer')
    print(' -l, --log\t\tOutput log file')
    today = time.strftime('%Y-%m-%d',time.localtime())
    LOGFILE = 'ESA_CS_{0}_sync_{1}.log'.format('SIR_SIN_L2',today)
    print('     Log file format: {}\n'.format(LOGFILE))

#-- Main program that calls esa_cryosat_sync()
def main():
    #-- Read the system arguments listed after the program
    long_options = ['help','year=','baseline=','directory=','bbox=','polygon=',
        'list','log','mode=','clobber']
    optlist,arglist = getopt.getopt(sys.argv[1:],'hY:B:LCM:l',long_options)

    #-- command line parameters
    YEAR = [2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020]
    BASELINE = 'D'
    DIRECTORY = os.getcwd()
    BBOX = None
    POLYGON = None
    LIST = False
    LOG = False
    #-- permissions mode of the local directories and files (number in octal)
    MODE = 0o775
    CLOBBER = False
    for opt, arg in optlist:
        if opt in ('-h','--help'):
            usage()
            sys.exit()
        elif opt in ("-Y","--year"):
            YEAR = [int(Y) for Y in arg.split(',')]
        elif opt in ("-B","--baseline"):
            BASELINE = arg.upper()
        elif opt in ("--directory"):
            DIRECTORY = os.path.expanduser(arg)
        elif opt in ("--bbox",):
            BBOX = [float(i) for i in arg.split(',')]
        elif opt in ("--polygon",):
            POLYGON = os.path.expanduser(arg)
        elif opt in ("-L","--list"):
            LIST = True
        elif opt in ("-l","--log"):
            LOG = True
        elif opt in ("-M","--mode"):
            MODE = int(arg, 8)
        elif opt in ("-C","--clobber"):
            CLOBBER = True

    if not arglist:
        #-- Input CryoSat-2 Product (sys.argv[0] is the python code)
        raise Exception('No CryoSat-2 Product Specified')

    #-- check internet connection before attempting to run program
    if check_connection():
        for PRODUCT in arglist:
            esa_cryosat_sync(PRODUCT, YEAR, BASELINE=BASELINE,
                DIRECTORY=DIRECTORY, BBOX=BBOX, POLYGON=POLYGON,
                LOG=LOG, LIST=LIST, MODE=MODE, CLOBBER=CLOBBER)

#-- run main program
if __name__ == '__main__':
    main()
