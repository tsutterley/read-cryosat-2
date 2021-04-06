#!/usr/bin/env python
u"""
esa_cryosat_ftp.py
Written by Tyler Sutterley (04/2021)

This program syncs Cryosat Elevation products
From the ESA Cryosat ftp dissemination server:
https://earth.esa.int/web/guest/-/how-to-access-cryosat-data-6842
https://earth.esa.int/web/guest/-/products-overview-6975

INPUTS:
    CryoSat-2 product to sync with ESA servers
        SIR_LRM_L2: CryoSat-2 Low-Resolution Mode
        SIR_SAR_L2: CryoSat-2 SAR Mode
        SIR_SIN_L2: CryoSat-2 SARin Mode

CALLING SEQUENCE:
    python esa_cryosat_ftp.py --baseline=C --user=<username> SIR_SIN_L2
    where <username> is your ESA data dissemination server username

COMMAND LINE OPTIONS:
    --help: list the command line options
    -Y X, --year X: years to sync separated by commas
    -B X, --baseline X: CryoSat-2 baseline to sync
    --user: username for CryoSat-2 FTP servers
    -N X, --netrc X: path to .netrc file for authentication
    --directory: working data directory (default: current working directory)
    --bbox X: Bounding box (lonmin,latmin,lonmax,latmax)
    --polygon X: Georeferenced file containing a set of polygons
    -M X, --mode X: Local permissions mode of the directories and files synced
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
    Updated 04/2021: set a default netrc file and check access
        default credentials from environmental variables
    Updated 10/2020: using argparse to set parameters
    Updated 07/2020: added netrc option for alternative authentication
    Updated 03/2020: add spatial subsetting to reduce files to sync
        increase ftplib timeout to prevent connection drops
    Updated 02/2020: convert from hard to soft tabulation
    Updated 08/2019: include baseline in regular expression patterns
    Updated 06/2018: using python3 compatible octal and input
    Updated 05/2018 for public release.
    Updated 05/2017: exception if ESA Cryosat-2 credentials weren't entered
        using os.makedirs to recursively create directories
        using getpass to enter server password securely (remove --password)
    Updated 04/2017: minor changes to check_connection function to use ftplib
        updated regular expression for months to include year of interest
    Updated 02/2017: switching username and password to login command
    Written 11/2016
"""
from __future__ import print_function

import sys
import re
import os
import io
import netrc
import getpass
import argparse
import builtins
import lxml.etree
import calendar, time
import ftplib, posixpath
import cryosat_toolkit.polygon
import shapely.geometry

#-- PURPOSE: check internet connection
def check_connection(USER, PASSWORD):
    #-- attempt to connect to ftp host for Cryosat-2 servers
    try:
        f = ftplib.FTP('science-pds.cryosat.esa.int')
        f.login(USER, PASSWORD)
        f.voidcmd("NOOP")
    except IOError:
        raise RuntimeError('Check internet connection')
    except ftplib.error_perm:
        raise RuntimeError('Check login credentials')
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
    regex_class = r'OFFL|NRT_|RPRO|TEST|LTA_'
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
    regex_baseline = r'({0})'.format(BASELINE) if BASELINE else r'(.*?)'
    #-- CRYOSAT LEVEL-2 PRODUCTS NAMING RULES
    #-- Mission Identifier
    #-- File Class
    #-- File Product
    #-- Validity Start Date and Time
    #-- Validity Stop Date and Time
    #-- Baseline Identifier
    #-- Version Number
    regex_pattern = r'({0})_({1})_({2})_({3})_({4})_{5}(\d+).{6}$'.format('CS',
        regex_class,regex_products[PRODUCT],START,STOP,regex_baseline,SUFFIX)
    return re.compile(regex_pattern, re.VERBOSE)

#-- PURPOSE: sync local Cryosat-2 files with ESA server
def esa_cryosat_ftp(PRODUCT, YEARS, BASELINE=None, DIRECTORY=None,
    USER='', PASSWORD='', BBOX=None, POLYGON=None, LOG=False, LIST=False,
    MODE=None, CLOBBER=False):
    #-- connect and login to ESA ftp server
    f = ftplib.FTP('science-pds.cryosat.esa.int', timeout=3600)
    f.login(USER, PASSWORD)
    #-- compile xml parser for lxml
    XMLparser = lxml.etree.XMLParser()

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

    #-- compile regular expression operator for years to sync
    regex_years = '|'.join('{0:d}'.format(y) for y in YEARS)
    R1 = re.compile(r'({0})'.format(regex_years), re.VERBOSE)
    #-- initial regular expression pattern for months of the year
    regex_months = '(' + '|'.join('{0:02d}'.format(m) for m in range(1,13)) + ')'

    #-- compile the regular expression operator to find CryoSat-2 files
    #-- spatially subset data using bounding box or polygon file
    if BBOX:
        #-- if using a bounding box to spatially subset data
        #-- only find header files to extract latitude and longitude coordinates
        R3 = compile_regex_pattern(PRODUCT, BASELINE, SUFFIX=r'(HDR)')
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
        R3 = compile_regex_pattern(PRODUCT, BASELINE, SUFFIX=r'(HDR)')
        #-- read shapefile, kml/kmz file or GeoJSON file
        fileBasename,fileExtension = os.path.splitext(POLYGON)
        #-- extract file name and subsetter indices lists
        match_object = re.match(r'(.*?)(\[(.*?)\])?$',POLYGON)
        f = os.path.expanduser(match_object.group(1))
        #-- read specific variables of interest
        v = match_object.group(3).split(',') if match_object.group(2) else None
        #-- get MultiPolygon object from input spatial file
        if fileExtension in ('.shp','.zip'):
            #-- if reading a shapefile or a zipped directory with a shapefile
            ZIP = (fileExtension == '.zip')
            m=cryosat_toolkit.polygon().from_shapefile(f,variables=v,zip=ZIP)
        elif fileExtension in ('.kml','.kmz'):
            #-- if reading a keyhole markup language (can be compressed kmz)
            KMZ = (fileExtension == '.kmz')
            m=cryosat_toolkit.polygon().from_kml(f,variables=v,kmz=KMZ)
        elif fileExtension in ('.json','.geojson'):
            #-- if reading a GeoJSON file
            m=cryosat_toolkit.polygon().from_geojson(f,variables=v)
        else:
            raise IOError('Unlisted polygon type ({0})'.format(fileExtension))
        #-- calculate the convex hull of the MultiPolygon object for subsetting
        poly_obj = m.convex_hull
        #-- Valid Polygon cannot have overlapping exterior or interior rings
        if (not poly_obj.is_valid):
            poly_obj = poly_obj.buffer(0)
    else:
        R3 = compile_regex_pattern(PRODUCT, BASELINE)

    #-- find remote yearly directories for PRODUCT within YEARS
    YRS = [R1.findall(Y).pop() for Y in f.nlst(PRODUCT) if R1.search(Y)]
    for Y in YRS:
        #-- compile regular expression operator for months in year to sync
        R2 = re.compile(posixpath.join(Y,regex_months), re.VERBOSE)
        #-- find remote monthly directories for PRODUCT within year
        MNS = [R2.findall(M).pop() for M in f.nlst(posixpath.join(PRODUCT,Y))
            if R2.search(M)]
        for M in MNS:
            #-- remote and local directory for data product of year and month
            remote_dir = posixpath.join(PRODUCT,Y,M)
            local_dir = os.path.join(DIRECTORY,PRODUCT,Y,M)
            #-- check if local directory exists and recursively create if not
            os.makedirs(local_dir,MODE) if not os.path.exists(local_dir) else None
            #-- get filenames from remote directory
            valid_lines = [fi for fi in f.nlst(remote_dir) if R3.search(fi)]
            #-- if spatially subsetting
            if BBOX or POLYGON:
                #-- for each header file
                for line in sorted(valid_lines):
                    #-- extract filename from regex object
                    f1 = R3.search(line).group(0)
                    remote_file = posixpath.join(remote_dir,f1)
                    local_file = os.path.join(local_dir,f1)
                    #-- extract information from filename
                    MI,CLASS,PRD,START,STOP,BSLN,VERS,SFX = R3.findall(f1).pop()
                    #-- read XML header file and check if intersecting
                    if parse_xml_file(f, remote_file, poly_obj, XMLparser):
                        #-- compile regular expression operator for times
                        R4 = compile_regex_pattern(PRODUCT, BASELINE,
                            START=START, STOP=STOP)
                        #-- will sync both the data and header files
                        subset=[f2 for f2 in f.nlst(remote_dir) if R4.search(f2)]
                        for subset_line in subset:
                            #-- extract filename from regex object
                            f2 = R4.search(subset_line).group(0)
                            remote_file = posixpath.join(remote_dir,f2)
                            local_file = os.path.join(local_dir,f2)
                            ftp_mirror_file(fid1,f,remote_file,local_file,
                                LIST,CLOBBER,MODE)
            else:
                #-- for each data and header file
                for line in sorted(valid_lines):
                    #-- extract filename from regex object
                    fi = R3.search(line).group(0)
                    remote_file = posixpath.join(remote_dir,fi)
                    local_file = os.path.join(local_dir,fi)
                    ftp_mirror_file(fid1,f,remote_file,local_file,
                        LIST,CLOBBER,MODE)

    #-- close the ftp connection
    f.quit()
    #-- close log file and set permissions level to MODE
    if LOG:
        fid1.close()
        os.chmod(os.path.join(DIRECTORY,LOGFILE), MODE)

#-- PURPOSE: pull and parse xml header file to check if intersecting subsetter
def parse_xml_file(ftp, remote_file, poly_obj, XMLparser):
    #-- copy remote file contents to BytesIO object
    fileobj = io.BytesIO()
    ftp.retrbinary('RETR {0}'.format(remote_file),fileobj.write)
    #-- rewind retrieved binary to start of file
    fileobj.seek(0)
    tree = lxml.etree.parse(fileobj,XMLparser)
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
def ftp_mirror_file(fid, ftp, remote_file, local_file, LIST, CLOBBER, MODE):
    #-- if file exists in file system: check if remote file is newer
    TEST = False
    OVERWRITE = ' (clobber)'
    #-- get last modified date of remote file and convert into unix time
    remote_mtime = get_mtime(ftp.sendcmd('MDTM {0}'.format(remote_file)))
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
        args=(posixpath.join('ftp://',ftp.host,remote_file),local_file,OVERWRITE)
        print('{0} -->\n\t{1}{2}\n'.format(*args), file=fid)
        #-- if executing copy command (not only printing the files)
        if not LIST:
            #-- copy remote file contents to local file
            with open(local_file, 'wb') as f:
                ftp.retrbinary('RETR {0}'.format(remote_file), f.write)
            #-- keep remote modification time of file and local access time
            os.utime(local_file, (os.stat(local_file).st_atime, remote_mtime))
            os.chmod(local_file, MODE)

#-- PURPOSE: returns the Unix timestamp value for a modification time
def get_mtime(collastmod, FORMAT="%Y%m%d%H%M%S"):
    lastmodtime = time.strptime(collastmod[4:], FORMAT)
    return calendar.timegm(lastmodtime)

#-- Main program that calls esa_cryosat_ftp()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Syncs Cryosat Elevation products
            from the ESA Cryosat ftp dissemination server
            """
    )
    #-- command line parameters
    parser.add_argument('product',
        metavar='PRODUCT', type=str, nargs='+',
        help='CryoSat-2 Product')
    #-- ESA ftp credentials
    parser.add_argument('--user','-U',
        type=str, default='',
        help='Username for CryoSat-2 FTP Login')
    parser.add_argument('--netrc','-N',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.path.join(os.path.expanduser('~'),'.netrc'),
        help='Path to .netrc file for authentication')
    #-- working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- CryoSat Baseline
    parser.add_argument('--baseline','-B',
        metavar='BASELINE', type=str,
        default='D', choices=['C','D'],
        help='CryoSat Baseline to sync')
    #-- temporal subsetting
    parser.add_argument('--year','-Y',
        type=int, nargs='+',
        help='Years to sync')
    #-- spatial subsetting
    parser.add_argument('--bbox','-b',
        type=float, nargs=4, metavar=('lon_min','lat_min','lon_max','lat_max'),
        help='Bounding box')
    parser.add_argument('--polygon','-p',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='Georeferenced file containing a set of polygons')
    #-- Output log file in form
    #-- ESA_CS_SIR_SIN_L2_sync_2002-04-01.log
    parser.add_argument('--log','-l',
        default=False, action='store_true',
        help='Output log file')
    #-- sync options
    parser.add_argument('--list','-L',
        default=False, action='store_true',
        help='Only print files that could be transferred')
    parser.add_argument('--clobber','-C',
        default=False, action='store_true',
        help='Overwrite existing data in transfer')
    #-- permissions mode of the directories and files synced (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files synced')
    args = parser.parse_args()

    #-- ESA CryoSat-2 FTP Server name
    HOST = 'science-pds.cryosat.esa.int'
    #-- get ESA CryoSat-2 FTP Server credentials
    if not args.user and not os.access(args.netrc,os.F_OK):
        #-- check that ESA CryoSat-2 FTP Server credentials were entered
        args.user=builtins.input('Username for {0}: '.format(HOST))
        #-- enter password securely from command-line
        PASSWORD=getpass.getpass('Password for {0}@{1}: '.format(args.user,HOST))
    elif not args.user and os.access(args.netrc,os.F_OK):
        args.user,_,PASSWORD = netrc.netrc(args.netrc).authenticators(HOST)
    else:
        #-- enter password securely from command-line
        PASSWORD=getpass.getpass('Password for {0}@{1}: '.format(args.user,HOST))

    #-- check internet connection before attempting to run program
    if check_connection(args.user,PASSWORD):
        for PRODUCT in args.product:
            esa_cryosat_ftp(PRODUCT, args.years, USER=args.user,
                PASSWORD=PASSWORD, BASELINE=args.baseline,
                DIRECTORY=args.directory, BBOX=args.bbox,
                POLYGON=args.polygon, LOG=args.log, LIST=args.list,
                CLOBBER=args.clobber, MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
