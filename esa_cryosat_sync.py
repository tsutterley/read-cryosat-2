#!/usr/bin/env python
u"""
esa_cryosat_sync.py
Written by Tyler Sutterley (09/2019)

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
	-M X, --mode=X: Local permissions mode of the directories and files synced
	--log: output log of files downloaded
	--list: print files to be transferred, but do not execute transfer
	--clobber: Overwrite existing data in transfer

PYTHON DEPENDENCIES:
	lxml: Pythonic XML and HTML processing library using libxml2/libxslt
		http://lxml.de/
		https://github.com/lxml/lxml
	future: Compatibility layer between Python 2 and Python 3
		(http://python-future.org/)

UPDATE HISTORY:
	Updated 09/2019: increase urlopen timeout. place regex compilations together
	Updated 08/2019: include baseline in regular expression patterns
		updated for https Cryosat-2 Science Server
		ftp program renamed esa_cryosat_ftp.py
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
import ssl
import shutil
import getopt
import getpass
import builtins
import posixpath
import lxml.etree
import calendar, time
if sys.version_info[0] == 2:
	import urllib2
else:
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
def compile_regex_pattern(PRODUCT, BASELINE):
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
	regex_products = {}
	regex_products['SIR_LRM_L2'] = 'SIR_LRM_2'
	regex_products['SIR_FDM_L2'] = 'SIR_FDM_2'
	regex_products['SIR_SIN_L2'] = 'SIR_SIN_2'
	regex_products['SIR_SID_L2'] = 'SIR_SID_2'
	regex_products['SIR_SAR_L2'] = 'SIR_SAR_2'
	regex_products['SIR_GDR_L2'] = 'SIR_GDR_2'
	regex_products['SIR_LRM_L2I'] = 'SIR_LRMI2'
	regex_products['SIR_SIN_L2I'] = 'SIR_SINI2'
	regex_products['SIR_SID_L2I'] = 'SIR_SIDI2'
	regex_products['SIR_SAR_L2I'] = 'SIR_SARI2'
	#-- Cryosat baseline Identifier
	regex_baseline = '({0})'.format(BASELINE) if BASELINE else '(.*?)'
	#-- Cryosat suffix
	suffix = '(DBL|HDR|nc)'
	#-- CRYOSAT LEVEL-2 PRODUCTS NAMING RULES
	#-- Mission Identifier
	#-- File Class
	#-- File Product
	#-- Validity Start Date and Time
	#-- Validity Stop Date and Time
	#-- Baseline Identifier
	#-- Version Number
	args = ('CS',regex_class,regex_products[PRODUCT],regex_baseline,suffix)
	regex_pattern = '({0})_({1})_({2})__(\d+T?\d+)_(\d+T?\d+)_{3}(\d+).{4}$'
	return re.compile(regex_pattern.format(*args), re.VERBOSE)

#-- PURPOSE: sync local Cryosat-2 files with ESA server
def esa_cryosat_sync(PRODUCT, YEARS, DIRECTORY=None, USER='', PASSWORD='',
	BASELINE=None, LOG=False, LIST=False, MODE=None, CLOBBER=False):

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

	#-- CryoSat-2 Science Server url [sic Cry0Sat2_data]
	HOST = posixpath.join('https://science-pds.cryosat.esa.int','Cry0Sat2_data')
	#-- compile HTML parser for lxml
	parser = lxml.etree.HTMLParser()
	#-- ssl context
	context = ssl.SSLContext()
	#-- compile regular expression operator for years to sync
	regex_years = '|'.join('{0:d}'.format(y) for y in YEARS)
	R1 = re.compile('({0})'.format(regex_years), re.VERBOSE)
	#-- regular expression pattern for months of the year
	regex_months = '|'.join('{0:02d}'.format(m) for m in range(1,13))
	R2 = re.compile('({0})'.format(regex_months), re.VERBOSE)
	#-- compile the regular expression operator to find CryoSat-2 files
	R3 = compile_regex_pattern(PRODUCT, BASELINE)

	#-- open connection with Cryosat-2 science server at remote directory
	request = urllib2.Request(url=posixpath.join(HOST,PRODUCT))
	response = urllib2.urlopen(request,timeout=60,context=context)
	tree = lxml.etree.parse(response, parser)
	#-- find remote yearly directories for PRODUCT within YEARS
	colnames = tree.xpath('//tr/td//a/text()')
	YRS = [d for i,d in enumerate(colnames) if R1.match(d)]
	for Y in YRS:
		#-- open connection with Cryosat-2 science server at remote directory
		request = urllib2.Request(url=posixpath.join(HOST,PRODUCT,Y))
		response = urllib2.urlopen(request,timeout=360,context=context)
		tree = lxml.etree.parse(response, parser)
		#-- find remote monthly directories for PRODUCT within year
		colnames = tree.xpath('//tr/td//a/text()')
		MNS = [d for i,d in enumerate(colnames) if R2.match(d)]
		for M in MNS:
			#-- local directory for data product of year and month
			local_dir = os.path.join(DIRECTORY,PRODUCT,Y,M)
			#-- check if local directory exists and recursively create if not
			os.makedirs(local_dir,MODE) if not os.path.exists(local_dir) else None
			#-- open connection with Cryosat-2 science server at remote directory
			request = urllib2.Request(url=posixpath.join(HOST,PRODUCT,Y,M))
			response = urllib2.urlopen(request,timeout=360,context=context)
			tree = lxml.etree.parse(response, parser)
			#-- find remote yearly directories for PRODUCT within YEARS
			colnames = tree.xpath('//tr/td//a/text()')
			collastmod = tree.xpath('//tr/td[@align="right"][1]/text()')
			valid_lines = [i for i,f in enumerate(colnames) if R3.match(f)]
			for i in valid_lines:
				#-- remote and local versions of the file
				remote_file = posixpath.join(HOST,PRODUCT,Y,M,colnames[i])
				local_file = os.path.join(local_dir,colnames[i])
				#-- get last modified date and convert into unix time
				LMD = time.strptime(collastmod[i].rstrip(),'%d-%b-%Y %H:%M')
				remote_mtime = calendar.timegm(LMD)
				http_pull_file(fid1, remote_file, remote_mtime, local_file,
					LIST, CLOBBER, MODE)

	#-- close log file and set permissions level to MODE
	if LOG:
		fid1.close()
		os.chmod(os.path.join(DIRECTORY,LOGFILE), MODE)

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
			resp = urllib2.urlopen(req,timeout=360,context=ssl.SSLContext())
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
	print(' -M X, --mode=X\t\tPermission mode of directories and files synced')
	print(' -L, --list\t\tOnly print files that are to be transferred')
	print(' -C, --clobber\t\tOverwrite existing data in transfer')
	print(' -l, --log\t\tOutput log file')
	today = time.strftime('%Y-%m-%d',time.localtime())
	LOGFILE = 'ESA_CS_{0}_sync_{1}.log'.format('SIR_SIN_L2',today)
	print('    Log file format: {}\n'.format(LOGFILE))

#-- Main program that calls esa_cryosat_sync()
def main():
	#-- Read the system arguments listed after the program
	long_options = ['help','year=','baseline=','directory=','list','log',
		'mode=','clobber']
	optlist,arglist = getopt.getopt(sys.argv[1:],'hY:B:LCM:l',long_options)

	#-- command line parameters
	YEAR = [2010,2011,2012,2013,2014,2015,2016,2017,2018,2019]
	BASELINE = 'D'
	DIRECTORY = os.getcwd()
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
		elif opt in ("-L","--list"):
			LIST = True
		elif opt in ("-l","--log"):
			LOG = True
		elif opt in ("-M","--mode"):
			MODE = int(arg, 8)
		elif opt in ("-C","--clobber"):
			CLOBBER = True

	if not arglist:
		#-- Input CryoSat Level-2 Product (sys.argv[0] is the python code)
		raise Exception('No CryoSat Level-2 Product Specified')

	#-- check internet connection before attempting to run program
	if check_connection():
		for PRODUCT in arglist:
			esa_cryosat_sync(PRODUCT,YEAR,DIRECTORY=DIRECTORY,BASELINE=BASELINE,
				LOG=LOG,LIST=LIST,MODE=MODE,CLOBBER=CLOBBER)

#-- run main program
if __name__ == '__main__':
	main()
