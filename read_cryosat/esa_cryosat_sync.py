#!/usr/bin/env python
u"""
esa_cryosat_sync.py
Written by Tyler Sutterley (05/2017)

This program syncs Cryosat Elevation products
From the ESA Cryosat dissemination server:
https://earth.esa.int/web/guest/-/how-to-access-cryosat-data-6842
https://earth.esa.int/web/guest/-/products-overview-6975

INPUTS:
	CryoSat-2 product to sync with ESA servers
		SIR_LRM_L2: CryoSat-2 Low-Resolution Mode
		SIR_SAR_L2: CryoSat-2 SAR Mode
		SIR_SIN_L2: CryoSat-2 SARin Mode

CALLING SEQUENCE:
	python esa_cryosat_sync.py --baseline=C --user=<username> SIR_SIN_L2
	where <username> is your ESA data dissemination server username

COMMAND LINE OPTIONS:
	--help: list the command line options
	-Y X, --year=X: years to sync separated by commas
	-B X, --baseline=X: CryoSat-2 baseline to sync
	--user: username for CryoSat-2 FTP servers
	--directory: working data directory (default: current working directory)
	-M X, --mode=X: Local permissions mode of the directories and files synced
	--log: output log of files downloaded
	--list: print files to be transferred, but do not execute transfer
	--clobber: Overwrite existing data in transfer

UPDATE HISTORY:
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
import getopt
import getpass
import calendar, time
import ftplib

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
def compile_regex_pattern(PRODUCT):
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
	#-- CRYOSAT LEVEL-2 PRODUCTS NAMING RULES
	#-- Mission Identifier
	#-- File Class
	#-- File Product
	#-- Validity Start Date and Time
	#-- Validity Stop Date and Time
	#-- Baseline Identifier
	#-- Version Number
	regex_pattern = '({0})_({1})_({2})__(\d+T?\d+)_(\d+T?\d+)_(.*?)(\d+).{3}' \
		.format('CS', regex_class, regex_products[PRODUCT], '(DBL|HDR)')
	return re.compile(regex_pattern, re.VERBOSE)

#-- PURPOSE: sync local Cryosat-2 files with ESA server
def esa_cryosat_sync(PRODUCT, YEARS, DIRECTORY=None, USER='', PASSWORD='',
	LOG=False, LIST=False, MODE=None, CLOBBER=False):
	#-- connect and login to ESA ftp server
	f = ftplib.FTP('science-pds.cryosat.esa.int')
	f.login(USER, PASSWORD)

	#-- create log file with list of synchronized files (or print to terminal)
	if LOG:
		#-- output to log file
		LOGDIR = os.path.join(DIRECTORY,'sync_logs.dir')
		#-- check if log directory exists and recursively create if not
		os.makedirs(LOGDIR,MODE) if not os.path.exists(LOGDIR) else None
		#-- format: PODAAC_sync_2002-04-01.log
		today = time.strftime('%Y-%m-%d',time.localtime())
		LOGFILE = 'ESA_CS_{0}_sync_{1}.log'.format(PRODUCT,today)
		fid1 = open(os.path.join(LOGDIR,LOGFILE),'w')
		print('ESA CryoSat-2 Sync Log ({0})'.format(today), file=fid1)
		print('PRODUCT={0}'.format(PRODUCT), file=fid1)
	else:
		#-- standard output (terminal output)
		fid1 = sys.stdout

	#-- compile regular expression operator for years to sync
	regex_years = '|'.join('{0:d}'.format(y) for y in YEARS)
	R1 = re.compile('({0})'.format(regex_years), re.VERBOSE)
	#-- initial regular expression pattern for months of the year
	regex_months = '|'.join('{0:02d}'.format(m) for m in range(1,13))

	#-- find remote yearly directories for PRODUCT within YEARS
	YRS = [R1.findall(Y).pop() for Y in f.nlst(PRODUCT) if R1.search(Y)]
	for Y in YRS:
		#-- compile regular expression operator for months in year to sync
		R2 = re.compile('{0}/({1})'.format(Y,regex_months), re.VERBOSE)
		#-- find remote monthly directories for PRODUCT within year
		MNS = [R2.findall(M).pop() for M in f.nlst('{0}/{1}'.format(PRODUCT,Y))
			if R2.search(M)]
		for M in MNS:
			#-- remote and local directory for data product of year and month
			remote_dir = '{0}/{1}/{2}/'.format(PRODUCT,Y,M)
			local_dir = os.path.join(DIRECTORY,PRODUCT,Y,M)
			#-- check if local directory exists and recursively create if not
			os.makedirs(local_dir,MODE) if not os.path.exists(local_dir) else None
			#-- compile the regular expression operator to find CryoSat-2 files
			R3 = compile_regex_pattern(PRODUCT)
			#-- get filenames from remote directory
			valid_lines = [fi for fi in f.nlst(remote_dir) if R3.search(fi)]
			for line in sorted(valid_lines):
				#-- extract filename from regex object
				fi = R3.search(line).group(0)
				remote_file = '{0}{1}'.format(remote_dir,fi)
				local_file = os.path.join(local_dir,fi)
				ftp_mirror_file(fid1,f,remote_file,local_file,LIST,CLOBBER,MODE)

	#-- close the ftp connection
	f.quit()
	#-- close log file and set permissions level to MODE
	if LOG:
		fid1.close()
		os.chmod(os.path.join(LOGDIR,LOGFILE), MODE)

#-- PURPOSE: pull file from a remote host checking if file exists locally
#-- and if the remote file is newer than the local file
def ftp_mirror_file(fid, ftp, remote_file, local_file, LIST, CLOBBER, MODE):
	#-- if file exists in file system: check if remote file is newer
	TEST = False
	OVERWRITE = ' (clobber)'
	#-- get last modified date of remote file and convert into unix time
	mdtm = ftp.sendcmd('MDTM {0}'.format(remote_file))
	remote_mtime = calendar.timegm(time.strptime(mdtm[4:],"%Y%m%d%H%M%S"))
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
		print('{0}{1}/{2} --> '.format('ftp://',ftp.host,remote_file),file=fid)
		print('\t{0}{1}\n'.format(local_file,OVERWRITE),file=fid)
		#-- if executing copy command (not only printing the files)
		if not LIST:
			#-- copy remote file contents to local file
			with open(local_file, 'wb') as f:
				ftp.retrbinary('RETR {0}'.format(remote_file), f.write)
			#-- keep remote modification time of file and local access time
			os.utime(local_file, (os.stat(local_file).st_atime, remote_mtime))
			os.chmod(local_file, MODE)

#-- PURPOSE: help module to describe the optional input parameters
def usage():
	print('\nHelp: {}'.format(os.path.basename(sys.argv[0])))
	print(' -Y X, --year=X\t\tYears to sync separated by commas')
	print(' -B X, --baseline=X\tCryoSat Baseline to run')
	print(' -U X, --user=X\t\tUsername for CryoSat-2 FTP servers')
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
	long_options = ['help','year=','baseline=','user=','directory=','list',
		'log','mode=','clobber']
	optlist,arglist = getopt.getopt(sys.argv[1:],'hY:B:U:LCM:l',long_options)

	#-- command line parameters
	years = [2010,2011,2012,2013,2014,2015,2016,2017,2018]
	BASELINE = 'C'
	USER = ''
	DIRECTORY = os.getcwd()
	LIST = False
	LOG = False
	#-- permissions mode of the local directories and files (number in octal)
	MODE = 0775
	CLOBBER = False
	for opt, arg in optlist:
		if opt in ('-h','--help'):
			usage()
			sys.exit()
		elif opt in ("-Y","--year"):
			years = [int(Y) for Y in arg.split(',')]
		elif opt in ("-B","--baseline"):
			BASELINE = arg.upper()
		elif opt in ("-U","--user"):
			USER = arg
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

	#-- check that ESA CryoSat-2 FTP Server credentials were entered
	if not USER:
		raise IOError('Please Enter your ESA CryoSat-2 FTP Server Username')
	#-- enter password securely from command-line
	HOST = 'science-pds.cryosat.esa.int'
	PASSWORD = getpass.getpass('Password for {0}@{1}: '.format(USER,HOST))

	#-- check internet connection before attempting to run program
	if check_connection(USER,PASSWORD):
		for PRODUCT in arglist:
			esa_cryosat_sync(PRODUCT, years, DIRECTORY=DIRECTORY, USER=USER,
				PASSWORD=PASSWORD,LOG=LOG,LIST=LIST,MODE=MODE,CLOBBER=CLOBBER)

#-- run main program
if __name__ == '__main__':
	main()
