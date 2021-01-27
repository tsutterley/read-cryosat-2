"""
utilities.py
Written by Tyler Sutterley (08/2020)
Download and management utilities for syncing time and auxiliary files

UPDATE HISTORY:
    Written 08/2020
"""
from __future__ import print_function

import sys
import os
import re
import io
import ssl
import ftplib
import shutil
import base64
import socket
import inspect
import hashlib
import posixpath
import lxml.etree
import calendar,time
if sys.version_info[0] == 2:
    from cookielib import CookieJar
    from urllib import urlencode
    import urllib2
else:
    from http.cookiejar import CookieJar
    from urllib.parse import urlencode
    import urllib.request as urllib2

def get_data_path(relpath):
    """
    Get the absolute path within a package from a relative path

    Arguments
    ---------
    relpath: relative path
    """
    #-- current file path
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    filepath = os.path.dirname(os.path.abspath(filename))
    if isinstance(relpath,list):
        #-- use *splat operator to extract from list
        return os.path.join(filepath,*relpath)
    elif isinstance(relpath,str):
        return os.path.join(filepath,relpath)

#-- PURPOSE: get the MD5 hash value of a file
def get_hash(local):
    """
    Get the MD5 hash value from a local file

    Arguments
    ---------
    local: path to file
    """
    #-- check if local file exists
    if os.access(os.path.expanduser(local),os.F_OK):
        #-- generate checksum hash for local file
        #-- open the local_file in binary read mode
        with open(os.path.expanduser(local), 'rb') as local_buffer:
            return hashlib.md5(local_buffer.read()).hexdigest()
    else:
        return ''

#-- PURPOSE: returns the Unix timestamp value for a formatted date string
def get_unix_time(time_string, format='%Y-%m-%d %H:%M:%S'):
    """
    Get the Unix timestamp value for a formatted date string

    Arguments
    ---------
    time_string: formatted time string to parse

    Keyword arguments
    -----------------
    format: format for input time string
    """
    try:
        parsed_time = time.strptime(time_string.rstrip(), format)
    except (TypeError, ValueError):
        return None
    else:
        return calendar.timegm(parsed_time)

#-- PURPOSE: list a directory on a ftp host
def ftp_list(HOST,username=None,password=None,timeout=None,basename=False,
    pattern=None,sort=False):
    """
    List a directory on a ftp host

    Arguments
    ---------
    HOST: remote ftp host path split as list

    Keyword arguments
    -----------------
    username: ftp server username
    password: ftp server password
    timeout: timeout in seconds for blocking operations
    basename: return the file or directory basename instead of the full path
    pattern: regular expression pattern for reducing list
    sort: sort output list

    Returns
    -------
    output: list of items in a directory
    mtimes: list of last modification times for items in the directory
    """
    #-- try to connect to ftp host
    try:
        ftp = ftplib.FTP(HOST[0],timeout=timeout)
        ftp.login(username,password)
        ftp.voidcmd("NOOP")
    except (socket.gaierror,IOError):
        raise RuntimeError('Unable to connect to {0}'.format(HOST[0]))
    except ftplib.error_perm:
        raise RuntimeError('Check login credentials')
    else:
        #-- list remote path
        output = ftp.nlst(posixpath.join(*HOST[1:]))
        #-- get last modified date of ftp files and convert into unix time
        mtimes = [None]*len(output)
        #-- iterate over each file in the list and get the modification time
        for i,f in enumerate(output):
            try:
                #-- try sending modification time command
                mdtm = ftp.sendcmd('MDTM {0}'.format(f))
            except ftplib.error_perm:
                #-- directories will return with an error
                pass
            else:
                #-- convert the modification time into unix time
                mtimes[i] = get_unix_time(time.strptime(mdtm[4:],
                    format="%Y%m%d%H%M%S"))
        #-- reduce to basenames
        if basename:
            output = [posixpath.basename(i) for i in output]
        #-- reduce using regular expression pattern
        if pattern:
            i = [i for i,f in enumerate(output) if re.search(pattern,f)]
            #-- reduce list of listed items and last modified times
            output = [output[indice] for indice in i]
            mtimes = [mtimes[indice] for indice in i]
        #-- sort the list
        if sort:
            i = [i for i,j in sorted(enumerate(output), key=lambda i: i[1])]
            #-- sort list of listed items and last modified times
            output = [output[indice] for indice in i]
            mtimes = [mtimes[indice] for indice in i]
        #-- close the ftp connection
        ftp.close()
        #-- return the list of items and last modified times
        return (output,mtimes)

#-- PURPOSE: download a file from a ftp host
def from_ftp(HOST,username=None,password=None,timeout=None,local=None,hash='',
    chunk=8192,verbose=False,mode=0o775):
    """
    Download a file from a ftp host

    Arguments
    ---------
    HOST: remote ftp host path split as list

    Keyword arguments
    -----------------
    username: ftp server username
    password: ftp server password
    timeout: timeout in seconds for blocking operations
    local: path to local file
    hash: MD5 hash of local file
    chunk: chunk size for transfer encoding
    verbose: print file transfer information
    mode: permissions mode of output local file

    Returns
    -------
    remote_buffer: BytesIO representation of file
    """
    #-- try downloading from ftp
    try:
        #-- try to connect to ftp host
        ftp = ftplib.FTP(HOST[0],timeout=timeout)
        ftp.login(username,password)
        ftp.voidcmd("NOOP")
    except (socket.gaierror,IOError):
        raise RuntimeError('Unable to connect to {0}'.format(HOST[0]))
    except ftplib.error_perm:
        raise RuntimeError('Check login credentials')
    else:
        #-- remote path
        ftp_remote_path = posixpath.join(*HOST[1:])
        #-- copy remote file contents to bytesIO object
        remote_buffer = io.BytesIO()
        ftp.retrbinary('RETR {0}'.format(ftp_remote_path),
            remote_buffer.write, blocksize=chunk)
        remote_buffer.seek(0)
        #-- save file basename with bytesIO object
        remote_buffer.filename = HOST[-1]
        #-- generate checksum hash for remote file
        remote_hash = hashlib.md5(remote_buffer.getvalue()).hexdigest()
        #-- compare checksums
        if local and (hash != remote_hash):
            #-- print file information
            if verbose:
                print('{0} -->\n\t{1}'.format(posixpath.join(*HOST),local))
            #-- store bytes to file using chunked transfer encoding
            remote_buffer.seek(0)
            with open(os.path.expanduser(local), 'wb') as f:
                shutil.copyfileobj(remote_buffer, f, chunk)
            #-- change the permissions mode
            os.chmod(local,mode)
        #-- close the ftp connection
        ftp.close()
        #-- return the bytesIO object
        remote_buffer.seek(0)
        return remote_buffer

#-- PURPOSE: download a file from a http host
def from_http(HOST,timeout=None,local=None,hash='',chunk=16384,
    verbose=False,mode=0o775):
    """
    Download a file from a http host

    Arguments
    ---------
    HOST: remote http host path split as list

    Keyword arguments
    -----------------
    timeout: timeout in seconds for blocking operations
    local: path to local file
    hash: MD5 hash of local file
    chunk: chunk size for transfer encoding
    verbose: print file transfer information
    mode: permissions mode of output local file

    Returns
    -------
    remote_buffer: BytesIO representation of file
    """
    #-- try downloading from http
    try:
        #-- Create and submit request.
        request = urllib2.Request(posixpath.join(*HOST))
        response = urllib2.urlopen(request,timeout=timeout,context=ssl.SSLContext())
    except (urllib2.HTTPError, urllib2.URLError):
        raise Exception('Download error from {0}'.format(posixpath.join(*HOST)))
    else:
        #-- copy remote file contents to bytesIO object
        remote_buffer = io.BytesIO(response.read())
        remote_buffer.seek(0)
        #-- save file basename with bytesIO object
        remote_buffer.filename = HOST[-1]
        #-- generate checksum hash for remote file
        remote_hash = hashlib.md5(remote_buffer.getvalue()).hexdigest()
        #-- compare checksums
        if local and (hash != remote_hash):
            #-- print file information
            if verbose:
                print('{0} -->\n\t{1}'.format(posixpath.join(*HOST),local))
            #-- store bytes to file using chunked transfer encoding
            remote_buffer.seek(0)
            with open(os.path.expanduser(local), 'wb') as f:
                shutil.copyfileobj(remote_buffer, f, chunk)
            #-- change the permissions mode
            os.chmod(local,mode)
        #-- return the bytesIO object
        remote_buffer.seek(0)
        return remote_buffer

#-- PURPOSE: create opener for downloading from ESA https server
def build_opener(context=ssl.SSLContext()):
    """
    build urllib opener for ESA CryoSat-2 Science Server

    Keyword arguments
    -----------------
    context: SSL context for opener object
    """
    #-- https://docs.python.org/3/howto/urllib2.html#id5
    handler = []
    #-- Create cookie jar for storing cookies
    cookie_jar = CookieJar()
    handler.append(urllib2.HTTPCookieProcessor(cookie_jar))
    #-- SSL context handler
    handler.append(urllib2.HTTPSHandler(context=context))
    #-- create "opener" (OpenerDirector instance)
    opener = urllib2.build_opener(*handler)
    #-- Now all calls to urllib2.urlopen use our opener.
    urllib2.install_opener(opener)
    #-- All calls to urllib2.urlopen will now use handler
    #-- Make sure not to include the protocol in with the URL, or
    #-- HTTPPasswordMgrWithDefaultRealm will be confused.