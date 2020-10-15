============
utilities.py
============

Download and management utilities for syncing time and auxiliary files

 - Can list a directory on a ftp host
 - Can download a file from a ftp or http host
 - Checks MD5 hashes between local and remote files

`Source code`__

.. __: https://github.com/tsutterley/read-cryosat-2/blob/master/cryosat_toolkit/utilities.py


General Methods
===============


.. method:: cryosat_toolkit.utilities.get_data_path(relpath)

    Get the absolute path within a package from a relative path

    Arguments:

        `relpath`: local relative path as list or string


.. method:: cryosat_toolkit.utilities.get_hash(local)

    Get the MD5 hash value from a local file

    Arguments:

        `local`: path to file


.. method:: cryosat_toolkit.utilities.get_unix_time(time_string, format='%Y-%m-%d %H:%M:%S')

    Get the Unix timestamp value for a formatted date string

    Arguments:

        `time_string`: formatted time string to parse

    Keyword arguments:

        `format`: format for input time string


.. method:: cryosat_toolkit.utilities.ftp_list(HOST,username,password,timeout=None,basename=False,pattern=None,sort=False)

    List a directory on a ftp host

    Arguments:

        `HOST`: remote ftp host path split as list

    Keyword arguments:

        `username`: ftp server username

        `password`: ftp server password

        `timeout`: timeout in seconds for blocking operations

        `basename`: return the file or directory basename instead of the full path

        `pattern`: regular expression pattern for reducing list

        `sort`: sort output list

    Returns:

        `output`: list of items in a directory

        `mtimes`: list of last modification times for items in the directory


.. method:: cryosat_toolkit.utilities.from_ftp(HOST,username,password,timeout=None,local=None,hash='',chunk=16384,verbose=False,mode=0o775)

    Download a file from a ftp host

    Arguments:

        `HOST`: remote ftp host path split as list

    Keyword arguments:

        `username`: ftp server username

        `password`: ftp server password

        `timeout`: timeout in seconds for blocking operations

        `local`: path to local file

        `hash`: MD5 hash of local file

        `chunk`: chunk size for transfer encoding

        `verbose`: print file transfer information

        `mode`: permissions mode of output local file


.. method:: cryosat_toolkit.utilities.from_http(HOST,timeout=None,local=None,hash='',chunk=16384,verbose=False,mode=0o775)

    Download a file from a http host

    Arguments:

        `HOST`: remote http host path split as list

    Keyword arguments:

        `timeout`: timeout in seconds for blocking operations

        `local`: path to local file

        `hash`: MD5 hash of local file

        `chunk`: chunk size for transfer encoding

        `verbose`: print file transfer information

        `mode`: permissions mode of output local file

.. method:: cryosat_toolkit.utilities.build_opener(context=ssl.SSLContext())

    build urllib opener for ESA CryoSat-2 Science Server

    Keyword arguments:

        `context`: SSL context for opener object
