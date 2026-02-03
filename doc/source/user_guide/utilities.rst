============
utilities.py
============

Download and management utilities for syncing time and auxiliary files

 - Can list a directory on a ftp host
 - Can download a file from a ftp or http host
 - Checks ``MD5`` or ``sha1`` hashes between local and remote files

`Source code`__

.. __: https://github.com/tsutterley/read-cryosat-2/blob/main/cryosat_toolkit/utilities.py


General Methods
===============

.. autofunction:: cryosat_toolkit.utilities.get_data_path

.. autofunction:: cryosat_toolkit.utilities.get_hash

.. autofunction:: cryosat_toolkit.utilities.url_split

.. autofunction:: cryosat_toolkit.utilities.convert_arg_line_to_args

.. autofunction:: cryosat_toolkit.utilities.get_unix_time

.. autofunction:: cryosat_toolkit.utilities.isoformat

.. autofunction:: cryosat_toolkit.utilities.even

.. autofunction:: cryosat_toolkit.utilities.ceil

.. autofunction:: cryosat_toolkit.utilities.copy

.. autofunction:: cryosat_toolkit.utilities.check_ftp_connection

.. autofunction:: cryosat_toolkit.utilities.ftp_list

.. autofunction:: cryosat_toolkit.utilities.from_ftp

.. autofunction:: cryosat_toolkit.utilities.check_connection

.. autofunction:: cryosat_toolkit.utilities.http_list

.. autofunction:: cryosat_toolkit.utilities.from_http

.. autofunction:: cryosat_toolkit.utilities.build_opener
