#
# Copyright 2008,2009 Free Software Foundation, Inc.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

# The presence of this file turns this directory into a Python package

'''
This is the GNU Radio SYMBOLMAPPING module. Place your Python package
description here (python/__init__.py).
'''
from __future__ import unicode_literals

# import swig generated symbols into the symbolmapping namespace
try:
    # this might fail if the module is python-only
    from .symbolmapping_swig import *
except ImportError:
    pass

try:
    from .symbolmapping_python import *
except ImportError:
    pass

# import any pure python here
#
from .symbol_constellation import generate_constellation, map_to_constellation
from .interleaver_indices import create_interleaver_indices