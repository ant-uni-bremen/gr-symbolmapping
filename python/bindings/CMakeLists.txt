# Copyright 2020 Free Software Foundation, Inc.
#
# This file is part of GNU Radio
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

GR_PYTHON_CHECK_MODULE_RAW(
    "pygccxml"
    "import pygccxml"
    PYGCCXML_FOUND
    )

include(GrPybind)

########################################################################
# Python Bindings
########################################################################

list(APPEND symbolmapping_python_files
    ../../lib/bit_interleaver.cc
    bit_interleaver_python.cc
    interleaver_python.cc
    symbol_demapper_cf_python.cc
    ../../lib/symbol_mapper.cc
    symbol_mapper_python.cc
    symbol_mapper_bc_python.cc
    python_bindings.cc)

GR_PYBIND_MAKE_OOT(symbolmapping
   ../..
   gr::symbolmapping
   "${symbolmapping_python_files}")

install(TARGETS symbolmapping_python DESTINATION ${GR_PYTHON_DIR}/symbolmapping COMPONENT pythonapi)
