# Copyright 2020 Free Software Foundation, Inc.
#
# This file is part of GNU Radio
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

if(Gnuradio_FOUND)
    GR_PYTHON_CHECK_MODULE_RAW(
        "pygccxml"
        "import pygccxml"
        PYGCCXML_FOUND
        )

    include(GrPybind)
endif()

########################################################################
# Python Bindings
########################################################################

list(APPEND symbolmapping_python_files
    ../../lib/bit_interleaver.cc
    bit_interleaver_python.cc
    ../../lib/symbol_mapper.cc
    symbol_mapper_python.cc
    python_bindings.cc)

if(Gnuradio_FOUND)
    list(APPEND symbolmapping_python_files
        interleaver_python.cc
        symbol_demapper_cf_python.cc
        symbol_mapper_bc_python.cc)

    GR_PYBIND_MAKE_OOT(symbolmapping
    ../..
    gr::symbolmapping
    "${symbolmapping_python_files}")
else()
    add_compile_definitions("SKIP_GNURADIO")
    pybind11_add_module(symbolmapping_python
        "${symbolmapping_python_files}"
    )

    target_include_directories(symbolmapping_python PUBLIC
        ${PYTHON_NUMPY_INCLUDE_DIR}
        ${CMAKE_SOURCE_DIR}/include
        ${PYBIND11_INCLUDE_DIR}
    )
    target_link_libraries(symbolmapping_python PRIVATE Volk::volk ${ORC_LIBRARIES})
endif()

install(TARGETS symbolmapping_python DESTINATION ${GR_PYTHON_DIR}/symbolmapping COMPONENT pythonapi)
