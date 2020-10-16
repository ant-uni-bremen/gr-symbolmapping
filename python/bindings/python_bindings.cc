/*
 * Copyright 2020 Johannes Demel
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 */

#include <pybind11/pybind11.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

namespace py = pybind11;

// Headers for binding functions
/**************************************/
/* The following comment block is used for
/* gr_modtool to insert function prototypes
/* Please do not delete
/**************************************/
// BINDING_FUNCTION_PROTOTYPES(
void bind_bit_interleaver(py::module& m);
void bind_symbol_mapper(py::module& m);

void bind_interleaver(py::module& m);
void bind_symbol_demapper_cf(py::module& m);
void bind_symbol_mapper_bc(py::module& m);
// ) END BINDING_FUNCTION_PROTOTYPES


// We need this hack because import_array() returns NULL
// for newer Python versions.
// This function is also necessary because it ensures access to the C API
// and removes a warning.
void* init_numpy()
{
    import_array();
    return NULL;
}

PYBIND11_MODULE(symbolmapping_python, m)
{
    // Initialize the numpy C API
    // (otherwise we will see segmentation faults)
    init_numpy();

    // Allow access to base block methods
#ifndef SKIP_GNURADIO
    py::module::import("gnuradio.gr");
#endif
    /**************************************/
    /* The following comment block is used for
    /* gr_modtool to insert binding function calls
    /* Please do not delete
    /**************************************/
    // BINDING_FUNCTION_CALLS(
    bind_bit_interleaver(m);
    bind_symbol_mapper(m);

#ifndef SKIP_GNURADIO
    bind_interleaver(m);
    bind_symbol_demapper_cf(m);
    bind_symbol_mapper_bc(m);
#endif
    // ) END BINDING_FUNCTION_CALLS
}
