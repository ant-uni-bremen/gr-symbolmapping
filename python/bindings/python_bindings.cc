
/*
 * Copyright 2020 Johannes Demel
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 */

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>


namespace py = pybind11;


void bind_bit_interleaver(py::module& m);
void bind_symbol_mapper(py::module& m);


PYBIND11_MODULE(symbolmapping_python, m) {
    bind_bit_interleaver(m);
    bind_symbol_mapper(m);

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif

}
