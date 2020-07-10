/*
 * Copyright 2020 Free Software Foundation, Inc.
 *
 * This file is part of GNU Radio
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 */

/***********************************************************************************/
/* This file is automatically generated using bindtool and can be manually edited  */
/* The following lines can be configured to regenerate this file during cmake      */
/* If manual edits are made, the following tags should be modified accordingly.    */
/* BINDTOOL_GEN_AUTOMATIC(0)                                                       */
/* BINDTOOL_USE_PYGCCXML(0)                                                        */
/* BINDTOOL_HEADER_FILE(symbol_demapper_cf.h)                                      */
/* BINDTOOL_HEADER_FILE_HASH(a0c72d4c039e745d75c78089c2be53c1)                     */
/***********************************************************************************/

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include <symbolmapping/symbol_demapper_cf.h>
// pydoc.h is automatically generated in the build directory
#include <symbol_demapper_cf_pydoc.h>

void bind_symbol_demapper_cf(py::module& m)
{

    using symbol_demapper_cf = ::gr::symbolmapping::symbol_demapper_cf;


    py::class_<symbol_demapper_cf,
               gr::sync_interpolator,
               std::shared_ptr<symbol_demapper_cf>>(
        m, "symbol_demapper_cf", D(symbol_demapper_cf))

        .def(py::init(&symbol_demapper_cf::make),
             py::arg("constellation_order"),
             py::arg("constellation_type"),
             py::arg("snr_tag_name") = std::string("snr"),
             D(symbol_demapper_cf, make))


        ;
}