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
/* BINDTOOL_HEADER_FILE(symbol_mapper_bc.h)                                        */
/* BINDTOOL_HEADER_FILE_HASH(f062549564e429bd1765ee8da408d007)                     */
/***********************************************************************************/

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include <symbolmapping/symbol_mapper_bc.h>
// pydoc.h is automatically generated in the build directory
#include <symbol_mapper_bc_pydoc.h>

void bind_symbol_mapper_bc(py::module& m)
{

    using symbol_mapper_bc = ::gr::symbolmapping::symbol_mapper_bc;


    py::class_<symbol_mapper_bc,
               gr::block,
               gr::basic_block,
               std::shared_ptr<symbol_mapper_bc>>(
        m, "symbol_mapper_bc", D(symbol_mapper_bc))

        .def(py::init(&symbol_mapper_bc::make),
             py::arg("constellation_order"),
             py::arg("constellation_type"),
             py::arg("is_packed"),
             D(symbol_mapper_bc, make))


        ;
}
