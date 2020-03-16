
/*
 * Copyright 2020 Johannes Demel
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 */

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include <cstdint>

#include <symbolmapping/symbol_mapper.h>


namespace py = pybind11;


void bind_symbol_mapper(py::module& m)
{
    m.def("lin2db", &lin2db);
    m.def("db2lin", &db2lin);
    py::class_<SymbolMapping>(m, "SymbolMapping")

        .def(py::init<unsigned, std::string>(),
             py::arg("constellation_order") = 2,
             py::arg("cstl_type") = std::string("GRAY"))
        .def("constellationOrder", &SymbolMapping::constellationOrder)
        .def("constellationSize", &SymbolMapping::constellationSize)
        .def("constellationType", &SymbolMapping::constellationType)
        .def("setConstellationOrder", &SymbolMapping::setConstellationOrder)
        .def("constellation",  [](SymbolMapping& self) {
                auto c = self.constellation();
                return py::array(c.size(), c.data());
            })

        .def("map_to_constellation",  [](SymbolMapping& self, 
                               const py::array_t<uint8_t, 
                                                 py::array::c_style | py::array::forcecast> array) {
            py::buffer_info inb = array.request();
            if(inb.ndim != 1){
                throw std::runtime_error("Only ONE-dimensional vectors allowed!");
            }
            auto result =py::array_t<fcmplx>(inb.size / self.constellationOrder());
            py::buffer_info resb = result.request();
            
            self.map_to_constellation((fcmplx*) resb.ptr, (uint8_t*) inb.ptr, inb.size);
            return result;
        })
        .def("calculate_ln_probabilities",  [](SymbolMapping& self, 
                const py::array_t<fcmplx, py::array::c_style | py::array::forcecast> array, 
                const float snr_db) {
            py::buffer_info inb = array.request();
            if(inb.ndim != 1){
                throw std::runtime_error("Only ONE-dimensional vectors allowed!");
            }
            auto result =py::array_t<float>(inb.size * self.constellationSize());
            py::buffer_info resb = result.request();
            
            self.calculate_ln_probabilities((float*) resb.ptr, 
                                            (fcmplx*) inb.ptr, 
                                            inb.size, snr_db);
            return result;
        })
        .def("demap_llrs",  [](SymbolMapping& self, 
                const py::array_t<fcmplx, py::array::c_style | py::array::forcecast> array, 
                const float snr_db) {
            py::buffer_info inb = array.request();
            if(inb.ndim != 1){
                throw std::runtime_error("Only ONE-dimensional vectors allowed!");
            }
            auto result =py::array_t<float>(inb.size * self.constellationOrder());
            py::buffer_info resb = result.request();
            
            self.demap_llrs((float*) resb.ptr, (fcmplx*) inb.ptr, 
                            inb.size, snr_db);
            return result;
        })

        ;
} 
