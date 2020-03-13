
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

#include <symbolmapping/bit_interleaver.h>


namespace py = pybind11;


void bind_bit_interleaver(py::module& m)
{
    py::class_<BitInterleaver>(m, "Interleaver")

        .def(py::init<std::vector<size_t> >())
        .def("interleaverLength", &BitInterleaver::interleaverLength)
        .def("interleaverIndices", &BitInterleaver::interleaverIndices)
        .def("deinterleaverIndices", &BitInterleaver::deinterleaverIndices)
        .def("interleave",  [](BitInterleaver& self, const py::array_t<float, py::array::c_style | py::array::forcecast> array) {
            py::buffer_info inb = array.request();
            if(inb.ndim != 1){
                throw std::runtime_error("Only ONE-dimensional vectors allowed!");
            }
            auto result =py::array_t<float>(inb.size);
            py::buffer_info resb = result.request();
            
            self.interleave<float>((float*) resb.ptr, (float*) inb.ptr);
            return result;
        })
        .def("interleave",  [](BitInterleaver& self, const py::array_t<uint8_t, py::array::c_style | py::array::forcecast> array) {
            py::buffer_info inb = array.request();
            if(inb.ndim != 1){
                throw std::runtime_error("Only ONE-dimensional vectors allowed!");
            }
            auto result =py::array_t<uint8_t>(inb.size);
            py::buffer_info resb = result.request();
            
            self.interleave<uint8_t>((uint8_t*) resb.ptr, (uint8_t*) inb.ptr);
            return result;
        })

        .def("deinterleave",  [](BitInterleaver& self, const py::array_t<float, py::array::c_style | py::array::forcecast> array) {
            py::buffer_info inb = array.request();
            if(inb.ndim != 1){
                throw std::runtime_error("Only ONE-dimensional vectors allowed!");
            }
            auto result =py::array_t<float>(inb.size);
            py::buffer_info resb = result.request();
            
            self.deinterleave<float>((float*) resb.ptr, (float*) inb.ptr);
            return result;
        })
        .def("deinterleave",  [](BitInterleaver& self, const py::array_t<uint8_t, py::array::c_style | py::array::forcecast> array) {
            py::buffer_info inb = array.request();
            if(inb.ndim != 1){
                throw std::runtime_error("Only ONE-dimensional vectors allowed!");
            }
            auto result =py::array_t<uint8_t>(inb.size);
            py::buffer_info resb = result.request();
            
            self.deinterleave<uint8_t>((uint8_t*) resb.ptr, (uint8_t*) inb.ptr);
            return result;
        })


        ;
} 
