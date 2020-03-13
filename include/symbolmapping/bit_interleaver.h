/* -*- c++ -*- */
/*
 * Copyright 2019 Johannes Demel.
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifndef SYMBOLMAPPING_BIT_INTERLEAVER_H
#define SYMBOLMAPPING_BIT_INTERLEAVER_H

#include <cstddef>
#include <vector>
#include <iostream>


class BitInterleaver
{
public:
    BitInterleaver(std::vector<size_t> interleaver_indices);
    ~BitInterleaver();

    size_t interleaverLength(){ return _interleaver_indices.size();};
    std::vector<size_t> interleaverIndices(){return _interleaver_indices;};
    std::vector<size_t> deinterleaverIndices(){return _deinterleaver_indices;};

    template <class T>
    void interleave(T* target, const T* src){
      for(auto idx : _interleaver_indices){
        *target++ = src[idx];
      }
    }

    template <class T>
    void deinterleave(T* target, const T* src){
      for(auto idx : _deinterleaver_indices){
        *target++ = src[idx];
      }
    }

private:
    std::vector<size_t> _interleaver_indices;
    std::vector<size_t> _deinterleaver_indices;

    void set_interleaver(const std::vector<size_t> &interleaver_indices);

};


#endif //SYMBOLMAPPING_BIT_INTERLEAVER_H
