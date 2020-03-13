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

#include <symbolmapping/bit_interleaver.h>
#include <cstring>
#include <chrono>
#include <cstdlib>
#include <algorithm>
#include <numeric>

BitInterleaver::BitInterleaver(std::vector<size_t> interleaver_indices)
{
  set_interleaver(interleaver_indices);
}

BitInterleaver::~BitInterleaver(){}

void BitInterleaver::set_interleaver(const std::vector<size_t> &interleaver_indices){
  _interleaver_indices = std::vector<size_t>(interleaver_indices);

  _deinterleaver_indices.resize(_interleaver_indices.size());
  std::iota(_deinterleaver_indices.begin(), _deinterleaver_indices.end(), 0);
  // sort indexes based on comparing values in v
  std::sort(_deinterleaver_indices.begin(), _deinterleaver_indices.end(),
            [&interleaver_indices](size_t i1, size_t i2) {
              return interleaver_indices[i1] < interleaver_indices[i2];
            });
}
