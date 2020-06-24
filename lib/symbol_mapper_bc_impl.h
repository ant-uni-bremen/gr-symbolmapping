/* -*- c++ -*- */
/*
 * Copyright 2020 Johannes Demel.
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

#ifndef INCLUDED_SYMBOLMAPPING_SYMBOL_MAPPER_BC_IMPL_H
#define INCLUDED_SYMBOLMAPPING_SYMBOL_MAPPER_BC_IMPL_H

#include <gnuradio/blocks/unpack_k_bits.h>
#include <symbolmapping/symbol_mapper.h>
#include <symbolmapping/symbol_mapper_bc.h>

#include <volk/volk_alloc.hh>

#include <memory>
#include <vector>

namespace gr {
namespace symbolmapping {

class symbol_mapper_bc_impl : public symbol_mapper_bc
{
private:
    std::unique_ptr<gr::blocks::kernel::unpack_k_bits> d_unpacker;
    std::unique_ptr<SymbolMapping> d_mapper;

    volk::vector<unsigned char> d_unpacked_bits;
    void update_unpacked_vector_length(const unsigned max_length);

    bool d_is_packed;

public:
    symbol_mapper_bc_impl(unsigned constellation_order,
                          std::string constellation_type,
                          bool is_packed);
    ~symbol_mapper_bc_impl();

    // Where all the action really happens
    void forecast(int noutput_items, gr_vector_int& ninput_items_required);

    int general_work(int noutput_items,
                     gr_vector_int& ninput_items,
                     gr_vector_const_void_star& input_items,
                     gr_vector_void_star& output_items);

    int fixed_rate_ninput_to_noutput(int ninput);
    int fixed_rate_noutput_to_ninput(int noutput);
};

} // namespace symbolmapping
} // namespace gr

#endif /* INCLUDED_SYMBOLMAPPING_SYMBOL_MAPPER_BC_IMPL_H */
