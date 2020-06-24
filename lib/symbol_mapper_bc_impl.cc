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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "symbol_mapper_bc_impl.h"
#include <gnuradio/io_signature.h>
#include <algorithm>

namespace gr {
namespace symbolmapping {

symbol_mapper_bc::sptr symbol_mapper_bc::make(unsigned constellation_order,
                                              std::string constellation_type,
                                              bool is_packed)
{
    return gnuradio::get_initial_sptr(
        new symbol_mapper_bc_impl(constellation_order, constellation_type, is_packed));
}


/*
 * The private constructor
 */
symbol_mapper_bc_impl::symbol_mapper_bc_impl(unsigned constellation_order,
                                             std::string constellation_type,
                                             bool is_packed)
    : gr::block("symbol_mapper_bc",
                gr::io_signature::make(1, 1, sizeof(unsigned char)),
                gr::io_signature::make(1, 1, sizeof(gr_complex))),
      d_unpacker(new gr::blocks::kernel::unpack_k_bits(8)),
      d_mapper(new SymbolMapping(constellation_order, constellation_type)),
      d_is_packed(is_packed)
{
    if (d_is_packed) {
        update_unpacked_vector_length(1024);
        set_relative_rate(8.0f / d_mapper->constellationOrder());
    } else {
        update_unpacked_vector_length(1);
        set_relative_rate(1.0f / d_mapper->constellationOrder());
    }

    set_fixed_rate(true);

    if (d_is_packed) {
        switch (constellation_order) {
        case 2:
            set_output_multiple(4);
            break;
        case 4:
            set_output_multiple(2);
            break;
        case 6:
            set_output_multiple(4);
            break;
        }
    }
}

/*
 * Our virtual destructor.
 */
symbol_mapper_bc_impl::~symbol_mapper_bc_impl() {}

int symbol_mapper_bc_impl::fixed_rate_ninput_to_noutput(int ninput)
{
    if (not d_is_packed) {
        return ninput / ((int)d_mapper->constellationOrder());
    } else {
        return 8 * ninput / ((int)d_mapper->constellationOrder());
    }
}

int symbol_mapper_bc_impl::fixed_rate_noutput_to_ninput(int noutput)
{
    if (not d_is_packed) {
        return noutput * ((int)d_mapper->constellationOrder());
    } else {
        return noutput * ((int)d_mapper->constellationOrder()) / 8;
    }
}

void symbol_mapper_bc_impl::forecast(int noutput_items,
                                     gr_vector_int& ninput_items_required)
{
    ninput_items_required[0] = std::max(1, fixed_rate_noutput_to_ninput(noutput_items));
}

void symbol_mapper_bc_impl::update_unpacked_vector_length(const unsigned max_length)
{
    if (d_unpacked_bits.size() < max_length) {
        d_unpacked_bits.resize(max_length);
    }
}

int symbol_mapper_bc_impl::general_work(int noutput_items,
                                        gr_vector_int& ninput_items,
                                        gr_vector_const_void_star& input_items,
                                        gr_vector_void_star& output_items)
{
    const unsigned char* in = (const unsigned char*)input_items[0];
    gr_complex* out = (gr_complex*)output_items[0];

    const unsigned num_bits =
        std::min(fixed_rate_noutput_to_ninput(noutput_items), ninput_items[0]);
    const int produced_items = fixed_rate_ninput_to_noutput(num_bits);
    // std::cout << noutput_items << ", " << produced_items << ", " << ninput_items[0] <<
    // ", " << num_bits << std::endl;
    if (d_is_packed) {
        update_unpacked_vector_length(8 * num_bits);
        d_unpacker->unpack(d_unpacked_bits.data(), in, num_bits);
        d_mapper->map_to_constellation(out, d_unpacked_bits.data(), 8 * num_bits);
    } else {
        d_mapper->map_to_constellation(out, in, num_bits);
    }

    // Tell runtime system how many input items we consumed on
    // each input stream.
    consume_each(num_bits);

    // Tell runtime system how many output items we produced.
    return produced_items;
}

} /* namespace symbolmapping */
} /* namespace gr */
