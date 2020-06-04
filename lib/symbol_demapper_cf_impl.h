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

#ifndef INCLUDED_SYMBOLMAPPING_SYMBOL_DEMAPPER_CF_IMPL_H
#define INCLUDED_SYMBOLMAPPING_SYMBOL_DEMAPPER_CF_IMPL_H

#include <symbolmapping/symbol_demapper_cf.h>
#include <symbolmapping/symbol_mapper.h>
#include <volk/volk_alloc.hh>

#include <memory>

namespace gr {
namespace symbolmapping {

class symbol_demapper_cf_impl : public symbol_demapper_cf
{
private:
    std::unique_ptr<SymbolMapping> d_mapper;

    float d_snr_db;
    volk::vector<float> d_snrs_lin;

    pmt::pmt_t d_tag_key;

    void demap_llrs(float* out, const gr_complex* in, const unsigned nitems);

    void update_snr(const float snr_db);
    void update_snr(const volk::vector<float>& snrs_lin);

    void handle_tag(const tag_t& tag);


public:
    symbol_demapper_cf_impl(unsigned constellation_order,
                            std::string constellation_type,
                            std::string snr_tag_name = std::string("snr"));
    ~symbol_demapper_cf_impl();

    // Where all the action really happens
    int work(int noutput_items,
             gr_vector_const_void_star& input_items,
             gr_vector_void_star& output_items);
};

} // namespace symbolmapping
} // namespace gr

#endif /* INCLUDED_SYMBOLMAPPING_SYMBOL_DEMAPPER_CF_IMPL_H */
