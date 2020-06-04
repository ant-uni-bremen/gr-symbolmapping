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

#include "symbol_demapper_cf_impl.h"
#include <gnuradio/io_signature.h>

namespace gr {
namespace symbolmapping {

symbol_demapper_cf::sptr symbol_demapper_cf::make(unsigned constellation_order,
                                                  std::string constellation_type,
                                                  std::string snr_tag_name)
{
    return gnuradio::get_initial_sptr(new symbol_demapper_cf_impl(
        constellation_order, constellation_type, snr_tag_name));
}


/*
 * The private constructor
 */
symbol_demapper_cf_impl::symbol_demapper_cf_impl(unsigned constellation_order,
                                                 std::string constellation_type,
                                                 std::string snr_tag_name)
    : gr::sync_interpolator("symbol_demapper_cf",
                            gr::io_signature::make(1, 1, sizeof(gr_complex)),
                            gr::io_signature::make(1, 1, sizeof(float)),
                            constellation_order),
      d_mapper(new SymbolMapping(constellation_order, constellation_type)),
      d_tag_key(pmt::string_to_symbol(snr_tag_name))
{
    update_snr(3.0);
}

/*
 * Our virtual destructor.
 */
symbol_demapper_cf_impl::~symbol_demapper_cf_impl() {}

void symbol_demapper_cf_impl::update_snr(const float snr_db)
{
    d_snrs_lin.resize(0);
    d_snr_db = snr_db;
}

void symbol_demapper_cf_impl::update_snr(const volk::vector<float>& snrs_lin)
{
    d_snrs_lin = snrs_lin;
}

void symbol_demapper_cf_impl::demap_llrs(float* out,
                                         const gr_complex* in,
                                         const unsigned nitems)
{
    // std::cout << "demap_llrs(.., nitems= " << nitems << ")\n";
    if (d_snrs_lin.size() < 1) {
        d_mapper->demap_llrs(out, in, nitems, d_snr_db);
    } else {
        const unsigned nin_frame = d_snrs_lin.size();
        const unsigned nout_frame = d_snrs_lin.size() * d_mapper->constellationOrder();
        const unsigned nframes = nitems / nin_frame;

        // std::cout << "demap_llrs(.., nitems= "
        //           << nitems
        //           << ", snr_lin.size= "
        //           << d_snrs_lin.size()
        //           << ", nframes= "
        //           << nframes
        //           << ")\n";

        for (unsigned i = 0; i < nframes; ++i) {
            d_mapper->demap_llrs_vec(out, in, d_snrs_lin.data(), nin_frame);
            in += nin_frame;
            out += nout_frame;
        }
        d_mapper->demap_llrs_vec(
            out, in, d_snrs_lin.data(), nitems - nframes * nin_frame);
    }
}

void symbol_demapper_cf_impl::handle_tag(const tag_t& tag)
{
    // std::cout << pmt::is_vector(tag.value) << ", " << pmt::is_f32vector(tag.value) <<
    // std::endl;
    if (pmt::is_f32vector(tag.value)) {
        // std::cout << "Vector SNR" << std::endl;
        auto vec = pmt::f32vector_elements(tag.value);
        auto snrs = volk::vector<float>(vec.begin(), vec.end());
        update_snr(snrs);
        // std::cout << "Vector SNR go" << std::endl;
    } else {
        // std::cout << "update_snr: " << pmt::to_float(tag.value) << "\t@" << tag.offset
        // << std::endl;
        update_snr(pmt::to_float(tag.value));
    }
}

int symbol_demapper_cf_impl::work(int noutput_items,
                                  gr_vector_const_void_star& input_items,
                                  gr_vector_void_star& output_items)
{
    const gr_complex* in = (const gr_complex*)input_items[0];
    float* out = (float*)output_items[0];
    // const auto tag_key = pmt::string_to_symbol("snr");
    std::vector<tag_t> tags;
    const uint64_t ninput_items = noutput_items / interpolation();
    get_tags_in_range(tags, 0, nitems_read(0), nitems_read(0) + ninput_items, d_tag_key);
    // std::cout << "nread: " << nitems_read(0) <<  ", nout: " << noutput_items <<
    // std::endl;
    if (tags.size() > 0) {
        uint64_t last_offset = nitems_read(0);
        for (auto& tag : tags) {
            const unsigned nin_pre = tag.offset - last_offset;
            this->demap_llrs(out, in, nin_pre);
            in += nin_pre;
            out += nin_pre * d_mapper->constellationOrder();
            last_offset = tag.offset;
            // tag update goes here!
            handle_tag(tag);
        }
        // any remaining tags items after last tag?
        const u_int64_t remainder = nitems_read(0) + ninput_items - last_offset;
        if (remainder) {
            this->demap_llrs(out, in, remainder);
        }
    } else {
        this->demap_llrs(out, in, noutput_items / d_mapper->constellationOrder());
    }

    // Tell runtime system how many output items we produced.
    return noutput_items;
}

} /* namespace symbolmapping */
} /* namespace gr */
