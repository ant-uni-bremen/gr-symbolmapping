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

#ifndef INCLUDED_SYMBOLMAPPING_SYMBOL_DEMAPPER_CF_H
#define INCLUDED_SYMBOLMAPPING_SYMBOL_DEMAPPER_CF_H

#include <gnuradio/sync_interpolator.h>
#include <symbolmapping/api.h>

namespace gr {
namespace symbolmapping {

/*!
 * \brief Compute soft bits, aka LLRs, for all bits per received symbol
 * \ingroup symbolmapping
 *
 */
class SYMBOLMAPPING_API symbol_demapper_cf : virtual public gr::sync_interpolator
{
public:
    typedef std::shared_ptr<symbol_demapper_cf> sptr;

    /*!
     * \brief Return a shared_ptr to a new instance of symbolmapping::symbol_demapper_cf.
     *
     * To avoid accidental use of raw pointers, symbolmapping::symbol_demapper_cf's
     * constructor is in a private implementation
     * class. symbolmapping::symbol_demapper_cf::make is the public interface for
     * creating new instances.
     */
    static sptr make(unsigned constellation_order,
                     std::string constellation_type,
                     std::string snr_tag_name = std::string("snr"));
};

} // namespace symbolmapping
} // namespace gr

#endif /* INCLUDED_SYMBOLMAPPING_SYMBOL_DEMAPPER_CF_H */
