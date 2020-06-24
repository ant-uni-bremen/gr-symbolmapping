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

#ifndef INCLUDED_SYMBOLMAPPING_SYMBOL_MAPPER_BC_H
#define INCLUDED_SYMBOLMAPPING_SYMBOL_MAPPER_BC_H

#include <gnuradio/block.h>
#include <symbolmapping/api.h>

namespace gr {
namespace symbolmapping {

/*!
 * \brief Map packed/unpacked bits to complex symbols
 * \ingroup symbolmapping
 *
 */
class SYMBOLMAPPING_API symbol_mapper_bc : virtual public gr::block
{
public:
    typedef std::shared_ptr<symbol_mapper_bc> sptr;

    /*!
     * \brief Return a shared_ptr to a new instance of symbolmapping::symbol_mapper_bc.
     *
     * To avoid accidental use of raw pointers, symbolmapping::symbol_mapper_bc's
     * constructor is in a private implementation
     * class. symbolmapping::symbol_mapper_bc::make is the public interface for
     * creating new instances.
     */
    static sptr
    make(unsigned constellation_order, std::string constellation_type, bool is_packed);
};

} // namespace symbolmapping
} // namespace gr

#endif /* INCLUDED_SYMBOLMAPPING_SYMBOL_MAPPER_BC_H */
