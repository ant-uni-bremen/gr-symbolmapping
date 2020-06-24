/* -*- c++ -*- */
/*
 * Copyright 2020 Johannes Demel.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_SYMBOLMAPPING_INTERLEAVER_H
#define INCLUDED_SYMBOLMAPPING_INTERLEAVER_H

#include <gnuradio/sync_block.h>
#include <symbolmapping/api.h>
#include <cstdint>

namespace gr {
namespace symbolmapping {

/*!
 * \brief Fully parameterizable interleaver
 * \ingroup symbolmapping
 *
 */
template <class T>
class SYMBOLMAPPING_API interleaver : virtual public gr::sync_block
{
public:
    typedef std::shared_ptr<interleaver<T>> sptr;

    /*!
     * \brief Return a shared_ptr to a new instance of symbolmapping::interleaver.
     *
     * To avoid accidental use of raw pointers, symbolmapping::interleaver's
     * constructor is in a private implementation
     * class. symbolmapping::interleaver::make is the public interface for
     * creating new instances.
     */
    static sptr
    make(std::vector<size_t> interleaver_indices, bool is_packed, bool interleave_mode);
};

typedef interleaver<uint8_t> interleaver_bb;
typedef interleaver<float> interleaver_ff;
} // namespace symbolmapping
} // namespace gr

#endif /* INCLUDED_SYMBOLMAPPING_INTERLEAVER_H */
