/* -*- c++ -*- */
/*
 * Copyright 2020 Johannes Demel.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_SYMBOLMAPPING_INTERLEAVER_IMPL_H
#define INCLUDED_SYMBOLMAPPING_INTERLEAVER_IMPL_H

#include <symbolmapping/interleaver.h>
#include <symbolmapping/bit_interleaver.h>

#include <vector>
#include <memory>
#include <gnuradio/blocks/pack_k_bits.h>
#include <gnuradio/blocks/unpack_k_bits.h>

namespace gr {
  namespace symbolmapping {

    template <class T>
    class interleaver_impl : public interleaver<T>
    {
     private:
      std::unique_ptr<gr::blocks::kernel::unpack_k_bits> d_unpacker;
      std::unique_ptr<gr::blocks::kernel::pack_k_bits> d_packer;
      std::unique_ptr<BitInterleaver> d_interleaver;

      bool d_is_packed;
      bool d_interleave_mode;  // True==interleave, False==deinterleave
      std::vector<uint8_t> d_unpacked_original;
      std::vector<uint8_t> d_unpacked_interleaved;
      
      void interleave_packed(T *out, const T *in, const unsigned nbytes_per_frame);
      void deinterleave_packed(T *out, const T *in, const unsigned nbytes_per_frame);

     public:
      interleaver_impl(std::vector<size_t> interleaver_indices, 
                       bool is_packed, bool interleave_mode);
      ~interleaver_impl();

      // Where all the action really happens
      int work(
              int noutput_items,
              gr_vector_const_void_star &input_items,
              gr_vector_void_star &output_items
      );
    };

  } // namespace symbolmapping
} // namespace gr

#endif /* INCLUDED_SYMBOLMAPPING_INTERLEAVER_IMPL_H */

