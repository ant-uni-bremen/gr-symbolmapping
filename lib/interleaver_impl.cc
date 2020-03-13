/* -*- c++ -*- */
/*
 * Copyright 2020 Johannes Demel.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include "interleaver_impl.h"
#include <type_traits>

namespace gr {
  namespace symbolmapping {

    template <class T>
    typename interleaver<T>::sptr
    interleaver<T>::make(std::vector<size_t> interleaver_indices, 
                         bool is_packed, bool interleave_mode)
    {
      return gnuradio::get_initial_sptr
        (new interleaver_impl<T>(interleaver_indices, is_packed, interleave_mode));
    }


    /*
     * The private constructor
     */
    template<class T>
    interleaver_impl<T>::interleaver_impl(std::vector<size_t> interleaver_indices, bool is_packed, bool interleave_mode)
      : gr::sync_block("interleaver",
              gr::io_signature::make(1, 1, sizeof(T)),
              gr::io_signature::make(1, 1, sizeof(T))),
        d_is_packed(is_packed),
        d_interleave_mode(interleave_mode),
        d_unpacker(new gr::blocks::kernel::unpack_k_bits(8)),
        d_packer(new gr::blocks::kernel::pack_k_bits(8)),
        d_interleaver(new BitInterleaver(interleaver_indices))
    {
      if(not std::is_same<T, uint8_t>::value){
        is_packed = false;
        d_is_packed = false;
      }
      if(is_packed and d_interleaver->interleaverLength() % 8){
        throw std::invalid_argument("Packed Interleaver requires 'interleaver_indices' to be a multiple of 8!");
      }
      d_unpacked_original.resize(d_interleaver->interleaverLength());
      d_unpacked_interleaved.resize(d_interleaver->interleaverLength());

      gr::sync_block::set_output_multiple(is_packed ? (d_interleaver->interleaverLength() / 8) : d_interleaver->interleaverLength());
    }

    /*
     * Our virtual destructor.
     */
    template<class T>
    interleaver_impl<T>::~interleaver_impl()
    {
    }

    template<class T>
    void 
    interleaver_impl<T>::interleave_packed(T *out, const T *in, const unsigned nbytes_per_frame){
      throw std::logic_error("Wrong type! packed interleaver is only available for uint8_t!");
    }

    template<>
    void 
    interleaver_impl<uint8_t>::interleave_packed(uint8_t *out, const uint8_t *in, const unsigned nbytes_per_frame)
    {
      d_unpacker->unpack(d_unpacked_original.data(), in, nbytes_per_frame);
      d_interleaver->interleave(d_unpacked_interleaved.data(),
                                d_unpacked_original.data());
      d_packer->pack(out, d_unpacked_interleaved.data(), nbytes_per_frame);
    }

    template<class T>
    void 
    interleaver_impl<T>::deinterleave_packed(T *out, const T *in, const unsigned nbytes_per_frame){
      throw std::logic_error("Wrong type! packed interleaver is only available for uint8_t!");
    }

    template<>
    void 
    interleaver_impl<uint8_t>::deinterleave_packed(uint8_t *out, const uint8_t *in, const unsigned nbytes_per_frame)
    {
      d_unpacker->unpack(d_unpacked_original.data(), in, nbytes_per_frame);
      d_interleaver->deinterleave(d_unpacked_interleaved.data(),
                                d_unpacked_original.data());
      d_packer->pack(out, d_unpacked_interleaved.data(), nbytes_per_frame);
    }

    template<class T>
    int
    interleaver_impl<T>::work(int noutput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items)
    {
      const T *in = (const T *) input_items[0];
      T *out = (T *) output_items[0];

      const unsigned nframes = (d_is_packed ? 8 : 1) * noutput_items / d_interleaver->interleaverLength();
      const unsigned nbytes_per_frame = d_interleaver->interleaverLength() / (d_is_packed ? 8 : 1);

      if(d_interleave_mode){
        if(d_is_packed){
          for(int i = 0; i < nframes; ++i){
            interleave_packed(out, in, nbytes_per_frame);

            in += nbytes_per_frame;
            out += nbytes_per_frame;
          }
        }
        else{
          for(int i = 0; i < nframes; ++i){
            d_interleaver->interleave(out, in);
            in += nbytes_per_frame;
            out += nbytes_per_frame;
          }
        }
      }
      else{
        if(d_is_packed){
          for(int i = 0; i < nframes; ++i){
            deinterleave_packed(out, in, nbytes_per_frame);

            in += nbytes_per_frame;
            out += nbytes_per_frame;
          }
        }
        else{
          for(int i = 0; i < nframes; ++i){
            d_interleaver->deinterleave(out, in);
            in += nbytes_per_frame;
            out += nbytes_per_frame;
          }
        }
      }

      noutput_items = nframes * nbytes_per_frame;

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

    template class interleaver<u_int8_t>;
    template class interleaver<float>;
  } /* namespace symbolmapping */
} /* namespace gr */

