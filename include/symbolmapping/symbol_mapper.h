/* -*- c++ -*- */
/*
 * Copyright 2020 Johannes Demel.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef SYMBOL_DEMAPPER_H
#define SYMBOL_DEMAPPER_H

#include <volk/volk_alloc.hh>
#include <complex>
#include <cstddef>
#include <string>
#include <vector>
#define fcmplx std::complex<float>


float lin2db(const float snr_lin);
float db2lin(const float snr_db);
float convert_snr_db2lin(const float snr_db);


class __attribute__((visibility("default"))) SymbolMapping
{
public:
    SymbolMapping(unsigned constellation_order,
                  std::string cstl_type = std::string("GRAY"));
    SymbolMapping() = default;
    ~SymbolMapping() = default;
    void setConstellationOrder(unsigned constellation_order);
    void setConstellationType(std::string cstl_type);

    unsigned constellationOrder() { return _constellation_order; };
    unsigned constellationSize() { return _constellation_size; };
    std::vector<fcmplx> constellation()
    {
        return std::vector<fcmplx>(_constellation.data(),
                                   _constellation.data() + _constellation.size());
    };
    std::string constellationType() { return _cstl_type; };

    // Use constellation as LUT to obtain num_bits / constellationOrder complex symbols
    void map_to_constellation(fcmplx* symbols,
                              const unsigned char* bits,
                              const unsigned num_bits);

    void calculate_ln_probabilities(float* ln_probs,
                                    const fcmplx* rx_symbols,
                                    const unsigned num_symbols,
                                    const float snr_db);
    void calculate_ln_probabilities_vec(float* ln_probs,
                                        const fcmplx* rx_symbols,
                                        const float* snr_lin,
                                        const unsigned num_symbols);
    void calculate_llrs(float* llrs, const float* ln_probs, const unsigned num_symbols);

    void demap_llrs(float* llrs,
                    const fcmplx* rx_symbols,
                    const unsigned num_symbols,
                    const float snr_db);
    void demap_llrs_vec(float* llrs,
                        const fcmplx* rx_symbols,
                        const float* snr_lin,
                        const unsigned num_symbols);


    size_t ln_prob_calculation_duration_ns() { return _ln_prob_calculation_duration; };
    size_t llr_calculation_duration_ns() { return _llr_calculation_duration; };
    size_t apllr_calculation_duration_ns() { return _apllr_calculation_duration; };
    size_t demap_llrs_duration_ns() { return _demap_llrs_duration; };

    void calculate_llrs_apriori(float* llrs,
                                const float* apllrs,
                                const float* ln_probs,
                                const unsigned num_symbols);

private:
    size_t _ln_prob_calculation_duration;
    size_t _llr_calculation_duration;
    size_t _demap_llrs_duration;
    size_t _apllr_calculation_duration;

    unsigned _constellation_order;
    unsigned _constellation_size;
    std::string _cstl_type;

    volk::vector<fcmplx> _constellation;

    volk::vector<float> _ln_prob_buffer;
    void update_ln_prob_buffer_size(const unsigned buffer_size);

    void generate_bpsk_gray_constellation();
    void generate_qpsk_gray_constellation();
    void generate_8psk_gray_constellation();
    void generate_16qam_constellation();
    void generate_16qam_gray_constellation();
    void generate_16qam_boronka_constellation();
    void generate_16qam_carson_constellation();
    void generate_64qam_gray_constellation();
    void generate_256qam_gray_constellation();

    void calculate_symbol_ln_probabilities(float* ln_probs,
                                           const fcmplx rx_symbol,
                                           const float scaling_factor);
    void calculate_symbol_llrs(float* llrs, const float* ln_probs);
    void calculate_symbol_llrs_apriori(float* llrs,
                                       const float* apllrs,
                                       const float* ln_probs);
    float calculate_apriori_sum(const float* apllrs, const unsigned bitmask);

    void demap_llrs_generic(float* llrs,
                            const fcmplx* rx_symbols,
                            const unsigned num_symbols,
                            const float snr_db);
    void demap_llrs_bpsk(float* llrs,
                         const fcmplx* rx_symbols,
                         const unsigned num_symbols,
                         const float snr_db);
    void demap_llrs_qpsk(float* llrs,
                         const fcmplx* rx_symbols,
                         const unsigned num_symbols,
                         const float snr_db);
    void demap_llrs_8psk(float* llrs,
                         const fcmplx* rx_symbols,
                         const unsigned num_symbols,
                         const float snr_db);
    void demap_llrs_16qam(float* llrs,
                          const fcmplx* rx_symbols,
                          const unsigned num_symbols,
                          const float snr_db);
    void demap_llrs_64qam(float* llrs,
                          const fcmplx* rx_symbols,
                          const unsigned num_symbols,
                          const float snr_db);
    void demap_llrs_256qam(float* llrs,
                           const fcmplx* rx_symbols,
                           const unsigned num_symbols,
                           const float snr_db);

    void demap_llrs_vec_generic(float* llrs,
                                const fcmplx* rx_symbols,
                                const float* snr_lin,
                                const unsigned num_symbols);
    void demap_llrs_vec_bpsk(float* llrs,
                             const fcmplx* rx_symbols,
                             const float* snr_lin,
                             const unsigned num_symbols);
    void demap_llrs_vec_qpsk(float* llrs,
                             const fcmplx* rx_symbols,
                             const float* snr_lin,
                             const unsigned num_symbols);
    void demap_llrs_vec_16qam(float* llrs,
                              const fcmplx* rx_symbols,
                              const float* snr_lin,
                              const unsigned num_symbols);
    void demap_llrs_vec_64qam(float* llrs,
                              const fcmplx* rx_symbols,
                              const float* snr_lin,
                              const unsigned num_symbols);
    void demap_llrs_vec_256qam(float* llrs,
                               const fcmplx* rx_symbols,
                               const float* snr_lin,
                               const unsigned num_symbols);

    template <int constellation_size>
    void calculate_symbol_ln_probabilities_t(float* ln_probs,
                                             const fcmplx rx_symbol,
                                             const float scaling_factor)
    {
        // Calculate exponents in LLR expression for all constellation points!
        for (const auto& c : _constellation) {
            // Calculate exponent: |y - x|^2 * SNR_lin
            *ln_probs++ = std::norm(rx_symbol - c) * scaling_factor;
        }
    }

    template <int constellation_size>
    void calculate_ln_probabilities_t(float* ln_probs,
                                      const fcmplx* rx_symbols,
                                      const unsigned num_symbols,
                                      const float snr_lin)
    {
        const float scaling_factor = -0.5f * snr_lin;
        for (unsigned i = 0; i < num_symbols; i++) {
            calculate_symbol_ln_probabilities_t<constellation_size>(
                ln_probs, rx_symbols[i], scaling_factor);
            ln_probs += constellation_size;
        }
    }

    template <int constellation_size>
    void calculate_ln_probabilities_vec_t(float* ln_probs,
                                          const fcmplx* rx_symbols,
                                          const float* snr_lin,
                                          const unsigned num_symbols)
    {
        for (unsigned i = 0; i < num_symbols; ++i) {
            const float scaling_factor = -0.5f * *snr_lin++;
            calculate_symbol_ln_probabilities_t<constellation_size>(
                ln_probs, rx_symbols[i], scaling_factor);
            ln_probs += constellation_size;
        }
    }
};


#endif // SYMBOL_DEMAPPER
