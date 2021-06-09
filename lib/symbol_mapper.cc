/*
 * Copyright 2020 Johannes Demel
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 */

#include <symbolmapping/symbol_mapper.h>
#include <volk/volk.h>
#include <bitset>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdexcept>


float lin2db(const float snr_lin) { return 10.0f * std::log10(snr_lin); }

float db2lin(const float snr_db) { return std::pow(10.0f, snr_db / 10.0f); }


float convert_snr_db2lin(const float snr_db) { return db2lin(snr_db); }

SymbolMapping::SymbolMapping(unsigned constellation_order, std::string cstl_type)
    : _constellation_order(constellation_order),
      _constellation_size(1 << constellation_order),
      _cstl_type(std::string("GRAY"))
{
    setConstellationType(cstl_type);
    setConstellationOrder(constellation_order);
}

// SymbolMapping::~SymbolMapping() {}

void SymbolMapping::update_ln_prob_buffer_size(const unsigned buffer_size)
{
    if (buffer_size > _ln_prob_buffer.size()) {
        _ln_prob_buffer.resize(buffer_size);
    }
}

void SymbolMapping::setConstellationType(std::string cstl_type)
{
    for (auto& c : cstl_type) {
        c = std::toupper(c);
    }
    std::string ctype("GRAY");
    if (_constellation_order == 4) {
        if (cstl_type == std::string("BORONKA")) {
            ctype = cstl_type;
        } else if (cstl_type == std::string("CARSON")) {
            ctype = cstl_type;
        }
    }
    _cstl_type = ctype;

    setConstellationOrder(_constellation_order);
}

void SymbolMapping::setConstellationOrder(unsigned constellation_order)
{
    _constellation_order = constellation_order;
    _constellation_size = 1 << constellation_order;
    update_ln_prob_buffer_size(1024 * _constellation_size);

    _constellation.resize(_constellation_size);
    switch (constellation_order) {
    case 1:
        generate_bpsk_gray_constellation();
        break;
    case 2:
        generate_qpsk_gray_constellation();
        break;
    case 3:
        generate_8psk_gray_constellation();
        break;
    case 4:
        generate_16qam_constellation();
        break;
    case 6:
        generate_64qam_gray_constellation();
        break;
    case 8:
        generate_256qam_gray_constellation();
        break;
    default:
        throw std::invalid_argument("Constellation not implemented!");
    }
}

void SymbolMapping::generate_bpsk_gray_constellation()
{
    // LTE constellation
    _constellation[0] = fcmplx(1.0f, 0.0f);
    _constellation[1] = fcmplx(-1.0f, 0.0f);
}

void SymbolMapping::generate_qpsk_gray_constellation()
{
    // LTE constellation
    float scale = 1.0f / std::sqrt(2.0f);
    _constellation[0] = fcmplx(1.0f * scale, 1.0f * scale);
    _constellation[1] = fcmplx(1.0f * scale, -1.0f * scale);
    _constellation[2] = fcmplx(-1.0f * scale, 1.0f * scale);
    _constellation[3] = fcmplx(-1.0f * scale, -1.0f * scale);
}

void SymbolMapping::generate_8psk_gray_constellation()
{
    // This Constellation is taken from the DVB-S2 standard
    // ETSI EN 302 307 V1.2.1 (2009-08), Sec. 5.4.2, p26
    // http://www.etsi.org/deliver/etsi_en/302300_302399/302307/01.02.01_60/en_302307v010201p.pdf
    float scale = 1.0f / std::sqrt(2.0f);
    _constellation[0] = fcmplx(1.0f * scale, 1.0f * scale);
    _constellation[1] = fcmplx(1.0f, 0.0f);
    _constellation[2] = fcmplx(-1.0f, 0.0f);
    _constellation[3] = fcmplx(-1.0f * scale, -1.0f * scale);
    _constellation[4] = fcmplx(0.0f, 1.0f);
    _constellation[5] = fcmplx(1.0f * scale, -1.0f * scale);
    _constellation[6] = fcmplx(-1.0f * scale, 1.0f * scale);
    _constellation[7] = fcmplx(0.0f * scale, -1.0f);
}

void SymbolMapping::generate_16qam_constellation()
{
    if (_cstl_type == std::string("BORONKA")) {
        generate_16qam_boronka_constellation();
    } else if (_cstl_type == std::string("CARSON")) {
        generate_16qam_carson_constellation();
    } else {
        generate_16qam_gray_constellation();
    }
}

void SymbolMapping::generate_16qam_gray_constellation()
{
    // LTE constellation
    float scale = 1.0f / std::sqrt(10.0f);
    for (unsigned i = 0; i < _constellation_size; ++i) {
        float b0 = 1.0f - 2.0f * float((i >> 3) & 0x01);
        float b1 = 1.0f - 2.0f * float((i >> 2) & 0x01);
        float b2 = 1.0f + 2.0f * float((i >> 1) & 0x01);
        float b3 = 1.0f + 2.0f * float(i & 0x01);
        _constellation[i] = fcmplx(b0 * b2 * scale, b1 * b3 * scale);
    }
}

void SymbolMapping::generate_16qam_boronka_constellation()
{
    float scale = 1.0f / std::sqrt(10.0f);
    _constellation[0] = scale * fcmplx(1.0f, 3.0f);
    _constellation[1] = scale * fcmplx(3.0f, -3.0f);
    _constellation[2] = scale * fcmplx(-1.0f, -3.0f);
    _constellation[3] = scale * fcmplx(1.0f, 1.0f);

    _constellation[4] = scale * fcmplx(-3.0f, 1.0f);
    _constellation[5] = scale * fcmplx(3.0f, 1.0f);
    _constellation[6] = scale * fcmplx(1.0f, -1.0f);
    _constellation[7] = scale * fcmplx(-3.0f, -1.0f);

    _constellation[8] = scale * fcmplx(3.0f, 3.0f);
    _constellation[9] = scale * fcmplx(-1.0f, -1.0f);
    _constellation[10] = scale * fcmplx(-1.0f, 1.0f);
    _constellation[11] = scale * fcmplx(-3.0f, -3.0f);

    _constellation[12] = scale * fcmplx(1.0f, -3.0f);
    _constellation[13] = scale * fcmplx(-1.0f, 3.0f);
    _constellation[14] = scale * fcmplx(-3.0f, 3.0f);
    _constellation[15] = scale * fcmplx(3.0f, -1.0f);
}

void SymbolMapping::generate_16qam_carson_constellation()
{
    float scale = 1.0f / std::sqrt(10.0f);
    _constellation[0] = scale * fcmplx(3.0f, 3.0f);
    _constellation[1] = scale * fcmplx(-3.0f, -3.0f);
    _constellation[2] = scale * fcmplx(-1.0f, 3.0f);
    _constellation[3] = scale * fcmplx(1.0f, -3.0f);

    _constellation[4] = scale * fcmplx(-3.0f, 1.0f);
    _constellation[5] = scale * fcmplx(3.0f, -1.0f);
    _constellation[6] = scale * fcmplx(1.0f, 1.0f);
    _constellation[7] = scale * fcmplx(-1.0f, -1.0f);

    _constellation[8] = scale * fcmplx(1.0f, -1.0f);
    _constellation[9] = scale * fcmplx(-1.0f, 1.0f);
    _constellation[10] = scale * fcmplx(-3.0f, -1.0f);
    _constellation[11] = scale * fcmplx(3.0f, 1.0f);

    _constellation[12] = scale * fcmplx(-1.0f, -3.0f);
    _constellation[13] = scale * fcmplx(1.0f, 3.0f);
    _constellation[14] = scale * fcmplx(3.0f, -3.0f);
    _constellation[15] = scale * fcmplx(-3.0f, 3.0f);
}

void SymbolMapping::generate_64qam_gray_constellation()
{
    // NR/LTE constellation
    float scale = 1.0f / std::sqrt(42.0f);
    float q_lut[] = { 3.0f, 1.0f, 5.0f, 7.0f };
    for (unsigned i = 0; i < _constellation_size; ++i) {
        float c0 = 1.0f - 2.0f * float((i >> 5) & 0x01);
        float c1 = 1.0f - 2.0f * float((i >> 4) & 0x01);
        int b2 = int((i >> 3) & 0x01);
        int b3 = int((i >> 2) & 0x01);
        int b4 = int((i >> 1) & 0x01);
        int b5 = int(i & 0x01);
        float c2 = q_lut[2 * b2 + b4];
        float c3 = q_lut[2 * b3 + b5];
        _constellation[i] = fcmplx(c0 * c2 * scale, c1 * c3 * scale);
    }
}

void SymbolMapping::generate_256qam_gray_constellation()
{
    // NR/LTE constellation
    float scale = 1.0f / std::sqrt(170.0f);
    // float q_lut[] = { 3.0f, 1.0f, 5.0f, 7.0f };
    for (unsigned i = 0; i < _constellation_size; ++i) {
        int b0 = int((i >> 7) & 0x01);
        int b1 = int((i >> 6) & 0x01);
        int b2 = int((i >> 5) & 0x01);
        int b3 = int((i >> 4) & 0x01);
        int b4 = int((i >> 3) & 0x01);
        int b5 = int((i >> 2) & 0x01);
        int b6 = int((i >> 1) & 0x01);
        int b7 = int(i & 0x01);

        float c0 = 1.0f - 2.0f * float(b0);
        float c1 = 1.0f - 2.0f * float(b1);
        float c2 = 1.0f - 2.0f * float(b2);
        float c3 = 1.0f - 2.0f * float(b3);
        float c4 = 1.0f - 2.0f * float(b4);
        float c5 = 1.0f - 2.0f * float(b5);
        float c6 = 1.0f - 2.0f * float(b6);
        float c7 = 1.0f - 2.0f * float(b7);

        float inphase = c0 * (8.0f - c2 * (4.0f - c4 * (2.0f - c6)));
        float quadrature = c1 * (8.0f - c3 * (4.0f - c5 * (2.0f - c7)));
        _constellation[i] = fcmplx(inphase * scale, quadrature * scale);
    }
}

void SymbolMapping::map_to_constellation(fcmplx* symbols,
                                         const unsigned char* bits,
                                         const unsigned num_bits)
{
    if (num_bits % _constellation_order) {
        std::runtime_error(
            "Number of provided bits and constellation order multiple do not match!");
    }
    const unsigned num_symbols = num_bits / _constellation_order;
    for (unsigned s = 0; s < num_symbols; ++s) {
        unsigned sym_byte = 0;
        for (unsigned b = _constellation_order; b > 0; --b) {
            sym_byte |= (*bits++ & 0x01) << (b - 1);
        }
        *symbols++ = _constellation[sym_byte];
    }
}

void SymbolMapping::calculate_ln_probabilities_vec(float* ln_probs,
                                                   const fcmplx* rx_symbols,
                                                   const float* snr_lin,
                                                   const unsigned num_symbols)
{
    switch (_constellation_size) {
    // case 2:
    //     calculate_ln_probabilities_t<2>(ln_probs, rx_symbols, snr_lin, num_symbols);
    //     break;
    // case 3:
    //     calculate_ln_probabilities_t<3>(ln_probs, rx_symbols, snr_lin, num_symbols);
    //     break;
    // case 4:
    //     calculate_ln_probabilities_t<4>(ln_probs, rx_symbols, snr_lin, num_symbols);
    //     break;
    // case 6:
    //     calculate_ln_probabilities_t<6>(ln_probs, rx_symbols, snr_lin, num_symbols);
    //     break;
    default:
        for (unsigned i = 0; i < num_symbols; i++) {
            const float scaling_factor = -0.5f * *snr_lin++;
            calculate_symbol_ln_probabilities(ln_probs, rx_symbols[i], scaling_factor);
            ln_probs += _constellation_size;
        }
    }
}

void SymbolMapping::calculate_ln_probabilities(float* ln_probs,
                                               const fcmplx* rx_symbols,
                                               const unsigned num_symbols,
                                               const float snr_db)
{
    /*
     * Assume white noise. Start with i.i.d. complex Gaussian noise:
     * p(x) = (1 / 2 * pi * sigma_n2) * exp(- norm(x - s) / 2 * sigma_n2)
     * Now cut of constant factors and transform this into the log domain:
     * L(x, s) = - norm(x - s) * SNR_lin
     * Assume E(x^2) = 1
     * SNR_lin = 1 / sigma_n2
     */

    const float snr_lin = convert_snr_db2lin(snr_db);
    switch (_constellation_size) {
    case 2:
        calculate_ln_probabilities_t<2>(ln_probs, rx_symbols, num_symbols, snr_lin);
        break;
    case 3:
        calculate_ln_probabilities_t<3>(ln_probs, rx_symbols, num_symbols, snr_lin);
        break;
    case 4:
        calculate_ln_probabilities_t<4>(ln_probs, rx_symbols, num_symbols, snr_lin);
        break;
    case 6:
        calculate_ln_probabilities_t<6>(ln_probs, rx_symbols, num_symbols, snr_lin);
        break;
    default:
        const float scaling_factor = -0.5f * snr_lin;
        for (unsigned i = 0; i < num_symbols; i++) {
            calculate_symbol_ln_probabilities(ln_probs, rx_symbols[i], scaling_factor);
            ln_probs += _constellation_size;
        }
    }
}


void SymbolMapping::calculate_symbol_ln_probabilities(float* ln_probs,
                                                      const fcmplx rx_symbol,
                                                      const float scaling_factor)
{
    // Calculate exponents in LLR expression for all constellation points!

    //    fcmplx val = rx_symbol;
    //    volk_32fc_x2_s32f_square_dist_scalar_mult_32f(ln_probs, &val, _constellation,
    //    scaling_factor, _constellation_size);

    for (const auto& c : _constellation) {
        // Calculate exponent: |y - x|^2 * SNR_lin
        *ln_probs++ = std::norm(rx_symbol - c) * scaling_factor;
    }
}

void SymbolMapping::calculate_llrs(float* llrs,
                                   const float* ln_probs,
                                   const unsigned num_symbols)
{
    for (unsigned i = 0; i < num_symbols; ++i) {
        calculate_symbol_llrs(llrs, ln_probs);
        llrs += _constellation_order;
        ln_probs += _constellation_size;
    }
}

void SymbolMapping::calculate_symbol_llrs(float* llrs, const float* ln_probs)
{
    // Loop through all bits per symbol
    int bitmask = _constellation_size;
    const float MIN_VAL = -1e16;
    for (unsigned l = 0; l < _constellation_order; ++l) {
        // Move to next lower bit
        bitmask >>= 1;

        // Determine minimum distance for bit l equal one / zero
        float dist0 = MIN_VAL; // minimum distances for b = 1 / 0
        float dist1 = MIN_VAL;

        for (unsigned k = 0; k < _constellation_size; ++k) {
            if (k & bitmask) {
                if (ln_probs[k] > dist1)
                    dist1 = ln_probs[k];
            } else {
                if (ln_probs[k] > dist0)
                    dist0 = ln_probs[k];
            }
        }

        // max-log-map approximation
        llrs[l] = dist0 - dist1;
    }
}

void SymbolMapping::calculate_llrs_apriori(float* llrs,
                                           const float* apllrs,
                                           const float* ln_probs,
                                           const unsigned num_symbols)
{
    for (unsigned i = 0; i < num_symbols; ++i) {
        calculate_symbol_llrs_apriori(llrs, apllrs, ln_probs);
        llrs += _constellation_order;
        apllrs += _constellation_order;
        ln_probs += _constellation_size;
    }
}

void SymbolMapping::calculate_symbol_llrs_apriori(float* llrs,
                                                  const float* apllrs,
                                                  const float* ln_probs)
{
    // Loop through all bits per symbol
    int bitmask = _constellation_size;
    const float MIN_VAL = -1e16;
    for (unsigned l = 0; l < _constellation_order; ++l) {
        // Move to next lower bit
        bitmask >>= 1;

        // Determine minimum distance for bit l equal one / zero
        float dist0 = MIN_VAL; // minimum distances for b = 1 / 0
        float dist1 = MIN_VAL;

        //      std::cout << "bit l=" << l << "\t\t" << apllrs[l] << std::endl;

        for (unsigned k = 0; k < _constellation_size; ++k) {
            const float app_factor = calculate_apriori_sum(apllrs, k);
            const float llr_factor = ln_probs[k] + app_factor;
            if (k & bitmask) {
                if (llr_factor > dist1)
                    dist1 = llr_factor;
            } else {
                if (llr_factor > dist0)
                    dist0 = llr_factor;
            }
            //        std::cout << k <<  "\tapp=" << app_factor << ",\tfactor=" <<
            //        llr_factor << ",\t0=" << dist0 << ",\t1=" << dist1 << std::endl;
        }

        // max-log-map approximation
        llrs[l] = (dist0 - dist1) - apllrs[l];
        //    llrs[l] = dist0 - dist1;
    }
}

float SymbolMapping::calculate_apriori_sum(const float* apllrs, const unsigned bitmask)
{
    float result = 0.0f;
    unsigned significant_bit = _constellation_size;
    for (unsigned i = 0; i < _constellation_order; ++i) {
        significant_bit >>= 1;
        if (bitmask & significant_bit) {
            result -= apllrs[i];
        }
    }
    return result;
}

void SymbolMapping::demap_llrs_vec(float* llrs,
                                   const fcmplx* rx_symbols,
                                   const float* snr_lin,
                                   const unsigned num_symbols)
{
    switch (_constellation_order) {
    case 1:
        demap_llrs_vec_bpsk(llrs, rx_symbols, snr_lin, num_symbols);
        break;
    case 2:
        demap_llrs_vec_qpsk(llrs, rx_symbols, snr_lin, num_symbols);
        break;
    case 4:
        demap_llrs_vec_16qam(llrs, rx_symbols, snr_lin, num_symbols);
        break;
    case 6:
        demap_llrs_vec_64qam(llrs, rx_symbols, snr_lin, num_symbols);
        break;
    case 8:
        demap_llrs_vec_256qam(llrs, rx_symbols, snr_lin, num_symbols);
        break;
    default:
        demap_llrs_vec_generic(llrs, rx_symbols, snr_lin, num_symbols);
    }
}

void SymbolMapping::demap_llrs(float* llrs,
                               const fcmplx* rx_symbols,
                               const unsigned num_symbols,
                               const float snr_db)
{
    const float snr_lin = convert_snr_db2lin(snr_db);
    switch (_constellation_order) {
    case 1:
        demap_llrs_bpsk(llrs, rx_symbols, num_symbols, snr_lin);
        break;
    case 2:
        demap_llrs_qpsk(llrs, rx_symbols, num_symbols, snr_lin);
        break;
    case 4:
        demap_llrs_16qam(llrs, rx_symbols, num_symbols, snr_lin);
        break;
    case 6:
        demap_llrs_64qam(llrs, rx_symbols, num_symbols, snr_lin);
        break;
    case 8:
        demap_llrs_256qam(llrs, rx_symbols, num_symbols, snr_lin);
        break;
    default:
        demap_llrs_generic(llrs, rx_symbols, num_symbols, snr_db);
    }
}

void SymbolMapping::demap_llrs_vec_generic(float* llrs,
                                           const fcmplx* rx_symbols,
                                           const float* snr_lin,
                                           const unsigned num_symbols)
{
    update_ln_prob_buffer_size(num_symbols * _constellation_size);
    calculate_ln_probabilities_vec(
        _ln_prob_buffer.data(), rx_symbols, snr_lin, num_symbols);
    calculate_llrs(llrs, _ln_prob_buffer.data(), num_symbols);
}

void SymbolMapping::demap_llrs_generic(float* llrs,
                                       const fcmplx* rx_symbols,
                                       const unsigned num_symbols,
                                       const float snr_db)
{
    update_ln_prob_buffer_size(num_symbols * _constellation_size);
    calculate_ln_probabilities(_ln_prob_buffer.data(), rx_symbols, num_symbols, snr_db);
    calculate_llrs(llrs, _ln_prob_buffer.data(), num_symbols);
}

void SymbolMapping::demap_llrs_vec_bpsk(float* llrs,
                                        const fcmplx* rx_symbols,
                                        const float* snr_lin,
                                        const unsigned num_symbols)
{
    const float scale = 2.0f;
    for (unsigned i = 0; i < num_symbols; ++i) {
        const float scaling_factor = scale * *snr_lin++;
        *llrs++ = scaling_factor * real(*rx_symbols++);
    }
}

void SymbolMapping::demap_llrs_bpsk(float* llrs,
                                    const fcmplx* rx_symbols,
                                    const unsigned num_symbols,
                                    const float snr_lin)
{
    // do not waste time with noise that does not carry information, i.e. only consider
    // real part.
    const float scaling_factor = 2.0f * snr_lin;
    for (unsigned i = 0; i < num_symbols; ++i) {
        *llrs++ = scaling_factor * real(*rx_symbols++);
    }
}

void SymbolMapping::demap_llrs_vec_qpsk(float* llrs,
                                        const fcmplx* rx_symbols,
                                        const float* snr_lin,
                                        const unsigned num_symbols)
{
    constexpr float scale = std::sqrt(2.0f);
    const float* rx_values = (const float*)rx_symbols;
    for (unsigned i = 0; i < num_symbols; ++i) {
        const float scaling_factor = scale * *snr_lin++;
        *llrs++ = scaling_factor * *rx_values++;
        *llrs++ = scaling_factor * *rx_values++;
    }
}

void SymbolMapping::demap_llrs_qpsk(float* llrs,
                                    const fcmplx* rx_symbols,
                                    const unsigned num_symbols,
                                    const float snr_lin)
{
    const float scaling_factor = std::sqrt(2.0f) * snr_lin;
    const float* rx_values = (const float*)rx_symbols;
    for (unsigned i = 0; i < 2 * num_symbols; ++i) {
        *llrs++ = scaling_factor * *rx_values++;
    }
}

void SymbolMapping::demap_llrs_8psk(float* llrs,
                                    const fcmplx* rx_symbols,
                                    const unsigned num_symbols,
                                    const float snr_lin)
{
    // Compare: Dirk Wuebben "Effiziente Detektionsverfahren fuer Multilayer-MIMO-Systeme"
    // Assume E_s=1, a=1
    //    const float eta = 4.0f / snr_lin; // 4 * a * sqrt(E_s) / N_0
    //    const float delta = 8.0f / snr_lin; // 8 * a^2 * E_s / N_0
    constexpr float decision_bound = 2.0f / std::sqrt(10.0f);
    const float scaling_factor = snr_lin * decision_bound;
    for (unsigned i = 0; i < num_symbols; ++i) {
        const fcmplx sym = *rx_symbols++;
        const float symi = real(sym);
        const float symq = imag(sym);
        *llrs++ = scaling_factor * symi;
        *llrs++ = scaling_factor * symq;
        *llrs++ = scaling_factor * (decision_bound - std::abs(symi));
        *llrs++ = scaling_factor * (decision_bound - std::abs(symq));
    }
}


void SymbolMapping::demap_llrs_vec_16qam(float* llrs,
                                         const fcmplx* rx_symbols,
                                         const float* snr_lin,
                                         const unsigned num_symbols)
{
    constexpr float decision_bound = 2.0f / std::sqrt(10.0f);
    constexpr float scale = decision_bound;
    //    const float* rx_values = (const float*) rx_symbols;
    for (unsigned i = 0; i < num_symbols; ++i) {
        const float scaling_factor = scale * *snr_lin++;
        const fcmplx sym = *rx_symbols++;
        const float symi = real(sym);
        const float symq = imag(sym);
        *llrs++ = scaling_factor * symi;
        *llrs++ = scaling_factor * symq;
        *llrs++ = scaling_factor * (decision_bound - std::abs(symi));
        *llrs++ = scaling_factor * (decision_bound - std::abs(symq));
    }
}

void SymbolMapping::demap_llrs_16qam(float* llrs,
                                     const fcmplx* rx_symbols,
                                     const unsigned num_symbols,
                                     const float snr_lin)
{
    // Compare: Dirk Wuebben "Effiziente Detektionsverfahren fuer Multilayer-MIMO-Systeme"
    // Assume E_s=1, a=1
    //    const float eta = 4.0f / snr_lin; // 4 * a * sqrt(E_s) / N_0
    //    const float delta = 8.0f / snr_lin; // 8 * a^2 * E_s / N_0
    constexpr float decision_bound = 2.0f / std::sqrt(10.0f);
    const float scaling_factor = snr_lin * decision_bound;
    //    const float* rx_values = (const float*) rx_symbols;
    for (unsigned i = 0; i < num_symbols; ++i) {
        const fcmplx sym = *rx_symbols++;
        const float symi = real(sym);
        const float symq = imag(sym);
        *llrs++ = scaling_factor * symi;
        *llrs++ = scaling_factor * symq;
        *llrs++ = scaling_factor * (decision_bound - std::abs(symi));
        *llrs++ = scaling_factor * (decision_bound - std::abs(symq));
    }
}

void SymbolMapping::demap_llrs_vec_64qam(float* llrs,
                                         const fcmplx* rx_symbols,
                                         const float* snr_lin,
                                         const unsigned num_symbols)
{
    constexpr float normalization_factor = 2.0f / std::sqrt(42.0f);
    constexpr float decision_bound1 = 2.0f * normalization_factor;
    constexpr float decision_bound2 = 1.0f * normalization_factor;

    for (unsigned i = 0; i < num_symbols; ++i) {
        const float scaling_factor = normalization_factor * *snr_lin++;
        const fcmplx sym = *rx_symbols++;
        const float symi = real(sym);
        const float symq = imag(sym);
        *llrs++ = scaling_factor * symi;
        *llrs++ = scaling_factor * symq;
        *llrs++ = scaling_factor * (decision_bound1 - std::abs(symi));
        *llrs++ = scaling_factor * (decision_bound1 - std::abs(symq));
        *llrs++ = scaling_factor *
                  (decision_bound2 - std::abs(std::abs(symi) - decision_bound1));
        *llrs++ = scaling_factor *
                  (decision_bound2 - std::abs(std::abs(symq) - decision_bound1));
    }
}

void SymbolMapping::demap_llrs_64qam(float* llrs,
                                     const fcmplx* rx_symbols,
                                     const unsigned num_symbols,
                                     const float snr_lin)
{
    constexpr float normalization_factor = 2.0f / std::sqrt(42.0f);
    constexpr float decision_bound1 = 2.0f * normalization_factor;
    constexpr float decision_bound2 = 1.0f * normalization_factor;
    const float scaling_factor = snr_lin * normalization_factor;

    for (unsigned i = 0; i < num_symbols; ++i) {
        const fcmplx sym = *rx_symbols++;
        const float symi = real(sym);
        const float symq = imag(sym);
        *llrs++ = scaling_factor * symi;
        *llrs++ = scaling_factor * symq;
        *llrs++ = scaling_factor * (decision_bound1 - std::abs(symi));
        *llrs++ = scaling_factor * (decision_bound1 - std::abs(symq));
        *llrs++ = scaling_factor *
                  (decision_bound2 - std::abs(std::abs(symi) - decision_bound1));
        *llrs++ = scaling_factor *
                  (decision_bound2 - std::abs(std::abs(symq) - decision_bound1));
    }
}

void SymbolMapping::demap_llrs_vec_256qam(float* llrs,
                                          const fcmplx* rx_symbols,
                                          const float* snr_lin,
                                          const unsigned num_symbols)
{
    constexpr float normalization_factor = 2.0f / std::sqrt(170.0f);
    constexpr float decision_bound1 = 4.0f * normalization_factor;
    constexpr float decision_bound2 = 2.0f * normalization_factor;
    constexpr float decision_bound3 = 1.0f * normalization_factor;

    for (unsigned i = 0; i < num_symbols; ++i) {
        const float scaling_factor = normalization_factor * *snr_lin++;
        const fcmplx sym = *rx_symbols++;
        const float symi = real(sym);
        const float symq = imag(sym);
        // layer 0
        *llrs++ = scaling_factor * symi;
        *llrs++ = scaling_factor * symq;
        // layer 1: bound @ 8
        *llrs++ = scaling_factor * (decision_bound1 - std::abs(symi));
        *llrs++ = scaling_factor * (decision_bound1 - std::abs(symq));
        // layer 2: bound @ 4
        *llrs++ = scaling_factor *
                  (decision_bound2 - std::abs(std::abs(symi) - decision_bound1));
        *llrs++ = scaling_factor *
                  (decision_bound2 - std::abs(std::abs(symq) - decision_bound1));
        // layer 3: bound @ 2
        *llrs++ = scaling_factor *
                  (decision_bound3 - std::abs(std::abs(std::abs(symi) - decision_bound1) -
                                              decision_bound2));
        *llrs++ = scaling_factor *
                  (decision_bound3 - std::abs(std::abs(std::abs(symq) - decision_bound1) -
                                              decision_bound2));
    }
}

void SymbolMapping::demap_llrs_256qam(float* llrs,
                                      const fcmplx* rx_symbols,
                                      const unsigned num_symbols,
                                      const float snr_lin)
{
    constexpr float normalization_factor = 2.0f / std::sqrt(170.0f);
    constexpr float decision_bound1 = 4.0f * normalization_factor;
    constexpr float decision_bound2 = 2.0f * normalization_factor;
    constexpr float decision_bound3 = 1.0f * normalization_factor;
    const float scaling_factor = snr_lin * normalization_factor;

    for (unsigned i = 0; i < num_symbols; ++i) {
        const fcmplx sym = *rx_symbols++;
        const float symi = real(sym);
        const float symq = imag(sym);
        // layer 0
        *llrs++ = scaling_factor * symi;
        *llrs++ = scaling_factor * symq;
        // layer 1: bound @ 8
        *llrs++ = scaling_factor * (decision_bound1 - std::abs(symi));
        *llrs++ = scaling_factor * (decision_bound1 - std::abs(symq));
        // layer 2: bound @ 4
        *llrs++ = scaling_factor *
                  (decision_bound2 - std::abs(std::abs(symi) - decision_bound1));
        *llrs++ = scaling_factor *
                  (decision_bound2 - std::abs(std::abs(symq) - decision_bound1));
        // layer 3: bound @ 2
        *llrs++ = scaling_factor *
                  (decision_bound3 - std::abs(std::abs(std::abs(symi) - decision_bound1) -
                                              decision_bound2));
        *llrs++ = scaling_factor *
                  (decision_bound3 - std::abs(std::abs(std::abs(symq) - decision_bound1) -
                                              decision_bound2));
    }
}
