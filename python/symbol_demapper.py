#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2020 Johannes Demel.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

import numpy as np
from symbol_constellation import map_to_constellation
import symbol_constellation as cstl


'''
LLR definition!
LLR = ln( P(b=0) / P(b=1) )
'''


def db2lin(snr_db):
    return 10 ** (snr_db / 10.)


def lin2db(snr_lin):
    return 10. * np.log10(snr_lin)


def get_symbol_prob(rx, constellation):
    p = np.zeros(len(constellation))
    for i, s in enumerate(constellation):
        p[i] = calculate_gaussian_probability(rx, s, 1.)
    return p


def calculate_log_probability_vector_lin(rx, constellation, snr_lin):
    '''
    Gaussian noise assumption!
    (1. / np.sqrt(2 * np.pi * sigma_sq)) * np.exp(-np.abs(rx - constellation)**2 / (2 * sigma_sq))

    Prune constant factors and move to LOG domain -->
    -.5 * snr_lin * np.abs(rx - constellation) ** 2
    '''
    pp = rx - constellation
    r = pp.real ** 2 + pp.imag ** 2
    return -1. * snr_lin * r


def calculate_log_probability_vector(rx, constellation, snr_db):
    snr_lin = db2lin(snr_db)
    return calculate_log_probability_vector_lin(rx, constellation, snr_lin)


def calculate_gaussian_probability(x, mu, variance):
    return (1. / np.sqrt(2 * np.pi * variance)) * np.exp(-.5 * (np.abs(x - mu) ** 2) / variance)


def calculate_qpsk_probs_to_llrs(probs):
    # compare qpsk bit map
    llr0 = np.sum(probs[2:4]) / np.sum(probs[0:2])
    llr1 = np.sum(probs[np.array((1, 3))]) / np.sum(probs[np.array((0, 2))])
    return llr0, llr1


def calculate_qpsk_log_probs_to_llrs(log_probs):
    # compare qpsk bit map
    llr0 = np.max(log_probs[2:4]) - np.max(log_probs[0:2])
    llr1 = np.max(log_probs[np.array((1, 3))]) - \
        np.max(log_probs[np.array((0, 2))])
    return llr0, llr1


def calculate_qpsk_llrs_messages(llrs, log_probs):
    # compare qpsk bit map
    res_llrs = np.zeros(2)
    res_llrs[0] = np.maximum(log_probs[2], log_probs[3] + llrs[1]) - \
        np.maximum(log_probs[0], log_probs[1] + llrs[1])
    res_llrs[1] = np.maximum(log_probs[1], log_probs[3] + llrs[0]) - \
        np.maximum(log_probs[0], log_probs[2] + llrs[0])
    return res_llrs * -1.


def calculate_qpsk_llrs(log_probs, llrs=None):
    n_symbols, constellation_size = np.shape(log_probs)
    constellation_order = int(np.log2(constellation_size))
    nc = n_symbols * constellation_order
    if llrs is None:
        llrs = np.zeros(nc, dtype=float)
    for i, l_prob in enumerate(log_probs):
        # llrs[i * constellation_order: (i + 1) * constellation_order] = calculate_qpsk_log_probs_to_llrs(l_prob)
        llrs[i * constellation_order: (i + 1) * constellation_order] = calculate_qpsk_llrs_messages(
            llrs[i * constellation_order: (i + 1) * constellation_order], l_prob)
    return llrs


def decide_bits(llrs):
    b = np.array(llrs) < 0.0
    return b.astype(int)


def qpsk_map_demap_chain():
    constellation_order = 2
    snr_db = 20
    constellation, bits_rep = cstl.generate_gray_constellation(
        constellation_order)
    print(constellation)
    print(bits_rep)
    app_llrs = np.zeros(constellation_order)
    for i in range(2):
        for j in range(2):
            bits = np.array((i, j))
            print(bits)
            s = map_to_constellation(bits, constellation)
            print(s)
            l_prob = calculate_log_probability_vector(s, constellation, snr_db)
            print(l_prob)
            llr0, llr1 = calculate_qpsk_llrs_messages(app_llrs, l_prob)
            print(llr0, llr1)
            hat_bits = decide_bits((llr0, llr1))
            print(hat_bits)
            print(np.all(bits == hat_bits))
            if not np.all(bits == hat_bits):
                raise RuntimeError('QPSK Mapping --> Demapping fails!')


def calculate_symbol_log_probabilities(symbols, constellation, snr_db):
    if isinstance(snr_db, np.ndarray):
        snr = np.tile(snr_db, int(
            np.ceil(symbols.size / snr_db.size)))[0:symbols.size]
    else:
        snr = np.full(symbols.size, db2lin(snr_db), dtype=np.float)
    log_probs = np.zeros((len(symbols), len(constellation)), dtype=float)
    for i, s in enumerate(symbols):
        l_prob = calculate_log_probability_vector_lin(s, constellation, snr[i])
        log_probs[i, :] = l_prob
    return log_probs


def unpack_bits(values, n_packed):
    n_vals = len(values)
    n_shift = np.arange(n_packed - 1, -1, -1)
    values = np.repeat(values, n_packed)
    values = np.right_shift(values, np.tile(n_shift, n_vals))
    values = np.bitwise_and(values, 1)
    return values


def calculate_app_llr_sum(llrs, n_bits, fix_pos):
    bits = unpack_bits(np.arange(2 ** n_bits, dtype=int), n_bits)
    bits = np.reshape(bits, (-1, n_bits))
    bits[:, fix_pos] = 0
    # print(bits)
    t = llrs.dot(bits.T)
    # print(t)
    return t


def calculate_16qam_llrs_messages(app_llrs, log_probs):
    llrs = np.zeros(4)
    for i in range(len(llrs)):
        app_sums = calculate_app_llr_sum(app_llrs, 4, i)
        # print(app_sums)
        l = log_probs + app_sums
        indices = np.arange(16)
        indices = np.reshape(indices, (2 ** (i + 1), -1))
        idx0 = indices[0::2].flatten()
        idx1 = indices[1::2].flatten()
        llrs[i] = np.max(l[idx1]) - np.max(l[idx0])
    return llrs * -1.


def calculate_16qam_llrs(log_probs, llrs=None):
    n_symbols, constellation_size = np.shape(log_probs)
    constellation_order = int(np.log2(constellation_size))
    nc = n_symbols * constellation_order
    # print(nc, n_symbols, constellation_order)
    if llrs is None:
        llrs = np.zeros(nc, dtype=float)
    for i, l_prob in enumerate(log_probs):
        # llrs[i * constellation_order: (i + 1) * constellation_order] = calculate_qpsk_log_probs_to_llrs(l_prob)
        llrs[i * constellation_order: (i + 1) * constellation_order] = calculate_16qam_llrs_messages(
            llrs[i * constellation_order: (i + 1) * constellation_order], l_prob)
    return llrs


def calculate_llrs(log_probs, llrs=None):
    n_symbols, constellation_size = np.shape(log_probs)
    constellation_order = int(np.log2(constellation_size))
    nc = n_symbols * constellation_order
    if llrs is None:
        llrs = np.zeros(nc, dtype=float)
    if constellation_order == 2:
        return calculate_qpsk_llrs(log_probs, llrs)
    elif constellation_order == 4:
        return calculate_16qam_llrs(log_probs, llrs)
    elif constellation_order == 6:
        return calculate_64qam_llrs(log_probs, llrs)
    else:
        raise NotImplementedError(
            'Constellation {:} NOT IMPLEMENTED!'.format(constellation_order))


def qam16_map_demap_chain():
    constellation_order = 4
    snr_db = 20
    constellation, bits_rep = cstl.generate_gray_constellation(
        constellation_order)
    print(constellation)
    print(bits_rep)
    app_llrs = np.zeros(constellation_order)
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for l in range(2):
                    bits = np.array((i, j, k, l))
                    # print(utils.pack_bits(bits, 4), bits)
                    app_llrs = 100. * (2. * bits - 1)
                    print(app_llrs)
                    s = map_to_constellation(bits, constellation)
                    print(s)
                    l_prob = calculate_log_probability_vector(
                        s, constellation, snr_db)
                    print(l_prob)
                    llrs = calculate_16qam_llrs_messages(app_llrs, l_prob)
                    # llrs = np.array(llrs)
                    print(llrs)
                    hat_bits = decide_bits(llrs)
                    print(hat_bits)
                    print(np.all(bits == hat_bits))
                    if not np.all(bits == hat_bits):
                        raise RuntimeError('16QAM Mapping --> Demapping fails!')


def calculate_64qam_llrs_messages(app_llrs, log_probs):
    llrs = np.zeros(6)
    for i in range(len(llrs)):
        app_l = np.sum(app_llrs[np.where(np.arange(len(app_llrs)) != i)])
        indices = np.arange(64)
        indices = np.reshape(indices, (2 ** (i + 1), -1))
        idx0 = indices[0::2].flatten()
        idx1 = indices[1::2].flatten()
        llrs[i] = np.max(log_probs[idx1] + app_l) - \
            np.max(log_probs[idx0] + app_l)
    return llrs * -1.


def calculate_64qam_llrs(log_probs, llrs=None):
    n_symbols, constellation_size = np.shape(log_probs)
    constellation_order = int(np.log2(constellation_size))
    nc = n_symbols * constellation_order
    if llrs is None:
        llrs = np.zeros(nc, dtype=float)
    for i, l_prob in enumerate(log_probs):
        llrs[i * constellation_order: (i + 1) * constellation_order] = calculate_64qam_llrs_messages(
            llrs[i * constellation_order: (i + 1) * constellation_order], l_prob)
    return llrs


def qam64_map_demap_chain():
    constellation_order = 6
    snr_db = 20
    constellation, bits_rep = cstl.generate_gray_constellation(
        constellation_order)
    print(constellation)
    print(bits_rep)
    app_llrs = np.zeros(constellation_order)
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for l in range(2):
                    for m in range(2):
                        for n in range(2):
                            bits = np.array((i, j, k, l, m, n))
                            print(bits)
                            s = map_to_constellation(bits, constellation)
                            print(s)
                            l_prob = calculate_log_probability_vector(
                                s, constellation, snr_db)
                            print(l_prob)
                            llrs = calculate_64qam_llrs_messages(
                                app_llrs, l_prob)
                            llrs = np.array(llrs)
                            print(llrs)
                            hat_bits = decide_bits(llrs)
                            print(hat_bits)
                            print(np.all(bits == hat_bits))
                            if not np.all(bits == hat_bits):
                                raise RuntimeError(
                                    '64QAM Mapping --> Demapping fails!')


def main():
    np.set_printoptions(precision=2)

    qpsk_map_demap_chain()
    qam16_map_demap_chain()
    qam64_map_demap_chain()

    try:
        calculate_llrs(np.array([[1, 4, 5, 6, 7, 8, 4, 3], ]))
    except NotImplementedError:
        pass
    # return
    rx = .6 + .7j
    constellation_order = 4
    snr_db = 20
    constellation, bits_rep = cstl.generate_gray_constellation(
        constellation_order)
    l_prob = calculate_log_probability_vector(rx, constellation, snr_db)
    # print(constellation)
    # print(bits_rep)
    # print(l_prob)
    llr0, llr1 = calculate_qpsk_log_probs_to_llrs(l_prob)
    calculate_qpsk_probs_to_llrs(np.exp(l_prob))
    # print(llr0)
    # print(llr1)

    bits = np.random.randint(0, 2, constellation_order * 300)
    symbols = map_to_constellation(bits, constellation)
    # print(symbols)
    # symbols += utils.generate_complex_noise_symbols(len(symbols), snr_db)
    log_probs = calculate_symbol_log_probabilities(
        symbols, constellation, snr_db)

    calculate_16qam_llrs(log_probs)


if __name__ == '__main__':
    main()
