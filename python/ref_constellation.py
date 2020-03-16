#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2020 Johannes Demel.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

import numpy as np


def get_bit_bins(constellation_order):
    n_points = 2 ** constellation_order
    bit_ints = np.arange(n_points)
    return [np.binary_repr(bit_ints[i], width=constellation_order)
                for i in range(n_points)]


def generate_bpsk_gray_constellation():
    return np.array([1, -1]), ['0', '1']


def generate_qpsk_gray_constellation():
    constellation_order = 2
    n_points = 2 ** constellation_order
    bit_bins = get_bit_bins(constellation_order)
    constellation = np.zeros(n_points, dtype=np.complex)
    for i in range(n_points):
        b = list(bit_bins[i])
        c = [1 - 2 * int(j) for j in b]
        constellation[i] = c[0] + 1j * c[1]
    constellation /= np.sqrt(2.)
    return constellation, bit_bins


def generate_8psk_gray_constellation():
    constellation_order = 3
    n_points = 2 ** constellation_order
    bit_bins = get_bit_bins(constellation_order)
    constellation = np.zeros(n_points, dtype=np.complex)
    scale = 1.0 / np.sqrt(2.0)
    constellation[0] = np.complex(  1.0 * scale,  1.0 * scale)
    constellation[1] = np.complex(  1.0        ,  0.0        )
    constellation[2] = np.complex( -1.0        ,  0.0        )
    constellation[3] = np.complex( -1.0 * scale, -1.0 * scale)
    constellation[4] = np.complex(  0.0        ,  1.0        )
    constellation[5] = np.complex(  1.0 * scale, -1.0 * scale)
    constellation[6] = np.complex( -1.0 * scale,  1.0 * scale)
    constellation[7] = np.complex(  0.0 * scale, -1.0        )
    return constellation, bit_bins

def generate_16qam_gray_constellation():
    constellation_order = 4
    n_points = 2 ** constellation_order
    bit_bins = get_bit_bins(constellation_order)
    constellation = np.zeros(n_points, dtype=np.complex)
    for i in range(n_points):
        b = list(bit_bins[i])
        c = [1 - 2 * int(j) for j in b]
        re = c[0] * (1 + 2 * int(b[2]))
        im = c[1] * (1 + 2 * int(b[3]))
        constellation[i] = re + 1j * im
    constellation /= np.sqrt(10.)
    return constellation, bit_bins


def generate_64qam_gray_constellation():
    constellation_order = 6
    n_points = 2 ** constellation_order
    bit_bins = get_bit_bins(constellation_order)
    constellation = np.zeros(n_points, dtype=np.complex)
    q_lut = np.array([3, 1, 5, 7])
    for i in range(n_points):
        b = list(bit_bins[i])
        c = [int(j) for j in b]
        re = (1 - 2 * c[0]) * q_lut[2 * c[2] + c[4]]
        im = (1 - 2 * c[1]) * q_lut[2 * c[3] + c[5]]
        constellation[i] = re + 1j * im
    constellation /= np.sqrt(42.)
    return constellation, bit_bins


def generate_gray_constellation(constellation_order):
    if constellation_order == 1:
        return generate_bpsk_gray_constellation()
    elif constellation_order == 2:
        return generate_qpsk_gray_constellation()
    elif constellation_order == 3:
        return generate_8psk_gray_constellation()
    elif constellation_order == 4:
        return generate_16qam_gray_constellation()
    elif constellation_order == 6:
        return generate_64qam_gray_constellation()
    else:
        raise NotImplementedError()


def generate_16qam_boronka_constellation():
    constellation_order = 4
    n_points = 2 ** constellation_order
    bit_bins = get_bit_bins(constellation_order)
    constellation = np.zeros(n_points, dtype=np.complex)
    constellation = np.array([ 1.+3.j,  3.-3.j, -1.-3.j,  1.+1.j,
                              -3.+1.j,  3.+1.j,  1.-1.j, -3.-1.j,
                               3.+3.j, -1.-1.j, -1.+1.j, -3.-3.j,
                               1.-3.j, -1.+3.j, -3.+3.j,  3.-1.j],
                             dtype=np.complex)
    constellation /= np.sqrt(10.)
    return constellation, bit_bins


def generate_16qam_carson_constellation():
    constellation_order = 4
    n_points = 2 ** constellation_order
    bit_bins = get_bit_bins(constellation_order)
    constellation = np.zeros(n_points, dtype=np.complex)
    constellation = np.array([ 3.+3.j, -3.-3.j, -1.+3.j,  1.-3.j,
                              -3.+1.j,  3.-1.j,  1.+1.j, -1.-1.j,
                               1.-1.j, -1.+1.j, -3.-1.j,  3.+1.j,
                              -1.-3.j,  1.+3.j,  3.-3.j, -3.+3.j],
                             dtype=np.complex)
    constellation /= np.sqrt(10.)
    return constellation, bit_bins


def pack_bits(bits, n_packed):
    b = np.reshape(bits, (-1, n_packed))
    vals = 2 ** np.arange(n_packed - 1, -1, -1)
    return np.sum(b * vals, axis=1)


def map_to_constellation(bits, constellation):
    constellation_order = int(np.log2(len(constellation)))
    points = pack_bits(bits, constellation_order)
    symbols = np.array([constellation[i] for i in points])
    return symbols
