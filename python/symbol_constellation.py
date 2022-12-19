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
    constellation = np.zeros(n_points, dtype=complex)
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
    constellation = np.zeros(n_points, dtype=complex)
    scale = 1.0 / np.sqrt(2.0)
    constellation[0] = complex(  1.0 * scale,  1.0 * scale)
    constellation[1] = complex(  1.0        ,  0.0        )
    constellation[2] = complex( -1.0        ,  0.0        )
    constellation[3] = complex( -1.0 * scale, -1.0 * scale)
    constellation[4] = complex(  0.0        ,  1.0        )
    constellation[5] = complex(  1.0 * scale, -1.0 * scale)
    constellation[6] = complex( -1.0 * scale,  1.0 * scale)
    constellation[7] = complex(  0.0 * scale, -1.0        )
    return constellation, bit_bins

def generate_16qam_gray_constellation():
    constellation_order = 4
    n_points = 2 ** constellation_order
    bit_bins = get_bit_bins(constellation_order)
    constellation = np.zeros(n_points, dtype=complex)
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
    constellation = np.zeros(n_points, dtype=complex)
    q_lut = np.array([3, 1, 5, 7])
    for i in range(n_points):
        b = list(bit_bins[i])
        c = [int(j) for j in b]
        re = (1 - 2 * c[0]) * q_lut[2 * c[2] + c[4]]
        im = (1 - 2 * c[1]) * q_lut[2 * c[3] + c[5]]
        constellation[i] = re + 1j * im
    constellation /= np.sqrt(42.)
    return constellation, bit_bins


def generate_256qam_gray_constellation():
    constellation_order = 8
    n_points = 2 ** constellation_order
    bit_bins = get_bit_bins(constellation_order)
    constellation = np.zeros(n_points, dtype=complex)
    for i in range(n_points):
        b = list(bit_bins[i])
        c = [int(j) for j in b]
        cz = 1. - 2. * np.array(c)
        re = cz[0] * (8 - cz[2] * (4 - cz[4] * (2 - cz[6])))
        im = cz[1] * (8 - cz[3] * (4 - cz[5] * (2 - cz[7])))
        constellation[i] = re + 1j * im
    constellation /= np.sqrt(170.)
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
    elif constellation_order == 8:
        return generate_256qam_gray_constellation()
    else:
        raise NotImplementedError()


def generate_16qam_boronka_constellation():
    constellation_order = 4
    bit_bins = get_bit_bins(constellation_order)
    constellation = np.array([ 1.+3.j,  3.-3.j, -1.-3.j,  1.+1.j,
                              -3.+1.j,  3.+1.j,  1.-1.j, -3.-1.j,
                               3.+3.j, -1.-1.j, -1.+1.j, -3.-3.j,
                               1.-3.j, -1.+3.j, -3.+3.j,  3.-1.j],
                             dtype=complex)
    constellation /= np.sqrt(10.)
    assert constellation.size == 2 ** constellation_order
    return constellation, bit_bins


def generate_16qam_carson_constellation():
    constellation_order = 4
    bit_bins = get_bit_bins(constellation_order)
    constellation = np.array([ 3.+3.j, -3.-3.j, -1.+3.j,  1.-3.j,
                              -3.+1.j,  3.-1.j,  1.+1.j, -1.-1.j,
                               1.-1.j, -1.+1.j, -3.-1.j,  3.+1.j,
                              -1.-3.j,  1.+3.j,  3.-3.j, -3.+3.j],
                             dtype=complex)
    constellation /= np.sqrt(10.)
    assert constellation.size == 2 ** constellation_order
    return constellation, bit_bins


def generate_constellation(constellation_order, constellation_type='GRAY'):
    constellation_type = constellation_type.upper()
    if constellation_order != 4:
        assert 'GRAY' in constellation_type
    else:
        constellation_type = constellation_type
        cstl_types = ['GRAY', 'CARSON', 'BORONKA']
        t = [i for i in cstl_types if constellation_type in i]
        if len(t) == 1:
            constellation_type = t[0]
        else:
            raise ValueError(f'Constellation type: "{constellation_type}" is not supported!')

    if 'GRAY' in constellation_type:
        return generate_gray_constellation(constellation_order)
    elif 'BORONKA' in constellation_type:
        return generate_16qam_boronka_constellation()
    elif 'CARSON' in constellation_type:
        return generate_16qam_carson_constellation()
    else:
        raise ValueError(f'Constellation type: "{constellation_type}" is not supported!')



def pack_bits(bits, n_packed):
    b = np.reshape(bits, (-1, n_packed))
    vals = 2 ** np.arange(n_packed - 1, -1, -1)
    return np.sum(b * vals, axis=1)


def map_to_constellation(bits, constellation):
    constellation_order = int(np.log2(len(constellation)))
    points = pack_bits(bits, constellation_order)
    symbols = np.array([constellation[i] for i in points])
    return symbols


def main():
    import matplotlib.pyplot as plt
    from latex_plot_magic import set_size
    latex_textwidth = 327.20668  # pt
    fig_size = set_size(latex_textwidth)
    fig = plt.figure(figsize= (fig_size[0], fig_size[0] * 3. / 4.))

    markersize = 40

    cstl, bits = generate_constellation(4)
    plt.scatter(cstl.real, cstl.imag, marker='o', s=markersize, label='16QAM')
    for c, b in zip(cstl, bits):
        pos = (c.real, c.imag)
        plt.annotate(b, pos, xytext=(-10, 5), textcoords='offset points',)

    cstl, bits = generate_constellation(2)
    plt.scatter(cstl.real, cstl.imag, marker='x', s=markersize, label='QPSK')
    for c, b in zip(cstl, bits):
        pos = (c.real, c.imag)
        plt.annotate(b, pos, xytext=(-5, 5), textcoords='offset points')

    # cstl, bits = generate_constellation(6)
    # plt.scatter(cstl.real, cstl.imag, marker='o', s=markersize, label='64QAM')
    # for c, b in zip(cstl, bits):
    #     pos = (c.real, c.imag)
    #     plt.annotate(b, pos, xytext=(-10, 5), textcoords='offset points',)

    plt.grid()
    plt.legend(fontsize='small')
    plt.xlabel('Inphase')
    plt.ylabel('Quadrature')
    # plt.ylim((-1., 1.11))
    ticks = plt.yticks()
    plt.yticks((-1., -.5, 0.0, .5, 1.))
    plt.tight_layout()
    # plt.savefig('constellation_qpsk_16qam.pgf')
    plt.show()


if __name__ == '__main__':
    main()
