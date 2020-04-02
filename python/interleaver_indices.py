#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2020 Johannes Demel.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

import numpy as np


'''
TS 36.212 Table 5.1.4-1
Turbo inter-column permutation pattern
< 0, 16, 8, 24, 4, 20, 12, 28, 2, 18, 10, 26, 6, 22, 14, 30, 1, 17, 9, 25, 5, 21, 13, 29, 3, 19, 11, 27, 7, 23, 15, 31 >

TS 36.212 Table 5.1.4-2
Convolutional inter-column permutation pattern
< 1, 17, 9, 25, 5, 21, 13, 29, 3, 19, 11, 27, 7, 23, 15, 31, 0, 16, 8, 24, 4, 20, 12, 28, 2, 18, 10, 26, 6, 22, 14, 30 >

TS 38.212 Section 5.4.1.1 describes Polar code interleaving + puncturing.
Since this is dependent on puncturing, we omit it here for now.

TS 38.212 Section 5.4.2.1 describes LDPC code interleaving + bit-selection
Since this is dependent on code-structure, we omit it here for now.

'''


permutation_pattern = {
    'turbo': np.array([0, 16, 8, 24, 4, 20, 12, 28, 2, 18, 10, 26, 6, 22, 14, 30, 1, 17, 9, 25, 5, 21, 13, 29, 3, 19, 11, 27, 7, 23, 15, 31], dtype=np.int32),
    'convolutional' : np.array([1, 17, 9, 25, 5, 21, 13, 29, 3, 19, 11, 27, 7, 23, 15, 31, 0, 16, 8, 24, 4, 20, 12, 28, 2, 18, 10, 26, 6, 22, 14, 30], dtype=np.int32)
}

def create_lte_interleaver_indices(frame_size, intl_type='convolutional'):
    column_len = 32
    row_len = int(np.ceil(frame_size / column_len))
    dummy_len = column_len * row_len - frame_size
    indices = np.concatenate((np.full(dummy_len, np.NaN), np.arange(frame_size)))
    indices = np.reshape(indices, (-1, column_len))
    indices = indices[:, permutation_pattern[intl_type]]
    indices = indices.T.flatten()
    indices = indices[np.logical_not(np.isnan(indices))]
    return indices.astype(np.uintp)


def create_interleaver_indices(frame_size, intl_type='convolutional'):
    assert frame_size % 8 == 0
    intl_type = intl_type.lower()
    if intl_type in 'random':
        return np.random.permutation(frame_size)
    assert intl_type in list(permutation_pattern.keys())
    return create_lte_interleaver_indices(frame_size, intl_type)
