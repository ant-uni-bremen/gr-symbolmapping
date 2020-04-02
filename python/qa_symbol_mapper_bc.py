#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2020 Johannes Demel.
#
# This is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
#
# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this software; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.
#

from gnuradio import gr, gr_unittest
from gnuradio import blocks
import symbolmapping_swig as symbolmapping
import numpy as np
from symbol_constellation import generate_gray_constellation
from symbol_constellation import map_to_constellation

# import os
# print('Blocked waiting for GDB attach (pid = %d)' % (os.getpid(),))
# input('Press Enter to continue: ')

class qa_symbol_mapper_bc(gr_unittest.TestCase):

    def setUp(self):
        self.tb = gr.top_block()

    def tearDown(self):
        self.tb = None

    def test_001_init(self):
        for co in (5,7,9):
            self.assertRaises((RuntimeError, TypeError), symbolmapping.symbol_mapper_bc,
                              co, "GRAY", False)

    def verify_constellation_unpacked(self, constellation_order,
                                      nbits, is_packed):
        constellation, bits_rep = generate_gray_constellation(constellation_order)

        data = np.random.randint(0, 2, constellation_order * nbits).astype(np.uint8)
        ref = map_to_constellation(data, constellation)

        if is_packed:
            data = np.packbits(data)

        # print(f'test: Order={constellation_order}, bits={nbits}/{data.size}, packed={is_packed}')

        mapper = symbolmapping.symbol_mapper_bc(constellation_order,
                                                "GRAY", is_packed)
        src = blocks.vector_source_b(data)
        snk = blocks.vector_sink_c()

        self.tb.connect(src, mapper, snk)
        self.tb.run()

        res = np.array(snk.data())
        self.assertComplexTuplesAlmostEqual(tuple(res), tuple(ref), 6)

    def test_002_unpacked_bpsk(self):
        nbits = 16 * 700
        constellation_order = 1
        self.verify_constellation_unpacked(constellation_order,
                                           nbits, False)

    def test_003_unpacked_qpsk(self):
        nbits = 16 * 700
        constellation_order = 2
        self.verify_constellation_unpacked(constellation_order,
                                           nbits, False)

    def test_004_unpacked_8psk(self):
        nbits = 16 * 700
        constellation_order = 3
        self.verify_constellation_unpacked(constellation_order,
                                           nbits, False)

    def test_005_unpacked_16qam(self):
        nbits = 16 * 700
        constellation_order = 4
        self.verify_constellation_unpacked(constellation_order,
                                           nbits, False)

    def test_006_unpacked_64qam(self):
        nbits = 16 * 700
        constellation_order = 6
        self.verify_constellation_unpacked(constellation_order,
                                           nbits, False)

    def test_007_packed_bpsk(self):
        nbits = 16 * 700
        constellation_order = 1
        self.verify_constellation_unpacked(constellation_order,
                                           nbits, True)

    def test_008_packed_qpsk(self):
        nbits = 16 * 700
        constellation_order = 2
        self.verify_constellation_unpacked(constellation_order,
                                           nbits, True)

    # def test_009_packed_8psk(self):
    #     nbits = 16 * 700
    #     constellation_order = 3
    #     self.verify_constellation_unpacked(constellation_order,
    #                                        nbits, True)

    def test_010_packed_16qam(self):
        nbits = 16 * 700
        constellation_order = 4
        self.verify_constellation_unpacked(constellation_order,
                                           nbits, True)

    # def test_011_packed_64qam(self):
    #     nbits = 16 * 700
    #     constellation_order = 6
    #     self.verify_constellation_unpacked(constellation_order,
    #                                        nbits, True)


if __name__ == '__main__':
    gr_unittest.run(qa_symbol_mapper_bc)
