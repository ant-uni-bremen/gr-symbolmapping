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
import pmt
import numpy as np
import symbolmapping_swig as symbolmapping
from ref_constellation import generate_gray_constellation
from ref_constellation import map_to_constellation
from symbol_demapper import calculate_symbol_log_probabilities
from symbol_demapper import calculate_llrs

class qa_symbol_demapper_cf(gr_unittest.TestCase):
    def setUp(self):
        self.tb = gr.top_block()

    def tearDown(self):
        self.tb = None

    def test_001_scalar_snr_qpsk(self):
        nbits = 16 * 100
        constellation_order = 2
        self.verify_scalar_snr(constellation_order, nbits)

    def test_002_scalar_snr_16qam(self):
        nbits = 16 * 100
        constellation_order = 4
        self.verify_scalar_snr(constellation_order, nbits)
    
    def test_003_scalar_snr_64qam(self):
        nbits = 16 * 100
        constellation_order = 6
        self.verify_scalar_snr(constellation_order, nbits)

    def verify_scalar_snr(self, constellation_order, nbits):
        constellation, bits_rep = generate_gray_constellation(constellation_order)
        data = np.random.randint(0, 2, constellation_order * nbits).astype(np.uint8)
        symbols = map_to_constellation(data, constellation)
        
        snr = 5.0
        snr_step = .5
        tags = []
        offsets = (0, 50, 400, 450, 800)
        last_offset = 0
        ref = np.array([], dtype=np.float32)
        for offset in offsets:
            t = gr.tag_utils.python_to_tag((offset, pmt.string_to_symbol("snr"), 
                                            pmt.from_float(snr)))
            tags.append(t)
            
            syms = symbols[last_offset:offset]
            # print(syms.size)
            ref_ln_probs = calculate_symbol_log_probabilities(syms,
                                                          constellation,
                                                          snr - snr_step)
            refs = calculate_llrs(ref_ln_probs)
            ref = np.concatenate((ref, refs))
            last_offset = offset
            snr += snr_step
            
        syms = symbols[last_offset:]
        ref_ln_probs = calculate_symbol_log_probabilities(syms,
                                                          constellation,
                                                          snr - snr_step)
        refs = calculate_llrs(ref_ln_probs)
        ref = np.concatenate((ref, refs))
        ref *= .5

        # print(f'test: Order={constellation_order}, bits={nbits}/{data.size}, packed={is_packed}')
        
        mapper = symbolmapping.symbol_demapper_cf(constellation_order, 
                                                  "GRAY")
        src = blocks.vector_source_c(symbols, False, 1, tags)
        snk = blocks.vector_sink_f()

        self.tb.connect(src, mapper, snk)
        self.tb.run()

        res = np.array(snk.data())
        if constellation_order > 3:
            self.assertFloatTuplesAlmostEqual(tuple(np.sign(res)), tuple(np.sign(ref)))
        else:
            # for i in range(ref.size - 10):
            #     if not np.abs(ref[i + 10] - res[i + 10]) < 1e-5:
            #         print(i, ref[i], res[i], np.abs(ref[i] - res[i]) < 1e-5)
            self.assertFloatTuplesAlmostEqual(tuple(res), tuple(ref), 5)

    def test_004_vector_snr_qpsk(self):
        nbits = 16 * 100
        constellation_order = 2
        self.verify_vector_snr(constellation_order, nbits)

    def test_005_vector_snr_16qam(self):
        nbits = 16 * 100
        constellation_order = 4
        self.verify_vector_snr(constellation_order, nbits)

    def test_006_vector_snr_64qam(self):
        nbits = 16 * 100
        constellation_order = 6
        self.verify_vector_snr(constellation_order, nbits)

    def verify_vector_snr(self, constellation_order, nbits):
        constellation, bits_rep = generate_gray_constellation(constellation_order)
        data = np.random.randint(0, 2, constellation_order * nbits).astype(np.uint8)
        symbols = map_to_constellation(data, constellation)
        
        snr = np.arange(32, dtype=np.float) + 1.
        snr_step = 2.5
        tags = []
        offsets = (0, 50, 400, 450, 800)
        last_offset = 0
        ref = np.array([], dtype=np.float32)
        for offset in offsets:
            valuevector = pmt.init_f32vector(snr.size, snr.tolist())
            t = gr.tag_utils.python_to_tag((offset, pmt.string_to_symbol("snr"), 
                                            valuevector))
            tags.append(t)
            
            syms = symbols[last_offset:offset]
            # print(syms.size)
            ref_ln_probs = calculate_symbol_log_probabilities(syms,
                                                          constellation,
                                                          snr - snr_step)
            refs = calculate_llrs(ref_ln_probs)
            ref = np.concatenate((ref, refs))
            last_offset = offset
            snr += snr_step
            
        syms = symbols[last_offset:]
        ref_ln_probs = calculate_symbol_log_probabilities(syms,
                                                          constellation,
                                                          snr - snr_step)
        refs = calculate_llrs(ref_ln_probs)
        ref = np.concatenate((ref, refs))
        ref *= .5

        # print(f'test: Order={constellation_order}, bits={nbits}/{data.size}, packed={is_packed}')
        
        mapper = symbolmapping.symbol_demapper_cf(constellation_order, 
                                                  "GRAY")
        src = blocks.vector_source_c(symbols, False, 1, tags)
        snk = blocks.vector_sink_f()

        self.tb.connect(src, mapper, snk)
        self.tb.run()

        res = np.array(snk.data())
        if constellation_order > 3:
            self.assertFloatTuplesAlmostEqual(tuple(np.sign(res)), tuple(np.sign(ref)))
        else:
            # for i in range(ref.size - 10):
            #     if not np.abs(ref[i + 10] - res[i + 10]) < 1e-5:
            #         print(i, ref[i], res[i], np.abs(ref[i] - res[i]) < 1e-5)
            self.assertFloatTuplesAlmostEqual(tuple(res), tuple(ref), 5)


if __name__ == '__main__':
    gr_unittest.run(qa_symbol_demapper_cf)
