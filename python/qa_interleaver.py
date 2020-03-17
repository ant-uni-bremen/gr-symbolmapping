#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2020 Johannes Demel.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

from gnuradio import gr, gr_unittest
from gnuradio import blocks
import symbolmapping_swig as symbolmapping
import numpy as np

class qa_interleaver(gr_unittest.TestCase):

    def setUp(self):
        self.tb = gr.top_block()

    def tearDown(self):
        self.tb = None

    def test_001_init(self):
        nbits = 16 * 10
        interleaver_indices = np.random.permutation(nbits)
        self.assertRaises(Exception, symbolmapping.interleaver_bb,
                          [interleaver_indices[0:13], True, True])

        interleaver = symbolmapping.interleaver_bb(interleaver_indices, True, True)
        self.assertEqual(interleaver.output_multiple(), nbits // 8)

        interleaveru = symbolmapping.interleaver_bb(interleaver_indices, False, True)
        self.assertEqual(interleaveru.output_multiple(), nbits)

    def test_002_unpacked(self):
        nframes = 5
        nbits = 16 * 7
        interleaver_indices = np.random.permutation(nbits)

        interleaver = symbolmapping.interleaver_bb(interleaver_indices, False, True)

        data = np.random.randint(0, 2, nbits * nframes)
        ref = np.array([], dtype=data.dtype)
        for f in np.reshape(data, (nframes, -1)):
            ref = np.concatenate((ref, f[interleaver_indices]))
        src = blocks.vector_source_b(data)
        snk = blocks.vector_sink_b()

        self.tb.connect(src, interleaver, snk)
        # set up fg
        self.tb.run()
        # # check data
        res = np.array(snk.data())
        self.assertTupleEqual(tuple(res), tuple(ref))

    def test_003_packed(self):
        nframes = 20
        nbits = 16 * 17
        interleaver_indices = np.random.permutation(nbits)

        interleaver = symbolmapping.interleaver_bb(interleaver_indices, True, True)

        data = np.random.randint(0, 2, nbits * nframes)
        ref = np.array([], dtype=data.dtype)
        for f in np.reshape(data, (nframes, -1)):
            ref = np.concatenate((ref, f[interleaver_indices]))
        datavec = np.packbits(data)
        refvec = np.packbits(ref)
        src = blocks.vector_source_b(datavec)
        snk = blocks.vector_sink_b()
        self.tb.connect(src, interleaver, snk)
        # set up fg
        self.tb.run()
        # check data
        resvec = np.array(snk.data())

        self.assertTupleEqual(tuple(resvec), tuple(refvec))

    def test_004_unpacked_float(self):
        nframes = 5
        nbits = 16 * 7
        interleaver_indices = np.random.permutation(nbits)

        interleaver = symbolmapping.interleaver_ff(interleaver_indices, False, True)

        data = np.random.normal(size=nbits * nframes).astype(np.float32)
        # data = np.random.randint(0, 2, nbits * nframes)
        ref = np.array([], dtype=data.dtype)
        for f in np.reshape(data, (nframes, -1)):
            ref = np.concatenate((ref, f[interleaver_indices]))
        src = blocks.vector_source_f(data)
        snk = blocks.vector_sink_f()

        self.tb.connect(src, interleaver, snk)
        # set up fg
        self.tb.run()
        # # check data
        res = np.array(snk.data())
        self.assertTupleEqual(tuple(res), tuple(ref))

    def test_005_deinterleave_unpacked_float(self):
        nframes = 5
        nbits = 16 * 17
        interleaver_indices = np.random.permutation(nbits)

        interleaver = symbolmapping.interleaver_ff(interleaver_indices, False, False)

        data = np.random.normal(size=nbits * nframes).astype(np.float32)
        # data = np.random.randint(0, 2, nbits * nframes)
        ref = np.array([], dtype=data.dtype)
        for f in np.reshape(data, (nframes, -1)):
            ref = np.concatenate((ref, f[interleaver_indices]))
        src = blocks.vector_source_f(ref)
        snk = blocks.vector_sink_f()

        self.tb.connect(src, interleaver, snk)
        # set up fg
        self.tb.run()
        # # check data
        res = np.array(snk.data())
        self.assertTupleEqual(tuple(res), tuple(data))


if __name__ == '__main__':
    gr_unittest.run(qa_interleaver)
