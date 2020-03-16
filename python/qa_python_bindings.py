#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2020 Johannes Demel.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

from gnuradio import gr, gr_unittest
import symbolmapping_python as symbolmapping
import numpy as np

from ref_constellation import generate_gray_constellation
from ref_constellation import generate_16qam_boronka_constellation
from ref_constellation import generate_16qam_carson_constellation
from ref_constellation import map_to_constellation
from symbol_demapper import calculate_symbol_log_probabilities
from symbol_demapper import calculate_llrs
from symbol_demapper import lin2db, db2lin


class qa_interleaver(gr_unittest.TestCase):
    def setUp(self):
        self._nbits = 128 - 8
        self._il_indices = np.random.permutation(self._nbits)
        self._il_indices = self._il_indices.astype(np.uintp)
        self._dil_indices = np.argsort(self._il_indices)
        self._uut = symbolmapping.Interleaver(self._il_indices)

    def tearDown(self):
        pass

    def test_001_initialize(self):
        ii = self._uut.interleaverIndices()
        di = self._uut.deinterleaverIndices()
        self.assertEqual(self._uut.interleaverLength(), self._nbits)
        self.assertTrue(np.all(self._il_indices == ii))
        self.assertTrue(np.all(self._dil_indices == di))

    def test_002_interleaver(self):
        bits = np.random.randint(0, 2, self._nbits).astype(np.uint8)
        ibits = self._uut.interleave(bits)
        self.assertEqual(ibits.dtype, np.uint8)
        pibits = bits[self._il_indices]
        self.assertTrue(np.all(ibits == pibits))
        self.assertTrue(np.all(self._uut.interleave(bits) == pibits))

        pfsyms = (1. - 2. * bits).astype(np.float32)
        isyms = self._uut.interleave(pfsyms)
        self.assertEqual(isyms.dtype, np.float32)
        pisyms = pfsyms[self._il_indices]
        self.assertTrue(np.all(isyms == pisyms))
        self.assertTrue(np.all(self._uut.interleave(pfsyms) == pisyms))

    def test_003_deinterleaver(self):
        bits = np.random.randint(0, 2, self._nbits).astype(np.uint8)
        pfsyms = (1. - 2. * bits).astype(np.float32)

        ibits = self._uut.deinterleave(bits)
        self.assertEqual(ibits.dtype, np.uint8)
        pibits = bits[self._dil_indices]
        self.assertTrue(np.all(ibits == pibits))
        self.assertTrue(np.all(self._uut.deinterleave(bits) == pibits))

        psyms = pfsyms[self._dil_indices]
        syms = self._uut.deinterleave(pfsyms)
        self.assertEqual(syms.dtype, np.float32)
        self.assertTrue(np.all(syms == psyms))
        self.assertTrue(np.all(self._uut.deinterleave(pfsyms) == syms))


class qa_constellation(gr_unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_001_generate_gray_valid(self):
        valid_constellation_orders = np.array([1, 2, 4, 6])
        for co in valid_constellation_orders:
            constellation, bits_rep = generate_gray_constellation(co)
            self.assertEqual(constellation.size, 2 ** co)
            self.assertEqual(constellation.size, len(bits_rep))

    def test_002_generate_gray_invalid(self):
        valid_constellation_orders = np.array([5, 7, 8])
        for co in valid_constellation_orders:
            self.assertRaises(NotImplementedError,
                              generate_gray_constellation, co)


class qa_SymbolMapping(gr_unittest.TestCase):
    def setUp(self):
        self._orders = np.array([1, 2, 3, 4, 6])
        self._precision = 1e-6

    def tearDown(self):
        pass

    def assertVectorAlmostEqual(self, a, b):
        self.assertTrue(np.all(np.abs(a - b) < self._precision))

    def test_001_constellation_ref(self):
        orders = np.array([1, 2, 4, 6])
        for co in self._orders:
            dm = symbolmapping.SymbolMapping(co)
            c, b = generate_gray_constellation(co)
            self.assertVectorAlmostEqual(c, dm.constellation())
            self.assertEqual(dm.constellationOrder(), co)

    def test_002_constellation_notimplemented(self):
        noi_orders = np.arange(12, dtype=int)
        for o in noi_orders:
            if o not in self._orders:
                self.assertRaises(ValueError, 
                                  symbolmapping.SymbolMapping, 
                                  o)

    def test_003_constellation_set(self):
        dm = symbolmapping.SymbolMapping(1)
        orders = np.array([1, 2, 4, 6])
        for co in orders:
            dm.setConstellationOrder(co)
            c, b = generate_gray_constellation(co)
            self.assertVectorAlmostEqual(c, dm.constellation())

    def test_004_mapping(self):
        for co in self._orders:
            dm = symbolmapping.SymbolMapping(co)
            bits = np.random.randint(0, 2, co * 500).astype(np.uint8)

            ref_syms = map_to_constellation(bits, dm.constellation())
            uut_syms = dm.map_to_constellation(bits)
            self.assertVectorAlmostEqual(uut_syms, ref_syms)

    def test_005_ln_prob_calculation(self):
        snr_db = float(3.0)
        for co in self._orders:
            dm = symbolmapping.SymbolMapping(co)
            bits = np.random.randint(0, 2, co * 500).astype(np.uint8)
            symbols = dm.map_to_constellation(bits)
            ln_probs = dm.calculate_ln_probabilities(symbols, snr_db)
            snr_lin = symbolmapping.db2lin(snr_db)
            c128cstl = dm.constellation().astype(np.complex128)
            ref_ln_probs = calculate_symbol_log_probabilities(symbols,
                                                            dm.constellation(),
                                                            snr_db).flatten()
            self.assertFloatTuplesAlmostEqual(ref_ln_probs * .5, ln_probs, 4)

    def test_006_llr_calculation(self):
        for co in (2, 4, 6):
            dm = symbolmapping.SymbolMapping(co)
            bits = np.random.randint(0, 2, co * 5).astype(np.uint8)
            symbols = dm.map_to_constellation(bits)

            ref_ln_probs = calculate_symbol_log_probabilities(symbols,
                                                            dm.constellation(),
                                                            float(3.0))
            ref_llrs = calculate_llrs(ref_ln_probs) * .5

            llrs_d = dm.demap_llrs(symbols, float(3.0))
            if co > 3:
                self.assertTrue(np.all(np.sign(ref_llrs) == np.sign(llrs_d)))
            else:
                self.assertTrue(np.all(np.abs(ref_llrs - llrs_d) < self._precision))

    def test_007_constellation_type(self):
        dm = symbolmapping.SymbolMapping(4, "BORonka")
        self.assertEqual(dm.constellationType(), "BORONKA")
        dm = symbolmapping.SymbolMapping(4, "carson")
        self.assertEqual(dm.constellationType(), "CARSON")
        dm = symbolmapping.SymbolMapping(2, "carson")
        self.assertEqual(dm.constellationType(), "GRAY")

    def test_008_constellation_boronka(self):
        dm = symbolmapping.SymbolMapping(4, "BORonka")
        self.assertEqual(dm.constellationType(), "BORONKA")
        res = dm.constellation()
        ref, bins = generate_16qam_boronka_constellation()
        self.assertTrue(np.all(np.abs(res - ref) < 1e-7))

    def test_009_constellation_carson(self):
        dm = symbolmapping.SymbolMapping(4, "carSon")
        self.assertEqual(dm.constellationType(), "CARSON")
        res = dm.constellation()
        ref, bins = generate_16qam_carson_constellation()
        self.assertTrue(np.all(np.abs(res - ref) < 1e-7))

    def test_010_db2lin_conversion(self):
        snrs = np.arange(0.0, 10., .5, dtype=np.float32)
        for s in snrs:
            ref = db2lin(s)
            sls = symbolmapping.db2lin(s)
            self.assertTrue(np.abs(ref - sls) < 1e-5)

    def test_011_lin2db_conversion(self):
        snrs = np.arange(0.0, 10., .5, dtype=np.float32)
        snrs = db2lin(snrs)
        for s in snrs:
            ref = lin2db(s)
            sls = symbolmapping.lin2db(s)
            self.assertTrue(np.abs(ref - sls) < 1e-5)

if __name__ == '__main__':
    gr_unittest.run(qa_interleaver)
    gr_unittest.run(qa_constellation)
    gr_unittest.run(qa_SymbolMapping)
