/* -*- c++ -*- */

#define SYMBOLMAPPING_API

%include "gnuradio.i"           // the common stuff

//load generated python docstrings
%include "symbolmapping_swig_doc.i"

%{
#include "symbolmapping/interleaver.h"
#include "symbolmapping/symbol_mapper_bc.h"
#include "symbolmapping/symbol_demapper_cf.h"
%}

%include "symbolmapping/interleaver.h"
GR_SWIG_BLOCK_MAGIC2_TMPL(symbolmapping, interleaver_bb, interleaver<std::uint8_t>);
GR_SWIG_BLOCK_MAGIC2_TMPL(symbolmapping, interleaver_ff, interleaver<float>);
%include "symbolmapping/symbol_mapper_bc.h"
GR_SWIG_BLOCK_MAGIC2(symbolmapping, symbol_mapper_bc);
%include "symbolmapping/symbol_demapper_cf.h"
GR_SWIG_BLOCK_MAGIC2(symbolmapping, symbol_demapper_cf);
