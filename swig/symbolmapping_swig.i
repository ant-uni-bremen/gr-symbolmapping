/* -*- c++ -*- */

#define SYMBOLMAPPING_API

%include "gnuradio.i"           // the common stuff

//load generated python docstrings
%include "symbolmapping_swig_doc.i"

%{
#include "symbolmapping/interleaver.h"
%}

%include "symbolmapping/interleaver.h"
GR_SWIG_BLOCK_MAGIC2_TMPL(symbolmapping, interleaver_bb, interleaver<std::uint8_t>);
GR_SWIG_BLOCK_MAGIC2_TMPL(symbolmapping, interleaver_ff, interleaver<float>);
