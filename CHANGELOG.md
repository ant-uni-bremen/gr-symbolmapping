# Release v1.0.0 2022-02-07

This is the first major release for `gr-symbolmapping`, a GNU Radio Out-Of-Tree (OOT) module.
The release targets GNU Radio 3.9 but CI indicates that GNU Radio 3.10 is supported as well.

The project includes a generic parameterizable interleaver and de-interleaver.
This enables arbitrary interleaver sequences.

The core of the project is a collection of standardized symbol constellations with optimized soft-demappers.
These soft-demappers employ well known approximations to boost throughput.
Currently, BPSK, QPSK, 16QAM, 64QAM, and 256QAM constellations are supported.
