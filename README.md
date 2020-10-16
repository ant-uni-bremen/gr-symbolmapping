# GNU Radio Symbolmapping


In digital communication, we rely on symbol mapping such as QPSK, 16QAM, etc.

There exists a set of common definitions for these mappings in across a wide variety of standards (WiFi, LTE, 5G NR, ...).
Unfortunately, GNU Radio does not implement this standard set.

Further, we want fast and efficient demappers that employ efficient approximations for these common mappings.

## Interleavers

In most systems, interleavers are employed right before symbol mapping in order to leverage full diversity in frequency-selective fading scenarios. Especially, multicarrier systems such as OFDM, GFDM or FBMC benefit from this approach.
Thus, we add bit interleavers as well. Sometimes this interleaver is called channel interleaver.

## Python interface

Finally, we need those functions in GNU Radio but also with a Python3 interface for use in simulations. Thus, we add a PyBind11 interface that offers just that.

## Prerequisites
We provide a `Dockerfile` to compile the Python interface only. This is only for test purposes. In case you want to compile and use the Python interface only, you need several packages.

### Debian/Ubuntu packages
For the python interface
- build-essential
- git
- cmake
- libvolk2-dev
- liborc-0.4-dev
- python3-dev
- python3-distutils
- python3-pybind11

### GNU Radio installation
In case you want to use the whole system with GNU Radio and Python interface, GNU Radio should already provide everything.

GNU Radio 3.9 is required to use this module!

## Usage

This is the symbolmapping-write-a-block package meant as a guide to building
out-of-tree packages. To use the symbolmapping blocks, the Python namespaces
is in 'symbolmapping', which is imported as:

    import symbolmapping

See the Doxygen documentation for details about the blocks available
in this package. A quick listing of the details can be found in Python
after importing by using:

    help(symbolmapping)


## Contributions

Pull Requests are highly welcome to add more standardized mappings aka constellations aka alphabets.
