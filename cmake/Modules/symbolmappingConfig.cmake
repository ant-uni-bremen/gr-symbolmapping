INCLUDE(FindPkgConfig)
PKG_CHECK_MODULES(PC_SYMBOLMAPPING symbolmapping)

FIND_PATH(
    SYMBOLMAPPING_INCLUDE_DIRS
    NAMES symbolmapping/api.h
    HINTS $ENV{SYMBOLMAPPING_DIR}/include
        ${PC_SYMBOLMAPPING_INCLUDEDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/include
          /usr/local/include
          /usr/include
)

FIND_LIBRARY(
    SYMBOLMAPPING_LIBRARIES
    NAMES gnuradio-symbolmapping
    HINTS $ENV{SYMBOLMAPPING_DIR}/lib
        ${PC_SYMBOLMAPPING_LIBDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/lib
          ${CMAKE_INSTALL_PREFIX}/lib64
          /usr/local/lib
          /usr/local/lib64
          /usr/lib
          /usr/lib64
          )

          include("${CMAKE_CURRENT_LIST_DIR}/gnuradio-symbolmappingTargets.cmake")

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(symbolmapping DEFAULT_MSG SYMBOLMAPPING_LIBRARIES SYMBOLMAPPING_INCLUDE_DIRS)
MARK_AS_ADVANCED(SYMBOLMAPPING_LIBRARIES SYMBOLMAPPING_INCLUDE_DIRS)
