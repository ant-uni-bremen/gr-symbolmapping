########################################################################
# Python Bindings
########################################################################

pybind11_add_module(symbolmapping_python
    ${CMAKE_SOURCE_DIR}/lib/bit_interleaver.cc
    ${CMAKE_SOURCE_DIR}/lib/symbol_mapper.cc
    bit_interleaver_python.cc
    symbol_mapper_python.cc
    python_bindings.cc
)

target_include_directories(symbolmapping_python PUBLIC
    ${PYTHON_NUMPY_INCLUDE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../../lib
    ${CMAKE_CURRENT_SOURCE_DIR}/../../../include
    ${PYBIND11_INCLUDE_DIR}
)
target_link_libraries(symbolmapping_python PUBLIC ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} gnuradio-runtime gnuradio-symbolmapping)

if(CMAKE_COMPILER_IS_GNUCC AND NOT APPLE)
    SET_TARGET_PROPERTIES(symbolmapping_python
    PROPERTIES LINK_FLAGS "-Wl,--no-as-needed")
endif()

install(TARGETS symbolmapping_python DESTINATION ${GR_PYTHON_DIR}/symbolmapping COMPONENT pythonapi)
