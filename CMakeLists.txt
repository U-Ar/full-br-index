CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(full-br-index C CXX)

SET(SDSL_INCLUDE "~/sdsl/include") #SDSL headeres
SET(SDSL_LIB "~/sdsl/lib") #SDSL lib

# copy & build using external makefile for pfpbuilder
FILE(COPY "${CMAKE_SOURCE_DIR}/src/pfpbuilder" DESTINATION "${CMAKE_BINARY_DIR}")
IF(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/__nothing.h)
    message(FATAL_ERROR "File \"${CMAKE_CURRENT_BINARY_DIR}/__nothing.h\" found, \
    this should never be created, remove!")
ENDIF()
ADD_CUSTOM_TARGET(
    PFP_builder_source ALL
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/__nothing.h
)
ADD_CUSTOM_COMMAND(
    OUTPUT
        ${CMAKE_CURRENT_BINARY_DIR}/__nothing.h  # fake! ensure we run!
	COMMAND make 
	WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/pfpbuilder"
	COMMENT "Build external PFPbuilder source."
)

# copy bri-build(Python script) to binary dir
FILE(COPY "${CMAKE_SOURCE_DIR}/src/bri-build" DESTINATION "${CMAKE_BINARY_DIR}")

# include&link sdsl's location
INCLUDE_DIRECTORIES(${SDSL_INCLUDE}) 
LINK_DIRECTORIES(${SDSL_LIB}) 

SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_FLAGS_DEBUG          "-O0 -g -ggdb")
SET(CMAKE_CXX_FLAGS_RELEASE        "-Ofast -march=native -DNDEBUG")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-Ofast -march=native -g -ggdb")

SET(CMAKE_C_STANDARD 99)
SET(CMAKE_C_FLAGS_DEBUG          "-O0 -g -ggdb")
SET(CMAKE_C_FLAGS_RELEASE        "-Ofast -march=native -DNDEBUG")
SET(CMAKE_C_FLAGS_RELWITHDEBINFO "-Ofast -march=native -g -ggdb")


INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/include")

SET(BRI_SOURCES 
	src/definitions.cpp
	src/sparse_sd_vector.cpp
	src/huffman_string.cpp
	src/br_index.cpp)

ADD_LIBRARY(brindex STATIC ${BRI_SOURCES})


# ------ add executables using br_index here ------
SET(PROJECT_EXECUTABLES
	bri-mem
	bri-space)

FOREACH( exe ${PROJECT_EXECUTABLES} )
	ADD_EXECUTABLE( ${exe} src/${exe}.cpp) # add executable source file to src directory
ENDFOREACH()

FOREACH( exe ${PROJECT_EXECUTABLES} )
	TARGET_LINK_LIBRARIES( ${exe} brindex )
    TARGET_LINK_LIBRARIES( ${exe} sdsl )
    TARGET_LINK_LIBRARIES( ${exe} divsufsort )
    TARGET_LINK_LIBRARIES( ${exe} divsufsort64 )
ENDFOREACH()
# -------------------------------------------------


# pfpbuild.cpp
ADD_EXECUTABLE(pfpbuild src/pfpbuild.cpp src/pfpbuilder/utils.c)
TARGET_LINK_LIBRARIES(pfpbuild brindex)
TARGET_LINK_LIBRARIES(pfpbuild sdsl)
TARGET_LINK_LIBRARIES(pfpbuild divsufsort)
TARGET_LINK_LIBRARIES(pfpbuild divsufsort64)


# unit tests
ADD_EXECUTABLE(run_tests
	test/run_tests.cpp 
	test/br_index_test.cpp
	test/sparse_sd_vector_test.cpp 
	test/huffman_string_test.cpp
	test/rle_string_test.cpp
	test/permuted_lcp_test.cpp
)
TARGET_LINK_LIBRARIES(run_tests PRIVATE brindex sdsl divsufsort divsufsort64)

TARGET_INCLUDE_DIRECTORIES(run_tests PRIVATE ${PROJECT_SOURCE_DIR}/external/iutest/include)


ADD_CUSTOM_TARGET(unit-tests
	COMMAND rm -rf test-tmp
	COMMAND mkdir test-tmp
	COMMAND ./run_tests
	COMMAND rm -rf test-tmp
	DEPENDS run_tests
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
	COMMENT "Execute run_tests."
)
