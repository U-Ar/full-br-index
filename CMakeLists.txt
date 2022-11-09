CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(full-br-index CXX)

FIND_PACKAGE(Git QUIET)


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


SET(SDSL_INCLUDE "~/sdsl/include") #SDSL headeres
SET(SDSL_LIB "~/sdsl/lib") #SDSL lib

INCLUDE_DIRECTORIES(${SDSL_INCLUDE}) 
LINK_DIRECTORIES(${SDSL_LIB}) 

SET(CMAKE_CXX_STANDARD 11)

SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g")
SET(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG -march=native")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -ggdb -Ofast -march=native")

# ADD_SUBDIRECTORY(test)
# INCLUDE_DIRECTORIES(src)

ADD_EXECUTABLE(bri-mem src/bri-mem.cpp src/br_index.cpp)
TARGET_LINK_LIBRARIES(bri-mem sdsl)
TARGET_LINK_LIBRARIES(bri-mem divsufsort)
TARGET_LINK_LIBRARIES(bri-mem divsufsort64)

ADD_EXECUTABLE(bri-space src/bri-space.cpp src/br_index.cpp)
TARGET_LINK_LIBRARIES(bri-space sdsl)
TARGET_LINK_LIBRARIES(bri-space divsufsort)
TARGET_LINK_LIBRARIES(bri-space divsufsort64)

ADD_EXECUTABLE(pfpbuild src/pfpbuild.cpp src/br_index.cpp src/pfpbuilder/utils.c)
TARGET_LINK_LIBRARIES(pfpbuild sdsl)
TARGET_LINK_LIBRARIES(pfpbuild divsufsort)
TARGET_LINK_LIBRARIES(pfpbuild divsufsort64)



add_executable(run_tests test/run_tests.cpp 
	test/sparse_sd_vector_test.cpp 
	test/huffman_string_test.cpp
	test/rle_string_test.cpp
	test/permuted_lcp_test.cpp
)
target_link_libraries(run_tests PRIVATE sdsl divsufsort divsufsort64)

target_include_directories(run_tests PRIVATE ${PROJECT_SOURCE_DIR}/external/iutest/include)


ADD_CUSTOM_TARGET(unit-tests
	COMMAND rm -rf test-tmp
	COMMAND mkdir test-tmp
	COMMAND ./run_tests
	COMMAND rm -rf test-tmp
	DEPENDS run_tests
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
	COMMENT "Execute run_tests."
)
