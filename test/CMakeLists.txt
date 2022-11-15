# ---- unit tests ----
ADD_EXECUTABLE(unit-test
	run_tests.cpp 
	br_index_test.cpp
	sparse_sd_vector_test.cpp 
	huffman_string_test.cpp
	rle_string_test.cpp
	permuted_lcp_test.cpp
)
TARGET_LINK_LIBRARIES(unit-test PRIVATE brindex sdsl divsufsort divsufsort64)

TARGET_INCLUDE_DIRECTORIES(unit-test PRIVATE ${PROJECT_SOURCE_DIR}/external/iutest/include)


# ---- integration tests ----

ADD_EXECUTABLE(itest1 itest1.cpp)
TARGET_LINK_LIBRARIES(itest1 PRIVATE brindex sdsl divsufsort divsufsort64)
ADD_EXECUTABLE(itest2 itest2.cpp)
TARGET_LINK_LIBRARIES(itest2 PRIVATE brindex sdsl divsufsort divsufsort64)

# copy dataset & testing Python script to build/test/dataset
FILE(COPY dataset DESTINATION ${CMAKE_BINARY_DIR}/test)
CONFIGURE_FILE(integration-test ${CMAKE_BINARY_DIR}/test/integration-test COPYONLY)
