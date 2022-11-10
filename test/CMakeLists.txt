# ---- unit tests ----
ADD_EXECUTABLE(run_tests
	run_tests.cpp 
	br_index_test.cpp
	sparse_sd_vector_test.cpp 
	huffman_string_test.cpp
	rle_string_test.cpp
	permuted_lcp_test.cpp
)
TARGET_LINK_LIBRARIES(run_tests PRIVATE brindex sdsl divsufsort divsufsort64)

TARGET_INCLUDE_DIRECTORIES(run_tests PRIVATE ${PROJECT_SOURCE_DIR}/external/iutest/include)


# ---- integration tests ----

# copy dataset to build/test/dataset
FILE(COPY dataset DESTINATION ${CMAKE_BINARY_DIR}/test)

