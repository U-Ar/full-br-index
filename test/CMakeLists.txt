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

ADD_EXECUTABLE(itest1 itest1.cpp)
TARGET_LINK_LIBRARIES(itest1 PRIVATE brindex sdsl divsufsort divsufsort64)
ADD_EXECUTABLE(itest2 itest1.cpp)
TARGET_LINK_LIBRARIES(itest2 PRIVATE brindex sdsl divsufsort divsufsort64)

# copy dataset & testing Python script to build/test/dataset
FILE(COPY dataset DESTINATION ${CMAKE_BINARY_DIR}/test)
FILE(COPY integration-test DESTINATION ${CMAKE_BINARY_DIR}/test)


