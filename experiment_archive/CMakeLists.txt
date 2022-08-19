CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(experiment_archive CXX)

FIND_PACKAGE(Git QUIET)

SET(SDSL_INCLUDE "~/sdsl/include") #SDSL headeres
SET(SDSL_LIB "~/sdsl/lib") #SDSL lib

INCLUDE_DIRECTORIES(${SDSL_INCLUDE}) 
LINK_DIRECTORIES(${SDSL_LIB}) 

SET(CMAKE_CXX_STANDARD 11)

SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g")
SET(CMAKE_CXX_FLAGS_RELEASE "-g -ggdb -Ofast -fstrict-aliasing -DNDEBUG -march=native")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -ggdb -Ofast -fstrict-aliasing -march=native")

INCLUDE_DIRECTORIES(src)

ADD_EXECUTABLE(bri-build src/bri-build.cpp)
TARGET_LINK_LIBRARIES(bri-build sdsl)
TARGET_LINK_LIBRARIES(bri-build divsufsort)
TARGET_LINK_LIBRARIES(bri-build divsufsort64)

ADD_EXECUTABLE(bri-count src/bri-count.cpp)
TARGET_LINK_LIBRARIES(bri-count sdsl)
TARGET_LINK_LIBRARIES(bri-count divsufsort)
TARGET_LINK_LIBRARIES(bri-count divsufsort64)

ADD_EXECUTABLE(bri-locate src/bri-locate.cpp)
TARGET_LINK_LIBRARIES(bri-locate sdsl)
TARGET_LINK_LIBRARIES(bri-locate divsufsort)
TARGET_LINK_LIBRARIES(bri-locate divsufsort64)

ADD_EXECUTABLE(bri-mem src/bri-mem.cpp)
TARGET_LINK_LIBRARIES(bri-mem sdsl)
TARGET_LINK_LIBRARIES(bri-mem divsufsort)
TARGET_LINK_LIBRARIES(bri-mem divsufsort64)

ADD_EXECUTABLE(bri-space src/bri-space.cpp)
TARGET_LINK_LIBRARIES(bri-space sdsl)
TARGET_LINK_LIBRARIES(bri-space divsufsort)
TARGET_LINK_LIBRARIES(bri-space divsufsort64)

ADD_EXECUTABLE(fixedtest src/fixedtest.cpp)
TARGET_LINK_LIBRARIES(fixedtest sdsl)
TARGET_LINK_LIBRARIES(fixedtest divsufsort)
TARGET_LINK_LIBRARIES(fixedtest divsufsort64)

ADD_EXECUTABLE(k_vs_time src/k_vs_time.cpp)
TARGET_LINK_LIBRARIES(k_vs_time sdsl)
TARGET_LINK_LIBRARIES(k_vs_time divsufsort)
TARGET_LINK_LIBRARIES(k_vs_time divsufsort64)

ADD_EXECUTABLE(k_vs_space src/k_vs_space.cpp)
TARGET_LINK_LIBRARIES(k_vs_space sdsl)
TARGET_LINK_LIBRARIES(k_vs_space divsufsort)
TARGET_LINK_LIBRARIES(k_vs_space divsufsort64)

ADD_EXECUTABLE(window_slide src/window_slide.cpp)
TARGET_LINK_LIBRARIES(window_slide sdsl)
TARGET_LINK_LIBRARIES(window_slide divsufsort)
TARGET_LINK_LIBRARIES(window_slide divsufsort64)

ADD_EXECUTABLE(x_treeop src/x_treeop.cpp)
TARGET_LINK_LIBRARIES(x_treeop sdsl)
TARGET_LINK_LIBRARIES(x_treeop divsufsort)
TARGET_LINK_LIBRARIES(x_treeop divsufsort64)
