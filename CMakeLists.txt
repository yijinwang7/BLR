cmake_minimum_required(VERSION 3.1)
project(BLRLibrary)

set(CMAKE_CXX_STANDARD 11)

# Find OpenMP
find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Use pkg-config to find fplll
find_package(PkgConfig)
pkg_check_modules(FPLLL REQUIRED fplll)

# Include directories for fplll (using variables from pkg_check_modules)
include_directories(${FPLLL_INCLUDE_DIRS})

# Link directories for fplll (using variables from pkg_check_modules)
link_directories(${FPLLL_LIBRARY_DIRS})

## Include Eigen headers (adjust path if necessary)
include_directories("/usr/local/include")

#include_directories("/home/oker/install/include/eigen3")

#add_executable(test_rand_perm bin/test_rand_perm.cpp randomPermutation.cpp utils.cpp)

# Add the executable
add_executable(c3
        src/blockLLLHy_standard.cpp
        src/blockLLLHy_standard.h
        src/blockLLLHy_BKZ.cpp
        src/blockLLLHy_BKZ.h
        src/randomPermutation.cpp
        src/randomPermutation.h
        src/utils.cpp
        include/BLR.h
        src/BLR.cpp
        test/test.cpp
        #test/test.cpp
)

# add_executable(fplll_test bin/fplll_test.cpp)

# Link OpenMP if it was found
if(OPENMP_FOUND)
    target_link_libraries(c3 PUBLIC OpenMP::OpenMP_CXX)
endif()

# Link fplll if it was found
if(FPLLL_FOUND)
    target_link_libraries(c3 PUBLIC ${FPLLL_LIBRARIES} gmp
            mpfr)
else()
    message(SEND_ERROR "fplll not found")
endif()

target_include_directories(c3 PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_include_directories(c3 PUBLIC ${PROJECT_SOURCE_DIR}/include)

#using ${PROJECT_SOURCE_DIR} in both cases ensures that your include directories are
#correctly and consistently located from the root of your project structure.
