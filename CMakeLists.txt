cmake_minimum_required(VERSION 3.10)
project(FMAlign2)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

option(DEBUG "Enable debugging" OFF)
option(M64 "Use 64-bit" OFF)

# Compiler config
set(CMAKE_CXX_FLAGS "-Wall -fopenmp -DLIBSAIS_OPENMP")
set(CMAKE_C_FLAGS "-fopenmp -DLIBSAIS_OPENMP")

if(DEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g -DDEBUG")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
endif()

if(M64)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -DM64")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m64 -DM64")
endif()

if(NOT WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lpthread -pthread -lrt")
endif()

include_directories(include ext/libsais)

# Main Sources
set(SOURCES
    main.cpp
    src/utils.cpp
    src/mem_finder.cpp
    src/sequence_split_align.cpp
    ext/SW/ssw.cpp
    ext/SW/ssw_cpp.cpp
)

if(NOT WIN32)
    list(APPEND SOURCES src/thread_pool.cpp src/thread_condition.cpp)
endif()

file(
  DOWNLOAD
  https://github.com/cpm-cmake/CPM.cmake/releases/download/v0.42.0/CPM.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake
  EXPECTED_HASH SHA256=2020b4fc42dba44817983e06342e682ecfc3d2f484a581f11cc5731fbe4dce8a
)
include(${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake)

# Libsais Config
CPMAddPackage(
  NAME libsais
  GITHUB_REPOSITORY IlyaGrebnov/libsais
  GIT_TAG v2.10.2
  OPTIONS
    "LIBSAIS_USE_OPENMP ON"
    "LIBSAIS_BUILD_SHARED_LIB OFF"
)

# SPOA Config
CPMAddPackage(
  NAME spoa
  GITHUB_REPOSITORY rvaser/spoa
  GIT_TAG 4.1.0
  OPTIONS
    "SPOA_BUILD_SHARED OFF"
    "SPOA_ENABLE_CPU_DISPATCH ON"
    "SPOA_OPTIMIZE_FOR_NATIVE ON"
)

# MAFFT Directories
set(MAFFT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/ext/mafft/mafft-7.526)
set(MAFFT_BINARIES ${CMAKE_SOURCE_DIR}/ext/mafft/linux/usr)

# MAFFT binary path (to check if mafft was compiled and exists)
set(MAFFT_BIN_CHECK ${MAFFT_BINARIES}/libexec/mafft/disttbfast)

# Create target for mafft
add_custom_command(
    OUTPUT ${MAFFT_BIN_CHECK}
    COMMAND make clean
    COMMAND make -j
    COMMAND make install
    COMMAND make clean
    COMMAND chmod 777 ${MAFFT_BIN_CHECK}
    WORKING_DIRECTORY ${MAFFT_SOURCE_DIR}/core
    COMMENT "Building MAFFT"
    DEPENDS ${MAFFT_SOURCE_DIR}/core/Makefile
)

add_custom_target(mafft DEPENDS ${MAFFT_BIN_CHECK})

add_executable(FMAlign2x ${SOURCES})
target_link_libraries(FMAlign2x libsais spoa)

# FMAlign2 depends on mafft
add_dependencies(FMAlign2x mafft spoa)

# Clean Rules
add_custom_target(clean-all
    COMMAND ${CMAKE_COMMAND} -E rm -rf FMAlign2x
    COMMAND ${CMAKE_COMMAND} -E rm -rf CMakeCache.txt CMakeFiles
    COMMAND ${CMAKE_COMMAND} -E rm -rf ${CMAKE_SOURCE_DIR}/ext/mafft/linux
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
