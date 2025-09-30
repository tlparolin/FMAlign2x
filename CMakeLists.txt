cmake_minimum_required(VERSION 3.10)
project(FMAlign2)

set(CMAKE_CXX_STANDARD 23)
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
    list(APPEND SOURCES src/thread_pool.cpp)
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

# WFA2 Config
CPMAddPackage(
  NAME wfa2
  GITHUB_REPOSITORY smarco/WFA2-lib
  GIT_TAG v2.3.5
  OPTIONS
    "OPENMP TRUE"
)

if(wfa2_ADDED)
    # INTERFACE to simplify include/link
    add_library(wfa2lib INTERFACE)
    target_include_directories(wfa2lib INTERFACE ${wfa2_SOURCE_DIR}/bindings/cpp)
    target_link_libraries(wfa2lib INTERFACE wfa2cpp)
endif()

find_package(OpenMP REQUIRED)

add_executable(FMAlign2x ${SOURCES})
if(OpenMP_CXX_FOUND)
    target_link_libraries(FMAlign2x PRIVATE OpenMP::OpenMP_CXX)
endif()
target_link_libraries(FMAlign2x PRIVATE libsais wfa2lib)

# Clean Rules
add_custom_target(clean-all
    COMMAND ${CMAKE_COMMAND} -E rm -rf FMAlign2x
    COMMAND ${CMAKE_COMMAND} -E rm -rf CMakeCache.txt CMakeFiles
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
