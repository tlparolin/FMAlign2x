cmake_minimum_required(VERSION 3.10)
project(FMAlign2)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

option(DEBUG "Enable debugging" 0)
option(M64 "Use 64-bit" 0)

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

# Libsais Config
if(M64)
    set(LIBSAIS_SRC ext/libsais/libsais64.c)
    set(LIBSAIS_HDR ext/libsais/libsais64.h)
else()
    set(LIBSAIS_SRC ext/libsais/libsais.c)
    set(LIBSAIS_HDR ext/libsais/libsais.h)
endif()

add_library(libsais STATIC ${LIBSAIS_SRC})

add_executable(FMAlign2 ${SOURCES})
target_link_libraries(FMAlign2 PRIVATE libsais)

# Clean Rules
add_custom_target(clean_all
    COMMAND ${CMAKE_COMMAND} -E rm -rf FMAlign2
    COMMAND ${CMAKE_COMMAND} -E rm -rf CMakeCache.txt CMakeFiles
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
