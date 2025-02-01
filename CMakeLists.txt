cmake_minimum_required(VERSION 3.10)
project(FMAlign2 C CXX)

set (CMAKE_CXX_STANDARD 20)
set (CMAKE_CXX_STANDARD_REQUIRED True)

# add seqan3 to search path
list(APPEND CMAKE_PREFIX_PATH "ext/seqan3")
# add the Sharg Parser to search path
list (APPEND CMAKE_PREFIX_PATH "ext/sharg-parser/cmake")

# require seqan3 with a version between >=3.0.0 and <4.0.0
find_package (seqan3 3.0 REQUIRED)
find_package (sharg 1.0 REQUIRED)

# Define source files
set(SRCS
    main.cpp
    src/utils.cpp
    src/mem_finder.cpp
    src/sequence_split_align.cpp
    ext/SW/ssw.cpp
    ext/SW/ssw_cpp.cpp
    src/gsacak.c
)

set_source_files_properties(src/gsacak.c PROPERTIES LANGUAGE C)

# Add platform-specific sources
if(NOT WIN32)
    list(APPEND SRCS
        src/thread_pool.cpp
        src/thread_condition.cpp
    )
endif()

# Include directories
include_directories(include ext/SW)

# Compiler flags
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O0 -g -DDEBUG")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O2 -Wno-interference-size")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")
endif()

# Platform-specific compiler flags
if(WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lpthread -pthread -lrt")
endif()

# Define executable
add_executable(FMAlign2 ${SRCS})

target_link_libraries(FMAlign2 seqan3::seqan3 sharg::sharg)

# Additional definitions
if(M64)
    add_definitions(-DM64)
endif()

# Clean target
add_custom_target(clean-all
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/FMAlign2
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/*.o
    COMMENT "Cleaning up binaries and object files"
)
