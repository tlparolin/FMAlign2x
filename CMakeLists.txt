cmake_minimum_required(VERSION 3.10)
cmake_policy(SET CMP0144 NEW)
cmake_policy(SET CMP0167 NEW)
project(FMAlign2)

# General configuration
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

# Debug or Release options
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_definitions(DEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
endif()

# Add M64, if necessary
option(M64 "Compile for 64-bit architecture" OFF)
if(M64)
    add_compile_definitions(M64)
endif()

# OpenMP, pthreads and other flags 
if(WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lpthread -pthread -lrt")
endif()

# Include directories
include_directories(include ext/SW)

# Source files
set(SOURCES
    main.cpp
    src/utils.cpp
    src/mem_finder.cpp
    src/sequence_split_align.cpp
    ext/SW/ssw.cpp
    ext/SW/ssw_cpp.cpp
)

# Add thread pool and condition
if(NOT WIN32)
    list(APPEND SOURCES src/thread_pool.cpp src/thread_condition.cpp)
endif()

# Add additional sources
set(ADDITIONAL_SOURCES src/gsacak.c)

# Find HPX
find_package(HPX REQUIRED)
if(HPX_FOUND)
    message(STATUS "HPX found: ${HPX_INCLUDE_DIRS}")
else()
    message(FATAL_ERROR "HPX not found")
endif()

# Add executable
add_executable(FMAlign2 ${SOURCES} ${ADDITIONAL_SOURCES})

# Link libraries
target_link_libraries(FMAlign2 PUBLIC HPX::hpx HPX::wrap_main HPX::iostreams_component)


# Clean all
add_custom_target(clean-all
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/FMAlign2 ${CMAKE_BINARY_DIR}/*.o
    COMMENT "Limpando arquivos gerados..."
)
