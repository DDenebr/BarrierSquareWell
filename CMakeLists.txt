cmake_minimum_required(VERSION 3.7...3.21)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

# Project name
project(
    EDMDSquareWellBarrier 
    VERSION 1.0
    DESCRIPTION "EDMD for square well barrier system"
    LANGUAGES CXX
)

# Build type
# Debug
# set(CMAKE_BUILD_TYPE "Debug")
if (CMAKE_BUILD_TYPE STREQUAL "")
    # Build type is not set eg. command was "cmake .."
    message(STATUS "  Diag: Build type was unspecified, set to Release")
    set(CMAKE_BUILD_TYPE Release)
else ()
    message(STATUS "  Diag: Build type specified as '${CMAKE_BUILD_TYPE}'")
endif ()

set(CMAKE_CXX_FLAGS "-std=c++17 -fopenmp -ltbb")

message(STATUS "${CMAKE_CXX_FLAGS}")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

# Include directory
include_directories(${PROJECT_SOURCE_DIR}/include)

# Source directory
add_subdirectory(src)

# Excutable 
add_executable(${PROJECT_NAME} ${PROJECT_SOURCE_DIR}/bin/main.cpp)

# Link library
target_link_libraries(${PROJECT_NAME} EDMD)

