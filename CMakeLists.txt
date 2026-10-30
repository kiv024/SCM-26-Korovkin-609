cmake_minimum_required(VERSION 3.24)

project(EqSolver VERSION 0.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(OpenMP REQUIRED)
add_definitions(-Wall -O3 -fopenmp)
set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")

include_directories("${PROJECT_SOURCE_DIR}/include")
file(GLOB CPP_SOURCES "source/*.cpp")
add_executable(EqSolver ${CPP_SOURCES} main.cpp)