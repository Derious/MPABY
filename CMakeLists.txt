cmake_minimum_required (VERSION 2.8.11)
project (emp-agmpc)
set(NAME "emp-agmpc")

set(CMAKE_CXX_COMPILER_WORKS TRUE)
set(CMAKE_C_COMPILER_WORKS TRUE)

find_path(CMAKE_FOLDER NAMES cmake/emp-tool-config.cmake)
include(${CMAKE_FOLDER}/cmake/emp-base.cmake)

FIND_PACKAGE(emp-ot REQUIRED)
INCLUDE_DIRECTORIES(${EMP-OT_INCLUDE_DIRS})

FIND_PACKAGE(OpenMP REQUIRED)
if (OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()

FIND_PACKAGE(Eigen3 REQUIRED)
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIRS})

install(FILES cmake/emp-agmpc-config.cmake DESTINATION cmake/)
install(DIRECTORY emp-agmpc DESTINATION include/)

ENABLE_TESTING()
ADD_SUBDIRECTORY(test/test_GC)
ADD_SUBDIRECTORY(test/test_GMW)
ADD_SUBDIRECTORY(test/test_CORE)

