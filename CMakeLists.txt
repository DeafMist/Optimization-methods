cmake_minimum_required(VERSION 3.17)

project(OptimizationMethods)

set(SOURCES zad3.cpp)

add_executable(OptimizationMethods ${SOURCES})

set_target_properties(
        OptimizationMethods PROPERTIES
        CXX_VERSION 17
        CXX_VERSION_REQUIRED ON
)