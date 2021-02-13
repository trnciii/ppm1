cmake_minimum_required(VERSION 3.1)
project(ppm1 CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-O3 -fopenmp")

add_executable(ppm1 ppm1.cpp)