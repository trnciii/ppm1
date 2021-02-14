cmake_minimum_required(VERSION 3.1)
project(ppm1 CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(STB_IMAGE REQUIRED)
include_directories(${STB_INCLUDE_DIRS})

if(STB_IMAGE_FOUND)
	set(CMAKE_CXX_STANDARD 17)
	set(CMAKE_CXX_FLAGS "-O3 -fopenmp")

	add_executable(ppm1 ppm1.cpp)
endif()