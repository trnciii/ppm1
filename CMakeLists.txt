cmake_minimum_required(VERSION 3.12)
project(ppm1 CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(STB_IMAGE REQUIRED)
include_directories(${STB_INCLUDE_DIRS})

find_package(OpenMP REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

if(APPLE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp-version=50")
	add_executable(ppm1 ppm1.cpp)
	target_link_libraries(ppm1 PRIVATE OpenMP::OpenMP_CXX)
else()
	add_executable(ppm1 ppm1.cpp)
endif()