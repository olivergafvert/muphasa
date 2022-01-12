# This file is an adaptation of its analog in Rivet
# see https://github.com/rivetTDA/rivet/

cmake_minimum_required (VERSION 3.1)
project (muphasa)
set (CMAKE_CXX_STANDARD 14)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-ftemplate-depth=1024 -Wall -Wextra -pedantic -fPIC")


find_package(Boost "1.58")

# So far, we haven’t been able to get OpenMP to work with OS X’s Clang, though this should be possible.
# See https://github.com/rivetTDA/rivet/issues/136
# For now, to avoid a distracting error message when running CMake in OS X, we use the following if statement.  
if(NOT APPLE)
	# Look for OpenMP and enable it if found.
	find_package(OpenMP)
	if (OPENMP_FOUND)
    		set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    		set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
        endif()
endif()


include_directories("${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/include" ${Boost_INCLUDE_DIR} )

add_definitions(-DMSGPACK_USE_BOOST)



add_executable (mph main.cpp)

