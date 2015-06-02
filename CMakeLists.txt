cmake_minimum_required(VERSION 3.2)                                                                                                  

set(PROJECT_NAME kallisto-export)

project(${PROJECT_NAME})

include(ExternalProject)

set(CMAKE_CXX_COMPILER g++)

set(CMAKE_VERBOSE_MAKEFILE on)

add_compile_options(-std=c++11)

# Set default ExternalProject root directory
SET_DIRECTORY_PROPERTIES(PROPERTIES EP_PREFIX ${CMAKE_BINARY_DIR}/external)

ExternalProject_Add(
kallisto
                URL https://github.com/churchill-lab/kallisto/archive/master.zip
                INSTALL_COMMAND ""
)

ExternalProject_Get_Property(kallisto source_dir)
set(KALLISTO_INCLUDE_DIR ${source_dir})

MESSAGE(STATUS "KALLISTO_INCLUDE_DIR=" ${KALLISTO_INCLUDE_DIR})

# Library
ExternalProject_Get_Property(kallisto binary_dir)
set(KALLISTO_LIBRARY_PATH ${binary_dir}/src/libkallisto_core.a)

MESSAGE(STATUS "KALLISTO_LIBRARY_PATH=" ${KALLISTO_LIBRARY_PATH})




add_executable(${PROJECT_NAME} main.cpp)
add_dependencies(${PROJECT_NAME} kallisto)
include_directories(${KALLISTO_INCLUDE_DIR}/src)
target_link_libraries(${PROJECT_NAME} ${KALLISTO_LIBRARY_PATH})

find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries( ${PROJECT_NAME} ${ZLIB_LIBRARIES} )
else()
    message(FATAL_ERROR "zlib not found. Required for to output files" )
endif( ZLIB_FOUND )

find_package(HDF5)
if(HDF5_FOUND)
    include_directories( ${HDF5_INCLUDE_DIR} )
    target_link_libraries( ${PROJECT_NAME} ${HDF5_LIBRARIES} )
else()
    message(FATAL_ERROR "HDF5 not found. Required to output files")
endif(HDF5_FOUND)

