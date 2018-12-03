cmake_minimum_required(VERSION 3.0)

project(Macroc VERSION 0.1 LANGUAGES C)

find_package(MPI REQUIRED)

find_library(
	PETSC_LIB petsc
	PATHS ENV PETSC_DIR
	PATH_SUFFIXES "$ENV{PETSC_ARCH}/lib" "lib"
	)

find_path(
	PETSC_INC_A petsc.h
	PATHS ENV PETSC_DIR
	PATH_SUFFIXES "include"
	)

find_path(
	PETSC_INC_B petscconf.h
	PATHS ENV PETSC_DIR
	PATH_SUFFIXES "$ENV{PETSC_ARCH}/include" "include"
	)

set(PETSC_INC ${PETSC_INC_A} ${PETSC_INC_B})

find_path(
	MICROPP_INC micro.hpp
	PATHS ENV MICROPP_DIR
	PATH_SUFFIXES "include"
	)

find_library(
	MICROPP_LIB micropp
	PATHS ENV MICROPP_DIR
	PATH_SUFFIXES "build"
	)

if(NOT MICROPP_LIB)
	  message(FATAL_ERROR "<libmicropp.a> library not found.")
endif()

if(NOT PETSC_LIB)
	message(FATAL_ERROR "<libpetsc.so> library not found.")
endif()

if(NOT MICROPP_INC)
	message(FATAL_ERROR "<include> directory of MicroPP was not found.")
endif()

SET(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ${C_COMPILE_FLAGS}")

include_directories(include ${CMAKE_BINARY_DIR} ${MICROPP_INC} ${PETSC_INC} ${MPI_INCLUDE_PATH})

file(GLOB SOURCES src/*.c)

link_libraries("-lm" "-lstdc++" ${PETSC_LIB} ${MPI_LIBRARIES} ${MICROPP_LIB})

add_executable(macroc ${SOURCES})

enable_testing()
add_subdirectory(tests)
