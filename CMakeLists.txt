cmake_minimum_required(VERSION 3.0)

project(Macroc VERSION 0.1 LANGUAGES C)

# MPI
find_package(MPI REQUIRED)
if (MPI_FOUND)
	#include_directories(SYSTEM ${MPI_INCLUDE_PATH})
	message("MPI C compiler    : " ${MPI_C_COMPILER})
	message("MPI Inc directory : " ${MPI_C_INCLUDE_PATH})
else (MPI_FOUND)
	message(SEND_ERROR "This application cannot compile without MPI")
endif (MPI_FOUND)

if(EXISTS "$ENV{PETSC_DIR}/include/" AND EXISTS "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/include/")
	set(PETSC_INC "$ENV{PETSC_DIR}/include" "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/include")
else()
	message(SEND_ERROR "PETSc includes not found")
endif()

if(EXISTS "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/libpetsc.so")
	set(PETSC_LIB "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/libpetsc.so")
else()
	message(SEND_ERROR "PETSc library not found")
endif()

if(EXISTS "$ENV{MICROPP_DIR}/include/")
	set(MICROPP_INC "$ENV{MICROPP_DIR}/include")
else()
	message(SEND_ERROR "MicroPP includes not found. Probably you should set MICROPP_DIR env variable.")
endif()

if(EXISTS "$ENV{MICROPP_DIR}/build_release/libmicropp.a")
	set(MICROPP_LIB "$ENV{MICROPP_DIR}/build_release/libmicropp.a")
else()
	message(SEND_ERROR "MicroPP library not found. Probably you should set MICROPP_DIR env variable.")
endif()


SET(C_COMPILE_FLAGS "-Wformat")
SET(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ${C_COMPILE_FLAGS}")

# Include Directories (for all targets)
include_directories(include ${CMAKE_BINARY_DIR} ${MICROPP_INC} ${PETSC_INC} ${MPI_INCLUDE_PATH})

# wildcard all the sources in src
file(GLOB SOURCES src/*.c)

link_libraries("-lm" "-lstdc++" ${PETSC_LIB} ${MPI_LIBRARIES} ${MICROPP_LIB})

add_executable(macroc ${SOURCES})
