cmake_minimum_required(VERSION 3.0)

project(Macroc VERSION 0.1 LANGUAGES C)

find_package(MPI REQUIRED)

set(PETSC_INC "$ENV{PETSC_DIR}/include" "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/include")
set(PETSC_LIB "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/libpetsc.so")
set(MICROPP_INC "$ENV{MICROPP_DIR}/include")
set(MICROPP_LIB "$ENV{MICROPP_DIR}/build/libmicropp.a")

SET(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ${C_COMPILE_FLAGS}")

include_directories(include ${CMAKE_BINARY_DIR} ${MICROPP_INC} ${PETSC_INC} ${MPI_INCLUDE_PATH})

file(GLOB SOURCES src/*.c)

link_libraries("-lm" "-lstdc++" ${PETSC_LIB} ${MPI_LIBRARIES} ${MICROPP_LIB})

add_executable(macroc ${SOURCES})

enable_testing()
add_subdirectory(tests)
