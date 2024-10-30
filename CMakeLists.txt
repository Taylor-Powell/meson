cmake_minimum_required(VERSION 3.21.3 FATAL_ERROR)

project(meson)

set(CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 20)

add_executable(meson
               src/main.cpp
               src/read_Ecm_dat.h
               src/read_Ecm_dat.cpp
               src/generic_funcs.h
               src/EvsQsq.h)

set(CMAKE_CXX_FLAGS "-O2")
set(CMAKE_C_FLAGS "-O2")
set(CMAKE_EXE_LINKER_FLAGS "-O2")

find_package(MPI REQUIRED)
message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")
target_link_libraries(meson MPI::MPI_CXX)