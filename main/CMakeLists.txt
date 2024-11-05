set(EXES
    read_Ecm
)

foreach(EXE ${EXES})
    add_executable(${EXE} ${EXE}.cpp)
    target_link_libraries(${EXE} MPI::MPI_CXX)
    target_link_libraries(${EXE} meson)
endforeach()