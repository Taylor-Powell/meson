set(EXES
    read_Ecm
    EvsQsq_out
    subHelicity_test
    Wigner_d_test
    Wigner_D_test
)

foreach(EXE ${EXES})
    add_executable(${EXE} ${EXE}.cpp)
    target_link_libraries(${EXE} MPI::MPI_CXX)
    target_link_libraries(${EXE} meson)
endforeach()