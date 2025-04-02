set(EXES
    subHelicity_test
    kinFactors_test
    EigenTest
)

foreach(EXE ${EXES})
    add_executable(${EXE} ${EXE}.cpp)
    target_link_libraries(${EXE} MPI::MPI_CXX)
    target_link_libraries(${EXE} meson)
    target_include_directories(${EXE} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../lib)
endforeach()