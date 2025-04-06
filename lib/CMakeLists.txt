add_library(meson STATIC)

target_include_directories(meson PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}
    /usr/include/eigen3)
target_link_libraries(meson MPI::MPI_CXX)
target_link_libraries(meson Eigen3::Eigen)

set(MESON_HDRS
    matelem.h
    generic_funcs.h
    rotations.h
    outputKinFactors.h
    Wigner.h
    errorHandling.h
    timing.h
)

install(FILES ${MESON_HDRS} DESTINATION include/meson)

target_sources(meson PUBLIC
    generic_funcs.cpp
    matelem.cpp
    rotations.cpp
    outputKinFactors.cpp
    Wigner.cpp
)

install(TARGETS meson)