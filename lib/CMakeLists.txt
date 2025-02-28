add_library(meson STATIC)

target_include_directories(meson PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(meson MPI::MPI_CXX)

set(MESON_HDRS
    EcmData.h
    matelem.h
    generic_funcs.h
    rotations.h
    readAndLoop.h
    Wigner.h
)

install(FILES ${MESON_HDRS} DESTINATION include/meson)

target_sources(meson PUBLIC
    generic_funcs.cpp
    EcmData.cpp
    matelem.cpp
    rotations.cpp
    readAndLoop.cpp
    Wigner.cpp
)

install(TARGETS meson)