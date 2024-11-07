add_library(meson STATIC)

target_include_directories(meson PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(meson MPI::MPI_CXX)

set(MESON_HDRS
    EvsQsq.h
    read_Ecm_dat.h
    generic_funcs.h
)

install(FILES ${MESON_HDRS} DESTINATION include/meson)

target_sources(meson PUBLIC
    EvsQsq.cpp
    generic_funcs.cpp
    read_Ecm_dat.cpp
)

install(TARGETS meson)