#Generate the local Makefile
if(CREATE_LOCAL_MAKEFILES)
    file(RELATIVE_PATH CISPBA_RELATIVE_PATH ${CISPBA_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
    configure_file(${CISPBA_SOURCE_DIR}/Scripts/Makefile.in "Wind/Makefile" @ONLY)
endif()

find_package(Eigen3 REQUIRED)
add_library(eigen INTERFACE IMPORTED)

add_cispba_executable(Wind main.cpp basic_math.cpp basic_math.h camera.cpp camera.h clock.h constants.cpp constants.h
        custom_output.h fps.cpp fps.h grid_data.cpp grid_data.h grid_data_matrix.h mac_grid.cpp mac_grid.h matrix.h
        open_gl_headers.h smoke_sim.cpp smoke_sim.h timer.h vec.cpp vec.h)

target_include_directories(Wind SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIR})
target_link_libraries(Wind partio eigen)