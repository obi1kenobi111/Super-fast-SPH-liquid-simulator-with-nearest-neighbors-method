
add_executable (SPH "SPH_Simulation.cpp")

find_package(GLUT REQUIRED)
include_directories(${PROJECT_SOURCE_DIR}/include_gl)


target_link_libraries (SPH solver ${PROJECT_SOURCE_DIR}/lib_gl)