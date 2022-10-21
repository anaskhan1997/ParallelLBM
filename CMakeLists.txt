cmake_minimum_required(VERSION 3.0)
project(fluidthing VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_compile_definitions(gpp9)

include_directories(include)

add_executable(LBM src/main.cpp src/solver.cpp)
#src/demo.cpp
cmake_policy(SET CMP0072 NEW)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
# add_dependencies(fluidthing OpenGL::OpenGL)
include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} )
target_link_libraries(LBM ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} )

#SET(CMAKE_CXX_FLAGS "-lGL -lGLU -lglut")
