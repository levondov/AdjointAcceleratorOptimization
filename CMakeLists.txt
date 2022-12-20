# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.12)

set( PROJECT_BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR} )
set( PROJECT_NAME "AAO" )

project (${PROJECT_NAME})

#python
find_package(Python COMPONENTS Interpreter Development)
include_directories(${Python_INCLUDE_DIRS})

# Include sub-projects.
add_subdirectory (AAO)