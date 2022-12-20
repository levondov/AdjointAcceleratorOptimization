# CMakeList.txt : CMake project for AdjointFTR, include source and define
# project specific logic here.
#

# CMakeList.txt : CMake project for AdjointFTR, include source and define
# project specific logic here.
#

set( PROJECT_FOLDER ${PROJECT_FOLDER}/Sim )

add_subdirectory(src)

# Add source to this project's executable.
#add_executable (AdjointFTR "AdjointFTR.cpp" "AdjointFTR.h")

#if (CMAKE_VERSION VERSION_GREATER 3.12)
#  set_property(TARGET AdjointFTR PROPERTY CXX_STANDARD 20)
#endif()

# TODO: Add tests and install targets if needed.
