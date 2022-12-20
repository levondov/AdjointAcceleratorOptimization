
set( PROJECT_FOLDER ${PROJECT_FOLDER}/AdjointFTR )

set( target AdjointFTR )

set( headers 
	AdjointFTR.h 
	)

set( sources 
	AdjointFTR.cpp 
	)

# Add source to this project's executable.
add_library ( ${target} SHARED ${sources} ${headers} )

TARGET_LINK_LIBRARIES( 
	${target} 
	${Python_LIBRARIES}
)  

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET AdjointFTR PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.