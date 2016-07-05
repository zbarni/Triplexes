
cmake_minimum_required (VERSION 2.6)
project (triplexator)

if (OPENMP_FOUND)
# We need a configuration header file.  First, configure this file, then
# add the target directory to the include path.
    configure_file (${triplexator_SOURCE_DIR}/config.h.in
            ${triplexator_BINARY_DIR}/config.h)

# Enable OpenMP.
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

elseif (OPENMP_FOUND)
    message(STATUS "Building triplexator without OpenMP support.")
endif (OPENMP_FOUND)

# Configure SHARED_LIBRARY define
if (SHAREDLIBRARY)
    set(SHARED_LIBRARY "1")
else ()
    set(SHARED_LIBRARY "0")
endif()

configure_file (${triplexator_SOURCE_DIR}/config.h.in
    ${triplexator_SOURCE_DIR}/config.h)

include_directories(${triplexator_BINARY_DIR})
include_directories(${triplexator_SOURCE_DIR})

if (SHAREDLIBRARY)
    message(STATUS "Building triplexator shared library instead of executable.")
    add_library(triplexator SHARED triplexator.cpp 
            triplexator.h
            triplex.h
            helper.h
            find_index_qgrams.h
            gardener.h
            edlib.cpp edlib.h
            triplex_alphabet.h
            triplex_functors.h
            triplex_pattern.h
            )
else()
    add_executable(triplexator triplexator.cpp 
            triplexator.h
            triplex.h
            helper.h
            find_index_qgrams.h
            gardener.h
            edlib.cpp edlib.h
            triplex_alphabet.h
            triplex_functors.h
            triplex_pattern.h
            )
endif()

# The boost component iostreams is optional and only be used to compress output files on request
if (Boost_IOSTREAMS_FOUND)
    set(Boost_USE_STATIC_LIBS ON)
    include_directories(${Boost_INCLUDE_DIRS})
    target_link_libraries(triplexator ${Boost_IOSTREAMS_LIBRARIES})
    message("Boost_IOSTREAMS_LIBRARY_RELEASE = ${Boost_IOSTREAMS_LIBRARY_RELEASE}")
endif (Boost_IOSTREAMS_FOUND)

if(UNIX AND NOT APPLE)
    target_link_libraries(triplexator rt)
endif()

if (SHAREDLIBRARY)
    install(TARGETS triplexator 
            RUNTIME DESTINATION bin 
            LIBRARY DESTINATION triplexator
            COMPONENT applications) 
else()
    install(TARGETS triplexator 
            RUNTIME DESTINATION bin 
            COMPONENT applications) 
endif()
