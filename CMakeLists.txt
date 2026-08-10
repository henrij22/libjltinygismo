
cmake_minimum_required(VERSION 3.18)

project(JlTinyGismo)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Defines CMAKE_INSTALL_INCLUDEDIR and friends. Without it the header install rules
# below expand to an absolute /<pkg>/... path and `cmake --install` tries to write to
# the filesystem root.
include(GNUInstallDirs)

find_package(JlCxx REQUIRED)
find_package(gismo REQUIRED)

include_directories( ${GISMO_INCLUDE_DIRS})
link_directories   ( ${GISMO_LIBRARY_DIR} )

include( ${GISMO_USE_FILE} )
set(CMAKE_MODULE_PATH ${GISMO_MODULE_PATH})
# include(gsConfig) #compiler configuration

get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${JlCxx_location}")

message(STATUS "Found JlCxx at ${JlCxx_location}")

option(JLTINYGISMO_BUILD_TESTS "Build the C++ unit tests (requires doctest, fetched if absent)" ON)

# A cross build cannot run the resulting test binary, and the BinaryBuilder sandbox that
# produces the JLL has no network for FetchContent to pull doctest with. Turn the tests
# off there rather than making every release recipe remember to pass the flag.
if(JLTINYGISMO_BUILD_TESTS AND CMAKE_CROSSCOMPILING)
  message(STATUS "Cross-compiling: not building the C++ tests")
  set(JLTINYGISMO_BUILD_TESTS OFF)
endif()

add_subdirectory(src)
# add_subdirectory(cmake)

if(JLTINYGISMO_BUILD_TESTS)
  enable_testing()
  add_subdirectory(test)
endif()
