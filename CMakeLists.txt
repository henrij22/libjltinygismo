# SPDX-FileCopyrightText: 2026 Henrik Jakob jakob@ibb.uni-stuttgart.de SPDX-License-Identifier:
# GPL-3.0-or-later

cmake_minimum_required(VERSION 3.18)

project(JlTinyGismo)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(JlCxx REQUIRED)
find_package(gismo  REQUIRED)

include_directories( ${GISMO_INCLUDE_DIRS})
link_directories   ( ${GISMO_LIBRARY_DIR} )


include( ${GISMO_USE_FILE} )
set(CMAKE_MODULE_PATH ${GISMO_MODULE_PATH})
include(gsConfig) #compiler configuration

get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${JlCxx_location}")

message(STATUS "Found JlCxx at ${JlCxx_location}")


add_subdirectory(src)
# add_subdirectory(cmake)
