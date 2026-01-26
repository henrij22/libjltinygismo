# SPDX-FileCopyrightText: 2026 Henrik Jakob jakob@ibb.uni-stuttgart.de SPDX-License-Identifier:
# GPL-3.0-or-later

add_subdirectory(jltinygismo)

set(JLTINYGISMO_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/jltinygismo)
set(JLTINYGISMO_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/jltinygismo)

file(GLOB_RECURSE JLTINYGISMO_HEADERS ${JLTINYGISMO_INCLUDE_DIR}/*.hh)

# Explicitly list the source files
set(JLTINYGISMO_SOURCES
    ${JLTINYGISMO_SOURCE_DIR}/gismo.cpp
    ${JLTINYGISMO_SOURCE_DIR}/basis/gsBasis.cpp
    ${JLTINYGISMO_SOURCE_DIR}/basis/bsplinebasis.cpp
    ${JLTINYGISMO_SOURCE_DIR}/basis/nurbsbasis.cpp
    ${JLTINYGISMO_SOURCE_DIR}/basis/knotvector.cpp
    ${JLTINYGISMO_SOURCE_DIR}/geometry/bspline.cpp
    ${JLTINYGISMO_SOURCE_DIR}/geometry/nurbs.cpp
    ${JLTINYGISMO_SOURCE_DIR}/geometry/tensorbspline.cpp
    ${JLTINYGISMO_SOURCE_DIR}/geometry/tensornurbs.cpp
    ${JLTINYGISMO_SOURCE_DIR}/geometry/gsGeometry.cpp
    ${JLTINYGISMO_SOURCE_DIR}/utility/utilities.cpp
    ${JLTINYGISMO_SOURCE_DIR}/utility/nurbsCreator.cpp
    ${JLTINYGISMO_SOURCE_DIR}/utility/filereader.cpp
    ${JLTINYGISMO_SOURCE_DIR}/utility/paraview.cpp
)


add_library(jltinygismo SHARED ${JLTINYGISMO_SOURCES})

set_target_properties(${JLCXX_TARGET} PROPERTIES PUBLIC_HEADER "${JLTINYGISMO_HEADERS}")

target_include_directories(jltinygismo PUBLIC "${CMAKE_SOURCE_DIR}/src")
target_include_directories(jltinygismo PUBLIC ${CMAKE_INSTALL_PREFIX}/include/gismo)

target_link_libraries(jltinygismo PUBLIC gismo JlCxx::cxxwrap_julia)
target_compile_features(jltinygismo PUBLIC cxx_std_20)

if(UNIX AND NOT APPLE)
   target_link_libraries(jltinygismo PUBLIC dl)
endif(UNIX AND NOT APPLE)

install(
  TARGETS jltinygismo
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION lib
)

