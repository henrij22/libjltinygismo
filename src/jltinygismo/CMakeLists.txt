add_subdirectory(basis)
add_subdirectory(geometry)
add_subdirectory(utility)

install(FILES helper.hh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/jltinygismo)
