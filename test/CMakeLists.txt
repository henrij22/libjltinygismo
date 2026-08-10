
include(FetchContent)

# Prefer an installed doctest if there is one, otherwise fetch it.
find_package(doctest QUIET)

if(NOT doctest_FOUND)
  message(STATUS "doctest not found, fetching it")
  FetchContent_Declare(
    doctest
    GIT_REPOSITORY https://github.com/doctest/doctest.git
    GIT_TAG v2.4.11
    GIT_SHALLOW TRUE
  )
  FetchContent_MakeAvailable(doctest)
  list(APPEND CMAKE_MODULE_PATH "${doctest_SOURCE_DIR}/scripts/cmake")
endif()

include(doctest)

add_executable(jltinygismo_tests
  test_helper.cpp
)

target_link_libraries(jltinygismo_tests PRIVATE
  jltinygismo
  doctest::doctest
)

target_compile_features(jltinygismo_tests PRIVATE cxx_std_20)

# Registers each TEST_CASE with CTest individually, so a failure names the case.
doctest_discover_tests(jltinygismo_tests)
