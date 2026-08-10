//
// Unit tests for the pure C++ helpers in jltinygismo/helper.hh.
//
// Only the helpers that do not touch the Julia runtime can be exercised here; everything
// that allocates or wraps a Julia array needs a live `jl_init`, and the bindings themselves
// only exist once the module is loaded from Julia. Those are covered by the Julia
// integration suite in julia/runtests.jl -- see the README.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <jltinygismo/helper.hh>

#include <gismo.h>

TEST_CASE("toGismoDir maps 1-based directions onto the G+Smo convention") {
  // 0 is "all directions", which G+Smo spells -1.
  CHECK(toGismoDir(0, 2, "f") == -1);
  CHECK(toGismoDir(0, 3, "f") == -1);

  // 1-based on the Julia side, 0-based in G+Smo.
  CHECK(toGismoDir(1, 2, "f") == 0);
  CHECK(toGismoDir(2, 2, "f") == 1);
  CHECK(toGismoDir(3, 3, "f") == 2);
}

TEST_CASE("toGismoDir rejects out-of-range directions instead of passing them through") {
  // These are the values that used to reach G+Smo unchecked and index out of bounds. -1 is
  // the important one: it reads as "all directions" on a basis but segfaults on a geometry,
  // so it must never make it past this helper.
  CHECK_THROWS_AS(toGismoDir(-1, 2, "f"), std::runtime_error);
  CHECK_THROWS_AS(toGismoDir(-2, 2, "f"), std::runtime_error);
  CHECK_THROWS_AS(toGismoDir(3, 2, "f"), std::runtime_error);
  CHECK_THROWS_AS(toGismoDir(2, 1, "f"), std::runtime_error);

  SUBCASE("the message names the function and the valid range") {
    try {
      toGismoDir(7, 2, "degreeElevate!");
      FAIL("expected toGismoDir to throw");
    } catch (const std::runtime_error& e) {
      const std::string message{e.what()};
      CHECK(message.find("degreeElevate!") != std::string::npos);
      CHECK(message.find("1:2") != std::string::npos);
      CHECK(message.find("7") != std::string::npos);
    }
  }
}

TEST_CASE("toGismoDirRequired has no all-directions mode") {
  // insertKnot! and friends act on exactly one direction. `0` must be rejected here: it
  // would otherwise reach G+Smo as -1 and index out of bounds.
  CHECK(toGismoDirRequired(1, 2, "f") == 0);
  CHECK(toGismoDirRequired(2, 2, "f") == 1);

  CHECK_THROWS_AS(toGismoDirRequired(0, 2, "f"), std::runtime_error);
  CHECK_THROWS_AS(toGismoDirRequired(-1, 2, "f"), std::runtime_error);
  CHECK_THROWS_AS(toGismoDirRequired(3, 2, "f"), std::runtime_error);

  SUBCASE("the message does not advertise 0") {
    try {
      toGismoDirRequired(0, 2, "insertKnot!");
      FAIL("expected toGismoDirRequired to throw");
    } catch (const std::runtime_error& e) {
      const std::string message{e.what()};
      CHECK(message.find("insertKnot!") != std::string::npos);
      CHECK(message.find("1:2") != std::string::npos);
      CHECK(message.find("all directions") == std::string::npos);
    }
  }
}

TEST_CASE("toGismoDir accepts every direction of a trivariate object") {
  for (int dir = 0; dir <= 3; ++dir)
    CHECK_NOTHROW(toGismoDir(dir, 3, "f"));
  CHECK_THROWS_AS(toGismoDir(4, 3, "f"), std::runtime_error);
}

TEST_CASE("incrementByOne shifts G+Smo's 0-based indices to Julia's 1-based ones") {
  SUBCASE("matrix") {
    gismo::gsMatrix<int> mat(2, 3);
    mat << 0, 1, 2, 3, 4, 5;

    incrementByOne(mat);

    gismo::gsMatrix<int> expected(2, 3);
    expected << 1, 2, 3, 4, 5, 6;
    CHECK(mat == expected);
  }

  SUBCASE("vector") {
    gismo::gsVector<int> vec(3);
    vec << 0, 7, 41;

    incrementByOne(vec);

    CHECK(vec[0] == 1);
    CHECK(vec[1] == 8);
    CHECK(vec[2] == 42);
  }

  SUBCASE("empty matrix is left alone") {
    gismo::gsMatrix<int> mat(0, 0);
    CHECK_NOTHROW(incrementByOne(mat));
    CHECK(mat.size() == 0);
  }
}
