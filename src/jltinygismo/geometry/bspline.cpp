#include <jltinygismo/helper.hh>

#include <jlcxx/module.hpp>

#include <gismo.h>
#include <gsNurbs/gsBSpline.h>

#include "registergeometries.hh"

void registerBSpline(jlcxx::Module& mod, jlcxx::TypeWrapper<gismo::gsGeometry<double>>& geo) {
  using jlcxx::arg;
  using JuliaMatrix = jlcxx::ArrayRef<double, 2>;
  using JuliaVector = jlcxx::ArrayRef<double, 1>;

  using BSpline = gismo::gsBSpline<>;

  auto spline = mod.add_type<BSpline>("BSpline", geo.dt());

  spline.constructor(
      [](const gismo::gsBSplineBasis<>& basis, JuliaVector coefs) { return new BSpline{basis, wrapVector(coefs)}; });

  spline.constructor(
      [](const gismo::gsBSplineBasis<>& basis, JuliaMatrix coefs) { return new BSpline{basis, wrapMatrix(coefs)}; });

  spline.method("numCoefs", [](const BSpline& spline) { return spline.numCoefs(); }, arg("spline"));

  // knots
  spline.method("knots", [](const BSpline& spline) { return spline.knots(); }, arg("spline"));

  spline.method(
      "insertKnot!", [](BSpline& spline, double knot) { spline.insertKnot(knot); }, arg("spline"), arg("knot"));

  spline.method(
      "uniformRefine!",
      [](BSpline& spline, int numKnots = 1, int mul = 1, const int dir = -1) {
        spline.uniformRefine(numKnots, mul, dir);
      },
      arg("spline"), arg("numKnots") = 1, arg("mul") = 1, arg("dir") = -1);

  spline.method(
      "uniformCoarsen!", [](BSpline& spline, int numKnots = 1) { spline.uniformCoarsen(numKnots); }, arg("spline"),
      arg("numKnots") = 1);

  // Degree
  spline.method("degree", [](const BSpline& spline) { return spline.degree(); }, arg("spline"));
  spline.method(
      "degreeElevate!", [](BSpline& spline, const int i = 1) { spline.degreeElevate(i); }, arg("spline"), arg("i") = 1);
  spline.method(
      "degreeReduce!", [](BSpline& spline, const int i = 1) { spline.degreeReduce(i); }, arg("spline"), arg("i") = 1);
  spline.method(
      "degreeIncrease!", [](BSpline& spline, const int i = 1) { spline.degreeIncrease(i); }, arg("spline"),
      arg("i") = 1);
  spline.method(
      "degreeDecrease!", [](BSpline& spline, const int i = 1) { spline.degreeDecrease(i); }, arg("spline"),
      arg("i") = 1);

  spline.method(
      "active!",
      [](const BSpline& basis, JuliaVector u, gismo::gsMatrix<int>& out) {
        basis.active_into(wrapVector(u), out);
        incrementByOne(out);
      },
      arg("basis"), arg("u"), arg("out"));

  // Basis
  spline.method(
      "basis", [](const BSpline& spline) -> const gismo::gsBSplineBasis<>& { return spline.basis(); }, arg("spline"));
  spline.method("boundary", [](const BSpline& spline, int c) { return spline.boundary(c); }, arg("spline"), arg("c"));

  /*
      enum corner { southwestfront = 1, southeastfront = 2, northwestfront = 3, northeastfront = 4,
                southwestback  = 5, southeastback  = 6, northwestback  = 7, northeastback  = 8,
                southwest      = 1, southeast      = 2, northwest      = 3, northeast      = 4
  };
  */
  spline.method(
      "coefAtCorner", [](const BSpline& spline, int c) { return gismo::gsVector<>{spline.coefAtCorner(c)}; },
      arg("spline"), arg("c"));
}

namespace jlcxx {

template <>
struct SuperType<gismo::gsBSpline<>>
{
  using type = gismo::gsGeometry<>;
};

} // namespace jlcxx