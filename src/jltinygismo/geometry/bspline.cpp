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

  spline.constructor([](const gismo::gsBSplineBasis<>& basis, JuliaVector coefs) {
    gsEigen::Map<const gsEigen::VectorXd> coefsMat(coefs.data(), coefs.size());
    return new BSpline{basis, coefsMat};
  });

  spline.constructor(
      [](const gismo::gsBSplineBasis<>& basis, JuliaMatrix coefs) { return new BSpline{basis, wrapMatrix(coefs)}; });

  spline.constructor(
      [](const gismo::gsBSplineBasis<>& basis, gismo::gsMatrix<> coefs) { return new BSpline{basis, coefs}; });

  // knots
  spline.method("knots", [](const BSpline& spline) { return spline.knots(); }, arg("spline"));

  spline.method("insertKnot", [](BSpline& spline, double knot) { spline.insertKnot(knot); }, arg("spline"),
                arg("knot"));

  spline.method(
      "uniformRefine",
      [](BSpline& spline, int numKnots = 1, int mul = 1, const int dir = -1) {
        spline.uniformRefine(numKnots, mul, dir);
      },
      arg("spline"), arg("numKnots") = 1, arg("mul") = 1, arg("dir") = -1);

  spline.method(
      "uniformCoarsen", [](BSpline& spline, int numKnots = 1) { spline.uniformCoarsen(numKnots); }, arg("spline"),
      arg("numKnots") = 1);

  // Degree
  spline.method("degree", [](const BSpline& spline) { return spline.degree(); }, arg("spline"));
  spline.method(
      "degreeElevate", [](BSpline& spline, const int i = 1) { spline.degreeElevate(i); }, arg("spline"), arg("i") = 1);
  spline.method(
      "degreeReduce", [](BSpline& spline, const int i = 1) { spline.degreeReduce(i); }, arg("spline"), arg("i") = 1);
  spline.method(
      "degreeIncrease", [](BSpline& spline, const int i = 1) { spline.degreeIncrease(i); }, arg("spline"),
      arg("i") = 1);
  spline.method(
      "degreeDecrease", [](BSpline& spline, const int i = 1) { spline.degreeDecrease(i); }, arg("spline"),
      arg("i") = 1);

  spline.method("active!", [](BSpline& basis, JuliaVector u, gismo::gsMatrix<int>& out) {
    basis.active_into(wrapVector(u), out);
  }, arg("basis"), arg("u"), arg("out"));

  spline.method("eval!",
                [](BSpline& basis, JuliaVector u, gismo::gsMatrix<>& out) { basis.eval_into(wrapVector(u), out); },
                arg("basis"), arg("u"), arg("out"));

  spline.method("deriv!",
                [](BSpline& spline, JuliaVector u, gismo::gsMatrix<>& out) { spline.deriv_into(wrapVector(u), out); },
                arg("spline"), arg("u"), arg("out"));

  spline.method("deriv2!",
                [](BSpline& spline, JuliaVector u, gismo::gsMatrix<>& out) { spline.deriv2_into(wrapVector(u), out); },
                arg("spline"), arg("u"), arg("out"));
  spline.method("_eval", [](BSpline& basis, JuliaVector u) { return basis.eval(wrapVector(u)); }, arg("basis"),
                arg("u"));

  spline.method("deriv", [](BSpline& spline, JuliaVector u) { return spline.deriv(wrapVector(u)); }, arg("spline"),
                arg("u"));

  spline.method("deriv2", [](BSpline& spline, JuliaVector u) { return spline.deriv2(wrapVector(u)); }, arg("spline"),
                arg("u"));

  // Basis
  spline.method("basis", [](BSpline& spline) -> gismo::gsBSplineBasis<>& { return spline.basis(); }, arg("spline"));
  spline.method("boundary", [](BSpline& spline, int c) { return spline.boundary(c); }, arg("spline"), arg("c"));

  spline.method("targetDim", &BSpline::targetDim);
  spline.method("coefDim", &BSpline::coefDim);
  spline.method("geoDim", &BSpline::geoDim);
  spline.method("parDim", &BSpline::parDim);

  spline.method("closestPointTo", [](BSpline& spline, JuliaVector pt, gismo::gsVector<>& result) {
    spline.closestPointTo(wrapVector(pt), result);
  }, arg("spline"), arg("pt"), arg("result"));

  spline.method("jacobian!", [](BSpline& spline, JuliaVector u, gismo::gsMatrix<>& out) {
    spline.jacobian_into(wrapVector(u), out);
  }, arg("spline"), arg("u"), arg("out"));

  spline.method("jacobian", [](BSpline& spline, JuliaVector u) { return spline.jacobian(wrapVector(u)); },
                arg("spline"), arg("u"));

  spline.method("hessian!", [](BSpline& spline, JuliaVector u, gismo::gsMatrix<>& out, int coord) {
    spline.hessian_into(wrapVector(u), out, coord);
  }, arg("spline"), arg("u"), arg("out"), arg("coord"));

  spline.method(
      "hessian", [](BSpline& spline, JuliaVector u, int coord = 1) { return spline.hessian(wrapVector(u), coord - 1); },
      arg("spline"), arg("u"), arg("coord") = 1);

  spline.method("coefAtCorner", [](BSpline& spline, int c) { return gismo::gsVector<>{spline.coefAtCorner(c)}; },
                arg("spline"), arg("c"));
}