#include <jltinygismo/helper.hh>

#include <jlcxx/module.hpp>

#include <gismo.h>
#include <gsCore/gsGeometry.h>

#include "registergeometries.hh"

jlcxx::TypeWrapper<gismo::gsGeometry<double>> registerGeometry(jlcxx::Module& mod) {
  using jlcxx::arg;

  using JuliaMatrix = jlcxx::ArrayRef<double, 2>;
  using JuliaVector = jlcxx::ArrayRef<double, 1>;

  using Geometry = gismo::gsGeometry<>;
  auto geometry  = mod.add_type<Geometry>("gsGeometry");

  {
    geometry.method("coefsSize", &Geometry::coefsSize, arg("spline"));
    geometry.method("coefs", [](const Geometry& spline) { return spline.coefs(); }, arg("spline"));

    geometry.method(
        "eval!",
        [](const Geometry& basis, JuliaVector u, gismo::gsMatrix<>& out) { basis.eval_into(wrapVector(u), out); },
        arg("basis"), arg("u"), arg("out"));

    geometry.method(
        "deriv!",
        [](const Geometry& geometry, JuliaVector u, gismo::gsMatrix<>& out) {
          geometry.deriv_into(wrapVector(u), out);
        },
        arg("geometry"), arg("u"), arg("out"));

    geometry.method(
        "deriv2!",
        [](const Geometry& geometry, JuliaVector u, gismo::gsMatrix<>& out) {
          geometry.deriv2_into(wrapVector(u), out);
        },
        arg("geometry"), arg("u"), arg("out"));
    geometry.method(
        "_eval", [](const Geometry& basis, JuliaVector u) { return basis.eval(wrapVector(u)); }, arg("basis"),
        arg("u"));

    geometry.method(
        "deriv", [](const Geometry& geometry, JuliaVector u) { return geometry.deriv(wrapVector(u)); }, arg("geometry"),
        arg("u"));

    geometry.method(
        "deriv2", [](const Geometry& geometry, JuliaVector u) { return geometry.deriv2(wrapVector(u)); },
        arg("geometry"), arg("u"));
  }
  {
    geometry.method(
        "eval!",
        [](const Geometry& basis, JuliaMatrix u, gismo::gsMatrix<>& out) { basis.eval_into(wrapMatrix(u), out); },
        arg("basis"), arg("u"), arg("out"));

    geometry.method(
        "deriv!",
        [](const Geometry& geometry, JuliaMatrix u, gismo::gsMatrix<>& out) {
          geometry.deriv_into(wrapMatrix(u), out);
        },
        arg("geometry"), arg("u"), arg("out"));

    geometry.method(
        "deriv2!",
        [](const Geometry& geometry, JuliaMatrix u, gismo::gsMatrix<>& out) {
          geometry.deriv2_into(wrapMatrix(u), out);
        },
        arg("geometry"), arg("u"), arg("out"));
    geometry.method(
        "_eval", [](const Geometry& basis, JuliaMatrix u) { return basis.eval(wrapMatrix(u)); }, arg("basis"),
        arg("u"));

    geometry.method(
        "deriv", [](const Geometry& geometry, JuliaMatrix u) { return geometry.deriv(wrapMatrix(u)); }, arg("geometry"),
        arg("u"));

    geometry.method(
        "deriv2", [](const Geometry& geometry, JuliaMatrix u) { return geometry.deriv2(wrapMatrix(u)); },
        arg("geometry"), arg("u"));
  }

  geometry.method(
      "closestPointTo",
      [](const Geometry& geometry, JuliaVector pt, gismo::gsVector<>& result) {
        geometry.closestPointTo(wrapVector(pt), result);
      },
      arg("geometry"), arg("pt"), arg("result"));

  geometry.method(
      "coefAtCorner", [](Geometry& geometry, int c) { return gismo::gsVector<>{geometry.coefAtCorner(c)}; },
      arg("geometry"), arg("c"));
  {
    geometry.method(
        "jacobian!",
        [](const Geometry& geometry, JuliaVector u, gismo::gsMatrix<>& out) {
          geometry.jacobian_into(wrapVector(u), out);
        },
        arg("geometry"), arg("u"), arg("out"));

    geometry.method(
        "jacobian", [](const Geometry& geometry, JuliaVector u) { return geometry.jacobian(wrapVector(u)); },
        arg("geometry"), arg("u"));

    geometry.method(
        "hessian!",
        [](const Geometry& geometry, JuliaVector u, gismo::gsMatrix<>& out, int coord) {
          geometry.hessian_into(wrapVector(u), out, coord - 1);
        },
        arg("geometry"), arg("u"), arg("out"), arg("coord"));

    geometry.method(
        "hessian",
        [](const Geometry& geometry, JuliaVector u, int coord) { return geometry.hessian(wrapVector(u), coord - 1); },
        arg("geometry"), arg("u"), arg("coord"));
  }
  {
    geometry.method(
        "jacobian!",
        [](const Geometry& geometry, JuliaMatrix u, gismo::gsMatrix<>& out) {
          geometry.jacobian_into(wrapMatrix(u), out);
        },
        arg("geometry"), arg("u"), arg("out"));

    geometry.method(
        "jacobian", [](const Geometry& geometry, JuliaMatrix u) { return geometry.jacobian(wrapMatrix(u)); },
        arg("geometry"), arg("u"));

    geometry.method(
        "hessian!",
        [](const Geometry& geometry, JuliaMatrix u, gismo::gsMatrix<>& out, int coord) {
          geometry.hessian_into(wrapMatrix(u), out, coord - 1);
        },
        arg("geometry"), arg("u"), arg("out"), arg("coord"));

    geometry.method(
        "hessian",
        [](const Geometry& geometry, JuliaMatrix u, int coord) { return geometry.hessian(wrapMatrix(u), coord - 1); },
        arg("geometry"), arg("u"), arg("coord"));
  }

  geometry.method("targetDim", &Geometry::targetDim);
  geometry.method("coefDim", &Geometry::coefDim);
  geometry.method("geoDim", &Geometry::geoDim);
  geometry.method("parDim", &Geometry::parDim);

  return geometry;
}