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
    geometry.method(
        "eval!", [](Geometry& basis, JuliaVector u, gismo::gsMatrix<>& out) { basis.eval_into(wrapVector(u), out); },
        arg("basis"), arg("u"), arg("out"));

    geometry.method("deriv!", [](Geometry& geometry, JuliaVector u, gismo::gsMatrix<>& out) {
      geometry.deriv_into(wrapVector(u), out);
    }, arg("geometry"), arg("u"), arg("out"));

    geometry.method("deriv2!", [](Geometry& geometry, JuliaVector u, gismo::gsMatrix<>& out) {
      geometry.deriv2_into(wrapVector(u), out);
    }, arg("geometry"), arg("u"), arg("out"));
    geometry.method("_eval", [](Geometry& basis, JuliaVector u) { return basis.eval(wrapVector(u)); },
                    arg("basis"), arg("u"));

    geometry.method("deriv", [](Geometry& geometry, JuliaVector u) { return geometry.deriv(wrapVector(u)); },
                    arg("geometry"), arg("u"));

    geometry.method("deriv2", [](Geometry& geometry, JuliaVector u) { return geometry.deriv2(wrapVector(u)); },
                    arg("geometry"), arg("u"));
  }
  {
    geometry.method(
        "eval!", [](Geometry& basis, JuliaMatrix u, gismo::gsMatrix<>& out) { basis.eval_into(wrapMatrix(u), out); },
        arg("basis"), arg("u"), arg("out"));

    geometry.method("deriv!", [](Geometry& geometry, JuliaMatrix u, gismo::gsMatrix<>& out) {
      geometry.deriv_into(wrapMatrix(u), out);
    }, arg("geometry"), arg("u"), arg("out"));

    geometry.method("deriv2!", [](Geometry& geometry, JuliaMatrix u, gismo::gsMatrix<>& out) {
      geometry.deriv2_into(wrapMatrix(u), out);
    }, arg("geometry"), arg("u"), arg("out"));
    geometry.method("_eval", [](Geometry& basis, JuliaMatrix u) { return basis.eval(wrapMatrix(u)); },
                    arg("basis"), arg("u"));

    geometry.method("deriv", [](Geometry& geometry, JuliaMatrix u) { return geometry.deriv(wrapMatrix(u)); },
                    arg("geometry"), arg("u"));

    geometry.method("deriv2", [](Geometry& geometry, JuliaMatrix u) { return geometry.deriv2(wrapMatrix(u)); },
                    arg("geometry"), arg("u"));
  }

  return geometry;
}