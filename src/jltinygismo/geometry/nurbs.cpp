#include <jltinygismo/helper.hh>

#include <jlcxx/module.hpp>

#include <gismo.h>
#include <gsNurbs/gsBSpline.h>

#include "registergeometries.hh"

void registerNurbs(jlcxx::Module& mod, jlcxx::TypeWrapper<gismo::gsGeometry<double>>& geo) {
  using jlcxx::arg;
  using JuliaMatrix = jlcxx::ArrayRef<double, 2>;
  using JuliaVector = jlcxx::ArrayRef<double, 1>;

  using Nurbs = gismo::gsNurbs<>;

  auto nurbs = mod.add_type<Nurbs>("Nurbs", geo.dt());

  nurbs.constructor(
      [](const gismo::gsNurbsBasis<>& basis, JuliaVector coefs) { return new Nurbs{basis, wrapVector(coefs)}; });

  nurbs.constructor(
      [](const gismo::gsNurbsBasis<>& basis, JuliaMatrix coefs) { return new Nurbs{basis, wrapMatrix(coefs)}; });

  // knots
  nurbs.method("knots", [](const Nurbs& nurbs) { return nurbs.knots(); }, arg("nurbs"));

  nurbs.method("weight", [](Nurbs& nurbs, int i) { return nurbs.weight(i - 1); }, arg("nurbs"), arg("i"));
  nurbs.method("weights", [](Nurbs& nurbs) { return nurbs.weights(); }, arg("nurbs"));

  nurbs.method("insertKnot", [](Nurbs& nurbs, double knot) { nurbs.insertKnot(knot); }, arg("nurbs"), arg("knot"));

  nurbs.method(
      "uniformRefine",
      [](Nurbs& nurbs, int numKnots = 1, int mul = 1, const int dir = -1) { nurbs.uniformRefine(numKnots, mul, dir); },
      arg("nurbs"), arg("numKnots") = 1, arg("mul") = 1, arg("dir") = -1);

  nurbs.method(
      "uniformCoarsen", [](Nurbs& nurbs, int numKnots = 1) { nurbs.uniformCoarsen(numKnots); }, arg("nurbs"),
      arg("numKnots") = 1);

  // Degree
  nurbs.method("degree", [](const Nurbs& nurbs) { return nurbs.degree(); }, arg("nurbs"));
  nurbs.method(
      "degreeElevate", [](Nurbs& nurbs, const int i = 1) { nurbs.degreeElevate(i); }, arg("nurbs"), arg("i") = 1);
  nurbs.method(
      "degreeReduce", [](Nurbs& nurbs, const int i = 1) { nurbs.degreeReduce(i); }, arg("nurbs"), arg("i") = 1);
  nurbs.method(
      "degreeIncrease", [](Nurbs& nurbs, const int i = 1) { nurbs.degreeIncrease(i); }, arg("nurbs"), arg("i") = 1);
  nurbs.method(
      "degreeDecrease", [](Nurbs& nurbs, const int i = 1) { nurbs.degreeDecrease(i); }, arg("nurbs"), arg("i") = 1);

  nurbs.method(
      "active!",
      [](Nurbs& basis, JuliaVector u, gismo::gsMatrix<int>& out) {
        basis.active_into(wrapVector(u), out);
        incrementByOne(out);
      },
      arg("basis"), arg("u"), arg("out"));

  // Basis
  nurbs.method("basis", [](Nurbs& nurbs) -> gismo::gsNurbsBasis<>& { return nurbs.basis(); }, arg("nurbs"));
  nurbs.method("boundary", [](Nurbs& nurbs, int c) { return nurbs.boundary(c); }, arg("nurbs"), arg("c"));

  /*
      enum corner { southwestfront = 1, southeastfront = 2, northwestfront = 3, northeastfront = 4,
                southwestback  = 5, southeastback  = 6, northwestback  = 7, northeastback  = 8,
                southwest      = 1, southeast      = 2, northwest      = 3, northeast      = 4
  };
  */
  nurbs.method(
      "coefAtCorner", [](Nurbs& nurbs, int c) { return gismo::gsVector<>{nurbs.coefAtCorner(c)}; }, arg("nurbs"),
      arg("c"));
}

namespace jlcxx {

template <>
struct SuperType<gismo::gsBSpline<>>
{
  using type = gismo::gsGeometry<>;
};

} // namespace jlcxx