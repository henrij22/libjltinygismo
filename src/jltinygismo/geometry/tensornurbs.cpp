#include <jltinygismo/helper.hh>

#include <jlcxx/module.hpp>

#include <gismo.h>
#include <gsNurbs/gsBSpline.h>

#include "registergeometries.hh"

struct WrapTensorNurbs
{
  template <typename T>
  struct nurbs_dimension;

  template <int n>
  struct nurbs_dimension<gismo::gsTensorNurbs<n>>
  {
    static constexpr int value = n;
  };

  template <typename NurbsT>
  void operator()(NurbsT&& nurbs) {
    using jlcxx::arg;
    using JuliaMatrix = jlcxx::ArrayRef<double, 2>;
    using JuliaVector = jlcxx::ArrayRef<double, 1>;

    using Nurbs     = typename NurbsT::type;
    constexpr int n = nurbs_dimension<Nurbs>::value;

    if constexpr (n == 1)
      nurbs.constructor([](const gismo::gsTensorNurbsBasis<n>& basis, JuliaVector coefs) {
        return new Nurbs{basis, wrapVector(coefs)};
      });

    nurbs.constructor([](const gismo::gsTensorNurbsBasis<n>& basis, JuliaMatrix coefs) {
      return new Nurbs{basis, wrapMatrix(coefs)};
    });

    nurbs.constructor(
        [](const gismo::gsTensorNurbsBasis<n>& basis, gismo::gsMatrix<> coefs) { return new Nurbs{basis, coefs}; });

    if constexpr (n == 2)
      nurbs.constructor([](const gismo::gsKnotVector<>& kv1, gismo::gsKnotVector<>& kv2, gismo::gsMatrix<> coefs) {
        return new Nurbs{kv1, kv2, coefs};
      });

    if constexpr (n == 3)
      nurbs.constructor([](const gismo::gsKnotVector<>& kv1, gismo::gsKnotVector<>& kv2, gismo::gsKnotVector<>& kv3,
                           gismo::gsMatrix<> coefs) { return new Nurbs{kv1, kv2, kv3, coefs}; });

    if constexpr (n == 4)
      nurbs.constructor([](const gismo::gsKnotVector<>& kv1, gismo::gsKnotVector<>& kv2, gismo::gsKnotVector<>& kv3,
                           gismo::gsKnotVector<>& kv4,
                           gismo::gsMatrix<> coefs) { return new Nurbs{kv1, kv2, kv3, kv4, coefs}; });

    // knots
    nurbs.method("knots", [](const Nurbs& nurbs, int i) { return nurbs.knots(i - 1); }, arg("nurbs"), arg("i"));

    // nurbs.method("weight", [](Nurbs& nurbs, int i) { return nurbs.weight(i - 1); }, arg("nurbs"), arg("i"));
    nurbs.method("weights", [](Nurbs& nurbs) { return nurbs.weights(); }, arg("nurbs"));

    nurbs.method(
        "insertKnot", [](Nurbs& nurbs, double knot, int dir, int mul = 1) { nurbs.insertKnot(knot, dir - 1, mul); },
        arg("nurbs"), arg("knot"), arg("dir"), arg("mul") = 1);

    nurbs.method(
        "uniformRefine",
        [](Nurbs& nurbs, int numKnots = 1, int mul = 1, const int dir = 0) {
          nurbs.uniformRefine(numKnots, mul, dir - 1);
        },
        arg("nurbs"), arg("numKnots") = 1, arg("mul") = 1, arg("dir") = 0);

    nurbs.method(
        "uniformCoarsen", [](Nurbs& nurbs, int numKnots = 1) { nurbs.uniformCoarsen(numKnots); }, arg("nurbs"),
        arg("numKnots") = 1);

    // Degree
    nurbs.method("degree", [](const Nurbs& nurbs, int i) { return nurbs.degree(i - 1); }, arg("nurbs"), arg("i"));
    nurbs.method(
        "degreeElevate", [](Nurbs& nurbs, int i = 1, int dir = 0) { nurbs.degreeElevate(i, dir - 1); }, arg("basis"),
        arg("i") = 1, arg("dir") = 0);
    nurbs.method(
        "degreeReduce", [](Nurbs& nurbs, int i = 1, int dir = 0) { nurbs.degreeReduce(i, dir - 1); }, arg("basis"),
        arg("i") = 1, arg("dir") = 0);
    nurbs.method(
        "degreeIncrease", [](Nurbs& nurbs, int i = 1, int dir = 0) { nurbs.degreeIncrease(i, dir - 1); }, arg("basis"),
        arg("i") = 1, arg("dir") = 0);
    nurbs.method(
        "degreeDecrease", [](Nurbs& nurbs, int i = 1, int dir = 0) { nurbs.degreeDecrease(i, dir - 1); }, arg("basis"),
        arg("i") = 1, arg("dir") = 0);

    nurbs.method(
        "active!",
        [](Nurbs& basis, JuliaVector u, gismo::gsMatrix<int>& out) {
          basis.active_into(wrapVector(u), out);
          incrementByOne(out);
        },
        arg("basis"), arg("u"), arg("out"));

    // Basis
    nurbs.method("basis", [](Nurbs& nurbs) -> gismo::gsTensorNurbsBasis<n>& { return nurbs.basis(); }, arg("nurbs"));
    nurbs.method("boundary", [](Nurbs& nurbs, int c) { return nurbs.boundary(c); }, arg("nurbs"), arg("c"));
  }
};

namespace jlcxx {

template <int Val, typename T>
struct BuildParameterList<gismo::gsTensorNurbs<Val, T>>
{
  // Using long to map to Julia Int64
  using type = ParameterList<std::integral_constant<long, Val>, T>;
};

template <int n>
struct SuperType<gismo::gsTensorNurbs<n>>
{
  using type = gismo::gsGeometry<>;
};

} // namespace jlcxx

void registerTensorNurbs(jlcxx::Module& mod, jlcxx::TypeWrapper<gismo::gsGeometry<double>>& gsGeometry) {
  mod.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>>>("TensorNurbs", gsGeometry.dt())
      .apply</*gismo::gsTensorNurbs<1>,*/ gismo::gsTensorNurbs<2>, gismo::gsTensorNurbs<3>>(WrapTensorNurbs());
}
