#include <jltinygismo/helper.hh>

#include <jlcxx/array.hpp>
#include <jlcxx/jlcxx.hpp>

#include <gsNurbs/gsKnotVector.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsEigen/Eigen>
#include <string>

#include "registerbasis.hh"

struct WrapTensorBSplineBasis
{
  template <typename T>
  struct spline_dimension;

  template <int n>
  struct spline_dimension<gismo::gsTensorBSplineBasis<n>>
  {
    static constexpr int value = n;
  };

  template <typename BasisT>
  void operator()(BasisT&& basis) {
    using jlcxx::arg;
    using JuliaMatrix = jlcxx::ArrayRef<double, 2>;
    using JuliaVector = jlcxx::ArrayRef<double, 1>;

    using Basis     = typename BasisT::type;
    constexpr int n = spline_dimension<Basis>::value;

    if constexpr (n == 1)
      basis.constructor([](const gismo::gsKnotVector<>& knotVector) { return new Basis{knotVector}; });
    else if constexpr (n == 2)
      basis.constructor([](const gismo::gsKnotVector<>& knotVector1, const gismo::gsKnotVector<>& knotVector2) {
        return new Basis{knotVector1, knotVector2};
      });

    // Knots
    basis.method("knots", [](const Basis& basis, int i = 1) { return basis.knots(i - 1); }, arg("basis"), arg("i") = 1);
    basis.method("knot", &Basis::knot, arg("basis"), arg("i"));
    basis.method("size", [](const Basis& basis) { return basis.size(); }, arg("basis"));

    /*
        enum side { west  = 1,  east = 2, south = 3, north= 4, front=5, back=6,
                stime = 7, etime = 8,
                left  = 1, right = 2, down  = 3, up   = 4 , none = 0, all = 0 };
    */
    basis.method(
        "numElements", [](const Basis& basis, short_t side = 0) { return basis.numElements(gismo::boxSide{side}); },
        arg("basis"), arg("side") = 0);

    if constexpr (n > 1)
      basis.method("numTotalElements", &Basis::numTotalElements);

    // Degree
    if constexpr (n == 1) {
      basis.method("degree", [](const Basis& basis) { return basis.degree(); }, arg("basis"));
      basis.method("order", [](const Basis& basis) { return basis.order(); }, arg("basis"));
    } else {
      basis.method("degree", [](const Basis& basis, int i) { return basis.degree(i); }, arg("basis"), arg("i"));
    }

    // Refinement
    basis.method(
        "insertKnot", [](Basis& basis, double knot, int mult = 1) { basis.insertKnot(knot, mult); }, arg("basis"),
        arg("knot"), arg("mult") = 1);

    basis.method(
        "removeKnot", [](Basis& basis, double knot, int mult = 1) { basis.removeKnot(knot, mult); }, arg("basis"),
        arg("knot"), arg("mult") = 1);

    basis.method(
        "insertKnots",
        [](Basis& basis, JuliaVector knots) {
          auto knots2 = std::vector<std::vector<double>>{std::vector<double>(knots.size())};
          std::copy(knots.begin(), knots.end(), knots2.front().begin());
          basis.insertKnots(knots2);
        },
        arg("basis"), arg("knots"));

    // Actives
    if constexpr (n == 1)
      basis.method("numActive", [](Basis& basis) { return basis.numActive(); }, arg("basis"), arg("u"));
    else {
      basis.method(
          "numActive", [](Basis& basis, JuliaVector u) { return basis.numActive(wrapVector(u)); }, arg("basis"),
          arg("u"));
      basis.method(
          "numActive!",
          [](Basis& basis, JuliaVector u, gismo::gsVector<int>& out) { basis.numActive_into(wrapVector(u), out); },
          arg("basis"), arg("u"));
      basis.method("component", [](const Basis& basis, int i) { return basis.component(i - 1); }, arg("i"));
      basis.method("component", [](Basis& basis, int i) { return basis.component(i - 1); }, arg("i"));
    }

    // First Derivative
    if constexpr (n == 1)
      basis.method(
          "deriv!",
          [](Basis& basis, JuliaVector u, JuliaVector coefs, gismo::gsMatrix<>& out) {
            basis.deriv_into(wrapVector(u), wrapVector(coefs), out);
          },
          arg("basis"), arg("u"), arg("coefs"), arg("out"));

    if constexpr (n == 1)
      basis.method(
          "deriv!",
          [](Basis& basis, JuliaVector u, JuliaMatrix coefs, gismo::gsMatrix<>& out) {
            basis.deriv_into(wrapVector(u), wrapMatrix(coefs), out);
          },
          arg("basis"), arg("u"), arg("coefs"), arg("out"));

    // Second Derivatives
    if constexpr (n == 1)
      basis.method(
          "deriv2!",
          [](Basis& basis, JuliaVector u, JuliaVector coefs, gismo::gsMatrix<>& out) {
            basis.deriv2_into(wrapVector(u), wrapVector(coefs), out);
          },
          arg("basis"), arg("u"), arg("coefs"), arg("out"));
    if constexpr (n == 1)
      basis.method(
          "deriv2!",
          [](Basis& basis, JuliaVector u, JuliaMatrix coefs, gismo::gsMatrix<>& out) {
            basis.deriv2_into(wrapVector(u), wrapMatrix(coefs), out);
          },
          arg("basis"), arg("u"), arg("coefs"), arg("out"));

    // Geometry
    // basis.method("makeGeometry", [](const Basis& basis, JuliaMatrix coefs) {
    //   gsEigen::Map<const gsEigen::MatrixXd> coefsMat(coefs.data(), coefs.size());
    //   return basis.makeGeometry(coefsMat);
    // });
  }
};

namespace jlcxx {

template <int Val, typename T>
struct BuildParameterList<gismo::gsTensorBSplineBasis<Val, T>>
{
  // Using long to map to Julia Int64
  using type = ParameterList<std::integral_constant<long, Val>, T>;
};

template <int n>
struct SuperType<gismo::gsTensorBSplineBasis<n>>
{
  using type = gismo::gsBasis<>;
};

} // namespace jlcxx

void regsiterBSplineBasis(jlcxx::Module& mod, jlcxx::TypeWrapper<gismo::gsBasis<double>>& gsBasis) {
  mod.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>>>("BSplineBasis", gsBasis.dt())
      .apply<gismo::gsTensorBSplineBasis<1>>(WrapTensorBSplineBasis());

  mod.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>>>("TensorBSplineBasis", gsBasis.dt())
      .apply<gismo::gsTensorBSplineBasis<2>, gismo::gsTensorBSplineBasis<3>>(WrapTensorBSplineBasis());
}
