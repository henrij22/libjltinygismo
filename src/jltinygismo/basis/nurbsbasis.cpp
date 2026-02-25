#include <jltinygismo/helper.hh>

#include <jlcxx/array.hpp>
#include <jlcxx/jlcxx.hpp>

#include <gsNurbs/gsKnotVector.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsEigen/Eigen>
#include <string>

#include "registerbasis.hh"

struct WrapTensorNURBSBasis
{
  template <typename T>
  struct spline_dimension;

  template <int n>
  struct spline_dimension<gismo::gsTensorNurbsBasis<n>>
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

    if constexpr (n == 1) {
      basis.constructor([](const gismo::gsKnotVector<>& knotVector) { return new Basis{knotVector}; });
      basis.constructor(
          [](const gismo::gsKnotVector<>& knotVector, JuliaVector w) { return new Basis{knotVector, wrapVector(w)}; });
      basis.constructor(
          [](const gismo::gsBSplineBasis<>& basis, JuliaVector w) { return new Basis{basis, wrapVector(w)}; });
    } else if constexpr (n == 2) {
      basis.constructor(
          [](const gismo::gsKnotVector<>& knotVector1, const gismo::gsKnotVector<>& knotVector2, JuliaMatrix w) {
            auto kvs = std::vector<gismo::gsKnotVector<>>{knotVector1, knotVector2};
            return new Basis{kvs, wrapMatrix(w)};
          });
    }

    // Knots
    if constexpr (n > 1) {
      basis.method("knots", [](const Basis& basis, int i) { return basis.knots(i - 1); }, arg("basis"), arg("i"));
      basis.method("knot", &Basis::knot, arg("basis"), arg("i"));
    } else
      basis.method("knots", [](const Basis& basis) { return basis.knots(); }, arg("basis"));

    basis.method("size", [](const Basis& basis) { return basis.size(); }, arg("basis"));

    if constexpr (n == 1) {
      // Refinement
      basis.method(
          "insertKnot!", [](Basis& basis, double knot, int mult = 1) { basis.insertKnot(knot, mult); }, arg("basis"),
          arg("knot"), arg("mult") = 1);
    }

    /*
        enum side { west  = 1,  east = 2, south = 3, north= 4, front=5, back=6,
                stime = 7, etime = 8,
                left  = 1, right = 2, down  = 3, up   = 4 , none = 0, all = 0 };
    */
    basis.method(
        "numElements", [](const Basis& basis, short_t side = 0) { return basis.numElements(gismo::boxSide{side}); },
        arg("basis"), arg("side") = 0);

    // Degree
    if constexpr (n == 1)
      basis.method("degree", [](const Basis& basis) { return basis.degree(); }, arg("basis"));
    else
      basis.method("degree", [](const Basis& basis, int i) { return basis.degree(i - 1); }, arg("basis"), arg("i"));

    // Actives
    if constexpr (n == 1)
      basis.method("numActive", [](const Basis& basis) { return basis.numActive(); }, arg("basis"), arg("u"));
    else {
      basis.method(
          "numActive", [](const Basis& basis, JuliaVector u) { return basis.numActive(wrapVector(u)); }, arg("basis"),
          arg("u"));
      basis.method(
          "numActive!",
          [](const Basis& basis, JuliaVector u, gismo::gsVector<int>& out) {
            basis.numActive_into(wrapVector(u), out);
          },
          arg("basis"), arg("u"));
    }

    // Geometry
    // basis.method("makeGeometry", [](const Basis& basis, JuliaMatrix coefs) {
    //   gsEigen::Map<const gsEigen::MatrixXd> coefsMat(coefs.data(), coefs.size());
    //   return basis.makeGeometry(coefsMat);
    // });
  }
};

namespace jlcxx {

template <int Val, typename T>
struct BuildParameterList<gismo::gsTensorNurbsBasis<Val, T>>
{
  // Using long to map to Julia Int64
  using type = ParameterList<std::integral_constant<int64_t, Val>, T>;
};

template <int n>
struct SuperType<gismo::gsTensorNurbsBasis<n>>
{
  using type = gismo::gsBasis<>;
};

} // namespace jlcxx

void registerNurbsBasis(jlcxx::Module& mod, jlcxx::TypeWrapper<gismo::gsBasis<double>>& gsBasis) {
  mod.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>>>("NurbsBasis", gsBasis.dt())
      .apply<gismo::gsTensorNurbsBasis<1>>(WrapTensorNURBSBasis());

  mod.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>>>("TensorNurbsBasis", gsBasis.dt())
      .apply<gismo::gsTensorNurbsBasis<2>, gismo::gsTensorNurbsBasis<3>>(WrapTensorNURBSBasis());

}
