#include <jltinygismo/helper.hh>

#include <jlcxx/module.hpp>

#include <gismo.h>
#include <gsHSplines/gsTHBSpline.h>
#include <gsHSplines/gsTHBSplineBasis.h>

#include "registergeometries.hh"

/// Wraps gsTHBSpline<d,double,Trunc>, the geometry counterpart of gsTHBSplineBasis.
///
/// As on the basis side, gsHBSpline is `gsTHBSpline<d,T,false>` rather than a class of its own,
/// so truncated and plain hierarchical geometries come out of this one functor.
struct WrapHierarchicalSpline
{
  template <typename T>
  struct h_traits;

  template <int n, bool Trunc>
  struct h_traits<gismo::gsTHBSpline<n, double, Trunc>>
  {
    static constexpr int dim      = n;
    static constexpr bool truncated = Trunc;
  };

  template <typename SplineT>
  void operator()(SplineT&& spline) {
    using jlcxx::arg;
    using JuliaMatrix = jlcxx::ArrayRef<double, 2>;
    using JuliaBoxes  = jlcxx::ArrayRef<int64_t, 1>;

    using Spline    = typename SplineT::type;
    constexpr int d = h_traits<Spline>::dim;
    using Basis     = gismo::gsTHBSplineBasis<d, double, h_traits<Spline>::truncated>;

    // ---------------------------------------------------------------- constructors

    spline.constructor(
        [](const Basis& basis, JuliaMatrix coefs) { return new Spline{basis, gismo::gsMatrix<>{wrapMatrix(coefs)}}; });

    spline.constructor([](const Basis& basis, gismo::gsMatrix<>& coefs) { return new Spline{basis, coefs}; });

    // Lift a tensor B-spline geometry into the hierarchical setting, unchanged but now refinable
    // locally.
    spline.constructor([](const gismo::gsTensorBSpline<d>& tensorSpline) { return new Spline{tensorSpline}; });

    // ---------------------------------------------------------------- queries

    spline.method("numCoefs", [](const Spline& spline) { return spline.coefsSize(); }, arg("spline"));

    spline.method("basis", [](const Spline& spline) -> Basis { return spline.basis(); }, arg("spline"));

    spline.method(
        "degree", [](const Spline& spline, int i) { return spline.basis().degree(i - 1); }, arg("spline"), arg("i"));

    spline.method("boundary", [](const Spline& spline, int c) { return spline.boundary(c); }, arg("spline"), arg("c"));

    spline.method(
        "coefAtCorner", [](const Spline& spline, int c) { return gismo::gsVector<>{spline.coefAtCorner(c)}; },
        arg("spline"), arg("c"));

    // ---------------------------------------------------------------- refinement

    // These go through gsGeometry, which refines the basis and carries the control points along,
    // so the geometry is unchanged as a map -- only its representation gets finer.
    spline.method(
        "refineElements!",
        [](Spline& spline, JuliaBoxes boxes) { spline.refineElements(toGismoBoxes<d>(boxes, "refineElements!")); },
        arg("spline"), arg("boxes"));

    spline.method(
        "unrefineElements!",
        [](Spline& spline, JuliaBoxes boxes) { spline.unrefineElements(toGismoBoxes<d>(boxes, "unrefineElements!")); },
        arg("spline"), arg("boxes"));

    spline.method(
        "uniformRefine!",
        [](Spline& spline, int numKnots = 1, int mul = 1) { spline.uniformRefine(numKnots, mul); }, arg("spline"),
        arg("numKnots") = 1, arg("mul") = 1);

    // Refines the whole domain to the finest level present, giving back an equivalent tensor
    // B-spline geometry.
    // gsTHBSpline::convertToBSpline refines *this to the finest level before copying the result
    // out, so it silently rewrites the geometry it is called on. That would make a non-bang
    // method mutate its argument, against the convention everywhere else here, so it runs on a
    // copy and leaves the caller's geometry alone.
    spline.method(
        "convertToBSpline",
        [](const Spline& spline) {
          Spline scratch{spline};
          gismo::gsTensorBSpline<d> result;
          scratch.convertToBSpline(result);
          return result;
        },
        arg("spline"));

    spline.method(
        "increaseMultiplicity!",
        [](Spline& spline, int level, int dir, double knotValue, int mult = 1) {
          checkLevel(level, spline.basis().numLevels(), "increaseMultiplicity!");
          spline.increaseMultiplicity(level - 1, toGismoDirRequired(dir, d, "increaseMultiplicity!"), knotValue, mult);
        },
        arg("spline"), arg("level"), arg("dir"), arg("knotValue"), arg("mult") = 1);
  }
};

namespace jlcxx {

template <int d, typename T, bool Trunc>
struct BuildParameterList<gismo::gsTHBSpline<d, T, Trunc>>
{
  // Using long to map to Julia Int64
  using type = ParameterList<std::integral_constant<int64_t, d>, T>;
};

template <int d, bool Trunc>
struct SuperType<gismo::gsTHBSpline<d, double, Trunc>>
{
  using type = gismo::gsGeometry<>;
};

} // namespace jlcxx

void registerHierarchicalSpline(jlcxx::Module& mod, jlcxx::TypeWrapper<gismo::gsGeometry<double>>& gsGeometry) {
  mod.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>>>("THBSpline", gsGeometry.dt())
      .apply<gismo::gsTHBSpline<1>, gismo::gsTHBSpline<2>, gismo::gsTHBSpline<3>>(WrapHierarchicalSpline());

  mod.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>>>("HBSpline", gsGeometry.dt())
      .apply<gismo::gsHBSpline<1>, gismo::gsHBSpline<2>, gismo::gsHBSpline<3>>(WrapHierarchicalSpline());
}
