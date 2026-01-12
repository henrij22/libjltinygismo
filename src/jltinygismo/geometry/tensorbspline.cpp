#include <jltinygismo/helper.hh>

#include <jlcxx/module.hpp>

#include <gismo.h>
#include <gsNurbs/gsBSpline.h>

#include "registergeometries.hh"

struct WrapTensorBSpline
{
  template <typename T>
  struct spline_dimension;

  template <int n>
  struct spline_dimension<gismo::gsTensorBSpline<n>>
  {
    static constexpr int value = n;
  };

  template <typename BSplineT>
  void operator()(BSplineT&& spline) {
    using jlcxx::arg;
    using JuliaMatrix = jlcxx::ArrayRef<double, 2>;
    using JuliaVector = jlcxx::ArrayRef<double, 1>;

    using BSpline   = typename BSplineT::type;
    constexpr int n = spline_dimension<BSpline>::value;

    if constexpr (n == 1)
      spline.constructor([](const gismo::gsTensorBSplineBasis<n>& basis, JuliaVector coefs) {
        return new BSpline{basis, wrapVector(coefs)};
      });

    spline.constructor([](const gismo::gsTensorBSplineBasis<n>& basis, JuliaMatrix coefs) {
      return new BSpline{basis, wrapMatrix(coefs)};
    });

    spline.constructor(
        [](const gismo::gsTensorBSplineBasis<n>& basis, gismo::gsMatrix<> coefs) { return new BSpline{basis, coefs}; });

    if constexpr (n == 2)
      spline.constructor([](const gismo::gsKnotVector<>& kv1, gismo::gsKnotVector<>& kv2, gismo::gsMatrix<> coefs) {
        return new BSpline{kv1, kv2, coefs};
      });

    if constexpr (n == 3)
      spline.constructor([](const gismo::gsKnotVector<>& kv1, gismo::gsKnotVector<>& kv2, gismo::gsKnotVector<>& kv3,
                            gismo::gsMatrix<> coefs) { return new BSpline{kv1, kv2, kv3, coefs}; });

    if constexpr (n == 4)
      spline.constructor([](const gismo::gsKnotVector<>& kv1, gismo::gsKnotVector<>& kv2, gismo::gsKnotVector<>& kv3,
                            gismo::gsKnotVector<>& kv4,
                            gismo::gsMatrix<> coefs) { return new BSpline{kv1, kv2, kv3, kv4, coefs}; });

    // corner is 4x3
    spline.constructor([](JuliaMatrix corner, const gismo::gsKnotVector<>& kv1, gismo::gsKnotVector<>& kv2) {
      return new BSpline{wrapMatrix(corner), kv1, kv2};
    });

    // knots
    spline.method("knots", [](const BSpline& spline, int i) { return spline.knots(i); });

    spline.method(
        "insertKnot", [](BSpline& spline, double knot, int dir, int mul = 1) { spline.insertKnot(knot, dir, mul); },
        arg("spline"), arg("knot"), arg("dir"), arg("mul") = 1);

    spline.method(
        "uniformRefine",
        [](BSpline& spline, int numKnots = 1, int mul = 1, const int dir = -1) {
          spline.uniformRefine(numKnots, mul, dir);
        },
        arg("spline"), arg("numKnots") = 1, arg("mul") = 1, arg("dir") = -1);

    spline.method(
        "uniformRefine", [](BSpline& spline, int numKnots = 1) { spline.uniformCoarsen(numKnots); }, arg("spline"),
        arg("numKnots") = 1);

    // Degree
    spline.method("degree", [](const BSpline& spline, int i) { return spline.degree(i); });
    spline.method(
        "degreeElevate", [](BSpline& spline, const int i = 1) { spline.degreeElevate(i); }, arg("basis"), arg("i") = 1);
    spline.method(
        "degreeReduce", [](BSpline& spline, const int i = 1) { spline.degreeReduce(i); }, arg("basis"), arg("i") = 1);
    spline.method(
        "degreeIncrease", [](BSpline& spline, const int i = 1) { spline.degreeIncrease(i); }, arg("basis"),
        arg("i") = 1);
    spline.method(
        "degreeDecrease", [](BSpline& spline, const int i = 1) { spline.degreeDecrease(i); }, arg("basis"),
        arg("i") = 1);

    spline.method("active!", [](BSpline& basis, JuliaVector u, gismo::gsMatrix<int>& out) {
      basis.active_into(wrapVector(u), out);
    });
    {
      spline.method("eval!",
                    [](BSpline& basis, JuliaVector u, gismo::gsMatrix<>& out) { basis.eval_into(wrapVector(u), out); });

      spline.method("deriv!", [](BSpline& spline, JuliaVector u, gismo::gsMatrix<>& out) {
        spline.deriv_into(wrapVector(u), out);
      });

      spline.method("deriv2!", [](BSpline& spline, JuliaVector u, gismo::gsMatrix<>& out) {
        spline.deriv2_into(wrapVector(u), out);
      });
      spline.method("_eval", [](BSpline& basis, JuliaVector u) { return basis.eval(wrapVector(u)); });

      spline.method("deriv", [](BSpline& spline, JuliaVector u) { return spline.deriv(wrapVector(u)); });

      spline.method("deriv2", [](BSpline& spline, JuliaVector u) { return spline.deriv2(wrapVector(u)); });
    }
    {
      spline.method("eval!",
                    [](BSpline& basis, JuliaMatrix u, gismo::gsMatrix<>& out) { basis.eval_into(wrapMatrix(u), out); });

      spline.method("deriv!", [](BSpline& spline, JuliaMatrix u, gismo::gsMatrix<>& out) {
        spline.deriv_into(wrapMatrix(u), out);
      });

      spline.method("deriv2!", [](BSpline& spline, JuliaMatrix u, gismo::gsMatrix<>& out) {
        spline.deriv2_into(wrapMatrix(u), out);
      });
      spline.method("_eval", [](BSpline& basis, JuliaMatrix u) { return basis.eval(wrapMatrix(u)); });

      spline.method("deriv", [](BSpline& spline, JuliaMatrix u) { return spline.deriv(wrapMatrix(u)); });

      spline.method("deriv2", [](BSpline& spline, JuliaMatrix u) { return spline.deriv2(wrapMatrix(u)); });
    }
    // Basis
    spline.method("basis", [](BSpline& spline) -> gismo::gsTensorBSplineBasis<n>& { return spline.basis(); });
    spline.method("boundary", [](BSpline& spline, int c) { return spline.boundary(c); });

    spline.method("targetDim", &BSpline::targetDim);
    spline.method("coefDim", &BSpline::coefDim);
    spline.method("geoDim", &BSpline::geoDim);
    spline.method("parDim", &BSpline::parDim);

    spline.method("closestPointTo", [](BSpline& spline, JuliaVector pt, gismo::gsVector<>& result) {
      spline.closestPointTo(wrapVector(pt), result);
    });

    spline.method("coefAtCorner", [](BSpline& spline, int c) { return gismo::gsVector<>{spline.coefAtCorner(c)}; });
    {
      spline.method("jacobian!", [](BSpline& spline, JuliaVector u, gismo::gsMatrix<>& out) {
        spline.jacobian_into(wrapVector(u), out);
      });

      spline.method("jacobian", [](BSpline& spline, JuliaVector u) { return spline.jacobian(wrapVector(u)); });

      spline.method("hessian!", [](BSpline& spline, JuliaVector u, gismo::gsMatrix<>& out, int coord) {
        spline.hessian_into(wrapVector(u), out, coord);
      });

      spline.method(
          "hessian", [](BSpline& spline, JuliaVector u, int coord = 0) { return spline.hessian(wrapVector(u), coord); },
          arg("spline"), arg("u"), arg("coord") = 0);
    }
    {
      spline.method("jacobian!", [](BSpline& spline, JuliaMatrix u, gismo::gsMatrix<>& out) {
        spline.jacobian_into(wrapMatrix(u), out);
      });

      spline.method("jacobian", [](BSpline& spline, JuliaMatrix u) { return spline.jacobian(wrapMatrix(u)); });

      spline.method("hessian!", [](BSpline& spline, JuliaMatrix u, gismo::gsMatrix<>& out, int coord) {
        spline.hessian_into(wrapMatrix(u), out, coord);
      });

      spline.method(
          "hessian", [](BSpline& spline, JuliaMatrix u, int coord = 0) { return spline.hessian(wrapMatrix(u), coord); },
          arg("spline"), arg("u"), arg("coord") = 0);
    }
  }
};

namespace jlcxx {

template <int Val, typename T>
struct BuildParameterList<gismo::gsTensorBSpline<Val, T>>
{
  // Using long to map to Julia Int64
  using type = ParameterList<std::integral_constant<long, Val>, T>;
};

} // namespace jlcxx

void registerTensorBSpline(jlcxx::Module& mod, jlcxx::TypeWrapper<gismo::gsGeometry<double>>& gsGeometry) {
  mod.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>>>("TensorBSpline", gsGeometry.dt())
      .apply<gismo::gsTensorBSpline<1>, gismo::gsTensorBSpline<2>, gismo::gsTensorBSpline<3>>(WrapTensorBSpline());
}
