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
    spline.method("knots", [](const BSpline& spline, int i) { return spline.knots(i - 1); }, arg("spline"), arg("i"));

    spline.method(
        "insertKnot", [](BSpline& spline, double knot, int dir, int mul = 1) { spline.insertKnot(knot, dir - 1, mul); },
        arg("spline"), arg("knot"), arg("dir"), arg("mul") = 1);

    spline.method(
        "uniformRefine",
        [](BSpline& spline, int numKnots = 1, int mul = 1, const int dir = 0) {
          spline.uniformRefine(numKnots, mul, dir - 1);
        },
        arg("spline"), arg("numKnots") = 1, arg("mul") = 1, arg("dir") = 0);

    spline.method(
        "uniformCoarsen", [](BSpline& spline, int numKnots = 1) { spline.uniformCoarsen(numKnots); }, arg("spline"),
        arg("numKnots") = 1);

    // Degree
    spline.method("degree", [](const BSpline& spline, int i) { return spline.degree(i - 1); }, arg("spline"), arg("i"));
    spline.method(
        "degreeElevate", [](BSpline& spline, int i = 1, int dir = 0) { spline.degreeElevate(i, dir - 1); },
        arg("basis"), arg("i") = 1, arg("dir") = 0);
    spline.method(
        "degreeReduce", [](BSpline& spline, int i = 1, int dir = 0) { spline.degreeReduce(i, dir - 1); }, arg("basis"),
        arg("i") = 1, arg("dir") = 0);
    spline.method(
        "degreeIncrease", [](BSpline& spline, int i = 1, int dir = 0) { spline.degreeIncrease(i, dir - 1); },
        arg("basis"), arg("i") = 1, arg("dir") = 0);
    spline.method(
        "degreeDecrease", [](BSpline& spline, int i = 1, int dir = 0) { spline.degreeDecrease(i, dir - 1); },
        arg("basis"), arg("i") = 1, arg("dir") = 0);

    spline.method(
        "active!",
        [](BSpline& basis, JuliaVector u, gismo::gsMatrix<int>& out) {
          basis.active_into(wrapVector(u), out);
          incrementByOne(out);
        },
        arg("basis"), arg("u"), arg("out"));
    {
      spline.method(
          "eval!", [](BSpline& basis, JuliaVector u, gismo::gsMatrix<>& out) { basis.eval_into(wrapVector(u), out); },
          arg("basis"), arg("u"), arg("out"));

      spline.method(
          "deriv!",
          [](BSpline& spline, JuliaVector u, gismo::gsMatrix<>& out) { spline.deriv_into(wrapVector(u), out); },
          arg("spline"), arg("u"), arg("out"));

      spline.method(
          "deriv2!",
          [](BSpline& spline, JuliaVector u, gismo::gsMatrix<>& out) { spline.deriv2_into(wrapVector(u), out); },
          arg("spline"), arg("u"), arg("out"));
      spline.method(
          "_eval", [](BSpline& basis, JuliaVector u) { return basis.eval(wrapVector(u)); }, arg("basis"), arg("u"));

      spline.method(
          "deriv", [](BSpline& spline, JuliaVector u) { return spline.deriv(wrapVector(u)); }, arg("spline"), arg("u"));

      spline.method(
          "deriv2", [](BSpline& spline, JuliaVector u) { return spline.deriv2(wrapVector(u)); }, arg("spline"),
          arg("u"));
    }
    {
      spline.method(
          "eval!", [](BSpline& basis, JuliaMatrix u, gismo::gsMatrix<>& out) { basis.eval_into(wrapMatrix(u), out); },
          arg("basis"), arg("u"), arg("out"));

      spline.method(
          "deriv!",
          [](BSpline& spline, JuliaMatrix u, gismo::gsMatrix<>& out) { spline.deriv_into(wrapMatrix(u), out); },
          arg("spline"), arg("u"), arg("out"));

      spline.method(
          "deriv2!",
          [](BSpline& spline, JuliaMatrix u, gismo::gsMatrix<>& out) { spline.deriv2_into(wrapMatrix(u), out); },
          arg("spline"), arg("u"), arg("out"));
      spline.method(
          "_eval", [](BSpline& basis, JuliaMatrix u) { return basis.eval(wrapMatrix(u)); }, arg("basis"), arg("u"));

      spline.method(
          "deriv", [](BSpline& spline, JuliaMatrix u) { return spline.deriv(wrapMatrix(u)); }, arg("spline"), arg("u"));

      spline.method(
          "deriv2", [](BSpline& spline, JuliaMatrix u) { return spline.deriv2(wrapMatrix(u)); }, arg("spline"),
          arg("u"));
    }
    // Basis
    spline.method(
        "basis", [](BSpline& spline) -> gismo::gsTensorBSplineBasis<n>& { return spline.basis(); }, arg("spline"));
    spline.method("boundary", [](BSpline& spline, int c) { return spline.boundary(c); }, arg("spline"), arg("c"));

    spline.method("targetDim", &BSpline::targetDim);
    spline.method("coefDim", &BSpline::coefDim);
    spline.method("geoDim", &BSpline::geoDim);
    spline.method("parDim", &BSpline::parDim);

    spline.method(
        "closestPointTo",
        [](BSpline& spline, JuliaVector pt, gismo::gsVector<>& result) {
          spline.closestPointTo(wrapVector(pt), result);
        },
        arg("spline"), arg("pt"), arg("result"));

    spline.method(
        "coefAtCorner", [](BSpline& spline, int c) { return gismo::gsVector<>{spline.coefAtCorner(c)}; }, arg("spline"),
        arg("c"));
    {
      spline.method(
          "jacobian!",
          [](BSpline& spline, JuliaVector u, gismo::gsMatrix<>& out) { spline.jacobian_into(wrapVector(u), out); },
          arg("spline"), arg("u"), arg("out"));

      spline.method(
          "jacobian", [](BSpline& spline, JuliaVector u) { return spline.jacobian(wrapVector(u)); }, arg("spline"),
          arg("u"));

      spline.method(
          "hessian!",
          [](BSpline& spline, JuliaVector u, gismo::gsMatrix<>& out, int coord) {
            spline.hessian_into(wrapVector(u), out, coord - 1);
          },
          arg("spline"), arg("u"), arg("out"), arg("coord"));

      spline.method(
          "hessian", [](BSpline& spline, JuliaVector u, int coord) { return spline.hessian(wrapVector(u), coord - 1); },
          arg("spline"), arg("u"), arg("coord"));
    }
    {
      spline.method(
          "jacobian!",
          [](BSpline& spline, JuliaMatrix u, gismo::gsMatrix<>& out) { spline.jacobian_into(wrapMatrix(u), out); },
          arg("spline"), arg("u"), arg("out"));

      spline.method(
          "jacobian", [](BSpline& spline, JuliaMatrix u) { return spline.jacobian(wrapMatrix(u)); }, arg("spline"),
          arg("u"));

      spline.method(
          "hessian!",
          [](BSpline& spline, JuliaMatrix u, gismo::gsMatrix<>& out, int coord) {
            spline.hessian_into(wrapMatrix(u), out, coord - 1);
          },
          arg("spline"), arg("u"), arg("out"), arg("coord"));

      spline.method(
          "hessian", [](BSpline& spline, JuliaMatrix u, int coord) { return spline.hessian(wrapMatrix(u), coord - 1); },
          arg("spline"), arg("u"), arg("coord"));
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
