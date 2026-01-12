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
    basis.method("knots", [](const Basis& basis, int i = 0) { return basis.knots(i); }, arg("basis"), arg("i") = 0);
    basis.method("knot", &Basis::knot);
    basis.method("size", [](const Basis& basis) { return basis.size(); });
    basis.method(
        "numElements", [](const Basis& basis, short_t side = 0) { return basis.numElements(side); }, arg("basis"),
        arg("side") = 0);

    // Degree
    if constexpr (n == 1) {
      basis.method("degree", [](const Basis& basis) { return basis.degree(); });
      basis.method("order", [](const Basis& basis) { return basis.order(); });
    } else {
      basis.method("degree", [](const Basis& basis, int i) { return basis.degree(i); });
    }

    basis.method(
        "degreeElevate", [](Basis& basis, const int i = 1) { basis.degreeElevate(i); }, arg("basis"), arg("i") = 1);
    basis.method(
        "degreeReduce", [](Basis& basis, const int i = 1) { basis.degreeReduce(i); }, arg("basis"), arg("i") = 1);
    basis.method(
        "degreeIncrease", [](Basis& basis, const int i = 1) { basis.degreeIncrease(i); }, arg("basis"), arg("i") = 1);
    basis.method(
        "degreeDecrease", [](Basis& basis, const int i = 1) { basis.degreeDecrease(i); }, arg("basis"), arg("i") = 1);
    basis.method(
        "elevateContinuity", [](Basis& basis, const int i = 1) { basis.elevateContinuity(i); }, arg("basis"),
        arg("i") = 1);
    basis.method(
        "reduceContinuity", [](Basis& basis, const int i = 1) { basis.reduceContinuity(i); }, arg("basis"),
        arg("i") = 1);

    // Refinement
    basis.method(
        "insertKnot", [](Basis& basis, double knot, int mult = 1) { basis.insertKnot(knot, mult); }, arg("basis"),
        arg("knot"), arg("mult") = 1);

    basis.method(
        "removeKnot", [](Basis& basis, double knot, int mult = 1) { basis.removeKnot(knot, mult); }, arg("basis"),
        arg("knot"), arg("mult") = 1);

    basis.method("insertKnots", [](Basis& basis, JuliaVector knots) {
      auto knots2 = std::vector<std::vector<double>>{std::vector<double>(knots.size())};
      std::copy(knots.begin(), knots.end(), knots2.front().begin());
      basis.insertKnots(knots2);
    });

    basis.method(
        "uniformRefine", [](Basis& basis, double numKnots = 1, int mul = 1) { basis.uniformRefine(numKnots, mul); },
        arg("basis"), arg("numKnots") = 1, arg("mul") = 1);

    basis.method(
        "uniformCoarsen", [](Basis& basis, double numKnots = 1) { basis.uniformCoarsen(numKnots); }, arg("basis"),
        arg("numKnots") = 1);

    basis.method(
        "uniformRefine_withCoefs",
        [](Basis& basis, JuliaMatrix coefs, int numKnots = 1, int mul = 1) {
          gsEigen::Map<const gsEigen::MatrixXd> coefsMatMap(coefs.data(), coefs.size());
          gismo::gsMatrix<double> coefsMat{coefsMatMap};

          basis.uniformRefine_withCoefs(coefsMat, numKnots, mul);
          assertSizeAndCopy(coefsMat, coefs);
        },
        arg("basis"), arg("coefs"), arg("numKnots") = 1, arg("mul") = 1);

    // Actives
    basis.method("active!",
                 [](Basis& basis, JuliaVector u, gismo::gsMatrix<int>& out) { basis.active_into(wrapVector(u), out); });

    basis.method("isActive", [](Basis& basis, int i, JuliaVector u) { return basis.isActive(i, wrapVector(u)); });

    // Eval
    basis.method("eval!",
                 [](Basis& basis, JuliaVector u, gismo::gsMatrix<>& out) { basis.eval_into(wrapVector(u), out); });
    basis.method("_eval", [](Basis& basis, JuliaVector u) { return basis.eval(wrapVector(u)); });

    basis.method("evalSingle!", [](Basis& basis, int i, JuliaVector u, gismo::gsMatrix<>& out) {
      basis.evalSingle_into(i, wrapVector(u), out);
    });

    basis.method("evalSingle", [](Basis& basis, int i, JuliaVector u) { return basis.evalSingle(i, wrapVector(u)); });

    basis.method("evalFunc!", [](Basis& basis, JuliaVector u, JuliaVector coefs, gismo::gsMatrix<>& out) {
      basis.evalFunc_into(wrapVector(u), wrapVector(coefs), out);
    });

    basis.method("evalFunc!", [](Basis& basis, JuliaVector u, JuliaMatrix coefs, gismo::gsMatrix<>& out) {
      basis.evalFunc_into(wrapVector(u), wrapMatrix(coefs), out);
    });

    basis.method("evalFunc", [](Basis& basis, JuliaVector u, JuliaVector coefs) {
      return basis.evalFunc(wrapVector(u), wrapVector(coefs));
    });

    basis.method("evalFunc", [](Basis& basis, JuliaVector u, JuliaMatrix coefs) {
      return basis.evalFunc(wrapVector(u), wrapMatrix(coefs));
    });

    // Eval Deriv
    basis.method("deriv!",
                 [](Basis& basis, JuliaVector u, gismo::gsMatrix<> out) { basis.deriv_into(wrapVector(u), out); });

    basis.method("derivSingle!", [](Basis& basis, int i, JuliaVector u, gismo::gsMatrix<>& out) {
      basis.derivSingle_into(i, wrapVector(u), out);
    });

    basis.method("derivSingle", [](Basis& basis, int i, JuliaVector u) { return basis.derivSingle(i, wrapVector(u)); });

    if constexpr (n == 1)
      basis.method("deriv!", [](Basis& basis, JuliaVector u, JuliaVector coefs, gismo::gsMatrix<>& out) {
        basis.deriv_into(wrapVector(u), wrapVector(coefs), out);
      });

    if constexpr (n == 1)
      basis.method("deriv!", [](Basis& basis, JuliaVector u, JuliaMatrix coefs, gismo::gsMatrix<>& out) {
        basis.deriv_into(wrapVector(u), wrapMatrix(coefs), out);
      });

    basis.method("derivFunc", [](Basis& basis, JuliaVector u, JuliaVector coefs) {
      return basis.derivFunc(wrapVector(u), wrapVector(coefs));
    });

    basis.method("derivFunc", [](Basis& basis, JuliaVector u, JuliaMatrix coefs) {
      return basis.derivFunc(wrapVector(u), wrapMatrix(coefs));
    });

    // Eval Second Derivative

    basis.method("deriv2!",
                 [](Basis& basis, JuliaVector u, gismo::gsMatrix<>& out) { basis.deriv2_into(wrapVector(u), out); });

    basis.method("deriv2Single!", [](Basis& basis, int i, JuliaVector u, gismo::gsMatrix<>& out) {
      basis.deriv2Single_into(i, wrapVector(u), out);
    });

    basis.method("deriv2Single",
                 [](Basis& basis, int i, JuliaVector u) { return basis.deriv2Single(i, wrapVector(u))(0, 0); });

    if constexpr (n == 1)
      basis.method("deriv2!", [](Basis& basis, JuliaVector u, JuliaVector coefs, gismo::gsMatrix<>& out) {
        basis.deriv2_into(wrapVector(u), wrapVector(coefs), out);
      });
    if constexpr (n == 1)
      basis.method("deriv2!", [](Basis& basis, JuliaVector u, JuliaMatrix coefs, gismo::gsMatrix<>& out) {
        basis.deriv2_into(wrapVector(u), wrapMatrix(coefs), out);
      });

    basis.method("deriv2Func", [](Basis& basis, JuliaVector u, JuliaVector coefs) {
      return basis.deriv2Func(wrapVector(u), wrapVector(coefs));
    });

    basis.method("deriv2Func", [](Basis& basis, JuliaVector u, JuliaMatrix coefs) {
      return basis.deriv2Func(wrapVector(u), wrapMatrix(coefs));
    });

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

} // namespace jlcxx

void regsiterBSplineBasis(jlcxx::Module& mod) {
  mod.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>>>("BSplineBasis")
      .apply<gismo::gsTensorBSplineBasis<1>>(WrapTensorBSplineBasis());

  mod.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>>>("TensorBSplineBasis")
      .apply<gismo::gsTensorBSplineBasis<2>, gismo::gsTensorBSplineBasis<3>>(WrapTensorBSplineBasis());
}
