#include <jltinygismo/helper.hh>

#include <jlcxx/module.hpp>

#include <gismo.h>

#include "registerbasis.hh"

jlcxx::TypeWrapper<gismo::gsBasis<double>> registerBasis(jlcxx::Module& mod) {
  using jlcxx::arg;

  using JuliaMatrix = jlcxx::ArrayRef<double, 2>;
  using JuliaVector = jlcxx::ArrayRef<double, 1>;

  using Basis = gismo::gsBasis<>;
  auto basis  = mod.add_type<Basis>("gsBasis");

  basis.method("elementIndex", [](Basis& basis, JuliaVector u) { return basis.elementIndex(wrapVector(u)) + 1; });
  basis.method("elementIndex", [](Basis& basis, JuliaMatrix u) { return basis.elementIndex(wrapMatrix(u)) + 1; });

  basis.method("elementInSupportOf", [](Basis& basis, int j) { return basis.elementInSupportOf(j - 1); }, arg("j"));

  // Actives
  basis.method(
      "active!",
      [](Basis& basis, JuliaVector u, gismo::gsMatrix<int>& out) {
        basis.active_into(wrapVector(u), out);
        incrementByOne(out);
      },
      arg("basis"), arg("u"), arg("out"));

  basis.method(
      "active!",
      [](Basis& basis, JuliaMatrix u, gismo::gsMatrix<int>& out) {
        basis.active_into(wrapMatrix(u), out);
        incrementByOne(out);
      },
      arg("basis"), arg("u"), arg("out"));

  basis.method(
      "isActive", [](Basis& basis, int i, JuliaVector u) { return basis.isActive(i - 1, wrapVector(u)); }, arg("basis"),
      arg("i"), arg("u"));

  // Degree
  basis.method(
      "degreeElevate", [](Basis& basis, const int i = 1, int dir = -1) { basis.degreeElevate(i); }, arg("basis"),
      arg("i") = 1, arg("dir") = -1);
  basis.method(
      "degreeReduce", [](Basis& basis, const int i = 1, int dir = -1) { basis.degreeReduce(i); }, arg("basis"),
      arg("i") = 1, arg("dir") = -1);
  basis.method(
      "degreeIncrease", [](Basis& basis, const int i = 1, int dir = -1) { basis.degreeIncrease(i); }, arg("basis"),
      arg("i") = 1, arg("dir") = -1);
  basis.method(
      "degreeDecrease", [](Basis& basis, const int i = 1, int dir = -1) { basis.degreeDecrease(i); }, arg("basis"),
      arg("i") = 1, arg("dir") = -1);
  basis.method(
      "elevateContinuity", [](Basis& basis, const int i = 1, int dir = -1) { basis.elevateContinuity(i); },
      arg("basis"), arg("i") = 1, arg("dir") = -1);
  basis.method(
      "reduceContinuity", [](Basis& basis, const int i = 1, int dir = -1) { basis.reduceContinuity(i); }, arg("basis"),
      arg("i") = 1, arg("dir") = -1);

  basis.method("setDegree", &Basis::setDegree, arg("i"));
  basis.method("setDegreePreservingMultiplicity", &Basis::setDegreePreservingMultiplicity, arg("i"));

  basis.method("reverse", &Basis::reverse);
  // Eval
  basis.method(
      "eval!", [](Basis& basis, JuliaVector u, gismo::gsMatrix<>& out) { basis.eval_into(wrapVector(u), out); },
      arg("basis"), arg("u"), arg("out"));
  basis.method("_eval", [](Basis& basis, JuliaVector u) { return basis.eval(wrapVector(u)); }, arg("basis"), arg("u"));

  basis.method(
      "evalSingle!",
      [](Basis& basis, int i, JuliaVector u, gismo::gsMatrix<>& out) {
        basis.evalSingle_into(i - 1, wrapVector(u), out);
      },
      arg("basis"), arg("i"), arg("u"), arg("out"));

  basis.method(
      "evalSingle", [](Basis& basis, int i, JuliaVector u) { return basis.evalSingle(i - 1, wrapVector(u)); },
      arg("basis"), arg("i"), arg("u"));

  basis.method(
      "evalFunc!",
      [](Basis& basis, JuliaVector u, JuliaVector coefs, gismo::gsMatrix<>& out) {
        basis.evalFunc_into(wrapVector(u), wrapVector(coefs), out);
      },
      arg("basis"), arg("u"), arg("coefs"), arg("out"));

  basis.method(
      "evalFunc!",
      [](Basis& basis, JuliaVector u, JuliaMatrix coefs, gismo::gsMatrix<>& out) {
        basis.evalFunc_into(wrapVector(u), wrapMatrix(coefs), out);
      },
      arg("basis"), arg("u"), arg("coefs"), arg("out"));

  basis.method(
      "evalFunc",
      [](Basis& basis, JuliaVector u, JuliaVector coefs) { return basis.evalFunc(wrapVector(u), wrapVector(coefs)); },
      arg("basis"), arg("u"), arg("coefs"));

  basis.method(
      "evalFunc",
      [](Basis& basis, JuliaVector u, JuliaMatrix coefs) { return basis.evalFunc(wrapVector(u), wrapMatrix(coefs)); },
      arg("basis"), arg("u"), arg("coefs"));

  // Eval Deriv
  basis.method(
      "deriv!", [](Basis& basis, JuliaVector u, gismo::gsMatrix<> out) { basis.deriv_into(wrapVector(u), out); },
      arg("basis"), arg("u"), arg("out"));

  basis.method(
      "derivSingle!",
      [](Basis& basis, int i, JuliaVector u, gismo::gsMatrix<>& out) {
        basis.derivSingle_into(i - 1, wrapVector(u), out);
      },
      arg("basis"), arg("i"), arg("u"), arg("out"));

  basis.method(
      "derivSingle", [](Basis& basis, int i, JuliaVector u) { return basis.derivSingle(i - 1, wrapVector(u)); },
      arg("basis"), arg("i"), arg("u"));

  basis.method(
      "derivFunc",
      [](Basis& basis, JuliaVector u, JuliaVector coefs) { return basis.derivFunc(wrapVector(u), wrapVector(coefs)); },
      arg("basis"), arg("u"), arg("coefs"));

  basis.method(
      "derivFunc",
      [](Basis& basis, JuliaVector u, JuliaMatrix coefs) { return basis.derivFunc(wrapVector(u), wrapMatrix(coefs)); },
      arg("basis"), arg("u"), arg("coefs"));

  // Eval Second Derivative

  basis.method(
      "deriv2!", [](Basis& basis, JuliaVector u, gismo::gsMatrix<>& out) { basis.deriv2_into(wrapVector(u), out); },
      arg("basis"), arg("u"), arg("out"));

  basis.method(
      "deriv2Single!",
      [](Basis& basis, int i, JuliaVector u, gismo::gsMatrix<>& out) {
        basis.deriv2Single_into(i - 1, wrapVector(u), out);
      },
      arg("basis"), arg("i"), arg("u"), arg("out"));

  basis.method(
      "deriv2Single", [](Basis& basis, int i, JuliaVector u) { return basis.deriv2Single(i - 1, wrapVector(u))(0, 0); },
      arg("basis"), arg("i"), arg("u"));

  basis.method(
      "deriv2Func",
      [](Basis& basis, JuliaVector u, JuliaVector coefs) { return basis.deriv2Func(wrapVector(u), wrapVector(coefs)); },
      arg("basis"), arg("u"), arg("coefs"));

  basis.method(
      "deriv2Func",
      [](Basis& basis, JuliaVector u, JuliaMatrix coefs) { return basis.deriv2Func(wrapVector(u), wrapMatrix(coefs)); },
      arg("basis"), arg("u"), arg("coefs"));

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

  // elements
  basis.method("knotSpans", [](Basis& basis) {
    jlcxx::Array<gismo::gsDomainIteratorWrapper<>> elements{};
    for (auto it = basis.domain()->beginAll(); it != basis.domain()->endAll(); ++it) {
      elements.push_back(it); // Construct wrapper from iterator
    }
    return elements;
  });

  return basis;
}