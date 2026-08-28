#include <jltinygismo/helper.hh>

#include <jlcxx/array.hpp>
#include <jlcxx/jlcxx.hpp>

#include <gsHSplines/gsTHBSplineBasis.h>
#include <gsDomain/gsHDomainIterator.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsEigen/Eigen>
#include <string>
#include <vector>

#include "registerbasis.hh"

/// Wraps gsTHBSplineBasis<d,double,Trunc>.
///
/// gsHBSplineBasis is not a class of its own -- G+Smo declares it as
/// `using gsHBSplineBasis = gsTHBSplineBasis<d,T,false>` -- so one functor covers both flavours.
/// They remain distinct C++ types and so map to distinct Julia types.
struct WrapHierarchicalBasis
{
  template <typename T>
  struct h_traits;

  template <int n, bool Trunc>
  struct h_traits<gismo::gsTHBSplineBasis<n, double, Trunc>>
  {
    static constexpr int dim = n;
  };

  template <typename BasisT>
  void operator()(BasisT&& basis) {
    using jlcxx::arg;
    using JuliaMatrix = jlcxx::ArrayRef<double, 2>;
    using JuliaVector = jlcxx::ArrayRef<double, 1>;
    using JuliaBoxes  = jlcxx::ArrayRef<int64_t, 1>;

    using Basis       = typename BasisT::type;
    constexpr int d   = h_traits<Basis>::dim;
    using TensorBasis = gismo::gsTensorBSplineBasis<d>;

    // ---------------------------------------------------------------- constructors

    // A hierarchical basis with a single level, covering the tensor basis it is built from.
    basis.constructor([](const TensorBasis& tensorBasis) { return new Basis{tensorBasis}; });

    // Pre-refined at construction, boxes given as parametric corner coordinates (d x 2k).
    basis.constructor([](const TensorBasis& tensorBasis, JuliaMatrix boxes) {
      return new Basis{tensorBasis, gismo::gsMatrix<>{wrapMatrix(boxes)}};
    });

    // Pre-refined at construction, boxes given as element indices; see toGismoBoxes.
    basis.constructor([](const TensorBasis& tensorBasis, JuliaBoxes boxes) {
      return new Basis{tensorBasis, toGismoBoxes<d>(boxes, "THBSplineBasis")};
    });

    // ---------------------------------------------------------------- queries

    basis.method("size", [](const Basis& basis) { return basis.size(); }, arg("basis"));
    basis.method("numElements", [](const Basis& basis) { return basis.numElements(); }, arg("basis"));
    basis.method("degree", [](const Basis& basis, int i) { return basis.degree(i - 1); }, arg("basis"), arg("i"));
    basis.method("maxDegree", [](const Basis& basis) { return basis.maxDegree(); }, arg("basis"));
    basis.method("minDegree", [](const Basis& basis) { return basis.minDegree(); }, arg("basis"));

    // Levels are 1-based here, running 1:numLevels(basis). C++ calls that finest level
    // maxLevel(), which is 0-based.
    basis.method("numLevels", [](const Basis& basis) { return basis.numLevels(); }, arg("basis"));

    basis.method("treeSize", [](const Basis& basis) { return basis.treeSize(); }, arg("basis"));

    // knotSpans, registered on gsBasis, hands Julia clones of a gsDomainIterator, and
    // gsHDomainIterator is not safely copyable upstream: m_curElement holds std::vector iterators
    // into its own m_breaks, and the defaulted copy constructor copies both, so a clone
    // dereferences the *source's* storage. The elements are therefore read out eagerly below,
    // from the one live iterator, and returned as owned data.
    basis.method(
        "knotSpans",
        [](const Basis&) -> jlcxx::Array<gismo::gsDomainIteratorWrapper<>> {
          throw std::runtime_error(
              "knotSpans: not available for hierarchical bases -- G+Smo's gsHDomainIterator "
              "cannot be safely copied, so the element handles would dangle. Use elementBoxes, "
              "which returns the element corners as a matrix.");
        },
        arg("basis"));

    // A 2d x numElements matrix: the first d rows are each element's lower corner, the last d
    // its upper corner, in parametric coordinates.
    basis.method(
        "elementBoxes",
        [](const Basis& basis) {
          gismo::gsMatrix<double> boxes(2 * d, basis.numElements());
          gismo::gsHDomainIterator<double, d> it(basis.tree(), basis);
          for (index_t e = 0; e != basis.numElements(); ++e, it.next()) {
            boxes.col(e).head(d) = it.lowerCorner();
            boxes.col(e).tail(d) = it.upperCorner();
          }
          return boxes;
        },
        arg("basis"));

    basis.method("levelOf", [](const Basis& basis, int i) { return basis.levelOf(i - 1) + 1; }, arg("basis"), arg("i"));

    basis.method(
        "getLevelAtPoint",
        [](const Basis& basis, JuliaVector u) {
          return basis.getLevelAtPoint(gismo::gsMatrix<>{wrapVector(u)}) + 1;
        },
        arg("basis"), arg("u"));

    basis.method(
        "levelAtCorner", [](const Basis& basis, int c) { return basis.levelAtCorner(gismo::boxCorner{c}) + 1; },
        arg("basis"), arg("c"));

    // Returned by value: tensorLevel hands back a reference into the hierarchical basis, which
    // Julia must not own.
    basis.method(
        "tensorLevel",
        [](const Basis& basis, int level) -> TensorBasis {
          checkLevel(level, basis.numLevels(), "tensorLevel");
          return basis.tensorLevel(level - 1);
        },
        arg("basis"), arg("level"));

    // ---------------------------------------------------------------- refinement

    // Boxes as parametric corner coordinates: a d x 2k matrix, each consecutive pair of columns
    // being the lower and upper corner of one box. refExt widens every box by that many cells.
    basis.method(
        "refine!",
        [](Basis& basis, JuliaMatrix boxes, int refExt = 0) {
          basis.refine(gismo::gsMatrix<>{wrapMatrix(boxes)}, refExt);
        },
        arg("basis"), arg("boxes"), arg("refExt") = 0);

    basis.method(
        "unrefine!",
        [](Basis& basis, JuliaMatrix boxes, int refExt = 0) {
          basis.unrefine(gismo::gsMatrix<>{wrapMatrix(boxes)}, refExt);
        },
        arg("basis"), arg("boxes"), arg("refExt") = 0);

    basis.method(
        "refineElements!",
        [](Basis& basis, JuliaBoxes boxes) { basis.refineElements(toGismoBoxes<d>(boxes, "refineElements!")); },
        arg("basis"), arg("boxes"));

    basis.method(
        "unrefineElements!",
        [](Basis& basis, JuliaBoxes boxes) { basis.unrefineElements(toGismoBoxes<d>(boxes, "unrefineElements!")); },
        arg("basis"), arg("boxes"));

    // The _withCoefs variants refine in place and return the matching coefficients. Those cannot
    // be updated in place: refinement changes the number of control points, so the result has a
    // different number of rows than the input.
    basis.method(
        "refineElements_withCoefs!",
        [](Basis& basis, JuliaMatrix coefs, JuliaBoxes boxes) {
          gismo::gsMatrix<double> coefsMat{wrapMatrix(coefs)};
          checkCoefRows(coefsMat.rows(), basis.size(), "refineElements_withCoefs!");
          basis.refineElements_withCoefs(coefsMat, toGismoBoxes<d>(boxes, "refineElements_withCoefs!"));
          return coefsMat;
        },
        arg("basis"), arg("coefs"), arg("boxes"));

    basis.method(
        "unrefineElements_withCoefs!",
        [](Basis& basis, JuliaMatrix coefs, JuliaBoxes boxes) {
          gismo::gsMatrix<double> coefsMat{wrapMatrix(coefs)};
          checkCoefRows(coefsMat.rows(), basis.size(), "unrefineElements_withCoefs!");
          basis.unrefineElements_withCoefs(coefsMat, toGismoBoxes<d>(boxes, "unrefineElements_withCoefs!"));
          return coefsMat;
        },
        arg("basis"), arg("coefs"), arg("boxes"));

    basis.method(
        "refine_withCoefs!",
        [](Basis& basis, JuliaMatrix coefs, JuliaMatrix boxes) {
          gismo::gsMatrix<double> coefsMat{wrapMatrix(coefs)};
          checkCoefRows(coefsMat.rows(), basis.size(), "refine_withCoefs!");
          basis.refine_withCoefs(coefsMat, gismo::gsMatrix<>{wrapMatrix(boxes)});
          return coefsMat;
        },
        arg("basis"), arg("coefs"), arg("boxes"));

    basis.method(
        "refineSide!",
        [](Basis& basis, int side, int level) {
          checkLevel(level, basis.numLevels() + 1, "refineSide!");
          basis.refineSide(gismo::boxSide{side}, level - 1);
        },
        arg("basis"), arg("side"), arg("level"));

    basis.method(
        "refineBasisFunction!", [](Basis& basis, int i) { basis.refineBasisFunction(i - 1); }, arg("basis"), arg("i"));

    basis.method(
        "increaseMultiplicity!",
        [](Basis& basis, int level, int dir, double knotValue, int mult = 1) {
          checkLevel(level, basis.numLevels(), "increaseMultiplicity!");
          basis.increaseMultiplicity(level - 1, toGismoDirRequired(dir, d, "increaseMultiplicity!"), knotValue, mult);
        },
        arg("basis"), arg("level"), arg("dir"), arg("knotValue"), arg("mult") = 1);
  }
};

namespace jlcxx {

template <int d, typename T, bool Trunc>
struct BuildParameterList<gismo::gsTHBSplineBasis<d, T, Trunc>>
{
  // Using long to map to Julia Int64
  using type = ParameterList<std::integral_constant<int64_t, d>, T>;
};

// Flat hierarchy, as for gsTensorBSplineBasis: gsHTensorBasis is not itself a Julia type. Must
// come before the add_type calls below -- a late SuperType specialization fails at run time, not
// at compile time.
template <int d, bool Trunc>
struct SuperType<gismo::gsTHBSplineBasis<d, double, Trunc>>
{
  using type = gismo::gsBasis<>;
};

} // namespace jlcxx

void registerHierarchicalBasis(jlcxx::Module& mod, jlcxx::TypeWrapper<gismo::gsBasis<double>>& gsBasis) {
  mod.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>>>("THBSplineBasis", gsBasis.dt())
      .apply<gismo::gsTHBSplineBasis<1>, gismo::gsTHBSplineBasis<2>, gismo::gsTHBSplineBasis<3>>(
          WrapHierarchicalBasis());

  mod.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>>>("HBSplineBasis", gsBasis.dt())
      .apply<gismo::gsHBSplineBasis<1>, gismo::gsHBSplineBasis<2>, gismo::gsHBSplineBasis<3>>(WrapHierarchicalBasis());
}
