#include <jltinygismo/helper.hh>
#include <jlcxx/jlcxx.hpp>

#include <gismo.h>
#include <gsEigen/Eigen>

struct WrapMatrix
{
  template <typename MatrixT>
  void operator()(MatrixT&& matrix) {
    using jlcxx::arg;
    using jlcxx::julia_base_type;

    using Matrix = typename MatrixT::type;

    matrix.constructor([](int rows, int cols) { return new Matrix(rows, cols); });

    matrix.method(
        "toMatrix", [](Matrix& matrix) { return jlcxx::make_julia_array(matrix.data(), matrix.rows(), matrix.cols()); },
        arg("matrix"));

    matrix.method(
        "toVector",
        [](Matrix& matrix) {
          if (matrix.cols() != 1)
            throw std::runtime_error("toVector: Matrix is not a vector!");
          return jlcxx::make_julia_array(matrix.data(), matrix.rows());
        },
        arg("matrix"));

    matrix.method(
        "toValue",
        [](Matrix& matrix) {
          if (matrix.cols() != 1 || matrix.rows() != 1)
            throw std::runtime_error("toValue: Matrix has more than one entry!");
          return matrix(0, 0);
        },
        arg("matrix"));

    // no bounds-check access
    matrix.method(
        "_value", [](Matrix& matrix, int i, int j) { return matrix(i - 1, j - 1); }, arg("matrix"), arg("i"), arg("j"));

    // bounds-check access
    matrix.method(
        "value",
        [](Matrix& matrix, int i, int j) {
          if (matrix.rows() < i || matrix.cols() < j || i < 1 || j < 1)
            throw std::runtime_error("value: Index " + std::to_string(i) + ", " + std::to_string(j) + " out-of-bounds");
          return matrix(i - 1, j - 1);
        },
        arg("matrix"), arg("i"), arg("j"));

    matrix.method("size", &Matrix::size);
    matrix.method("rows", &Matrix::rows);
    matrix.method("cols", &Matrix::cols);
  }
};

struct WrapVector
{
  template <typename VectorT>
  void operator()(VectorT&& vector) {
    using jlcxx::arg;
    using jlcxx::julia_base_type;

    using Vector = typename VectorT::type;

    vector.constructor([](int rows) { return new Vector(rows); });

    vector.method(
        "toVector", [](Vector& vector) { return jlcxx::make_julia_array(vector.data(), vector.rows()); },
        arg("vector"));

    vector.method(
        "toValue",
        [](Vector& vector) {
          if (vector.rows() != 1)
            throw std::runtime_error("toValue: Vector has more than one entry!");
          return vector(0);
        },
        arg("vector"));

    // no bounds-check access
    vector.method("_value", [](Vector& vector, int i) { return vector(i - 1); }, arg("vector"), arg("i"));

    // bounds-check access
    vector.method(
        "value",
        [](Vector& vector, int i) {
          if (vector.rows() < i)
            throw std::runtime_error("value: Index " + std::to_string(i) + " out-of-bounds");
          return vector(i - 1);
        },
        arg("vector"), arg("i"));

    vector.method("size", &Vector::size);
    vector.method("rows", &Vector::rows);
    vector.method("cols", &Vector::cols);
  }
};

namespace jlcxx {

template <class T, int _Rows, int _Cols, int _Options>
struct BuildParameterList<gismo::gsMatrix<T, _Rows, _Cols, _Options>>
{
  // Using long to map to Julia Int64
  using type = ParameterList<T, std::integral_constant<long, _Rows>, std::integral_constant<long, _Cols>,
                             std::integral_constant<long, _Options>>;
};

template <class T, int _Rows, int _Options>
struct BuildParameterList<gismo::gsVector<T, _Rows, _Options>>
{
  // Using long to map to Julia Int64
  using type = ParameterList<T, std::integral_constant<long, _Rows>, std::integral_constant<long, _Options>>;
};

} // namespace jlcxx

void registerUtilities(jlcxx::Module& mod) {
  using jlcxx::arg;
  using jlcxx::julia_base_type;

  mod.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>>>("gsMatrix")
      .apply<gismo::gsMatrix<double>, gismo::gsMatrix<int>, gismo::gsMatrix<long>>(WrapMatrix());

  mod.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>>>("gsVector")
      .apply<gismo::gsVector<double>, gismo::gsVector<int>, gismo::gsVector<long>>(WrapVector());

  using KnotSpan = gismo::gsDomainIteratorWrapper<>;
  auto knotSpan  = mod.add_type<KnotSpan>("KnotSpan");

  knotSpan.method("centerPoint", &KnotSpan::centerPoint);
  knotSpan.method("upperCorner", &KnotSpan::upperCorner);
  knotSpan.method("lowerCorner", &KnotSpan::lowerCorner);
}
