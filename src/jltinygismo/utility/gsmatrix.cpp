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

    matrix.method("toMatrix",
                  [](Matrix& matrix) { return jlcxx::make_julia_array(matrix.data(), matrix.rows(), matrix.cols()); });

    matrix.method("toVector", [](Matrix& matrix) {
      if (matrix.cols() != 1)
        throw std::runtime_error("toVector: Matrix is not a vector!");
      return jlcxx::make_julia_array(matrix.data(), matrix.rows());
    });

    matrix.method("toValue", [](Matrix& matrix) {
      if (matrix.cols() != 1 || matrix.rows() != 1)
        throw std::runtime_error("toValue: Matrix has more than one entry!");
      return matrix(0, 0);
    });
  }
};

struct WrapVector
{
  template <typename VectorT>
  void operator()(VectorT&& matrix) {
    using jlcxx::arg;
    using jlcxx::julia_base_type;

    using Vector = typename VectorT::type;

    matrix.constructor([](int rows) { return new Vector(rows); });

    matrix.method("toVector", [](Vector& vector) { return jlcxx::make_julia_array(vector.data(), vector.rows()); });

    matrix.method("toValue", [](Vector& vector) {
      if (vector.rows() != 1)
        throw std::runtime_error("toValue: Vector has more than one entry!");
      return vector(0);
    });
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

void registerGsMatrix(jlcxx::Module& mod) {
  using jlcxx::arg;
  using jlcxx::julia_base_type;

  mod.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>>>("gsMatrix")
      .apply<gismo::gsMatrix<double>, gismo::gsMatrix<int>, gismo::gsMatrix<long>>(WrapMatrix());

  mod.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>>>("gsVector")
      .apply<gismo::gsVector<double>, gismo::gsVector<int>, gismo::gsVector<long>>(WrapMatrix());
}
