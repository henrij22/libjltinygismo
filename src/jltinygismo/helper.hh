#include <jlcxx/const_array.hpp>
#include <jlcxx/jlcxx.hpp>

#include <gsCore/gsDebug.h>
#include <gsNurbs/gsKnotVector.h>
#include <gsEigen/Eigen>
#include <string>

template <typename Scalar>
inline void assertSizeAndCopy(const gismo::gsMatrix<Scalar>& fromMat, jlcxx::ArrayRef<Scalar, 2> out) {
  if (fromMat.rows() * fromMat.cols() != static_cast<int>(out.size()))
    throw std::runtime_error("Output size mismatch, should be " + std::to_string(fromMat.rows()) + "," + std::to_string(fromMat.cols()));
  std::copy(fromMat.data(), fromMat.data() + fromMat.size(), out.data());
}

template <typename Scalar>
inline auto rowsAndCols(jlcxx::ArrayRef<Scalar, 2> mat) {
  return std::array<size_t, 2>{mat.m_array->dimsize[0], mat.m_array->dimsize[1]};
}

template <typename Scalar>
inline auto wrapVector(jlcxx::ArrayRef<Scalar, 1> mat) {
  return gsEigen::Map<const gsEigen::VectorX<Scalar>>(mat.data(), mat.size());
}

template <typename Scalar>
inline auto wrapMatrix(jlcxx::ArrayRef<Scalar, 2> mat) {
  auto [rows, cols] = rowsAndCols(mat);
  return gsEigen::Map<const gsEigen::MatrixX<Scalar>>(mat.data(), rows, cols);
}
