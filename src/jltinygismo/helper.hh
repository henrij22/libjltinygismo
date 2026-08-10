#include <jlcxx/const_array.hpp>
#include <jlcxx/jlcxx.hpp>

#include <gsCore/gsDebug.h>
#include <gsNurbs/gsKnotVector.h>
#include <gsEigen/Eigen>
#include <string>
#include <algorithm>

template <typename Scalar>
inline void assertSizeAndCopy(const gismo::gsMatrix<Scalar>& fromMat, jlcxx::ArrayRef<Scalar, 2> out) {
  if (fromMat.rows() * fromMat.cols() != static_cast<int>(out.size()))
    throw std::runtime_error("Output size mismatch, should be " + std::to_string(fromMat.rows()) + "," +
                             std::to_string(fromMat.cols()));
  std::copy(fromMat.data(), fromMat.data() + fromMat.size(), out.data());
}

/// Copy \a mat into a freshly allocated, Julia-owned 2-D array.
///
/// jlcxx::make_julia_array only wraps the C++ buffer, so the Julia array stays valid just as
/// long as the gsMatrix behind it. Nearly every call site hands us a temporary --
/// `toMatrix(coefs(geo))` boxes a gsMatrix that becomes garbage the moment the call returns
/// -- and the view then silently reads freed memory. Copying is the only safe option here;
/// the bang methods remain available when the copy is not wanted.
template <typename Scalar>
inline jlcxx::ArrayRef<Scalar, 2> copyToJuliaArray(const gismo::gsMatrix<Scalar>& mat) {
  jl_value_t* arrayType = jlcxx::apply_array_type(jlcxx::julia_type<Scalar>(), 2);
  jl_array_t* array     = jl_alloc_array_2d(arrayType, mat.rows(), mat.cols());
  // No allocation happens between here and the return, so `array` needs no extra rooting.
  std::copy(mat.data(), mat.data() + mat.size(), jlcxx::jlcxx_array_data<Scalar>(array));
  return jlcxx::ArrayRef<Scalar, 2>(array);
}

/// Copy \a size elements starting at \a data into a freshly allocated, Julia-owned vector.
/// See copyToJuliaArray for why this copies.
template <typename Scalar>
inline jlcxx::ArrayRef<Scalar, 1> copyToJuliaVector(const Scalar* data, std::size_t size) {
  jl_value_t* arrayType = jlcxx::apply_array_type(jlcxx::julia_type<Scalar>(), 1);
  jl_array_t* array     = jl_alloc_array_1d(arrayType, size);
  std::copy(data, data + size, jlcxx::jlcxx_array_data<Scalar>(array));
  return jlcxx::ArrayRef<Scalar, 1>(array);
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

/// Translate a Julia-side parametric direction to the G+Smo convention.
///
/// On the Julia side directions are 1-based and `0` means "all directions"; G+Smo uses
/// 0-based directions with `-1` for "all". Anything outside that range is rejected here:
/// G+Smo indexes into a fixed-size direction array without bounds checking, so passing e.g.
/// `-1` straight through segfaults the process instead of raising an error.
inline short_t toGismoDir(int dir, short_t parDim, const char* fname) {
  if (dir < 0 || dir > static_cast<int>(parDim))
    throw std::runtime_error(std::string{fname} + ": direction must be 0 (= all directions) or in 1:" +
                             std::to_string(parDim) + ", got " + std::to_string(dir));
  return static_cast<short_t>(dir - 1);
}

/// As toGismoDir, but for operations that act on exactly one direction and have no
/// "all directions" mode -- inserting a knot at a given parameter, for instance. Passing
/// `0` here would reach G+Smo as `-1` and index out of bounds just like any other invalid
/// direction, so it is rejected too.
inline short_t toGismoDirRequired(int dir, short_t parDim, const char* fname) {
  if (dir < 1 || dir > static_cast<int>(parDim))
    throw std::runtime_error(std::string{fname} + ": direction must be in 1:" + std::to_string(parDim) + ", got " +
                             std::to_string(dir));
  return static_cast<short_t>(dir - 1);
}

template <typename Scalar>
inline void incrementByOne(gismo::gsMatrix<Scalar>& mat) {
  std::for_each(mat.reshaped().begin(), mat.reshaped().end(), [](auto& i) { i += 1; });
}
template <typename Scalar>
inline void incrementByOne(gismo::gsVector<Scalar>& vec) {
  std::for_each(vec.begin(), vec.end(), [](auto& i) { i += 1; });
}