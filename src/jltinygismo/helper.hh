#include <jlcxx/const_array.hpp>
#include <jlcxx/jlcxx.hpp>

#include <gsCore/gsDebug.h>
#include <gsNurbs/gsKnotVector.h>
#include <gsEigen/Eigen>
#include <string>
#include <algorithm>
#include <vector>

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
/// Reject a level outside the levels a hierarchical basis has. Levels are 1-based here, so they
/// run `1:numLevels`. G+Smo indexes its level containers unchecked, so an out-of-range level
/// would reach it as an out-of-bounds read rather than an error.
inline void checkLevel(int level, index_t numLevels, const char* fname) {
  if (level < 1 || level > static_cast<int>(numLevels))
    throw std::runtime_error(std::string{fname} + ": level must be in 1:" + std::to_string(numLevels) + ", got " +
                             std::to_string(level));
}

/// Guard the coefficient array of a *_withCoefs refinement against a basis it does not belong to.
inline void checkCoefRows(index_t rows, index_t basisSize, const char* fname) {
  if (rows != basisSize)
    throw std::runtime_error(std::string{fname} + ": coefs has " + std::to_string(rows) +
                             " rows but the basis has " + std::to_string(basisSize) + " functions");
}

/// Translate Julia-side element boxes into the flat index vector G+Smo's refineElements expects.
///
/// Both are flat arrays of `2d+1` entries per box, `[level, lower..., upper...]`. G+Smo uses a
/// 0-based level and 0-based span indices on that level's grid, upper corner *exclusive*;
/// Julia-side both are 1-based and the corners *inclusive*, so a box reads as the cell range
/// `lower:upper`. The upper corner is therefore the one entry not shifted -- an inclusive
/// 1-based bound is already the exclusive 0-based one.
template <int d>
inline std::vector<index_t> toGismoBoxes(jlcxx::ArrayRef<int64_t, 1> boxes, const char* fname) {
  constexpr std::size_t stride = 2 * d + 1;
  const std::string prefix{fname};

  if (boxes.size() % stride != 0)
    throw std::runtime_error(prefix + ": box array length must be a multiple of " + std::to_string(stride) + " for a " +
                             std::to_string(d) + "-dimensional basis, got " + std::to_string(boxes.size()));

  std::vector<index_t> out;
  out.reserve(boxes.size());

  for (std::size_t box = 0; box != boxes.size() / stride; ++box) {
    const int64_t* entry = boxes.data() + box * stride;

    if (entry[0] < 1)
      throw std::runtime_error(prefix + ": level must be >= 1, got " + std::to_string(entry[0]));
    out.push_back(static_cast<index_t>(entry[0] - 1));

    for (std::size_t i = 0; i != d; ++i) {
      const int64_t lower = entry[1 + i];
      const int64_t upper = entry[1 + d + i];
      if (lower < 1)
        throw std::runtime_error(prefix + ": lower corner must be >= 1, got " + std::to_string(lower));
      if (upper < lower)
        throw std::runtime_error(prefix + ": upper corner (" + std::to_string(upper) +
                                 ") must not be below the lower corner (" + std::to_string(lower) + ")");
    }

    for (std::size_t i = 0; i != d; ++i)
      out.push_back(static_cast<index_t>(entry[1 + i] - 1));
    for (std::size_t i = 0; i != d; ++i)
      out.push_back(static_cast<index_t>(entry[1 + d + i]));
  }

  return out;
}
