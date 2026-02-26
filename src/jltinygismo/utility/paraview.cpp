#include <jltinygismo/helper.hh>
#include <jlcxx/jlcxx.hpp>
#include <gismo.h>

// register paraview export functions
void registerParaviewFunctions(jlcxx::Module& mod) {
  using jlcxx::arg;
  using JuliaMatrix = jlcxx::ArrayRef<double, 2>;

  // ------------------------------------------------------------
  // gsGeometry exports
  mod.method(
      "writeParaview",
      [](const gismo::gsGeometry<double>& geo, const std::string& fn, unsigned npts = 1000, bool mesh = false,
         bool ctrlNet = false) { gismo::gsWriteParaview(geo, fn, npts, mesh, ctrlNet); },
      arg("geo"), arg("fn"), arg("npts") = 1000, arg("mesh") = false, arg("ctrlNet") = false);

  // ------------------------------------------------------------
  // gsBasis exports
  mod.method(
      "writeParaview",
      [](const gismo::gsBasis<double>& basis, const std::string& fn, unsigned npts = 1000, bool mesh = false) {
        gismo::gsWriteParaview(basis, fn, npts, mesh);
      },
      arg("basis"), arg("fn"), arg("npts") = 1000, arg("mesh") = false);

  // gsBasis export with indices
  mod.method(
      "writeParaview",
      [](const gismo::gsBasis<double>& basis, const std::vector<index_t>& indices, const std::string& fn,
         unsigned npts = 1000, bool mesh = false) { gismo::gsWriteParaview(basis, indices, fn, npts, mesh); },
      arg("basis"), arg("indices"), arg("fn"), arg("npts") = 1000, arg("mesh") = false);

  // Single basis function export
  mod.method(
      "writeParaviewBasisFnct",
      [](int i, const gismo::gsBasis<double>& basis, const std::string& fn, unsigned npts = 1000) {
        gismo::gsWriteParaview_basisFnct(i, basis, fn, npts);
      },
      arg("i"), arg("basis"), arg("fn"), arg("npts") = 1000);

  // ------------------------------------------------------------
  // Point exports
  mod.method(
      "writeParaviewPoints",
      [](JuliaMatrix points, const std::string& fn) {
        gismo::gsWriteParaviewPoints(gismo::gsMatrix<>(wrapMatrix(points)), fn);
      },
      arg("points"), arg("fn"));

  mod.method(
      "writeParaviewPoints",
      [](JuliaMatrix X, JuliaMatrix Y, const std::string& fn) {
        gismo::gsWriteParaviewPoints(gismo::gsMatrix<>(wrapMatrix(X)), gismo::gsMatrix<>(wrapMatrix(Y)), fn);
      },
      arg("X"), arg("Y"), arg("fn"));

  mod.method(
      "writeParaviewPoints",
      [](JuliaMatrix X, JuliaMatrix Y, JuliaMatrix Z, const std::string& fn) {
        gismo::gsWriteParaviewPoints(gismo::gsMatrix<>(wrapMatrix(X)), gismo::gsMatrix<>(wrapMatrix(Y)),
                                     gismo::gsMatrix<>(wrapMatrix(Z)), fn);
      },
      arg("X"), arg("Y"), arg("Z"), arg("fn"));

  mod.method(
      "writeParaviewPoints",
      [](JuliaMatrix X, JuliaMatrix Y, JuliaMatrix Z, JuliaMatrix V, const std::string& fn) {
        gismo::gsWriteParaviewPoints(gismo::gsMatrix<>(wrapMatrix(X)), gismo::gsMatrix<>(wrapMatrix(Y)),
                                     gismo::gsMatrix<>(wrapMatrix(Z)), gismo::gsMatrix<>(wrapMatrix(V)), fn);
      },
      arg("X"), arg("Y"), arg("Z"), arg("V"), arg("fn"));
}
