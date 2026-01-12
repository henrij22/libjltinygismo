#include <jltinygismo/basis/registerbasis.hh>
#include <jltinygismo/geometry/registergeometries.hh>
#include <jltinygismo/utility/registerutilities.hh>

#include <jlcxx/module.hpp>

#include <gismo.h>
#include <gsCore/gsBasis.h>

JLCXX_MODULE define_julia_module(jlcxx::Module& mod) {
  using jlcxx::arg;
  using jlcxx::julia_base_type;

  std::vector<double> knots{0, 0, 0, 0.5, 1, 1, 1};
  gismo::gsKnotVector<> knotVector{knots};

  auto basis = gismo::gsBSplineBasis<>{knotVector};

  registerGsMatrix(mod);

  registerKnotVector(mod);
  regsiterBSplineBasis(mod);

  auto gsGeometry = registerGeometry(mod);
  registerTensorBSpline(mod, gsGeometry);
  registerBSpline(mod, gsGeometry);

  registerNurbsCreatorFunctions(mod);

  // gismo::gsNurbsCreator<>::BSplineRe

  // auto geo = gismo::gsBSpline<>{};
  // auto geo2 = gismo::gsTensorBSpline<2>{};
}