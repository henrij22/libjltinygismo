#include <jltinygismo/basis/registerbasis.hh>
#include <jltinygismo/geometry/registergeometries.hh>
#include <jltinygismo/utility/registerutilities.hh>

#include <jlcxx/module.hpp>

#include <gismo.h>
#include <gsCore/gsBasis.h>

JLCXX_MODULE define_julia_module(jlcxx::Module& mod) {
  registerGsMatrix(mod);

  registerKnotVector(mod);
  regsiterBSplineBasis(mod);

  auto gsGeometry = registerGeometry(mod);
  registerTensorBSpline(mod, gsGeometry);
  registerBSpline(mod, gsGeometry);

  registerNurbsCreatorFunctions(mod);
}