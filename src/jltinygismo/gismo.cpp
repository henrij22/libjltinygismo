#include <jltinygismo/basis/registerbasis.hh>
#include <jltinygismo/geometry/registergeometries.hh>
#include <jltinygismo/utility/registerutilities.hh>

#include <jlcxx/module.hpp>

#include <gismo.h>
#include <gsCore/gsBasis.h>

JLCXX_MODULE define_julia_module(jlcxx::Module& mod) {
  registerGsMatrix(mod);

  registerKnotVector(mod);
  auto gsBasis = registerBasis(mod);
  regsiterBSplineBasis(mod, gsBasis);
  registerNurbsBasis(mod, gsBasis);

  auto gsGeometry = registerGeometry(mod);
  registerTensorBSpline(mod, gsGeometry);
  registerBSpline(mod, gsGeometry);
  registerNurbs(mod, gsGeometry);
  registerTensorNurbs(mod, gsGeometry);

  registerNurbsCreatorFunctions(mod);
  registerFileReaderFunctions(mod);
}