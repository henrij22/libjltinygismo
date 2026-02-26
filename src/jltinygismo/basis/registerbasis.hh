#include <jlcxx/jlcxx.hpp>
#include <gismo.h>

jlcxx::TypeWrapper<gismo::gsBasis<double>> registerBasis(jlcxx::Module& mod);
void registerKnotVector(jlcxx::Module& mod);
void regsiterBSplineBasis(jlcxx::Module& mod, jlcxx::TypeWrapper<gismo::gsBasis<double>>& gsBasis);
void registerNurbsBasis(jlcxx::Module& mod, jlcxx::TypeWrapper<gismo::gsBasis<double>>& gsBasis);
