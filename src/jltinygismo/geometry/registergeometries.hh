#include <jlcxx/jlcxx.hpp>
#include <gismo.h>

jlcxx::TypeWrapper<gismo::gsGeometry<double>> registerGeometry(jlcxx::Module& mod);
void registerBSpline(jlcxx::Module& mod, jlcxx::TypeWrapper<gismo::gsGeometry<double>>& gsGeometry);
void registerTensorBSpline(jlcxx::Module& mod, jlcxx::TypeWrapper<gismo::gsGeometry<double>>& gsGeometry);