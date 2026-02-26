#include <jltinygismo/helper.hh>
#include <jlcxx/jlcxx.hpp>

#include <gismo.h>

template <typename Type>
void registerFileReaderFunctions_IMPL(jlcxx::Module& mod) {
  using jlcxx::arg;

  mod.method(
      "readFile",
      [](jlcxx::SingletonType<Type>, const std::string& filename) {
        gismo::gsFileData<> data{filename};
        return *data.getAnyFirst<Type>();
      },
      arg("type"), arg("filename"));
}

void registerFileReaderFunctions(jlcxx::Module& mod) {
  using jlcxx::arg;

  registerFileReaderFunctions_IMPL<gismo::gsTensorBSplineBasis<1>>(mod);
  registerFileReaderFunctions_IMPL<gismo::gsTensorBSplineBasis<2>>(mod);
  registerFileReaderFunctions_IMPL<gismo::gsTensorBSplineBasis<3>>(mod);

  registerFileReaderFunctions_IMPL<gismo::gsTensorNurbsBasis<1>>(mod);
  registerFileReaderFunctions_IMPL<gismo::gsTensorNurbsBasis<2>>(mod);
  registerFileReaderFunctions_IMPL<gismo::gsTensorNurbsBasis<3>>(mod);

  registerFileReaderFunctions_IMPL<gismo::gsBSpline<>>(mod);
  registerFileReaderFunctions_IMPL<gismo::gsTensorBSpline<1>>(mod);
  registerFileReaderFunctions_IMPL<gismo::gsTensorBSpline<2>>(mod);
  registerFileReaderFunctions_IMPL<gismo::gsTensorBSpline<3>>(mod);

  registerFileReaderFunctions_IMPL<gismo::gsNurbs<>>(mod);
  registerFileReaderFunctions_IMPL<gismo::gsTensorNurbs<2>>(mod);
  registerFileReaderFunctions_IMPL<gismo::gsTensorNurbs<3>>(mod);
}