#include <jltinygismo/helper.hh>
#include <jlcxx/jlcxx.hpp>

#include <gismo.h>

template <typename Type>
void registerFileReaderFunctions_IMPL(jlcxx::Module& mod) {
  using jlcxx::arg;

  mod.method(
      "readFile",
      [](jlcxx::SingletonType<Type>, const std::string& filename) {
        if (!gismo::gsFileManager::fileExists(filename))
          throw std::runtime_error("readFile: no such file: " + filename);

        gismo::gsFileData<> data{filename};
        auto obj = data.getAnyFirst<Type>();
        if (!obj)
          throw std::runtime_error("readFile: " + filename + " contains no object of the requested type");
        return *obj;
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

  registerFileReaderFunctions_IMPL<gismo::gsTHBSplineBasis<1>>(mod);
  registerFileReaderFunctions_IMPL<gismo::gsTHBSplineBasis<2>>(mod);
  registerFileReaderFunctions_IMPL<gismo::gsTHBSplineBasis<3>>(mod);

  registerFileReaderFunctions_IMPL<gismo::gsHBSplineBasis<1>>(mod);
  registerFileReaderFunctions_IMPL<gismo::gsHBSplineBasis<2>>(mod);
  registerFileReaderFunctions_IMPL<gismo::gsHBSplineBasis<3>>(mod);

  registerFileReaderFunctions_IMPL<gismo::gsTHBSpline<1>>(mod);
  registerFileReaderFunctions_IMPL<gismo::gsTHBSpline<2>>(mod);
  registerFileReaderFunctions_IMPL<gismo::gsTHBSpline<3>>(mod);

  registerFileReaderFunctions_IMPL<gismo::gsHBSpline<1>>(mod);
  registerFileReaderFunctions_IMPL<gismo::gsHBSpline<2>>(mod);
  registerFileReaderFunctions_IMPL<gismo::gsHBSpline<3>>(mod);
}