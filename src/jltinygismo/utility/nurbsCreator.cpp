#include <jlcxx/jlcxx.hpp>
#include <gismo.h>

// register as free functions
void registerNurbsCreatorFunctions(jlcxx::Module& mod) {
  using jlcxx::arg;

  mod.method(
      "createBSplineRectangle",
      [](double low_x = 0, double low_y = 0, double upp_x = 1, double upp_y = 1, double turndeg = 0) {
        return *gismo::gsNurbsCreator<>::BSplineRectangle(low_x, low_y, upp_x, upp_y, turndeg);
      },
      arg("low_x") = 0, arg("low_y") = 0, arg("upp_x") = 1, arg("upp_y") = 1, arg("turndeg") = 0);

  // ------------------------------------------------------------
  // BSpline Trapezium (parameter version)
  mod.method(
      "createBSplineTrapezium",
      [](double Lbot = 1, double Ltop = 0.5, double H = 1, double d = 0, double turndeg = 0) {
        return *gismo::gsNurbsCreator<>::BSplineTrapezium(Lbot, Ltop, H, d, turndeg);
      },
      arg("Lbot") = 1, arg("Ltop") = 0.5, arg("H") = 1, arg("d") = 0, arg("turndeg") = 0);

  // BSpline Trapezium (corner points)
  mod.method(
      "createBSplineTrapezium",
      [](double Ax, double Ay, double Bx, double By, double Cx, double Cy, double Dx, double Dy, double turndeg = 0) {
        return *gismo::gsNurbsCreator<>::BSplineTrapezium(Ax, Ay, Bx, By, Cx, Cy, Dx, Dy, turndeg);
      },
      arg("Ax"), arg("Ay"), arg("Bx"), arg("By"), arg("Cx"), arg("Cy"), arg("Dx"), arg("Dy"), arg("turndeg") = 0);

  // ------------------------------------------------------------
  // // Nurbs Arc Trapezium (parameter version)
  // mod.method(
  //     "createNurbsArcTrapezium",
  //     [](double Lbot = 1, double Ltop = 0.5, double H = 1, double d = 0, double turndeg = 0) {
  //       return *gismo::gsNurbsCreator<>::NurbsArcTrapezium(Lbot, Ltop, H, d, turndeg);
  //     },
  //     arg("Lbot") = 1, arg("Ltop") = 0.5, arg("H") = 1, arg("d") = 0, arg("turndeg") = 0);

  // // Nurbs Arc Trapezium (corner points)
  // mod.method(
  //     "createNurbsArcTrapezium",
  //     [](double Ax, double Ay, double Bx, double By, double Cx, double Cy, double Dx, double Dy, double turndeg = 0) {
  //       return *gismo::gsNurbsCreator<>::NurbsArcTrapezium(Ax, Ay, Bx, By, Cx, Cy, Dx, Dy, turndeg);
  //     },
  //     arg("Ax"), arg("Ay"), arg("Bx"), arg("By"), arg("Cx"), arg("Cy"), arg("Dx"), arg("Dy"), arg("turndeg") = 0);

  // ------------------------------------------------------------
  // BSpline shapes
  mod.method(
      "createBSplineAmoeba",
      [](double r = 1, double x = 0, double y = 0) { return *gismo::gsNurbsCreator<>::BSplineAmoeba(r, x, y); },
      arg("r") = 1, arg("x") = 0, arg("y") = 0);

  mod.method(
      "createBSplineAmoebaBig",
      [](double r = 1, double x = 0, double y = 0) { return *gismo::gsNurbsCreator<>::BSplineAmoebaBig(r, x, y); },
      arg("r") = 1, arg("x") = 0, arg("y") = 0);

  mod.method(
      "createBSplineAustria",
      [](double r = 1, double x = 0, double y = 0) { return *gismo::gsNurbsCreator<>::BSplineAustria(r, x, y); },
      arg("r") = 1, arg("x") = 0, arg("y") = 0);

  mod.method(
      "createBSplineFish",
      [](double r = 1, double x = 0, double y = 0) { return *gismo::gsNurbsCreator<>::BSplineFish(r, x, y); },
      arg("r") = 1, arg("x") = 0, arg("y") = 0);

  mod.method(
      "createBSplineAmoeba3degree",
      [](double r = 1, double x = 0, double y = 0) { return *gismo::gsNurbsCreator<>::BSplineAmoeba3degree(r, x, y); },
      arg("r") = 1, arg("x") = 0, arg("y") = 0);

  // ------------------------------------------------------------
  // NURBS quarter plate with hole
  mod.method("createNurbsQrtPlateWHoleC0", []() { return *gismo::gsNurbsCreator<>::NurbsQrtPlateWHoleC0(); });

  // ------------------------------------------------------------
  // Triangle
  mod.method(
      "createBSplineTriangle",
      [](double H = 1, double W = 1) { return *gismo::gsNurbsCreator<>::BSplineTriangle(H, W); }, arg("H") = 1,
      arg("W") = 1);
}