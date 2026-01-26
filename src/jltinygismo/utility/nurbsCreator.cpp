#include <jlcxx/jlcxx.hpp>
#include <gismo.h>

// register as free functions
void registerNurbsCreatorFunctions(jlcxx::Module& mod) {
  using jlcxx::arg;

  mod.method(
      "createBSplineUnitInterval", [](short_t deg) { return *gismo::gsNurbsCreator<>::BSplineUnitInterval(deg); },
      arg("deg"));

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
  // Nurbs Arc Trapezium (parameter version)
  mod.method(
      "createNurbsArcTrapezium",
      [](double Lbot = 1, double Ltop = 0.5, double H = 1, double d = 0, double turndeg = 0) {
        return *gismo::gsNurbsCreator<>::NurbsArcTrapezium(Lbot, Ltop, H, d, turndeg);
      },
      arg("Lbot") = 1, arg("Ltop") = 0.5, arg("H") = 1, arg("d") = 0, arg("turndeg") = 0);

  // Nurbs Arc Trapezium (corner points)
  mod.method(
      "createNurbsArcTrapezium",
      [](double Ax, double Ay, double Bx, double By, double Cx, double Cy, double Dx, double Dy, double turndeg = 0) {
        return *gismo::gsNurbsCreator<>::NurbsArcTrapezium(Ax, Ay, Bx, By, Cx, Cy, Dx, Dy, turndeg);
      },
      arg("Ax"), arg("Ay"), arg("Bx"), arg("By"), arg("Cx"), arg("Cy"), arg("Dx"), arg("Dy"), arg("turndeg") = 0);

  // ------------------------------------------------------------
  // Basic B-spline and NURBS primitives
  mod.method(
      "createBSplineSquare",
      [](double r = 1, double x = 0, double y = 0) { return *gismo::gsNurbsCreator<>::BSplineSquare(r, x, y); },
      arg("r") = 1, arg("x") = 0, arg("y") = 0);

  mod.method(
      "createBSplineCube",
      [](double r = 1, double x = 0, double y = 0, double z = 0) {
        return *gismo::gsNurbsCreator<>::BSplineCube(r, x, y, z);
      },
      arg("r") = 1, arg("x") = 0, arg("y") = 0, arg("z") = 0);

  mod.method(
      "createBSplineHalfCube",
      [](double r = 1, double x = 0, double y = 0, double z = 0) {
        return *gismo::gsNurbsCreator<>::BSplineHalfCube(r, x, y, z);
      },
      arg("r") = 1, arg("x") = 0, arg("y") = 0, arg("z") = 0);

  mod.method(
      "createNurbsCube",
      [](double r = 1, double x = 0, double y = 0, double z = 0) {
        return *gismo::gsNurbsCreator<>::NurbsCube(r, x, y, z);
      },
      arg("r") = 1, arg("x") = 0, arg("y") = 0, arg("z") = 0);

  // ------------------------------------------------------------
  // Annulus and circular primitives
  mod.method(
      "createNurbsQuarterAnnulus",
      [](double r0 = 1, double r1 = 2) { return *gismo::gsNurbsCreator<>::NurbsQuarterAnnulus(r0, r1); }, arg("r0") = 1,
      arg("r1") = 2);

  mod.method(
      "createNurbsAnnulus", [](double r0 = 1, double r1 = 2) { return *gismo::gsNurbsCreator<>::NurbsAnnulus(r0, r1); },
      arg("r0") = 1, arg("r1") = 2);

  mod.method(
      "createNurbsSphere",
      [](double r = 1, double x = 0, double y = 0, double z = 0) {
        return *gismo::gsNurbsCreator<>::NurbsSphere(r, x, y, z);
      },
      arg("r") = 1, arg("x") = 0, arg("y") = 0, arg("z") = 0);

  mod.method(
      "createNurbsCircle",
      [](double r = 1, double x = 0, double y = 0) { return *gismo::gsNurbsCreator<>::NurbsCircle(r, x, y); },
      arg("r") = 1, arg("x") = 0, arg("y") = 0);

  mod.method(
      "createBSplineFatCircle",
      [](double r = 1, double x = 0, double y = 0) { return *gismo::gsNurbsCreator<>::BSplineFatCircle(r, x, y); },
      arg("r") = 1, arg("x") = 0, arg("y") = 0);

  mod.method(
      "createBSplineFatDisk",
      [](double r = 1, double x = 0, double y = 0) { return *gismo::gsNurbsCreator<>::BSplineFatDisk(r, x, y); },
      arg("r") = 1, arg("x") = 0, arg("y") = 0);

  // ------------------------------------------------------------
  // Curves and segments
  mod.method(
      "createNurbsCurve1",
      [](double r = 1, double x = 0, double y = 0) { return *gismo::gsNurbsCreator<>::NurbsCurve1(r, x, y); },
      arg("r") = 1, arg("x") = 0, arg("y") = 0);

  mod.method(
      "createNurbsCurve2",
      [](double r = 1, double x = 0, double y = 0) { return *gismo::gsNurbsCreator<>::NurbsCurve2(r, x, y); },
      arg("r") = 1, arg("x") = 0, arg("y") = 0);

  mod.method(
      "createNurbsBean",
      [](double r = 1, double x = 0, double y = 0) { return *gismo::gsNurbsCreator<>::NurbsBean(r, x, y); },
      arg("r") = 1, arg("x") = 0, arg("y") = 0);

  mod.method(
      "createBSplineE",
      [](double r = 1, double x = 0, double y = 0) { return *gismo::gsNurbsCreator<>::BSplineE(r, x, y); },
      arg("r") = 1, arg("x") = 0, arg("y") = 0);

  mod.method(
      "createNurbsAmoebaFull",
      [](double r = 1, double x = 0, double y = 0) { return *gismo::gsNurbsCreator<>::NurbsAmoebaFull(r, x, y); },
      arg("r") = 1, arg("x") = 0, arg("y") = 0);

  mod.method(
      "createBSplineSegment",
      [](double u0 = 0, double u1 = 1) { return *gismo::gsNurbsCreator<>::BSplineSegment(u0, u1); }, arg("u0") = 0,
      arg("u1") = 1);

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

// register transformation helpers (no multipatch variants)
void registerNurbsTransformFunctions(jlcxx::Module& mod) {
  using jlcxx::arg;

  // Rotation helpers
  mod.method(
      "rotate2D",
      [](const gismo::gsTensorBSpline<2, double>& geo, double turndeg = 0, double Tx = 0, double Ty = 0) {
        return *gismo::gsNurbsCreator<>::rotate2D(geo, turndeg, Tx, Ty);
      },
      arg("geo"), arg("turndeg") = 0, arg("Tx") = 0, arg("Ty") = 0);

  mod.method(
      "rotate2D!",
      [](gismo::gsGeometry<double>& geo, double turndeg = 0, double Tx = 0, double Ty = 0) {
        gismo::gsNurbsCreator<>::rotate2D(geo, turndeg, Tx, Ty);
      },
      arg("geo"), arg("turndeg") = 0, arg("Tx") = 0, arg("Ty") = 0);

  mod.method(
      "rotate3D!",
      [](gismo::gsGeometry<double>& geo, double phi_z = 0, double phi_y = 0, double phi_x = 0) {
        gismo::gsNurbsCreator<>::rotate3D(geo, phi_z, phi_y, phi_x);
      },
      arg("geo"), arg("phi_z") = 0, arg("phi_y") = 0, arg("phi_x") = 0);

  // Shifting helpers
  mod.method(
      "shift2D",
      [](const gismo::gsTensorBSpline<2, double>& geo, double dx = 0, double dy = 0, double dz = 0) {
        return *gismo::gsNurbsCreator<>::shift2D(geo, dx, dy, dz);
      },
      arg("geo"), arg("dx") = 0, arg("dy") = 0, arg("dz") = 0);

  mod.method(
      "shift2D!",
      [](gismo::gsGeometry<double>& geo, double dx = 0, double dy = 0, double dz = 0) {
        gismo::gsNurbsCreator<>::shift2D(geo, dx, dy, dz);
      },
      arg("geo"), arg("dx") = 0, arg("dy") = 0, arg("dz") = 0);

  // Mirroring helpers
  mod.method(
      "mirror2D",
      [](gismo::gsTensorBSpline<2, double>& geo, bool axis) { return *gismo::gsNurbsCreator<>::mirror2D(geo, axis); },
      arg("geo"), arg("axis"));

  mod.method(
      "mirror2D!", [](gismo::gsGeometry<double>& geo, bool axis) { gismo::gsNurbsCreator<>::mirror2D(geo, axis); },
      arg("geo"), arg("axis"));

  // Scaling helpers
  mod.method(
      "scale2D",
      [](const gismo::gsTensorBSpline<2, double>& geo, double factor = 1.0) {
        return *gismo::gsNurbsCreator<>::scale2D(geo, factor);
      },
      arg("geo"), arg("factor") = 1.0);

  mod.method(
      "scale2DVec",
      [](const gismo::gsTensorBSpline<2, double>& geo, const std::vector<double>& factors) {
        return *gismo::gsNurbsCreator<>::scale2D(geo, factors);
      },
      arg("geo"), arg("factors"));

  mod.method(
      "scale2D!",
      [](gismo::gsGeometry<double>& geo, double factor = 1.0) { gismo::gsNurbsCreator<>::scale2D(geo, factor); },
      arg("geo"), arg("factor") = 1.0);

  mod.method(
      "scale2DVec!",
      [](gismo::gsGeometry<double>& geo, const std::vector<double>& factors) {
        gismo::gsNurbsCreator<>::scale2D(geo, factors);
      },
      arg("geo"), arg("factors"));
}