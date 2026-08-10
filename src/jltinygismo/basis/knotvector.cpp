#include <jltinygismo/helper.hh>

#include <jlcxx/const_array.hpp>
#include <jlcxx/array.hpp>
#include <jlcxx/jlcxx.hpp>

#include <gsNurbs/gsKnotVector.h>

void registerKnotVector(jlcxx::Module& mod) {
  using jlcxx::arg;
  using jlcxx::julia_base_type;
  using JuliaVector = jlcxx::ArrayRef<double, 1>;

  using KnotVector = gismo::gsKnotVector<>;

  auto kv = mod.add_type<KnotVector>("KnotVector");

  kv.constructor([](JuliaVector knots) {
    auto knotVec = std::vector<double>{knots.begin(), knots.end()};
    return new KnotVector{knotVec};
  });

  kv.method("size", [](const KnotVector& kv) { return kv.size(); }, arg("kv"));
  kv.method("uSize", [](const KnotVector& kv) { return kv.uSize(); }, arg("kv"));
  kv.method("numElements", [](const KnotVector& kv) { return kv.numElements(); }, arg("kv"));
  kv.method("unique", [](const KnotVector& kv) { return kv.unique(); }, arg("kv"));
  kv.method("multiplicity", [](const KnotVector& kv, double u) { return kv.multiplicity(u); }, arg("kv"), arg("u"));
  kv.method("multiplicities", [](const KnotVector& kv) { return kv.multiplicities(); }, arg("kv"));
  kv.method(
      "knotContainer",
      [](const KnotVector& knotVector) { return copyToJuliaVector(knotVector.data(), knotVector.size()); },
      arg("kv"));

  // gsKnotVector::degree is a const getter, not a mutator -- the old "degree!" binding
  // called it and discarded the result, so it always returned nothing.
  kv.method("degree", [](const KnotVector& kv) { return kv.degree(); }, arg("kv"));
  kv.method("degreeIncrease!", [](KnotVector& kv, int i = 1) { kv.degreeIncrease(i); }, arg("kv"), arg("i") = 1);
  kv.method(
      "degreeDecrease!",
      [](KnotVector& kv, int i = 1, bool updateInterior = false) { kv.degreeDecrease(i, updateInterior); }, arg("kv"),
      arg("i") = 1, arg("updateInterior") = false);

  kv.method("degreeElevate!", [](KnotVector& kv, int i = 1) { kv.degreeElevate(i); }, arg("kv"), arg("i") = 1);

  kv.method(
      "uniformRefine!", [](KnotVector& kv, int numKnots = 1, int mult = 1) { kv.uniformRefine(numKnots, mult); },
      arg("kv"), arg("numKnots") = 1, arg("mult") = 1);

  // The index used to be bound as `bool`, which silently collapsed every index onto the
  // first two abscissae. Bind the no-argument overload (all abscissae) and a properly
  // 1-based indexed one.
  kv.method("greville", [](const KnotVector& kv) { return kv.greville(); }, arg("kv"));
  kv.method(
      "greville",
      [](const KnotVector& kv, int i) {
        if (i < 1 || i > kv.size() - kv.degree() - 1)
          throw std::runtime_error("greville: index " + std::to_string(i) + " out of range 1:" +
                                   std::to_string(kv.size() - kv.degree() - 1));
        return kv.greville(i - 1);
      },
      arg("kv"), arg("i"));
  kv.method(
      "greville!", [](const KnotVector& kv, gismo::gsMatrix<>& out) { kv.greville_into(out); }, arg("kv"), arg("out"));
}