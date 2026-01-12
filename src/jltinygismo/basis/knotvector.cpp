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

  kv.method("size", &KnotVector::size);
  kv.method("uSize", &KnotVector::uSize);
  kv.method("numElements", &KnotVector::numElements);
  kv.method("unique", &KnotVector::unique);
  kv.method("multiplicities", &KnotVector::multiplicities);
  kv.method("knotContainer",
            [](const KnotVector& knotVector) { return jlcxx::make_const_array(knotVector.data(), knotVector.size()); });
}