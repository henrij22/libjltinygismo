module GismoTest
    path_to_lib = "/home/USADR/ac136110/dev/libjltinygismo/out/build/debug/lib/"

    using CxxWrap
    @wrapmodule(() -> joinpath(path_to_lib, "libjltinygismo"))

    BSplineBasis(args...) = BSplineBasis{1}(args...)
end

using .GismoTest
using BenchmarkTools
vec = [0, 0, 0, 0.5, 1, 1, 1]
kv = GismoTest.KnotVector(vec)

GismoTest.unique(kv)

basis = GismoTest.BSplineBasis(kv)
GismoTest.size(basis)

geo = GismoTest.BSpline(basis, rand(4))

geo = GismoTest.BSpline(basis, rand(4, 2))

GismoTest.boundary(geo, 1)[]

evals = GismoTest.gsMatrix{Float64}()
GismoTest.eval!(geo, [0.4], evals)
GismoTest.toMatrix(evals)

basis2 = GismoTest.TensorBSplineBasis{2}(kv, kv )
GismoTest.eval!(basis2, [0.4, 0.4], evals)
GismoTest.toMatrix(evals)

geo2 = GismoTest.TensorBSpline{2}(basis2, randn(16,3))
GismoTest.eval!(geo2, [0.4, 0.4], evals)
GismoTest.toMatrix(evals)

corners = [
    0.0 0.0 0.0
    1.0 0.0 0.0
    1.0 1.0 0.0
    0.0 1.0 0.0
    ]
geo3 = GismoTest.TensorBSpline{2}(corners, kv, kv)
GismoTest.eval!(geo3, [0.4, 0.4], evals)
GismoTest.toMatrix(evals)

vec = GismoTest.gsVector{Float64}()
GismoTest.closestPointTo(geo3, [0.5, 0.6, 0.5], vec)
GismoTest.toVector(vec)

GismoTest.toVector(GismoTest.coefAtCorner(geo3, 4))

geo4 = GismoTest.createBSplineRectangle()
GismoTest.eval!(geo4, [0.4, 0.4], evals)
GismoTest.toMatrix(evals)


boundary = GismoTest.boundary(geo4, 1)
GismoTest.eval!(boundary[], [0.4], evals)


nothing
# GismoTest.degreeIncrease(basis, 1)
# GismoTest.elevateContinuity(basis, 1)
# GismoTest.degreeIncrease(basis, 1)

# knotsNew = GismoTest.knots(basis)
# GismoTest.knotContainer(knotsNew)

# actives = GismoTest.gsMatrixInt()
# GismoTest.active!(basis, [0.4], actives)
# GismoTest.toJulia(actives)

# # GismoTest.insertKnot(basis, 0.8, 4)

# Test: Memory is shared
# evals = GismoTest.gsMatrix()
# GismoTest.eval!(basis, [0.4], evals)
# m = GismoTest.toJulia(evals)
# m[2] = 1.0
# m2 = GismoTest.toJulia(evals)
# @assert m2[2] ≈ 1.0


# GismoTest.evalSingle!(basis, 1, [0.4], evals)
# @benchmark GismoTest.evalSingle!($basis, 1, $[0.4], $evals)
# GismoTest.toJulia(evals)

# GismoTest.derivSingle!(basis, 1, [0.4], evals)
# GismoTest.toJulia(evals)

# GismoTest.evalFunc!(basis, [0.4], randn(4), evals)
# GismoTest.toJulia(evals)

# GismoTest.deriv!(basis, [0.4], randn(4), evals)
# GismoTest.toJulia(evals)

