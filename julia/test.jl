module GismoTest
    path_to_lib = "/home/USADR/ac136110/dev/libjltinygismo/out/build/debug/lib/"

    using CxxWrap
    @wrapmodule(() -> joinpath(path_to_lib, "libjltinygismo"))


    function __init__()
        return @initcxx
    end

    BSplineBasis(args...) = BSplineBasis{1}(args...)

end

using .GismoTest
using BenchmarkTools
vec = [0, 0, 0, 0.5, 1, 1, 1]
kv = GismoTest.KnotVector(vec)

GismoTest.degreeElevate!(kv)
# using BenchmarkTools
# @benchmark GismoTest.knotContainer($kv)

GismoTest.toMatrix(GismoTest.greville(kv))

GismoTest.unique(kv)

basis = GismoTest.BSplineBasis(kv)
GismoTest.size(basis)
GismoTest.elementIndex(basis, [0.4])


ks = GismoTest.knotSpans(basis)[1]

GismoTest.lowerCorner(ks)

basis = GismoTest.TensorBSplineBasis{2}(kv, kv)

boundaryIndices = GismoTest.boundary(basis, 2)
GismoTest.toMatrix(boundaryIndices)

# GismoTest.writeParaview(basis, "basis")

methods(GismoTest.BSplineBasis{Int64(2)})
# basis = GismoTest.TensorBSplineBasis{2}(kv, kv)
# GismoTest.size(basis)
# GismoTest.elementIndex(basis, [0.4, 0.8])

filename = "julia/geometry/scordelis_lo.xml"
basis = GismoTest.readFile(GismoTest.TensorNurbsBasis{2}, filename)
geo2 = GismoTest.readFile(GismoTest.TensorNurbs{2}, filename)
weights = GismoTest.weights(geo2)
GismoTest.toVector(weights)
coefs = GismoTest.coefs(geo2)
GismoTest.toMatrix(coefs)
coefs = GismoTest.coefsSize(geo2)

geo = GismoTest.BSpline(basis, rand(4))

geo = GismoTest.BSpline(basis, rand(4, 2))

GismoTest.boundary(geo, 1)[]

evals = GismoTest.gsMatrix{Float64}()
GismoTest.eval!(geo, [0.4], evals)
GismoTest.toMatrix(evals)

basis2 = GismoTest.TensorBSplineBasis{2}(kv, kv)
GismoTest.eval!(basis2, [0.4, 0.4], evals)
GismoTest.toMatrix(evals)

actives = GismoTest.gsMatrix{Int32}()
GismoTest.active!(basis2, [0.1, 0.1], actives)
GismoTest.toMatrix(actives)

geo2 = GismoTest.TensorBSpline{2}(basis2, randn(16, 3))
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

GismoTest._value(vec, 2)

GismoTest.toVector(GismoTest.coefAtCorner(geo3, 4))

geo4 = GismoTest.createBSplineRectangle()
GismoTest.eval!(geo4, [0.4, 0.4], evals)
GismoTest.toMatrix(evals)

GismoTest.toMatrix(evals)
GismoTest.value(evals, 2, 1)

boundary = GismoTest.boundary(geo4, 1)
GismoTest.eval!(boundary[], [0.4], evals)


basis = GismoTest.BSplineBasis(kv)

result = GismoTest.gsMatrix{Float64}()
GismoTest.eval!(basis, [0.4], result)
GismoTest.toMatrix(result)

coefs = randn(4, 1)
spline = GismoTest.BSpline(basis, coefs)

result = GismoTest._eval(basis, [0.4])

result = GismoTest.gsMatrix{Float64}()
using BenchmarkTools
@benchmark GismoTest.evalFunc!($basis, $[0.4], $coefs, $result)

nothing
GismoTest.degreeIncrease!(basis, 1)
GismoTest.elevateContinuity(basis, 1)
GismoTest.degreeIncrease(basis, 1)

knotsNew = GismoTest.knots(basis)
GismoTest.knotContainer(knotsNew)

actives = GismoTest.gsMatrix{Int32}()
GismoTest.active!(basis, [0.4], actives)
GismoTest.toVector(actives)

# GismoTest.insertKnot(basis, 0.8, 4)
# 
# Test: Memory is shared
evals = GismoTest.gsMatrix{Float64}()
GismoTest.eval!(basis, [0.4], evals)
m = GismoTest.toMatrix(evals)
m[2] = 1.0
m2 = GismoTest.toJulia(evals)
@assert m2[2] ≈ 1.0


GismoTest.evalSingle!(basis, 1, [0.4], evals)
@benchmark GismoTest.evalSingle!($basis, 1, $[0.4], $evals)
GismoTest.toJulia(evals)

GismoTest.derivSingle!(basis, 1, [0.4], evals)
GismoTest.toJulia(evals)

GismoTest.evalFunc!(basis, [0.4], randn(4), evals)
GismoTest.toJulia(evals)

GismoTest.deriv!(basis, [0.4], randn(4), evals)
GismoTest.toJulia(evals)
