# Integration test suite for the libjltinygismo bindings.
#
# Usage:  julia --project=julia julia/runtests.jl <path_to_lib_dir>
#
# These tests load the freshly built shared library directly, so they exercise the actual
# registration code rather than a released TinyGismo_jll. They are the real coverage for the
# bindings -- the C++ unit tests in test/ can only reach helpers that do not touch the Julia
# runtime.

using Test

const LIBDIR = if !isempty(ARGS)
    ARGS[1]
elseif haskey(ENV, "JLTINYGISMO_LIB")
    ENV["JLTINYGISMO_LIB"]
else
    error("Usage: julia --project=julia julia/runtests.jl <path_to_lib_dir>")
end

isdir(LIBDIR) || error("Library directory not found: $LIBDIR")

module G
using CxxWrap
@wrapmodule(() -> joinpath(Main.LIBDIR, "libjltinygismo"))
__init__() = @initcxx

# Mirrors the convenience constructors TinyGismo.jl defines on top of the raw bindings.
BSplineBasis(args...) = BSplineBasis{1}(args...)
NurbsBasis(args...) = NurbsBasis{1}(args...)
gsMatrix(args...) = gsMatrix{Float64}(args...)
gsVector(args...) = gsVector{Float64}(args...)
end

using .G

# Shorthands used throughout
mat(x) = G.toMatrix(x)
vec_(x) = G.toVector(x)
evalat(obj, u) = (out = G.gsMatrix(); G.eval!(obj, u, out); mat(out))

const QUADRATIC = [0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0]

@testset "libjltinygismo" begin

    @testset "KnotVector" begin
        kv = G.KnotVector(QUADRATIC)
        @test G.size(kv) == 7
        @test G.uSize(kv) == 3
        @test G.numElements(kv) == 2
        @test G.knotContainer(kv) ≈ QUADRATIC
        @test G.unique(kv) ≈ [0.0, 0.5, 1.0]
        @test G.multiplicities(kv) == [3, 1, 3]

        @testset "degree is a getter" begin
            # Used to be bound as `degree!`, which called the const getter and threw the
            # result away, so it always returned `nothing`.
            @test G.degree(kv) == 2
            @test G.degree(G.KnotVector([0.0, 0.0, 1.0, 1.0])) == 1
        end

        @testset "multiplicity of a single knot" begin
            kv2 = G.KnotVector([0.0, 0.0, 0.0, 0.25, 0.5, 0.5, 1.0, 1.0, 1.0])
            @test G.multiplicity(kv2, 0.0) == 3
            @test G.multiplicity(kv2, 0.5) == 2
            @test G.multiplicity(kv2, 0.25) == 1
            @test G.multiplicity(kv2, 0.3) == 0   # not a knot
        end

        @testset "greville" begin
            kv2 = G.KnotVector([0.0, 0.0, 0.0, 0.25, 0.5, 0.5, 1.0, 1.0, 1.0])
            expected = [0.0, 0.125, 0.375, 0.5, 0.75, 1.0]

            # The no-argument overload returns every abscissa as a 1 x n row.
            @test vec_(G.greville(kv2)) ≈ expected

            out = G.gsMatrix()
            G.greville!(kv2, out)
            @test vec_(out) ≈ expected

            # The index used to be bound as a Bool, collapsing every index onto the first
            # two abscissae.
            @test G.greville(kv2, 1) ≈ expected[1]
            @test G.greville(kv2, 3) ≈ expected[3]
            @test_throws Exception G.greville(kv2, 0)
            @test_throws Exception G.greville(kv2, 99)
        end

        @testset "refinement and degree operations" begin
            kv2 = G.KnotVector(QUADRATIC)
            G.uniformRefine!(kv2)
            @test G.numElements(kv2) == 4

            kv3 = G.KnotVector(QUADRATIC)
            G.degreeElevate!(kv3)
            @test G.size(kv3) == 10   # every unique knot gains one multiplicity

            kv4 = G.KnotVector(QUADRATIC)
            G.degreeIncrease!(kv4)
            @test G.size(kv4) == 9    # only the end knots repeat
        end
    end

    @testset "BSplineBasis" begin
        kv = G.KnotVector(QUADRATIC)
        basis = G.BSplineBasis(kv)

        @test G.size(basis) == 4
        @test G.degree(basis) == 2
        @test G.order(basis) == 3
        @test G.numElements(basis) == 2
        @test G.numActive(basis) == 3
        @test G.knotContainer(G.knots(basis)) ≈ QUADRATIC

        @testset "evaluation" begin
            N = evalat(basis, [0.4])
            @test size(N) == (3, 1)
            @test sum(N) ≈ 1.0                       # partition of unity

            actives = G.gsMatrix{Int32}()
            G.active!(basis, [0.4], actives)
            @test vec_(actives) == Int32[1, 2, 3]    # 1-based

            @test G.isActive(basis, 1, [0.4])
            @test !G.isActive(basis, 4, [0.4])
            @test G.elementIndex(basis, [0.4]) == 1
            @test G.elementIndex(basis, [0.6]) == 2

            # evalSingle agrees with the corresponding entry of eval!
            for i in 1:3
                @test mat(G.evalSingle(basis, i, [0.4]))[1] ≈ N[i]
            end
        end

        @testset "derivatives" begin
            dN = G.gsMatrix()
            G.deriv!(basis, [0.4], dN)
            @test sum(mat(dN)) ≈ 0.0 atol = 1e-12    # differentiated partition of unity

            for i in 1:3
                @test mat(G.derivSingle(basis, i, [0.4]))[1] ≈ mat(dN)[i]
            end
            @test G.deriv2Single(basis, 2, [0.4]) isa Float64
        end

        @testset "evalFunc matches manual assembly" begin
            coefs = [1.0, 2.0, 0.0, -1.0]
            N = vec_(G._eval(basis, [0.4]))
            actives = G.gsMatrix{Int32}()
            G.active!(basis, [0.4], actives)
            idx = vec_(actives)

            manual = sum(N[k] * coefs[idx[k]] for k in eachindex(idx))
            @test mat(G.evalFunc(basis, [0.4], coefs))[1] ≈ manual
        end

        @testset "refinement" begin
            b = G.BSplineBasis(G.KnotVector(QUADRATIC))
            G.uniformRefine!(b)
            @test G.numElements(b) == 4
            @test G.size(b) == 6
            G.uniformCoarsen!(b)
            @test G.numElements(b) == 2

            b2 = G.BSplineBasis(G.KnotVector(QUADRATIC))
            G.insertKnot!(b2, 0.25)
            @test G.knotContainer(G.knots(b2)) ≈ [0, 0, 0, 0.25, 0.5, 1, 1, 1]
            G.insertKnot!(b2, 0.75, 2)
            @test G.multiplicity(G.knots(b2), 0.75) == 2
            G.removeKnot!(b2, 0.25)
            @test G.multiplicity(G.knots(b2), 0.25) == 0
        end

        @testset "uniformRefine_withCoefs! returns the refined coefficients" begin
            # Used to write back in place, which cannot work: refinement adds control
            # points. It also built its Eigen map from the wrong dimensions.
            kv2 = G.KnotVector(QUADRATIC)
            cp = [0.0 0.0; 1.0 1.5; 2.0 -0.5; 3.0 1.0]

            before = evalat(G.BSpline(G.BSplineBasis(kv2), cp), [0.3])

            b = G.BSplineBasis(kv2)
            newcoefs = mat(G.uniformRefine_withCoefs!(b, copy(cp)))

            @test size(newcoefs) == (6, 2)
            @test G.size(b) == 6
            # The refined representation is the same curve.
            @test evalat(G.BSpline(b, newcoefs), [0.3]) ≈ before

            @test_throws Exception G.uniformRefine_withCoefs!(G.BSplineBasis(kv2), zeros(3, 2))
        end

        @testset "continuity operations take no direction" begin
            b = G.BSplineBasis(G.KnotVector(QUADRATIC))
            G.elevateContinuity!(b, 1)
            @test G.numElements(b) == 1      # the single interior knot is removed
            # The bogus three-argument method is gone.
            @test_throws MethodError G.elevateContinuity!(b, 1, -1)
        end
    end

    @testset "TensorBSplineBasis" begin
        kvu = G.KnotVector([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])   # quadratic, 3 functions
        kvv = G.KnotVector([0.0, 0.0, 1.0, 1.0])             # linear, 2 functions
        tb = G.TensorBSplineBasis{2}(kvu, kvv)

        @test G.size(tb) == 6
        @test G.degree(tb, 1) == 2
        @test G.degree(tb, 2) == 1
        @test G.numActive(tb) == 6                            # (2+1) * (1+1)
        @test G.size(G.component(tb, 1)) == 3
        @test G.size(G.component(tb, 2)) == 2
        @test sum(evalat(tb, [0.5, 0.5])) ≈ 1.0

        @testset "trivariate constructor" begin
            # gsTensorBSplineBasis{3} used to expose only the default constructor.
            tb3 = G.TensorBSplineBasis{3}(kvu, kvv, kvv)
            @test G.size(tb3) == 12
            @test G.numActive(tb3) == 12
            @test (G.degree(tb3, 1), G.degree(tb3, 2), G.degree(tb3, 3)) == (2, 1, 1)
            @test sum(evalat(tb3, [0.5, 0.5, 0.5])) ≈ 1.0
        end

        @testset "numElements second argument is a box side" begin
            b = G.basis(G.createBSplineRectangle())
            G.uniformRefine!(b, 2)          # 3 x 3 elements
            @test G.numElements(b) == 9
            @test G.numElements(G.component(b, 1)) == 3
            @test G.numElements(G.component(b, 2)) == 3
            @test all(G.numElements(b, s) == 3 for s in 1:4)
        end
    end

    @testset "TensorNurbsBasis" begin
        kvu = G.KnotVector([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
        kvv = G.KnotVector([0.0, 0.0, 1.0, 1.0])
        w = ones(6, 1)
        tnb = G.TensorNurbsBasis{2}(kvu, kvv, w)
        @test G.size(tnb) == 6
        @test G.numActive(tnb) == 6

        tnb3 = G.TensorNurbsBasis{3}(kvu, kvv, kvv, ones(12, 1))
        @test G.size(tnb3) == 12
        @test G.numActive(tnb3) == 12
    end

    @testset "BSpline geometry" begin
        kv = G.KnotVector(QUADRATIC)
        cp = [0.0 0.0; 1.0 1.5; 2.0 -0.5; 3.0 1.0]
        curve = G.BSpline(G.BSplineBasis(kv), cp)

        @test G.parDim(curve) == 1
        @test G.targetDim(curve) == 2
        @test G.geoDim(curve) == 2
        @test G.coefDim(curve) == 2
        @test G.coefsSize(curve) == 4
        @test G.numCoefs(curve) == 4
        @test mat(G.coefs(curve)) ≈ cp
        @test vec_(G.coefAtCorner(curve, 1)) ≈ cp[1, :]

        @test evalat(curve, [0.0]) ≈ reshape(cp[1, :], 2, 1)
        @test evalat(curve, [1.0]) ≈ reshape(cp[end, :], 2, 1)

        @testset "derivative agrees with finite differences" begin
            u, h = 0.3, 1e-6
            fd = (vec(evalat(curve, [u + h])) .- vec(evalat(curve, [u - h]))) ./ (2h)
            @test vec(mat(G.deriv(curve, [u]))) ≈ fd atol = 1e-6
            @test vec(mat(G.jacobian(curve, [u]))) ≈ fd atol = 1e-6
        end

        @testset "closestPointTo" begin
            par = G.gsVector()
            onCurve = vec(evalat(curve, [0.42]))
            G.closestPointTo(curve, onCurve, par)
            @test vec(evalat(curve, vec_(par))) ≈ onCurve atol = 1e-6
        end

        @testset "refinement preserves the curve" begin
            c = G.BSpline(G.BSplineBasis(kv), cp)
            before = evalat(c, [0.3])
            G.uniformRefine!(c)
            @test G.coefsSize(c) == 6
            @test evalat(c, [0.3]) ≈ before

            G.degreeElevate!(c)
            @test G.degree(c) == 3
            @test evalat(c, [0.3]) ≈ before
        end

        @testset "scalar-valued spline" begin
            scalar = G.BSpline(G.BSplineBasis(kv), [0.0, 1.0, -1.0, 0.5])
            @test G.parDim(scalar) == 1
            @test G.targetDim(scalar) == 1
        end
    end

    @testset "Nurbs geometry" begin
        circle = G.createNurbsCircle(1.0)

        @testset "gsGeometry methods are reachable" begin
            # gsNurbs used to be missing its SuperType specialization, so every inherited
            # gsGeometry method failed with an upcast error.
            @test G.parDim(circle) == 1
            @test G.targetDim(circle) == 2
            @test G.coefsSize(circle) == 9
            @test size(mat(G.coefs(circle))) == (9, 2)
        end

        @testset "it really is a circle" begin
            for u in 0.0:0.1:1.0
                @test sqrt(sum(abs2, evalat(circle, [u]))) ≈ 1.0
            end
        end

        weights = vec_(G.weights(circle))
        @test length(weights) == 9
        @test G.weight(circle, 2) ≈ weights[2]
        @test weights[2] ≈ sqrt(2) / 2
    end

    @testset "Tensor geometries" begin
        rect = G.createBSplineRectangle(0.0, 0.0, 2.0, 1.0)
        @test G.parDim(rect) == 2
        @test G.targetDim(rect) == 2
        @test evalat(rect, [0.5, 0.5]) ≈ [1.0; 0.5;;]

        @testset "affine map has a constant Jacobian and vanishing Hessian" begin
            @test mat(G.jacobian(rect, [0.5, 0.5])) ≈ [2.0 0.0; 0.0 1.0]
            @test mat(G.hessian(rect, [0.5, 0.5], 1)) ≈ zeros(2, 2) atol = 1e-12
        end

        @testset "boundary" begin
            west = G.boundary(rect, 1)[]
            @test G.parDim(west) == 1
            @test mat(G.coefs(west)) ≈ [0.0 0.0; 0.0 0.5; 0.0 1.0]
            @test vec_(G.boundary(G.basis(rect), 1)) == Int32[1, 4, 7]
        end

        @testset "corner constructor" begin
            kvlin = G.KnotVector([0.0, 0.0, 1.0, 1.0])
            corners = [0.0 0.0 0.0; 4.0 0.0 0.0; 3.0 2.0 0.0; 1.0 2.0 0.0]
            patch = G.TensorBSpline{2}(corners, kvlin, kvlin)
            @test G.targetDim(patch) == 3
            @test vec(evalat(patch, [0.0, 0.0])) ≈ corners[1, :]
            @test vec(evalat(patch, [1.0, 0.0])) ≈ corners[2, :]
            @test vec(evalat(patch, [1.0, 1.0])) ≈ corners[3, :]
            @test vec(evalat(patch, [0.0, 1.0])) ≈ corners[4, :]
        end

        @testset "quarter annulus is exactly circular" begin
            ann = G.createNurbsQuarterAnnulus(1.0, 2.0)
            for v in 0.0:0.125:1.0
                @test sqrt(sum(abs2, evalat(ann, [0.0, v]))) ≈ 1.0   # inner edge
                @test sqrt(sum(abs2, evalat(ann, [1.0, v]))) ≈ 2.0   # outer edge
            end
        end

        @testset "trivariate" begin
            cube = G.createBSplineCube(1.0)
            @test G.parDim(cube) == 3
            @test G.targetDim(cube) == 3
            @test G.coefsSize(cube) == 8
            @test vec(evalat(cube, [0.5, 0.5, 0.5])) ≈ [0.5, 0.5, 0.5]
        end
    end

    @testset "Direction arguments" begin
        # `dir` is 1-based with 0 meaning "all". Out-of-range values used to be forwarded
        # to G+Smo, which indexes without bounds checks and segfaults the process.
        @testset "valid directions" begin
            r = G.createBSplineRectangle()
            G.degreeElevate!(r, 1, 1)
            @test (G.degree(r, 1), G.degree(r, 2)) == (3, 2)

            r2 = G.createBSplineRectangle()
            G.degreeElevate!(r2, 1, 2)
            @test (G.degree(r2, 1), G.degree(r2, 2)) == (2, 3)

            r3 = G.createBSplineRectangle()
            G.degreeElevate!(r3, 1, 0)
            @test (G.degree(r3, 1), G.degree(r3, 2)) == (3, 3)

            r4 = G.createBSplineRectangle()
            G.uniformRefine!(r4, 1, 1, 1)
            b = G.basis(r4)
            @test (G.numElements(G.component(b, 1)), G.numElements(G.component(b, 2))) == (2, 1)
        end

        @testset "invalid directions raise instead of crashing" begin
            for f in (G.degreeElevate!, G.degreeIncrease!, G.degreeDecrease!, G.degreeReduce!)
                @test_throws Exception f(G.createBSplineRectangle(), 1, -1)
                @test_throws Exception f(G.createBSplineRectangle(), 1, 3)
            end
            @test_throws Exception G.uniformRefine!(G.createBSplineRectangle(), 1, 1, -1)
            @test_throws Exception G.insertKnot!(G.createBSplineRectangle(), 0.5, 3)
            # A univariate geometry only has direction 1.
            @test_throws Exception G.uniformRefine!(G.createBSplineUnitInterval(2), 1, 1, 2)
        end
    end

    @testset "gsMatrix / gsVector" begin
        m = G.gsMatrix(3, 2)
        @test G.rows(m) == 3
        @test G.cols(m) == 2
        @test G.size(m) == 6

        for j in 1:2, i in 1:3
            G.setValue!(m, i, j, 10i + j)
        end
        @test mat(m) == [11.0 12.0; 21.0 22.0; 31.0 32.0]
        @test G.value(m, 2, 1) == 21.0
        @test_throws Exception G.value(m, 9, 1)
        @test_throws Exception G.setValue!(m, 0, 1, 1.0)
        @test_throws Exception G.setValue!(m, 9, 1, 1.0)

        v = G.gsVector(3)
        for i in 1:3
            G.setValue!(v, i, float(i^2))
        end
        @test vec_(v) == [1.0, 4.0, 9.0]
        @test G.value(v, 3) == 9.0

        @testset "toVector accepts rows and columns" begin
            col = G.gsMatrix(3, 1)
            @test length(vec_(col)) == 3
            kv = G.KnotVector([0.0, 0.0, 1.0, 1.0])
            @test length(vec_(G.greville(kv))) == 2   # a 1 x n row
            @test_throws Exception vec_(G.gsMatrix(2, 2))
        end

        @testset "a gsMatrix can be handed to a constructor" begin
            # Only possible now that setValue! exists.
            kvu = G.KnotVector([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
            kvv = G.KnotVector([0.0, 0.0, 1.0, 1.0])
            cp = [0.0 0.0; 0.5 0.4; 1.0 0.0; 0.0 1.0; 0.5 1.4; 1.0 1.0]
            cm = G.gsMatrix(6, 2)
            for j in 1:2, i in 1:6
                G.setValue!(cm, i, j, cp[i, j])
            end
            fromKnots = G.TensorBSpline{2}(kvu, kvv, cm)
            fromBasis = G.TensorBSpline{2}(G.TensorBSplineBasis{2}(kvu, kvv), cp)
            @test evalat(fromKnots, [0.5, 0.5]) ≈ evalat(fromBasis, [0.5, 0.5])
        end
    end

    @testset "Returned arrays own their memory" begin
        # toMatrix/toVector/knotContainer used to hand back a view into a C++ buffer. When
        # the source was a temporary -- which it is in every one of these expressions -- the
        # collector could free it and the array silently read freed memory.
        kv = G.KnotVector(QUADRATIC)
        basis = G.BSplineBasis(kv)
        rect = G.createBSplineRectangle(0.0, 0.0, 2.0, 1.0)
        ann = G.createNurbsQuarterAnnulus(1.0, 2.0)

        results = [
            ("knotContainer(knots(basis))", () -> G.knotContainer(G.knots(basis)), QUADRATIC),
            ("toMatrix(coefs(rect))", () -> mat(G.coefs(rect))[1:3, :], [0.0 0.0; 1.0 0.0; 2.0 0.0]),
            ("toVector(weights(ann))", () -> vec_(G.weights(ann)), [1, 1, √2 / 2, √2 / 2, 1, 1]),
            ("toMatrix(_eval(basis,u))", () -> mat(G._eval(basis, [0.4])), [0.04, 0.64, 0.32]),
            ("toMatrix(jacobian(rect,u))", () -> mat(G.jacobian(rect, [0.5, 0.5])), [2.0 0.0; 0.0 1.0]),
            ("toVector(greville(kv))", () -> vec_(G.greville(kv)), [0.0, 0.25, 0.75, 1.0]),
        ]

        for (name, f, expected) in results
            @testset "$name" begin
                got = f()
                GC.gc()
                GC.gc()
                @test vec(got) ≈ vec(expected)
            end
        end
    end

    @testset "File I/O" begin
        xml = joinpath(@__DIR__, "geometry", "scordelis_lo.xml")
        if isfile(xml)
            geo = G.readFile(G.TensorNurbs{2}, xml)
            @test G.parDim(geo) == 2
            @test G.targetDim(geo) == 3
            basis = G.readFile(G.TensorNurbsBasis{2}, xml)
            @test G.size(basis) > 0
        else
            @info "skipping file I/O tests, $xml not found"
        end

        @test_throws Exception G.readFile(G.TensorNurbs{2}, "definitely_not_a_file.xml")
    end

    @testset "Geometry factories" begin
        # parDim / targetDim of each factory, which the TinyGismo.jl docstrings quote.
        expectations = [
            (G.createBSplineUnitInterval(2), 1, 1),
            (G.createBSplineRectangle(), 2, 2),
            (G.createBSplineSquare(), 2, 2),
            (G.createBSplineTriangle(), 2, 2),
            (G.createBSplineSegment(), 1, 1),
            (G.createBSplineCube(), 3, 3),
            (G.createBSplineHalfCube(), 3, 3),
            (G.createNurbsCube(), 3, 3),
            (G.createNurbsSphere(), 2, 3),     # a surface, despite the name
            (G.createBSplineFatCircle(), 1, 2),  # a curve
            (G.createBSplineFatDisk(), 2, 2),  # a surface
            (G.createNurbsQuarterAnnulus(), 2, 2),
            (G.createNurbsAnnulus(), 2, 2),
            (G.createNurbsCircle(), 1, 2),
        ]
        for (geo, parDim, targetDim) in expectations
            @test G.parDim(geo) == parDim
            @test G.targetDim(geo) == targetDim
            @test G.coefsSize(geo) > 0
        end
    end

    @testset "Knot spans" begin
        basis = G.BSplineBasis(G.KnotVector(QUADRATIC))
        spans = G.knotSpans(basis)
        @test length(spans) == 2
        @test vec_(G.lowerCorner(spans[1])) ≈ [0.0]
        @test vec_(G.upperCorner(spans[1])) ≈ [0.5]
        @test vec_(G.centerPoint(spans[1])) ≈ [0.25]

        @testset "tensor spans tile the domain" begin
            tb = G.basis(G.createBSplineRectangle())
            G.uniformRefine!(tb)
            spans2 = G.knotSpans(tb)
            @test length(spans2) == 4
            area = sum(spans2) do span
                lo, up = vec_(G.lowerCorner(span)), vec_(G.upperCorner(span))
                prod(up .- lo)
            end
            @test area ≈ 1.0    # the spans partition the unit parameter square
        end

        @testset "elementInSupportOf returns an element inside the support" begin
            b = G.BSplineBasis(G.KnotVector([0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0]))
            for i in 1:G.size(b)
                span = mat(G.elementInSupportOf(b, i))
                @test length(span) == 2
                @test span[1] < span[2]
                @test G.isActive(b, i, [(span[1] + span[2]) / 2])
            end
        end
    end
end
