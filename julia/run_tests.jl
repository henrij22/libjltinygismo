# Test runner for jlTinyGismo
# Usage: julia run_tests.jl <path_to_lib>

if length(ARGS) < 1
    error("Usage: julia run_tests.jl <path_to_lib>")
end

lib_path = ARGS[1]

module GismoTest
    using CxxWrap
    @wrapmodule(() -> joinpath(ARGS[1], "libjltinygismo"))

    function __init__()
        return @initcxx
    end

    BSplineBasis(args...) = BSplineBasis{1}(args...)
end

using .GismoTest

# Run basic tests
try
    # Test 1: KnotVector
    vec = [0, 0, 0, 0.5, 1, 1, 1]
    kv = GismoTest.KnotVector(vec)
    unique_kv = GismoTest.unique(kv)
    println("✓ Test 1: KnotVector creation passed")
    
    # Test 2: BSplineBasis
    basis = GismoTest.BSplineBasis(kv)
    size_basis = GismoTest.size(basis)
    println("✓ Test 2: BSplineBasis creation passed")
    
    # Test 3: TensorBSplineBasis
    basis2d = GismoTest.TensorBSplineBasis{2}(kv, kv)
    println("✓ Test 3: TensorBSplineBasis creation passed")
    
    # Test 4: File reading (if file exists)
    filename = "julia/geometry/scordelis_lo.xml"
    if isfile(filename)
        basis_nurbs = GismoTest.readFile(GismoTest.TensorNurbsBasis{2}, filename)
        geo_nurbs = GismoTest.readFile(GismoTest.TensorNurbs{2}, filename)
        println("✓ Test 4: File reading passed")
    else
        println("⊘ Test 4: Skipped (file not found)")
    end
    
    # Test 5: BSpline geometry
    geo = GismoTest.BSpline(basis, rand(4))
    evals = GismoTest.gsMatrix{Float64}()
    GismoTest.eval!(geo, [0.4], evals)
    println("✓ Test 5: BSpline evaluation passed")
    
    # Test 6: TensorBSpline
    corners = [
        0.0 0.0 0.0
        1.0 0.0 0.0
        1.0 1.0 0.0
        0.0 1.0 0.0
    ]
    geo3 = GismoTest.TensorBSpline{2}(corners, kv, kv)
    GismoTest.eval!(geo3, [0.4, 0.4], evals)
    println("✓ Test 6: TensorBSpline evaluation passed")
    
    # Test 7: Rectangle creation
    geo4 = GismoTest.createBSplineRectangle()
    GismoTest.eval!(geo4, [0.4, 0.4], evals)
    println("✓ Test 7: BSpline rectangle creation passed")
    
    println("\n✓✓✓ All tests passed! ✓✓✓")
catch e
    println("✗ Test failed with error:")
    println(e)
    rethrow()
end
