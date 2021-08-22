using Test

using StaticArrays
import CDDLib
import QHull
import VoronoiDelaunay
include("../src/DynamicFoam.jl")

function _test_grid(n, lib)
    points = [SVector{2,Float64}(x, y) for x in -n:n for y in -n:n]
    return DynamicFoam.Foam(points, lib, DynamicFoam.Circumcenter())
end

function test_grid_0(lib)
    foam = _test_grid(0, lib)
    @test length(foam.centers) == 0
end

function test_grid_1(lib)
    foam = _test_grid(1, lib)
    @test size(foam.simplices) == (3, 8)
    @test length(foam.centers) == 8
    @test all(c -> all(abs.(c) .≈ 0.5), foam.centers)
    @test size(foam.voronoi_edges) == (3, 8)
    @test 8 <= count(iszero, foam.voronoi_edges) <= 18
    @test size(foam.active_edges) == (3, 8)
    @test 8 <= count(iszero, foam.active_edges) <= 20
    @test 0 <= length(foam.knots) <= 2
    for knot in foam.knots
        @test sort(getindex.(knot, 2)) == 1:8
    end
    if length(foam.knots) == 2
        for i in 1:8
            @test sort(foam.facet_knot[:,i]) == 0:2
            for j in 1:3
                if iszero(foam.facet_knot[j,i])
                    @test foam.knot_dist[j,i] == 1
                else
                    @test foam.knot_dist[j,i] == 0
                end
            end
        end
    end
    @test length(foam.num_catched) == length(foam.knots)
end

# Test failing with VoronoiDelaunay
# See https://github.com/JuliaGeometry/VoronoiDelaunay.jl/issues/55
function test_issue_55(lib)
    points = [
        SVector(0.0, 0.0),
        SVector(-0.8, -0.4),
        SVector( 0.8, -0.2),
        SVector( 0.0,  0.6),
    ]
    foam = DynamicFoam.Foam(points, lib, DynamicFoam.Circumcenter())
    @test size(foam.voronoi_edges) == (3, 3)
end

using Random
function test_random_points(algo)
    Random.seed!(0)
    points = DynamicFoam.random_points(10, SVector(-1.0, -1.0), SVector(1.0, 1.0), algo)
    @test points ≈ [
        [ 0.0,          0.0],
        [ 0.670087346,  0.153600819],
        [ 0.666462902, -0.503069459],
        [ 0.977731941,  0.407310927],
        [ 0.344192411,  0.514585898],
        [ 0.235859624, -0.292309593],
        [-0.359376319, -0.169391956],
        [ 0.164568279, -0.664826577],
        [-0.359350030, -0.928460908],
        [-0.240282715,  0.438854948],
    ]
end

function test_center()
    points = [
        SVector(-1,  0, -1/√2),
        SVector( 1,  0, -1/√2),
        SVector( 0, -1,  1/√2),
        SVector( 0,  1,  1/√2),
    ]
    @test DynamicFoam.center(points, DynamicFoam.Circumcenter()) ≈ zeros(3) atol=1e-12
    @test DynamicFoam.center(points, DynamicFoam.Centroid()) ≈ zeros(3) atol=1e-12
end

LIBRARIES = [
    CDDLib.Library(:float),
    VoronoiDelaunay.DelaunayTessellation2D,
    QHull.Library(),
    MiniQhull.delaunay,
]

@testset "Test issue 55 $lib" for lib in LIBRARIES
    if lib != VoronoiDelaunay.DelaunayTessellation2D
        test_issue_55(lib)
    end
end

@testset "Grid $lib" for lib in LIBRARIES
    @testset "0" begin
        if lib != VoronoiDelaunay.DelaunayTessellation2D && !isa(lib, QHull.Library) && lib != MiniQhull.delaunay
            test_grid_0(lib)
        end
    end
    @testset "1" begin
        test_grid_1(lib)
    end
end

@testset "random_points" begin
    test_random_points(DynamicFoam.NN.GridTree)
end

@testset "center" begin
    test_center()
end
