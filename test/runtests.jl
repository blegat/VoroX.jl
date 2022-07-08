using Test

using StaticArrays
import CDDLib
import QHull
import VoronoiDelaunay
import MiniQhull
import HyperVoronoiDelaunay
include("../src/VoroX.jl")

function _test_grid(n, lib)
    points = [SVector{2,Float64}(x, y) for x in -n:n for y in -n:n]
    return VoroX.Foam(points, lib, HyperVoronoiDelaunay.NonPeriodic(), VoroX.Circumcenter())
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
    @test 8 <= count(iszero, foam.voronoi_edges) <= 22
    @test size(foam.active_edges) == (3, 8)
    @test 8 <= count(iszero, foam.active_edges) <= 24
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
    foam = VoroX.Foam(points, lib, HyperVoronoiDelaunay.NonPeriodic(), VoroX.Circumcenter())
    @test size(foam.voronoi_edges) == (3, 3)
end

using Random
function test_random_points(lib)
    Random.seed!(0)
    points = VoroX.random_points(10, SVector(-1.0, -1.0), SVector(1.0, 1.0), lib)
    @test points ≈ [
        [0.0,           0.0],
        [0.3800315515,  0.030874917],
        [0.2960086941, -0.358495838],
        [-0.134591813, -0.341442487],
        [0.1603910181,  0.455165450],
        [-0.529081597,  0.070297047],
        [0.2486058990, -0.917995138],
        [-0.355779492, -0.664489558],
        [-0.857012677, -0.621931050],
        [-0.731416925, -0.238293904],
    ]
end

function test_center()
    points = [
        SVector(-1,  0, -1/√2),
        SVector( 1,  0, -1/√2),
        SVector( 0, -1,  1/√2),
        SVector( 0,  1,  1/√2),
    ]
    @test VoroX.center(points, VoroX.Circumcenter()) ≈ zeros(3) atol=1e-12
    @test VoroX.center(points, VoroX.Centroid()) ≈ zeros(3) atol=1e-12
end

function hascol(A, cols)
    for col in cols
        sort!(col)
        for i in 1:size(A, 2)
            if col == sort(A[:, i])
                return true
            end
        end
    end
    return false
end

function test_periodic(algo)
    points = [
        SVector( 1/4,  1/2),
        SVector(-1/2, -1/2),
        SVector( 1/2, -1/2),
    ]
    d = HyperVoronoiDelaunay.delaunay(points, algo, HyperVoronoiDelaunay.Periodic(SVector(2.0, 2.0)))
    @test size(d) == (3, 6)
    @test hascol(d, [[13, 14, 15]])
    @test hascol(d, [[4, 14, 15], [13, 23, 24]])
    @test hascol(d, [[13, 15, 17], [10, 12, 14]])
    @test hascol(d, [[4, 15, 17], [1, 12, 14], [13, 24, 26]])
    @test hascol(d, [[10, 13, 23], [13, 16, 26], [1, 4, 14]])
    @test hascol(d, [[10, 13, 14], [13, 16, 17]])
end

LIBRARIES = VoroX.LIBRARIES

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
    test_random_points(VoroX.NN.GridTree)
end

@testset "center" begin
    test_center()
end

@testset "periodic" for lib in LIBRARIES
    test_periodic(lib)
end
