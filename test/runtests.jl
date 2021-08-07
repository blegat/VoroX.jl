using Test
using CDDLib
lib = CDDLib.Library(:float)
#lib = Polyhedra.DefaultLibrary{Float64}()
include("../src/DynamicFoam.jl")

function _test_grid(n)
    points = [SVector{2,Float64}(x, y) for x in -n:n for y in -n:n]
    return DynamicFoam.Foam(points, lib, DynamicFoam.Circumcenter())
end

function test_grid_0()
    foam = _test_grid(0)
    @test length(foam.centers) == 0
end

function test_grid_1()
    foam = _test_grid(1)
    @test size(foam.simplices) == (3, 8)
    @test length(foam.centers) == 8
    @test all(c -> all(abs.(c) .≈ 0.5), foam.centers)
    @test size(foam.voronoi_edges) == (3, 8)
    @test count(iszero, foam.voronoi_edges) == 8
    @test size(foam.active_edges) == (3, 8)
    @test count(iszero, foam.active_edges) == 8
    @test length(foam.knots) == 2
    for knot in foam.knots
        @test sort(getindex.(knot, 2)) == 1:8
    end
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

@testset "Grid" begin
    test_grid_0()
    test_grid_1()
end

@testset "random_points" begin
    test_random_points(DynamicFoam.NN.GridTree)
end
