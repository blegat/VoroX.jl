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

@testset "Grid" begin
    test_grid_0()
end
