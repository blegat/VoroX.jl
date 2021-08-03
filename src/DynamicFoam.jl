module DynamicFoam

using LinearAlgebra
using Polyhedra, GLMakie, StaticArrays

include("center.jl")
include("nearest_neighbors.jl")
include("foam.jl")

end
