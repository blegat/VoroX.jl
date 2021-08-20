module DynamicFoam

using LinearAlgebra
using GLMakie, StaticArrays
import ProgressMeter
import ColorSchemes

include("center.jl")
include("nearest_neighbors.jl")
include("rand.jl")
include("foam.jl")
include("motion.jl")
include("plot.jl")
#include("plotting.jl")


#function main(K, a, b, delaunay_algo, centering=Barycenter(), scatterfun=scatter!, args...)
#    points = random_points(K, a, b, args...)
#    foam = Foam(points, delaunay_algo, centering)
#    return viz(foam, K, a, b, delaunay_algo, centering, scatterfun, args...)
#end

end
