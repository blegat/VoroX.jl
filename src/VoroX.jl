module VoroX

export foam

using LinearAlgebra
using GLMakie, StaticArrays
import ProgressMeter
import ColorSchemes

include("center.jl")
include("nearest_neighbors.jl")
include("rand.jl")
include("periodic.jl")
include("delaunay.jl")
include("foam.jl")
include("motion.jl")
include("plot.jl")

import MiniQhull, QHull, CDDLib, VoronoiDelaunay
const LIBRARIES = [MiniQhull.delaunay, QHull.Library(), CDDLib.Library(:float), VoronoiDelaunay.DelaunayTessellation2D]

end
