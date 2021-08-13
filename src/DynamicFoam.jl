module DynamicFoam

using LinearAlgebra
using GLMakie, StaticArrays

include("center.jl")
include("nearest_neighbors.jl")
include("rand.jl")
include("foam.jl")
#include("plotting.jl")

function markersize_range(::typeof(scatter!), r)
    return LinRange(r, 100r, 20), 10r
end
function markersize_range(::typeof(meshscatter!), r)
    return LinRange(r / 1000, r / 10, 20), r / 100
end

function viz(foam, a, b, scatterfun=scatter!)
    fig = Figure()
    #s = Axis(fig[1, 1])
    s = LScene(fig[1, 1])
    r, start = markersize_range(scatterfun, norm(b - a))
    point_size = Slider(fig[2, 1], range = r, startvalue = start)
    center_size = Slider(fig[3, 1], range = r, startvalue = start / 2)
    knot_width = Slider(fig[4, 1], range = LinRange(0.1, 10, 20), startvalue = 1)
    dist_decay = Slider(fig[5, 1], range = LinRange(0, 1, 100), startvalue = 0.9)
    for simplex in 1:size(foam.simplices, 2)
        for i in 1:size(foam.simplices, 1)
            for j in 1:(i-1)
                from = foam.simplices[i, simplex]
                to = foam.simplices[j, simplex]
                a = foam.points[from]
                b = foam.points[to]
                lines!(s, [a, b])
            end
        end
    end
    for facet in CartesianIndices(foam.simplices)
        d = foam.knot_dist[facet]
        next = foam.voronoi_edges[facet]
        if d > 0 && !iszero(next)
            a = foam.centers[facet[2]]
            b = foam.centers[next[2]]
            lw = @lift($(knot_width.value) * $(dist_decay.value)^d)
            lines!(s, [a, b], color=:yellow, linewidth=lw)
        end
    end
    for knot in foam.knots
        c = foam.centers[getindex.(knot, 2)]
        push!(c, c[1])
        lines!(s, getindex.(c, 1), getindex.(c, 2), color=:orange, linewidth = knot_width.value)
    end
    scatterfun(s, foam.points, markersize = point_size.value, color=:red)
    scatterfun(s, foam.centers, markersize = center_size.value, color=:blue)
    return fig
end

function main(K, a, b, delaunay_algo, centering=Barycenter(), scatterfun=scatter!, args...)
    points = random_points(K, a, b, args...)
    foam = Foam(points, delaunay_algo, centering)
    viz(foam, a, b, scatterfun)
end

end
