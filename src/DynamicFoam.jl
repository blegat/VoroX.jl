module DynamicFoam

using LinearAlgebra
using GLMakie, StaticArrays

include("center.jl")
include("nearest_neighbors.jl")
include("rand.jl")
include("foam.jl")
include("motion.jl")
#include("plotting.jl")

function markersize_range(::typeof(scatter!), r)
    return LinRange(r, 100r, 20), 10r
end
function markersize_range(::typeof(meshscatter!), r)
    return LinRange(r / 1000, r / 10, 20), r / 100
end

function _flat(points::Vector{SVector{N,T}}) where {N,T}
    #ntuple(i -> getindex.(points, i), Val(N))
    m = Matrix{T}(undef, N, length(points))
    for i in eachindex(points)
        m[:, i] = points[i]
    end
    return m
end

function viz(foam::Foam{N,T}, a, b, algo, centering, scatterfun=scatter!) where {N,T}
    fig = Figure(resolution = (1200, 800))
    #s = Axis(fig[1, 1])
    height = Slider(fig, range = 100:1000, startvalue = 400)
    r, start = markersize_range(scatterfun, norm(b - a))
    point_size = Slider(fig, range = r, startvalue = start)
    center_size = Slider(fig, range = r, startvalue = start / 2)
    knot_width = Slider(fig, range = LinRange(0.1, 10, 20), startvalue = 1)
    dist_decay = Slider(fig, range = LinRange(0, 1, 100), startvalue = 0.9)
    energy = Slider(fig, range = LinRange(0, 10, 100), startvalue = 0.1)
    scale = Slider(fig, range = LinRange(0, 10, 100), startvalue = 0.1)
    equilibration = Toggle(fig, active = true)
    contractive = Toggle(fig, active = true)
    expansive = Toggle(fig, active = true)
    toggles = GridLayout()
    toggles[:h] = [equilibration, contractive, expansive]
    menu = GridLayout()
    menu[:v] = [height, point_size, center_size, knot_width, dist_decay, energy, scale, toggles]
    obs_foam = Node(foam)
    t = Node(0)
    dt = 1
    delaunay_edges = Node(zeros(T, 0, 0))
    function draw(foam)
    end
    s = LScene(fig, height = height.value)
    text!(s, @lift(string($t)), position = a + 0.99 * (b - a))
    points = SVector{N,T}[]
    for simplex in 1:size(foam.simplices, 2)
        for i in 1:size(foam.simplices, 1)
            for j in 1:(i-1)
                from = foam.simplices[i, simplex]
                to = foam.simplices[j, simplex]
                a = foam.points[from]
                b = foam.points[to]
                push!(points, a)
                push!(points, b)
            end
        end
    end
    delaunay_edges[] = _flat(points)
    empty!(points)
    linesegments!(s, delaunay_edges)
    dists = Int[]
    out_knot_edges = Dict()
    #out_knot_edges = Node(zeros(T, 0, 0))
    #dists = Node(Int[])
    for facet in CartesianIndices(foam.simplices)
        d = foam.knot_dist[facet]
        next = foam.voronoi_edges[facet]
        if d > 0 && !iszero(next)
            a = foam.centers[facet[2]]
            b = foam.centers[next[2]]
            #push!(out_knot_lw, @lift($(knot_width.value) * $(dist_decay.value)^d))
            #push!(out_knot_lw, knot_width.value[] * dist_decay.value[]^d)
            push!(points, a)
            push!(points, b)
            if haskey(out_knot_edges, d)
                push!(out_knot_edges[d], a)
                push!(out_knot_edges[d], b)
            else
                out_knot_edges[d] = [a, b]
            end
        end
    end
    out_knot_edges_obs = Dict()
    for (d, v) in out_knot_edges
        m = _flat(out_knot_edges[d])
        if haskey(out_knot_edges_obs, d)
            out_knot_edges_obs[d][] = m
        else
            lw = @lift($(knot_width.value) * $(dist_decay.value)^d)
            mobs = Node(m)
            linesegments!(s, mobs, color=:yellow, linewidth=lw)
            out_knot_edges_obs[d] = mobs
        end
    end
    #@out_knot_edges_obs[] = _flat(points)
    #out_knot_lw = @lift($(knot_width.value) * $(dist_decay.value).^$dists)
    #@show length(out_knot_lw)
    #linesegments!(s, out_knot_edges)
    #linesegments!(s, out_knot_edges, color=:yellow, linewidth=out_knot_lw)
    for knot in foam.knots
        c = foam.centers[getindex.(knot, 2)]
        push!(c, c[1])
        lines!(s, getindex.(c, 1), getindex.(c, 2), color=:orange, linewidth = knot_width.value)
    end
    scatterfun(s, foam.points, markersize = point_size.value, color=:red)
    scatterfun(s, foam.centers, markersize = center_size.value, color=:blue)
    fig[1, 1] = menu
    fig[1, 2] = s
    draw(foam)
    on(draw, obs_foam)
    on(events(fig).keyboardbutton) do event
        if event.action in (Keyboard.press, Keyboard.repeat) &&
            (event.key == Keyboard.right || event.key == Keyboard.down)
            t[] += dt
            obs_foam[] = integrate(obs_foam[], dt, algo, centering, scale.value[], energy.value[], equilibration.active[], contractive.active[], expansive.active[])
        end
        # Let the event reach other listeners
        return Consume(false)
    end
    trim!(fig.layout)
    return fig
end

function main(K, a, b, delaunay_algo, centering=Barycenter(), scatterfun=scatter!, args...)
    points = random_points(K, a, b, args...)
    foam = Foam(points, delaunay_algo, centering)
    return viz(foam, a, b, delaunay_algo, centering, scatterfun)
end

end
