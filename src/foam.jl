using Polyhedra
import VoronoiDelaunay

# `CartesianIndex(id, simplex_id)` where
#    id::Int # which facet of the simplex from 1 to `N+1` ?
#    simplex_id::Int # id of center in `centers` of the `N`-dimensional simplex
const Facet = CartesianIndex{2}

struct Foam{N,T}
    # Delaunay nodes or Voronoi region/sites
    points::Vector{SVector{N,T}}
    # `simplices[id, :]` contains the of the `points` in the simplex `id`.
    simplices::Matrix{Int}
    # Delaunay centers or Voronoi vertices
    centers::Vector{SVector{N,T}}
    # Maps a `f1::Facet` to the unique `f2::Facet` such that
    # their points are the same.
    voronoi_edges::Matrix{Facet}
    # Maps a `f1::Facet` such that `f2::Facet` is taken after `f1`
    # in the knots dynamic.
    active_edges::Matrix{Facet}
    # The cycles in the union of cycles formed by `active_edges`.
    knots::Vector{Vector{Facet}}
    # Maps a facet to the index in `knots` of the knot the dynamics leads it to.
    facet_knot::Matrix{Int}
    # Maps a facet to the length of the path util it reaches the knot or zero
    # if it is part of the knot.
    knot_dist::Matrix{Int}
    # Number of facets out of the knot that converge to the knot.
    num_catched::Vector{Int}
end

function facet_dir(points, simplices, centers, from::Facet, to::Facet)
    c1 = centers[from[2]]
    c2 = centers[to[2]]
    if c1 ≈ c2
        # Happens for instance if the points of the
        # two simplices belong to the same hypersphere.
        #return normalize(points[simplices[to]] - points[simplices[from]])
        c1 = center(points[simplices[:, from[2]]], Centroid())
        c2 = center(points[simplices[:, to[2]]], Centroid())
    end
    return normalize(c2 - c1)
end

function Foam(points::Vector{SVector{N,T}}, simplices::Matrix{Int}, centering) where {N,T}
    centers = SVector{N,T}[center(points[simplices[:, i]], centering) for i in 1:size(simplices, 2)]
    voronoi_edges = zeros(Facet, size(simplices)...)
    active = Dict{Vector{Int},Facet}()
    for facet in CartesianIndices(simplices)
        I = sort(setdiff(1:(N+1), facet[1]))
        pidx = simplices[I, facet[2]]
        if haskey(active, pidx)
            facet2 = active[pidx]
            @assert iszero(voronoi_edges[facet2])
            voronoi_edges[facet] = facet2
            voronoi_edges[facet2] = facet
        else
            active[pidx] = facet
        end
    end
    active_edges = map(CartesianIndices(simplices)) do facet
        mirror = voronoi_edges[facet]
        best = zero(Facet)
        best_dot = -one(T)
        if iszero(mirror)
            return best
        end
        simplex = mirror[2]
        dir = facet_dir(points, simplices, centers, facet, mirror)
        for i in 1:(N+1)
            next = Facet(i, simplex)
            if next != mirror
                next_next = voronoi_edges[next]
                if !iszero(next_next)
                    next_dir = facet_dir(points, simplices, centers, mirror, next_next)
                    cur_dot = dot(dir, next_dir)
                    if cur_dot > best_dot
                        best = next
                        best_dot = cur_dot
                    end
                end
            end
        end
        return best
    end
    visited = falses(size(simplices)...)
    facet_knot = zeros(Int, size(simplices)...)
    knot_dist = zeros(Int, size(simplices)...)
    knots = Vector{Facet}[]
    num_catched = Int[]
    for facet in CartesianIndices(simplices)
        if !visited[facet]
            knot = [facet]
            knot_index = Dict(facet => 1)
            cur = active_edges[facet]
            while !iszero(cur) && !haskey(knot_index, cur) && !visited[cur]
                push!(knot, cur)
                knot_index[cur] = length(knot)
                cur = active_edges[cur]
            end
            if iszero(cur)
                start_knot = length(knot) + 1
                target_knot = 0
            elseif visited[cur]
                start_knot = length(knot) + 1 + knot_dist[cur]
                target_knot = facet_knot[cur]
            else
                start_knot = knot_index[cur]
                push!(knots, knot[start_knot:end])
                push!(num_catched, 0)
                target_knot = length(knots)
            end
            for i in eachindex(knot)
                visited[knot[i]] = true
                facet_knot[knot[i]] = target_knot
                if i < start_knot
                    knot_dist[knot[i]] = start_knot - i
                    if !iszero(target_knot)
                        num_catched[target_knot] += 1
                    end
                end
            end
        end
    end
    return Foam(
        points,
        simplices,
        centers,
        voronoi_edges,
        active_edges,
        knots,
        facet_knot,
        knot_dist,
        num_catched,
    )
end

function Foam(points::Vector{SVector{N,T}}, algo::Polyhedra.Library, args...) where {N,T}
    lifted = [SVector(p..., norm(p)^2) for p in points]
    p = polyhedron(vrep(lifted), algo)
    _simplices = Vector{Int}[]
    function add_Δ(vr, vv, vertices_idx, idxmap)
        @assert length(vertices_idx) == N + 1
        push!(_simplices, map(1:(N+1)) do j
            vi = vertices_idx[j]
            @assert get(vr, vi) == vv[vi.value]
            idxmap[vi.value]
        end)
    end
    for hi in Polyhedra.Indices{T,Polyhedra.halfspacetype(p)}(p)
        h = get(p, hi)
        if Polyhedra._neg(h.a[end])
            vertices_idx = incidentpointindices(p, hi)
            @assert length(vertices_idx) >= N + 1
            if length(vertices_idx) == N + 1
                add_Δ(p, lifted, vertices_idx, eachindex(lifted))
            else
                idx = map(vertices_idx) do vi
                    @assert get(p, vi) == lifted[vi.value]
                    vi.value
                end
                vv = points[idx]
                vr = polyhedron(vrep(vv), algo)
                for Δ in triangulation_indices(vr)
                    add_Δ(vr, vv, Δ, idx)
                end
            end
        end
    end
    simplices = Matrix{Int}(undef, N+1, length(_simplices))
    for (i, Δ) in enumerate(_simplices)
        simplices[:, i] = Δ
    end
    return Foam(points, simplices, args...)
end

function Foam(points::Vector{SVector{2,T}}, algo::Type{<:VoronoiDelaunay.DelaunayTessellation2D}, args...) where {T}
    tess = VoronoiDelaunay.DelaunayTessellation(length(points))
    # VoronoiDelaunay currently needs the points to be between 1 + ε and 2 - 2ε
    a = reduce((a, b) -> min.(a, b), points, init=SVector(Inf, Inf))
    b = reduce((a, b) -> max.(a, b), points, init=SVector(-Inf, -Inf))
    width = b .- a
    scaled = map(points) do p
        # Multipltiply by 1.0001 to be sure to be in [1 + ε, 2 - 2ε]
        x = (p .- a) ./ (width * (1 + 1e-4)) .+ (1 + 1e-6)
        @assert all(x) do c
            1 + eps(Float64) < c < 2 - 2eps(Float64)
        end
        return VoronoiDelaunay.Point2D(x...)
    end
    # Should construct before as `tess` modifies `scaled`.
    back = Dict(p => i for (i, p) in enumerate(scaled))
    push!(tess, scaled)
    Δs = VoronoiDelaunay.DelaunayTriangle{VoronoiDelaunay.GeometricalPredicates.Point2D}[]
    # Cannot collect as it does not implement length nor IteratorSize
    for Δ in tess
        push!(Δs, Δ)
    end
    simplices = Matrix{Int}(undef, 3, length(Δs))
    for (i, Δ) in enumerate(Δs)
        simplices[1, i] = back[VoronoiDelaunay.geta(Δ)]
        simplices[2, i] = back[VoronoiDelaunay.getb(Δ)]
        simplices[3, i] = back[VoronoiDelaunay.getc(Δ)]
    end
    return Foam(points, simplices, args...)
end
