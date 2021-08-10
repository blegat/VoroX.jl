using Polyhedra

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
    # Maps a `f1::Facet` the unique `f2::Facet` such that
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
end

function facet_dir(points, simplices, centers, from::Facet, to::Facet)
    c1 = centers[from[2]]
    c2 = centers[to[2]]
    if c1 ≈ c2
        # Happens for instance if the points of the
        # two simplices belong to the same hypersphere.
        #return normalize(points[simplices[to]] - points[simplices[from]])
        c1 = center(points[simplices[:, from[2]]], Barycenter())
        c2 = center(points[simplices[:, to[2]]], Barycenter())
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
                target_knot = length(knots)
            end
            for i in eachindex(knot)
                visited[knot[i]] = true
                facet_knot[knot[i]] = target_knot
                if i < start_knot
                    knot_dist[knot[i]] = start_knot - i
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
