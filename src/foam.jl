# `CartesianIndex(id, simplex_id)` where
#    id::Int # which facet of the simplex from 1 to `N+1` ?
#    simplex_id::Int # id of center in `centers` of the `N`-dimensional simplex
const Facet = CartesianIndex{2}

struct Foam{N,T}
    # Delaunay nodes or Voronoi region/sites
    points::Union{PeriodicVector{N,T},Vector{SVector{N,T}}}
    # `simplices[:, id]` contains the of the `points` in the simplex `id`.
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

# Coordinate to the shift it is at
function coord_shift(a, min, max)
    return a < min ? -1 : (a > max ? 1 : 0)
end
function point_shift(point::SVector, min_coords::SVector, max_coords::SVector)
    return coord_shift.(point, min_coords, max_coords)
end

function Foam(points::Union{Vector{SVector{N,T}},PeriodicVector{N,T}}, simplices::Matrix{Int}, centering, args...) where {N,T}
    centers = SVector{N,T}[]
    for i in 1:size(simplices, 2)
        c = center(points[simplices[:, i]], centering)
        if points isa PeriodicVector
            shift = point_shift(c, args...)
            # If the centers is outside the box,
            # we consider a different shift of the simplex so that
            # the center is inside the box
            if !all(iszero, shift)
                for j in 1:size(simplices, 1)
                    id = simplices[j, i]
                    new_shift = tuple((id_shift(id, points) .- shift)...)
                    @assert all(k -> -1 <= k <= 1, new_shift)
                    simplices[j, i] = shift_range(new_shift, length(points.points))[index(points, id)]
                end
                c = center(points[simplices[:, i]], centering)
                shift = point_shift(c, args...)
                @assert all(iszero, shift)
            end
        end
        push!(centers, c)
    end
    centers = SVector{N,T}[center(points[simplices[:, i]], centering) for i in 1:size(simplices, 2)]
    voronoi_edges = zeros(Facet, size(simplices)...)
    active = Dict{Vector{Int},Facet}()
    for facet in CartesianIndices(simplices)
        I = sort(index.(Ref(points), setdiff(1:(N+1), facet[1])))
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
    @assert all(CartesianIndices(simplices)) do facet
        !iszero(knot_dist[facet]) || any(knots) do knot
            facet in knot
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

periodify(points, ::NonPeriodic) = points
periodify(points, p::Periodic) = PeriodicVector(points, p.period)

function Foam(points::Vector{SVector{N,T}}, algo, periodic::Union{NonPeriodic,Periodic}, args...) where {N,T}
    simplices = delaunay(points, algo, periodic)
    return Foam(periodify(points, periodic), convert(Matrix{Int}, simplices), args...)
end
