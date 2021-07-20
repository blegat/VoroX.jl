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
    return Foam(
        points,
        simplices,
        centers,
        voronoi_edges,
        active_edges,
    )
end

function Foam(points::Vector{SVector{N,T}}, algo::Polyhedra.Library, args...) where {N,T}
    lifted = [SVector(p..., norm(p)) for p in points]
    p = polyhedron(vrep(lifted), algo)
    simplices_h = filter(collect(Polyhedra.Indices{T,Polyhedra.halfspacetype(p)}(p))) do hi
        h = get(p, hi)
        return Polyhedra._neg(h.a[end])
    end
    simplices = Matrix{Int}(undef, N+1, length(simplices_h))
    for (i, hi) in enumerate(simplices_h)
        vertices_idx = incidentpointindices(p, hi)
        @assert length(vertices_idx) == N + 1
        for j in 1:(N+1)
            vi = vertices_idx[j]
            @assert get(p, vi) == lifted[vi.value]
            simplices[j, i] = vi.value
        end
    end
    return Foam(points, simplices, args...)
end
