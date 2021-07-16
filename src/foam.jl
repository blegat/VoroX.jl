struct Foam{N,T}
    # Delaunay nodes or Voronoi region/sites
    points::Vector{SVector{N,T}}
    # `simplices[id, :]` contains the of the `points` in the simplex `id`.
    simplices::Matrix{Int}
end

function Foam(points::Vector{SVector{N,T}}, algo::Polyhedra.Library) where {N,T}
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
    return Foam(points, simplices)
end
