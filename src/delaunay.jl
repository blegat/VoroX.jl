using Polyhedra
import VoronoiDelaunay
import MiniQhull

struct NonPeriodic end
struct Periodic{N,T}
    period::SVector{N,T}
end

function delaunay(points::Vector{SVector{N,T}}, algo::Polyhedra.Library, ::NonPeriodic) where {N,T}
    lifted = [SVector(p..., norm(p)^2) for p in points]
    p = polyhedron(vrep(lifted), algo)
    _simplices = Vector{Int}[]
    function add_Δ(vr, vv, vertices_idx, idxmap)
        @assert length(vertices_idx) == N + 1
        push!(_simplices, map(1:(N+1)) do j
            vi = vertices_idx[j]
            v = get(vr, vi)
            if v == vv[vi.value]
                i = vi.value
            else
                i = findfirst(isequal(v), vv)
            end
            idxmap[i]
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
    return simplices
end

function delaunay(points::Vector{SVector{2,T}}, algo::Type{<:VoronoiDelaunay.DelaunayTessellation2D}, ::NonPeriodic) where {T}
    tess = VoronoiDelaunay.DelaunayTessellation(length(points))
    # VoronoiDelaunay currently needs the points to be between 1 + ε and 2 - 2ε
    a = reduce((a, b) -> min.(a, b), points, init=SVector(Inf, Inf))
    b = reduce((a, b) -> max.(a, b), points, init=SVector(-Inf, -Inf))
    width = b .- a
    scaled = map(points) do p
        # Multipltiply by 1.0001 to be sure to be in [1 + ε, 2 - 2ε]
        x = (p .- a) ./ (width * (1 + 1e-4)) .+ (1 + 1e-6)
        for c in x
            if !(1 + eps(Float64) < c < 2 - 2eps(Float64))
                error("Point $p was mapped to $c for which the coordonate $x is not in the interval [$(1 + eps(Float64)), $(2 - 2eps(Float64)))]")
            end
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
    return simplices
end

function delaunay(points::Vector{SVector{N,T}}, algo::typeof(MiniQhull.delaunay), ::NonPeriodic) where {N,T}
    return MiniQhull.delaunay(points)
end

function delaunay(points::Vector{SVector{N,T}}, algo, p::Periodic{N,T}) where {N,T}
    pp = all_shift(points, p.period)
    simplices = delaunay(all_shift(points, p.period), algo, NonPeriodic())
    selected = Dict{Tuple{Vector{eltype(simplices)},Vector{NTuple{N,Int}}},Int}()
#    function score(simplex)
#        count(id_shift.(simplex, length(points), Val(N))) do shift
#            @assert all(i -> -1 <= i <= 1, shift)
#            all(iszero, shift)
#        end
#    end
    for j in 1:size(simplices, 2)
        s = simplices[:, j]
        sort!(s)
        id = mod1.(s, length(points))
        id_shifts = id_shift.(s, length(points), Val(N))
        if !any(shift -> all(iszero, shift), id_shifts)
            continue
        end
        root = id_shifts[1]
        id_shifts = map(id_shifts) do shift
            shift .- root
        end
        key = (id, id_shifts)
        if !haskey(selected, key)
            selected[key] = j
        end
    end
    return simplices[:, collect(values(selected))]
end
