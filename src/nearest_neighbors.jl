module NN

struct GridTree{N,T}
    a::SVector{N,T}
    b::SVector{N,T}
    grid::Array{N,Int}
end
function index(g::GridTree, p::SVector)
    return Int.(size(g.grid) .* (p-g.a) ./ (g.b-g.a))
end

widen(i::Int) = (i - 2):(i + 2)
function radius(g::InRadius{<:GridTree}, p::SVector) where {N,T}
    i = index(g.inner, p)
    radius = typemax(T)
    for I in Iterators.product(widen.(i)...)
        j = g.inner.grid[I]
        if !iszero(j)
            radius = min(radius, norm(g.points[j] - p))
        end
    end
    return radius
end
function Base.push!(g::InRadius{<:GridTree}, p::SVector)
    push!(g.points, p)
    g.inner.grid[index(g.inner, p)...] = length(p.points)
    return
end

function radius(t::InRadius{<:NNTree}, p::SVector)
    return norm(nn(t.inner, x) - x)
end
function Base.push!(g::InRadius{NN}, p::SVector) where {NN<:NNTree}
    push!(g.points, p)
    g.inner = NN(g.points, Euclidean())
    return
end

mutable struct InRadius{NN, N, T}
    inner::NN
    points::Vector{SVector{N,T}}
    r::T
    function InRadius(::Type{GridTree}, data, r, a, b)
        length = r / √2
        n = ceil.((b .- a) ./ length)
        grid = zeros(Int, n...)
        obj = InRadius(GridTree(a, b, grid), data, r)
        for p in data
            push!(obj, p)
        end
        return obj
    end
    function InRadius(::Type{T}, data, r, a, b) where {T<:NNTree}
        return InRadius(T(data, Euclidean()), data, r)
    end
end

function Base.in(r::InRadius, p)
    return radius(r.inner, p) < r
end
Base.length(r::InRadius) = length(r.points)
Base.getindex(r::InRadius, i) = r.points[i]

end
