struct PeriodicVector{N,T}
    points::Vector{SVector{N,T}}
    period::SVector{N,T}
end

index(::Vector, i) = i
index(v::PeriodicVector, i) = mod1(i, length(v.points))
function Base.getindex(v::PeriodicVector{N}, i::Integer) where {N}
    return shift_point(v.points[index(v, i)], id_shift(i, length(v.points), Val(N)), v.period)
end
Base.getindex(v::PeriodicVector, i::Vector) = getindex.(Ref(v), i)

function shift_point(point::SVector{N,T}, shift::NTuple{N,Int}, period::SVector{N,T}) where {N,T}
    return point .+ shift .* period
end

function all_shift(points::Vector{SVector{N,T}}, p::SVector{N,T}) where {N,T}
    n = length(points)
    ps = Vector{SVector{N,T}}(undef, 3^N * length(points))
    for shift in Iterators.product(ntuple(_ -> -1:1, Val(3))...)
        ps[shift_range(shift, n)] = shift_point.(points, Ref(shift), Ref(p))
    end
    return ps
end

shift_offset(shift::NTuple{0,Int}) = 0
function shift_offset(shift::NTuple{N,Int}) where {N}
    return 1 + first(shift) + 3shift_offset(Base.tail(shift))
end
shift_range(shift, len) = len * shift_offset(shift) .+ (1:len)

id_shift(id, ::Tuple{}) = tuple()
function id_shift(id, t::Tuple)
    return ((id % 3) - 1), id_shift(div(id, 3), Base.tail(t))...
end
function id_shift(id, n, l::Val)
    return id_shift(div(id - 1, n), ntuple(_ -> nothing, l))
end
