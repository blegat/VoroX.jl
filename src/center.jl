abstract type AbstractSimplexCenter end

struct Circumcenter <: AbstractSimplexCenter end
function center(points::Vector{SVector{2,T}}, ::Circumcenter) where T
    @assert length(points) == 3
    a, b, c = points
    ad = a[1] * a[1] + a[2] * a[2]
    bd = b[1] * b[1] + b[2] * b[2]
    cd = c[1] * c[1] + c[2] * c[2]
    D = 2.0 * (a[1] * (b[2] - c[2]) + b[1] * (c[2] - a[2]) + c[1] * (a[2] - b[2]))
    return inv(D) * SVector(
        ad*(b[2]-c[2]) + bd*(c[2]-a[2]) + cd*(a[2]-b[2]),
        ad*(c[1]-b[1]) + bd*(a[1]-c[1]) + cd*(b[1]-a[1]),
    )
end

struct Barycenter <: AbstractSimplexCenter end
function center(points::Vector, ::Barycenter) where T
    return sum(points) / length(points)
end
