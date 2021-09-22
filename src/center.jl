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

# From https://en.wikipedia.org/wiki/Tetrahedron#Circumcenter
# Doesn't seem to work
#function center(points::Vector{SVector{3,T}}, ::Circumcenter) where T
#    @assert length(points) == 4
#    a, b, c, d = points
#    A = [(b - a)'; (c - a)'; (d - a)']::SMatrix{3,3,T}
#    display(A)
#    display(A.^2)
#    display(sum(A.^2, dims=2))
#    B = vec(sum(A.^2, dims=2))::SVector{3,T}
#    B = SVector(norm(b - a)^2, norm(c - a)^2, norm(d - a)^2)::SVector{3,T}
#    display(B)
#    return (inv(A) * B) / 2
#end

function _det(A::SMatrix{4,3,T}) where T
    o = SVector(one(T), one(T), one(T), one(T))
    return det([A o]::SMatrix{4,4,T})
end

# From pp. 6-7 of https://people.sc.fsu.edu/~jburkardt/presentations/cg_lab_tetrahedrons.pdf
function center(points::Vector{SVector{3,T}}, ::Circumcenter) where T
    @assert length(points) == 4
    a, b, c, d = points
    o = SVector(one(T), one(T), one(T), one(T))
    α = _det([a'; b'; c'; d'])
    na = sum(a.^2)
    nb = sum(b.^2)
    nc = sum(c.^2)
    nd = sum(d.^2)
    Dx = _det(@SMatrix([
        na a[2] a[3]
        nb b[2] b[3]
        nc c[2] c[3]
        nd d[2] d[3]
    ]))
    Dy = _det(@SMatrix([
        na a[1] a[3]
        nb b[1] b[3]
        nc c[1] c[3]
        nd d[1] d[3]
    ]))
    Dz = _det(@SMatrix([
        na a[1] a[2]
        nb b[1] b[2]
        nc c[1] c[2]
        nd d[1] d[2]
    ]))
    return SVector(Dx, Dy, Dz) / (2 * α)
end

struct Centroid <: AbstractSimplexCenter end
function center(points, ::Centroid) where T
    return sum(points) / length(points)
end
