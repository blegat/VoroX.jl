function gradient(foam::Foam{N,T}, scale, energy, equilibration::Bool, contractive::Bool, expansive::Bool) where {N,T}
    ∇ = [zero(SVector{N,T}) for i in eachindex(foam.points)]
    for simplex in 1:size(foam.simplices, 2)
        catchment = 0.0
        for id in 1:size(foam.simplices, 1)
            facet = Facet(id, simplex)
            if iszero(foam.knot_dist[facet])
                k = foam.facet_knot[facet]
                l = length(foam.knots[k])
                catchment += (l + foam.num_catched[k]) / l
            end
        end
        catchment /= size(foam.simplices, 1)
        pidxs = foam.simplices[:,simplex]
        ps = foam.points[pidxs]
        c = center(ps, Centroid())
        for (pidx, p) in zip(pidxs, ps)
            Δ = c - p
            h = 0.0
            if equilibration
                h += 1 - scale / norm(Δ)
            end
            if contractive
                h += 1 - scale / norm(Δ) * (1 - catchment)
            end
            if expansive
                h += 1 - scale / norm(Δ) * catchment
            end
            ∇[pidx] += energy * (c - foam.points[pidx]) * h
        end
    end
    return ∇
end

function integrate(foam, dt, algo, centering, args...)
    g = gradient(foam, args...)
    points = foam.points .+ dt * g
    return Foam(points, algo, centering)
end
