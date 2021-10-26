function homothety(Δ, catchment, scale, energy, equilibration::Bool, contractive::Bool, expansive::Bool)
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
    return energy * Δ * h
end
function simplex_catchment(foam, simplex)
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
    return catchment
end
function gradient(foam::Foam{N,T}, edge_scale::Bool, args...) where {N,T}
    ∇ = [zero(SVector{N,T}) for i in eachindex(_points(foam.points))]
    for simplex in 1:size(foam.simplices, 2)
        catchment = simplex_catchment(foam, simplex)
        pidxs = index.(Ref(foam.points), foam.simplices[:,simplex])
        ps = foam.points[pidxs]
        c = center(ps, Centroid())
        for (pidx, p) in zip(pidxs, ps)
            if edge_scale
                # Similar to https://github.com/weigert/DynamicFoam
                for (qidx, q) in zip(pidxs, ps)
                    if qidx != pidx
                        ∇[pidx] += homothety(q - p, catchment, args...)
                    end
                end
            else
                ∇[pidx] += homothety(c - p, catchment, args...)
            end
        end
    end
    return ∇
end

function integrate(foam, dt, algo, periodic, centering, args...)
    g = gradient(foam, args...)
    points = _points(foam.points) .+ dt * g
    return Foam(points, algo, periodic, centering)
end
