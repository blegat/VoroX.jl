# Sample Set of Points in Region
function random_points(K::Int, a::SVector, b::SVector, algo)
    r = norm(b - a) / √(π * K)

    # Setup Grid
    set = InRadius(algo, [(a + b) / 2], r, a, b)

    # Now Generate New Samples by Iterating over the Set
    tries = 0

    K -= 1
    while K > 0
        if tries > length(set)
            break
        end

        # Sample Uniformly from Set
        n = rand(1:length(set))

        # Generate Sample Surrounding It
        len = r + r * rand(Float64)
        dir = randn(typeof(a))
        dir /= norm(dir)

        # New Sample Position, Grid Index
        p = set[n] + len * dir
        if !all(a .< p .< b)
            continue
        end

        if p in set
            tries += 1
        else
            push!(set, p)
            K -= 1
        end
    end
    return set.points
end