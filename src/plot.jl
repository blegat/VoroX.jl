function markersize_range(::typeof(scatter!), r)
    return LinRange(r, 100r, 20), 10r
end
function markersize_range(::typeof(meshscatter!), r)
    return LinRange(r / 10000, r / 100, 20), r / 200
end

function _flat(points::Vector{SVector{N,T}}) where {N,T}
    #ntuple(i -> getindex.(points, i), Val(N))
    m = Matrix{T}(undef, N, length(points))
    for i in eachindex(points)
        m[:, i] = points[i]
    end
    return m
end

const HEIGHT = 1
const POINT_SIZE = 2
const CENTER_SIZE = 3
const EDGE_WIDTH = 4
const CIRCUIT_WIDTH = 5
const DECAY = 6

const WIDTH = 350

const ENERGY = 1
const SCALE = 2

function main(K, min_coords::SVector{N,T}, max_coords::SVector{N,T}, algo, centering, scatterfun=scatter!, args...) where {N,T}
    set_theme!(backgroundcolor = ColorSchemes.tableau_red_black[end])
    fig = Figure(resolution = (1200, 800))
    r, start = markersize_range(scatterfun, norm(max_coords - min_coords))

    display_sl = labelslidergrid!(
        fig,
        ["Height", "Point size", "Center size", "Edge width", "Circuit width", "Decay"],
        [100:2000, r, r, LinRange(0, 10, 20), LinRange(0.1, 10, 20), LinRange(0, 1, 100)],
        formats = [x -> "$(round(x, digits = 5))" for _ in 1:6],
        width = WIDTH,
        tellheight = true
    )

    set_close_to!(display_sl.sliders[1], 800)
    set_close_to!(display_sl.sliders[2], start)
    set_close_to!(display_sl.sliders[3], start / 2)
    set_close_to!(display_sl.sliders[4], 1)
    set_close_to!(display_sl.sliders[5], 1)
    set_close_to!(display_sl.sliders[6], 0.9)

    dynamic_sl = labelslidergrid!(
        fig,
        ["Energy", "Scale"],
        [LinRange(0, 0.01, 1000), LinRange(0, 10, 100)],
        formats = [x -> "$(round(x, digits = 5))" for _ in 1:2],
        width = WIDTH,
        tellheight = true
    )

    set_close_to!(dynamic_sl.sliders[ENERGY], 0.0005)
    set_close_to!(dynamic_sl.sliders[SCALE], 0.5)

    equilibration = Toggle(fig, active = true)
    contractive = Toggle(fig, active = true)
    expansive = Toggle(fig, active = true)
    toggles = GridLayout()
    toggles[1, 1] = Label(fig, "Equilibration")
    toggles[1, 2] = equilibration
    toggles[2, 1] = Label(fig, "Contractive")
    toggles[2, 2] = contractive
    toggles[3, 1] = Label(fig, "Expansive")
    toggles[3, 2] = expansive

    num_points = labelslider!(fig, "Number of points", N:1000, format = x -> string(x), width=WIDTH, tellheight = true)
    set_close_to!(num_points.slider, K)
    resample_button = Button(fig, label = @lift(string("Resample ", $(num_points.slider.value), " points")))

    function resample(n)
        points = random_points(num_points.slider.value[], min_coords, max_coords, args...)
        obs_foam[] = Foam(points, algo, centering)
    end

    on(resample, resample_button.clicks)

    record_frames = labelslidergrid!(
        fig,
        ["Framerate", "Length"],
        [1:30, 1:300],
        formats = [x -> "$x fps", x -> "$x frames"],
        width = WIDTH,
        tellheight = true
    )
    set_close_to!(record_frames.sliders[1], 10)
    set_close_to!(record_frames.sliders[2], 100)
    duration = Label(fig, @lift(string("Duration: ", $(record_frames.sliders[2].value) / $(record_frames.sliders[1].value), " seconds.")))
    record_button = Button(fig, label = "Record")

    t = Node(1)
    dt = 1
    time = Label(fig, @lift(string($t, "th frame.")))
    menu = GridLayout()
    menu[:v] = [time, display_sl.layout, dynamic_sl.layout, toggles, num_points.layout, resample_button, record_frames.layout, duration, record_button]

    obs_foam = Node{Foam{N,T}}()
    delaunay_edges = Node(zeros(T, 0, 0))
    out_knot_edges_obs = Dict()
    knot_obs = []
    foam_points = Node(SVector{N,T}[])
    foam_centers = Node(SVector{N,T}[])
    s = LScene(fig, height = display_sl.sliders[HEIGHT].value, scenekw = (show_axis = false,))

    function draw(foam)
        points = SVector{N,T}[]
        for simplex in 1:size(foam.simplices, 2)
            for i in 1:size(foam.simplices, 1)
                for j in 1:(i-1)
                    from = foam.simplices[i, simplex]
                    to = foam.simplices[j, simplex]
                    a = foam.points[from]
                    b = foam.points[to]
                    push!(points, a)
                    push!(points, b)
                end
            end
        end
        delaunay_edges[] = _flat(points)

        dists = Int[]
        out_knot_edges = Dict()
        for facet in CartesianIndices(foam.simplices)
            d = foam.knot_dist[facet]
            next = foam.voronoi_edges[facet]
            if d > 0 && !iszero(next)
                a = foam.centers[facet[2]]
                b = foam.centers[next[2]]
                if haskey(out_knot_edges, d)
                    push!(out_knot_edges[d], a)
                    push!(out_knot_edges[d], b)
                else
                    out_knot_edges[d] = [a, b]
                end
            end
        end
        for d in keys(out_knot_edges_obs)
            if haskey(out_knot_edges, d)
                out_knot_edges_obs[d][] = _flat(out_knot_edges[d])
                delete!(out_knot_edges, d)
            else
                out_knot_edges_obs[d][] = zeros(T, N, 0)
            end
        end
        for d in keys(out_knot_edges)
            @assert !haskey(out_knot_edges_obs, d)
            lw = @lift($(display_sl.sliders[CIRCUIT_WIDTH].value) * $(display_sl.sliders[DECAY].value)^d)
            col = @lift(ColorSchemes.haline[$(display_sl.sliders[DECAY].value)^d])
            mobs = Node(_flat(out_knot_edges[d]))
            linesegments!(s, mobs, color=col, linewidth=lw)
            out_knot_edges_obs[d] = mobs
        end

        i = 0
        for knot in foam.knots
            for facet in knot
                @assert iszero(foam.knot_dist[facet])
            end
            c = foam.centers[getindex.(knot, 2)]
            push!(c, c[1])
            m = _flat(c)
            i += 1
            if i > length(knot_obs)
                obs = Node(m)
                lines!(s, obs, color=ColorSchemes.speed[0.5], linewidth = display_sl.sliders[CIRCUIT_WIDTH].value)
                push!(knot_obs, obs)
            else
                knot_obs[i][] = m
            end
        end
        for j in (i+1):length(knot_obs)
            knot_obs[j][] = zeros(T, N, 0)
        end

        foam_points[] = foam.points
        foam_centers[] = foam.centers
    end

    fig[1, 1] = menu
    fig[1, 2] = s
    on(draw, obs_foam)
    function integ()
        t[] += dt
        obs_foam[] = integrate(obs_foam[], dt, algo, centering, dynamic_sl.sliders[SCALE].value[], dynamic_sl.sliders[ENERGY].value[], equilibration.active[], contractive.active[], expansive.active[])
    end
    on(events(fig).keyboardbutton) do event
        if event.action in (Keyboard.press, Keyboard.repeat) &&
            (event.key == Keyboard.right || event.key == Keyboard.down)
            integ()
        end
        # Let the event reach other listeners
        return Consume(false)
    end

    on(record_button.clicks) do n
        i = 1
        while isfile("recording_$i.mp4")
            i += 1
        end
        n = record_frames.sliders[2].value[]
        p = ProgressMeter.Progress(n, 1)
        path = joinpath(pwd(), "recording_$i.mp4")
        @info("Starting recording of $n frames to be recorded in `$path`.")
        record(fig, path, 1:n; framerate = record_frames.sliders[1].value[]) do j
            integ()
            ProgressMeter.next!(p)
        end
        @info("Recording completed")
    end

    resample(0)
    linesegments!(s, delaunay_edges, linewidth = display_sl.sliders[EDGE_WIDTH].value)
    scatterfun(s, foam_points, markersize = display_sl.sliders[POINT_SIZE].value, color=ColorSchemes.inferno[0.5])
    scatterfun(s, foam_centers, markersize = display_sl.sliders[CENTER_SIZE].value, color=ColorSchemes.roma[0.65])


    trim!(fig.layout)
    return fig
end
