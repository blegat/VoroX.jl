function markersize_range(::typeof(scatter!), r)
    return LinRange(r / 10, 40r, 40), 1r
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

function shrink(v::VRepresentation, scaling)
    c = center(points(v), Centroid())
    return vrep(map(points(v)) do p
        c + (p - c) * scaling
    end)
end

const HEIGHT = 1
const POINT_SIZE = 2
const CENTER_SIZE = 3
const EDGE_WIDTH = 4
const CIRCUIT_WIDTH = 5
const DECAY = 6
const TRANSPARENCY = 7

const WIDTH = 350

const ENERGY = 1
const SCALE = 2

function clip_point(p, dir, min_coords, max_coords::SVector{N}) where {N}
    for i in 1:N
        if p[i] < min_coords[i]
            p = p + dir * (min_coords[i] - p[i]) / dir[i]
        end
        if p[i] > max_coords[i]
            p = p + dir * (max_coords[i] - p[i]) / dir[i]
        end
    end
    return p
end

function _shift(f, p::SVector{N}, min_coords::SVector{N}, max_coords::SVector{N}) where {N}
    return ntuple(Val(N)) do i
        if f(p[i], min_coords[i])
            return 1
        elseif f(max_coords[i], p[i])
            return -1
        else
            return 0
        end
    end
end

_points(v::Vector) = v
_points(v::PeriodicVector) = v.points

function clipped_edge(from, to, min_coords, max_coords::SVector{N}) where {N}
    shift = _shift(<, from, min_coords, max_coords)
    period = max_coords - min_coords
    from = shift_point(from, shift, period)
    to = shift_point(to, shift, period)
    out = [from]
    for i in 1:3
        if all(iszero, _shift(<, to, min_coords, max_coords))
            break
        end
        mid = clip_point(to, from - to, min_coords, max_coords)
        push!(out, mid)
        shift = _shift(≈, mid, min_coords, max_coords)
        from = shift_point(mid, shift, period)
        to = shift_point(to, shift, period)
        push!(out, from)
    end
    if !all(iszero, _shift(<, to, min_coords, max_coords))
        error("$from, $to")
    end
    push!(out, to)
    return out
end

function main(K, min_coords::SVector{N,T}, max_coords::SVector{N,T}, scatterfun=scatter!, args...) where {N,T}
    _periodic(is_periodic) = is_periodic ? Periodic(max_coords - min_coords) : NonPeriodic()
    set_theme!(backgroundcolor = ColorSchemes.tableau_red_black[end])
    fig = Figure(resolution = (1200, 1200))
    r, start = markersize_range(scatterfun, norm(max_coords - min_coords))

    display_sl = labelslidergrid!(
        fig,
        ["Height", "Point size", "Center size", "Edge width", "Circuit width", "Decay", "Transparency"],
        [100:2000, r, r, LinRange(0, 1, 20), LinRange(0.1, 10, 20), LinRange(0, 1, 100), LinRange(0.0, 1.0, 21)],
        formats = [[x -> "$(round(x, digits = 5))" for _ in 1:6]; [x -> "$(round(x, digits = 2))"]],
        width = WIDTH,
        tellheight = true
    )

    set_close_to!(display_sl.sliders[1], 800)
    set_close_to!(display_sl.sliders[2], start)
    set_close_to!(display_sl.sliders[3], start / 2)
    set_close_to!(display_sl.sliders[4], 0.5)
    set_close_to!(display_sl.sliders[5], 1)
    set_close_to!(display_sl.sliders[6], 0.9)
    set_close_to!(display_sl.sliders[TRANSPARENCY], 0.4)

    dynamic_sl = labelslidergrid!(
        fig,
        ["Energy", "Scale"],
        [LinRange(0, 0.01, 1000), LinRange(0, 1, 100)],
        formats = [x -> "$(round(x, digits = 5))" for _ in 1:2],
        width = WIDTH,
        tellheight = true
    )

    set_close_to!(dynamic_sl.sliders[ENERGY], 0.0005)
    set_close_to!(dynamic_sl.sliders[SCALE], 0.5)

    edge_scale = Toggle(fig, active = true)
    equilibration = Toggle(fig, active = true)
    contractive = Toggle(fig, active = false)
    expansive = Toggle(fig, active = true)
    periodic = Toggle(fig, active = true)
    toggles = GridLayout()
    toggles[1, 1] = Label(fig, "Edge scale")
    toggles[1, 2] = edge_scale
    toggles[1, 3] = Label(fig, "Equilibration")
    toggles[1, 4] = equilibration
    toggles[2, 1] = Label(fig, "Contractive")
    toggles[2, 2] = contractive
    toggles[2, 3] = Label(fig, "Expansive")
    toggles[2, 4] = expansive
    toggles[3, 1] = Label(fig, "Periodic")
    toggles[3, 2] = periodic
    on(periodic.active) do a
        obs_foam[] = Foam(_points(obs_foam[].points), current_library(), _periodic(a), current_centering())
    end

    lib_layout = GridLayout()
    lib = Menu(fig, options = ["MiniQhull", "Qhull", "CDDLib", "VoronoiDelaunay"])
    lib.i_selected = 1
    lib_layout[:h] = [Label(fig, "Library"), lib]
    current_library() = LIBRARIES[lib.i_selected[]]
    on(lib.selection) do s
        obs_foam[] = Foam(_points(obs_foam[].points), current_library(), _periodic(periodic.active[]), current_centering())
    end

    centering_layout = GridLayout()
    centering = Menu(fig, options = ["Centroid", "Circumcenter"])
    centering.i_selected = 1
    centering_layout[:h] = [Label(fig, "Centering"), centering]
    function current_centering()
        if centering.i_selected[] == 1
            return Centroid()
        else
            @assert centering.i_selected[] == 2
            return Circumcenter()
        end
    end
    on(centering.selection) do s
        obs_foam[] = Foam(_points(obs_foam[].points), current_library(), _periodic(periodic.active[]), current_centering())
    end

    num_points = labelslider!(fig, "Number of points", N:1000, format = x -> string(x), width=WIDTH, tellheight = true)
    set_close_to!(num_points.slider, K)
    resample_button = Button(fig, label = @lift(string("Resample ", $(num_points.slider.value), " points")))

    function resample(n)
        points = random_points(num_points.slider.value[], min_coords, max_coords, args...)
        t[] = 1
        obs_foam[] = Foam(points, current_library(), _periodic(periodic.active[]), current_centering())
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
    duration_fps = GridLayout()
    duration = Label(fig, @lift(string("Duration: ", round($(record_frames.sliders[2].value) / $(record_frames.sliders[1].value), digits=2), " seconds.")))
    #cur_fps = Node(NaN)
    #cur_fps_label = Label(fig, @lift(isnan($cur_fps) ? "" : string(round($cur_fps, digits=2), " fps.")))
    #duration_fps[:h] = [duration, cur_fps_label]
    record_layout = GridLayout()
    record_button = Button(fig, label = "⬤", labelcolor="red")
    live_toggle = Toggle(fig, active = false)
    record_layout[:h] = [record_button, Label(fig, "Live"), live_toggle]

    t = Node(1)
    dt = 1
    last_frame_time = Node(NaN)
    time = Label(fig, @lift(string($t, "th frame.", isnan($last_frame_time) ? "" : string(" ", round(inv($last_frame_time), digits=2), " fps."))))
    menu = GridLayout()
    menu[:v] = [time, display_sl.layout, dynamic_sl.layout, toggles, lib_layout, centering_layout, num_points.layout, resample_button, record_frames.layout, duration, record_layout]

    obs_foam = Node{Foam{N,T}}()
    delaunay_edges = Node(zeros(T, 0, 0))
    out_knot_edges_obs = Dict()
    knot_obs = []
    foam_points = Node(SVector{N,T}[])
    foam_centers = Node(SVector{N,T}[])
    foam_cells = []
    s = LScene(fig, height = display_sl.sliders[HEIGHT].value, scenekw = (show_axis = false,))

    function draw(foam)
        delaunay_edges_points = SVector{N,T}[]
        for simplex in 1:size(foam.simplices, 2)
            for i in 1:size(foam.simplices, 1)
                for j in 1:(i-1)
                    from = foam.simplices[i, simplex]
                    to = foam.simplices[j, simplex]
                    append!(delaunay_edges_points, clipped_edge(foam.points[from], foam.points[to], min_coords, max_coords))
                end
            end
        end
        delaunay_edges[] = _flat(delaunay_edges_points)

        dists = Int[]
        out_knot_edges = Dict()
        for facet in CartesianIndices(foam.simplices)
            d = foam.knot_dist[facet]
            next = foam.voronoi_edges[facet]
            if d > 0 && !iszero(next)
                append!(
                    get!(out_knot_edges, d, SVector{N,T}[]),
                    clipped_edge(foam.centers[facet[2]], foam.centers[next[2]], min_coords, max_coords),
                )
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
            d = SVector{N,T}[]
            for i in eachindex(c)
                append!(d, clipped_edge(c[i], c[mod1(i + 1, length(c))], min_coords, max_coords))
            end
            m = _flat(d)
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

        foam_points[] = _points(foam.points)
        foam_centers[] = foam.centers

        for cell in foam_cells
            delete!(s, cell)
        end
        empty!(foam_cells)
        if display_sl.sliders[TRANSPARENCY].value[] > 0
            if foam.points isa PeriodicVector
                cube = reduce(*, [HalfSpace(SVector(one(T)), max_coords[i]) ∩ HalfSpace(SVector(-one(T)), -min_coords[i]) for i in 1:N])
            else
                cube = nothing
            end
            cell_points = Tuple{Vector{SVector{N,T}},Float64}[]
            for i in 1:size(foam.simplices, 2)
                simplex = foam.simplices[:, i]
                ps = foam.points[simplex]
                catchment = simplex_catchment(foam, i)
                push!(cell_points, (ps, catchment))
                if foam.points isa PeriodicVector
                    period = max_coords - min_coords
                    for shift in Iterators.product(ntuple(_ -> -1:1, Val(N))...)
                        push!(cell_points, ([shift_point(p, map(-, shift), period) for p in ps], catchment))
                    end
                end
            end
            cell_polyhedra = map(cell_points) do (cell_p, catchment)
                vr = vrep(cell_p)
                #p = polyhedron(shrink(vr, 0.8))
                p = polyhedron(vr)
                if cube !== nothing
                    p = p ∩ cube
                end
                return p, catchment
            end
            if foam.points isa PeriodicVector
                filter!(cell_polyhedra) do (p, _)
                    Polyhedra.computevrep!(p)
                    !isempty(p)
                end
            end
            #volumes = Polyhedra.volume.(cell_polyhedra)
            #min_volume, max_volume = extrema(volumes)
            for i in eachindex(cell_polyhedra)
                cell_polyhedron, catchment = cell_polyhedra[i]
                ratio = sqrt(min(catchment, 1))
                #ratio = (volumes[i] - min_volume) / (max_volume - min_volume)
                #scheme = ColorSchemes.hsv
                scheme = ColorSchemes.coolwarm
                push!(foam_cells, mesh!(s, Polyhedra.Mesh(cell_polyhedron), color = (scheme[ratio], @lift($(display_sl.sliders[TRANSPARENCY].value) * (1 - ratio) + ratio))))
            end
        end
    end

    fig[1, 1] = menu
    fig[1, 2] = s
    on(draw, obs_foam)
    function _integ()
        t[] += dt
        obs_foam[] = integrate(obs_foam[], dt, current_library(), _periodic(periodic.active[]), current_centering(), edge_scale.active[], dynamic_sl.sliders[SCALE].value[], dynamic_sl.sliders[ENERGY].value[], equilibration.active[], contractive.active[], expansive.active[])
    end
    function integ()
        cur_time = @elapsed _integ()
        @info("New frame computed in $cur_time. Possible fps: $(inv(cur_time))")
        last_frame_time[] = cur_time
        return cur_time
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

    live = true
    running = false
    on(live_toggle.active) do active
        if active
            live = true
            if running
                @info("Live animation already running")
            else
                running = true
                @async begin
                    @info("Starting live animation")
                    while live
                        cur_time = integ()
                        target_fps = record_frames.sliders[1].value[]
                        target_time = inv(target_fps)
                        if target_time > cur_time
                            @info("Sleeping for $(target_time - cur_time).")
                            sleep(target_time - cur_time)
                        end
                    end
                    @info("Stopping live animation")
                    running = false
                end
            end
        else
            live = false
        end
    end

    return fig
end


function foam(K::Int, N::Int, args...)
    main(K, SVector(ntuple(_ -> -1.0, Val(N))), SVector(ntuple(_ -> 1.0, Val(N))))
end
