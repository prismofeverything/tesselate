module tesselate

using Random
using Distributions

function each(m)
    x, y = axes(m)
    eachindex(view(m, x, y))
end

# grid

unit_size = 40

function generate_world(bounds)
    Dict(
        "states" => ones(Int8, bounds[1], bounds[2]),
        "links" => Vector{Matrix{Int8}}(),
        "from" => Dict{Vector{Int8}, Vector{Vector{Int8}}}()
    )
end

function world_size(world)
    size(world["states"])
end

function set_state(world, location, state)
    world["states"][location[1], location[2]] = state
    world
end

function get_state(world, location)
    world["states"][location[1], location[2]]
end

function state_indexes(world)
    each(world["states"])
end

function has_link(world, a, b)
    if haskey(world["from"], a)
        b in world["from"][a]
    else
        false
    end
end

function add_from(world, a, b)
    if !haskey(world["from"], a)
        world["from"][a] = []
    end
    push!(world["from"][a], b)
end

function add_link(world, a, b)
    if !has_link(world, a, b)
        link = reduce(hcat, [a, b])
        push!(world["links"], link)
        add_from(world, a, b)
        add_from(world, b, a)
    end
    world
end

function link_count(world, a)
    if haskey(world["from"], a)
        size(world["from"][a])
    else
        0
    end
end

EMPTY = 0
SUBSTRATE = 1
MEMBRANE = 2
ENZYME = 3
REPAIR = 4

symbols = Dict(
    EMPTY => " ",
    SUBSTRATE => "◯",
    MEMBRANE => "⧇",
    ENZYME => "*",
    REPAIR => "△"
)

colors = Dict(
    EMPTY => [0, 0, 0],
    SUBSTRATE => [0.6, 0.6, 0.6],
    MEMBRANE => [0.8, 0.7, 0.1],
    ENZYME => [0.7, 0.1, 0.6],
    REPAIR => [0.1, 0.7, 0.4]
)

function random_location(bounds)
    [rand(1:bounds[1]), rand(1:bounds[2])]
end

function initialize_world(bounds, counts)
    world = generate_world(bounds)
    for (state, count) in counts
        for n in 1:count
            location = random_location(bounds)
            while get_state(world, location) != 1
                println(location)
                location = random_location(bounds)
            end
            set_state(world, location, state)
        end
    end
    world
end

# dynamics

# O -> A
# [] -> B
# * -> F
# △ -> Phi
# all we are missing is a △ to "repair" our * from [], and also to have [] act on * to produce △

function mod_index(x, d)
    r = mod1(x, d)
    if r == 0
        r = d
    end
    r
end

function mod_location(bounds, location)
    [
        mod_index(location[1], bounds[1]),
        mod_index(location[2], bounds[2]),
    ]
end

function adjacent_locations(world, location)
    bounds = world_size(world)
    mod_bounds = l -> mod_location(bounds, l)
    adjacent = [
        [location[1], location[2]+1],
        [location[1], location[2]-1],
        [location[1]+1, location[2]],
        [location[1]-1, location[2]]
    ]

    map(mod_bounds, adjacent)
end

function no_action(world, location)
    world
end

function move_element(world, location)
    bounds = world_size(world)
    state = get_state(world, location)
    adjacent = adjacent_locations(world, location)
    order = shuffle(adjacent)
    for space in order
        if get_state(world, space) == EMPTY
            world = set_state(world, location, EMPTY)
            world = set_state(world, space, state)
            break
        elseif state == SUBSTRATE && get_state(world, space) == MEMBRANE
            direction = space - location
            beyond = mod_location(bounds, space + (direction * 2))
            if get_state(world, beyond) == EMPTY
                world = set_state(world, location, EMPTY)
                world = set_state(world, beyond, state)
                break
            end
        end
    end
    world
end

function move_membrane(world, location)
    if link_count(world, location) == 0
        world = move_element(world, location)
    end
    world
end

function consolidate_actions(dynamics)
    consolidate = Dict()

    for (state, actions) in dynamics
        action_map = Dict(
            "total_propensity" => 0,
            "propensities" => [],
            "actions" => []
        )

        for possibility in actions
            propensity = possibility["propensity"]
            action_map["total_propensity"] += propensity
            push!(action_map["propensities"], propensity)
            push!(action_map["actions"], possibility["action"])
        end

        consolidate[state] = action_map
    end
    consolidate
end

function generate_dynamics()
    dynamics = Dict(
        EMPTY => [
            Dict(
                "propensity" => 1,
                "action" => no_action
            )
        ],
        SUBSTRATE => [
            Dict(
                "propensity" => 1,
                "action" => no_action
            ),
            Dict(
                "propensity" => 1,
                "action" => move_element
            )
        ],
        MEMBRANE => [
            Dict(
                "propensity" => 1,
                "action" => no_action
            ),
            Dict(
                "propensity" => 1,
                "action" => move_membrane
            )
        ],
        ENZYME => [
            Dict(
                "propensity" => 1,
                "action" => no_action
            ),
            Dict(
                "propensity" => 1,
                "action" => move_element
            )
        ],
        REPAIR => [
            Dict(
                "propensity" => 1,
                "action" => no_action
            ),
            Dict(
                "propensity" => 1,
                "action" => move_element
            )
        ],
    )

    consolidate_actions(dynamics)
end

function choose_action(possibilities)
    total = possibilities["total_propensity"]
    choice = rand(Uniform(0, total))
    index = 1
    action = possibilities["actions"][index]
    propensity = possibilities["propensities"][index]
    while choice > propensity
        choice -= propensity
        index += 1
        action = possibilities["actions"][index]
        propensity = possibilities["propensities"][index]
    end
    action
end

function generate_actions(world, dynamics)
    actions = []
    for index in state_indexes(world)
        state = world["states"][index]
        possibilities = dynamics[state]
        action = choose_action(possibilities)
        push!(actions, ([i for i in Tuple(index)], action))
    end
    shuffle(actions)
end

function apply_action(world, location_action)
    location, action = location_action
    action(world, location)
end

function apply_actions(world, actions)
    for action in actions
        world = apply_action(world, action)
    end
    world
end

# drawing

using Cairo

grid_size = [13, 13]

surface_size = [
    grid_size[1] * unit_size,
    grid_size[2] * unit_size
]

function initialize_cairo(surface)
    cairo_surface = CairoRGBSurface(
        surface[1],
        surface[2])
    cairo_surface, CairoContext(cairo_surface)
end

function set_background(cairo, surface)
    save(cairo)
    set_source_rgb(cairo, 0.8, 0.8, 0.8)
    rectangle(cairo, 0.0, 0.0, surface[1], surface[2])
    fill(cairo)
    restore(cairo)
end

function draw_nothing(cairo, unit, location, color)
    
end

function draw_substrate(cairo, unit, location, color)
    radius = unit / 3
    set_source_rgb(cairo, color[1], color[2], color[3])
    arc(cairo, location[1], location[2], radius, 0, 2*pi)
    stroke(cairo)
end

function draw_membrane(cairo, unit, location, color)
    radius = unit / 3
    set_source_rgb(cairo, color[1], color[2], color[3])
    rectangle(cairo, location[1] - radius, location[2] - radius, radius * 2, radius * 2)
    stroke(cairo)
end

function draw_enzyme(cairo, unit, location, color)
    radius = unit / 3    
    diagonal = sin(pi / 4) * radius

    move_to(cairo, location[1] - radius, location[2])
    line_to(cairo, location[1] + radius, location[2])
    close_path(cairo)

    move_to(cairo, location[1], location[2] - radius)
    line_to(cairo, location[1], location[2] + radius)
    close_path(cairo)

    move_to(cairo, location[1] - diagonal, location[2] - diagonal)
    line_to(cairo, location[1] + diagonal, location[2] + diagonal)
    close_path(cairo)

    move_to(cairo, location[1] - diagonal, location[2] + diagonal)
    line_to(cairo, location[1] + diagonal, location[2] - diagonal)
    close_path(cairo)

    set_source_rgb(cairo, color[1], color[2], color[3])
    stroke(cairo)
end

function draw_repair(cairo, unit, location, color)
    radius = unit / 3    
    a = [sin(2*pi/3), cos(2*pi/3)] * radius + location
    b = [sin(4*pi/3), cos(4*pi/3)] * radius + location
    c = [sin(6*pi/3), cos(6*pi/3)] * radius + location

    move_to(cairo, a[1], a[2])
    line_to(cairo, b[1], b[2])
    line_to(cairo, c[1], c[2])
    line_to(cairo, a[1], a[2])
    close_path(cairo)

    set_source_rgb(cairo, color[1], color[2], color[3])
    stroke(cairo)
end

function draw_link(cairo, unit, a, b, color)
    set_source_rgb(cairo, color[1], color[2], color[3])

    radius = unit / 3
    width = unit / 6
    vertical = a[1] == b[1]
    horizontal = a[2] == b[2]

    if vertical
        top = a[2]
        if b[2] < top
            top = b[2]
        end
        rectangle(cairo, a[1] - width, top + radius, width * 2, unit - (radius * 2))
        fill(cairo)
    elseif horizontal
        left = a[1]
        if b[1] < left
            left = b[1]
        end
        rectangle(cairo, left + radius, a[2] - width, unit - (radius * 2), width * 2)
        fill(cairo)
    else
        
    end
end

draw_symbols = Dict(
    EMPTY => draw_nothing,
    SUBSTRATE => draw_substrate,
    MEMBRANE => draw_membrane,
    ENZYME => draw_enzyme,
    REPAIR => draw_repair
)

function find_location(index, unit)
    location = (index - [1, 0]) * unit
    location += [unit / 2, -unit / 2]
end

function draw_world(cairo, surface_size, world, unit, symbols, colors)
    set_background(cairo, surface_size)
    set_line_width(cairo, 5.0)
    set_line_cap(cairo, Cairo.CAIRO_LINE_CAP_ROUND)

    for index in state_indexes(world)
        state = world["states"][index]
        if state > EMPTY
            location = find_location([i for i in Tuple(index)], unit)
            color = colors[state]
            draw = draw_symbols[state]
            draw(cairo, unit, location, color)
        end
    end

    for link in world["links"]
        a = find_location(link[:, 1], unit)
        b = find_location(link[:, 2], unit)
        draw_link(cairo, unit, a, b, colors[0])
    end
end

function test_world()
    world = generate_world(grid_size)
    set_state(world, [3, 4], MEMBRANE)
    set_state(world, [3, 5], MEMBRANE)
    set_state(world, [5, 7], MEMBRANE)
    set_state(world, [6, 7], MEMBRANE)
    set_state(world, [11, 8], ENZYME)
    set_state(world, [8, 12], REPAIR)
    add_link(world, [3, 4], [3, 5])
    add_link(world, [5, 7], [6, 7])

    println(world)
    world
end

cairo_surface, cairo = initialize_cairo(surface_size)

function run_simulation(bounds, counts)
    world = initialize_world(bounds, counts)
    dynamics = generate_dynamics()

    draw_world(cairo, surface_size, world, unit_size, symbols, colors)
    write_to_png(cairo_surface, "out/frames/tesselate-1.png")

    actions = generate_actions(world, dynamics)
    world = apply_actions(world, actions)

    draw_world(cairo, surface_size, world, unit_size, symbols, colors)
    write_to_png(cairo_surface, "out/frames/tesselate-2.png")
end

run_simulation(
    [13, 13],
    Dict(
        EMPTY => 21,
        MEMBRANE => 8,
        ENZYME => 5,
        REPAIR => 3
    )
)

end # module tesselate
