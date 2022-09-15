module tesselate

function each(m)
    x, y = axes(m)
    eachindex(view(m, x, y))
end

# grid

unit_size = 40
grid_size = [30, 30]

function generate_world(grid)
    Dict(
        "states" => ones(Int8, grid[1], grid[2]),
        "links" => Vector{Matrix{Int8}}()
    )
end

function set_state(world, location, state)
    world["states"][location[1], location[2]] = state
    world
end

function add_link(world, a, b)
    link = reduce(hcat, [a, b])
    push!(world["links"], link)
    world
end

symbols = Dict(
    0 => " ",
    1 => "◯",
    2 => "⧇",
    3 => "*",
    4 => "△"
)

colors = Dict(
    0 => [0, 0, 0],
    1 => [0.6, 0.6, 0.6],
    2 => [0.8, 0.7, 0.1],
    3 => [0.7, 0.1, 0.6],
    4 => [0.1, 0.7, 0.4]
)

function initialize_world()
    world = generate_world(grid_size)
    set_state(world, [3, 4], 2)
    set_state(world, [3, 5], 2)
    set_state(world, [5, 7], 2)
    set_state(world, [6, 7], 2)
    set_state(world, [11, 8], 3)
    set_state(world, [8, 14], 4)
    add_link(world, [3, 4], [3, 5])
    add_link(world, [5, 7], [6, 7])

    println(world)
    world
end

# dynamics

# O -> A
# [] -> B
# * -> F
# △ -> Phi
# all we are missing is a △ to "repair" our * from [], and also to have [] act on * to produce △

function mod_location(bounds, location)
    [
        mod1(location[1], bounds[1]),
        mod1(location[2], bounds[2]),
    ]
end

function adjacent_locations(world, location)
    bounds = 
    [[location[1], location[2]]]
end

function no_action(world, location)
    world
end

function move_substrate(world, location)
    
end

function generate_dynamics()
    Dict(
        0 => Dict(
            
        ),
        1 => [
            Dict(
                "propensity" => 1,
                "action" => no_action
            ),
            Dict(
                "propensity" => 1,
                "action" => move_substrate
            )
        ],
        2 => Dict(

        ),
        3 => Dict(

        ),
        4 => Dict(

        )
    )
end

function generate_actions(world, dynamics)
    for index in each(world["states"])
        state = world["states"][index]
        possibilities = dynamics[state]
    end    
end

# drawing

using Cairo

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
    0 => draw_nothing,
    1 => draw_substrate,
    2 => draw_membrane,
    3 => draw_enzyme,
    4 => draw_repair
)

function find_location(index, unit)
    location = (index - [1, 0]) * unit
    location += [unit / 2, -unit / 2]
end

function draw_world(cairo, world, unit, symbols, colors)
    set_line_width(cairo, 5.0)
    set_line_cap(cairo, Cairo.CAIRO_LINE_CAP_ROUND)

    for index in each(world["states"])
        state = world["states"][index]
        if state > 0
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

cairo_surface, cairo = initialize_cairo(surface_size)
set_background(cairo, surface_size)

world = initialize_world()
draw_world(cairo, world, unit_size, symbols, colors)

write_to_png(cairo_surface, "tesselate.png")

end # module tesselate
