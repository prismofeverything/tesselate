module tesselate

function each(m)
    x, y = axes(m)
    eachindex(view(m, x, y))
end

# grid

unit_size = 40
grid_size = [30, 30]

function generate_world(grid)
    ones(Int8, grid[1], grid[2])
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

world = generate_world(grid_size)
world[3, 4] = 2
world[11, 8] = 3
world[8, 14] = 4

println(world)

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

function draw_link(cairo, unit, location, color)
    
end

draw_symbols = Dict(
    0 => draw_nothing,
    1 => draw_substrate,
    2 => draw_membrane,
    3 => draw_enzyme,
    4 => draw_repair
)

function draw_world(cairo, world, unit, symbols, colors)
    for index in each(world)
        state = world[index]
        if state > 0
            location = ([i for i in Tuple(index)] - [1, 0]) * unit
            location += [unit / 2, -unit / 2]
            color = colors[state]
            draw = draw_symbols[state]
            draw(cairo, unit, location, color)
        end
    end
end

cairo_surface, cairo = initialize_cairo(surface_size)
set_background(cairo, surface_size)
draw_world(cairo, world, unit_size, symbols, colors)

write_to_png(cairo_surface, "tesselate.png")

end # module tesselate
