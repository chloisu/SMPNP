# 
# Copyright (C) 2025 - 2025 by the authors of the SMPNP code.
# 
# This file is part of SMPNP
# 


"""
simple 1D grid
"""
function generate_grid_1d(parameters)
    # first, we normalize the geometry inputs
    r = parameters.pore_radius / parameters.pore_length
    X = collect(0:0.001:r)
    return simplexgrid(X)
end

"""
this function generates the 2D grid of the simulation with the pore geometry and the refinement towards the pore wall.
"""
function generate_grid_2d(parameters)
    # first, we normalize the geometry inputs
    r = parameters.pore_radius / parameters.pore_length
    l = parameters.pore_length / parameters.pore_length
    h = parameters.reservoir_height / parameters.pore_length
    # we will generate a bounding box that goes from 0 to 3l with the middle l section to be the pore
    # for the radius, we will have the tube at the top with a total height of l
    builder = SimplexGridBuilder(; Generator=Triangulate)
    cellregion!(builder, 1) # this is our inside domain
    p1 = point!(builder, 0, 0)
    p2 = point!(builder, l, 0)
    p3 = point!(builder, l, h - r)
    p4 = point!(builder, 2 * l, h - r)
    p5 = point!(builder, 2 * l, 0)
    p6 = point!(builder, 3 * l, 0)
    p7 = point!(builder, 3 * l, h)
    p8 = point!(builder, 0, h)
    # left reservoir
    facetregion!(builder, 1)
    facet!(builder, p8, p1)
    # bottom reservoir
    facetregion!(builder, 2)
    facet!(builder, p1, p2)
    facet!(builder, p2, p3)
    facet!(builder, p4, p5)
    facet!(builder, p5, p6)
    # pore wall
    facetregion!(builder, 3)
    facet!(builder, p3, p4)
    # right reservoir
    facetregion!(builder, 4)
    facet!(builder, p6, p7)
    # symmetry line
    facetregion!(builder, 5)
    facet!(builder, p7, p8)

    """
    this function returns a bool telling us which cells need refinement.
    """
    function unsuitable(x1, y1, x2, y2, x3, y3, area)
        bary = [(x1 + x2 + x3) / 3, (y1 + y2 + y3) / 3]
        needs_refinement = 0
        # towards the pore wall below
        min_x = l
        max_x = 2 * l
        rf_x = max(min_x, min(bary[1], max_x))
        refinement_center = [rf_x, h - r]
        dist = norm(bary - refinement_center)
        if area > 0.0005 * dist
            needs_refinement = 1
        end
        return needs_refinement
    end
    options!(builder; unsuitable=unsuitable)
    grid = simplexgrid(builder)
    println(num_cells(grid))
    return grid
end