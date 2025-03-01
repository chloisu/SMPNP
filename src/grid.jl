# 
# Copyright (C) 2025 - 2025 by the authors of the SMPNP code.
# 
# This file is part of SMPNP
# 
using ExtendableGrids

"""
This file implements the different grid types that we cover.
For each grid type, we implement the function to generate the particular 
grid. At the bottom of the file, we register the different grids available
in dictionary. The user can then select the given grid
"""

"""
equally spaced grid in one direction
"""

@option struct EquallySpaced1dParams
    pore_radius::Float64 = 10e-9 # pore radius in [meter]
    pore_length::Float64 = 50e-9 # pore length in [meter]
end

function get_equally_spaced_1d(parameters)
    # first, we normalize the geometry inputs
    r = parameters.grid.equally_spaced_1d.pore_radius / parameters.grid.equally_spaced_1d.pore_length
    X = collect(0:0.001:r)
    return simplexgrid(X)
end



"""
two-dimensional grid of a slit connected to two reservoirs with a symmetry 
in the middle of the slit. The mesh is refined towards the slit wall in the
"""
@option struct NanoSlitWithResevoirs2DParams
    pore_radius::Float64 = 10e-9 # pore radius in [meter]
    pore_length::Float64 = 50e-9 # pore length in [meter]
    reservoir_height::Float64 = 80e-9 # reservoir height [meter]
end

function get_nano_slit_with_reservoirs_2d(parameters)
    # first, we normalize the geometry inputs
    r = parameters.grid.nano_slit_with_reservoirs_2d.pore_radius / parameters.grid.nano_slit_with_reservoirs_2d.pore_length
    l = parameters.grid.nano_slit_with_reservoirs_2d.pore_length / parameters.grid.nano_slit_with_reservoirs_2d.pore_length
    h = parameters.grid.nano_slit_with_reservoirs_2d.reservoir_height / parameters.grid.nano_slit_with_reservoirs_2d.pore_length
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
    return grid
end


"""
This is the dictionary that will store all the different grids and return the proper one
"""
# Dictionary mapping keys to function calls
grids_dict = Dict(
    "equally_spaced_1d" => get_equally_spaced_1d,
    "nano_slit_with_reservoirs_2d" => get_nano_slit_with_reservoirs_2d,
)


"""
Generic function that calls the correct function from the dictionary
This function was generated with the help of ChatGPT (OpenAI).
"""
function get_grid(grid_type::String, args...)
    if haskey(grids_dict, grid_type)
        return grids_dict[grid_type](args...)
    else
        error("Function '$grid_type' not found in dictionary")
    end
end


"""
Here we define the struct of Grid Parameters which contains all the parameters of the different grid types.
"""
@option struct GridParameters
    equally_spaced_1d::EquallySpaced1dParams = EquallySpaced1dParams()
    nano_slit_with_reservoirs_2d::NanoSlitWithResevoirs2DParams = NanoSlitWithResevoirs2DParams()
end
