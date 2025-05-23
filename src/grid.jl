# 
# Copyright (C) 2025 - 2025 by the authors of the SMPNP code.
# 
# This file is part of SMPNP
# 

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
end

function validate_equally_spaced_1d_parameters(parameters)
    params = parameters.grid.equally_spaced_1d
    params.pore_radius > 0.0 || throw(AssertionError("The pore radius must be larger than 0."))

end

function get_equally_spaced_1d(parameters)
    # first, we normalize the geometry inputs
    r = parameters.grid.equally_spaced_1d.pore_radius / parameters.non_dim.L_REF
    X = collect(0:0.001:r)
    return simplexgrid(X)
end



"""
two-dimensional grid of a slit connected to two reservoirs with a symmetry 
in the middle of the slit. The mesh is refined towards the slit wall in the
"""
@option mutable struct NanoSlitWithResevoirs2DParams
    pore_radius::Float64 = 10e-9 # pore radius in [meter]
    pore_length::Float64 = 50e-9 # pore length in [meter]
    reservoir_height::Float64 = 80e-9 # reservoir height [meter]
end

function validate_nano_slit_with_reservoirs_2d(parameters)
    params = parameters.grid.nano_slit_with_reservoirs_2d
    params.pore_radius > 0.0 || throw(AssertionError("The pore radius must be larger than 0."))
    params.pore_length > 0.0 || throw(AssertionError("The pore length must be larger than 0."))
    params.reservoir_height > 0.0 || throw(AssertionError("The reservoir height must be larger than 0."))
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
Aligned NanoSlit
"""
@option mutable struct AlignedNanoSlitWithResevoirs2DParams
    pore_radius::Float64 = 10e-9 # pore radius in [meter]
    pore_length::Float64 = 50e-9 # pore length in [meter]
    reservoir_height::Float64 = 80e-9 # reservoir height [meter]
end

function validate_aligned_nano_slit_with_reservoirs_2d(parameters)
    params = parameters.grid.aligned_nano_slit_with_reservoirs_2d
    params.pore_radius > 0.0 || throw(AssertionError("The pore radius must be larger than 0."))
    params.pore_length > 0.0 || throw(AssertionError("The pore length must be larger than 0."))
    params.reservoir_height > 0.0 || throw(AssertionError("The reservoir height must be larger than 0."))
end

function get_aligned_nano_slit_with_reservoirs_2d(parameters)
    # first, we normalize the geometry inputs
    params = parameters.grid.aligned_nano_slit_with_reservoirs_2d
    r = params.pore_radius / parameters.non_dim.L_REF
    l = params.pore_length / parameters.non_dim.L_REF
    h = params.reservoir_height / parameters.non_dim.L_REF
    # we will generate a bounding box that goes from 0 to 3l with the middle l section to be the pore
    # for the radius, we will have the tube at the top with a total height of l
    hmin = 0.005
    hmax = 0.1
    Xleft = geomspace(0.0, l, hmax, hmin)
    Xmid = collect(l:hmin:2*l)
    Xright = geomspace(2 * l, 3 * l, hmin, hmax)
    hbottom = 0.05
    hwall = 0.005
    hcenterline = 0.01
    Ybot = geomspace(0.0, h - r, hmax, hwall)
    Ytop = geomspace(h - r, h, hwall, hcenterline)
    X = glue(Xleft, glue(Xmid, Xright))
    Y = glue(Ybot, Ytop)
    parent = simplexgrid(X, Y)
    rect!(parent, [l, 0], [2 * l, h - r]; region=2, bregion=1)
    grid = subgrid(parent, [1])
    bfacemask!(grid, [l, h - r], [2 * l, h - r], 5)
    return grid
end

"""
Aligned Channel with two Unstructured Reservoirs
"""
@option mutable struct AlignedNanoSlitWithUnstructuredResevoirs2DParams
    pore_radius::Float64 = 10e-9 # pore radius in [meter]
    pore_length::Float64 = 50e-9 # pore length in [meter]
    reservoir_height::Float64 = 80e-9 # reservoir height [meter]
    hwall_normal::Float64 = 0.0002 # cell size normal to the wall [non-dimensional units]
    hwall_tangential::Float64 = 0.002 # cell size tangential to the wall [non-dimensional units]
    hcenterline::Float64 = 0.01 # cell size at the centerline in fraction of non-dimensional radius 
    #(i.e. 0.1 means that the cell will have a size of one-tenth of the radius.)
    maxvolume::Float64 = 0.0001 # maximum cell volume in non-dimensional units for the simplex grid builder
end

function validate_aligned_nano_slit_with_unstructured_reservoirs_2d(parameters)
    params = parameters.grid.aligned_nano_slit_with_unstructured_reservoirs_2d
    params.pore_radius > 0.0 || throw(AssertionError("The pore radius must be larger than 0."))
    params.pore_length > 0.0 || throw(AssertionError("The pore length must be larger than 0."))
    params.reservoir_height > 0.0 || throw(AssertionError("The reservoir height must be larger than 0."))
    params.hwall_normal > 0.0 || throw(AssertionError("The wall normal cell height must be larger than 0."))
    params.hwall_tangential > 0.0 || throw(AssertionError("The wall tangential cell width must be larger than 0."))
    params.hcenterline > 0.0 || throw(AssertionError("The centerline cell height must be larger than 0."))
    params.maxvolume > 0.0 || throw(AssertionError("The maximum cell volume must be larger than 0."))
end

function get_aligned_nano_slit_with_unstructured_reservoirs_2d(parameters)
    # first, we normalize the geometry inputs
    params = parameters.grid.aligned_nano_slit_with_unstructured_reservoirs_2d
    r = params.pore_radius / parameters.non_dim.L_REF
    l = params.pore_length / parameters.non_dim.L_REF
    h = params.reservoir_height / parameters.non_dim.L_REF
    # we first generate the grid for the nano channel only
    hwall_normal = params.hwall_normal
    hwall_tangential = params.hwall_tangential
    hcenterline = params.hcenterline * r
    Xchannel = collect(l:hwall_tangential:2*l)
    Ychannel = geomspace(h - r, h, hwall_normal, hcenterline)
    channel_grid = simplexgrid(Xchannel, Ychannel)
    # next, we generate the grid builder for the left reservoir
    left_reservoir_builder = SimplexGridBuilder(; Generator=Triangulate)
    cellregion!(left_reservoir_builder, 1) # this is our inside domain
    p1 = point!(left_reservoir_builder, 0, 0)
    p2 = point!(left_reservoir_builder, l, 0)
    p3 = point!(left_reservoir_builder, l, h - r)
    p4 = point!(left_reservoir_builder, l, h)
    p5 = point!(left_reservoir_builder, 0, h)
    # we also add the boundaries
    facetregion!(left_reservoir_builder, 5)
    facet!(left_reservoir_builder, p5, p1)
    facetregion!(left_reservoir_builder, 6)
    facet!(left_reservoir_builder, p1, p2)
    facet!(left_reservoir_builder, p2, p3)
    facetregion!(left_reservoir_builder, 4)
    facet!(left_reservoir_builder, p3, p4)
    facetregion!(left_reservoir_builder, 3)
    facet!(left_reservoir_builder, p4, p5)
    bregions!(left_reservoir_builder, channel_grid, 4)
    left_reservoir = simplexgrid(left_reservoir_builder, maxvolume=params.maxvolume)
    # we glue the left reservoir to the channel
    left_plus_channel = glue(left_reservoir, channel_grid)
    # finally, we do the same for the right reservior
    # next, we generate the grid builder for the left reservoir
    right_reservoir_builder = SimplexGridBuilder(; Generator=Triangulate)
    cellregion!(right_reservoir_builder, 1) # this is our inside domain
    p1r = point!(right_reservoir_builder, 3 * l, 0)
    p2r = point!(right_reservoir_builder, 2 * l, 0)
    p3r = point!(right_reservoir_builder, 2 * l, h - r)
    p4r = point!(right_reservoir_builder, 2 * l, h)
    p5r = point!(right_reservoir_builder, 3 * l, h)
    # we also add the boundaries
    facetregion!(right_reservoir_builder, 9)
    facet!(right_reservoir_builder, p5r, p1r)
    facetregion!(right_reservoir_builder, 6)
    facet!(right_reservoir_builder, p1r, p2r)
    facet!(right_reservoir_builder, p2r, p3r)
    facetregion!(right_reservoir_builder, 2)
    facet!(right_reservoir_builder, p3r, p4r)
    facetregion!(right_reservoir_builder, 3)
    facet!(right_reservoir_builder, p4r, p5r)
    bregions!(right_reservoir_builder, channel_grid, 2)
    right_reservoir = simplexgrid(right_reservoir_builder, maxvolume=params.maxvolume)
    # we glue the right reservoir to the channel + left
    grid = glue(left_plus_channel, right_reservoir)
    return grid
end

"""
Partially Aligned Channel with Rest Unstructured and Unstructured Reservoirs
"""
@option mutable struct PartiallyAlignedNanoSlitWithRestUnstructuredAndUnstructuredResevoirs2DParams
    pore_radius::Float64 = 10e-9 # pore radius in [meter]
    pore_length::Float64 = 50e-9 # pore length in [meter]
    reservoir_height::Float64 = 80e-9 # reservoir height [meter]
    width_of_structured_region::Float64 = 20e-9 # width of the structured region. [meter] 
    # The region will be centered around the middle of the x-axis (i.e. middle of nanopore).
    hwall_normal::Float64 = 0.0002 # cell size normal to the wall [non-dimensional units]
    hwall_tangential::Float64 = 0.002 # cell size tangential to the wall [non-dimensional units]
    hcenterline::Float64 = 0.01 # cell size at the centerline in fraction of non-dimensional radius 
    #(i.e. 0.1 means that the cell will have a size of one-tenth of the radius.)
    maxvolume::Float64 = 0.0001 # maximum cell volume in non-dimensional units for the simplex grid builder
end

function validate_partially_aligned_nano_slit_with_rest_unstructured_and_unstructured_reservoirs_2d(parameters)
    params = parameters.grid.partially_aligned_nano_slit_with_rest_unstructured_and_unstructured_reservoirs_2d
    params.pore_radius > 0.0 || throw(AssertionError("The pore radius must be larger than 0."))
    params.pore_length > 0.0 || throw(AssertionError("The pore length must be larger than 0."))
    params.reservoir_height > 0.0 || throw(AssertionError("The reservoir height must be larger than 0."))
    params.width_of_structured_region > 0.0 || throw(AssertionError("The structured region width must be larger than 0."))
    params.hwall_normal > 0.0 || throw(AssertionError("The wall normal cell height must be larger than 0."))
    params.hwall_tangential > 0.0 || throw(AssertionError("The wall tangential cell width must be larger than 0."))
    params.hcenterline > 0.0 || throw(AssertionError("The centerline cell height must be larger than 0."))
    params.maxvolume > 0.0 || throw(AssertionError("The maximum cell volume must be larger than 0."))
end

function get_partially_aligned_nano_slit_with_rest_unstructured_and_unstructured_reservoirs_2d(parameters)
    # first, we normalize the geometry inputs
    params = parameters.grid.partially_aligned_nano_slit_with_rest_unstructured_and_unstructured_reservoirs_2d
    r = params.pore_radius / parameters.non_dim.L_REF
    l = params.pore_length / parameters.non_dim.L_REF
    h = params.reservoir_height / parameters.non_dim.L_REF
    w = params.width_of_structured_region / parameters.non_dim.L_REF
    # we first generate the grid for the nano channel only
    hwall_normal = params.hwall_normal
    hwall_tangential = params.hwall_tangential
    hcenterline = params.hcenterline * r
    Xchannel = collect(1.5*l-0.5*w:hwall_tangential:1.5*l+0.5*w)
    Ychannel = geomspace(h - r, h, hwall_normal, hcenterline)
    channel_grid = simplexgrid(Xchannel, Ychannel)
    # helper function for the refinement toward the wall 
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
    # next, we generate the grid builder for the left reservoir
    left_reservoir_builder = SimplexGridBuilder(; Generator=Triangulate)
    cellregion!(left_reservoir_builder, 1) # this is our inside domain
    p1 = point!(left_reservoir_builder, 0, 0)
    p2 = point!(left_reservoir_builder, l, 0)
    p3 = point!(left_reservoir_builder, l, h - r)
    p4 = point!(left_reservoir_builder, 1.5 * l - 0.5 * w, h - r)
    p5 = point!(left_reservoir_builder, 1.5 * l - 0.5 * w, h)
    p6 = point!(left_reservoir_builder, 0, h)
    # we also add the boundaries
    facetregion!(left_reservoir_builder, 5)
    facet!(left_reservoir_builder, p6, p1)
    facetregion!(left_reservoir_builder, 6)
    facet!(left_reservoir_builder, p1, p2)
    facet!(left_reservoir_builder, p2, p3)
    facetregion!(left_reservoir_builder, 1)
    facet!(left_reservoir_builder, p3, p4)
    facetregion!(left_reservoir_builder, 4)
    facet!(left_reservoir_builder, p4, p5)
    facetregion!(left_reservoir_builder, 3)
    facet!(left_reservoir_builder, p5, p6)
    bregions!(left_reservoir_builder, channel_grid, 4)
    left_reservoir = simplexgrid(left_reservoir_builder, unsuitable=unsuitable, maxvolume=params.maxvolume)
    # we glue the left reservoir to the channel
    left_plus_channel = glue(left_reservoir, channel_grid)
    # finally, we do the same for the right reservior
    # next, we generate the grid builder for the left reservoir
    right_reservoir_builder = SimplexGridBuilder(; Generator=Triangulate)
    cellregion!(right_reservoir_builder, 1) # this is our inside domain
    p1r = point!(right_reservoir_builder, 3 * l, 0)
    p2r = point!(right_reservoir_builder, 2 * l, 0)
    p3r = point!(right_reservoir_builder, 2 * l, h - r)
    p4r = point!(right_reservoir_builder, 1.5 * l + 0.5 * w, h - r)
    p5r = point!(right_reservoir_builder, 1.5 * l + 0.5 * w, h)
    p6r = point!(right_reservoir_builder, 3 * l, h)
    # we also add the boundaries
    facetregion!(right_reservoir_builder, 9)
    facet!(right_reservoir_builder, p6r, p1r)
    facetregion!(right_reservoir_builder, 6)
    facet!(right_reservoir_builder, p1r, p2r)
    facet!(right_reservoir_builder, p2r, p3r)
    facetregion!(right_reservoir_builder, 1)
    facet!(right_reservoir_builder, p3r, p4r)
    facetregion!(right_reservoir_builder, 2)
    facet!(right_reservoir_builder, p4r, p5r)
    facetregion!(right_reservoir_builder, 3)
    facet!(right_reservoir_builder, p5r, p6r)
    bregions!(right_reservoir_builder, channel_grid, 2)
    right_reservoir = simplexgrid(right_reservoir_builder, unsuitable=unsuitable, maxvolume=params.maxvolume)
    # we glue the right reservoir to the channel + left
    grid = glue(left_plus_channel, right_reservoir)
    return grid
end

"""
This is the dictionary that will store all the different grids and return the proper one
"""
# Dictionary mapping keys to function calls
grids_dict = Dict(
    "equally_spaced_1d" => get_equally_spaced_1d,
    "nano_slit_with_reservoirs_2d" => get_nano_slit_with_reservoirs_2d,
    "aligned_nano_slit_with_reservoirs_2d" => get_aligned_nano_slit_with_reservoirs_2d,
    "aligned_nano_slit_with_unstructured_reservoirs_2d" => get_aligned_nano_slit_with_unstructured_reservoirs_2d,
    "partially_aligned_nano_slit_with_rest_unstructured_and_unstructured_reservoirs_2d" => get_partially_aligned_nano_slit_with_rest_unstructured_and_unstructured_reservoirs_2d)

# Dictionary mapping keys to function calls
validation_dict = Dict(
    "equally_spaced_1d" => validate_equally_spaced_1d_parameters,
    "nano_slit_with_reservoirs_2d" => validate_nano_slit_with_reservoirs_2d,
    "aligned_nano_slit_with_reservoirs_2d" => validate_aligned_nano_slit_with_reservoirs_2d,
    "partially_aligned_nano_slit_with_rest_unstructured_and_unstructured_reservoirs_2d" => validate_partially_aligned_nano_slit_with_rest_unstructured_and_unstructured_reservoirs_2d)

"""
asdf
This function was generated with the help of ChatGPT (OpenAI).
"""
function validate_grid_type(mode::String)
    # Check if the mode is a key in the dictionary
    mode in keys(grids_dict) || throw(NotImplementedError("Not Implemented Grid type: $mode. The grid_type specified in the parameter file is not yet implemented. Currently available values are: " * string(collect(keys(grids_dict)))))
end

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
Generic function that calls the correct function from the dictionary
This function was generated with the help of ChatGPT (OpenAI).
"""
function validate_grid_parameters(grid_type::String, args...)
    if haskey(grids_dict, grid_type)
        validation_dict[grid_type](args...)
    else
        error("Validation Function for grid type '$grid_type' not found in validation dictionary.")
    end
end


"""
Here we define the struct of Grid Parameters which contains all the parameters of the different grid types.
"""
@option struct GridParameters
    equally_spaced_1d::EquallySpaced1dParams = EquallySpaced1dParams()
    nano_slit_with_reservoirs_2d::NanoSlitWithResevoirs2DParams = NanoSlitWithResevoirs2DParams()
    aligned_nano_slit_with_reservoirs_2d::AlignedNanoSlitWithResevoirs2DParams = AlignedNanoSlitWithResevoirs2DParams()
    aligned_nano_slit_with_unstructured_reservoirs_2d::AlignedNanoSlitWithUnstructuredResevoirs2DParams = AlignedNanoSlitWithUnstructuredResevoirs2DParams()
    partially_aligned_nano_slit_with_rest_unstructured_and_unstructured_reservoirs_2d::PartiallyAlignedNanoSlitWithRestUnstructuredAndUnstructuredResevoirs2DParams = PartiallyAlignedNanoSlitWithRestUnstructuredAndUnstructuredResevoirs2DParams()
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