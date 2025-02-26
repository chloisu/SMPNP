
ENV["JULIA_NUM_THREADS"] = 1

using Printf
using VoronoiFVM
using ExtendableGrids
using ExtendableSparse: ILUZeroPreconBuilder
using GridVisualize
using LinearSolve
using ILUZero
using Triangulate
using SimplexGridFactory
using LinearAlgebra
using YAML
using Configurations
include("constants.jl")
include("analytical_solutions.jl")
include("parameters.jl")  # Include the parameters file
include("chemical_potential.jl")
#using .ParametersModule  # Use the module
#using Configurations
#@option struct Parameters
#    case::String = "test_case" # name of the case 
#    D_ref::Float64 = 5e-9 # reference diffusivity [m^2/s]
#    pore_radius::Float64 = 5e-9 # pore radius [meter]
#    pore_length::Float64 = 50e-9 # pore length [meter]
#    #cation_size::Float64 = 1e-10 # cation size [meter]
#    #anion_size::Float64 = 1e-10 # cation size [meter]
#    D_anion::Float64 = 5e-9 # anion diffusivity [m^2/s]
#    D_cation::Float64 = 5e-9 # cation diffusivity [m^2/s]
#    surface_potential::Float64 = 0.01 # surface potential [V]
#    #c_0::Float64 = 1 # initial concentration [mol/L]
#    #max_iter::Int = 20 # maximum nonlinear Newton iterations 
#    #dt::Float64 = 1e-6 # timestep
#    #dt_max::Float64 = 2e-1 # maximum timestep
#    #t_final::Float64 = 1 # final time
#end

function load_config_file(filename)
    # Define a struct for your configuration
    file = YAML.load_file(filename; dicttype=Dict{String,Any})
    parameters = from_dict(Parameters, file) # parse the 
    return parameters
end

function generate_grid_2d(parameters)
    # first, we normalize the geometry inputs
    r = parameters.pore_radius / parameters.pore_radius
    l = parameters.pore_length / parameters.pore_radius
    # we will generate a bounding box that goes from 0 to 3l with the middle l section to be the pore
    # for the radius, we will have the tube at the top with a total height of l
    builder = SimplexGridBuilder(; Generator=Triangulate)
    l = 5 * r
    cellregion!(builder, 1)
    p1 = point!(builder, 0, 0)
    p2 = point!(builder, l, 0)
    p3 = point!(builder, l, l - r)
    p4 = point!(builder, 2 * l, l - r)
    p5 = point!(builder, 2 * l, l)
    #p6 = point!(builder, 3 * l, 0)
    #p7 = point!(builder, 3 * l, l)
    p8 = point!(builder, 0, l)
    # left reservoir
    facetregion!(builder, 1)
    facet!(builder, p8, p1)
    # bottom reservoir
    facetregion!(builder, 2)
    facet!(builder, p1, p2)
    facet!(builder, p2, p3)
    #facet!(builder, p5, p6)
    # pore wall
    facetregion!(builder, 3)
    facet!(builder, p3, p4)
    facet!(builder, p4, p5)
    # right reservoir
    #facetregion!(builder, 4)
    #facet!(builder, p6, p7)
    # symmetry line
    facetregion!(builder, 4)
    facet!(builder, p5, p8)


    function unsuitable(x1, y1, x2, y2, x3, y3, area)
        bary = [(x1 + x2 + x3) / 3, (y1 + y2 + y3) / 3]
        needs_refinement = 0
        # towards the pore wall below
        min_x = l
        max_x = 2 * l
        rf_x = max(min_x, min(bary[1], max_x))
        refinement_center = [rf_x, l - r]
        dist = norm(bary - refinement_center)
        if area > 0.1 * dist
            needs_refinement = 1
        end
        # towards the right wall
        min_y = l - r
        max_y = l
        rf_y = max(min_y, min(bary[2], max_y))
        refinement_center = [2 * l, rf_y]
        dist = norm(bary - refinement_center)
        if area > 0.1 * dist
            needs_refinement = 1
        end

        return needs_refinement
    end
    options!(builder; unsuitable=unsuitable)
    grid = simplexgrid(builder)
    println(num_cells(grid))
    return grid
end

function generate_grid_1d(parameters)
    # first, we normalize the geometry inputs
    r = parameters.pore_radius / parameters.pore_radius
    X = collect(0:0.001:r)
    return simplexgrid(X)
end

function initial_timestep_initial_and_boundary_conditions_2d!(sys, parameters)
    # we add the boundary conditions
    # we scale the surface potential to the correct value
    surface_potential_norm = parameters.surface_potential * E_CHARGE * BETA
    boundary_dirichlet!(sys, 1, 1, 0)
    boundary_dirichlet!(sys, 1, 3, surface_potential_norm)
    # we set the inital potential distribution to zero
    inival = unknowns(sys)
    inival .= 0
end

function initial_timestep_initial_and_boundary_conditions_1d!(sys, parameters)
    # we add the boundary conditions
    # we scale the surface potential to the correct value
    surface_potential_norm = parameters.surface_potential * E_CHARGE * BETA
    boundary_dirichlet!(sys, 1, 1, surface_potential_norm)
    boundary_dirichlet!(sys, 1, 2, 0)
    # we set the inital potential distribution to zero
    inival = unknowns(sys)
    inival .= 0
end

function time_dependent_initial_and_boundary_conditions_2d!(sys, initial_potential, parameters)
    # potential boundary conditions
    surface_potential_norm = parameters.surface_potential * E_CHARGE * BETA
    boundary_dirichlet!(sys, 1, 1, 0)
    boundary_dirichlet!(sys, 1, 3, surface_potential_norm)
    # anion boundary conditions
    boundary_dirichlet!(sys, 2, 1, 1.0)
    #boundary_neumann!(sys, 2, 2, 0.0)
    # cation bounda ry conditions
    boundary_dirichlet!(sys, 3, 1, 1.0)
    #boundary_neumann!(sys, 3, 2, 0.0)
    U = unknowns(sys)
    U[1, :] = initial_potential
    U[2, :] .= 1.0
    U[3, :] .= 1.0
end
function time_dependent_initial_and_boundary_conditions_1d!(sys, initial_potential, parameters)
    # potential boundary conditions
    surface_potential_norm = parameters.surface_potential * E_CHARGE * BETA
    boundary_dirichlet!(sys, 1, 1, surface_potential_norm)
    boundary_dirichlet!(sys, 1, 2, 0)
    # anion boundary conditions
    boundary_dirichlet!(sys, 2, 2, 1.0)
    #boundary_neumann!(sys, 2, 2, 0.0)
    # cation bounda ry conditions
    boundary_dirichlet!(sys, 3, 2, 1.0)
    #boundary_neumann!(sys, 3, 2, 0.0)
    U = unknowns(sys)
    U[1, :] = initial_potential
    U[2, :] .= 1.0
    U[3, :] .= 1.0
end


"""
    get_initial_timestep_system(config)

returns a VoronoiFVM.System to solve the initial poisson 
equation for the potential distribution.
"""
function get_initial_timestep_system(grid, parameters)
    # we setup the physics for the poisson system only
    # we compute the prefactor for the Poisson equation
    r = parameters.pore_radius
    prefactor = 1.0 / (4 * pi * r^2 * L_B * MOL_PER_LITER_TO_PER_CUBIC_METER)
    # now we are ready to define the physics of this problem.
    physics = VoronoiFVM.Physics(;
        flux=function (f, u, edge, data)
            f[1] = prefactor * (u[1, 1] - u[1, 2])
            return nothing
        end
    )
    sys = VoronoiFVM.System(grid, physics; is_linear=false, unknown_storage=:sparse, assembly=:edgewise)
    # we add the potential 'species'
    enable_species!(sys, 1, [1])
    return sys
end


function get_time_dependent_system(grid, parameters=:nothing)
    # this system describes the phyiscs of the time-dependent PNP equations
    # the index of the different species are 
    # 1 = potential $\phi$
    # 2 = anion concentration $c_a$
    # 3 = cation concentration $c_c$

    # we compute the prefactor for the Poisson equation
    r = parameters.pore_radius
    prefactor = 1.0 / (4 * pi * r^2 * L_B * MOL_PER_LITER_TO_PER_CUBIC_METER)
    # we scale the diffusivities of the two concentration equations
    D_a_norm = parameters.D_anion / parameters.D_ref
    D_c_norm = parameters.D_cation / parameters.D_ref

    physics = VoronoiFVM.Physics(
        ; reaction=function (f, u, node, data)
            # source term of the poisson equation
            f[1] = -(Z_ANION * u[2] + Z_CATION * u[3])
            return nothing
        end,
        flux=function (f, u, edge, data)
            # compute the chemical potentials for the two cells
            μ_anion1 = μ(u[2, 1], u[3, 1], Z_ANION, u[1, 1], parameters)
            μ_anion2 = μ(u[2, 2], u[3, 2], Z_ANION, u[1, 2], parameters)

            μ_cation1 = μ(u[3, 1], u[2, 1], Z_CATION, u[1, 1], parameters)
            μ_cation2 = μ(u[3, 2], u[2, 2], Z_CATION, u[1, 2], parameters)
            # potential flux is simple poisson
            f[1] = prefactor * (u[1, 1] - u[1, 2])
            # anion flux is diffusion + migration
            c_anion_interface = 0.5 * (u[2, 1] + u[2, 2])
            f[2] = D_a_norm * c_anion_interface * (μ_anion1 - μ_anion2)#(((u[2, 1] - u[2, 2]) + c_anion_interface * Z_ANION * (u[1, 1] - u[1, 2])))
            # cation flux is diffusion + migration
            c_cation_interface = 0.5 * (u[3, 1] + u[3, 2])
            f[3] = D_c_norm * c_cation_interface * (μ_cation1 - μ_cation2)#((u[3, 1] - u[3, 2]) + c_cation_interface * Z_CATION * (u[1, 1] - u[1, 2]))
            return nothing
        end,
        storage=function (f, u, node, data)
            # time derivative for concentrations
            f[2] = u[2]
            f[3] = u[3]
            return nothing
        end
    )

    sys = VoronoiFVM.System(grid, physics; is_linear=false, unknown_storage=:sparse, assembly=:edgewise)
    enable_species!(sys, 1, [1]) # add potential
    enable_species!(sys, 2, [1]) # add anion
    enable_species!(sys, 3, [1]) # add cation
    return sys
end



function main(;
    n=10, Plotter=nothing, verbose=false, unknown_storage=:sparse,
    method_linear=nothing, assembly=:edgewise
)
    # load the parameter file 
    parameters = load_config_file("test.yml")

    # generate grid
    grid = generate_grid_1d(parameters)
    p = gridplot(grid; Plotter, size=(3000, 1000))


    sys = get_initial_timestep_system(grid, parameters)
    initial_timestep_initial_and_boundary_conditions_1d!(sys, parameters)

    control = VoronoiFVM.NewtonControl()
    control.verbose = verbose
    control.reltol_linear = 1.0e-8
    control.abstol = 1e-6
    control.method_linear = method_linear
    control.maxiters = 20

    tstep = 1e-6
    time = 0
    U = solve(sys; control, tstep)
    initial_potential = U[1, :]
    if any(isnan.(U))
        error("Initial potential contains NaN values!")
    end
    p = GridVisualizer(;
        Plotter,
        layout=(3, 1),
        clear=true)
    #scalarplot!(p[1, 1], grid, initial_potential, clear=true, show=true)

    # now solve the time-dependent once
    control2 = VoronoiFVM.NewtonControl()
    #control2.verbose = verbose
    control2.reltol_linear = 1.0e-8
    control2.method_linear = method_linear
    sys2 = get_time_dependent_system(grid, parameters)
    time_dependent_initial_and_boundary_conditions_1d!(sys2, initial_potential, parameters)
    inival = unknowns(sys2)
    inival[1, :] = initial_potential
    inival[2, :] .= 1.0
    inival[3, :] .= 1.0
    # time loop
    t_plot = 0.0
    while time < 20
        if time > t_plot
            #scalarplot!(p[1, 1], grid, U[1, :], xlimits=(0, 0.1), clear=false, show=true)
            #scalarplot!(p[2, 1], grid, U[2, :], xlimits=(0, 0.1), clear=false, show=true)
            #scalarplot!(p[3, 1], grid, U[3, :], xlimits=(0, 0.1), clear=false, show=true)
            t_plot = t_plot + 0.01
        end
        try
            print("Solving timestep at time: ")
            println(time)
            U = solve(sys2; inival, control, tstep)
            time = time + tstep
            inival .= U
            tstep *= 2
            tstep = min(tstep, 0.2)
        catch error
            tstep *= 0.5
            print("Repeating timestep at time: ")
            println(time)
            print("with timestep: ")
            println(tstep)
            # this means we didn't converge, so we will decrease the timestep
        end

    end
    surface_potential_norm = parameters.surface_potential * E_CHARGE * BETA
    f = potential_pb_1d(grid, surface_potential_norm, parameters.pore_radius)
    ca, cc = concentrations_pb_1d(grid, surface_potential_norm, parameters.pore_radius)
    scalarplot!(p[1, 1], grid, U[1, :], show=true)
    scalarplot!(p[1, 1], grid, f, clear=false, show=true, linestyle=:dash, color=(1, 0, 0))
    scalarplot!(p[2, 1], grid, U[2, :], clear=false, show=true)
    scalarplot!(p[2, 1], grid, ca, clear=false, show=true, linestyle=:dash, color=(1, 0, 0))
    scalarplot!(p[3, 1], grid, U[3, :], clear=false, show=true)
    scalarplot!(p[3, 1], grid, cc, clear=false, show=true, linestyle=:dash, color=(1, 0, 0))
    return reveal(p)

end

GC.gc()  # Force garbage collection
using GLMakie
p = main(Plotter=GLMakie, verbose=true)
GLMakie.save(joinpath(".", "out.jpg"), p)  #hide