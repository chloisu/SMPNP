# 
# Copyright (C) 2025 - 2025 by the authors of the SMPNP code.
# 
# This file is part of SMPNP
# 

ENV["JULIA_NUM_THREADS"] = 1
using Printf
using VoronoiFVM
using ExtendableGrids
using ExtendableSparse
using GridVisualize
using LinearSolve
using ILUZero
using Triangulate
using SimplexGridFactory
using LinearAlgebra
using YAML
using Configurations
using Metal
using ArgParse

include("constants.jl")
include("parameters.jl")
include("chemical_potential.jl")
include("grid.jl")
include("utils.jl")


function initial_timestep_initial_and_boundary_conditions_2d!(sys, parameters)
    # we add the boundary conditions
    # we scale the surface potential to the correct value
    surface_potential_norm = parameters.surface_potential * E_CHARGE * BETA
    boundary_dirichlet!(sys, 1, 1, 0) # left reservoir 
    #boundary_dirichlet!(sys, 1, 4, 0) # right reservoir
    boundary_dirichlet!(sys, 1, 3, surface_potential_norm) # wall
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
    boundary_dirichlet!(sys, 1, 1, 0) # left reservoir 
    #boundary_dirichlet!(sys, 1, 4, 0) # right reservoir
    boundary_dirichlet!(sys, 1, 3, surface_potential_norm) # wall
    # anion boundary conditions
    boundary_dirichlet!(sys, 2, 1, 1.0) #left reservoir
    #boundary_dirichlet!(sys, 2, 4, 1.1) # right reservoir
    # cation bounda ry conditions
    boundary_dirichlet!(sys, 3, 1, 1.0)#left reservoir
    boundary_dirichlet!(sys, 3, 4, 1.1)# right reservoir

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
    L = parameters.pore_length
    prefactor = 1.0 / (4 * pi * L^2 * L_B * MOL_PER_LITER_TO_PER_CUBIC_METER)
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
    L = parameters.pore_length
    prefactor = 1.0 / (4 * pi * L^2 * L_B * MOL_PER_LITER_TO_PER_CUBIC_METER)
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

"""
function to parse the command line arguments of the code.
This function was generated with the help of ChatGPT (OpenAI).
"""
function parse_args()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--parameters", "-p"
        help = "Path to the parameter .yml file."
        arg_type = String
        required = true
        #"--verbose", "-v"
        #help = "Enable verbose mode"
        #action = :store_true  # Boolean flag
        #"--output", "-o"
        #help = "Output file path"
        #arg_type = String
        #default = "output.txt"
    end
    # Make a copy of ARGS to avoid modifying it directly
    args = copy(ARGS)

    # If the first argument is present but does not start with "--", assume it is --parameters
    if !isempty(args) && !startswith(args[1], "--")
        args = vcat(["--parameters"], args)  # Prepend `--parameters`
    end

    parsed_args = ArgParse.parse_args(args, s)
    return parsed_args
end

function main(;
    n=10, Plotter=nothing, verbose=false, unknown_storage=:sparse,
    method_linear=nothing, assembly=:edgewise
)
    # the very first thing we do is parse the command line arguments
    args = parse_args()
    # load the parameter file 
    parameters = load_config_file(args["parameters"])
    # next, we want to generate the output directory
    output_dir = create_output_directory(parameters)
    # Save parameters to YAML
    save_parameters_to_yaml(parameters, output_dir)
    # generate grid
    grid = get_grid(parameters.grid_type, parameters)

    #grid = generate_grid_2d(parameters)
    p = gridplot(grid; Plotter, size=(3000, 1000))

    sys = get_initial_timestep_system(grid, parameters)
    initial_timestep_initial_and_boundary_conditions_2d!(sys, parameters)

    control = VoronoiFVM.NewtonControl()
    control.verbose = verbose
    control.reltol_linear = 1.0e-8
    control.abstol = 1e-6
    control.method_linear = method_linear
    control.maxiters = 5

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
    time_dependent_initial_and_boundary_conditions_2d!(sys2, initial_potential, parameters)
    inival = unknowns(sys2)
    inival[1, :] = initial_potential
    inival[2, :] .= 1.0
    inival[3, :] .= 1.0
    # time loop
    t_plot = 0.0
    while time < 1000
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
            tstep = min(tstep, 20)
        catch error
            tstep *= 0.5
            print("Repeating timestep at time: ")
            println(time)
            print("with timestep: ")
            println(tstep)
            # this means we didn't converge, so we will decrease the timestep
        end

    end
    nf = nodeflux(sys2, U)
    surface_potential_norm = parameters.surface_potential * E_CHARGE * BETA
    #f = potential_pb_1d(grid, surface_potential_norm, parameters.pore_radius)
    #ca, cc = concentrations_pb_1d(grid, surface_potential_norm, parameters.pore_radius)
    scalarplot!(p[1, 1], grid, U[1, :], show=true, xlimits=(5, 5.1), ylimits=(4, 4.1))
    vectorplot!(p[1, 1], grid, nf[:, 1, :]; clear=false, vscale=1.5, xlimits=(5, 5.1), ylimits=(4, 4.1))

    #scalarplot!(p[1, 1], grid, f, clear=false, show=true, linestyle=:dash, color=(1, 0, 0))

    scalarplot!(p[2, 1], grid, U[2, :], clear=false, show=true, xlimits=(5, 5.1), ylimits=(4, 4.1))
    vectorplot!(p[2, 1], grid, nf[:, 2, :]; clear=false, vscale=1.5, xlimits=(5, 5.1), ylimits=(4, 4.1))

    #scalarplot!(p[2, 1], grid, ca, clear=false, show=true, linestyle=:dash, color=(1, 0, 0))
    scalarplot!(p[3, 1], grid, U[3, :], clear=false, show=true, xlimits=(5, 5.1), ylimits=(4, 4.1))
    vectorplot!(p[3, 1], grid, nf[:, 3, :]; clear=false, vscale=1.5, xlimits=(5, 5.1), ylimits=(4, 4.1))
    #scalarplot!(p[3, 1], grid, cc, clear=false, show=true, linestyle=:dash, color=(1, 0, 0))
    println("Norm of J- flux")
    println(sum(sum(abs.(nf[:, 2, :]))))

    println("Norm of J+ flux")
    println(sum(sum(abs.(nf[:, 3, :]))))

    writeVTK("out", grid, phi=U[1, :], cminus=U[2, :], cplus=U[3, :], nminus=nf[:, 2, :], nplus=nf[:, 3, :])
    return 0

end



main(verbose=true)#method_linear=KrylovJL_GMRES(precs=BlockPreconBuilder(precs=LinearSolvePreconBuilder(UMFPACKFactorization())))
