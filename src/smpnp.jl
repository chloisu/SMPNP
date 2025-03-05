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
include("config.jl")
include("nondimensionalization.jl")
include("grid.jl")
include("utils.jl")
include("boundary_conditions.jl")
include("initial_conditions.jl")
include("solver_config.jl")
include("time_config.jl")
include("parameters.jl")
include("chemical_potential.jl")
include("physics.jl")

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

"""
this is a small custom function to check if we are done with the time loop of the simulation.
That depends on whether or not we solve for a specific final time or for the steady state.
"""
function check_if_stay_in_time_loop(current_time, U_current_timestep, U_previous_timestep, parameters)
    # check if we have even started
    if current_time < parameters.time_parameters.initial_timestep
        return true
    end
    # else we check the conditions
    if parameters.time_parameters.solve_to_steady_state
        # check if we have reached the steady state
        if norm(U_current_timestep - U_previous_timestep, Inf) > parameters.time_parameters.steady_state_tol
            return true
        else
            # we are done
            return false
        end
    else
        # have we reached the final time
        if current_time < parameters.time_parameters.final_time
            return true
        else
            # we are done
            return false
        end
    end
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
    grid = get_grid(parameters.grid_type, parameters)#generate_grid_2d(parameters)#

    p = gridplot(grid; Plotter, size=(3000, 1000))

    sys = get_initial_timestep_system(grid, parameters)

    initial_timestep_initial_and_boundary_conditions_2d!(sys, parameters)
    #apply_dirichlet_for_initial_timestep!(sys, parameters)
    #apply_initial_conditions_for_intial_timestep!(sys, parameters)

    #control = VoronoiFVM.SolverControl()
    control = VoronoiFVM.NewtonControl()
    control.verbose = verbose
    control.reltol_linear = 1.0e-8
    control.abstol = 1e-6
    control.method_linear = method_linear
    control.maxiters = 5

    #setup_solver_control!(control, parameters)
    # solve the initial condition
    Δt = 1e-6
    time = 0
    U = solve(sys; control, Δt)
    initial_potential = U[1, :]
    writeVTK("out", grid, phi=U[1, :])

    if any(isnan.(U))
        error("Initial potential contains NaN values!")
    end
    p = GridVisualizer(;
        Plotter,
        layout=(3, 1),
        clear=true)

    # setup the time parameters correclty
    #setup_time_parameters!(parameters)

    sys2 = get_time_dependent_system(grid, parameters)
    #apply_dirichlet_for_time_dependent!(time_dependent_system, parameters)
    time_dependent_initial_and_boundary_conditions_2d!(sys2, initial_potential, parameters)
    #
    #U_current_timestep = apply_initial_condition_for_time_dependent!(time_dependent_system, parameters, initial_potential)
    #U_current_timestep = unknowns(time_dependent_system)
    #U_current_timestep[1, :] = initial_potential
    #U_current_timestep[2, :] .= 1.0
    #U_current_timestep[3, :] .= 1.0
    #U_previous_timestep = zeros(size(U_current_timestep))
    #U_previous_timestep .= U_current_timestep
    #time_dependent_initial_and_boundary_conditions_2d!(time_dependent_system, initial_potential, parameters)
    inival = unknowns(sys2)
    inival[1, :] = initial_potential
    inival[2, :] .= 1.0
    inival[3, :] .= 1.0
    # time loop
    t_plot = 0.0
    Δt = 1e-6#parameters.time_parameters.initial_timestep
    while time < 1000#check_if_stay_in_time_loop(time, U_current_timestep, U_previous_timestep, parameters)
        if time > t_plot
            #scalarplot!(p[1, 1], grid, U[1, :], xlimits=(0, 0.1), clear=false, show=true)
            #scalarplot!(p[2, 1], grid, U[2, :], xlimits=(0, 0.1), clear=false, show=true)
            #scalarplot!(p[3, 1], grid, U[3, :], xlimits=(0, 0.1), clear=false, show=true)
            t_plot = t_plot + 0.01
        end
        try
            print("Solving timestep at time: ")
            println(time)
            U = solve(sys2; inival, control, tstep=Δt)
            time = time + Δt
            inival .= U
            Δt *= 2
            Δt = min(Δt, 20)
        catch error
            Δt *= 0.5
            print("Repeating timestep at time: ")
            println(time)
            print("with timestep: ")
            println(Δt)
            # this means we didn't converge, so we will decrease the timestep
        end

    end
    nf = nodeflux(time_dependent_system, U)
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
