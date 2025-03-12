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
using GLMakie
using JLD2

include("constants.jl")
include("config.jl")
include("nondimensionalization.jl")
include("grid.jl")
include("utils.jl")
include("boundary_conditions.jl")
include("initial_conditions.jl")
include("solver_config.jl")
include("time_config.jl")
include("pnp_config.jl")
include("parameters.jl")
include("chemical_potential.jl")
include("physics.jl")

function initial_timestep_initial_and_boundary_conditions_2d!(sys, parameters)
    # we add the boundary conditions
    # we scale the surface potential to the correct value
    surface_potential_norm = 11.8#parameters.surface_potential * E_CHARGE * BETA
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
    surface_potential_norm = 11.8#parameters.surface_potential * E_CHARGE * BETA
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

function smpnp()
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
    # save the grid for postprocessing
    filename = joinpath(output_dir, "grid.jld2")
    @save filename grid
    # TODO: make the grid output nicer
    vis = GridVisualizer(Plotter=GLMakie)
    gridplot!(vis, grid; size=(3000, 1000))
    GLMakie.save(joinpath(output_dir, "grid.png"), reveal(vis))
    # setup the system for the initial timestep
    sys = get_initial_timestep_system(grid, parameters)
    # specfiy boundary + initial conditions for the intial timestep
    apply_dirichlet_for_initial_timestep!(sys, parameters)
    apply_initial_conditions_for_intial_timestep!(sys, parameters)
    # setup the Newton solver
    control = VoronoiFVM.SolverControl()
    setup_solver_control!(control, parameters)
    # solve the initial condition
    U = solve(sys; control)
    # write the initial potential distribution as VTK
    @assert PHI_EQ == 1 # make sure that the potential species is the first species
    initial_potential = U[1, :]
    # solve initial timestep
    nf = nodeflux(sys, U)
    filename = joinpath(output_dir, @sprintf("%.5f.jld2", 0.0))
    @save filename U
    #writeVTK(joinpath(output_dir, filename), grid, phi=U[1, :], gradphi=nf[:, 1, :])
    # check if there are any nan values in the initial solution
    if any(isnan.(U))
        error("Initial potential contains NaN values!")
    end
    # setup the time parameters correclty
    setup_time_parameters!(parameters) # non dimensionalize some of the input parameters
    time = 0.0
    Δt = parameters.time_parameters.initial_timestep
    t_plot = parameters.time_parameters.plot_time_interval

    # setup the time-dependent system
    sys2 = get_time_dependent_system(grid, parameters)
    # apply boundary and initial conditions
    apply_dirichlet_for_time_dependent!(sys2, parameters)
    U = apply_initial_condition_for_time_dependent!(sys2, parameters, initial_potential)
    # setup helper solution vector
    U_previous_timestep = similar(U)
    U_previous_timestep .= U
    # solve the system for postprocessing
    filename = joinpath(output_dir, "system.jld2")
    system_state = VoronoiFVM.SystemState(sys2)
    @save filename system_state
    # time loop
    while time < 100#check_if_stay_in_time_loop(time, U, U_previous_timestep, parameters)
        try
            print("Solving timestep at time: ")
            println(time)
            U = solve(sys2; inival=U_previous_timestep, control=control, tstep=Δt)
            time = time + Δt
            U_previous_timestep .= U
            Δt *= 2
        catch error
            # this means we didn't converge, so we will decrease the timestep
            Δt *= 0.5
            print("Repeating timestep at time: ")
            println(time)
            print("with timestep: ")
            println(Δt)
        end
        if time ≈ t_plot
            nf = nodeflux(sys2, U)
            filename = joinpath(output_dir, @sprintf("%.5f.jld2", time))
            #writeVTK(joinpath(output_dir, filename), grid, phi=U[1, :], cminus=U[2, :], cplus=U[3, :], nminus=nf[:, 2, :], nplus=nf[:, 3, :])
            @save filename U
            # next plot in
            t_plot = t_plot + parameters.time_parameters.plot_time_interval
        end
        tf = TestFunctionFactory(sys2)
        T = testfunction(tf, [1], [3])
        I = integrate(sys2, T, U)
        println("Current Flux")
        println(I)
        # decide on new timestep
        Δt = min(min(Δt, parameters.time_parameters.max_timestep), t_plot - time)

    end
    # and we are done
    return grid, U

end

# Only call main() if the script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    smpnp()
end
#method_linear=KrylovJL_GMRES(precs=BlockPreconBuilder(precs=LinearSolvePreconBuilder(UMFPACKFactorization())))
