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
using ArgParse
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
function check_if_stay_in_time_loop(current_time, U_current_timestep, U_previous_timestep, parameters, in_time_loop)
    # check if we have even started
    if !in_time_loop
        print("Exit from (1)")
        return true
    end
    # else we check the conditions
    if parameters.time_parameters.solve_to_steady_state
        # check if we have reached the steady state
        current_tol = maximum(abs.(U_current_timestep .- U_previous_timestep)) / (maximum(abs.(U_previous_timestep)) + eps())
        print("Current tolerance")
        println(current_tol)
        print("Tolerance")
        println(parameters.time_parameters.steady_state_tol)
        if current_time < parameters.time_parameters.minimum_steady_state_time || current_tol > parameters.time_parameters.steady_state_tol
            return true
        else
            # we are done
            print("Exit from (2)")
            return false
        end
    else
        # have we reached the final time
        if current_time < parameters.time_parameters.final_time
            return true
        else
            # we are done
            print("Exit from (3)")
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
    # Copy the provided input yml file as input.yml
    # and save all parameters including default to parameters.yml
    copy_parameters_yaml(args["parameters"], output_dir)
    save_parameters_to_yaml(parameters, output_dir)
    # generate grid
    grid = get_grid(parameters.grid_type, parameters)#generate_grid_2d(parameters)#
    # save the grid for postprocessing
    filename = joinpath(output_dir, "grid")
    writeVTK(filename * ".vtu", grid)
    @save filename * ".jld2" grid
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
    filename = joinpath(output_dir, @sprintf("%010.*f", TIME_WRITE_PRECISION, 0.0))
    writeVTK(filename * ".vtu", grid, phi=U[1, :])
    @save filename * ".jld2" U
    # check if there are any nan values in the initial solution
    if any(isnan.(U))
        error("Initial potential contains NaN values!")
    end
    # setup the time parameters correclty
    setup_time_parameters!(parameters) # non dimensionalize some of the input parameters
    time = 0.0
    Δt = parameters.time_parameters.initial_timestep
    t_plot = parameters.time_parameters.plot_time_interval
    # flush output
    flush(stdout)
    # setup the time-dependent system
    sys2 = get_time_dependent_system(grid, parameters)
    # apply boundary and initial conditions
    apply_dirichlet_for_time_dependent!(sys2, parameters)
    U = apply_initial_condition_for_time_dependent!(sys2, parameters, initial_potential)
    # setup helper solution vector
    U_previous_timestep = similar(U)
    U_previous_timestep .= U
    U_old = deepcopy(U)
    # solve the system for postprocessing
    filename = joinpath(output_dir, "system.jld2")
    system_state = VoronoiFVM.SystemState(sys2)
    @save filename system_state
    # time loop
    in_time_loop = false
    # time loop via simple call to solve
    while check_if_stay_in_time_loop(time, U, U_old, parameters, in_time_loop)
        in_time_loop = true
        try
            print("Solving timestep at time: ")
            println(time)
            U = solve(sys2; inival=U_previous_timestep, control=control, tstep=Δt)
            time = time + Δt
            U_old .= U_previous_timestep
            U_previous_timestep .= U
            Δt *= 2
        catch error
            # this means we didn't converge, so we will decrease the timestep
            Δt *= 0.5
            print("Repeating timestep at time: ")
            println(time)
            print("with timestep: ")
            println(Δt)
            in_time_loop = false
            continue
        end
        if time ≈ t_plot
            filename = joinpath(output_dir, @sprintf("%010.*f", TIME_WRITE_PRECISION, time))
            writeVTK(filename * ".vtu", grid, phi=U[1, :], cminus=U[2, :], cplus=U[3, :])
            @save filename * ".jld2" U
            # next plot in
            t_plot = t_plot + parameters.time_parameters.plot_time_interval
        end
        # decide on new timestep
        Δt = min(min(Δt, parameters.time_parameters.max_timestep), t_plot - time)
        # flush output 
        flush(stdout)
    end
    # and we are done
    generate_pvd(output_dir, "out.pvd")
    return grid, U
end

# Only call main() if the script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    smpnp()
end
#method_linear=KrylovJL_GMRES(precs=BlockPreconBuilder(precs=LinearSolvePreconBuilder(UMFPACKFactorization())))
