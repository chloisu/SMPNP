"""
Here we define the struct of Grid Parameters which contains all the parameters of the different grid types.
This is a small wrapper to make the solver arguments of VoronoiFVM.SolverControl() available through the 
parameter file.
"""
@option mutable struct TimeParameters
    solve_to_steady_state::Bool = true
    steady_state_tol::Float64 = 1e-14
    final_time::Float64 = 1000 # measured in [s]
    max_timestep::Float64 = (50e-9)^2 / 5e-9 #  measure in [s]
    initial_timestep::Float64 = 1e-12 # measured in [s]
    plot_time_interval::Float64 = 1e-8# measured in [s]
end

function validate_time_parameters(params)
    params.steady_state_tol > 0 || throw(ArgumentError("The steady-state tolerance must be larger than 0."))
    params.initial_timestep <= params.max_timestep || throw(ArgumentError("The initial timestep must be <=  the maximum timestep."))
end

"""
take the parameters from the parameter file and use them to configure the solver
"""
function setup_time_parameters!(parameters)
    params = parameters.time_parameters
    L_REF = parameters.non_dim.L_REF
    D_REF = parameters.non_dim.D_REF
    T_REF = L_REF^2 / D_REF
    params.final_time = params.final_time / T_REF
    params.max_timestep = params.max_timestep / T_REF
    params.initial_timestep = 1e-6#params.initial_timestep / T_REF
    params.plot_time_interval = params.plot_time_interval / T_REF
end
