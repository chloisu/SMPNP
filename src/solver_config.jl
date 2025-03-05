"""
Here we define the struct of Grid Parameters which contains all the parameters of the different grid types.
This is a small wrapper to make the solver arguments of VoronoiFVM.SolverControl() available through the 
parameter file.
"""
@option struct SolverParameters
    verbose::Union{Bool,String} = false
    abstol::Float64 = 1.0e-10
    reltol::Float64 = 1.0e-10
    maxiters::Int = 100
    reltol_linear::Float64 = 1.0e-4
    abstol_linear::Float64 = 1.0e-8
    method_linear::Union{Nothing,LinearSolve.SciMLLinearSolveAlgorithm} = nothing
end

"""
take the parameters from the parameter file and use them to configure the solver
"""
function setup_solver_control!(control, parameters)
    params = parameters.solver_parameters
    control.verbose::Union{Bool,String} = params.verbose
    control.abstol::Float64 = params.abstol
    control.reltol::Float64 = params.reltol
    control.maxiters::Int = params.maxiters
    control.reltol_linear::Float64 = params.reltol_linear
    control.abstol_linear::Float64 = params.abstol_linear
    control.method_linear::Union{Nothing,LinearSolve.SciMLLinearSolveAlgorithm} = params.method_linear
end
