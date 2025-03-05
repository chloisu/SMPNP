# 
# Copyright (C) 2025 - 2025 by the authors of the SMPNP code.
# 
# This file is part of SMPNP
# 

@option struct Parameters
    case::String = "h_20nm_r10nm_l50nmh20h" # name of the case 
    # the type of the PNP system. There is three different options to choose from
    # PNP : the classical Poisson-Nernst-Planck system without taking the size of the ions into account.
    # MPNP: the size of the ions is taken into account but with a single constant size parameters (i.e. both ions are assumed to have the same size.)
    # SMPNP: each of the two ions has a differnt size
    PNP_mode::String = "SMPNP" |> validate_PNP_mode
    grid_type::String = "equally_spaced_1d"
    grid::GridParameters = GridParameters()
    cation_size::Float64 = 5e-10 # cation size [meter]
    anion_size::Float64 = 5e-10 # cation size [meter]
    solvent_molecule_size::Float64 = 5e-10
    D_anion::Float64 = 5e-9 # anion diffusivity [m^2/s]
    D_cation::Float64 = 5e-9 # cation diffusivity [m^2/s]
    surface_potential::Float64 = 0.3# surface potential [V]
    output_directory::String = ""
    # the following section describes the references scales
    # used for the non-dimensionalization of the equations
    # based on these quantities, we can then also run the post-processing
    # steps to convert back any non-dimensional results to 
    # dimensional units.
    non_dim::NonDimensionalization = NonDimensionalization()
    boundary_conditions::BoundaryConditions = BoundaryConditions()
    initial_conditions::InitialConditions = InitialConditions()
    solver_parameters::SolverParameters = SolverParameters()
    time_parameters::TimeParameters = TimeParameters()
    #L_REF::Float64 = pore_length # reference lengthscale in [meter]
    #C_REF::Float64 = 1 # reference concentration in [mol / L]
    #D_REF::Float64 = 5e-9 # reference diffusivity [m^2/s]
    #T_REF::Float64 = L_REF^2 / D_REF
    #PHI_REF::Float64 = K_B * T / E_CHARGE
    # the following section allows the customization of the linear and 
    # non-linear solvers. The arguments provided here are going to be passed 
    # on to the Solver Control object of VoronoiFVM.
    #max_iter::Int = 20 # maximum nonlinear Newton iterations 
    #dt::Float64 = 1e-6 # timestep
    #dt_max::Float64 = 2e-1 # maximum timestep
    #t_final::Float64 = 1 # final time
end

"""
    load_config_file(filename)

loads the parameter file with the given filename and parses the parameters in the file.
Any parameters not present in the input file will get their default value according
to the definition above.
"""
function load_config_file(filename)
    # Define a struct for your configuration
    file = YAML.load_file(filename; dicttype=Dict{String,Any})
    parameters = from_dict(Parameters, file) # parse input parameters
    return parameters
end
