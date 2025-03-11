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
    pnp_type::String = "SMPNP"
    species_parameters::SpeciesParameters = SpeciesParameters()
    grid_type::String = "equally_spaced_1d"
    grid::GridParameters = GridParameters()
    #pore_radius::Float64 = 10e-9 # pore radius [meter]
    #pore_length::Float64 = 50e-9 # pore length [meter]
    #reservoir_height::Float64 = 80e-9
    #D_REF::Float64 = 5e-9
    #cation_size::Float64 = 5e-10 # cation size [meter]
    #anion_size::Float64 = 5e-10 # cation size [meter]
    ##solvent_molecule_size::Float64 = 5e-10
    #D_anion::Float64 = 5e-9 # anion diffusivity [m^2/s]
    #D_cation::Float64 = 5e-9 # cation diffusivity [m^2/s]
    #surface_potential::Float64 = 0.3# surface potential [V]
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
    if check_parameters(parameters)
        return parameters
    else
        print("error")
    end
end

"""
This function verifies the provided input parameters. It makes sure that the input arguments
    provided are valid. It also provided error messages if they are not.
"""
function check_parameters(parameters)
    # Validate the fields after parsing
    validate_pnp_type(parameters.pnp_type)
    validate_species_parameters(parameters.species_parameters, parameters.pnp_type)
    validate_grid_type(parameters.grid_type)
    validate_grid_parameters(parameters.grid_type, parameters)
    validate_non_dim(parameters.non_dim)
    validate_dirichlet_boundary_conditions(parameters.boundary_conditions.dirichlet)
    validate_initial_conditions(parameters.initial_conditions)
    validate_time_parameters(parameters.time_parameters)
    return true
end
