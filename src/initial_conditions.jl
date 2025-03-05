# 
# Copyright (C) 2025 - 2025 by the authors of the SMPNP code.
# 
# This file is part of SMPNP
# 


@option struct InitialConditions
    species::Vector{Int} = Vector{Int}(1:NUM_SPECIES) # initial species
    values::Vector{Float64} = DEFAULT_INITIAL_CONDITIONS # initial conditions to apply
end

function apply_initial_conditions_for_intial_timestep!(sys, parameters)
    initial_condition_parameters = parameters.initial_conditions
    # we only want the boundary conditions for the PHI_EQ
    # so we first check if there is an initial conditions specified
    if PHI_EQ ∈ initial_condition_parameters.species
        # we know it's in there but we don't know the exact index, let get it
        idx = findfirst(==(PHI_EQ), initial_condition_parameters.species)
        # we can now set the initial condition value value
        inival = unknowns(sys)
        inival .= initial_condition_parameters.values[idx]
    else
        throw(ParameterError("You need to specify either all or no initial condition in the parameter file. You can currently not specify initial conditions only for a subset of species."))
    end
end

function apply_initial_condition_for_time_dependent!(sys, parameters, potential_from_initial_timestep)
    initial_condition_parameters = parameters.initial_conditions
    inival = unknowns(sys)
    # we loop over all the different species that we have
    if length(initial_condition_parameters.species) != NUM_SPECIES
        throw(ParameterError("You got this error becaues the number of species for which initial conditions are defined is different from  NUM_SPECIES in config.jl. You need to specify either all or no initial condition in the parameter file. You can currently not specify initial conditions only for a subset of species."))
    end
    for species in initial_condition_parameters.species
        # if we get to the PHI_EQ, we will not use the initial condition from the parameter
        # file again but rather specify the result of the initial "timestep" to be the initial condition.
        if species == PHI_EQ
            inival[PHI_EQ, :] = potential_from_initial_timestep
        else
            inival[species, :] .= initial_condition_parameters.values[species]
        end
    end
    return inival
end