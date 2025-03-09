# 
# Copyright (C) 2025 - 2025 by the authors of the SMPNP code.
# 
# This file is part of SMPNP
# 


@option struct DirichletBoundary
    species::Vector{Int} = Vector{Int}()
    boundary_indices::Vector{Vector{Int}} = Vector{Vector{Int}}() # indicies to apply the boundary conditions
    boundary_values::Vector{Vector{Float64}} = Vector{Vector{Float64}}() # values to apply at these boundaries
end

function validate_dirichlet_boundary_conditions(params)
    for species in params.species
        if species < 1
            throw(ArgumentError("The species $species is not a valid species. Indicies must be positive and taken from the SPECIES_VECTOR in config.jl"))
        end
        # it must be in the list of default species
        if species ∉ SPECIES_VECTOR
            throw(ArgumentError("The species $species cannot be found in the SPECIES_VECTOR of config.jl."))
        end
    end
    # TODO : add validation for the boundary indices that they are valid
    # simply check that boundary_values have same shape as boundary_indices
    size(params.boundary_indices) == size(params.boundary_values) || throw(AssertionError("The size of the indices and values vectors have to be the same."))
end

function apply_dirichlet_for_initial_timestep!(sys, parameters)
    dirichlet_parameters = parameters.boundary_conditions.dirichlet
    # we only want the boundary conditions for the PHI_EQ
    # so we first check if there are any boundary conditions
    if PHI_EQ ∈ dirichlet_parameters.species
        # know we need to apply the correct boundary conditions,
        # so we loop over the boundary_indices and boundary values
        for (id, value) in zip(dirichlet_parameters.boundary_indices[PHI_EQ], dirichlet_parameters.boundary_values[PHI_EQ])
            boundary_dirichlet!(sys, PHI_EQ, id, value)
        end
    end
end


function apply_dirichlet_for_time_dependent!(sys, parameters)
    dirichlet_parameters = parameters.boundary_conditions.dirichlet
    # we loop over all the different species that we have
    for species in dirichlet_parameters.species
        # know we need to apply the correct boundary conditions,
        # so we loop over the boundary_indices and boundary values
        println("Boundary conditions")

        for (id, value) in zip(dirichlet_parameters.boundary_indices[species], dirichlet_parameters.boundary_values[species])
            print("Species " * string(species))
            print("ID" * string(id))
            print("value " * string(value))
            println("")
            boundary_dirichlet!(sys, species, id, value)
        end
    end
end


"""
Here we define the struct of Grid Parameters which contains all the parameters of the different grid types.
"""
@option struct BoundaryConditions
    dirichlet::DirichletBoundary = DirichletBoundary()
end
