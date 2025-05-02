# 
# Copyright (C) 2025 - 2025 by the authors of the SMPNP code.
# 
# This file is part of SMPNP
# 

# configuration
ALLOWED_PNP_MODES = ["PNP", "MPNP", "SMPNP"]
"""
function to validate the PNP mode is a valid PNP mode.
This function was generated with the help of ChatGPT (OpenAI).
"""
function validate_pnp_type(mode::String)
    mode in ALLOWED_PNP_MODES || throw(ArgumentError("You choose an invalid PNP_mode throug the parameter file. Please choose one of the available modes. Available Modes: $mode. Allowed values: $(join(ALLOWED_PNP_MODES, ", "))"))
end

"""
Here we define the struct of Species parameters with all the information 
about the different species.
"""
@option struct SpeciesParameters
    species::Vector{Int} = Vector{Int}(1:NUM_CONCENTRATION_SPECIES) # species without potential
    solvent_size::Float64 = 0.0
    ion_sizes::Vector{Float64} = zeros(Float64, NUM_CONCENTRATION_SPECIES)# ion sizes
    diffusivities::Vector{Float64} = zeros(Float64, NUM_CONCENTRATION_SPECIES) # the diffusivities are required, so we don't provide a default value 
    epsilon_r::Float64 = 80 # the relative permittivity
end

"""
asdf
This function was generated with the help of ChatGPT (OpenAI).
"""
function validate_species_parameters(species_params, pnp_type)
    are_sets_equal(species_params.species, CONCENTRATION_SPECIES_VECTOR) || throw(ArgumentError("The species vector must contain all the elements of the reference species vector in config.jl and have the same length (i.e. contain each species only once). Check your provided vector in the input file."))
    length(species_params.ion_sizes) == length(species_params.species) || throw(ArgumentError("The vector of ion sizes must have the same length as the vector of species."))
    for (ion_size, species) in zip(species_params.ion_sizes, species_params.species)
        if pnp_type == "PNP"
            species_params.solvent_size == 0.0 || throw(ArgumentError("The solvent size parameter will not have an effect since the pnp_type is chosen as 'PNP'."))
            ion_size == 0.0 || throw(ArgumentError("The ion size parameter will not have an effect since the pnp_type is chosen as 'PNP'."))
        elseif pnp_type == "SMPNP"
            species_params.solvent_size > 0.0 || throw(ArgumentError("The solvent size should be positive."))
            ion_size > 0.0 || throw(ArgumentError("The ion size should be positive."))
        end
    end

    length(species_params.diffusivities) == length(species_params.species) || throw(ArgumentError("The vector of diffusivities must have the same length as the vector of species."))
    for (diffusivity, species) in zip(species_params.diffusivities, species_params.species)
        diffusivity > 0.0 || throw(ArgumentError("The diffusivity size of species $species cannot be negative."))
    end

    species_params.epsilon_r > 0.0 || throw(ArgumentError("The relative permittivity cannot be negative or zero."))
end

