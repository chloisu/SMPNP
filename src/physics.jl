# 
# Copyright (C) 2025 - 2025 by the authors of the SMPNP code.
# 
# This file is part of SMPNP
# 

"""
    get_initial_timestep_system(grid, parameters)

returns a VoronoiFVM.System to solve the initial poisson 
equation for the potential distribution.
"""
function get_initial_timestep_system(grid, parameters)
    # we setup the physics for the poisson system only
    # we compute the prefactor for the Poisson equation
    L = parameters.non_dim.L_REF
    reference_permittivity = parameters.non_dim.EPSILON_REF * EPSILON_VAC
    non_dim_permittivity = parameters.species_parameters.epsilon_r * EPSILON_VAC / reference_permittivity
    additional_non_dim_prefactor = reference_permittivity / L^2 / parameters.non_dim.C_REF * non_dim_permittivity
    prefactor = POISSON_PHYS_PREFACTOR * additional_non_dim_prefactor
    # now we are ready to define the physics of this problem.
    physics = VoronoiFVM.Physics(;
        flux=function (f, u, edge, data)
            f[1] = prefactor * (u[1, 1] - u[1, 2])
            return nothing
        end
    )
    sys = VoronoiFVM.System(grid, physics; is_linear=false, unknown_storage=:sparse, assembly=:edgewise)
    # we add the potential 'species'
    enable_species!(sys, 1, [1])
    return sys
end



function get_time_dependent_system(grid, parameters=:nothing)
    # this system describes the phyiscs of the time-dependent PNP equations
    # the index of the different species are 
    # 1 = potential $\phi$
    # 2 = anion concentration $c_a$
    # 3 = cation concentration $c_c$

    # we compute the prefactor for the Poisson equation
    # we compute the prefactor in two parts, the constant one and the additional nondim parameters
    L = parameters.non_dim.L_REF
    reference_permittivity = parameters.non_dim.EPSILON_REF * EPSILON_VAC
    non_dim_permittivity = parameters.species_parameters.epsilon_r * EPSILON_VAC / reference_permittivity
    additional_non_dim_prefactor = reference_permittivity / L^2 / parameters.non_dim.C_REF * non_dim_permittivity
    prefactor = POISSON_PHYS_PREFACTOR * additional_non_dim_prefactor
    #1.0 / (4 * pi * L^2 * L_B * parameters.non_dim.C_REF)
    # we scale the diffusivities of the two concentration equations
    D_a_norm = parameters.species_parameters.diffusivities[ANION_EQ-1] / parameters.non_dim.D_REF
    D_c_norm = parameters.species_parameters.diffusivities[CATION_EQ-1] / parameters.non_dim.D_REF

    physics = VoronoiFVM.Physics(
        ; reaction=function (f, u, node, data)
            # source term of the poisson equation
            f[1] = -(Z_ANION * u[2] + Z_CATION * u[3])
            return nothing
        end,
        flux=function (f, u, edge, data)
            # compute the chemical potentials for the two cells
            μ_anion1 = μ(u[2, 1], u[3, 1], Z_ANION, u[1, 1], parameters)
            μ_anion2 = μ(u[2, 2], u[3, 2], Z_ANION, u[1, 2], parameters)

            μ_cation1 = μ(u[3, 1], u[2, 1], Z_CATION, u[1, 1], parameters)
            μ_cation2 = μ(u[3, 2], u[2, 2], Z_CATION, u[1, 2], parameters)
            # potential flux is simple poisson
            f[1] = prefactor * (u[1, 1] - u[1, 2])
            # anion flux is diffusion + migration
            c_anion_interface = 0.5 * (u[2, 1] + u[2, 2])
            f[2] = D_a_norm * c_anion_interface * (μ_anion1 - μ_anion2)#(((u[2, 1] - u[2, 2]) + c_anion_interface * Z_ANION * (u[1, 1] - u[1, 2])))
            # cation flux is diffusion + migration
            c_cation_interface = 0.5 * (u[3, 1] + u[3, 2])
            f[3] = D_c_norm * c_cation_interface * (μ_cation1 - μ_cation2)#((u[3, 1] - u[3, 2]) + c_cation_interface * Z_CATION * (u[1, 1] - u[1, 2]))
            return nothing
        end,
        storage=function (f, u, node, data)
            # time derivative for concentrations
            f[2] = u[2]
            f[3] = u[3]
            return nothing
        end
    )

    sys = VoronoiFVM.System(grid, physics; is_linear=false, unknown_storage=:sparse, assembly=:edgewise)
    enable_species!(sys, 1, [1]) # add potential
    enable_species!(sys, 2, [1]) # add anion
    enable_species!(sys, 3, [1]) # add cation
    return sys
end


function get_diffusion_flux_system(grid, parameters=:nothing)
    # this system describes the phyiscs of the time-dependent PNP equations
    # the index of the different species are 
    # 1 = potential $\phi$
    # 2 = anion concentration $c_a$
    # 3 = cation concentration $c_c$

    # we compute the prefactor for the Poisson equation
    L = parameters.non_dim.L_REF
    prefactor = 1.0 / (4 * pi * L^2 * L_B * parameters.non_dim.C_REF)
    # we scale the diffusivities of the two concentration equations
    D_a_norm = parameters.species_parameters.diffusivities[ANION_EQ-1] / parameters.non_dim.D_REF
    D_c_norm = parameters.species_parameters.diffusivities[CATION_EQ-1] / parameters.non_dim.D_REF

    physics = VoronoiFVM.Physics(
        ; reaction=function (f, u, node, data)
            # source term of the poisson equation
            f[1] = -(Z_ANION * u[2] + Z_CATION * u[3])
            return nothing
        end,
        flux=function (f, u, edge, data)
            # compute the chemical potentials for the two cells
            μ_anion1 = μ_due_to_diff(u[2, 1], u[3, 1], Z_ANION, u[1, 1], parameters)
            μ_anion2 = μ_due_to_diff(u[2, 2], u[3, 2], Z_ANION, u[1, 2], parameters)

            μ_cation1 = μ_due_to_diff(u[3, 1], u[2, 1], Z_CATION, u[1, 1], parameters)
            μ_cation2 = μ_due_to_diff(u[3, 2], u[2, 2], Z_CATION, u[1, 2], parameters)
            # potential flux is simple poisson
            f[1] = prefactor * (u[1, 1] - u[1, 2])
            # anion flux is diffusion + migration
            c_anion_interface = 0.5 * (u[2, 1] + u[2, 2])
            f[2] = D_a_norm * c_anion_interface * (μ_anion1 - μ_anion2)#(((u[2, 1] - u[2, 2]) + c_anion_interface * Z_ANION * (u[1, 1] - u[1, 2])))
            # cation flux is diffusion + migration
            c_cation_interface = 0.5 * (u[3, 1] + u[3, 2])
            f[3] = D_c_norm * c_cation_interface * (μ_cation1 - μ_cation2)#((u[3, 1] - u[3, 2]) + c_cation_interface * Z_CATION * (u[1, 1] - u[1, 2]))
            return nothing
        end,
        storage=function (f, u, node, data)
            # time derivative for concentrations
            f[2] = u[2]
            f[3] = u[3]
            return nothing
        end
    )

    sys = VoronoiFVM.System(grid, physics; is_linear=false, unknown_storage=:sparse, assembly=:edgewise)
    enable_species!(sys, 1, [1]) # add potential
    enable_species!(sys, 2, [1]) # add anion
    enable_species!(sys, 3, [1]) # add cation
    return sys
end

function get_size_flux_system(grid, parameters=:nothing)
    # this system describes the phyiscs of the time-dependent PNP equations
    # the index of the different species are 
    # 1 = potential $\phi$
    # 2 = anion concentration $c_a$
    # 3 = cation concentration $c_c$

    # we compute the prefactor for the Poisson equation
    L = parameters.non_dim.L_REF
    prefactor = 1.0 / (4 * pi * L^2 * L_B * parameters.non_dim.C_REF)
    # we scale the diffusivities of the two concentration equations
    D_a_norm = parameters.species_parameters.diffusivities[ANION_EQ-1] / parameters.non_dim.D_REF
    D_c_norm = parameters.species_parameters.diffusivities[CATION_EQ-1] / parameters.non_dim.D_REF

    physics = VoronoiFVM.Physics(
        ; reaction=function (f, u, node, data)
            # source term of the poisson equation
            f[1] = -(Z_ANION * u[2] + Z_CATION * u[3])
            return nothing
        end,
        flux=function (f, u, edge, data)
            # compute the chemical potentials for the two cells
            μ_anion1 = μ_due_to_finite_size(u[2, 1], u[3, 1], Z_ANION, u[1, 1], parameters)
            μ_anion2 = μ_due_to_finite_size(u[2, 2], u[3, 2], Z_ANION, u[1, 2], parameters)

            μ_cation1 = μ_due_to_finite_size(u[3, 1], u[2, 1], Z_CATION, u[1, 1], parameters)
            μ_cation2 = μ_due_to_finite_size(u[3, 2], u[2, 2], Z_CATION, u[1, 2], parameters)
            # potential flux is simple poisson
            f[1] = prefactor * (u[1, 1] - u[1, 2])
            # anion flux is diffusion + migration
            c_anion_interface = 0.5 * (u[2, 1] + u[2, 2])
            f[2] = D_a_norm * c_anion_interface * (μ_anion1 - μ_anion2)#(((u[2, 1] - u[2, 2]) + c_anion_interface * Z_ANION * (u[1, 1] - u[1, 2])))
            # cation flux is diffusion + migration
            c_cation_interface = 0.5 * (u[3, 1] + u[3, 2])
            f[3] = D_c_norm * c_cation_interface * (μ_cation1 - μ_cation2)#((u[3, 1] - u[3, 2]) + c_cation_interface * Z_CATION * (u[1, 1] - u[1, 2]))
            return nothing
        end,
        storage=function (f, u, node, data)
            # time derivative for concentrations
            f[2] = u[2]
            f[3] = u[3]
            return nothing
        end
    )

    sys = VoronoiFVM.System(grid, physics; is_linear=false, unknown_storage=:sparse, assembly=:edgewise)
    enable_species!(sys, 1, [1]) # add potential
    enable_species!(sys, 2, [1]) # add anion
    enable_species!(sys, 3, [1]) # add cation
    return sys
end

function get_potential_flux_system(grid, parameters=:nothing)
    # this system describes the phyiscs of the time-dependent PNP equations
    # the index of the different species are 
    # 1 = potential $\phi$
    # 2 = anion concentration $c_a$
    # 3 = cation concentration $c_c$

    # we compute the prefactor for the Poisson equation
    L = parameters.non_dim.L_REF
    prefactor = 1.0 / (4 * pi * L^2 * L_B * parameters.non_dim.C_REF)
    # we scale the diffusivities of the two concentration equations
    D_a_norm = parameters.species_parameters.diffusivities[ANION_EQ-1] / parameters.non_dim.D_REF
    D_c_norm = parameters.species_parameters.diffusivities[CATION_EQ-1] / parameters.non_dim.D_REF

    physics = VoronoiFVM.Physics(
        ; reaction=function (f, u, node, data)
            # source term of the poisson equation
            f[1] = -(Z_ANION * u[2] + Z_CATION * u[3])
            return nothing
        end,
        flux=function (f, u, edge, data)
            # compute the chemical potentials for the two cells
            μ_anion1 = μ_due_to_potential(u[2, 1], u[3, 1], Z_ANION, u[1, 1], parameters)
            μ_anion2 = μ_due_to_potential(u[2, 2], u[3, 2], Z_ANION, u[1, 2], parameters)

            μ_cation1 = μ_due_to_potential(u[3, 1], u[2, 1], Z_CATION, u[1, 1], parameters)
            μ_cation2 = μ_due_to_potential(u[3, 2], u[2, 2], Z_CATION, u[1, 2], parameters)
            # potential flux is simple poisson
            f[1] = prefactor * (u[1, 1] - u[1, 2])
            # anion flux is diffusion + migration
            c_anion_interface = 0.5 * (u[2, 1] + u[2, 2])
            f[2] = D_a_norm * c_anion_interface * (μ_anion1 - μ_anion2)#(((u[2, 1] - u[2, 2]) + c_anion_interface * Z_ANION * (u[1, 1] - u[1, 2])))
            # cation flux is diffusion + migration
            c_cation_interface = 0.5 * (u[3, 1] + u[3, 2])
            f[3] = D_c_norm * c_cation_interface * (μ_cation1 - μ_cation2)#((u[3, 1] - u[3, 2]) + c_cation_interface * Z_CATION * (u[1, 1] - u[1, 2]))
            return nothing
        end,
        storage=function (f, u, node, data)
            # time derivative for concentrations
            f[2] = u[2]
            f[3] = u[3]
            return nothing
        end
    )

    sys = VoronoiFVM.System(grid, physics; is_linear=false, unknown_storage=:sparse, assembly=:edgewise)
    enable_species!(sys, 1, [1]) # add potential
    enable_species!(sys, 2, [1]) # add anion
    enable_species!(sys, 3, [1]) # add cation
    return sys
end

