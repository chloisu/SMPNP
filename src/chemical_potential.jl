"""
this function returns the chemical potential for a given species 

"""
function μ(c1, c2, z1, ϕ, parameters::Parameters)

    if parameters.pnp_type == "PNP"
        μ1 = rlog(c1) + z1 * ϕ
    end
    #if parameters.PNP_mode == "MPNP"
    #
    #end
    #
    if parameters.pnp_type == "SMPNP"
        m0 = parameters.species_parameters.solvent_size^3 * MOL_PER_LITER_TO_PER_CUBIC_METER
        if z1 == Z_ANION
            m1 = parameters.species_parameters.ion_sizes[ANION_EQ-1]^3 * MOL_PER_LITER_TO_PER_CUBIC_METER
            m2 = parameters.species_parameters.ion_sizes[CATION_EQ-1]^3 * MOL_PER_LITER_TO_PER_CUBIC_METER
        end
        if z1 == Z_CATION
            m1 = parameters.species_parameters.ion_sizes[CATION_EQ-1]^3 * MOL_PER_LITER_TO_PER_CUBIC_METER
            m2 = parameters.species_parameters.ion_sizes[ANION_EQ-1]^3 * MOL_PER_LITER_TO_PER_CUBIC_METER
        end
        sum = 1 - c1 * m1 - c2 * m2
        μ1 = m1 / m0 * (rlog(c1 * m1) - rlog(sum)) + z1 * ϕ
    end
    return μ1
end

function μ_due_to_diff(c1, c2, z1, ϕ, parameters::Parameters)
    if parameters.pnp_type == "SMPNP"
        m0 = parameters.species_parameters.solvent_size^3 * MOL_PER_LITER_TO_PER_CUBIC_METER
        if z1 == Z_ANION
            m1 = parameters.species_parameters.ion_sizes[ANION_EQ-1]^3 * MOL_PER_LITER_TO_PER_CUBIC_METER
            m2 = parameters.species_parameters.ion_sizes[CATION_EQ-1]^3 * MOL_PER_LITER_TO_PER_CUBIC_METER
        end
        if z1 == Z_CATION
            m1 = parameters.species_parameters.ion_sizes[CATION_EQ-1]^3 * MOL_PER_LITER_TO_PER_CUBIC_METER
            m2 = parameters.species_parameters.ion_sizes[ANION_EQ-1]^3 * MOL_PER_LITER_TO_PER_CUBIC_METER
        end
        sum = 1 - c1 * m1 - c2 * m2
        μ1 = m1 / m0 * (rlog(c1 * m1))
    end
    return μ1
end

function μ_due_to_potential(c1, c2, z1, ϕ, parameters::Parameters)
    if parameters.pnp_type == "SMPNP"
        m0 = parameters.species_parameters.solvent_size^3 * MOL_PER_LITER_TO_PER_CUBIC_METER
        if z1 == Z_ANION
            m1 = parameters.species_parameters.ion_sizes[ANION_EQ-1]^3 * MOL_PER_LITER_TO_PER_CUBIC_METER
            m2 = parameters.species_parameters.ion_sizes[CATION_EQ-1]^3 * MOL_PER_LITER_TO_PER_CUBIC_METER
        end
        if z1 == Z_CATION
            m1 = parameters.species_parameters.ion_sizes[CATION_EQ-1]^3 * MOL_PER_LITER_TO_PER_CUBIC_METER
            m2 = parameters.species_parameters.ion_sizes[ANION_EQ-1]^3 * MOL_PER_LITER_TO_PER_CUBIC_METER
        end
        sum = 1 - c1 * m1 - c2 * m2
        μ1 = z1 * ϕ
    end
    return μ1
end

function μ_due_to_finite_size(c1, c2, z1, ϕ, parameters::Parameters)
    if parameters.pnp_type == "SMPNP"
        m0 = parameters.species_parameters.solvent_size^3 * MOL_PER_LITER_TO_PER_CUBIC_METER
        if z1 == Z_ANION
            m1 = parameters.species_parameters.ion_sizes[ANION_EQ-1]^3 * MOL_PER_LITER_TO_PER_CUBIC_METER
            m2 = parameters.species_parameters.ion_sizes[CATION_EQ-1]^3 * MOL_PER_LITER_TO_PER_CUBIC_METER
        end
        if z1 == Z_CATION
            m1 = parameters.species_parameters.ion_sizes[CATION_EQ-1]^3 * MOL_PER_LITER_TO_PER_CUBIC_METER
            m2 = parameters.species_parameters.ion_sizes[ANION_EQ-1]^3 * MOL_PER_LITER_TO_PER_CUBIC_METER
        end
        sum = 1 - c1 * m1 - c2 * m2
        μ1 = m1 / m0 * (-rlog(sum))
    end
    return μ1
end