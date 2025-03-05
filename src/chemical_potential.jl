"""
this function returns the chemical potential for a given species 

"""
function μ(c1, c2, z1, ϕ, parameters::Parameters)

    if parameters.PNP_mode == "PNP"
        μ1 = rlog(c1) + z1 * ϕ
    end
    #if parameters.PNP_mode == "MPNP"
    #
    #end
    #
    if parameters.PNP_mode == "SMPNP"
        m0 = parameters.solvent_molecule_size^3 * MOL_PER_LITER_TO_PER_CUBIC_METER
        if z1 == Z_ANION
            m1 = parameters.anion_size^3 * MOL_PER_LITER_TO_PER_CUBIC_METER
            m2 = parameters.cation_size^3 * MOL_PER_LITER_TO_PER_CUBIC_METER
        end
        if z1 == Z_CATION
            m1 = parameters.cation_size^3 * MOL_PER_LITER_TO_PER_CUBIC_METER
            m2 = parameters.anion_size^3 * MOL_PER_LITER_TO_PER_CUBIC_METER
        end
        sum = 1 - c1 * m1 - c2 * m2
        μ1 = m1 / m0 * (rlog(c1 * m1) - rlog(sum)) + z1 * ϕ
    end
    return μ1
end