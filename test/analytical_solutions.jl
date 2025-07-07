using ExtendableGrids
using GridVisualize
include("../src/constants.jl")


function grid1d(; n=50)
    X = collect(0:(1/n):1)
    return g = simplexgrid(X)
end

function potential_pb_1d(grid, wall_potential, r)
    g = grid
    lambda_d = sqrt(EPSILON_VAC * 80 * K_B * T / E_CHARGE^2 / 2 / MOL_PER_LITER_TO_PER_CUBIC_METER) / r
    a = exp(wall_potential / 2) + 1
    b = exp(wall_potential / 2) - 1
    return map(x -> 2 * log((a + b * exp(-x / lambda_d)) / (a - b * exp(-x / lambda_d))), g)
end

"""
Analytical potential expression taken from 
Hunter, R. J. (2013). Zeta potential in colloid science: 
principles and applications (Vol. 2). Academic press.
Kappa is compute from equation (2.3.9) and normalized with the radius here.
The potential expression is then given by (2.3.16)
"""
function potential_pb_1d_v2(grid, wall_potential, r, C_REF=MOL_PER_LITER_TO_PER_CUBIC_METER)
    g = grid
    kappa = 1.0 / (sqrt(EPSILON_VAC * 80 * K_B * T / E_CHARGE^2 / 2 / C_REF) / r)
    print("Kappa is " * string(kappa))
    return map(x -> 4 * atanh(tanh(wall_potential / 4) * exp(-x * kappa)), g)
end

function concentrations_pb_1d(grid, wall_potential, r, C_REF=MOL_PER_LITER_TO_PER_CUBIC_METER)
    potential_fun = potential_pb_1d_v2(grid, wall_potential, r, C_REF)
    ca = map(x -> exp(-Z_ANION * x), potential_fun)
    cc = map(x -> exp(-Z_CATION * x), potential_fun)
    return ca, cc
end

function concentrations_smpnp_1d(potential_array, aa, ac, a0, ca_centerline, cc_centerline, ka=nothing, kc=nothing)
    # we first make sure that in case the user provides ka and kc that both are provided
    if (ka !== nothing && kc === nothing)
        throw(AssertionError("If ka is provided, kc must be provided too."))
    end
    if (kc !== nothing && ka === nothing)
        throw(AssertionError("If kc is provided, ka must be provided too."))
    end
    # we also check that both ka and kc come from one of the implemented cases
    ka ∈ (nothing, 1, 2) || throw(AssertionError("The functionality for the provided ka is not implemented yet."))
    kc ∈ (nothing, 1, 2) || throw(AssertionError("The functionality for the provided ka is not implemented yet."))
    # at this point both ka and kc are provided or not provided. So for the case that both are provided,
    # we assert that the provided ion radii match the provided k's.
    #if (kc !== nothing)
    #    isapprox(aa^3 / a0^3, ka) || throw(AssertionError("The provided ka is not close to the ka computed from the provided aa and a0"))
    #    isapprox(ac^3 / a0^3, kc) || throw(AssertionError("The provided kc is not close to the kc computed from the provided ac and a0"))
    #end
    ma = aa^3 * MOL_PER_LITER_TO_PER_CUBIC_METER
    mc = ac^3 * MOL_PER_LITER_TO_PER_CUBIC_METER
    m0 = a0^3 * MOL_PER_LITER_TO_PER_CUBIC_METER

    ζa_centerline = ca_centerline * ma
    ζc_centerline = cc_centerline * mc

    Ca = log(ζa_centerline / (1 - ζa_centerline - ζc_centerline)^(ma / m0)) + Z_ANION * potential_array[end]
    Cc = log(ζc_centerline / (1 - ζa_centerline - ζc_centerline)^(mc / m0)) + Z_CATION * potential_array[end]

    Ka = exp.((Ca .+ potential_array))
    Kc = exp.((Cc .- potential_array))

    # now we have to distinguish the differnt cases.
    # case with equal ion sizes and equal to solvent size
    # this is also the default case if not ka and kc are provided
    if (ka === nothing || (ka == 1 && kc == 1))
        S = 1 ./ (Ka .+ Kc .+ 1)
        ζa = Ka .* S
        ζc = Kc .* S
    end
    # case with anion larger than cation. cation equal to solvent size
    if (ka == 2 && kc == 1)
        S = (-(Kc .+ 1) + sqrt.((Kc .+ 1) .^ 2 + 4 .* Ka)) ./ (2 .* Ka)
        ζa = Ka .* S .^ ka
        ζc = Kc .* S .^ kc
    end
    ca = ζa / ma
    cc = ζc / mc
    return ca, cc
end


