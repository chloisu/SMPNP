using ExtendableGrids
using GridVisualize
include("../src/constants.jl")


function grid1d(; n=50)
    X = collect(0:(1/n):1)
    return g = simplexgrid(X)
end

function potential_pb_1d(grid, wall_potential, r)
    g = grid
    lambda_d = LAMBDA_D / r
    a = exp(wall_potential / 2) + 1
    b = exp(wall_potential / 2) - 1
    return map(x -> 2 * log((a + b * exp(-x / lambda_d)) / (a - b * exp(-x / lambda_d))), g)
end

function concentrations_pb_1d(grid, wall_potential, r)
    potential_fun = potential_pb_1d(grid, wall_potential, r)
    ca = map(x -> exp(-Z_ANION * x), potential_fun)
    cc = map(x -> exp(-Z_CATION * x), potential_fun)
    return ca, cc
end

function concentrations_smpnp_1d(potential_array, aa, ac, a0, ca_centerline, cc_centerline)
    ma = aa^3 * MOL_PER_LITER_TO_PER_CUBIC_METER
    mc = ac^3 * MOL_PER_LITER_TO_PER_CUBIC_METER
    m0 = a0^3 * MOL_PER_LITER_TO_PER_CUBIC_METER

    xa_centerline = ca_centerline * ma
    xc_centerline = cc_centerline * mc

    Ca = ma / m0 * log(xa_centerline / (1 - xa_centerline - xc_centerline)) + Z_ANION * potential_array[end]
    Cc = mc / m0 * log(xc_centerline / (1 - xa_centerline - xc_centerline)) + Z_CATION * potential_array[end]

    Aa = exp.((Ca .+ potential_array) * ma / m0)
    Ac = exp.((Cc .- potential_array) * mc / m0)

    ca = 1.0 ./ ma .* (Aa - Aa .* Ac ./ (1 .+ Ac)) ./ (1 .+ Aa .- Aa .* Ac ./ (1 .+ Ac))
    cc = 1.0 ./ mc .* (Ac - Aa .* Ac ./ (1 .+ Aa)) ./ (1 .+ Ac .- Aa .* Ac ./ (1 .+ Aa))
    return ca, cc
end