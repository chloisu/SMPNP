using Printf
using VoronoiFVM
using ExtendableGrids
using ExtendableSparse: ILUZeroPreconBuilder
using GridVisualize
using LinearSolve
using ILUZero
using Triangulate
using SimplexGridFactory
using LinearAlgebra

function generate_grid()
    builder = SimplexGridBuilder(; Generator=Triangulate)
    cellregion!(builder, 1)
    maxvolume!(builder, 0.01)
    regionpoint!(builder, 0.5, 0.5)
    p1 = point!(builder, 0, 0)
    p2 = point!(builder, 1, 0)
    p3 = point!(builder, 1, 0.4)
    p4 = point!(builder, 3, 0.4)
    p5 = point!(builder, 3, 0)
    p6 = point!(builder, 4, 0)
    p7 = point!(builder, 4, 0.5)
    p8 = point!(builder, 0, 0.5)
    # left reservoir
    facetregion!(builder, 1)
    facet!(builder, p8, p1)
    # bottom reservoir
    facetregion!(builder, 2)
    facet!(builder, p1, p2)
    facet!(builder, p5, p6)
    # pore wall
    facetregion!(builder, 3)
    facet!(builder, p2, p3)
    facet!(builder, p3, p4)
    facet!(builder, p4, p5)
    # right reservoir
    facetregion!(builder, 4)
    facet!(builder, p6, p7)
    # symmetry line
    facetregion!(builder, 5)
    facet!(builder, p7, p8)


    function unsuitable(x1, y1, x2, y2, x3, y3, area)
        bary = [(x1 + x2 + x3) / 3, (y1 + y2 + y3) / 3]
        min_x = 1
        max_x = 3
        rf_x = max(min_x, min(bary[1], max_x))
        refinement_center = [rf_x, 0.4]
        dist = norm(bary - refinement_center)
        if area > 0.001 * dist
            return 1
        else
            return 0
        end
    end
    options!(builder; unsuitable=unsuitable)
    grid = simplexgrid(builder)
    print(num_cells(grid))
    return grid
end


function main(;
    n=10, Plotter=nothing, verbose=false, unknown_storage=:sparse,
    method_linear=nothing, assembly=:edgewise
)
    h = 1.0 / convert(Float64, n)
    X = collect(0.0:h:1.0)
    Y = collect(0.0:h:1.0)

    grid = generate_grid()

    physics = VoronoiFVM.Physics(;
        reaction=function (f, u, node, data)
            f[1] = 0
            return nothing
        end, flux=function (f, u, edge, data)
            f[1] = 0.001 * (u[1, 1] - u[1, 2])
            return nothing
        end, source=function (f, node, data)
            x1 = node[1] - 0.5
            x2 = node[2] - 0.5
            f[1] = 0.0#exp(-20.0 * (x1^2 + x2^2))
            return nothing
        end, storage=function (f, u, node, data)
            f[1] = u[1]
            return nothing
        end
    )
    sys = VoronoiFVM.System(grid, physics; unknown_storage, assembly=assembly)
    enable_species!(sys, 1, [1])

    boundary_dirichlet!(sys, 1, 1, 1)
    boundary_dirichlet!(sys, 1, 4, 0)

    inival = unknowns(sys)
    inival .= 0.5

    control = VoronoiFVM.NewtonControl()
    control.verbose = verbose
    control.reltol_linear = 1.0e-5
    control.method_linear = method_linear
    tstep = 1
    time = 0
    U = solve(sys; inival, control, tstep)
    p = GridVisualizer(;
        Plotter,
        resolution=(600, 600),
        clear=true)

    while time < 100
        time = time + tstep
        U = solve(sys; inival, control, tstep)
        u15 = U[15]
        inival .= U
        tstep *= 1.0
    end
    scalarplot!(p[1, 1], grid, U[1, :], clear=true, show=true)
    return reveal(p)
end

using GLMakie
p = main(Plotter=GLMakie, verbose=true)
GLMakie.save(joinpath(".", "out.png"), p)  #hide