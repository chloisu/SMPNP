using ExtendableGrids
using GridVisualize
include("constants.jl")


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

#function plotting_func1d(; Plotter=default_plotter(), kwargs...)
#    g, f = potential_pb_1d(0.4, 5e-9)
#    return scalarplot(g, f; Plotter=Plotter, resolution=(500, 300), kwargs...)
#end

#using GLMakie
#p = plotting_func1d(Plotter=GLMakie, verbose=true)
#GLMakie.save(joinpath(".", "db.jpg"), p)  #hide