include("../src/smpnp.jl")
include("analytical_solutions.jl")
using CairoMakie
# This setup was generated with the help of ChatGPT (OpenAI).
# now that we have main available, lets run it and get the steady state solution
global ARGS = ["verification_simple_pnp_different_epsilon_r.yml"]  # Simulate command-line arguments
# first we compute the numerical solution 
grid, U = Base.invokelatest(smpnp)  # Call main() after modifying ARGS
phi_numerical = U[1, :]
ca_numerical = U[2, :]
cc_numerical = U[3, :]
# second, we load the analytical solutions
phi_analytical = potential_pb_1d_v2(grid, 0.3, 10e-9; epsilon_r=31)
ca_analytical, cc_analytical = concentrations_pb_1d(grid, 0.3, 10e-9; epsilon_r=31)

p = GridVisualizer(;
    Plotter=CairoMakie,
    layout=(3, 1))
scalarplot!(p[1, 1], grid, phi_analytical, clear=false, show=true)
scalarplot!(p[1, 1], grid, phi_numerical, clear=false, show=true, linestyle=:dash, color=(1, 0, 0))

scalarplot!(p[2, 1], grid, ca_analytical, clear=false, show=true)
scalarplot!(p[2, 1], grid, ca_numerical, clear=false, show=true, linestyle=:dash, color=(1, 0, 0))

scalarplot!(p[3, 1], grid, cc_analytical, clear=false, show=true)
scalarplot!(p[3, 1], grid, cc_numerical, clear=false, show=true, linestyle=:dash, color=(1, 0, 0))

CairoMakie.save(joinpath(".", "verification_simple_pnp_different_epsilon_r.pdf"), reveal(p))  #hide