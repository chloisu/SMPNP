include("../src/smpnp.jl")
include("analytical_solutions.jl")
include("helpers.jl")
using CairoMakie
# This setup was generated with the help of ChatGPT (OpenAI).
# now that we have main available, lets run it and get the steady state solution
global ARGS = ["verification_simple_pnp_dimensional_residual_check.yml"]  # Simulate command-line arguments
# first we compute the numerical solution 
grid, U = Base.invokelatest(smpnp)  # Call main() after modifying ARGS
phi_numerical = U[1, :]
ca_numerical = U[2, :]
cc_numerical = U[3, :]
# second, we load the analytical solutions
phi_analytical = potential_pb_1d_v2(grid, 0.3, 10e-9)
ca_analytical, cc_analytical = concentrations_pb_1d(grid, 0.3, 10e-9)

phi_analytical_dimensional = K_B * T / E_CHARGE * phi_analytical
ca_analytical_dimensional = MOL_PER_LITER_TO_PER_CUBIC_METER * ca_analytical
cc_analytical_dimensional = MOL_PER_LITER_TO_PER_CUBIC_METER * cc_analytical

potential_eq_lhs = second_derivative(phi_analytical_dimensional, 0.00001 * 10e-9)
potential_eq_rhs = -E_CHARGE / EPSILON_VAC / 80 * (cc_analytical_dimensional - ca_analytical_dimensional)

relative_error_potential_eq_residual = abs.(potential_eq_lhs[2:end-1] - potential_eq_rhs[2:end-1]) ./ abs.(potential_eq_lhs[2:end-1])
x = range(0, stop=1, length=length(relative_error_potential_eq_residual))

# Create figure and axis
fig = Figure()
ax = Axis(fig[1, 1]; limits=((0, 0.5), (0, maximum(relative_error_potential_eq_residual[1:floor(Int, length(relative_error_potential_eq_residual) / 2)]))))  # Set x-limits onlyt x-axis limits
lines!(ax, x, relative_error_potential_eq_residual)
# Display the figure
display(fig)
CairoMakie.save(joinpath(".", "verification_simple_pnp_dimensional_residual_check.pdf"), fig)
