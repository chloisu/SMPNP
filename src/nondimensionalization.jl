# 
# Copyright (C) 2025 - 2025 by the authors of the SMPNP code.
# 
# This file is part of SMPNP
# 

"""
Here we define the struct of NonDimensionalization variables to nondimensionalize the system.
"""
@option struct NonDimensionalization
    L_REF::Float64 = 50e-9 # reference length [meter]
    C_REF::Float64 = MOL_PER_LITER_TO_PER_CUBIC_METER # reference concentration [1/m^3]
    D_REF::Float64 = 5e-9 # reference diffusivity [meter^2/second]
    PHI_REF::Float64 = K_B * T / E_CHARGE # reference potential [V]
    EPSILON_REF::Float64 = 80 * EPSILON_VAC
end

function validate_non_dim(non_dim_params)
    non_dim_params.L_REF > 0.0 || throw(AssertionError("The reference length L_REF needs to be chosen larger than 0.0"))
    non_dim_params.C_REF > 0.0 || throw(AssertionError("The reference concentration C_REF needs to be chosen larger than 0.0"))
    non_dim_params.D_REF > 0.0 || throw(AssertionError("The reference diffusion coefficient D_REF needs to be chosen larger than 0.0"))
    non_dim_params.PHI_REF > 0.0 || throw(AssertionError("The reference potential PHI_REF needs to be chosen larger than 0.0"))
    non_dim_params.EPSILON_REF > 0.0 || throw(AssertionError("The reference permittivity EPSILON_REF needs to be chosen larger than 0.0"))
end