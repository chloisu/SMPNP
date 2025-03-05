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
    C_REF::Float64 = 1.0 # reference concentration [mol/liter]
    D_REF::Float64 = 5e-9 # reference diffusivity [meter^2/second]
    PHI_REF::Float64 = K_B * T / E_CHARGE # reference potential [V]
end
