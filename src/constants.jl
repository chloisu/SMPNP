# 
# Copyright (C) 2025 - 2025 by the authors of the SMPNP code.
# 
# This file is part of SMPNP
# 

# Physical constants
const NM = 1e-9 # nanometer conversion factor
const AO = 1e-10 # angstrom conversion factor
const K_B = 1.380649e-23 # Boltzmann constant [J/K]
const N_A = 6.02214076e23 # Avogadro constant
const Z_CATION = 1 # cation charge sign
const Z_ANION = -1 # cation charge sign
const E_CHARGE = 1.6020e-19 # Elementary charge [C]
const EPSILON_VAC = 8.8541878188e-12 # vacuum permittivity [F/m]
# Environtment constants
const T = 293.15 #Temperature [K]
const BETA = 1.0 / K_B / T # [1/J]
# physical prefactor of Poisson equation
const POISSON_PHYS_PREFACTOR = K_B * T / E_CHARGE^2

# conversion constants
MOL_PER_LITER_TO_PER_CUBIC_METER = N_A * 1e3 #from mol/L to m^-3

