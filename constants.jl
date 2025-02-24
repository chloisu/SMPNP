
# Physical constants
const NM = 1e-9 # nanometer conversion factor
const AO = 1e-10 # angstrom conversion factor
const K_B = 1.380649e-23 # Boltzmann constant [J/K]
const R = 8.31446261815324 # Molar gas constant [J/(K⋅mol)]
const N_A = 6.02214076e23 # Avogadro constant
const F = 9.64853321233100184e4 # Faraday constant [C/mol]
const Z_CATION = 1 # cation charge sign
const Z_ANION = -1 # cation charge sign
const E_CHARGE = 1.6020e-19 # Elementary charge [C]
const EPSILON_VAC = 8.8541878188e-12 # vacuum permittivity [F/m]
# Material constants
const EPSILON_R = 80 # Relative dielectric constant

# Environtment constants
const T = 293.15 #Temperature [K]
const L_B = 7.5e-10 # Bjerrum length [m]
const BETA = 1.0 / K_B / T # [1/J]
# conversion constants
MOL_PER_LITER_TO_PER_CUBIC_METER = N_A * 1e3 #from mol/L to m^-3

# double layer size for 1 Molar
const LAMBDA_D = sqrt(EPSILON_R * EPSILON_VAC * K_B * T / 2 / E_CHARGE^2 / MOL_PER_LITER_TO_PER_CUBIC_METER)
