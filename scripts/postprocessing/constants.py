PHI_EQ = 1
ANION_EQ = 2
CATION_EQ = 3

NM = 1e-9 # nanometer conversion factor
AO = 1e-10 # angstrom conversion factor
K_B = 1.380649e-23 # Boltzmann constant [J/K]
N_A = 6.02214076e23 # Avogadro constant
Z_CATION = 1 # cation charge sign
Z_ANION = -1 # cation charge sign
E_CHARGE = 1.6020e-19 # Elementary charge [C]
EPSILON_VAC = 8.8541878188e-12 # vacuum permittivity [F/m]
# Environtment constants
T = 293.15 #Temperature [K]
BETA = 1.0 / K_B / T # [1/J]
# physical prefactor of Poisson equation
POISSON_PHYS_PREFACTOR = K_B * T / E_CHARGE**2

# conversion constants
MOL_PER_LITER_TO_PER_CUBIC_METER = N_A * 1e3 #from mol/L to m^-3