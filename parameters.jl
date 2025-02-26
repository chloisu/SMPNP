
using Configurations

@option struct Parameters
    case::String = "test_case" # name of the case 
    # the type of the PNP system. There is three different options to choose from
    # PNP : the classical Poisson-Nernst-Planck system without taking the size of the ions into account.
    # MPNP: the size of the ions is taken into account but with a single constant size parameters (i.e. both ions are assumed to have the same size.)
    # SMPNP: each of the two ions has a differnt size
    PNP_mode::String = "PNP" |> validate_PNP_mode
    D_ref::Float64 = 5e-9 # reference diffusivity [m^2/s]
    pore_radius::Float64 = 5e-9 # pore radius [meter]
    pore_length::Float64 = 10e-9 # pore length [meter]
    cation_size::Float64 = 5e-10 # cation size [meter]
    anion_size::Float64 = 5e-10 # cation size [meter]
    solvent_molecule_size::Float64 = 5e-10
    D_anion::Float64 = 5e-9 # anion diffusivity [m^2/s]
    D_cation::Float64 = 5e-9 # cation diffusivity [m^2/s]
    surface_potential::Float64 = 0.1 # surface potential [V]
    #c_0::Float64 = 1 # initial concentration [mol/L]
    #max_iter::Int = 20 # maximum nonlinear Newton iterations 
    #dt::Float64 = 1e-6 # timestep
    #dt_max::Float64 = 2e-1 # maximum timestep
    #t_final::Float64 = 1 # final time
end