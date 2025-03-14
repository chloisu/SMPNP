# 
# Copyright (C) 2025 - 2025 by the authors of the SMPNP code.
# 
# This file is part of SMPNP
# 

"""
This enum is to specify the solver which primary is assigned to what 
index in the solver. This also allows us to quickly switch this order
if necessary.
"""
const PHI_EQ = 1
const ANION_EQ = 2
const CATION_EQ = 3

const NUM_SPECIES = 3
const NUM_CONCENTRATION_SPECIES = 2
const DEFAULT_INITIAL_CONDITIONS = [0.0, 1.0, 1.0]

const SPECIES_VECTOR = Vector{Int}(1:NUM_SPECIES)
const CONCENTRATION_SPECIES_VECTOR = SPECIES_VECTOR[2:end]


# write out 
TIME_WRITE_PRECISION = 5



