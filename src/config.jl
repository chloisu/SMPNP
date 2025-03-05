# 
# Copyright (C) 2025 - 2025 by the authors of the SMPNP code.
# 
# This file is part of SMPNP
# 

# configuration
ALLOWED_PNP_MODES = ["PNP", "MPNP", "SMPNP"]
"""
function to validate the PNP mode is a valid PNP mode.
This function was generated with the help of ChatGPT (OpenAI).
"""
function validate_PNP_mode(mode::String)
    mode in ALLOWED_PNP_MODES || throw(ArgumentError("Invalid PNP_mode: $mode. Allowed values: $(join(ALLOWED_PNP_MODES, ", "))"))
    return mode
end

"""
This enum is to specify the solver which primary is assigned to what 
index in the solver. This also allows us to quickly switch this order
if necessary.
"""
const PHI_EQ = 1
const ANION_EQ = 2
const CATION_EQ = 3
const NUM_SPECIES = 3

const DEFAULT_INITIAL_CONDITIONS = [0.0, 1.0, 1.0]
