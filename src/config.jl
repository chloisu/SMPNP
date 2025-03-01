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
    rlog(u; eps=1.0e-20)

Regularized logarithm. Smooth linear continuation for `x<eps`.
This means we can calculate a "logarithm"  of a small negative number.
"""
function rlog(x; eps=1.0e-20)
    if x < eps
        return log(eps) + (x - eps) / eps
    else
        return log(x)
    end
end
"""
rexp(x;trunc=500.0)

Regularized exponential. Linear continuation for `x>trunc`,  
returns 1/rexp(-x) for `x<-trunc`.
"""
function rexp(x; trunc=500.0)
    if x < -trunc
        1.0 / rexp(-x; trunc)
    elseif x <= trunc
        exp(x)
    else
        exp(trunc) * (x - trunc + 1)
    end
end
