# 
# Copyright (C) 2025 - 2025 by the authors of the SMPNP code.
# 
# This file is part of SMPNP
# 

# in this file, we implement the post processing routines that are configured
# using an input .yml file similar to the simulation runs. The post processing 
# options are quite limited at this point but could evolve over apply_dirichlet_for_time_dependent

# as a given initial state, we assume that the simulation run preceeding the postrpocesing
# step has produced outpufiles of the patter *.out where * is the time in non-dimensional diffusion
# times.


@option struct ExtractSlice
    case::String = ""# name of the line
    slice::Any = "" #specify the slice (required)
end

@option struct Integration
    case::String = ""# name of the integration
    boundary_id::Int = 0    # boundary_id for which we want to know the flux.
end

@option struct QuantityOverTime
    case::String

end

"""
this is the main function to postprocess the solution. It get's invoked if the parameter
file provided is a postprocessing parameter file or if the command line argument -postprocess
is provided.
"""
function smpnp_postprocess()
    # we do some checks on the solution 
    # we do all the extract lines parts
    # we do all the integration parts
    # we do all the quantitiy over time parts
end
