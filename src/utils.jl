using Configurations
using YAML

""" 
This function was generated with the help of ChatGPT (OpenAI).
"""
function create_output_directory(params)
    output_dir = params.output_directory

    if isempty(output_dir)
        # Generate a default unique output directory
        base_name = "output_"
        counter = 1
        while isdir(@sprintf("%s%03d", base_name, counter))
            counter += 1
        end
        output_dir = @sprintf("%s%03d", base_name, counter)
    end

    # Create the directory
    mkpath(output_dir)
    return output_dir
end

""" 
This function was generated with the help of ChatGPT (OpenAI).
"""
function save_parameters_to_yaml(params, output_dir)
    output_file = joinpath(output_dir, "parameters.yml")
    d = to_dict(params; include_defaults=true, exclude_nothing=false)
    YAML.write_file(output_file, d)
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


""" 
Custom Error messages
"""
struct InputError <: Exception
    msg::String
end

function Base.showerror(io::IO, e::InputError)
    print(io, "InputError: ", e.msg)
end

struct ParameterError <: Exception
    msg::String
end

function Base.showerror(io::IO, e::ParameterError)
    print(io, "ParameterError: ", e.msg)
end