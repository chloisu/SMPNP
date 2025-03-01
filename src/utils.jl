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

