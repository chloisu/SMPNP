configfile: "snakemake_config.yml"

from pathlib import Path


raw_data = config["RAW_DATA"]
# Get a list of folder names in 'raw_data/' directory
cases = [p.name for p in Path(raw_data).iterdir() if p.is_dir()]

postprocessing_folder = config["POSTPROCESSING_DATA"]

rule all:
    input:
       expand("{folder}/{case}/with_additional_fields.vtu",
       folder=postprocessing_folder,
       case=cases) + 
       expand("{folder}/{case}/extract_pore_center_radius_line.csv",
       folder=postprocessing_folder,
       case=cases) +
       expand("{folder}/{case}/extract_inlet_line.csv",
       folder=postprocessing_folder,
       case=cases)

rule compute_additional_flux_fields:
    input:
        yml=lambda wildcards: f"{raw_data}/{wildcards.case}/parameters.yml",
        pvd=lambda wildcards: f"{raw_data}/{wildcards.case}/out.pvd"
    output:
        vtu="{postprocessing_folder}/{case}/with_additional_fields.vtu"
    shell:
        """
        export PATH=/Applications/ParaView-5.13.2.app/Contents/bin:$PATH
        pvpython --venv /Users/nik/.pyenv/versions/paraview_3.10_env ./scripts/postprocessing/paraview_compute_fluxes.py {input.yml} {input.pvd} {output.vtu}
        """

rule extract_pore_center_radius_line:
    input:
        yml=lambda wildcards: f"{raw_data}/{wildcards.case}/parameters.yml",
        vtu=lambda wildcards: f"{postprocessing_folder}/{wildcards.case}/with_additional_fields.vtu"
    output:
        csv="{postprocessing_folder}/{case}/extract_pore_center_radius_line.csv"
    shell:
        """
        export PATH=/Applications/ParaView-5.13.2.app/Contents/bin:$PATH
        pvpython --venv /Users/nik/.pyenv/versions/paraview_3.10_env ./scripts/postprocessing/paraview_extract_pore_center_radius_line.py {input.yml} {input.vtu} {output.csv}
        """

rule extract_inlet_line:
    input:
        yml=lambda wildcards: f"{raw_data}/{wildcards.case}/parameters.yml",
        vtu=lambda wildcards: f"{postprocessing_folder}/{wildcards.case}/with_additional_fields.vtu"
    output:
        csv="{postprocessing_folder}/{case}/extract_inlet_line.csv"
    shell:
        """
        export PATH=/Applications/ParaView-5.13.2.app/Contents/bin:$PATH
        pvpython --venv /Users/nik/.pyenv/versions/paraview_3.10_env ./scripts/postprocessing/paraview_extract_inlet_line.py {input.yml} {input.vtu} {output.csv}
        """
