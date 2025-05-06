configfile: "snakemake_config.yml"

cases = ["20250506_r_10nm_h_20nm_l_50nm_phiT_50"]
raw_data = config["RAW_DATA"]
postprocessing_folder = config["POSTPROCESSING_DATA"]

rule all:
    input:
       expand("{folder}/{case}/with_additional_fields.vtu",
       folder=postprocessing_folder,
       case=cases)

rule process_case:
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
