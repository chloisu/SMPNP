import os
from pathlib import Path

configfile: "papers/2026/Overlapping-Electrical-Double-Layer-in-Nanoporous-Carbon-Separator-Coating-Enhances-Rate-Capability-of-Lithium-Ion-Batteries/snakemake_config.yml"

# --- CONFIGURATION AND SETUP ---
raw_data              = config["RAW_DATA"]
postprocessing_folder = config["POSTPROCESSING_DATA"]
scripts_dir           = config["SCRIPTS_DIR"]

# --- PORTABLE pvpython HANDLE ---
# Resolution order: explicit config value  ->  environment variable  ->  bare
# `pvpython` found on PATH. The container image exports PVPYTHON and
# PVPYTHON_ARGS, so nothing ParaView-specific needs to be committed here.
#
# On a local macOS machine, either export the variables before running, e.g.
#   export PVPYTHON="/Applications/ParaView-5.13.2.app/Contents/bin/pvpython"
#   export PVPYTHON_ARGS="--venv /Users/nik/.pyenv/versions/paraview_3.10_env"
# or add PVPYTHON / PVPYTHON_ARGS to snakemake_config.yml.
PVPYTHON      = config.get("PVPYTHON",      os.environ.get("PVPYTHON", "pvpython"))
PVPYTHON_ARGS = config.get("PVPYTHON_ARGS", os.environ.get("PVPYTHON_ARGS", ""))

# --- WILDCARDS ---
# Get a list of folder names in the 'raw_data/' directory.
cases = [p.name for p in Path(raw_data).iterdir() if p.is_dir()]

# Define the different extraction types to avoid repetition.
EXTRACTIONS = [
    "pore_center_radius_line",
    "inlet_line",
    "along_x_center_line",
]

# --- RULES ---
rule all:
    input:
        expand("{folder}/{case}/with_additional_fields.vtu",
               folder=postprocessing_folder,
               case=cases),
        expand("{folder}/{case}/extract_{extraction}.csv",
               folder=postprocessing_folder,
               case=cases,
               extraction=EXTRACTIONS)

rule compute_additional_flux_fields:
    input:
        yml=f"{raw_data}/{{case}}/parameters.yml",
        pvd=f"{raw_data}/{{case}}/out.pvd",
        script=f"{scripts_dir}/postprocessing/paraview_compute_fluxes.py"
    output:
        vtu="{postprocessing_folder}/{case}/with_additional_fields.vtu"
    params:
        pvpython=PVPYTHON,
        pvargs=PVPYTHON_ARGS,
    shell:
        # Portable: pvpython is resolved from the environment/config, and runs
        # headlessly via PVPYTHON_ARGS (--force-offscreen-rendering) inside the
        # container. No hardcoded application paths.
        "{params.pvpython} {params.pvargs} {input.script} {input.yml} {input.pvd} {output.vtu}"

rule extract_line:
    input:
        yml=f"{raw_data}/{{case}}/parameters.yml",
        vtu="{postprocessing_folder}/{case}/with_additional_fields.vtu",
        script=lambda wildcards: f"{scripts_dir}/postprocessing/paraview_extract_{wildcards.extraction}.py"
    output:
        csv="{postprocessing_folder}/{case}/extract_{extraction}.csv"
    params:
        pvpython=PVPYTHON,
        pvargs=PVPYTHON_ARGS,
    shell:
        "{params.pvpython} {params.pvargs} {input.script} {input.yml} {input.vtu} {output.csv}"
