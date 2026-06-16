# syntax=docker/dockerfile:1.6
#
# Container for a mixed Julia + Python (Snakemake) codebase.
# Works directly with Docker, and can be consumed by Apptainer via either:
#   apptainer build app.sif docker-daemon://<image>:<tag>     (after `docker build`)
#   apptainer build app.sif docker://<registry>/<image>:<tag> (after `docker push`)
# Or build natively with Apptainer using the accompanying apptainer.def.
#
# Multi-arch: builds natively on linux/amd64 and linux/arm64. On Apple Silicon
# this means a fast native arm64 build. For x86-based HPC deployment, pass
# `--platform=linux/amd64` explicitly.

FROM ubuntu:22.04

ARG JULIA_VERSION=1.11.3
ARG MINIFORGE_VERSION=24.11.3-0
ARG PYTHON_VERSION=3.11

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    TZ=Etc/UTC

# ---------------------------------------------------------------------------
# System packages
#   - OpenGL / X11 / xkbcommon — needed by GLMakie and Gmsh GUI.
#   - Xvfb — used during the Julia precompile so GLMakie can initialise a
#     virtual GL context without a real display.
#   - libgomp — OpenMP runtime pulled in by several SciML / sparse JLLs.
# Note: libquadmath0 is intentionally NOT installed; it's x86-only and the
# Julia JLLs bundle their own compiler-support libraries.
# ---------------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates curl wget bzip2 xz-utils tar unzip git \
        build-essential \
        locales tzdata \
        xvfb \
        libgl1 libglu1-mesa libglx-mesa0 libegl1 libgles2 \
        libglvnd0 libglx0 libopengl0 \
        libxrandr2 libxinerama1 libxcursor1 libxi6 libxext6 \
        libxrender1 libxfixes3 libxft2 libxt6 libxkbcommon0 libxkbcommon-x11-0 \
        libxcb1 libxcb-glx0 libxcb-render0 libxcb-shm0 libxcb-xkb1 \
        libsm6 libice6 \
        libgomp1 libfontconfig1 libfreetype6 \
        libpng16-16 libjpeg-turbo8 libtiff5 \
 && locale-gen en_US.UTF-8 \
 && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# Julia (matches the version recorded in Manifest.toml). Multi-arch download.
# ---------------------------------------------------------------------------
ENV JULIA_PATH=/opt/julia
ENV PATH=${JULIA_PATH}/bin:${PATH}

RUN ARCH=$(uname -m) \
 && case "$ARCH" in \
        x86_64)  JL_PATH_ARCH=x64;     JL_FILE_ARCH=x86_64  ;; \
        aarch64) JL_PATH_ARCH=aarch64; JL_FILE_ARCH=aarch64 ;; \
        *) echo "Unsupported architecture: $ARCH" >&2; exit 1 ;; \
    esac \
 && JULIA_MM="${JULIA_VERSION%.*}" \
 && mkdir -p ${JULIA_PATH} \
 && curl -fsSL "https://julialang-s3.julialang.org/bin/linux/${JL_PATH_ARCH}/${JULIA_MM}/julia-${JULIA_VERSION}-linux-${JL_FILE_ARCH}.tar.gz" \
      | tar -xz -C ${JULIA_PATH} --strip-components=1 \
 && julia --version

# ---------------------------------------------------------------------------
# Miniforge (conda + mamba). Multi-arch installer.
# ---------------------------------------------------------------------------
ENV CONDA_DIR=/opt/conda
ENV PATH=${CONDA_DIR}/bin:${PATH}

RUN ARCH=$(uname -m) \
 && curl -fsSL -o /tmp/miniforge.sh \
        "https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}/Miniforge3-${MINIFORGE_VERSION}-Linux-${ARCH}.sh" \
 && bash /tmp/miniforge.sh -b -p ${CONDA_DIR} \
 && rm /tmp/miniforge.sh \
 && conda config --system --set always_yes yes \
 && conda config --system --set auto_update_conda false \
 && conda clean -afy

# ---------------------------------------------------------------------------
# Python environment (conda env "app", pip-installed from requirements.txt)
# ---------------------------------------------------------------------------
COPY requirements.txt /tmp/requirements.txt

RUN mamba create -n app -y -c conda-forge python=${PYTHON_VERSION} pip \
 && /opt/conda/envs/app/bin/pip install --no-cache-dir -r /tmp/requirements.txt \
 && conda clean -afy

# Put the "app" env first on PATH so `python` resolves to it by default.
ENV PATH=/opt/conda/envs/app/bin:${CONDA_DIR}/bin:${JULIA_PATH}/bin:$PATH \
    CONDA_DEFAULT_ENV=app

# ---------------------------------------------------------------------------
# Julia project: copy Project.toml + Manifest.toml and instantiate.
# Precompile is best-effort: Metal.jl is macOS-only and will fail on Linux;
# Pkg.precompile continues past per-package errors, which we swallow here.
# Xvfb gives GLMakie a virtual display for its precompile-time GL context.
# ---------------------------------------------------------------------------
ENV JULIA_DEPOT_PATH=/opt/julia-depot \
    JULIA_PROJECT=/opt/project \
    JULIA_NUM_THREADS=auto

RUN mkdir -p ${JULIA_DEPOT_PATH} ${JULIA_PROJECT}
COPY Project.toml Manifest.toml ${JULIA_PROJECT}/

RUN cd ${JULIA_PROJECT} \
 && xvfb-run -a julia --project=. -e ' \
        using Pkg; \
        Pkg.instantiate(); \
        try \
            Pkg.precompile(); \
        catch e \
            @warn "Pkg.precompile finished with errors (expected for Metal.jl on Linux)" exception=(e, catch_backtrace()); \
        end' \
 && chmod -R a+rX ${JULIA_DEPOT_PATH} ${JULIA_PROJECT}

# ---------------------------------------------------------------------------
# Runtime depot path: prepend the user's writable ~/.julia so that the
# pre-instantiated /opt/julia-depot stays effectively read-only in Apptainer
# (where the container runs as the host user, not root).
# ---------------------------------------------------------------------------
ENV JULIA_DEPOT_PATH=":/opt/julia-depot"

WORKDIR /work

CMD ["/bin/bash", "-l"]
