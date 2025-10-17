# ======================
#  STAGE 1: Base with R and Conda
# ======================
FROM continuumio/miniconda3:23.3.1-0 AS r_base
LABEL authors="nicholas.sunderland@bristol.ac.uk"

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH=/opt/conda/bin:$PATH

# Install system dependencies
USER root
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    git \
    curl \
    zlib1g-dev \
    libbz2-dev \
    libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Install extra tools
RUN conda install -y -c conda-forge \
    mamba \
    && conda clean --all --yes

# Install AWS CLI
RUN mamba install -y -c conda-forge \
    awscli \
    && conda clean --all --yes

# Preinstall base R to reduce build time for envs
RUN mamba install -y -c conda-forge \
    r-base=4.4.3 \
    r-essentials \
    && conda clean --all --yes

# ===========================
#  STAGE 2: Hermes Container
# ===========================
FROM r_base AS hermes

# Create environments from YAML
COPY environment_hermes.yml /tmp/environment_hermes.yml
RUN mamba env create -f /tmp/environment_hermes.yml && \
    mamba clean --all --yes

COPY environment_gwas2vcf.yml /tmp/environment_gwas2vcf.yml
RUN mamba env create -f /tmp/environment_gwas2vcf.yml && \
    mamba clean --all --yes

# Activate GWAS2VCF env and reinstall vgraph from GitHub
SHELL ["bash", "-c"]
RUN source /opt/conda/etc/profile.d/conda.sh && \
    conda activate gwas2vcf && \
    pip uninstall -y vgraph || true && \
    pip install --no-cache-dir --force-reinstall \
        git+https://github.com/bioinformed/vgraph@5866e6a686ce98e8d02f6f97f6ab025c89f123b9#egg=vgraph

# Get gwas2vcf repository
RUN git clone --branch 1.4.3 https://github.com/MRCIEU/gwas2vcf.git /opt/gwas2vcf
RUN chmod +x /opt/gwas2vcf/main.py

# Copy the QC scripts and entrypoint
COPY run.sh /opt/run.sh
COPY scripts /opt/scripts

WORKDIR /opt
RUN chmod +x run.sh

ENTRYPOINT ["/opt/run.sh"]



# ====================== # Build & Push Commands # ======================
# docker buildx create --use
# push
# docker buildx build --platform linux/amd64 --progress=plain --tag nicksunderland/hermes_docker:latest --file Dockerfile --push . 2>&1 | tee docker_build.log
# or load locally
# docker buildx build --platform linux/amd64 --progress=plain --tag nicksunderland/hermes_docker:latest --file Dockerfile --load . 2>&1 | tee docker_build.log
