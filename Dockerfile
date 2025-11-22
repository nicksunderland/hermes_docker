# ======================
#  STAGE 1: Base with R and Conda
# ======================
FROM continuumio/miniconda3:23.3.1-0 AS r_base
LABEL authors="nicholas.sunderland@bristol.ac.uk"

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH=/opt/conda/bin:$PATH
ENV PYTHONUNBUFFERED=1

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
        git+https://github.com/bioinformed/vgraph@5866e6a686ce98e8d02f6f97f6ab025c89f123b9#egg=vgraph && \
    conda deactivate

# Activate hermes env and install my fixed version of ldscr
RUN source /opt/conda/etc/profile.d/conda.sh && \
    conda activate hermes && \
    git clone https://github.com/nicksunderland/ldscr.git /tmp/ldscr && \
    Rscript -e "install.packages(c('checkmate','corrr','gdata','ggtext'), repos='https://cloud.r-project.org')" && \
    Rscript -e "install.packages('/tmp/ldscr', repos = NULL, type = 'source')" && \
    rm -rf /tmp/ldscr && \
    conda deactivate

# Get gwas2vcf repository
RUN git clone --branch 1.4.3 https://github.com/MRCIEU/gwas2vcf.git /gwas2vcf && \
    chmod +x /gwas2vcf/main.py

# Patch vcf.py to add null-allele check
RUN python3 - <<EOF
fname = "/gwas2vcf/vcf.py"
lines = []

with open(fname) as f:
    for line in f:
        lines.append(line)
        if "record.id = Vcf.remove_illegal_chars(result.dbsnpid)" in line:
            # Detect leading whitespace of the line
            indent = line[:len(line) - len(line.lstrip())]
            inner_indent = indent + "    "  # 1 level deeper for nested block

            # Append patch using the detected indentation
            lines.append(f"{indent}# NS added - for some reason there can be null alleles\n")
            lines.append(f"{indent}if not result.ref or not result.alt:\n")
            lines.append(f"{inner_indent}logging.warning(f\"Skipping variant with null allele: chrom={{result.chrom}}, pos={{result.pos}}, ref={{result.ref}}, alt={{result.alt}}, dbsnpid={{result.dbsnpid}}\")\n")
            lines.append(f"{inner_indent}continue\n")

with open(fname, "w") as f:
    f.writelines(lines)
EOF

# check patch compiles
RUN python3 -m py_compile /gwas2vcf/vcf.py

# Copy the QC scripts and entrypoint
COPY run.sh /run.sh
COPY scripts /scripts

RUN chmod +x /run.sh

ENTRYPOINT ["/run.sh"]



# ====================== # Build & Push Commands # ======================
# docker buildx create --use

# all in one:
# docker buildx build --platform linux/amd64 --progress=plain --tag nicksunderland/hermes_docker:latest --file Dockerfile --cache-from type=registry,ref=nicksunderland/hermes_docker:latest --cache-to type=inline --push . 2>&1 | tee docker_build.log

# clean local and pull:
# docker rmi nicksunderland/hermes_docker:latest
# docker image prune -f
# docker pull --platform linux/amd64 nicksunderland/hermes_docker:latest