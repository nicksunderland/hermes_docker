#!/usr/bin/env bash
set -euo pipefail

show_help() {
  cat << EOF
Usage: $(basename "$0") <data-path> <json-path> <resources-dir> <file-guid>

Description:
  This script runs the Hermes QC job using the specified input data and configuration.

Arguments:
  <data-path>        Path to the input data (e.g., ~/Desktop/data.tsv, or, s3://bucket/path/data.csv)
  <json-path>        Path to the JSON config file
  <resources-dir>    Local or S3 path to the resources directory
  <file-guid>        Unique output file GUID or identifier

Options:
  -h, --help         Show this help message and exit

Examples:
  $(basename "$0") s3://my-bucket/input.tsv.gz config.json ./resources abc123
  $(basename "$0") ~/Desktop/data.tsv config.json ./resources abc123
EOF
}

# Check if user asked for help
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  show_help
  exit 0
fi

# Ensure 4 arguments are provided
if [[ $# -lt 4 ]]; then
  echo "Error: Missing required arguments."
  echo "Use --help for usage information."
  exit 1
fi

# Arguments passed to the script
DATA_PATH="$1"
JSON_PATH="$2"
RESOURCES_DIR="$3"
FILE_GUID="$4"

echo "=== RUN INPUT ==="
echo "GWAS file path:       $DATA_PATH"
echo "JSON file path:       $JSON_PATH"
echo "Resources directory:  $RESOURCES_DIR"
echo "File output / guid:   $FILE_GUID"
echo "================="


#########################################
# UTILITY FUNCTIONS
#########################################
activate_conda() {
    local env="$1"
    set +u
    source /opt/conda/etc/profile.d/conda.sh
    conda activate "$env"
    set -u
}
download_or_copy() {
    local src="$1"
    local dest="$2"
    if [[ "$src" == s3://* ]]; then
        # S3 - Must download to /tmp
        aws s3 cp "$src" "$dest" --no-progress
    else
        cp "$src" "$dest"
    fi
}
setup_resource() {
    local file_name="$1"
    if [[ "$RESOURCES_DIR" == s3://* ]]; then
        local local_dest="/tmp/resources/$file_name"
        mkdir -p "/tmp/resources"
        download_or_copy "$RESOURCES_DIR/$file_name" "$local_dest"
        echo "$local_dest"
    else
        echo "${RESOURCES_DIR}/${file_name}"
    fi
}
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}
validate_output() {
    local file="$1"
    local step="$2"
    if [[ ! -s "$file" ]]; then 
        log "CRITICAL ERROR: ${step} failed to produce output."
        echo "Missing or empty file: $file"
        exit 1
    else
        log "SUCCESS: ${step} verified output: $file"
    fi
}
ensure_unzipped() {
    local gz_path="$1"     # The path to the .gz file (e.g., /resources/ref.fasta.gz)
    local unz_path="$2"    # The path where we expect the unzipped file (e.g., /resources/ref.fasta)
    # 1. If the unzipped file already exists at the source (e.g. mounted), use it
    if [[ -f "$unz_path" ]]; then
        echo "$unz_path"
        return
    fi
    # 2. If unzipped doesn't exist, we must unzip the GZ file to /tmp/scratch
    local scratch_dest
    scratch_dest="/tmp/$(basename "$unz_path")"
    # Check if we already unzipped it in a previous run/step
    if [[ -f "$scratch_dest" ]]; then
        echo "$scratch_dest"
        return
    fi
    log "Unzipping $(basename "$gz_path") to /tmp..." >&2
    # Decompress to stdout > destination to handle read-only source mounts
    gzip -dc "$gz_path" > "$scratch_dest"
    echo "$scratch_dest"
}

#########################################
# SETUP FILES
#########################################
LOCAL_GWAS_FILE="/tmp/input_gwas.tsv.gz"
LOCAL_OUTPUT="/tmp/output"
mkdir -p "$LOCAL_OUTPUT"

# -------------------------------------------------------
# 1. DEFINE RESOURCE PATHS
# -------------------------------------------------------
# Reference and resource files
EAF_PATH=$(setup_resource "all_afreq_b38.tsv.gz")
FASTA_38_PATH=$(setup_resource "Homo_sapiens_assembly38_nochr.fasta.gz")
FASTA_38_UNZ_TARGET=$(setup_resource "Homo_sapiens_assembly38_nochr.fasta") # might not exist yet
DBSNP_PATH=$(setup_resource "dbSNP157_b38_clean.vcf.gz")

# -------------------------------------------------------
# 2. HANDLE GENOME BUILD
# -------------------------------------------------------
log "Getting GWAS build"
BUILD=$(python -c "
import json, os

json_input = '$JSON_PATH'

try:
    # Try reading as a file (local dev case)
    with open(json_input) as f:
        data = json.load(f)
except (FileNotFoundError, OSError, IsADirectoryError, json.JSONDecodeError):
    # Otherwise, assume it's a raw JSON string (server case)
    data = json.loads(json_input)

print(data.get('referenceGenome', ''))
")

# Check we know the build
if [[ -z "$BUILD" ]]; then
    echo 'ERROR: Could not extract referenceGenome from JSON_PATH'
    exit 1
else
    log "Detected build: $BUILD"
fi

# update paths based on build
if [[ "$BUILD" == "Hg19" ]]; then
  CHAIN_PATH=$(setup_resource "hg19ToHg38.over.chain.gz")
  FASTA_37_PATH=$(setup_resource "human_g1k_v37.fasta.gz")
  FASTA_37_UNZ_TARGET=$(setup_resource "human_g1k_v37.fasta") # might not exist yet
fi

# -------------------------------------------------------
# 3. GET INPUT DATA
# -------------------------------------------------------
log "Getting GWAS input"
download_or_copy "$DATA_PATH" "$LOCAL_GWAS_FILE"


#########################################
# STEP 1: QC PARSE
#########################################
activate_conda hermes
log "Running QC Step 1"
OUTPUT_PARSED_GWAS="${LOCAL_OUTPUT}/parsed_gwas.tsv"
OUTPUT_GWAS2VCF_JSON="${LOCAL_OUTPUT}/gwas2vcf.json"
OUTPUT_STEP1_SUMMARY="${LOCAL_OUTPUT}/step1_summary.tsv"

Rscript /scripts/01_parse_data.R \
  "$LOCAL_GWAS_FILE" \
  "$JSON_PATH" \
  "$OUTPUT_PARSED_GWAS" \
  "$OUTPUT_GWAS2VCF_JSON" \
  "$OUTPUT_STEP1_SUMMARY"
  
validate_output "$OUTPUT_PARSED_GWAS" "Step 1 (QC Parse)"
validate_output "$OUTPUT_GWAS2VCF_JSON" "Step 1 (JSON Generation)"


#########################################
# STEP 2: gwas2vcf
#########################################
activate_conda gwas2vcf
log "Running Step 2 gwas2vcf"
OUTPUT_VCF="${LOCAL_OUTPUT}/gwas.vcf.gz"
OUTPUT_STEP2_SUMMARY="${LOCAL_OUTPUT}/step2_summary.tsv"

REF=$(ensure_unzipped "$FASTA_38_PATH" "$FASTA_38_UNZ_TARGET")
if [[ "$BUILD" == "Hg19" ]]; then
    REF=$(ensure_unzipped "$FASTA_37_PATH" "$FASTA_37_UNZ_TARGET")
fi

python /gwas2vcf/main.py \
    --data "$OUTPUT_PARSED_GWAS" \
    --json "$OUTPUT_GWAS2VCF_JSON" \
    --id "gwas_qc_pipeline" \
    --ref "$REF" \
    --out "$OUTPUT_VCF" \
    2>&1 | tee "$OUTPUT_STEP2_SUMMARY"
    
validate_output "$OUTPUT_VCF" "Step 2 (gwas2vcf conversion)"


#########################################
# STEP 3: Liftover & annotate
#########################################
activate_conda hermes
log "Running Step 3 bcftools liftover & annotate"
DBSNP="$DBSNP_PATH"
CORES=$(nproc)
OUTPUT_VCF_B38="${LOCAL_OUTPUT}/gwas_b38.vcf.gz"
OUTPUT_STEP3_SUMMARY="${LOCAL_OUTPUT}/step3_summary.tsv"

if [[ "$BUILD" == "Hg19" ]]; then

    REF_37=$(ensure_unzipped "$FASTA_37_PATH" "$FASTA_37_UNZ_TARGET")
    REF_38=$(ensure_unzipped "$FASTA_38_PATH" "$FASTA_38_UNZ_TARGET")

    bcftools +liftover "$OUTPUT_VCF" \
        -- -s "$REF_37" \
        -f "$REF_38" \
        -c "$CHAIN_PATH" \
        2>> "$OUTPUT_STEP3_SUMMARY" \
        | bcftools view -Oz -o "$OUTPUT_VCF_B38.tmp.gz" -

    bcftools sort -Oz -o "$OUTPUT_VCF_B38.tmp.sorted.gz" "$OUTPUT_VCF_B38.tmp.gz"
    mv "$OUTPUT_VCF_B38.tmp.sorted.gz" "$OUTPUT_VCF_B38.tmp.gz"
else
    cp "$OUTPUT_VCF" "$OUTPUT_VCF_B38.tmp.gz"
    : > "$OUTPUT_STEP3_SUMMARY"
fi

bcftools index -f --csi "$OUTPUT_VCF_B38.tmp.gz"

# Annotate with dbSNP
bcftools annotate -a "$DBSNP" -c ID -O z --threads "$CORES" \
    "$OUTPUT_VCF_B38.tmp.gz" -o "$OUTPUT_VCF_B38"

# Cleanup
rm -f "$OUTPUT_VCF_B38.tmp.gz" "$OUTPUT_VCF_B38.tmp.gz.csi"

bcftools index -f --csi "$OUTPUT_VCF_B38"

validate_output "$OUTPUT_VCF_B38" "Step 3 (Liftover & Annotation)"


#########################################
# STEP 4: Reporting
#########################################
log "Running QC Step 4 reporting"
EAF_REF="$EAF_PATH"
OUTPUT_CLEAN_GWAS="${LOCAL_OUTPUT}/gwas_b38_clean.tsv.gz"
OUTPUT_REPORT="${LOCAL_OUTPUT}/gwas_report.html"

Rscript /scripts/04_qc_report.R \
    "$OUTPUT_VCF_B38" \
    "$JSON_PATH" \
    "$OUTPUT_STEP1_SUMMARY" \
    "$OUTPUT_STEP2_SUMMARY" \
    "$OUTPUT_STEP3_SUMMARY" \
    "$EAF_REF" \
    "/scripts/04_qc_report.Rmd" \
    "$OUTPUT_CLEAN_GWAS" \
    "$OUTPUT_REPORT"
    
validate_output "$OUTPUT_CLEAN_GWAS" "Step 4 (Cleaned GWAS TSV)"
validate_output "$OUTPUT_REPORT" "Step 4 (HTML Report)"


#########################################
# STEP 5: Upload QC outputs to S3
#########################################
log "Uploading QC outputs"

S3_BUCKET="hermes-qc"

# Ensure outputs exist
if [[ ! -f "$OUTPUT_CLEAN_GWAS" ]]; then
    echo "ERROR: file not found: $OUTPUT_CLEAN_GWAS"
    exit 1
fi

# If FILE_GUID looks like a directory, use local mode; otherwise assume S3 upload
if [[ "$FILE_GUID" == /* || "$FILE_GUID" == ./* || -d "$FILE_GUID" ]]; then
    log "Detected local output path — copying locally"
    mkdir -p "$FILE_GUID/images" "$FILE_GUID/tables"

    cp "$OUTPUT_CLEAN_GWAS" "$FILE_GUID/tables/"
    for file in "$LOCAL_OUTPUT/gwas_report.html" "$LOCAL_OUTPUT/eaf_plot.png" "$LOCAL_OUTPUT/pz_plot.png" "$LOCAL_OUTPUT/qq_plot.png"; do
        if [[ ! -f "$file" ]]; then
            echo "ERROR: file not found: $file"
            exit 1
        fi
        cp "$file" "$FILE_GUID/images/"
    done
else
    log "Detected S3 output mode — uploading to hermes-qc bucket under GUID ${FILE_GUID}"

    # Upload clean GWAS TSV to tables/
    aws s3 cp "$OUTPUT_CLEAN_GWAS" \
        "s3://${S3_BUCKET}/tables/${FILE_GUID}/$(basename "$OUTPUT_CLEAN_GWAS")" \
        --content-type "text/tab-separated-values" \
        --no-progress

    # Upload HTML + PNG files to images/
    for file in "$LOCAL_OUTPUT/gwas_report.html" "$LOCAL_OUTPUT/eaf_plot.png" "$LOCAL_OUTPUT/pz_plot.png" "$LOCAL_OUTPUT/qq_plot.png"; do
        if [[ ! -f "$file" ]]; then
            echo "ERROR: file not found: $file"
            exit 1
        fi
        content_type="image/png"
        [[ "$file" == *.html ]] && content_type="text/html"

        aws s3 cp "$file" \
            "s3://${S3_BUCKET}/images/${FILE_GUID}/$(basename "$file")" \
            --content-type "$content_type" \
            --no-progress
    done
fi

log "QC pipeline completed successfully"