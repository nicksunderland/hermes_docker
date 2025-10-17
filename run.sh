#!/usr/bin/env bash
set -euo pipefail

# Arguments passed to the script
S3_DATA_PATH="$1"
JSON_PATH="$2"
RESOURCES_DIR="$3"
FILE_GUID="$4"

echo "=== RUN INPUT ==="
echo "S3_PATH file path:    $S3_DATA_PATH"
echo "JSON file path:       $JSON_PATH"
echo "Resources directory:  $RESOURCES_DIR"
echo "File guid:            $FILE_GUID"
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
        aws s3 cp "$src" "$dest" --no-progress
    else
        cp "$src" "$dest"
    fi
}
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

#########################################
# SETUP FILES
#########################################
LOCAL_GWAS_FILE="/tmp/input_gwas.tsv.gz"
LOCAL_RESOURCES="/tmp/resources"
LOCAL_OUTPUT="/tmp/output"
mkdir -p "$LOCAL_RESOURCES" "$LOCAL_OUTPUT"

# Load genome build from JSON file (or string)
BUILD=$(python -c "
import json, os

json_input = os.environ.get('JSON_PATH')

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
fi

# Reference and resource files
EAF="all_afreq_b38.tsv.gz"
FASTA_37="human_g1k_v37.fasta.gz"
FASTA_37_UNZ="human_g1k_v37.fasta"
FASTA_38="Homo_sapiens_assembly38_nochr.fasta.gz"
FASTA_38_UNZ="Homo_sapiens_assembly38_nochr.fasta"
CHAIN_37_38="hg19ToHg38.over.chain.gz"
DBSNP_38="dbSNP157_b38_clean.vcf.gz"

# Download or copy input files
log "Getting GWAS input"
download_or_copy "$S3_DATA_PATH" "$LOCAL_GWAS_FILE"

log "Getting resource files"
download_or_copy "$RESOURCES_DIR/$EAF" "$LOCAL_RESOURCES/$EAF"
download_or_copy "$RESOURCES_DIR/$FASTA_38" "$LOCAL_RESOURCES/$FASTA_38"
download_or_copy "$RESOURCES_DIR/$FASTA_38.fai" "$LOCAL_RESOURCES/${FASTA_38}.fai"
download_or_copy "$RESOURCES_DIR/$DBSNP_38" "$LOCAL_RESOURCES/$DBSNP_38"
download_or_copy "$RESOURCES_DIR/$DBSNP_38.tbi" "$LOCAL_RESOURCES/${DBSNP_38}.tbi"

if [[ "$BUILD" == "Hg19" ]]; then
    log "Getting Hg19 resources for liftover"
    download_or_copy "$RESOURCES_DIR/$CHAIN_37_38" "$LOCAL_RESOURCES/$CHAIN_37_38"
    download_or_copy "$RESOURCES_DIR/$FASTA_37" "$LOCAL_RESOURCES/$FASTA_37"
    download_or_copy "$RESOURCES_DIR/$FASTA_37.fai" "$LOCAL_RESOURCES/${FASTA_37}.fai"
fi

# Uncompress FASTA files
gzip -dk "$LOCAL_RESOURCES/$FASTA_38"
[[ "$BUILD" == "Hg19" ]] && gzip -dk "$LOCAL_RESOURCES/$FASTA_37"


#########################################
# STEP 1: QC PARSE
#########################################
activate_conda hermes
log "Running QC Step 1"
OUTPUT_PARSED_GWAS="${LOCAL_OUTPUT}/parsed_gwas.tsv"
OUTPUT_GWAS2VCF_JSON="${LOCAL_OUTPUT}/gwas2vcf.json"
OUTPUT_STEP1_SUMMARY="${LOCAL_OUTPUT}/step1_summary.tsv"

Rscript /opt/scripts/01_parse_data.R \
  "$LOCAL_GWAS_FILE" \
  "$JSON_PATH" \
  "$OUTPUT_PARSED_GWAS" \
  "$OUTPUT_GWAS2VCF_JSON" \
  "$OUTPUT_STEP1_SUMMARY"


#########################################
# STEP 2: gwas2vcf
#########################################
activate_conda gwas2vcf
log "Running Step 2 gwas2vcf"
OUTPUT_VCF="${LOCAL_OUTPUT}/gwas.vcf.gz"
OUTPUT_STEP2_SUMMARY="${LOCAL_OUTPUT}/step2_summary.tsv"

REF="$LOCAL_RESOURCES/$FASTA_38_UNZ"
[[ "$BUILD" == "Hg19" ]] && REF="$LOCAL_RESOURCES/$FASTA_37_UNZ"

python /opt/gwas2vcf/main.py \
    --data "$OUTPUT_PARSED_GWAS" \
    --json "$OUTPUT_GWAS2VCF_JSON" \
    --id "foo" \
    --ref "$REF" \
    --out "$OUTPUT_VCF" \
    > "$OUTPUT_STEP2_SUMMARY" 2>&1


#########################################
# STEP 3: Liftover & annotate
#########################################
activate_conda hermes
log "Running Step 3 bcftools liftover & annotate"
DBSNP="$LOCAL_RESOURCES/$DBSNP_38"
CORES=$(nproc)
OUTPUT_VCF_B38="${LOCAL_OUTPUT}/gwas_b38.vcf.gz"
OUTPUT_STEP3_SUMMARY="${LOCAL_OUTPUT}/step3_summary.tsv"

if [[ "$BUILD" == "Hg19" ]]; then
    bcftools +liftover "$OUTPUT_VCF" \
        -- -s "$LOCAL_RESOURCES/$FASTA_37_UNZ" \
        -f "$LOCAL_RESOURCES/$FASTA_38_UNZ" \
        -c "$LOCAL_RESOURCES/$CHAIN_37_38" \
        2>> "$OUTPUT_STEP3_SUMMARY" \
        | bcftools view -Oz -o "$OUTPUT_VCF_B38.tmp.gz" -
    bcftools sort -Oz -o "$OUTPUT_VCF_B38.tmp.sorted.gz" "$OUTPUT_VCF_B38.tmp.gz"
    mv "$OUTPUT_VCF_B38.tmp.sorted.gz" "$OUTPUT_VCF_B38.tmp.gz"
else
    cp "$OUTPUT_VCF" "$OUTPUT_VCF_B38.tmp.gz"
fi

bcftools index -f --csi "$OUTPUT_VCF_B38.tmp.gz"

# Annotate with dbSNP
bcftools annotate -a "$DBSNP" -c ID -O z --threads "$CORES" \
    "$OUTPUT_VCF_B38.tmp.gz" -o "$OUTPUT_VCF_B38"

# Cleanup
rm -f "$OUTPUT_VCF_B38.tmp.gz" "$OUTPUT_VCF_B38.tmp.gz.csi"

bcftools index -f --csi "$OUTPUT_VCF_B38"
log "Step 3 completed"


#########################################
# STEP 4: Reporting
#########################################
log "Running QC Step 4 reporting"
EAF_REF="$LOCAL_RESOURCES/$EAF"
OUTPUT_CLEAN_GWAS="${LOCAL_OUTPUT}/gwas_b38_clean.tsv.gz"
OUTPUT_REPORT="${LOCAL_OUTPUT}/gwas_report.pdf"

Rscript /opt/scripts/04_qc_report.R \
    "$OUTPUT_VCF_B38" \
    "$JSON_PATH" \
    "$OUTPUT_STEP1_SUMMARY" \
    "$OUTPUT_STEP2_SUMMARY" \
    "$OUTPUT_STEP3_SUMMARY" \
    "$EAF_REF" \
    "/opt/scripts/04_qc_report.Rmd" \
    "$OUTPUT_CLEAN_GWAS" \
    "$OUTPUT_REPORT"


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
if [[ -d "$FILE_GUID" ]]; then
    log "Detected local output path — copying locally"
    mkdir -p "$FILE_GUID/images" "$FILE_GUID/tables"

    cp "$OUTPUT_CLEAN_GWAS" "$FILE_GUID/tables/"
    for file in "$LOCAL_OUTPUT/gwas_qc.html" "$LOCAL_OUTPUT/eaf_plot.png" "$LOCAL_OUTPUT/pz_plot.png" "$LOCAL_OUTPUT/qq_plot.png"; do
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
    for file in "$LOCAL_OUTPUT/gwas_qc.html" "$LOCAL_OUTPUT/eaf_plot.png" "$LOCAL_OUTPUT/pz_plot.png" "$LOCAL_OUTPUT/qq_plot.png"; do
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