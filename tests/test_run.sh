# server command ???
#ContainerProperties:
  #  Image: nicksunderland/hermes_docker:latest
#  Command:
#    - '/bin/bash'
#    - '/opt/run.sh'
#    - 'Ref::s3-path'           # $1 → GWAS input (S3 path)
#    - 'Ref::config-json'       # $2 → JSON file or inline JSON string (the gwas meta-data)
#    - 'Ref::s3-resources-dir'  # $3 → reference resources dir (S3 hermes repository directory)
#    - 'Ref::file-guid'         # $4 → GUID of file

# run test 20,000 variants
DATA_DIR="$HOME/git/hermes_docker/tests/data"
RES_DIR="$HOME/git/hermes_docker/tests/resources"
OUT_DIR="$HOME/git/hermes_docker/tests/output"
TMP_DIR="$HOME/git/hermes_docker/tests/tmp"
mkdir -p "$TMP_DIR" "$OUT_DIR"
docker run --rm \
  --platform linux/amd64 \
  -v "$DATA_DIR":/data \
  -v "$RES_DIR":/resources \
  -v "$OUT_DIR":/output \
  -v "$TMP_DIR":/tmp \
  nicksunderland/hermes_docker:latest \
  /data/test_data.tsv.gz \
  /data/config.json \
  /resources \
  /output/test_output
rm -rf "$TMP_DIR"

# Run AOU
DATA_DIR="$HOME/git/hermes_docker/tests/data"
RES_DIR="$HOME/git/hermes_docker/tests/resources"
OUT_DIR="$HOME/git/hermes_docker/tests/output"
TMP_DIR="$HOME/git/hermes_docker/tests/tmp"
mkdir -p "$TMP_DIR" "$OUT_DIR"
docker run --rm \
  --platform linux/amd64 \
  -v "$DATA_DIR":/data \
  -v "$RES_DIR":/resources \
  -v "$OUT_DIR":/output \
  -v "$TMP_DIR":/tmp \
  nicksunderland/hermes_docker:latest \
  /data/aou_eur_p1_combo.tsv.gz \
  /data/config_hermes_aou.json \
  /resources \
  /output/aou
rm -rf "$TMP_DIR"

# Run UKBB
DATA_DIR="$HOME/git/hermes_docker/tests/data"
RES_DIR="$HOME/git/hermes_docker/tests/resources"
OUT_DIR="$HOME/git/hermes_docker/tests/output"
TMP_DIR="$HOME/git/hermes_docker/tests/tmp"
mkdir -p "$TMP_DIR" "$OUT_DIR"
docker run --rm \
  --platform linux/amd64 \
  -v "$DATA_DIR":/data \
  -v "$RES_DIR":/resources \
  -v "$OUT_DIR":/output \
  -v "$TMP_DIR":/tmp \
  nicksunderland/hermes_docker:latest \
  /data/GWAS_UKB_Pheno1_PREV_MIX_EUR_BOTH_20250114_CJ.tsv.gz \
  /data/config_hermes_ukbb.json \
  /resources \
  /output/ukbb
rm -rf "$TMP_DIR"