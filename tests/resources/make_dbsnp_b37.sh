#!/bin/bash
# Build dbSNP build 157 for GRCh37/Hg19 with simple chromosome renaming
# Warning: large download

set -euo pipefail

# configurable paths
TARGET_DIR="/Users/xx20081/git/hermes_docker/tests/resources"
PREFIX="GCF_000001405.25"
BCFTOOLS="${BCFTOOLS:-bcftools}"

mkdir -p "$TARGET_DIR"
cd "$TARGET_DIR"

echo "📥 Downloading dbSNP b157 GRCh37 (${PREFIX}.*)..."
wget -c "https://ftp.ncbi.nih.gov/snp/archive/b157/VCF/${PREFIX}.gz"
wget -c "https://ftp.ncbi.nih.gov/snp/archive/b157/VCF/${PREFIX}.gz.tbi"
wget -c "https://ftp.ncbi.nih.gov/snp/archive/b157/VCF/${PREFIX}.gz.md5"
wget -c "https://ftp.ncbi.nih.gov/snp/archive/b157/VCF/${PREFIX}.gz.tbi.md5"

echo "🔍 Verifying MD5 checksums..."
md5sum -c "${PREFIX}.gz.md5"
md5sum -c "${PREFIX}.gz.tbi.md5"

echo "🧬 Creating GRCh37 chromosome rename map..."
cat > hg19_rename_chrom_names.tsv <<'EOF'
NC_000001.10    1
NC_000002.11    2
NC_000003.11    3
NC_000004.11    4
NC_000005.9     5
NC_000006.11    6
NC_000007.13    7
NC_000008.10    8
NC_000009.11    9
NC_000010.10    10
NC_000011.9     11
NC_000012.11    12
NC_000013.10    13
NC_000014.8     14
NC_000015.9     15
NC_000016.9     16
NC_000017.10    17
NC_000018.9     18
NC_000019.9     19
NC_000020.10    20
NC_000021.8     21
NC_000022.10    22
NC_000023.10    X
NC_000024.9     Y
NC_012920.1     MT
EOF

echo "🔧 Renaming chromosome names and compressing..."
$BCFTOOLS annotate \
  --rename-chrs hg19_rename_chrom_names.tsv \
  --output-type z \
  --output dbSNP157_b37_clean.vcf.gz \
  "${PREFIX}.gz"

echo "🧩 Indexing cleaned dbSNP file..."
$BCFTOOLS index -t dbSNP157_b37_clean.vcf.gz

PRIMARY_CHROMS=$(bcftools view -h dbSNP157_b37_clean.vcf.gz \
  | grep '^##contig=' \
  | sed -E 's/.*ID=([^,>]+).*/\1/' \
  | grep -E '^(?:[1-9]|1[0-9]|2[0-2]|X|Y|MT)$')

echo "🧬 Splitting dbSNP primary chromosomes only..."
while read -r chr; do
  echo "→ chr${chr}"
  bcftools view \
    -r "$chr" \
    -O b \
    -o "dbSNP157_b37_chr${chr}.bcf" \
    dbSNP157_b37_clean.vcf.gz

  bcftools index -f "dbSNP157_b37_chr${chr}.bcf"
done <<< "$PRIMARY_CHROMS"

echo "✅ Done. Outputs:"
ls -lh dbSNP157_b37_chr*
echo "Place dbSNP157_b37_chr*.bcf(.csi) files into tests/resources/ for the pipeline."
