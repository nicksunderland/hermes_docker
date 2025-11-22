#!/bin/bash

# Safer bash behavior
set -euo pipefail

# --- CONFIGURATION ---
# Define the working directory here
WORK_DIR="$HOME/resources"

# Create directory if it doesn't exist and enter it
mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

echo "📂 Working in: $(pwd)"

# Check if samtools is actually available
if command -v samtools &> /dev/null; then
    echo "✅ Samtools found. Indexing will be performed."
    RUN_SAMTOOLS=true
else
    echo "⚠️  WARNING: 'samtools' command not found!"
    echo "   - Files will be downloaded and 'chr' prefixes stripped."
    echo "   - Indexing (.fai) and Dictionary (.dict) creation will be SKIPPED."
    RUN_SAMTOOLS=false
fi

# --- STEP 1: DOWNLOAD B37 (Legacy) ---
echo "⬇️  Checking/Downloading Build 37 files..."

if [ ! -f "human_g1k_v37.fasta" ]; then
    # Check if the gz exists, if not download it
    if [ ! -f "human_g1k_v37.fasta.gz" ]; then
        wget -c http://fileserve.mrcieu.ac.uk/ref/2.8/b37/human_g1k_v37.fasta.gz
    fi
    echo "📦 Decompressing human_g1k_v37.fasta.gz..."
    gzip -d human_g1k_v37.fasta.gz
else
    echo "   human_g1k_v37.fasta already exists. Skipping."
fi

# Download indices for b37 if missing
[ ! -f "human_g1k_v37.fasta.fai" ] && wget -c http://fileserve.mrcieu.ac.uk/ref/2.8/b37/human_g1k_v37.fasta.fai
[ ! -f "human_g1k_v37.dict" ] && wget -c http://fileserve.mrcieu.ac.uk/ref/2.8/b37/human_g1k_v37.dict

# --- STEP 2: DOWNLOAD HG38 ---
echo "⬇️  Checking/Downloading HG38 files..."

HG38_URL="https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0"
IN_FASTA="Homo_sapiens_assembly38.fasta"

[ ! -f "$IN_FASTA" ] && wget -c "${HG38_URL}/${IN_FASTA}"
[ ! -f "${IN_FASTA}.fai" ] && wget -c "${HG38_URL}/${IN_FASTA}.fai"
[ ! -f "${IN_FASTA%.*}.dict" ] && wget -c "${HG38_URL}/${IN_FASTA%.*}.dict"

# --- STEP 3: CONVERT HG38 (STRIP CHR) ---
OUT_FASTA="Homo_sapiens_assembly38_nochr.fasta"

if [ ! -f "$OUT_FASTA" ]; then
    echo "🔧 Stripping 'chr' prefixes from $IN_FASTA..."
    # The awk command replaces '>chr' with '>' in header lines only
    awk '/^>/{sub(/^>chr/,">",$0)}1' "$IN_FASTA" > "$OUT_FASTA"
else
    echo "   $OUT_FASTA already exists. Skipping awk conversion."
fi

# --- STEP 4: INDEX NEW FASTA (if samtools installed - recommended)
if [ "$RUN_SAMTOOLS" = true ]; then
    echo "🧩 Indexing new 'no-chr' FASTA with samtools..."

    if [ ! -f "${OUT_FASTA}.fai" ]; then
        samtools faidx "$OUT_FASTA"
    fi

    if [ ! -f "${OUT_FASTA%.*}.dict" ]; then
        echo "📖 Creating sequence dictionary..."
        samtools dict "$OUT_FASTA" -o "${OUT_FASTA%.*}.dict"
    fi
else
    echo "🛑 Samtools not installed. Skipping .fai and .dict generation for $OUT_FASTA."
fi

echo "✅ All tasks complete."
ls -lh "$OUT_FASTA"* 2>/dev/null || true