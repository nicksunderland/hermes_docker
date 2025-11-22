#!/bin/bash

# Create folder 'resources' if it doesn't exist, and go into it
mkdir -p resources
cd resources

echo "📂 Downloading chain files into 'resources'..."

# Download files (flag -c allows resuming if interrupted)
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
wget -c http://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh37_to_GRCh38.chain.gz
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHs1.over.chain.gz

echo "✅ Done! Files are in $(pwd)"
ls -lh