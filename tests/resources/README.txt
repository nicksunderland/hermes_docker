This resources folder needs to contain the files:

FASTA files:
Homo_sapiens_assembly38_nochr.fasta.fai
Homo_sapiens_assembly38_nochr.fasta.gz
human_g1k_v37.fasta.fai
human_g1k_v37.fasta.gz

dbSNP VCF files:
dbSNP157_b38_clean.vcf.gz
dbSNP157_b38_clean.vcf.gz.tbi

Liftover chain file:
hg19ToHg38.over.chain.gz

Ancestry-specific allele frequency reference file:
all_afreq_b38.tsv.gz

Either, I can try to transfer the files, or you can attempt to create them using the scripts below:
make_afreq.R
make_chain.sh
make_dbsnp.sh
make_fasta.sh
