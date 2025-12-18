# Build a 1KG Phase 3 allele-frequency reference on GRCh37/Hg19
# Source data: plink2 "all_phase3" (GRCh37) V2 array from https://www.cog-genomics.org/plink/2.0/resources#phase3_1kg
# Before running:
#   1) Download and unzip the GRCh37 bundle; set `ref` below to the prefix of the .pgen/.pvar/.psam (without extension).
#   2) Ensure plink2 is on PATH or update `plink2` to the full binary path.
# Output:
#   all_afreq_b37.tsv.gz with columns: rsid, chr, pos_b37, ref, alt, <ancestry>_afreq

library(data.table)

plink2     <- "/usr/local/plink/plink2"
ref        <- "/Users/xx20081/Documents/local_data/genome_reference/1kg.p3_grc37_all/all_phase3"
output_dir <- "/Users/xx20081/git/hermes_docker/tests/resources"

ancestries <- c("EUR","AFR","SAS","EAS","AMR")

# Create ancestry-specific subsets (keeps files small and speeds up freq calc)
for (a in ancestries) {
  cmd <- paste(
    plink2,
    "--pfile vzs", ref,
    "--keep-if", paste0("SuperPop==", a),
    "--rm-dup force-first",
    "--set-missing-var-ids '@:#_$r_$a'",
    "--new-id-max-allele-len 662",
    "--make-pgen vzs",
    "--out", paste0(tolower(a), "_1kg_p3_b37")
  )
  system(cmd)
}

data <- data.table()

for (a in ancestries) {
  message("📊 Calculating frequencies for ", a)
  tmp <- tempfile()
  cmd <- paste(
    plink2,
    "--pfile vzs", file.path("/Users/xx20081/git/hermes_docker", paste0(tolower(a), "_1kg_p3_b37")),
    "--freq", "cols=chrom,pos,ref,alt,altfreq",
    "--out", tmp
  )
  system(cmd)
  freq <- fread(paste0(tmp, ".afreq"))
  unlink(paste0(tmp, ".afreq"))

  # Expand multi-allelic rows
  multi <- freq[grepl(",", ALT_FREQS), ]
  freq  <- freq[!grepl(",", ALT_FREQS), ][, ALT_FREQS := as.numeric(ALT_FREQS)]
  multi <- multi[
    , .(
      ALT = unlist(strsplit(ALT, ",")),
      ALT_FREQS = as.numeric(unlist(strsplit(ALT_FREQS, ",")))
    ),
    by = c("ID", "#CHROM", "POS", "REF")
  ]
  freq <- rbind(freq, multi, use.names = TRUE, fill = TRUE)
  freq <- freq[!is.na(ALT_FREQS) & ALT_FREQS >= 0 & ALT_FREQS <= 1]

  setnames(freq,
           c("ID", "#CHROM", "POS", "REF", "ALT", "ALT_FREQS"),
           c("rsid","chr","pos_b37","ref","alt",paste0(tolower(a), "_afreq")))

  if (nrow(data) == 0) {
    data <- freq
  } else {
    data <- merge(data, freq, by = c("rsid","chr","pos_b37","ref","alt"), all = TRUE)
  }
}

fwrite(data, file.path(output_dir, "all_afreq_b37.tsv.gz"), sep = "\t")
