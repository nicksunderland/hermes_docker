# first download the all_hg38.pgen / pvar / psam from the plink2 website
# https://www.cog-genomics.org/plink/2.0/resources#phase3_1kg

library(data.table)
plink2="/usr/local/plink/plink2" # plink binary
ref="path_to_downloaded_file/all_hg38" # plink pgen/pvar/psam prefix
output_dir="/path_to_resources_output_directory/resources"

ancestries = c("EUR","AFR","SAS","EAS","AMR")

for (a in ancestries) {

  cmd <- paste(
    plink2,
    "--pfile vzs", ref,
    "--keep-if", paste0("SuperPop==", a),
    "--rm-dup force-first",
    "--set-missing-var-ids '@:#_$r_$a'",
    "--new-id-max-allele-len 662",
    "--make-pgen vzs",
    "--out", paste0(tolower(a), "_1kg_p3_b38")
  )
  system(cmd)

}

# rsid    chr     pos_b37 ref     alt     afr_afreq       eur_afreq       amr_afreq       eas_afreq       sas_afreq
data = data.table()

#loop ancestries again getting the allele frequencies, rsid, chr, pos, alt, afr_afreq       eur_afreq       amr_afreq       eas_afreq       sas_afreq
for (a in ancestries) {
  message("📊 Calculating frequencies for ", a)
  tmp <- tempfile()
  cmd <- paste(
    plink2,
    "--pfile vzs", file.path(dirname(ref), paste0(tolower(a), "_1kg_p3_b38")),
    "--freq", "cols=chrom,pos,ref,alt,altfreq",
    "--out", tmp
  )
  system(cmd)
  freq <- fread(paste0(tmp, ".afreq"))
  unlink(paste0(tmp, ".afreq"))
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

  setnames(freq, c("ID", "#CHROM", "POS", "REF", "ALT", "ALT_FREQS"), c("rsid","chr","pos_b38","ref","alt",paste0(tolower(a), "_afreq")))
  if (nrow(data)==0) {
    data = freq
  } else {
    data = merge(data, freq, by=c("rsid","chr","pos_b38","ref","alt"), all=TRUE)
  }

}
fwrite(data, file.path(output_dir, "all_afreq_b38.tsv.gz"), sep="\t")