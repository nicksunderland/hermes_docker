#!/usr/bin/env Rscript

#### DESCRIPTION #######################################
# This script performs the GWAS QC step (standalone version)
# Usage:
#   Rscript gwas_qc_step2.R <vcf> <meta_file> <step1_tbl> <step2_tbl> <step3_tbl> <eaf_ref> <report_tmp> <clean_tsv> <report>

#### SET INPUT #########################################
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 9) {
  stop("Usage: Rscript gwas_qc_step2.R <vcf> <meta_file> <step1_tbl> <step2_tbl> <step3_tbl> <eaf_ref> <report_tmp> <clean_tsv> <report>")
}

vcf        <- args[1]
meta_file  <- args[2]
step1_tbl  <- args[3]
step2_tbl  <- args[4]
step3_tbl  <- args[5]
eaf_ref    <- args[6]
report_tmp <- args[7]
clean_tsv  <- args[8]
report     <- args[9]

# testing ====
if (FALSE) {
  vcf        = '/myriadfs/projects/mol_cardio/Projects/p051-h3dcm/output/tmp_objects/parsed_gwas/dcm/eur/mixed/biovu-dcm-eur-mixed.vcf.gz'
  meta_file  = '/myriadfs/projects/mol_cardio/data/hermes_cohorts/biovu/biovu_h3pheno4_topmed_eur_mixed.json'
  step1_tbl  = '/myriadfs/projects/mol_cardio/Projects/p051-h3dcm/output/tmp_objects/parsed_gwas/dcm/eur/mixed/biovu-dcm-eur-mixed_step1_summary.tsv'
  step2_tbl  = '/myriadfs/projects/mol_cardio/Projects/p051-h3dcm/output/tmp_objects/parsed_gwas/dcm/eur/mixed/biovu-dcm-eur-mixed-step2_summary.txt'
  eaf_ref    = '/myriadfs/projects/mol_cardio/data/genome_references/1kg_v3/all_afreq_b38.tsv.gz'
  clean_gwas = ''
  report     = '/myriadfs/projects/mol_cardio/Projects/p051-h3dcm/output/figures/qc_reports/foo.pdf'
  report_tmp = '/myriadfs/projects/mol_cardio/Projects/p051-h3dcm/scripts/gwas_qc.Rmd'
}


# requirements ====
library(data.table)
library(jsonlite)
library(vcfR)
library(hexbin)
library(ggplot2)
library(ldscr)
library(stringr)


# read gwas ====
cat(sprintf("Reading GWAS file: %s\n", basename(meta_file)))
meta <- fromJSON(meta_file)
cat(sprintf("Reading GWAS file: %s\n", basename(vcf)))
vcf_dat <- read.vcfR(vcf)

# extract from vcf object ====
cat("converting VCF to data.table\n")
gwas <- vcfR2tidy(vcf_dat,
                  single_frame  = TRUE,
                  info_types    = TRUE,
                  format_fields = c("ES","SE","LP","SS","NC"),
                  format_types  = TRUE,
                  alleles       = FALSE,
                  gt_column_prepend = "")
gwas <- as.data.table(gwas$dat)
cols <- c(rsid="ID", chr="CHROM", bp_b38="POS", oa="REF", ea="ALT", eaf="AF", beta="ES", se="SE", nlog10p="LP", n="SS", ncase="NC")
gwas[, setdiff(names(gwas), unname(cols)) := NULL]
setnames(gwas, unname(cols), names(cols))
setcolorder(gwas, names(cols))
rm(vcf_dat)
gc()

# fix data
cat("fixing data types\n")
gwas[, chr := as.character(chr)]
gwas[, bp_b38 := as.integer(bp_b38)]
gwas[, eaf := as.numeric(eaf)]
gwas[, beta := as.numeric(beta)]
gwas[, se := as.numeric(se)]
gwas[, nlog10p := as.numeric(nlog10p)]
gwas[, p := 2 * pnorm(-abs(beta / se))]
gwas[, c("nlog10p") := NULL]
setcolorder(gwas, "p", after = "se")
gwas[, n := as.integer(n)]
gwas[, ncase := as.integer(ncase)]


#=============================================================================
# formatting allele frequency reference
#=============================================================================
cat("reading allele frequency reference\n")
freq <- data.table::fread(eaf_ref)
col_fn_list <- list(
  rsid      = as.character,
  chr       = as.character,
  pos_b38   = as.integer,
  ref       = as.character,
  alt       = as.character,
  afr_afreq = as.numeric,
  eur_afreq = as.numeric,
  amr_afreq = as.numeric,
  eas_afreq = as.numeric,
  sas_afreq = as.numeric
)

ancestry_map <- list(
  eur = "eur",
  eu = "eur",
  af = "afr",
  afr = "afr",
  aa = "afr",
  ssaf = "afr",
  ea = "eas",
  eas = "eas",
  hs = "amr",
  amr = "amr",
  mixed = "amr",
  sa = "sas",
  sas = "sas",
  gme = "sas",
  mid = "sas"
)
ancestry_parsed <- ancestry_map[[tolower(meta$ancestry)]]

freq[ , names(col_fn_list) := Map(function(fn, col) fn(col), col_fn_list[names(col_fn_list)], .SD), .SDcols = names(col_fn_list)]
freq <- freq[, .SD, .SDcols = c(names(freq)[!grepl("afreq", names(freq))], paste0(ancestry_parsed, "_afreq"))]
setnames(freq, paste0(ancestry_parsed, "_afreq"), "afreq")
freq <- freq[!is.na(afreq)]


#=============================================================================
# allele frequency analysis
#=============================================================================
cat("analysing allele frequency difference\n")

# get the allele frequencies for this ancestry
join_vec  <- stats::setNames(c("chr", "pos_b38", "ref", "alt"), c("chr", "bp_b38", "oa", "ea"))
rjoin_vec <- stats::setNames(c("chr", "pos_b38", "ref", "alt"), c("chr", "bp_b38", "ea", "oa"))
gwas[freq, `:=`(reference_id    = i.rsid,
                reference_afreq = i.afreq),     on = join_vec]
gwas[freq, `:=`(reference_id    = i.rsid,
                reference_afreq = 1 - i.afreq), on = rjoin_vec]

# absolute frequency difference cohort vs reference
freq_diff_thresh <- 0.2
frq_diff_labs <- stats::setNames(c(F,T), paste0(c("freq_diff_lt","freq_diff_gt"), freq_diff_thresh))
gwas[, freq_diff := factor(abs(eaf - reference_afreq) > freq_diff_thresh,
                           levels = unname(frq_diff_labs),
                           labels = names(frq_diff_labs))]

# print frequency difference counts
freq_diff_summary <- data.table(table(gwas$freq_diff, useNA = "always"))
freq_diff_summary[, pct := sprintf("%.1f%%", 100*(N / nrow(gwas)))]

# print summary to console
print(freq_diff_summary)

# plot the frequency differences
eaf_plot <- ggplot(mapping = aes(x = reference_afreq, y = eaf)) +
  geom_hex(data = gwas[freq_diff == paste0("freq_diff_lt", freq_diff_thresh)], bins = 100, fill = "royalblue", color = NA) +
  geom_point(data = gwas[freq_diff == paste0("freq_diff_gt", freq_diff_thresh)], aes(color = freq_diff)) +
  geom_abline(slope = 1, intercept =  freq_diff_thresh, linetype = "dashed", color = "red") +
  geom_abline(slope = 1, intercept = -1 * freq_diff_thresh, linetype = "dashed", color = "red") +
  scale_color_manual(values = stats::setNames(c("royalblue", "firebrick", "black"), levels(gwas$freq_diff))) +
  labs(x = paste0("Reference [", toupper(ancestry_parsed), "] allele frequency"), y = "Cohort allele frequency", fill = "Freq. difference") +
  lims(x = c(0,1), y = c(0,1)) +
  theme_minimal(base_size = 18) +
  theme(legend.position = "top",
        legend.title = element_blank())


#=============================================================================
# PZ plot
#=============================================================================
cat("assessing for analytical issues\n")

pz_dat <- gwas[, list(observed = -log10(p),
                      expected = -log10(2 * pnorm(abs(beta) / se, lower.tail=FALSE)))]

# generate PZ plot
pz_plot <- ggplot(pz_dat, aes(x = expected, y = observed)) +
  geom_hex(bins = 1000, color="darkblue") +
  geom_abline(slope=1, intercept=0, color="darkred", linetype = "dotted") +
  guides(fill = "none") +
  labs(x     = "P.ztest (-log10)",
       y     = "P (-log10)",
       fill  = NULL) +
  theme_minimal(base_size = 18) +
  theme(aspect.ratio = 1)


#=============================================================================
# population stratification analysis 1
#=============================================================================
cat("assessing for population stratification issues\n")

# get Z score and run LDSC
ldscr_dat <- gwas[, list(SNP = reference_id,
                         A1  = ea,
                         A2  = oa,
                         Z   = beta / se,
                         N   = n)]

# try to run LDSC
tryCatch({
  ldsc_res <- ldsc_h2(munged_sumstats = ldscr_dat, ancestry = toupper(ancestry_parsed))
},
error = function(e) {
  ldsc_res <<- list(intercept = NA_real_)
})


#=============================================================================
# population stratification analysis 2
#=============================================================================
cat("calculating lambda-GC\n")

# data for the QQ plot
gwas[, `:=`(chisq          = qchisq(log(p), df = 1, lower.tail = FALSE, log.p = TRUE))] # https://stackoverflow.com/questions/40144267/calculating-miniscule-numbers-for-chi-squared-distribution-numerical-precisio
gwas[, `:=`(lambda         = median(chisq) / qchisq(0.5, 1))]
gwas[, `:=`(adj_chisq_gc   = chisq/lambda)]
gwas[, `:=`(adj_p_gc       = pchisq(adj_chisq_gc, 1, lower.tail=FALSE))]
gwas[, `:=`(adj_se_gc      = se * sqrt(lambda[1]))]
gwas[, `:=`(adj_chisq_ldsc = chisq/ldsc_res$intercept)]
gwas[, `:=`(adj_p_ldsc     = pchisq(adj_chisq_ldsc, 1, lower.tail=FALSE))]
gwas[, `:=`(adj_se_ldsc    = se * sqrt(ldsc_res$intercept))]

qq_data <- rbind(
  # unadjusted
  gwas[, .(stat     = "Unadjusted",
           value    = NA_real_,
           observed = -log10(sort(p, decreasing=FALSE, na.last=TRUE)),
           expected = -log10(ppoints(.N)),
           clower   = -log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:.N, shape2 = .N:1)),
           cupper   = -log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:.N, shape2 = .N:1)))],
  # GC lambda adjusted
  gwas[, .(stat     = "GC lambda",
           value    = lambda[1],
           observed = -log10(sort(adj_p_gc, decreasing=FALSE, na.last=TRUE)),
           expected = -log10(ppoints(.N)),
           clower   = NA_real_,
           cupper   = NA_real_)],
  # LDSC intercept adjusted
  gwas[, .(stat     = "LDSC intercept",
           value    = ldsc_res$intercept,
           observed = -log10(sort(adj_p_ldsc, decreasing=FALSE, na.last=TRUE)),
           expected = -log10(ppoints(.N)),
           clower   = NA_real_,
           cupper   = NA_real_)]
)
setorder(qq_data, expected)

# generate QQ-plot
# axis labels
log10Pe <- expression(paste("Expected -log"[10], plain(P)))
log10Po <- expression(paste("Observed -log"[10], plain(P)))

# lambda labels
labels <- qq_data[stat != "Unadjusted",
                  list(value    = value[1],
                       expected = 2.0,
                       observed = .GRP-1 + 0.25), by = "stat"]
labels[, label := sprintf("%s = %.3f", stat, value)]

# plot
# colors
adjustment <- NULL
point_colors <- if (is.null(adjustment)) {
                  c(Unadjusted = "blue", `GC lambda` = "darkgray", `LDSC intercept` = "black")
                } else if (adjustment=="lambda") {
                  c(Unadjusted = "darkgray", `GC lambda` = "blue", `LDSC intercept` = "black")
                } else if (adjustment=="ldsc") {
                  c(Unadjusted = "darkgray", `GC lambda` = "black", `LDSC intercept` = "blue")
                }

qq_plot <- qq_data |>
  ggplot(aes(x = expected, y = observed, fill = stat)) +
  geom_hex(bins = 1000) +
  #geom_point(size = 0.5) +
  geom_ribbon(aes(x = expected, ymin = clower, ymax = cupper), alpha = 0.1, color="transparent") +
  geom_abline(slope=1, intercept=0, color="darkred", linetype = "dotted") +
  geom_text(data = labels, aes(x=expected, y=observed, label=label), hjust = 0, color="black", show.legend = FALSE) +
  scale_fill_manual(values = point_colors) +
  labs(x = log10Pe,
       y = log10Po,
       color = "Adjustment") +
  theme_minimal(base_size = 18) +
  theme(aspect.ratio    = 1,
        legend.position = "top")


#=============================================================================
# extract then save the clean GWAS and summary table
#=============================================================================

# if using adjustment with lambda GC or LDSC intercept, switch the P and SE columns
if (!is.null(adjustment)) {
  if (adjustment=="lambda") {
    gwas[, c("se", "p") := .(adj_se_gc, adj_p_gc)]
  } else if (adjustment=="ldsc") {
    gwas[, c("se", "p") := .(adj_se_ldsc, adj_p_ldsc)]
  }
}

# sorted allele ID column
gwas[, varid := paste0(chr, ":", bp_b38, "_", pmin(ea, oa), "_", pmax(ea, oa))]
gwas[, rsid := fcoalesce(rsid, reference_id)]
data.table::setcolorder(gwas, c("varid", "rsid"))

# clean
gwas[, c("reference_id", "freq_diff",
         "chisq", "lambda", "adj_chisq_gc", "adj_p_gc", "adj_se_gc",
         "adj_chisq_ldsc", "adj_p_ldsc", "adj_se_ldsc") := NULL]

# save gwas
cat("saving clean GWAS file\n")
fwrite(gwas, clean_tsv, sep = "\t")

# get summary data for report from step1 log
summary   <- fread(step1_tbl)

# get harmonisation data for report from step2 log
step2_lines <- readLines(step2_tbl)


# function to extract a numeric stat by name
extract_stat <- function(stat_name){
  line <- step2_lines[str_detect(step2_lines, paste0(stat_name, ":"))]
  if(length(line) == 0) return(NA_integer_)
  as.integer(str_match(line, paste0(stat_name, ": (\\d+)"))[,2])
}

# get gwas2vcf details
log_info <- list()
log_info$total_variants <- extract_stat("Total variants")
log_info$variants_harmonised <- extract_stat("Variants harmonised")
log_info$variants_discarded <- extract_stat("Variants discarded during harmonisation")
log_info$alleles_switched <- extract_stat("Alleles switched")
log_info$normalised_variants <- extract_stat("Normalised variants")
log_info$skipped_variants <- extract_stat("Skipped")
harm_summary <- as.data.table(log_info)

# liftover log
step3_lines <- readLines(step3_tbl)
line <- step3_lines[str_detect(step3_lines, "Lines\\s+total/swapped/reference added/rejected:")]
if(length(line) == 0){
  liftover_stats <- data.table(total=NA_integer_, swapped=NA_integer_,
                               reference_added=NA_integer_, rejected=NA_integer_)
} else {
  # Extract the four numbers
  numbers <- str_match(line, "(\\d+)/(\\d+)/(\\d+)/(\\d+)")[1,2:5]
  liftover_stats <- data.table(
    total = as.integer(numbers[1]),
    swapped = as.integer(numbers[2]),
    reference_added = as.integer(numbers[3]),
    rejected = as.integer(numbers[4])
  )
}

# save png files
eaf_plot_fp <- "/tmp/output/eaf_plot.png"
qq_plot_fp  <- "/tmp/output/qq_plot.png"
pz_plot_fp  <- "/tmp/output/pz_plot.png"
ggsave(eaf_plot_fp, eaf_plot, width = 6, height = 6, dpi = 300, bg = "white")
ggsave(qq_plot_fp,  qq_plot,  width = 6, height = 6, dpi = 300, bg = "white")
ggsave(pz_plot_fp,  pz_plot,  width = 6, height = 6, dpi = 300, bg = "white")

# create report
dir.create(dirname(report), showWarnings = FALSE, recursive = TRUE)
rmarkdown::render(input             = report_tmp,
                  knit_root_dir     = dirname(report),
                  intermediates_dir = dirname(report),
                  output_dir        = dirname(report),
                  output_file       = basename(report),
                  envir             = globalenv())