#!/usr/bin/env Rscript
# GWAS QC Step 1: Missingness & basic parsing
# Usage: Rscript gwas_qc_step1.R <gwas_file> <meta_file> <parsed_gwas> <gwas2vcf_json> <step1_summary>

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 5){
  stop("Usage: Rscript gwas_qc_step1.R <gwas_file> <meta_file> <parsed_gwas> <gwas2vcf_json> <step1_summary>")
}
gwas_file     <- args[1]
meta_file     <- args[2]
parsed_gwas   <- args[3]
gwas2vcf_json <- args[4]
step1_summary <- args[5]


# testing ====
if (FALSE) {
  gwas_file     = '/Users/xx20081/git/hermes_docker/tests/data/aou_eur_p1_combo.tsv.gz'
  meta_file     = '/Users/xx20081/git/hermes_docker/tests/data/config_hermes_aou.json'
  parsed_gwas   = '/Users/xx20081/git/hermes_docker/tests/output/test_data_parsed.tsv.gz'
  gwas2vcf_json = '/Users/xx20081/git/hermes_docker/tests/output/test_gwas2vcf.json'
  step1_summary = '/Users/xx20081/git/hermes_docker/tests/output/test_summary1.txt'
}


# requirements ====
library(data.table)
library(jsonlite)


# read gwas & meta-data ====
cat(sprintf("Reading GWAS meta-data file: %s\n", basename(meta_file)))
meta <- fromJSON(meta_file)
cat(sprintf("Reading GWAS file: %s\n", basename(gwas_file)))
gwas <- fread(gwas_file, select = unname(unlist(meta$column_map)), col.names = names(meta$column_map))


# summary table to be iteratively filled ====
summary <- data.table()


#=============================================================================
# check column data missingness
#=============================================================================
cat("checking data missingness\n")

req_cols  <- c("rsid", "chr", "pos", "alt", "reference", "eaf", "beta", "stdErr", "pval", "n", "ncase", "ncontrol", "imputed", "info")
miss_cols <- setdiff(req_cols, names(gwas))
if (length(miss_cols) > 0) {
  cat("[!] missing columns: ", paste0(miss_cols, collapse = ", "), "\n")
  for (col in miss_cols) {
    gwas[, (col) := NA]
  }
}

num_na_summary <- gwas[, lapply(.SD, function(col) c(num = .N, num_na = sum(is.na(col)), pct_na = sprintf("%.1f%%", 100 * (sum(is.na(col)) / .N))))]
summary <- rbind(summary,
                 data.table(column   = names(gwas),
                            type     = sapply(gwas, typeof),
                            num      = as.integer(sapply(num_na_summary, `[[`, 1)),
                            num_na   = as.integer(sapply(num_na_summary, `[[`, 2)),
                            pct_na   = sapply(num_na_summary, `[[`, 3)))
print(summary)


#=============================================================================
# check input data validity
#=============================================================================
cat("checking data validity\n")

# numeric checker
check_num <- function(x, 
                      lb = -Inf, ub = Inf,        
                      lb_incl = FALSE, ub_incl = FALSE,
                      is_int = FALSE) {
  val <- suppressWarnings(as.numeric(x))
  res <- !is.na(val) & is.finite(val)
  if (lb_incl) res <- res & val >= lb
  else         res <- res & val >  lb
  if (ub_incl) res <- res & val <= ub
  else         res <- res & val <  ub
  if (is_int) {
    res <- res & (abs(val - round(val)) < .Machine$double.eps^0.5)
  }
  return(res)
}

# per column functions that return true if that row is valid
col_fn_list <- list(
  rsid      = function(x) grepl("^rs[0-9]+|^[0-9]+:[0-9]+.*", x),
  chr       = function(x) x %in% c(as.character(1:22), "X", "Y"),
  pos       = function(x) check_num(x, lb=0, lb_incl=F, is_int=T),
  alt       = function(x) grepl("^[ACTG]+$", x),
  reference = function(x) grepl("^[ACTG]+$", x),
  eaf       = function(x) check_num(x, lb=0, ub=1),
  beta      = function(x) check_num(x, lb=-20, ub=20),
  stdErr    = function(x) check_num(x, lb=0),
  pval      = function(x) check_num(x, lb=0, ub=1, lb_incl=T, ub_incl=T),
  n         = function(x) check_num(x, lb=0, is_int=T),
  ncase     = function(x) check_num(x, lb=0, is_int=T),
  ncontrol  = function(x) check_num(x, lb=0, is_int=T),
  imputed   = function(x) !is.na(x) & x %in% list(TRUE, FALSE, 0, 1),
  info      = function(x) check_num(x, lb=0, ub=1, lb_incl=T, ub_incl = T)
)

# apply the functions to the corresponding columns and create summary table
valid_summary <- gwas[, Map(function(fn, col) {
  valid <- sum(fn(col), na.rm = TRUE)
  valid_pct <- sprintf("%.1f%%", 100 * (valid / .N))
  c(valid = valid, valid_pct = valid_pct)
}, col_fn_list[names(.SD)], .SD)]

# add to main summary table
summary <- cbind(summary,
                 data.table(valid     = as.integer(sapply(valid_summary, `[[`, 1)),
                            valid_pct = sapply(valid_summary, `[[`, 2)))
print(summary)


#=============================================================================
# attempt recoding / data fixes
#=============================================================================
cat("formatting data\n")

# (re)code chromosome column for fasta
cat("\t - unique chr:", paste0(unique(gwas$chr), collapse = ", "), "\n")
gwas[, chr := sub("(?i)^chr", "", chr)]
gwas[, chr := toupper(as.character(chr))]
ok_chrom <- c(as.character(1:22), "X", "Y")
gwas[, chr := fcase(chr %in% ok_chrom, chr,
                    chr == "23",       "X",
                    chr == "24",       "Y",
                    default = NA_character_)]

#=============================================================================
# Pre-processing: Handle censored N values (e.g., "<=40")
#=============================================================================
cat("handling censored values in sample size columns\n")

n_cols <- c("n", "ncase", "ncontrol")
for (col in n_cols) {
  if (col %in% names(gwas)) {
    if (is.character(gwas[[col]])) {
      cat(sprintf("\t - cleaning column '%s' (removing '<=' etc)\n", col))
      gwas[, (col) := gsub("[^0-9.]", "", get(col))]
    }
  }
}

# parse base position and n-sample columns to integer
integer_cols <- c("pos", "n", "ncase", "ncontrol")
gwas[, (integer_cols) := lapply(.SD, as.integer), .SDcols = integer_cols]

# if total N is NA, then set to N in the meta-data
gwas[is.na(n), n := meta$totalSampleSize]
gwas[is.na(ncase), ncase := meta$cases]
gwas[is.na(ncontrol), ncontrol := meta$totalSampleSize - meta$cases]

# apply integer filters
for (col in integer_cols) {
    gwas[!col_fn_list[[col]](get(col)), (col) := NA_integer_]
}

# parse frequency, beta, se, p, info columns to numeric
numeric_cols <- c("eaf", "beta", "stdErr", "pval", "info")
if ("log10p" %in% names(gwas)) {
  numeric_cols <- c(numeric_cols, "log10p")
}
gwas[, (numeric_cols) := lapply(.SD, as.numeric), .SDcols = numeric_cols]

# detect 0–100 frequencies and rescale only if clearly percentage-based
eaf_vals <- gwas$eaf[is.finite(gwas$eaf)]
if (length(eaf_vals)) {
  pct_gt1_lt100  <- mean(eaf_vals > 1 & eaf_vals <= 100)
  if (pct_gt1_lt100 > 0.9) {
    cat(sprintf("[!] EAF likely 0–100 (%.1f%% of finite values); scaling to 0–1\n", 100 * pct_gt1_lt100))
    gwas[, eaf := eaf / 100]
  } else {
    cat(sprintf("[i] EAF appears 0–1 (%.1f%% ≤1); leaving as-is\n", 100 * pct_le1))
  }
}

# sometime -log10P might be provided - need to recode
if (range(gwas$pval[!is.na(gwas$pval) & is.finite(gwas$pval)])[2] > 1) {
  cat("[!] -log10(P) values detected, max value in pval column =", range(gwas$pval[!is.na(gwas$pval) & is.finite(gwas$pval)])[2], "- recalculating P-values from beta and se\n")
  gwas[, pval := 2 * pnorm(-abs( beta / stdErr ))]
}

# apply numeric filters
for (col in numeric_cols) {
    gwas[!col_fn_list[[col]](get(col)), (col) := NA_real_]
}

# recode pvalue=0 to minimum machine precision and remove p>1 or p<0
gwas[, pval := fcase(pval == 0, .Machine$double.xmin,
                     pval > 0 & pval <= 1, pval,
                     default = NA_real_)]

# alleles must be characters
allele_cols <- c("alt", "reference")
gwas[, (allele_cols) := lapply(.SD, function(x) toupper(as.character(x))), .SDcols = allele_cols]
for (col in allele_cols) {
    gwas[!col_fn_list[[col]](get(col)), (col) := NA_character_]
}

# find indels and and number to summary
indel_idx <- which(nchar(gwas$alt) > 1 | nchar(gwas$reference) > 1)
summary[grep("^(alt|reference)$", column), indels := length(indel_idx)]

# summarise the counts during this recoding process
postfix_na_summary <- gwas[, lapply(.SD, function(col) c(postfix_valid = sum(!is.na(col)), pct_na = sprintf("%.1f%%", 100 * (sum(!is.na(col)) / .N))))]

# add to overall summary table
summary <- cbind(summary,
                 data.table(postfix_valid     = as.integer(sapply(postfix_na_summary, `[[`, 1)),
                            postfix_valid_pct = sapply(postfix_na_summary, `[[`, 2)))

# remove invalid rows (those containing NAs)
essential_cols  <- c("chr", "pos", "alt", "reference", "eaf", "beta", "stdErr", "pval", "n", "ncase", "ncontrol")
na_variants <- gwas[!stats::complete.cases(gwas[, mget(essential_cols)]) ]
gwas        <- gwas[ stats::complete.cases(gwas[, mget(essential_cols)]) ]

# print updated summary to console
print(summary)
data.table::fwrite(summary, step1_summary, sep = "\t")

# save gwas
cat("saving parsed GWAS file\n")
data.table::fwrite(gwas, parsed_gwas, sep = "\t")

# save column formatting for gwas2vcf
json_data <- list(
  chr_col       = which(names(gwas) == "chr") - 1,
  pos_col       = which(names(gwas) == "pos") - 1,
  snp_col       = which(names(gwas) == "rsid") - 1,
  ea_col        = which(names(gwas) == "alt") - 1,
  oa_col        = which(names(gwas) == "reference") - 1,
  beta_col      = which(names(gwas) == "beta") - 1,
  se_col        = which(names(gwas) == "stdErr") - 1,
  ncase_col     = which(names(gwas) == "ncase") - 1,
  ncontrol_col  = which(names(gwas) == "ncontrol") - 1,
  pval_col      = which(names(gwas) == "pval") - 1,
  eaf_col       = which(names(gwas) == "eaf") - 1,
  delimiter     = "\t",
  header        = TRUE,
  build         = fcase(
    grepl("hg[-_.]?19|grch[-_.]?37|b(uild)?[-_.]?37|ncbi[-_.]?37", meta$referenceGenome, ignore.case = TRUE), "GRCh37",
    grepl("hg[-_.]?38|grch[-_.]?38|b(uild)?[-_.]?38|ncbi[-_.]?38", meta$referenceGenome, ignore.case = TRUE), "GRCh38",
    default = NA_character_
  )
)

# check we created a valid file
if (is.na(json_data$build)) {
  stop(paste0("ERROR: Could not determine Genome Build from metadata string: '", meta$referenceGenome, "'"))
}

missing_indices <- names(json_data)[sapply(json_data, length) == 0]
if (length(missing_indices) > 0) {
  stop(paste0("ERROR: The following required columns were missing from the parsed GWAS data object when creating the JSON map: ", 
              paste(missing_indices, collapse = ", ")))
}

# write JSON to disk
jsonlite::write_json(json_data, path = gwas2vcf_json, auto_unbox = TRUE, pretty = TRUE)
