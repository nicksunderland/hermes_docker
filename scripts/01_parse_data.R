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
  gwas_file     = '/Users/xx20081/git/hermes_docker/tests/data/test_data.tsv.gz'
  meta_file     = '/Users/xx20081/git/hermes_docker/tests/data/config.json'
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

# per column functions that return true if that row is valid
col_fn_list <- list(
  rsid      = function(x) grepl("^rs[0-9]+|^[0-9]+:[0-9]+.*", x),
  chr       = function(x) x %in% c(as.character(1:22), "X", "Y"),
  pos       = function(x) is.integer(x) & x > 0,
  alt       = function(x) grepl("^[ACTG]+$", x),
  reference = function(x) grepl("^[ACTG]+$", x),
  eaf       = function(x) is.numeric(x) & x > 0 & x < 1,
  beta      = function(x) is.numeric(x) & is.finite(x) & abs(x) < 20,
  stdErr    = function(x) is.numeric(x) & x > 0 & is.finite(x),
  pval      = function(x) is.numeric(x) & x >= 0 & x <=1,
  n         = function(x) is.integer(x) & x > 0,
  ncase     = function(x) is.integer(x) & x > 0,
  ncontrol  = function(x) is.integer(x) & x > 0,
  imputed   = function(x) is.logical(x) & !is.na(x),
  info      = function(x) is.numeric(x) & x >= 0 & x <=1
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

# find indels and and number to summary, remove indels if flag specified
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
  build         = fcase(grepl("hg19|b37", meta$referenceGenome, ignore.case=T), "GRCh37",
                        grepl("hg38|b38", meta$referenceGenome, ignore.case=T), "GRCh38",
                        default = NA_character_)
)

# write JSON to disk
jsonlite::write_json(json_data, path = gwas2vcf_json, auto_unbox = TRUE, pretty = TRUE)
