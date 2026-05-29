# ============================================================
#
# Goal:
# Create 7 cleaned cohort files only.
#
# Input:
# ../data/fraga.rds
# ../data/johansson.rds
# ../data/kalyakulina.rds
# ../data/kurushima.rds
# ../data/mahdi.rds
# ../data/quinn.rds
# ../data/xu.rds
#
# Outputs:
# ../results/cleaned_cohorts/fraga_cleaned.rds
# ../results/cleaned_cohorts/johansson_cleaned.rds
# ../results/cleaned_cohorts/kalyakulina_cleaned.rds
# ../results/cleaned_cohorts/kurushima_cleaned.rds
# ../results/cleaned_cohorts/mahdi_cleaned.rds
# ../results/cleaned_cohorts/quinn_cleaned.rds
# ../results/cleaned_cohorts/xu_cleaned.rds
#
# Also saves:
# ../results/common_cpgs.rds
# ../results/qc_summary.csv
#
# Notes:
# - Rows are CpGs / cytosines.
# - Columns are samples.
# - Sample annotations are in colanns(x).
# - CpG annotations are in rowanns(x).
# ============================================================


# ------------------------------------------------------------
# 0. Package
# ------------------------------------------------------------

library(annmatrix)


# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

data_dir <- file.path("..", "data")
results_dir <- file.path("..", "results")
cleaned_dir <- file.path(results_dir, "cleaned_cohorts")

if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

if (!dir.exists(cleaned_dir)) {
  dir.create(cleaned_dir)
}

files <- file.path(
  data_dir,
  c(
    "fraga.rds",
    "johansson.rds",
    "kalyakulina.rds",
    "kurushima.rds",
    "mahdi.rds",
    "quinn.rds",
    "xu.rds"
  )
)

if (!all(file.exists(files))) {
  print(files)
  stop("Some data files were not found. Make sure you run this script from the scripts/ directory.")
}


# ------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------

status <- function(...) {
  message(format(Sys.time(), "%H:%M:%S"), " | ", ...)
}


get_cpg_ids <- function(x) {
  as.character(rowanns(x, "id"))
}


is_true <- function(x) {
  x %in% c(TRUE, "TRUE", "True", "true", "T", "t", 1, "1")
}


check_required_variables <- function(x, cohort) {
  
  needed_sample_vars <- c(
    "id",
    "age",
    "qc_iac",
    "qc_detection",
    "qc_badsex",
    "qc_badsnp"
  )
  
  missing_sample_vars <- setdiff(needed_sample_vars, colnames(colanns(x)))
  
  if (length(missing_sample_vars) > 0) {
    stop(
      paste(
        "Missing sample annotation variables in",
        cohort,
        ":",
        paste(missing_sample_vars, collapse = ", ")
      )
    )
  }
  
  if (!"id" %in% colnames(rowanns(x))) {
    stop(paste("Missing CpG ID variable 'id' in row annotations for", cohort))
  }
}


filter_samples_qc <- function(x) {
  
  ca <- colanns(x)
  
  keep <- rep(TRUE, ncol(x))
  
  # QC rules from the task:
  # qc_iac: remove samples with qc_iac < -3
  # qc_detection: remove samples with qc_detection < 0.95
  # qc_badsex: remove TRUE
  # qc_badsnp: remove TRUE
  #
  # NA means the QC type was not applicable, so we keep NA.
  
  keep <- keep & (is.na(ca$qc_iac) | ca$qc_iac >= -3)
  keep <- keep & (is.na(ca$qc_detection) | ca$qc_detection >= 0.95)
  keep <- keep & (is.na(ca$qc_badsex) | !is_true(ca$qc_badsex))
  keep <- keep & (is.na(ca$qc_badsnp) | !is_true(ca$qc_badsnp))
  
  list(
    data = x[, keep, drop = FALSE],
    samples_before = ncol(x),
    samples_after = sum(keep),
    samples_removed = sum(!keep)
  )
}


# ------------------------------------------------------------
# 3. First pass: find CpGs common to all 7 datasets
# ------------------------------------------------------------

status("FIRST PASS started: finding CpGs common to all datasets.")

common_cpgs <- NULL

for (i in seq_along(files)) {
  
  file <- files[i]
  cohort <- sub("\\.rds$", "", basename(file))
  
  status("[", i, "/", length(files), "] Reading CpG IDs from: ", cohort)
  
  x <- readRDS(file)
  
  check_required_variables(x, cohort)
  
  ids <- get_cpg_ids(x)
  
  if (is.null(common_cpgs)) {
    common_cpgs <- ids
  } else {
    common_cpgs <- intersect(common_cpgs, ids)
  }
  
  status("Current common CpGs after ", cohort, ": ", length(common_cpgs))
  
  rm(x, ids)
  gc()
}

if (length(common_cpgs) == 0) {
  stop("No common CpGs found across datasets.")
}

saveRDS(
  common_cpgs,
  file = file.path(results_dir, "common_cpgs.rds")
)

status("FIRST PASS finished. Final number of common CpGs: ", length(common_cpgs))
status("Saved common CpGs to: ", file.path(results_dir, "common_cpgs.rds"))


# ------------------------------------------------------------
# 4. Second pass: clean each cohort and save separately
# ------------------------------------------------------------

status("SECOND PASS started: cleaning each cohort separately.")

qc_summary <- data.frame(
  cohort = character(),
  samples_before = integer(),
  samples_after_qc = integer(),
  samples_removed_qc = integer(),
  samples_after_age_filter = integer(),
  samples_removed_age_filter = integer(),
  cpgs_common = integer(),
  output_file = character(),
  stringsAsFactors = FALSE
)

for (i in seq_along(files)) {
  
  file <- files[i]
  cohort <- sub("\\.rds$", "", basename(file))
  
  status("[", i, "/", length(files), "] Cleaning cohort: ", cohort)
  
  # ----------------------------------------------------------
  # 4.1 Load
  # ----------------------------------------------------------
  
  status(cohort, ": loading raw file.")
  
  x <- readRDS(file)
  
  check_required_variables(x, cohort)
  
  status(
    cohort, ": loaded. Raw dimensions: ",
    nrow(x), " CpGs x ", ncol(x), " samples."
  )
  
  # Add cohort label to sample annotations
  colanns(x, "cohort") <- cohort
  
  # Make sample names unique across all future work
  colnames(x) <- paste0(cohort, "_", colnames(x))
  
  # ----------------------------------------------------------
  # 4.2 Sample QC
  # ----------------------------------------------------------
  
  status(cohort, ": applying sample QC.")
  
  qc_res <- filter_samples_qc(x)
  x <- qc_res$data
  
  status(
    cohort, ": QC finished. Samples: ",
    qc_res$samples_before, " -> ", qc_res$samples_after,
    " | removed: ", qc_res$samples_removed
  )
  
  # ----------------------------------------------------------
  # 4.3 Keep only common CpGs
  # ----------------------------------------------------------
  
  status(cohort, ": matching and keeping common CpGs.")
  
  ids <- get_cpg_ids(x)
  idx <- match(common_cpgs, ids)
  
  if (any(is.na(idx))) {
    stop(paste("CpG matching failed for cohort:", cohort))
  }
  
  x <- x[idx, , drop = FALSE]
  
  # Make rownames consistent
  rownames(x) <- common_cpgs
  
  status(cohort, ": common CpG filter finished. CpGs kept: ", nrow(x))
  
  # ----------------------------------------------------------
  # 4.4 Remove samples without age
  # ----------------------------------------------------------
  
  status(cohort, ": checking age values.")
  
  age_numeric <- as.numeric(colanns(x, "age"))
  keep_age <- !is.na(age_numeric)
  
  samples_before_age_filter <- ncol(x)
  
  x <- x[, keep_age, drop = FALSE]
  
  samples_after_age_filter <- ncol(x)
  samples_removed_age_filter <- samples_before_age_filter - samples_after_age_filter
  
  # Recalculate after filtering
  age_numeric <- as.numeric(colanns(x, "age"))
  
  # Keep original numeric age and transformed modelling age.
  # We use log1p(age) = log(age + 1), because some cohorts may contain age 0.
  # Later predictions can be transformed back with expm1(prediction).
  colanns(x, "age_numeric") <- age_numeric
  colanns(x, "age_for_model") <- log1p(age_numeric)
  
  status(
    cohort, ": age filter finished. Samples: ",
    samples_before_age_filter, " -> ", samples_after_age_filter,
    " | removed: ", samples_removed_age_filter
  )
  
  # ----------------------------------------------------------
  # 4.5 Save cleaned cohort
  # ----------------------------------------------------------
  
  output_file <- file.path(cleaned_dir, paste0(cohort, "_cleaned.rds"))
  
  status(cohort, ": saving cleaned file.")
  
  saveRDS(x, output_file)
  
  status(cohort, ": saved to ", output_file)
  
  # ----------------------------------------------------------
  # 4.6 Update summary
  # ----------------------------------------------------------
  
  qc_summary <- rbind(
    qc_summary,
    data.frame(
      cohort = cohort,
      samples_before = qc_res$samples_before,
      samples_after_qc = qc_res$samples_after,
      samples_removed_qc = qc_res$samples_removed,
      samples_after_age_filter = samples_after_age_filter,
      samples_removed_age_filter = samples_removed_age_filter,
      cpgs_common = nrow(x),
      output_file = output_file,
      stringsAsFactors = FALSE
    )
  )
  
  rm(x, ids, idx, age_numeric, keep_age)
  gc()
  
  status("[", i, "/", length(files), "] Finished cohort: ", cohort)
}


# ------------------------------------------------------------
# 5. Save QC summary
# ------------------------------------------------------------

write.csv(
  qc_summary,
  file = file.path(results_dir, "qc_summary.csv"),
  row.names = FALSE
)

status("Cleaning finished.")
status("Cleaned cohort files saved in: ", cleaned_dir)
status("QC summary saved to: ", file.path(results_dir, "qc_summary.csv"))

print(qc_summary)