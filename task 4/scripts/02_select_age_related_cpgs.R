# ============================================================
# 02_select_age_related_cpgs.R
#
# Run from the scripts/ directory:
# source("02_select_age_related_cpgs.R")
#
# Goal:
# Select CpGs most associated with age.
#
# Input:
# ../results/cleaned_cohorts/*_cleaned.rds
#
# Outputs:
# ../results/age_cpg_selection/age_cpg_scores.rds
# ../results/age_cpg_selection/top_1000_age_cpgs.rds
# ../results/age_cpg_selection/top_5000_age_cpgs.rds
# ../results/age_cpg_selection/top_10000_age_cpgs.rds
# ../results/age_cpg_selection/top_age_cpg_scores.csv
# ============================================================


# ------------------------------------------------------------
# 0. Package
# ------------------------------------------------------------

library(annmatrix)


# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

cleaned_dir <- file.path("..", "results", "cleaned_cohorts")
selection_dir <- file.path("..", "results", "age_cpg_selection")

if (!dir.exists(selection_dir)) {
  dir.create(selection_dir)
}

cleaned_files <- file.path(
  cleaned_dir,
  c(
    "fraga_cleaned.rds",
    "johansson_cleaned.rds",
    "kalyakulina_cleaned.rds",
    "kurushima_cleaned.rds",
    "mahdi_cleaned.rds",
    "quinn_cleaned.rds",
    "xu_cleaned.rds"
  )
)

if (!all(file.exists(cleaned_files))) {
  print(cleaned_files)
  stop("Some cleaned cohort files were not found.")
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


safe_cor_with_age <- function(mat, age) {
  
  # mat:
  # rows = CpGs
  # columns = samples
  #
  # age:
  # numeric vector, length = number of samples
  
  age <- as.numeric(age)
  
  if (any(is.na(age))) {
    stop("Age vector contains NA values.")
  }
  
  # Center age
  age_centered <- age - mean(age)
  age_sd <- sqrt(sum(age_centered^2))
  
  if (age_sd == 0) {
    stop("Age has zero variability in this cohort.")
  }
  
  # CpG means
  cpg_means <- rowMeans(mat, na.rm = TRUE)
  
  # Center CpGs
  mat_centered <- mat - cpg_means
  
  # Numerator of correlation
  numerator <- as.vector(mat_centered %*% age_centered)
  
  # Denominator
  cpg_sd <- sqrt(rowSums(mat_centered^2, na.rm = TRUE))
  denominator <- cpg_sd * age_sd
  
  cor_values <- numerator / denominator
  
  # Handle CpGs with zero variance or numerical issues
  cor_values[!is.finite(cor_values)] <- NA
  
  cor_values
}


# ------------------------------------------------------------
# 3. Main selection settings
# ------------------------------------------------------------

top_n_values <- c(1000, 5000, 10000)

# If memory becomes a problem, set this lower.
# Usually 369543 CpGs x samples is okay file-by-file.
chunk_size <- 50000


# ------------------------------------------------------------
# 4. Process cohorts one by one
# ------------------------------------------------------------

status("Starting age-related CpG selection.")

all_scores <- NULL
cohort_summary <- data.frame(
  cohort = character(),
  samples = integer(),
  cpgs = integer(),
  age_min = numeric(),
  age_median = numeric(),
  age_max = numeric(),
  stringsAsFactors = FALSE
)

for (i in seq_along(cleaned_files)) {
  
  file <- cleaned_files[i]
  cohort <- sub("_cleaned\\.rds$", "", basename(file))
  
  status("[", i, "/", length(cleaned_files), "] Loading cleaned cohort: ", cohort)
  
  x <- readRDS(file)
  
  if (!"age_for_model" %in% colnames(colanns(x))) {
    stop(paste("age_for_model not found in:", cohort))
  }
  
  if (!"age_numeric" %in% colnames(colanns(x))) {
    stop(paste("age_numeric not found in:", cohort))
  }
  
  cpg_ids <- get_cpg_ids(x)
  age <- as.numeric(colanns(x, "age_for_model"))
  age_original <- as.numeric(colanns(x, "age_numeric"))
  
  status(
    cohort, ": dimensions = ",
    nrow(x), " CpGs x ", ncol(x), " samples."
  )
  
  status(
    cohort, ": age range = ",
    min(age_original), " to ", max(age_original), " years."
  )
  
  mat <- as.matrix(x)
  
  # ----------------------------------------------------------
  # Calculate correlations in chunks
  # ----------------------------------------------------------
  
  status(cohort, ": calculating CpG-age correlations.")
  
  cor_values <- numeric(nrow(mat))
  cor_values[] <- NA_real_
  
  starts <- seq(1, nrow(mat), by = chunk_size)
  
  for (s in starts) {
    
    e <- min(s + chunk_size - 1, nrow(mat))
    
    chunk_idx <- s:e
    
    cor_values[chunk_idx] <- safe_cor_with_age(
      mat = mat[chunk_idx, , drop = FALSE],
      age = age
    )
    
    status(
      cohort, ": correlation progress ",
      e, "/", nrow(mat), " CpGs."
    )
  }
  
  cohort_scores <- data.frame(
    cpg_id = cpg_ids,
    cohort = cohort,
    correlation = cor_values,
    abs_correlation = abs(cor_values),
    stringsAsFactors = FALSE
  )
  
  if (is.null(all_scores)) {
    all_scores <- cohort_scores
  } else {
    all_scores <- rbind(all_scores, cohort_scores)
  }
  
  cohort_summary <- rbind(
    cohort_summary,
    data.frame(
      cohort = cohort,
      samples = ncol(x),
      cpgs = nrow(x),
      age_min = min(age_original),
      age_median = median(age_original),
      age_max = max(age_original),
      stringsAsFactors = FALSE
    )
  )
  
  rm(x, mat, cpg_ids, age, age_original, cor_values, cohort_scores)
  gc()
  
  status("[", i, "/", length(cleaned_files), "] Finished cohort: ", cohort)
}


# ------------------------------------------------------------
# 5. Combine cohort-level scores
# ------------------------------------------------------------

status("Combining CpG scores across cohorts.")

combined_scores <- aggregate(
  cbind(abs_correlation, correlation) ~ cpg_id,
  data = all_scores,
  FUN = function(z) mean(z, na.rm = TRUE)
)

colnames(combined_scores) <- c(
  "cpg_id",
  "mean_abs_correlation",
  "mean_correlation"
)

# Count in how many cohorts the CpG had a finite correlation
finite_counts <- aggregate(
  is.finite(all_scores$correlation),
  by = list(cpg_id = all_scores$cpg_id),
  FUN = sum
)

colnames(finite_counts) <- c("cpg_id", "n_cohorts_finite")

combined_scores <- merge(
  combined_scores,
  finite_counts,
  by = "cpg_id",
  all.x = TRUE
)

combined_scores <- combined_scores[order(-combined_scores$mean_abs_correlation), ]

status("Combined scores created.")


# ------------------------------------------------------------
# 6. Save results
# ------------------------------------------------------------

saveRDS(
  combined_scores,
  file = file.path(selection_dir, "age_cpg_scores.rds")
)

write.csv(
  head(combined_scores, 50000),
  file = file.path(selection_dir, "top_age_cpg_scores.csv"),
  row.names = FALSE
)

write.csv(
  cohort_summary,
  file = file.path(selection_dir, "cohort_age_summary.csv"),
  row.names = FALSE
)

for (top_n in top_n_values) {
  
  selected_cpgs <- head(combined_scores$cpg_id, top_n)
  
  saveRDS(
    selected_cpgs,
    file = file.path(selection_dir, paste0("top_", top_n, "_age_cpgs.rds"))
  )
  
  status("Saved top ", top_n, " CpGs.")
}


# ------------------------------------------------------------
# 7. Final summary
# ------------------------------------------------------------

status("Age-related CpG selection finished.")
status("Results saved in: ", selection_dir)

status("Top 10 CpGs by mean absolute correlation:")
print(head(combined_scores, 10))

status("Cohort age summary:")
print(cohort_summary)