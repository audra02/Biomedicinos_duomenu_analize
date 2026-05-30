# ============================================================
# Audra.R
#
# Run from the task4/scripts/ directory:
# source("Audra.R")
#
# Goal:
# Create and compare two epigenetic age clocks.
# These are not existing epigenetic clocks. They are trained from the
# provided methylation cohorts only.
#
# Clock 1: m_value_ridge_calibrated
# - beta methylation values are transformed to M-values
# - features are scaled using training data only
# - a closed-form ridge regression is fitted on log1p(age)
# - final predictions are corrected using a train-only quadratic calibration
#
# Clock 2: signature_kernel
# - creates a signed age-methylation signature from age-related CpGs
# - increasing and decreasing CpGs are combined into one directional score
# - age is predicted by a one-dimensional kernel smoother
# - the kernel prediction is blended with a cubic signature curve
# - final predictions are corrected using a train-only quadratic calibration
#
# Input:
# ../results/cleaned_cohorts/*_cleaned.rds
# ../results/age_cpg_selection/age_cpg_scores.rds
#
# Outputs:
# ../results/model_comparison/model_predictions_1.csv
# ../results/model_comparison/model_metrics_1.csv
# ../results/model_comparison/trained_models_1.rds
# ../results/model_comparison/ridge_tuning_results_1.csv
# ../results/model_comparison/signature_tuning_results_1.csv
# ../results/model_comparison/prediction_plot_*_1.png
# ============================================================


# ------------------------------------------------------------
# 0. Package
# ------------------------------------------------------------

library(annmatrix)


# ------------------------------------------------------------
# 1. Paths and settings
# ------------------------------------------------------------

cleaned_dir <- file.path("..", "results", "cleaned_cohorts")
selection_dir <- file.path("..", "results", "age_cpg_selection")
model_dir <- file.path("..", "results", "model_comparison")

if (!dir.exists(model_dir)) {
  dir.create(model_dir, recursive = TRUE)
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
  stop("Some cleaned cohort files were not found. Run prepare_data.R first.")
}

score_file <- file.path(selection_dir, "age_cpg_scores.rds")

if (!file.exists(score_file)) {
  stop("Could not find age_cpg_scores.rds. Run 02_select_age_related_cpgs.R first.")
}

set.seed(123)

train_fraction <- 0.80
internal_train_fraction <- 0.80

# The maximum number of CpGs loaded into memory. Both clocks tune within this set.
max_loaded_cpgs <- 10000

# Ridge clock tuning.
ridge_top_n_candidates <- c(500, 1000, 2000, 5000, 10000)
ridge_lambda_candidates <- c(0.1, 1, 10, 100, 1000)

# Signature clock tuning.
signature_top_n_candidates <- c(100, 250, 500, 1000, 2000, 5000)
signature_k_candidates <- c(15, 25, 50, 75, 100)
signature_blend_candidates <- c(0.00, 0.25, 0.50, 0.75, 1.00)

# Small value used when beta values are converted to M-values.
m_value_eps <- 1e-5


# ------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------

status <- function(...) {
  message(format(Sys.time(), "%H:%M:%S"), " | ", ...)
}


get_cpg_ids <- function(x) {
  as.character(rowanns(x, "id"))
}


inverse_age_transform <- function(predicted_log_age) {
  pmax(expm1(predicted_log_age), 0)
}


cap_age_predictions <- function(predicted_age, training_age) {
  upper <- max(training_age) + 5
  pmin(pmax(predicted_age, 0), upper)
}


calculate_metrics <- function(actual_age, predicted_age) {
  ok <- is.finite(actual_age) & is.finite(predicted_age)
  actual_age <- actual_age[ok]
  predicted_age <- predicted_age[ok]
  
  mae <- mean(abs(actual_age - predicted_age))
  rmse <- sqrt(mean((actual_age - predicted_age)^2))
  
  if (length(unique(actual_age)) > 1 && length(unique(predicted_age)) > 1) {
    cor_value <- cor(actual_age, predicted_age)
  } else {
    cor_value <- NA_real_
  }
  
  data.frame(
    n = length(actual_age),
    MAE = mae,
    RMSE = rmse,
    correlation = cor_value
  )
}


make_train_test_split <- function(metadata, train_fraction = 0.80) {
  train_index <- logical(nrow(metadata))
  
  for (cohort in unique(metadata$cohort)) {
    idx <- which(metadata$cohort == cohort)
    n_train <- floor(length(idx) * train_fraction)
    train_index[sample(idx, n_train)] <- TRUE
  }
  
  list(
    train = which(train_index),
    test = which(!train_index)
  )
}


evaluate_prediction_groups <- function(pred_df, model_name) {
  out <- data.frame()
  
  groups <- list(
    all_test = rep(TRUE, nrow(pred_df)),
    non_quinn_test = pred_df$cohort != "quinn",
    quinn_test = pred_df$cohort == "quinn"
  )
  
  for (group_name in names(groups)) {
    idx <- groups[[group_name]]
    
    if (sum(idx) == 0) {
      next
    }
    
    m <- calculate_metrics(
      actual_age = pred_df$actual_age[idx],
      predicted_age = pred_df$predicted_age[idx]
    )
    
    m$model <- model_name
    m$test_group <- group_name
    out <- rbind(out, m)
  }
  
  out[, c("model", "test_group", "n", "MAE", "RMSE", "correlation")]
}


plot_predictions_base <- function(pred_df, model_name, output_file) {
  png(output_file, width = 900, height = 700)
  
  plot(
    pred_df$actual_age,
    pred_df$predicted_age,
    xlab = "Actual age",
    ylab = "Predicted age",
    main = model_name,
    pch = 19,
    cex = 0.7
  )
  
  abline(0, 1, lwd = 2)
  dev.off()
}


looks_like_beta_values <- function(x) {
  vals <- x[is.finite(x)]
  
  if (length(vals) == 0) {
    return(FALSE)
  }
  
  min(vals) >= -0.01 && max(vals) <= 1.01
}


beta_to_m_values <- function(x, eps = 1e-5) {
  x <- pmin(pmax(x, eps), 1 - eps)
  log2(x / (1 - x))
}


impute_nonfinite_by_train_mean <- function(x_train, x_test) {
  feature_means <- colMeans(x_train, na.rm = TRUE)
  feature_means[!is.finite(feature_means)] <- 0
  
  bad_train <- which(!is.finite(x_train), arr.ind = TRUE)
  if (nrow(bad_train) > 0) {
    x_train[bad_train] <- feature_means[bad_train[, 2]]
  }
  
  bad_test <- which(!is.finite(x_test), arr.ind = TRUE)
  if (nrow(bad_test) > 0) {
    x_test[bad_test] <- feature_means[bad_test[, 2]]
  }
  
  list(
    x_train = x_train,
    x_test = x_test,
    imputation_means = feature_means
  )
}


preprocess_ridge_train_test <- function(x_train_raw, x_test_raw) {
  use_m_values <- looks_like_beta_values(x_train_raw)
  
  if (use_m_values) {
    x_train <- beta_to_m_values(x_train_raw, eps = m_value_eps)
    x_test <- beta_to_m_values(x_test_raw, eps = m_value_eps)
  } else {
    x_train <- x_train_raw
    x_test <- x_test_raw
  }
  
  imputed <- impute_nonfinite_by_train_mean(x_train, x_test)
  x_train <- imputed$x_train
  x_test <- imputed$x_test
  
  feature_means <- colMeans(x_train)
  feature_sds <- apply(x_train, 2, sd)
  feature_sds[!is.finite(feature_sds) | feature_sds == 0] <- 1
  
  x_train_scaled <- sweep(x_train, 2, feature_means, "-")
  x_train_scaled <- sweep(x_train_scaled, 2, feature_sds, "/")
  
  x_test_scaled <- sweep(x_test, 2, feature_means, "-")
  x_test_scaled <- sweep(x_test_scaled, 2, feature_sds, "/")
  
  list(
    x_train = x_train_scaled,
    x_test = x_test_scaled,
    use_m_values = use_m_values,
    imputation_means = imputed$imputation_means,
    feature_means = feature_means,
    feature_sds = feature_sds
  )
}


preprocess_signature_train_test <- function(x_train_raw, x_test_raw) {
  # The signature clock intentionally uses robust scaling on beta values.
  # This makes it different from the M-value ridge clock.
  x_train <- x_train_raw
  x_test <- x_test_raw
  
  imputed <- impute_nonfinite_by_train_mean(x_train, x_test)
  x_train <- imputed$x_train
  x_test <- imputed$x_test
  
  feature_centers <- apply(x_train, 2, median)
  feature_scales <- apply(x_train, 2, IQR)
  fallback_sds <- apply(x_train, 2, sd)
  
  use_fallback <- !is.finite(feature_scales) | feature_scales == 0
  feature_scales[use_fallback] <- fallback_sds[use_fallback]
  feature_scales[!is.finite(feature_scales) | feature_scales == 0] <- 1
  
  x_train_scaled <- sweep(x_train, 2, feature_centers, "-")
  x_train_scaled <- sweep(x_train_scaled, 2, feature_scales, "/")
  
  x_test_scaled <- sweep(x_test, 2, feature_centers, "-")
  x_test_scaled <- sweep(x_test_scaled, 2, feature_scales, "/")
  
  list(
    x_train = x_train_scaled,
    x_test = x_test_scaled,
    imputation_means = imputed$imputation_means,
    feature_centers = feature_centers,
    feature_scales = feature_scales
  )
}


fit_ridge_dual <- function(x_train, y_train, lambda) {
  y_mean <- mean(y_train)
  y_centered <- y_train - y_mean
  
  k_mat <- tcrossprod(x_train)
  diag(k_mat) <- diag(k_mat) + lambda
  
  alpha <- solve(k_mat, y_centered)
  beta <- as.vector(crossprod(x_train, alpha))
  
  list(
    intercept = y_mean,
    beta = beta,
    lambda = lambda
  )
}


predict_ridge_dual <- function(model, x_new) {
  as.numeric(model$intercept + x_new %*% model$beta)
}


fit_age_calibration <- function(predicted_age, actual_age) {
  df <- data.frame(
    actual_age = actual_age,
    predicted_age = predicted_age
  )
  
  df <- df[is.finite(df$actual_age) & is.finite(df$predicted_age), ]
  
  if (nrow(df) < 20 || length(unique(df$predicted_age)) < 5) {
    return(NULL)
  }
  
  lm(actual_age ~ predicted_age + I(predicted_age^2), data = df)
}


apply_age_calibration <- function(calibration_model, predicted_age, training_age) {
  if (is.null(calibration_model)) {
    return(cap_age_predictions(predicted_age, training_age))
  }
  
  df <- data.frame(predicted_age = predicted_age)
  calibrated <- as.numeric(predict(calibration_model, newdata = df))
  cap_age_predictions(calibrated, training_age)
}


make_signature_weights <- function(cpg_info) {
  corr <- cpg_info$mean_correlation
  corr[!is.finite(corr)] <- 0
  
  weights <- sign(corr) * sqrt(abs(corr))
  
  if (sum(abs(weights)) == 0) {
    weights[] <- 1
  }
  
  weights / sum(abs(weights))
}


calculate_signature_score <- function(x_scaled, weights) {
  as.numeric(x_scaled %*% weights)
}


predict_signature_kernel <- function(train_score, train_y, test_score, k) {
  pred <- numeric(length(test_score))
  k <- min(k, length(train_score))
  
  for (i in seq_along(test_score)) {
    distances <- abs(train_score - test_score[i])
    neighbor_idx <- order(distances)[seq_len(k)]
    neighbor_dist <- distances[neighbor_idx]
    
    bandwidth <- max(neighbor_dist)
    if (!is.finite(bandwidth) || bandwidth == 0) {
      bandwidth <- 1e-8
    }
    
    weights <- exp(-0.5 * (neighbor_dist / bandwidth)^2)
    weights[!is.finite(weights)] <- 0
    
    if (sum(weights) == 0) {
      weights[] <- 1
    }
    
    pred[i] <- sum(weights * train_y[neighbor_idx]) / sum(weights)
  }
  
  pred
}


fit_signature_curve <- function(train_score, train_y) {
  df <- data.frame(
    y = train_y,
    score = train_score
  )
  
  lm(y ~ score + I(score^2) + I(score^3), data = df)
}


predict_signature_curve <- function(model, new_score) {
  as.numeric(predict(model, newdata = data.frame(score = new_score)))
}


# ------------------------------------------------------------
# 3. Load CpG score table
# ------------------------------------------------------------

status("Loading age-related CpG scores.")

cpg_scores <- readRDS(score_file)

needed_score_cols <- c("cpg_id", "mean_abs_correlation", "mean_correlation")

if (!all(needed_score_cols %in% colnames(cpg_scores))) {
  stop("age_cpg_scores.rds does not have the expected columns.")
}

cpg_scores <- cpg_scores[order(-cpg_scores$mean_abs_correlation), ]
cpg_scores <- cpg_scores[is.finite(cpg_scores$mean_abs_correlation), ]

if (nrow(cpg_scores) < max_loaded_cpgs) {
  max_loaded_cpgs <- nrow(cpg_scores)
}

selected_cpgs <- head(cpg_scores$cpg_id, max_loaded_cpgs)
selected_cpg_info <- cpg_scores[match(selected_cpgs, cpg_scores$cpg_id), ]

status("CpGs loaded for clocks: ", length(selected_cpgs))

ridge_top_n_candidates <- ridge_top_n_candidates[ridge_top_n_candidates <= length(selected_cpgs)]
signature_top_n_candidates <- signature_top_n_candidates[signature_top_n_candidates <= length(selected_cpgs)]


# ------------------------------------------------------------
# 4. Load cleaned cohorts and extract selected CpGs
# ------------------------------------------------------------

status("Loading cleaned cohorts and extracting selected CpGs.")

x_list <- list()
metadata_list <- list()

for (i in seq_along(cleaned_files)) {
  file <- cleaned_files[i]
  cohort <- sub("_cleaned\\.rds$", "", basename(file))
  
  status("[", i, "/", length(cleaned_files), "] Loading: ", cohort)
  
  x <- readRDS(file)
  cpg_ids <- get_cpg_ids(x)
  idx <- match(selected_cpgs, cpg_ids)
  
  if (any(is.na(idx))) {
    stop("Some selected CpGs were not found in cohort: ", cohort)
  }
  
  x_sub <- x[idx, , drop = FALSE]
  mat <- t(as.matrix(x_sub))
  colnames(mat) <- selected_cpgs
  
  meta <- as.data.frame(colanns(x_sub))
  meta$sample_id <- colnames(x_sub)
  
  needed_cols <- c("sample_id", "cohort", "age_numeric", "age_for_model")
  
  if (!all(needed_cols %in% colnames(meta))) {
    stop("Missing required metadata columns in cohort: ", cohort)
  }
  
  meta <- meta[, needed_cols, drop = FALSE]
  meta$age_numeric <- as.numeric(as.character(meta$age_numeric))
  meta$age_for_model <- as.numeric(as.character(meta$age_for_model))
  
  if (any(is.na(meta$age_numeric))) {
    stop("NA values found in age_numeric in cohort: ", cohort)
  }
  
  if (any(is.na(meta$age_for_model))) {
    stop("NA values found in age_for_model in cohort: ", cohort)
  }
  
  x_list[[cohort]] <- mat
  metadata_list[[cohort]] <- meta
  
  status(
    cohort, ": extracted matrix = ",
    nrow(mat), " samples x ", ncol(mat), " CpGs."
  )
  
  rm(x, x_sub, mat, meta, cpg_ids, idx)
  gc()
}

X <- do.call(rbind, x_list)
metadata <- do.call(rbind, metadata_list)
rownames(metadata) <- metadata$sample_id
metadata <- metadata[rownames(X), , drop = FALSE]

y <- as.numeric(metadata$age_for_model)
actual_age <- as.numeric(metadata$age_numeric)

if (any(is.na(y)) || any(is.na(actual_age))) {
  stop("Age values contain NA after data loading.")
}

if (length(unique(y)) <= 1) {
  stop("Model target is constant. Check age annotations.")
}

status("Final modelling matrix:")
status("Samples: ", nrow(X))
status("CpGs: ", ncol(X))
status("Age range: ", min(actual_age), " to ", max(actual_age), " years.")
status("Samples per cohort:")
print(table(metadata$cohort))


# ------------------------------------------------------------
# 5. Train/test split
# ------------------------------------------------------------

status("Creating stratified train/test split by cohort.")

split <- make_train_test_split(metadata, train_fraction = train_fraction)

train_idx <- split$train
test_idx <- split$test

status("Training samples: ", length(train_idx))
status("Testing samples: ", length(test_idx))
status("Test samples per cohort:")
print(table(metadata$cohort[test_idx]))

X_train_all_raw <- X[train_idx, , drop = FALSE]
X_test_all_raw <- X[test_idx, , drop = FALSE]

y_train <- y[train_idx]
y_test <- y[test_idx]
actual_train_age <- actual_age[train_idx]

test_metadata <- metadata[test_idx, , drop = FALSE]

internal_meta <- metadata[train_idx, , drop = FALSE]
internal_split <- make_train_test_split(
  metadata = internal_meta,
  train_fraction = internal_train_fraction
)

subtrain_local <- internal_split$train
valid_local <- internal_split$test

y_subtrain <- y_train[subtrain_local]
y_valid <- y_train[valid_local]
actual_subtrain_age <- actual_train_age[subtrain_local]
actual_valid_age <- actual_train_age[valid_local]

status("Internal tuning split:")
status("Subtrain samples: ", length(subtrain_local))
status("Validation samples: ", length(valid_local))


# ------------------------------------------------------------
# 6. Clock 1 tuning: M-value ridge with calibration
# ------------------------------------------------------------

status("CLOCK 1 tuning started: m_value_ridge_calibrated.")

ridge_tuning_results <- data.frame()

for (top_n in ridge_top_n_candidates) {
  status("CLOCK 1: tuning top_n = ", top_n)
  
  feature_idx <- seq_len(top_n)
  
  prep <- preprocess_ridge_train_test(
    x_train_raw = X_train_all_raw[subtrain_local, feature_idx, drop = FALSE],
    x_test_raw = X_train_all_raw[valid_local, feature_idx, drop = FALSE]
  )
  
  for (lambda in ridge_lambda_candidates) {
    ridge_model_tmp <- fit_ridge_dual(
      x_train = prep$x_train,
      y_train = y_subtrain,
      lambda = lambda
    )
    
    pred_subtrain_log <- predict_ridge_dual(ridge_model_tmp, prep$x_train)
    pred_valid_log <- predict_ridge_dual(ridge_model_tmp, prep$x_test)
    
    pred_subtrain_age <- inverse_age_transform(pred_subtrain_log)
    pred_valid_age_raw <- inverse_age_transform(pred_valid_log)
    
    calibration_tmp <- fit_age_calibration(
      predicted_age = pred_subtrain_age,
      actual_age = actual_subtrain_age
    )
    
    pred_valid_age <- apply_age_calibration(
      calibration_model = calibration_tmp,
      predicted_age = pred_valid_age_raw,
      training_age = actual_subtrain_age
    )
    
    metrics <- calculate_metrics(actual_valid_age, pred_valid_age)
    
    ridge_tuning_results <- rbind(
      ridge_tuning_results,
      data.frame(
        top_n = top_n,
        lambda = lambda,
        validation_MAE = metrics$MAE,
        validation_RMSE = metrics$RMSE,
        validation_correlation = metrics$correlation,
        stringsAsFactors = FALSE
      )
    )
    
    rm(ridge_model_tmp, pred_subtrain_log, pred_valid_log)
    gc()
  }
  
  rm(prep)
  gc()
}

ridge_tuning_results <- ridge_tuning_results[order(ridge_tuning_results$validation_MAE), ]

write.csv(
  ridge_tuning_results,
  file = file.path(model_dir, "ridge_tuning_results_1.csv"),
  row.names = FALSE
)

best_ridge_top_n <- ridge_tuning_results$top_n[1]
best_ridge_lambda <- ridge_tuning_results$lambda[1]

status(
  "CLOCK 1 best settings: top_n = ", best_ridge_top_n,
  ", lambda = ", best_ridge_lambda
)


# ------------------------------------------------------------
# 7. Clock 1 final fit and prediction
# ------------------------------------------------------------

status("CLOCK 1 final fitting started.")

ridge_feature_idx <- seq_len(best_ridge_top_n)
ridge_cpgs <- selected_cpgs[ridge_feature_idx]

ridge_prep_final <- preprocess_ridge_train_test(
  x_train_raw = X_train_all_raw[, ridge_feature_idx, drop = FALSE],
  x_test_raw = X_test_all_raw[, ridge_feature_idx, drop = FALSE]
)

ridge_model_final <- fit_ridge_dual(
  x_train = ridge_prep_final$x_train,
  y_train = y_train,
  lambda = best_ridge_lambda
)

ridge_train_pred_log <- predict_ridge_dual(
  ridge_model_final,
  ridge_prep_final$x_train
)

ridge_test_pred_log <- predict_ridge_dual(
  ridge_model_final,
  ridge_prep_final$x_test
)

ridge_train_pred_age <- inverse_age_transform(ridge_train_pred_log)
ridge_test_pred_age_raw <- inverse_age_transform(ridge_test_pred_log)

ridge_calibration_model <- fit_age_calibration(
  predicted_age = ridge_train_pred_age,
  actual_age = actual_train_age
)

ridge_test_pred_age <- apply_age_calibration(
  calibration_model = ridge_calibration_model,
  predicted_age = ridge_test_pred_age_raw,
  training_age = actual_train_age
)

ridge_predictions <- data.frame(
  sample_id = rownames(X_test_all_raw),
  cohort = test_metadata$cohort,
  actual_age = as.numeric(test_metadata$age_numeric),
  predicted_age = ridge_test_pred_age,
  predicted_log_age = ridge_test_pred_log,
  model = "m_value_ridge_calibrated",
  stringsAsFactors = FALSE
)

status("CLOCK 1 finished: m_value_ridge_calibrated.")


# ------------------------------------------------------------
# 8. Clock 2 tuning: signed signature kernel
# ------------------------------------------------------------

status("CLOCK 2 tuning started: signature_kernel.")

signature_tuning_results <- data.frame()

for (top_n in signature_top_n_candidates) {
  status("CLOCK 2: tuning top_n = ", top_n)
  
  feature_idx <- seq_len(top_n)
  cpg_info_tmp <- selected_cpg_info[feature_idx, , drop = FALSE]
  weights_tmp <- make_signature_weights(cpg_info_tmp)
  
  prep <- preprocess_signature_train_test(
    x_train_raw = X_train_all_raw[subtrain_local, feature_idx, drop = FALSE],
    x_test_raw = X_train_all_raw[valid_local, feature_idx, drop = FALSE]
  )
  
  score_subtrain <- calculate_signature_score(prep$x_train, weights_tmp)
  score_valid <- calculate_signature_score(prep$x_test, weights_tmp)
  
  curve_model_tmp <- fit_signature_curve(score_subtrain, y_subtrain)
  pred_curve_valid_log <- predict_signature_curve(curve_model_tmp, score_valid)
  
  for (k in signature_k_candidates) {
    if (k >= length(score_subtrain)) {
      next
    }
    
    pred_kernel_valid_log <- predict_signature_kernel(
      train_score = score_subtrain,
      train_y = y_subtrain,
      test_score = score_valid,
      k = k
    )
    
    for (blend in signature_blend_candidates) {
      pred_valid_log <- blend * pred_kernel_valid_log +
        (1 - blend) * pred_curve_valid_log
      
      pred_subtrain_kernel_log <- predict_signature_kernel(
        train_score = score_subtrain,
        train_y = y_subtrain,
        test_score = score_subtrain,
        k = k
      )
      
      pred_curve_subtrain_log <- predict_signature_curve(
        curve_model_tmp,
        score_subtrain
      )
      
      pred_subtrain_log <- blend * pred_subtrain_kernel_log +
        (1 - blend) * pred_curve_subtrain_log
      
      pred_subtrain_age <- inverse_age_transform(pred_subtrain_log)
      pred_valid_age_raw <- inverse_age_transform(pred_valid_log)
      
      calibration_tmp <- fit_age_calibration(
        predicted_age = pred_subtrain_age,
        actual_age = actual_subtrain_age
      )
      
      pred_valid_age <- apply_age_calibration(
        calibration_model = calibration_tmp,
        predicted_age = pred_valid_age_raw,
        training_age = actual_subtrain_age
      )
      
      metrics <- calculate_metrics(actual_valid_age, pred_valid_age)
      
      signature_tuning_results <- rbind(
        signature_tuning_results,
        data.frame(
          top_n = top_n,
          k = k,
          blend = blend,
          validation_MAE = metrics$MAE,
          validation_RMSE = metrics$RMSE,
          validation_correlation = metrics$correlation,
          stringsAsFactors = FALSE
        )
      )
    }
    
    rm(pred_kernel_valid_log)
    gc()
  }
  
  rm(prep, weights_tmp, score_subtrain, score_valid, curve_model_tmp)
  gc()
}

signature_tuning_results <- signature_tuning_results[order(signature_tuning_results$validation_MAE), ]

write.csv(
  signature_tuning_results,
  file = file.path(model_dir, "signature_tuning_results_1.csv"),
  row.names = FALSE
)

best_signature_top_n <- signature_tuning_results$top_n[1]
best_signature_k <- signature_tuning_results$k[1]
best_signature_blend <- signature_tuning_results$blend[1]

status(
  "CLOCK 2 best settings: top_n = ", best_signature_top_n,
  ", k = ", best_signature_k,
  ", blend = ", best_signature_blend
)


# ------------------------------------------------------------
# 9. Clock 2 final fit and prediction
# ------------------------------------------------------------

status("CLOCK 2 final fitting started.")

signature_feature_idx <- seq_len(best_signature_top_n)
signature_cpgs <- selected_cpgs[signature_feature_idx]
signature_cpg_info <- selected_cpg_info[signature_feature_idx, , drop = FALSE]
signature_weights <- make_signature_weights(signature_cpg_info)

signature_prep_final <- preprocess_signature_train_test(
  x_train_raw = X_train_all_raw[, signature_feature_idx, drop = FALSE],
  x_test_raw = X_test_all_raw[, signature_feature_idx, drop = FALSE]
)

signature_train_score <- calculate_signature_score(
  signature_prep_final$x_train,
  signature_weights
)

signature_test_score <- calculate_signature_score(
  signature_prep_final$x_test,
  signature_weights
)

signature_curve_model <- fit_signature_curve(signature_train_score, y_train)

signature_train_kernel_log <- predict_signature_kernel(
  train_score = signature_train_score,
  train_y = y_train,
  test_score = signature_train_score,
  k = best_signature_k
)

signature_test_kernel_log <- predict_signature_kernel(
  train_score = signature_train_score,
  train_y = y_train,
  test_score = signature_test_score,
  k = best_signature_k
)

signature_train_curve_log <- predict_signature_curve(
  signature_curve_model,
  signature_train_score
)

signature_test_curve_log <- predict_signature_curve(
  signature_curve_model,
  signature_test_score
)

signature_train_pred_log <- best_signature_blend * signature_train_kernel_log +
  (1 - best_signature_blend) * signature_train_curve_log

signature_test_pred_log <- best_signature_blend * signature_test_kernel_log +
  (1 - best_signature_blend) * signature_test_curve_log

signature_train_pred_age <- inverse_age_transform(signature_train_pred_log)
signature_test_pred_age_raw <- inverse_age_transform(signature_test_pred_log)

signature_calibration_model <- fit_age_calibration(
  predicted_age = signature_train_pred_age,
  actual_age = actual_train_age
)

signature_test_pred_age <- apply_age_calibration(
  calibration_model = signature_calibration_model,
  predicted_age = signature_test_pred_age_raw,
  training_age = actual_train_age
)

signature_predictions <- data.frame(
  sample_id = rownames(X_test_all_raw),
  cohort = test_metadata$cohort,
  actual_age = as.numeric(test_metadata$age_numeric),
  predicted_age = signature_test_pred_age,
  predicted_log_age = signature_test_pred_log,
  model = "signature_kernel",
  stringsAsFactors = FALSE
)

status("CLOCK 2 finished: signature_kernel.")


# ------------------------------------------------------------
# 10. Combine predictions and calculate metrics
# ------------------------------------------------------------

status("Combining model predictions and calculating metrics.")

all_predictions <- rbind(
  ridge_predictions,
  signature_predictions
)

all_metrics <- data.frame()

for (model_name in unique(all_predictions$model)) {
  pred_df <- all_predictions[all_predictions$model == model_name, ]
  model_metrics <- evaluate_prediction_groups(
    pred_df = pred_df,
    model_name = model_name
  )
  
  all_metrics <- rbind(all_metrics, model_metrics)
}

all_metrics <- all_metrics[order(all_metrics$test_group, all_metrics$MAE), ]

write.csv(
  all_predictions,
  file = file.path(model_dir, "model_predictions_1.csv"),
  row.names = FALSE
)

write.csv(
  all_metrics,
  file = file.path(model_dir, "model_metrics_1.csv"),
  row.names = FALSE
)

status("Saved predictions to: ", file.path(model_dir, "model_predictions_1.csv"))
status("Saved metrics to: ", file.path(model_dir, "model_metrics_1.csv"))


# ------------------------------------------------------------
# 11. Save trained models and preprocessing objects
# ------------------------------------------------------------

trained_models <- list(
  selected_cpgs_loaded = selected_cpgs,
  selected_cpg_info_loaded = selected_cpg_info,
  train_idx = train_idx,
  test_idx = test_idx,
  ridge_clock = list(
    model_name = "m_value_ridge_calibrated",
    cpgs = ridge_cpgs,
    best_top_n = best_ridge_top_n,
    best_lambda = best_ridge_lambda,
    model = ridge_model_final,
    preprocessing = ridge_prep_final[c(
      "use_m_values",
      "imputation_means",
      "feature_means",
      "feature_sds"
    )],
    calibration_model = ridge_calibration_model,
    tuning_results = ridge_tuning_results
  ),
  signature_clock = list(
    model_name = "signature_kernel",
    cpgs = signature_cpgs,
    cpg_info = signature_cpg_info,
    weights = signature_weights,
    best_top_n = best_signature_top_n,
    best_k = best_signature_k,
    best_blend = best_signature_blend,
    train_score = signature_train_score,
    train_y = y_train,
    curve_model = signature_curve_model,
    preprocessing = signature_prep_final[c(
      "imputation_means",
      "feature_centers",
      "feature_scales"
    )],
    calibration_model = signature_calibration_model,
    tuning_results = signature_tuning_results
  )
)

saveRDS(
  trained_models,
  file = file.path(model_dir, "trained_models_1.rds")
)

status("Saved trained models to: ", file.path(model_dir, "trained_models_1.rds"))


# ------------------------------------------------------------
# 12. Save prediction plots
# ------------------------------------------------------------

status("Saving prediction plots.")

for (model_name in unique(all_predictions$model)) {
  pred_df <- all_predictions[all_predictions$model == model_name, ]
  
  output_file <- file.path(
    model_dir,
    paste0("prediction_plot_", model_name, "_1.png")
  )
  
  plot_predictions_base(
    pred_df = pred_df,
    model_name = model_name,
    output_file = output_file
  )
}

status("Prediction plots saved.")


# ------------------------------------------------------------
# 13. Final summary
# ------------------------------------------------------------

status("Model comparison finished.")
status("Results saved in: ", model_dir)

status("Model metrics:")
print(all_metrics)

status("Best model on all_test by MAE:")
print(all_metrics[all_metrics$test_group == "all_test", ][1, ])
