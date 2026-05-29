# ============================================================
# justinas.r
#
# Run from the scripts/ directory:
# source("justinas.r")
#
# Goal:
# Compare age prediction models using methods from course/task:
# 1. PCA + linear regression
# 2. kNN regression after PCA
# 3. Random forest regression
# 4. PCA + linear regression without Quinn as sensitivity analysis
#
# Input:
# ../results/cleaned_cohorts/*_cleaned.rds
# ../results/age_cpg_selection/top_5000_age_cpgs.rds
#
# Outputs:
# ../results/model_comparison/model_metrics.csv
# ../results/model_comparison/model_predictions.csv
# ../results/model_comparison/trained_models.rds
# ../results/model_comparison/prediction_plot_*.png
# ============================================================


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

needed_packages <- c("annmatrix", "ranger", "FNN")

for (pkg in needed_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(annmatrix)
library(ranger)
library(FNN)


# ------------------------------------------------------------
# 1. Paths and settings
# ------------------------------------------------------------

cleaned_dir <- file.path("..", "results", "cleaned_cohorts")
selection_dir <- file.path("..", "results", "age_cpg_selection")
model_dir <- file.path("..", "results", "model_comparison")

if (!dir.exists(model_dir)) {
  dir.create(model_dir)
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

selected_cpg_file <- file.path(selection_dir, "top_5000_age_cpgs.rds")

if (!file.exists(selected_cpg_file)) {
  stop("Could not find selected CpG file: ", selected_cpg_file)
}

set.seed(123)

train_fraction <- 0.80

# PCA + linear regression and kNN use this many CpGs.
model_top_n <- 5000

# Random forest is slower, so use fewer CpGs.
# If slow, reduce to 500.
rf_top_n <- 1000

# PCA settings
linear_pcs <- 100

# kNN tuning settings
pca_candidate_pcs <- c(20, 50, 100, 150)
knn_candidate_k <- c(3, 5, 7, 9, 15, 25)


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
  # We used log1p(age), so transform back with expm1().
  # pmax prevents impossible negative age predictions.
  pmax(expm1(predicted_log_age), 0)
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


scale_train_test <- function(x_train, x_test) {
  
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
    means = feature_means,
    sds = feature_sds
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


# ------------------------------------------------------------
# 3. Load selected CpGs
# ------------------------------------------------------------

status("Loading selected CpGs.")

selected_cpgs <- readRDS(selected_cpg_file)
selected_cpgs <- head(selected_cpgs, model_top_n)

rf_cpgs <- head(selected_cpgs, rf_top_n)

status("Selected CpGs for PCA linear and kNN/PCA: ", length(selected_cpgs))
status("Selected CpGs for random forest: ", length(rf_cpgs))


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
    stop("Some selected CpGs not found in cohort: ", cohort)
  }
  
  x_sub <- x[idx, , drop = FALSE]
  
  # Convert to samples x CpGs for modelling.
  mat <- t(as.matrix(x_sub))
  colnames(mat) <- selected_cpgs
  
  # Important metadata fix:
  # different cohorts may have different annotation columns,
  # so keep only columns needed for modelling.
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
    stop("NA values found in age_numeric after conversion in cohort: ", cohort)
  }
  
  if (any(is.na(meta$age_for_model))) {
    stop("NA values found in age_for_model after conversion in cohort: ", cohort)
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

# Align metadata to matrix rows.
metadata <- metadata[rownames(X), , drop = FALSE]

if (any(is.na(metadata$age_for_model))) {
  stop("age_for_model contains NA after metadata alignment.")
}

y <- as.numeric(metadata$age_for_model)
actual_age <- as.numeric(metadata$age_numeric)

if (length(unique(y)) <= 1) {
  print(summary(y))
  print(table(metadata$cohort))
  stop("Model target y is constant. Check metadata extraction.")
}

status("Final modelling matrix:")
status("Samples: ", nrow(X))
status("CpGs: ", ncol(X))
status("Age range: ", min(actual_age), " to ", max(actual_age), " years.")
status("Unique transformed age values: ", length(unique(y)))
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

X_train_raw <- X[train_idx, , drop = FALSE]
X_test_raw <- X[test_idx, , drop = FALSE]

y_train <- y[train_idx]
y_test <- y[test_idx]

test_metadata <- metadata[test_idx, , drop = FALSE]

status("Training age_for_model summary:")
print(summary(y_train))
status("Unique training y values: ", length(unique(y_train)))

if (length(unique(y_train)) <= 1) {
  stop("Training target y_train is constant. Model cannot be trained.")
}


# ------------------------------------------------------------
# 6. Scale features using training data only
# ------------------------------------------------------------

status("Scaling features using training data only.")

scaled <- scale_train_test(X_train_raw, X_test_raw)

X_train <- scaled$x_train
X_test <- scaled$x_test

rm(X_train_raw, X_test_raw)
gc()


# ------------------------------------------------------------
# 7. Model 1: PCA + linear regression, all cohorts
# ------------------------------------------------------------

status("MODEL 1 started: PCA + linear regression, all cohorts.")

max_linear_pcs <- min(linear_pcs, nrow(X_train) - 1, ncol(X_train))

status("MODEL 1: fitting PCA with ", max_linear_pcs, " PCs.")

pca_linear_model <- prcomp(
  X_train,
  center = FALSE,
  scale. = FALSE,
  rank. = max_linear_pcs
)

train_scores_linear <- as.data.frame(pca_linear_model$x)
train_scores_linear$age_for_model <- y_train

linear_model <- lm(age_for_model ~ ., data = train_scores_linear)

test_scores_linear <- as.data.frame(
  predict(pca_linear_model, newdata = X_test)
)

linear_pred_log <- as.numeric(
  predict(linear_model, newdata = test_scores_linear)
)

linear_pred_age <- inverse_age_transform(linear_pred_log)

linear_predictions <- data.frame(
  sample_id = rownames(X_test),
  cohort = test_metadata$cohort,
  actual_age = as.numeric(test_metadata$age_numeric),
  predicted_age = linear_pred_age,
  predicted_log_age = linear_pred_log,
  model = "pca_linear_all",
  stringsAsFactors = FALSE
)

rm(train_scores_linear, test_scores_linear)
gc()

status("MODEL 1 finished: PCA + linear regression, all cohorts.")


# ------------------------------------------------------------
# 8. Model 2: kNN regression after PCA, all cohorts
# ------------------------------------------------------------

status("MODEL 2 started: kNN regression after PCA, all cohorts.")

# Internal validation split within training data.
internal_meta <- metadata[train_idx, , drop = FALSE]
internal_split <- make_train_test_split(internal_meta, train_fraction = 0.80)

subtrain_local <- internal_split$train
valid_local <- internal_split$test

X_subtrain <- X_train[subtrain_local, , drop = FALSE]
X_valid <- X_train[valid_local, , drop = FALSE]

y_subtrain <- y_train[subtrain_local]
y_valid <- y_train[valid_local]

max_pcs <- min(max(pca_candidate_pcs), nrow(X_subtrain) - 1, ncol(X_subtrain))

status("MODEL 2: fitting PCA for tuning. Max PCs: ", max_pcs)

pca_tune <- prcomp(
  X_subtrain,
  center = FALSE,
  scale. = FALSE,
  rank. = max_pcs
)

subtrain_scores <- pca_tune$x
valid_scores_all <- predict(pca_tune, newdata = X_valid)

tuning_results <- data.frame()

status("MODEL 2: tuning k and number of PCs.")

for (pcs in pca_candidate_pcs) {
  
  if (pcs > ncol(subtrain_scores)) {
    next
  }
  
  for (k in knn_candidate_k) {
    
    if (k >= nrow(subtrain_scores)) {
      next
    }
    
    pred_valid_log <- FNN::knn.reg(
      train = subtrain_scores[, seq_len(pcs), drop = FALSE],
      test = valid_scores_all[, seq_len(pcs), drop = FALSE],
      y = y_subtrain,
      k = k
    )$pred
    
    pred_valid_age <- inverse_age_transform(pred_valid_log)
    actual_valid_age <- inverse_age_transform(y_valid)
    
    mae <- mean(abs(actual_valid_age - pred_valid_age))
    
    tuning_results <- rbind(
      tuning_results,
      data.frame(
        pcs = pcs,
        k = k,
        validation_MAE = mae
      )
    )
  }
}

tuning_results <- tuning_results[order(tuning_results$validation_MAE), ]

best_pcs <- tuning_results$pcs[1]
best_k <- tuning_results$k[1]

status("MODEL 2: best settings: PCs = ", best_pcs, ", k = ", best_k)

write.csv(
  tuning_results,
  file = file.path(model_dir, "knn_pca_tuning_results.csv"),
  row.names = FALSE
)

status("MODEL 2: fitting final PCA.")

pca_knn_model <- prcomp(
  X_train,
  center = FALSE,
  scale. = FALSE,
  rank. = best_pcs
)

train_scores_knn <- pca_knn_model$x[, seq_len(best_pcs), drop = FALSE]
test_scores_knn <- predict(pca_knn_model, newdata = X_test)[, seq_len(best_pcs), drop = FALSE]

knn_pred_log <- FNN::knn.reg(
  train = train_scores_knn,
  test = test_scores_knn,
  y = y_train,
  k = best_k
)$pred

knn_pred_age <- inverse_age_transform(knn_pred_log)

knn_predictions <- data.frame(
  sample_id = rownames(X_test),
  cohort = test_metadata$cohort,
  actual_age = as.numeric(test_metadata$age_numeric),
  predicted_age = knn_pred_age,
  predicted_log_age = knn_pred_log,
  model = "knn_pca_all",
  stringsAsFactors = FALSE
)

rm(
  X_subtrain,
  X_valid,
  subtrain_scores,
  valid_scores_all,
  train_scores_knn,
  test_scores_knn
)
gc()

status("MODEL 2 finished: kNN regression after PCA, all cohorts.")


# ------------------------------------------------------------
# 9. Model 3: Random forest regression, all cohorts
# ------------------------------------------------------------

status("MODEL 3 started: Random forest regression, all cohorts.")

rf_features <- rf_cpgs

X_train_rf <- X_train[, rf_features, drop = FALSE]
X_test_rf <- X_test[, rf_features, drop = FALSE]

rf_train_df <- data.frame(
  age_for_model = y_train,
  X_train_rf,
  check.names = FALSE
)

rf_threads <- max(1, parallel::detectCores() - 1)

status("MODEL 3: fitting random forest with ", ncol(X_train_rf), " CpGs.")
status("MODEL 3: using ", rf_threads, " threads.")

rf_model <- ranger(
  formula = age_for_model ~ .,
  data = rf_train_df,
  num.trees = 500,
  mtry = floor(sqrt(ncol(X_train_rf))),
  min.node.size = 5,
  importance = "none",
  num.threads = rf_threads,
  seed = 123
)

rf_pred_log <- predict(
  rf_model,
  data = data.frame(X_test_rf, check.names = FALSE)
)$predictions

rf_pred_age <- inverse_age_transform(rf_pred_log)

rf_predictions <- data.frame(
  sample_id = rownames(X_test),
  cohort = test_metadata$cohort,
  actual_age = as.numeric(test_metadata$age_numeric),
  predicted_age = rf_pred_age,
  predicted_log_age = rf_pred_log,
  model = "random_forest_all",
  stringsAsFactors = FALSE
)

rm(rf_train_df, X_train_rf, X_test_rf)
gc()

status("MODEL 3 finished: Random forest regression, all cohorts.")


# ------------------------------------------------------------
# 10. Sensitivity model: PCA + linear regression without Quinn
# ------------------------------------------------------------

status("SENSITIVITY MODEL started: PCA + linear regression without Quinn.")

non_quinn_train <- metadata$cohort[train_idx] != "quinn"

if (sum(non_quinn_train) < 20) {
  stop("Too few non-Quinn training samples.")
}

X_train_no_quinn <- X_train[non_quinn_train, , drop = FALSE]
y_train_no_quinn <- y_train[non_quinn_train]

max_linear_no_quinn_pcs <- min(
  linear_pcs,
  nrow(X_train_no_quinn) - 1,
  ncol(X_train_no_quinn)
)

status("SENSITIVITY MODEL: fitting PCA with ", max_linear_no_quinn_pcs, " PCs.")

pca_linear_no_quinn_model <- prcomp(
  X_train_no_quinn,
  center = FALSE,
  scale. = FALSE,
  rank. = max_linear_no_quinn_pcs
)

train_scores_no_quinn <- as.data.frame(pca_linear_no_quinn_model$x)
train_scores_no_quinn$age_for_model <- y_train_no_quinn

linear_no_quinn_model <- lm(age_for_model ~ ., data = train_scores_no_quinn)

test_scores_no_quinn <- as.data.frame(
  predict(pca_linear_no_quinn_model, newdata = X_test)
)

linear_no_quinn_pred_log <- as.numeric(
  predict(linear_no_quinn_model, newdata = test_scores_no_quinn)
)

linear_no_quinn_pred_age <- inverse_age_transform(linear_no_quinn_pred_log)

linear_no_quinn_predictions <- data.frame(
  sample_id = rownames(X_test),
  cohort = test_metadata$cohort,
  actual_age = as.numeric(test_metadata$age_numeric),
  predicted_age = linear_no_quinn_pred_age,
  predicted_log_age = linear_no_quinn_pred_log,
  model = "pca_linear_without_quinn",
  stringsAsFactors = FALSE
)

rm(
  X_train_no_quinn,
  y_train_no_quinn,
  train_scores_no_quinn,
  test_scores_no_quinn
)
gc()

status("SENSITIVITY MODEL finished: PCA + linear regression without Quinn.")


# ------------------------------------------------------------
# 11. Combine predictions and calculate metrics
# ------------------------------------------------------------

status("Combining predictions and calculating metrics.")

all_predictions <- rbind(
  linear_predictions,
  knn_predictions,
  rf_predictions,
  linear_no_quinn_predictions
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
  file = file.path(model_dir, "model_predictions.csv"),
  row.names = FALSE
)

write.csv(
  all_metrics,
  file = file.path(model_dir, "model_metrics.csv"),
  row.names = FALSE
)

status("Saved predictions to: ", file.path(model_dir, "model_predictions.csv"))
status("Saved metrics to: ", file.path(model_dir, "model_metrics.csv"))


# ------------------------------------------------------------
# 12. Save models
# ------------------------------------------------------------

trained_models <- list(
  selected_cpgs = selected_cpgs,
  rf_cpgs = rf_cpgs,
  scaling_means = scaled$means,
  scaling_sds = scaled$sds,
  train_idx = train_idx,
  test_idx = test_idx,
  pca_linear_all = list(
    pca = pca_linear_model,
    model = linear_model,
    pcs = max_linear_pcs
  ),
  knn_pca_all = list(
    pca = pca_knn_model,
    best_pcs = best_pcs,
    best_k = best_k,
    tuning_results = tuning_results
  ),
  random_forest_all = rf_model,
  pca_linear_without_quinn = list(
    pca = pca_linear_no_quinn_model,
    model = linear_no_quinn_model,
    pcs = max_linear_no_quinn_pcs
  )
)

saveRDS(
  trained_models,
  file = file.path(model_dir, "trained_models.rds")
)

status("Saved trained models to: ", file.path(model_dir, "trained_models.rds"))


# ------------------------------------------------------------
# 13. Save prediction plots
# ------------------------------------------------------------

status("Saving prediction plots.")

for (model_name in unique(all_predictions$model)) {
  
  pred_df <- all_predictions[all_predictions$model == model_name, ]
  
  output_file <- file.path(
    model_dir,
    paste0("prediction_plot_", model_name, ".png")
  )
  
  plot_predictions_base(
    pred_df = pred_df,
    model_name = model_name,
    output_file = output_file
  )
}

status("Prediction plots saved.")


# ------------------------------------------------------------
# 14. Final summary
# ------------------------------------------------------------

status("Model comparison finished.")
status("Results saved in: ", model_dir)

status("Model metrics:")
print(all_metrics)

status("Best model on all_test by MAE:")
print(all_metrics[all_metrics$test_group == "all_test", ][1, ])

status("Best model on non_quinn_test by MAE:")
print(all_metrics[all_metrics$test_group == "non_quinn_test", ][1, ])