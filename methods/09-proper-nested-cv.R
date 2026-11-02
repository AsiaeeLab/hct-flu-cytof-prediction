# ==============================================================================
# Proper Nested Cross-Validation (k tuned inside outer fold)
# ==============================================================================
# This implements statistically correct nested CV where k selection happens
# INSIDE the outer fold, not on the outer test set.

library(caret)

#' Proper nested CV with k selection inside outer fold
#'
#' This function implements the correct statistical procedure:
#' - Outer loop: Performance estimation
#' - Inner loop: Tune BOTH k and model hyperparameters
#' - Test set touched only once per outer fold
#'
#' @param mini_data Cell-level feature matrix
#' @param p_id Vector mapping cells to patient IDs
#' @param patient_outcomes List with y, x, and optionally covariates
#' @param k_range Vector of k values to try
#' @param clustering_method Clustering method name
#' @param prediction_method Prediction method name
#' @param outer_folds Pre-generated outer CV folds (for consistency across methods)
#' @param lambda_grid Grid of lambda values (for elasticnet)
#' @param mtry_grid Grid of mtry values (for RF)
#' @param nodesize_grid Grid of nodesize values (for RF)
#' @param alpha Elastic net alpha
#' @param scale_features Whether to scale features before clustering
#' @param verbose Print progress
#' @return List with results including k chosen per fold
proper_nested_cv <- function(mini_data, p_id, patient_outcomes,
                             k_range,
                             clustering_method, prediction_method,
                             outer_folds,
                             lambda_grid = 10^seq(-3, 1, length.out = 20),
                             mtry_grid = NULL,
                             nodesize_grid = c(1, 5, 10),
                             alpha = 0.5,
                             scale_features = TRUE,
                             verbose = FALSE) {

  n_outer_folds <- length(outer_folds)
  unique_p_ids <- names(patient_outcomes$y)

  # Storage for results
  results <- list(
    test_errors = numeric(n_outer_folds),
    train_errors = numeric(n_outer_folds),
    k_selected = integer(n_outer_folds),
    test_predictions = list(),
    test_true = list(),
    # NEW: Store inner CV curves per fold
    inner_cv_training = list(),      # inner_cv_training[[fold]] = training error vs k
    inner_cv_validation = list()     # inner_cv_validation[[fold]] = validation error vs k
  )

  # Outer CV loop
  for (fold_id in 1:n_outer_folds) {
    if (verbose) cat("    Outer fold", fold_id, "/", n_outer_folds, "...")

    # Split patients
    test_p_ids <- unique_p_ids[outer_folds[[fold_id]]]
    train_p_ids <- setdiff(unique_p_ids, test_p_ids)

    # Split cells
    train_idx <- which(p_id %in% train_p_ids)
    test_idx <- which(p_id %in% test_p_ids)

    train_data_raw <- mini_data[train_idx, , drop = FALSE]
    test_data_raw <- mini_data[test_idx, , drop = FALSE]
    train_p_id <- p_id[train_idx]
    test_p_id <- p_id[test_idx]

    # Feature scaling (inside fold!)
    if (scale_features) {
      sc <- scale(train_data_raw)
      train_data <- sc
      test_data <- scale(test_data_raw,
                        center = attr(sc, "scaled:center"),
                        scale = attr(sc, "scaled:scale"))
    } else {
      train_data <- train_data_raw
      test_data <- test_data_raw
    }

    # === INNER LOOP: Tune k and hyperparameters ===
    best <- tune_k_inside_outer(
      train_data, train_p_id, patient_outcomes,
      k_range, clustering_method, prediction_method,
      lambda_grid, mtry_grid, nodesize_grid,
      alpha, n_inner_folds = 5
    )

    best_k <- best$k
    best_hyper <- best$hyper
    results$k_selected[fold_id] <- best_k

    # Store inner CV curves (both training and validation)
    results$inner_cv_training[[fold_id]] <- best$inner_cv_training_errors
    results$inner_cv_validation[[fold_id]] <- best$inner_cv_validation_errors

    # === Refit on full outer training set with best k ===
    # Cluster training data
    cluster_result <- perform_clustering(train_data, best_k, method = clustering_method)
    train_clusters <- cluster_result$clusters

    # Build training features
    train_covariates <- if (!is.null(patient_outcomes$covariates)) {
      patient_outcomes$covariates[train_p_ids, , drop = FALSE]
    } else NULL

    if (!is.null(train_covariates)) {
      train_features <- build_patient_features_with_covariates(
        train_p_id, train_clusters,
        patient_outcomes$x[train_p_ids],
        covariates = train_covariates
      )
    } else {
      train_features <- build_patient_features(
        train_p_id, train_clusters,
        patient_outcomes$x[train_p_ids]
      )
    }

    # Assign test cells to training clusters
    test_clusters <- assign_to_clusters(best_k, train_data, train_clusters, test_data)

    # Build test features
    test_features_raw <- as.data.frame.matrix(table(test_p_id, test_clusters))
    test_features_raw <- log(test_features_raw + 1)

    # Align columns
    non_cluster_cols <- c("baseline")
    if (!is.null(train_covariates)) {
      non_cluster_cols <- c(non_cluster_cols, colnames(train_covariates))
    }
    train_cluster_cols <- setdiff(colnames(train_features), non_cluster_cols)

    test_features_aligned <- test_features_raw
    missing_clusters <- setdiff(train_cluster_cols, colnames(test_features_aligned))
    for (col in missing_clusters) {
      test_features_aligned[[col]] <- 0
    }
    test_features_aligned <- test_features_aligned[, train_cluster_cols, drop = FALSE]

    # Add baseline and covariates to test
    test_patient_ids <- rownames(test_features_aligned)
    test_baseline <- patient_outcomes$x[test_patient_ids]
    test_features <- cbind(test_features_aligned, baseline = test_baseline)

    if (!is.null(train_covariates)) {
      test_covariates <- patient_outcomes$covariates[test_patient_ids, , drop = FALSE]
      test_features <- cbind(test_features, test_covariates)
    }

    # Convert to matrix
    train_features <- as.matrix(train_features)
    test_features <- as.matrix(test_features)

    # Train final model with best hyperparameters
    train_y <- patient_outcomes$y[rownames(train_features)]
    test_y <- patient_outcomes$y[rownames(test_features)]

    if (prediction_method == "elasticnet") {
      model <- glmnet(train_features, train_y,
                     alpha = best_hyper$alpha,
                     lambda = best_hyper$lambda)
      test_pred <- predict(model, newx = test_features, s = best_hyper$lambda)
      train_pred <- predict(model, newx = train_features, s = best_hyper$lambda)
    } else if (prediction_method == "random_forest") {
      train_df <- as.data.frame(train_features)
      train_df$y <- train_y
      model <- randomForest(y ~ ., data = train_df,
                           ntree = best_hyper$ntree,
                           mtry = best_hyper$mtry,
                           nodesize = best_hyper$nodesize)
      test_pred <- predict(model, newdata = as.data.frame(test_features))
      train_pred <- predict(model, newdata = train_df)
    }

    # Store results
    results$test_errors[fold_id] <- mean((test_pred - test_y)^2)
    results$train_errors[fold_id] <- mean((train_pred - train_y)^2)
    results$test_predictions[[fold_id]] <- as.numeric(test_pred)
    results$test_true[[fold_id]] <- as.numeric(test_y)

    if (verbose) {
      cat(" k=", best_k, " Test:", round(results$test_errors[fold_id], 4),
          " Train:", round(results$train_errors[fold_id], 4), "\n", sep = "")
    }
  }

  # Aggregate results
  return(list(
    mean_test_error = mean(results$test_errors),
    se_test_error = sd(results$test_errors) / sqrt(n_outer_folds),
    mean_train_error = mean(results$train_errors),
    se_train_error = sd(results$train_errors) / sqrt(n_outer_folds),
    k_selected = results$k_selected,
    test_errors_per_fold = results$test_errors,
    train_errors_per_fold = results$train_errors,
    test_predictions = results$test_predictions,
    test_true = results$test_true,
    # Inner CV curves for visualization (what happened inside each outer fold)
    inner_cv_training = results$inner_cv_training,      # Training error vs k per fold
    inner_cv_validation = results$inner_cv_validation   # Validation error vs k per fold
  ))
}


#' Baseline-only control model (no clusters)
#'
#' Trains a model using ONLY baseline and covariates (no cellular features).
#' This provides a control to measure how much the cluster features actually help.
#'
#' @param patient_outcomes List with y, x, and optionally covariates
#' @param outer_folds Pre-generated outer CV folds
#' @param prediction_method Prediction method name
#' @param lambda_grid Grid of lambda values (for elasticnet)
#' @param mtry_grid Grid of mtry values (for RF)
#' @param nodesize_grid Grid of nodesize values (for RF)
#' @param alpha Elastic net alpha
#' @param verbose Print progress
#' @return List with CV results
baseline_only_cv <- function(patient_outcomes,
                             outer_folds,
                             prediction_method = "elasticnet",
                             lambda_grid = 10^seq(-3, 1, length.out = 20),
                             mtry_grid = NULL,
                             nodesize_grid = c(1, 5, 10),
                             alpha = 0.5,
                             verbose = FALSE) {

  n_outer_folds <- length(outer_folds)
  unique_p_ids <- names(patient_outcomes$y)

  # Storage
  results <- list(
    test_errors = numeric(n_outer_folds),
    train_errors = numeric(n_outer_folds),
    test_predictions = list(),
    test_true = list()
  )

  for (fold_id in 1:n_outer_folds) {
    if (verbose) cat("    Outer fold", fold_id, "/", n_outer_folds, "...")

    # Split patients
    test_p_ids <- unique_p_ids[outer_folds[[fold_id]]]
    train_p_ids <- setdiff(unique_p_ids, test_p_ids)

    # Build features: baseline + covariates ONLY (no clusters!)
    train_baseline <- patient_outcomes$x[train_p_ids]
    test_baseline <- patient_outcomes$x[test_p_ids]

    if (!is.null(patient_outcomes$covariates)) {
      train_cov <- patient_outcomes$covariates[train_p_ids, , drop = FALSE]
      test_cov <- patient_outcomes$covariates[test_p_ids, , drop = FALSE]

      train_features <- cbind(baseline = train_baseline, train_cov)
      test_features <- cbind(baseline = test_baseline, test_cov)
    } else {
      train_features <- matrix(train_baseline, ncol = 1)
      colnames(train_features) <- "baseline"
      rownames(train_features) <- train_p_ids

      test_features <- matrix(test_baseline, ncol = 1)
      colnames(test_features) <- "baseline"
      rownames(test_features) <- test_p_ids
    }

    train_y <- patient_outcomes$y[train_p_ids]
    test_y <- patient_outcomes$y[test_p_ids]

    # Tune hyperparameters using inner CV
    if (prediction_method == "elasticnet") {
      gs <- grid_search_elasticnet(train_features, train_y, alpha = alpha,
                                   lambda_grid = lambda_grid, n_inner_folds = 5)
      best_lambda <- gs$best_lambda

      # Train with best lambda
      model <- glmnet(train_features, train_y, alpha = alpha, lambda = best_lambda)
      test_pred <- predict(model, newx = test_features, s = best_lambda)
      train_pred <- predict(model, newx = train_features, s = best_lambda)

    } else if (prediction_method == "random_forest") {
      gs <- grid_search_rf(train_features, train_y, mtry_grid = mtry_grid,
                          nodesize_grid = nodesize_grid, ntree = 500, n_inner_folds = 3)

      train_df <- as.data.frame(train_features)
      train_df$y <- train_y
      model <- randomForest(y ~ ., data = train_df, ntree = 500,
                           mtry = gs$best_mtry, nodesize = gs$best_nodesize)
      test_pred <- predict(model, newdata = as.data.frame(test_features))
      train_pred <- predict(model, newdata = train_df)
    }

    # Store results
    results$test_errors[fold_id] <- mean((test_pred - test_y)^2)
    results$train_errors[fold_id] <- mean((train_pred - train_y)^2)
    results$test_predictions[[fold_id]] <- as.numeric(test_pred)
    results$test_true[[fold_id]] <- as.numeric(test_y)

    if (verbose) {
      cat(" Test:", round(results$test_errors[fold_id], 4),
          " Train:", round(results$train_errors[fold_id], 4), "\n")
    }
  }

  return(list(
    mean_test_error = mean(results$test_errors),
    se_test_error = sd(results$test_errors) / sqrt(n_outer_folds),
    mean_train_error = mean(results$train_errors),
    se_train_error = sd(results$train_errors) / sqrt(n_outer_folds),
    test_errors_per_fold = results$test_errors,
    train_errors_per_fold = results$train_errors,
    test_predictions = results$test_predictions,
    test_true = results$test_true
  ))
}


#' Proper nested CV with subsampling
#'
#' Runs proper nested CV with multiple cell subsamples for robustness.
#' Key improvements:
#' - k tuned inside outer fold (unbiased)
#' - Patient-stratified cell subsampling
#' - Fixed outer folds across all k values
#' - Feature scaling inside folds
#'
#' @param full_data Full cell-level feature matrix
#' @param full_p_id Vector mapping all cells to patient IDs
#' @param outcome_obj Outcome data frame
#' @param meta_obj Metadata data frame with covariates (optional)
#' @param outcome_col_1 Baseline outcome column
#' @param outcome_col_2 Target outcome column
#' @param clustering_methods Vector of clustering methods to test
#' @param prediction_methods Vector of prediction methods to test
#' @param k_range Vector of k values to try
#' @param n_cells Target number of cells per subsample
#' @param n_subsamples Number of subsampling iterations
#' @param n_outer_folds Number of outer CV folds
#' @param lambda_grid Lambda grid for elasticnet
#' @param alpha Elasticnet alpha
#' @param include_dose Include dose covariate
#' @param include_age Include age covariate
#' @param include_sex Include sex covariate
#' @param include_time Include time covariate
#' @param scale_features Whether to scale features before clustering
#' @param patient_stratified_sampling Whether to sample cells per patient
#' @param verbose Print progress
#' @return Comprehensive results structure
proper_nested_cv_with_subsampling <- function(
    full_data, full_p_id, outcome_obj,
    meta_obj = NULL,
    outcome_col_1, outcome_col_2,
    clustering_methods = c("kmeans", "flowsom"),
    prediction_methods = c("elasticnet", "random_forest"),
    k_range = seq(10, 50, by = 10),
    n_cells = 10000,
    n_subsamples = 3,
    n_outer_folds = 5,
    lambda_grid = 10^seq(-3, 1, length.out = 20),
    alpha = 0.5,
    include_dose = TRUE,
    include_age = TRUE,
    include_sex = TRUE,
    include_time = FALSE,
    scale_features = TRUE,
    patient_stratified_sampling = TRUE,
    outcome_type = "fold_change",  # "fold_change" or "absolute"
    n_cores = 1,  # NEW: number of cores for parallel processing
    verbose = TRUE) {

  if (verbose) {
    cat("\n", rep("=", 70), "\n", sep = "")
    cat("PROPER NESTED CV WITH SUBSAMPLING\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("Parameters:\n")
    cat("  - n_cells per subsample:", n_cells, "\n")
    cat("  - n_subsamples:", n_subsamples, "\n")
    cat("  - k_range:", paste(k_range, collapse = ", "), "\n")
    cat("  - Clustering methods:", paste(clustering_methods, collapse = ", "), "\n")
    cat("  - Prediction methods:", paste(prediction_methods, collapse = ", "), "\n")
    cat("  - Outer CV folds:", n_outer_folds, "\n")
    cat("  - Patient-stratified sampling:", patient_stratified_sampling, "\n")
    cat("  - Feature scaling:", scale_features, "\n")
    cat("  - Parallel cores:", n_cores, "\n")
    cat("  - Outcome type:", outcome_type, "\n")
    if (outcome_type == "absolute") {
      cat("    (predicting log(v3) with log(v1) as feature)\n")
    } else {
      cat("    (predicting log(v3/v1) fold-change)\n")
    }
    covars <- c()
    if (include_dose) covars <- c(covars, "dose")
    if (include_age) covars <- c(covars, "age")
    if (include_sex) covars <- c(covars, "sex")
    if (include_time) covars <- c(covars, "time")
    cat("  - Covariates:", paste(covars, collapse = ", "), "\n\n")
  }

  # Prepare outcomes (once, before subsampling)
  if (!is.null(meta_obj)) {
    if (outcome_type == "absolute") {
      patient_outcomes <- prepare_outcome_with_covariates_absolute(
        outcome_obj, meta_obj, outcome_col_1, outcome_col_2, full_p_id,
        include_dose = include_dose,
        include_age = include_age,
        include_sex = include_sex,
        include_time = include_time
      )
    } else {
      patient_outcomes <- prepare_outcome_with_covariates(
        outcome_obj, meta_obj, outcome_col_1, outcome_col_2, full_p_id,
        include_dose = include_dose,
        include_age = include_age,
        include_sex = include_sex,
        include_time = include_time
      )
    }
  } else {
    if (outcome_type == "absolute") {
      patient_outcomes <- prepare_outcome_absolute(outcome_obj, outcome_col_1, outcome_col_2, full_p_id)
    } else {
      patient_outcomes <- prepare_outcome(outcome_obj, outcome_col_1, outcome_col_2, full_p_id)
    }
  }

  # Initialize results
  results <- list()
  method_combinations <- expand.grid(
    clustering = clustering_methods,
    prediction = prediction_methods,
    stringsAsFactors = FALSE
  )

  # Function to process one subsample (for parallelization)
  process_subsample <- function(subsample_idx) {
    # Patient-stratified cell subsampling
    set.seed(2025 + subsample_idx)
    if (patient_stratified_sampling) {
      by_pid <- split(seq_len(nrow(full_data)), full_p_id)
      target_per_patient <- ceiling(n_cells / length(by_pid))
      sample_idx <- unlist(lapply(by_pid, function(ix) {
        sample(ix, size = min(length(ix), target_per_patient))
      }))
    } else {
      sample_idx <- sample(nrow(full_data), n_cells)
    }

    mini_data <- full_data[sample_idx, , drop = FALSE]
    p_id <- full_p_id[sample_idx]

    # Create outer folds ONCE for this subsample (re-used across all methods)
    set.seed(2025 + subsample_idx)  # Same seed = same folds
    unique_p_ids <- names(patient_outcomes$y)
    outer_folds <- createFolds(seq_along(unique_p_ids), k = n_outer_folds, list = TRUE)

    # Results for this subsample
    subsample_results <- list()

    # Test each method combination
    for (i in 1:nrow(method_combinations)) {
      clustering_method <- method_combinations$clustering[i]
      prediction_method <- method_combinations$prediction[i]
      method_name <- paste(clustering_method, prediction_method, sep = "_")

      # Run proper nested CV
      cv_result <- proper_nested_cv(
        mini_data, p_id, patient_outcomes,
        k_range = k_range,
        clustering_method = clustering_method,
        prediction_method = prediction_method,
        outer_folds = outer_folds,
        lambda_grid = lambda_grid,
        mtry_grid = NULL,
        nodesize_grid = c(1, 5, 10),
        alpha = alpha,
        scale_features = scale_features,
        verbose = FALSE  # Disable verbose in parallel
      )

      subsample_results[[method_name]] <- cv_result
    }

    # === BASELINE-ONLY CONTROL ===
    # Run baseline-only model for each prediction method (no clustering needed!)
    for (pred_method in prediction_methods) {
      method_name <- paste0("baseline_only_", pred_method)

      baseline_result <- baseline_only_cv(
        patient_outcomes = patient_outcomes,
        outer_folds = outer_folds,
        prediction_method = pred_method,
        lambda_grid = lambda_grid,
        mtry_grid = NULL,
        nodesize_grid = c(1, 5, 10),
        alpha = alpha,
        verbose = FALSE  # Disable verbose in parallel
      )

      subsample_results[[method_name]] <- baseline_result
    }

    return(subsample_results)
  }

  # Run subsamples (in parallel if n_cores > 1)
  if (n_cores > 1) {
    if (verbose) {
      cat("Running", n_subsamples, "subsamples in parallel using", n_cores, "cores...\n")
    }

    # Load required packages on each worker
    library(parallel)
    library(doParallel)
    library(foreach)

    # Setup parallel backend
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)

    # Load required functions on each worker
    clusterEvalQ(cl, {
      library(here)
      source(here::here("methods", "00-paths.R"))
      source(here::here("methods", "02-data-preparation.R"))
      source(here::here("methods", "03-clustering-methods.R"))
      source(here::here("methods", "04-prediction-models.R"))
      source(here::here("methods", "07-nested-cv-pipeline.R"))
      source(here::here("methods", "08-multi-outcome-covariates.R"))
      # Note: Don't source 09 (current file) to avoid recursion
    })

    # Export necessary objects to workers (including functions from this file)
    clusterExport(cl, c("full_data", "full_p_id", "patient_outcomes", "n_cells",
                        "patient_stratified_sampling", "n_outer_folds",
                        "method_combinations", "k_range", "lambda_grid", "alpha",
                        "scale_features", "prediction_methods",
                        "proper_nested_cv", "baseline_only_cv"),
                  envir = environment())

    # Run in parallel
    subsample_results_list <- foreach(subsample_idx = 1:n_subsamples,
                                      .packages = c("glmnet", "randomForest", "caret", "FlowSOM", "flowCore")) %dopar% {
      process_subsample(subsample_idx)
    }

    # Stop cluster
    stopCluster(cl)

    if (verbose) cat("✓ Parallel processing complete\n")

  } else {
    # Sequential execution (original behavior)
    subsample_results_list <- list()
    for (subsample_idx in 1:n_subsamples) {
      if (verbose) {
        cat("\n", rep("-", 70), "\n", sep = "")
        cat("SUBSAMPLE ", subsample_idx, " / ", n_subsamples, "\n", sep = "")
        cat(rep("-", 70), "\n", sep = "")
      }

      subsample_results_list[[subsample_idx]] <- process_subsample(subsample_idx)

      if (verbose) {
        cat("✓ Subsample", subsample_idx, "complete\n")
      }
    }
  }

  # Reorganize results from list of subsamples to list of methods
  for (subsample_idx in 1:n_subsamples) {
    subsample_res <- subsample_results_list[[subsample_idx]]
    for (method_name in names(subsample_res)) {
      if (is.null(results[[method_name]])) {
        results[[method_name]] <- list()
      }
      results[[method_name]][[paste0("subsample", subsample_idx)]] <- subsample_res[[method_name]]
    }
  }

  # Aggregate across subsamples
  aggregated <- aggregate_proper_nested_cv(results, clustering_methods,
                                          prediction_methods, n_subsamples,
                                          include_baseline = TRUE)

  return(list(
    raw_results = results,
    aggregated = aggregated
  ))
}


#' Aggregate proper nested CV results across subsamples
#'
#' @param results Raw results from proper_nested_cv_with_subsampling
#' @param clustering_methods Vector of clustering methods
#' @param prediction_methods Vector of prediction methods
#' @param n_subsamples Number of subsamples
#' @return Aggregated results
aggregate_proper_nested_cv <- function(results, clustering_methods,
                                      prediction_methods, n_subsamples,
                                      include_baseline = TRUE) {

  method_combinations <- expand.grid(
    clustering = clustering_methods,
    prediction = prediction_methods,
    stringsAsFactors = FALSE
  )
  method_combinations$method_name <- paste(method_combinations$clustering,
                                          method_combinations$prediction, sep = "_")

  aggregated <- list()

  # Aggregate cluster-based methods
  for (i in 1:nrow(method_combinations)) {
    method_name <- method_combinations$method_name[i]
    method_data <- results[[method_name]]

    # Collect metrics across subsamples
    test_errors <- sapply(1:n_subsamples, function(s) {
      method_data[[paste0("subsample", s)]]$mean_test_error
    })
    train_errors <- sapply(1:n_subsamples, function(s) {
      method_data[[paste0("subsample", s)]]$mean_train_error
    })

    # Collect k selections
    k_selections <- lapply(1:n_subsamples, function(s) {
      method_data[[paste0("subsample", s)]]$k_selected
    })

    aggregated[[method_name]] <- list(
      clustering = method_combinations$clustering[i],
      prediction = method_combinations$prediction[i],
      mean_test_error = mean(test_errors),
      se_test_error = sd(test_errors) / sqrt(n_subsamples),
      mean_train_error = mean(train_errors),
      se_train_error = sd(train_errors) / sqrt(n_subsamples),
      k_selected_per_subsample = k_selections,
      k_mode = as.integer(names(sort(table(unlist(k_selections)), decreasing = TRUE)[1])),
      is_baseline_only = FALSE
    )
  }

  # Aggregate baseline-only methods
  if (include_baseline) {
    for (pred_method in prediction_methods) {
      method_name <- paste0("baseline_only_", pred_method)
      method_data <- results[[method_name]]

      if (!is.null(method_data)) {
        test_errors <- sapply(1:n_subsamples, function(s) {
          method_data[[paste0("subsample", s)]]$mean_test_error
        })
        train_errors <- sapply(1:n_subsamples, function(s) {
          method_data[[paste0("subsample", s)]]$mean_train_error
        })

        aggregated[[method_name]] <- list(
          clustering = "baseline_only",
          prediction = pred_method,
          mean_test_error = mean(test_errors),
          se_test_error = sd(test_errors) / sqrt(n_subsamples),
          mean_train_error = mean(train_errors),
          se_train_error = sd(train_errors) / sqrt(n_subsamples),
          k_selected_per_subsample = NULL,
          k_mode = NA,
          is_baseline_only = TRUE
        )
      }
    }
  }

  return(aggregated)
}

#' Diagnostic: Evaluate all k values (for visualization only)
#'
#' This function evaluates performance for each k value across CV folds.
#' IMPORTANT: Use this ONLY for visualization/diagnostics, NOT for model selection!
#' The proper nested CV above is the statistically correct approach.
#'
#' @param mini_data Cell-level feature matrix
#' @param p_id Vector mapping cells to patient IDs
#' @param patient_outcomes List with y, x, and optionally covariates
#' @param k_range Vector of k values to try
#' @param clustering_method Clustering method name
#' @param prediction_method Prediction method name
#' @param cv_folds Pre-generated CV folds
#' @param lambda_grid Grid of lambda values (for elasticnet)
#' @param alpha Elastic net alpha
#' @param scale_features Whether to scale features before clustering
#' @return List with mean_test_error and mean_train_error vectors (one per k)
diagnostic_evaluate_all_k <- function(mini_data, p_id, patient_outcomes,
                                     k_range,
                                     clustering_method, prediction_method,
                                     cv_folds,
                                     lambda_grid = 10^seq(-3, 1, length.out = 20),
                                     alpha = 0.5,
                                     scale_features = TRUE) {

  n_folds <- length(cv_folds)
  unique_p_ids <- names(patient_outcomes$y)

  # Storage: matrix with rows = folds, cols = k values
  test_errors <- matrix(NA, nrow = n_folds, ncol = length(k_range))
  train_errors <- matrix(NA, nrow = n_folds, ncol = length(k_range))

  for (k_idx in seq_along(k_range)) {
    k <- k_range[k_idx]

    for (fold_id in 1:n_folds) {
      # Split patients
      test_p_ids <- unique_p_ids[cv_folds[[fold_id]]]
      train_p_ids <- setdiff(unique_p_ids, test_p_ids)

      # Split cells
      train_idx <- which(p_id %in% train_p_ids)
      test_idx <- which(p_id %in% test_p_ids)

      train_data_raw <- mini_data[train_idx, , drop = FALSE]
      test_data_raw <- mini_data[test_idx, , drop = FALSE]
      train_p_id <- p_id[train_idx]
      test_p_id <- p_id[test_idx]

      # Feature scaling
      if (scale_features) {
        sc <- scale(train_data_raw)
        train_data <- sc
        test_data <- scale(test_data_raw,
                          center = attr(sc, "scaled:center"),
                          scale = attr(sc, "scaled:scale"))
      } else {
        train_data <- train_data_raw
        test_data <- test_data_raw
      }

      # Cluster with this k
      cluster_result <- perform_clustering(train_data, k, method = clustering_method)
      train_clusters <- cluster_result$clusters

      # Build training features (old API - just baseline, no covariates)
      train_features <- build_patient_features(
        p_id = train_p_id,
        clusters = train_clusters,
        baseline_outcome = patient_outcomes$x[train_p_ids]
      )

      # Add covariates if present
      if (!is.null(patient_outcomes$covariates)) {
        train_cov <- patient_outcomes$covariates[rownames(train_features), , drop = FALSE]
        train_features <- cbind(train_features, train_cov)
      }

      # Assign test cells to clusters
      test_clusters <- assign_to_clusters(k, train_data, train_clusters, test_data)

      test_features <- build_patient_features(
        p_id = test_p_id,
        clusters = test_clusters,
        baseline_outcome = patient_outcomes$x[test_p_ids]
      )

      # Add covariates if present
      if (!is.null(patient_outcomes$covariates)) {
        test_cov <- patient_outcomes$covariates[rownames(test_features), , drop = FALSE]
        test_features <- cbind(test_features, test_cov)
      }

      # Align features (convert to data.frame first since align_features expects data.frame)
      train_colnames <- colnames(train_features)
      test_features <- as.data.frame(test_features)
      test_features <- align_features(test_features, train_colnames)
      test_features <- as.matrix(test_features)  # Convert back to matrix for prediction

      # Fit model (uses default hyperparameters for speed)
      # elasticnet: uses CV to select lambda internally
      # RF: uses default mtry and nodesize
      model_result <- fit_prediction_model(
        X_train = train_features,
        y_train = patient_outcomes$y[train_p_ids],
        X_test = test_features,
        method = prediction_method
      )

      # Calculate errors
      train_mse <- mean((model_result$train_pred - patient_outcomes$y[train_p_ids])^2)
      test_mse <- mean((model_result$test_pred - patient_outcomes$y[test_p_ids])^2)

      test_errors[fold_id, k_idx] <- test_mse
      train_errors[fold_id, k_idx] <- train_mse
    }
  }

  return(list(
    k_range = k_range,
    mean_test_error = colMeans(test_errors),
    se_test_error = apply(test_errors, 2, sd) / sqrt(n_folds),
    mean_train_error = colMeans(train_errors),
    se_train_error = apply(train_errors, 2, sd) / sqrt(n_folds)
  ))
}
