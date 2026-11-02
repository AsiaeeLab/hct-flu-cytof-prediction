# ==============================================================================
# Nested Cross-Validation Pipeline with Subsampling
# ==============================================================================
# Full nested CV: Outer loop for k selection, inner loop for hyperparameter tuning
# With multiple subsamples for robustness

library(caret)

#' Tune k inside an outer fold
#'
#' This function tunes both k and model hyperparameters using ONLY the training
#' data from a single outer fold. This prevents optimistic bias from selecting
#' k based on outer test performance.
#'
#' @param train_data Training cell-level features
#' @param train_p_id Training patient IDs for cells
#' @param patient_outcomes List with y, x, and optionally covariates
#' @param k_range Vector of k values to try
#' @param clustering_method Clustering method name
#' @param prediction_method Prediction method name
#' @param lambda_grid Grid of lambda values (for elasticnet)
#' @param mtry_grid Grid of mtry values (for RF)
#' @param nodesize_grid Grid of nodesize values (for RF)
#' @param alpha Elastic net alpha
#' @param n_inner_folds Number of inner CV folds for hyperparameter tuning
#' @return List with best k, inner CV error, and best hyperparameters
tune_k_inside_outer <- function(train_data, train_p_id, patient_outcomes,
                                k_range, clustering_method, prediction_method,
                                lambda_grid = 10^seq(-3, 1, length.out = 20),
                                mtry_grid = NULL,
                                nodesize_grid = c(1, 5, 10),
                                alpha = 0.5,
                                n_inner_folds = 5) {

  best <- list(k = NA, inner_cv_error = Inf, hyper = NULL)

  # Track inner CV validation AND training errors for ALL k values
  inner_cv_validation_errors <- numeric(length(k_range))
  inner_cv_training_errors <- numeric(length(k_range))
  names(inner_cv_validation_errors) <- k_range
  names(inner_cv_training_errors) <- k_range

  for (k_idx in seq_along(k_range)) {
    k <- k_range[k_idx]

    # Cluster on training cells only
    cl <- perform_clustering(train_data, k, method = clustering_method)
    tr_clusters <- cl$clusters

    # Get training patient IDs that have outcomes
    train_p_ids_with_outcome <- unique(train_p_id)
    train_p_ids_with_outcome <- train_p_ids_with_outcome[train_p_ids_with_outcome %in% names(patient_outcomes$y)]

    # Build patient features (with covariates if present)
    if (!is.null(patient_outcomes$covariates)) {
      covars <- patient_outcomes$covariates[train_p_ids_with_outcome, , drop = FALSE]
      Xtr <- build_patient_features_with_covariates(
        train_p_id, tr_clusters,
        patient_outcomes$x[train_p_ids_with_outcome],
        covariates = covars
      )
    } else {
      Xtr <- build_patient_features(
        train_p_id, tr_clusters,
        patient_outcomes$x[train_p_ids_with_outcome]
      )
    }

    ytr <- patient_outcomes$y[rownames(Xtr)]

    # Inner CV to tune hyperparameters for this k
    if (prediction_method == "elasticnet") {
      gs <- grid_search_elasticnet(Xtr, ytr, alpha = alpha,
                                   lambda_grid = lambda_grid,
                                   n_inner_folds = n_inner_folds)
      inner_validation_err <- min(gs$cv_errors, na.rm = TRUE)
      hyper <- list(lambda = gs$best_lambda, alpha = alpha)

      # Compute training error with best hyperparameters
      model_train <- glmnet(Xtr, ytr, alpha = hyper$alpha, lambda = hyper$lambda)
      train_pred <- predict(model_train, newx = Xtr, s = hyper$lambda)
      inner_training_err <- mean((train_pred - ytr)^2)

    } else if (prediction_method == "random_forest") {
      gs <- grid_search_rf(Xtr, ytr, mtry_grid = mtry_grid,
                          nodesize_grid = nodesize_grid, ntree = 500,
                          n_inner_folds = max(3, floor(n_inner_folds / 2)))
      inner_validation_err <- min(gs$grid_results$cv_error, na.rm = TRUE)
      hyper <- list(mtry = gs$best_mtry, nodesize = gs$best_nodesize, ntree = 500)

      # Compute training error with best hyperparameters
      train_df <- as.data.frame(Xtr)
      train_df$y <- ytr
      model_train <- randomForest(y ~ ., data = train_df,
                                  ntree = hyper$ntree,
                                  mtry = hyper$mtry,
                                  nodesize = hyper$nodesize)
      train_pred <- predict(model_train, newdata = train_df)
      inner_training_err <- mean((train_pred - ytr)^2)
    }

    # Store both errors for this k
    inner_cv_validation_errors[k_idx] <- inner_validation_err
    inner_cv_training_errors[k_idx] <- inner_training_err

    # Track best k based on validation error
    if (inner_validation_err < best$inner_cv_error) {
      best <- list(k = k, inner_cv_error = inner_validation_err, hyper = hyper)
    }
  }

  # Return best k AND all inner CV errors (both training and validation)
  best$inner_cv_validation_errors <- inner_cv_validation_errors
  best$inner_cv_training_errors <- inner_cv_training_errors
  return(best)
}


#' Nested cross-validation for a single k value
#'
#' Performs rigorous nested CV:
#' - Outer loop: 5-fold CV for performance evaluation
#' - Inner loop: Hyperparameter tuning on training folds only
#'
#' @param mini_data Cell-level feature matrix
#' @param p_id Vector mapping cells to patient IDs
#' @param patient_outcomes List with y and x outcomes (from prepare_outcome)
#' @param k Number of clusters
#' @param clustering_method Clustering method name
#' @param prediction_method Prediction method name
#' @param n_outer_folds Number of outer CV folds
#' @param lambda_grid Grid of lambda values (for elasticnet)
#' @param mtry_grid Grid of mtry values (for RF)
#' @param nodesize_grid Grid of nodesize values (for RF)
#' @param alpha Elastic net alpha (fixed at 0.5)
#' @param verbose Print progress
#' @return List with CV results and best hyperparameters
nested_cv_single_k <- function(mini_data, p_id, patient_outcomes, k,
                               clustering_method, prediction_method,
                               n_outer_folds = 5,
                               lambda_grid = 10^seq(-3, 1, length.out = 20),
                               mtry_grid = NULL,
                               nodesize_grid = c(1, 5, 10),
                               alpha = 0.5,
                               verbose = FALSE) {

  # Filter data to valid patients
  valid_indices <- p_id %in% names(patient_outcomes$y)
  mini_data <- mini_data[valid_indices, ]
  p_id <- p_id[valid_indices]

  # Create outer CV folds (patient-level)
  unique_p_ids <- names(patient_outcomes$y)
  outer_folds <- createFolds(unique_p_ids, k = n_outer_folds, list = TRUE)

  # Storage for results
  outer_fold_test_errors <- numeric(n_outer_folds)
  outer_fold_train_errors <- numeric(n_outer_folds)
  outer_fold_hyperparams <- list()

  # Outer CV loop
  for (fold_id in 1:n_outer_folds) {
    if (verbose) cat("    Outer fold", fold_id, "/", n_outer_folds, "...")

    # Split data
    test_p_ids <- unique_p_ids[outer_folds[[fold_id]]]
    train_p_ids <- setdiff(unique_p_ids, test_p_ids)

    train_idx <- which(p_id %in% train_p_ids)
    test_idx <- which(p_id %in% test_p_ids)

    train_data <- mini_data[train_idx, ]
    test_data <- mini_data[test_idx, ]
    train_p_id <- p_id[train_idx]
    test_p_id <- p_id[test_idx]

    # Cluster training data only
    cluster_result <- perform_clustering(train_data, k, method = clustering_method)
    train_clusters <- cluster_result$clusters

    # Build training features (with covariates if available)
    train_covariates <- if (!is.null(patient_outcomes$covariates)) {
      patient_outcomes$covariates[train_p_ids, , drop = FALSE]
    } else {
      NULL
    }

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
    test_clusters <- assign_to_clusters(k, train_data, train_clusters, test_data)

    # Build test features
    test_features_raw <- as.data.frame.matrix(table(test_p_id, test_clusters))
    test_features_raw <- log(test_features_raw + 1)

    # Determine which columns are clusters (not baseline or covariates)
    non_cluster_cols <- c("baseline")
    if (!is.null(train_covariates)) {
      non_cluster_cols <- c(non_cluster_cols, colnames(train_covariates))
    }
    train_cluster_cols <- setdiff(colnames(train_features), non_cluster_cols)

    test_features_aligned <- align_features(test_features_raw, train_cluster_cols)

    # Add baseline and covariates
    test_patient_ids <- rownames(test_features_aligned)
    test_features <- as.matrix(test_features_aligned)
    test_features <- cbind(test_features, baseline = patient_outcomes$x[test_patient_ids])

    if (!is.null(train_covariates)) {
      test_covariates <- patient_outcomes$covariates[test_patient_ids, , drop = FALSE]
      test_features <- cbind(test_features, test_covariates)
    }

    # INNER CV: Tune hyperparameters on training data
    # Extract y values matching the train_features rows
    train_y <- patient_outcomes$y[rownames(train_features)]

    if (prediction_method == "elasticnet") {
      # Grid search for lambda
      grid_result <- grid_search_elasticnet(train_features,
                                            train_y,
                                            alpha = alpha,
                                            lambda_grid = lambda_grid,
                                            n_inner_folds = 5)
      best_lambda <- grid_result$best_lambda
      outer_fold_hyperparams[[fold_id]] <- list(lambda = best_lambda, alpha = alpha)

      # Train with best lambda on full training set
      model <- glmnet(train_features, train_y,
                     alpha = alpha, lambda = best_lambda, intercept = TRUE)
      test_pred <- as.numeric(predict(model, newx = test_features, s = best_lambda))

    } else if (prediction_method == "random_forest") {
      # Grid search for RF hyperparameters
      grid_result <- grid_search_rf(train_features,
                                    train_y,
                                    mtry_grid = mtry_grid,
                                    nodesize_grid = nodesize_grid,
                                    ntree = 500,
                                    n_inner_folds = 3)
      best_mtry <- grid_result$best_mtry
      best_nodesize <- grid_result$best_nodesize
      outer_fold_hyperparams[[fold_id]] <- list(mtry = best_mtry,
                                                nodesize = best_nodesize,
                                                ntree = 500)

      # Train with best hyperparameters on full training set
      train_df <- as.data.frame(train_features)
      train_df$y <- train_y
      model <- randomForest(y ~ ., data = train_df,
                           ntree = 500, mtry = best_mtry,
                           nodesize = best_nodesize)
      test_pred <- predict(model, newdata = as.data.frame(test_features))
    }

    # Calculate test error for this fold
    test_y <- patient_outcomes$y[rownames(test_features)]
    outer_fold_test_errors[fold_id] <- mean((test_pred - test_y)^2)

    # Calculate training error for this fold
    train_y <- patient_outcomes$y[rownames(train_features)]
    if (prediction_method == "elasticnet") {
      train_pred <- predict(model, newx = train_features, s = best_lambda)
    } else {
      train_pred <- predict(model, newdata = as.data.frame(train_features))
    }
    outer_fold_train_errors[fold_id] <- mean((train_pred - train_y)^2)

    if (verbose) cat(" Test:", round(outer_fold_test_errors[fold_id], 4),
                     " Train:", round(outer_fold_train_errors[fold_id], 4), "\n")
  }

  # Average results across outer folds
  mean_test_error <- mean(outer_fold_test_errors)
  se_test_error <- sd(outer_fold_test_errors) / sqrt(n_outer_folds)
  mean_train_error <- mean(outer_fold_train_errors)
  se_train_error <- sd(outer_fold_train_errors) / sqrt(n_outer_folds)

  # Average hyperparameters
  if (prediction_method == "elasticnet") {
    avg_lambda <- mean(sapply(outer_fold_hyperparams, function(x) x$lambda))
    best_hyperparams <- list(lambda = avg_lambda, alpha = alpha)
  } else if (prediction_method == "random_forest") {
    avg_mtry <- round(mean(sapply(outer_fold_hyperparams, function(x) x$mtry)))
    avg_nodesize <- round(mean(sapply(outer_fold_hyperparams, function(x) x$nodesize)))
    best_hyperparams <- list(mtry = avg_mtry, nodesize = avg_nodesize, ntree = 500)
  }

  return(list(
    mean_test_error = mean_test_error,
    se_test_error = se_test_error,
    mean_train_error = mean_train_error,
    se_train_error = se_train_error,
    outer_fold_test_errors = outer_fold_test_errors,
    outer_fold_train_errors = outer_fold_train_errors,
    best_hyperparams = best_hyperparams,
    hyperparams_per_fold = outer_fold_hyperparams
  ))
}


#' Full nested CV with subsampling
#'
#' Complete pipeline:
#' - Subsample cells n_subsamples times
#' - For each subsample, test all k values with nested CV
#' - Aggregate results across subsamples
#'
#' @param full_data Full cell-level feature matrix (all 18M cells or subset)
#' @param full_p_id Vector mapping all cells to patient IDs
#' @param outcome_obj Outcome data frame
#' @param meta_obj Metadata data frame with covariates (optional)
#' @param outcome_col_1 Baseline outcome column
#' @param outcome_col_2 Target outcome column
#' @param clustering_methods Vector of clustering methods to test
#' @param prediction_methods Vector of prediction methods to test
#' @param k_range Vector of k values to test
#' @param n_cells Number of cells to subsample
#' @param n_subsamples Number of subsampling iterations
#' @param n_outer_folds Number of outer CV folds
#' @param lambda_grid Lambda grid for elasticnet
#' @param alpha Elasticnet alpha (fixed)
#' @param include_dose Include dose covariate
#' @param include_age Include age covariate
#' @param include_sex Include sex covariate
#' @param verbose Print progress
#' @return Comprehensive results structure
nested_cv_with_subsampling <- function(full_data, full_p_id, outcome_obj,
                                       meta_obj = NULL,
                                       outcome_col_1, outcome_col_2,
                                       clustering_methods = c("hierarchical", "kmeans"),
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
                                       verbose = TRUE) {

  if (verbose) {
    cat("\n", rep("=", 70), "\n", sep = "")
    cat("NESTED CV WITH SUBSAMPLING\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("Parameters:\n")
    cat("  - n_cells per subsample:", n_cells, "\n")
    cat("  - n_subsamples:", n_subsamples, "\n")
    cat("  - k_range:", paste(k_range, collapse = ", "), "\n")
    cat("  - Clustering methods:", paste(clustering_methods, collapse = ", "), "\n")
    cat("  - Prediction methods:", paste(prediction_methods, collapse = ", "), "\n")
    cat("  - Outer CV folds:", n_outer_folds, "\n")
    cat("  - Covariates: ")
    covars <- c()
    if (include_dose) covars <- c(covars, "dose")
    if (include_age) covars <- c(covars, "age")
    if (include_sex) covars <- c(covars, "sex")
    if (include_time) covars <- c(covars, "time_post_transplant")
    cat(paste(covars, collapse = ", "), "\n\n")
  }

  # Prepare outcomes (once, before subsampling)
  if (!is.null(meta_obj)) {
    patient_outcomes <- prepare_outcome_with_covariates(
      outcome_obj, meta_obj, outcome_col_1, outcome_col_2, full_p_id,
      include_dose = include_dose,
      include_age = include_age,
      include_sex = include_sex,
      include_time = include_time
    )
  } else {
    patient_outcomes <- prepare_outcome(outcome_obj, outcome_col_1, outcome_col_2, full_p_id)
  }

  # Initialize results structure
  results <- list()

  # Loop over subsamples
  for (subsample_idx in 1:n_subsamples) {
    if (verbose) {
      cat("\n", rep("-", 70), "\n", sep = "")
      cat("SUBSAMPLE ", subsample_idx, " / ", n_subsamples, "\n", sep = "")
      cat(rep("-", 70), "\n", sep = "")
    }

    # Subsample cells
    set.seed(123 + subsample_idx)
    if (nrow(full_data) > n_cells) {
      sample_idx <- sample(nrow(full_data), n_cells)
      mini_data <- full_data[sample_idx, ]
      p_id <- full_p_id[sample_idx]
    } else {
      mini_data <- full_data
      p_id <- full_p_id
    }

    if (verbose) cat("Sampled", nrow(mini_data), "cells\n\n")

    # Loop over method combinations
    for (clust_method in clustering_methods) {
      for (pred_method in prediction_methods) {
        method_name <- paste(clust_method, pred_method, sep = "_")

        if (verbose) {
          cat("Method:", method_name, "\n")
        }

        # Storage for this method
        k_results <- list()

        # Loop over k values
        for (k_idx in seq_along(k_range)) {
          k <- k_range[k_idx]
          if (verbose) cat("  k =", k, "\n")

          # Nested CV for this k
          cv_result <- nested_cv_single_k(
            mini_data, p_id, patient_outcomes, k,
            clustering_method = clust_method,
            prediction_method = pred_method,
            n_outer_folds = n_outer_folds,
            lambda_grid = lambda_grid,
            alpha = alpha,
            verbose = verbose
          )

          k_results[[paste0("k", k)]] <- cv_result
        }

        # Store results for this method and subsample
        if (is.null(results[[method_name]])) {
          results[[method_name]] <- list()
        }
        results[[method_name]][[paste0("subsample", subsample_idx)]] <- k_results
      }
    }
  }

  # Aggregate results across subsamples
  aggregated <- aggregate_subsample_results(results, k_range,
                                           clustering_methods, prediction_methods,
                                           n_subsamples)

  return(list(
    raw_results = results,
    aggregated = aggregated,
    parameters = list(
      k_range = k_range,
      n_cells = n_cells,
      n_subsamples = n_subsamples,
      clustering_methods = clustering_methods,
      prediction_methods = prediction_methods,
      n_outer_folds = n_outer_folds,
      alpha = alpha
    )
  ))
}


#' Aggregate results across subsamples
#'
#' @param results Raw results from nested_cv_with_subsampling
#' @param k_range Vector of k values
#' @param clustering_methods Vector of clustering methods
#' @param prediction_methods Vector of prediction methods
#' @param n_subsamples Number of subsamples
#' @return Aggregated results
aggregate_subsample_results <- function(results, k_range, clustering_methods,
                                        prediction_methods, n_subsamples) {

  method_combinations <- expand.grid(
    clustering = clustering_methods,
    prediction = prediction_methods,
    stringsAsFactors = FALSE
  )
  method_combinations$method_name <- paste(method_combinations$clustering,
                                          method_combinations$prediction, sep = "_")

  aggregated <- list()

  for (i in 1:nrow(method_combinations)) {
    method_name <- method_combinations$method_name[i]
    method_data <- results[[method_name]]

    # Extract errors for each k across subsamples
    mean_test_errors <- matrix(NA, nrow = n_subsamples, ncol = length(k_range))
    se_test_errors <- matrix(NA, nrow = n_subsamples, ncol = length(k_range))
    mean_train_errors <- matrix(NA, nrow = n_subsamples, ncol = length(k_range))
    se_train_errors <- matrix(NA, nrow = n_subsamples, ncol = length(k_range))

    for (subsample_idx in 1:n_subsamples) {
      subsample_data <- method_data[[paste0("subsample", subsample_idx)]]

      for (k_idx in seq_along(k_range)) {
        k <- k_range[k_idx]
        k_data <- subsample_data[[paste0("k", k)]]

        mean_test_errors[subsample_idx, k_idx] <- k_data$mean_test_error
        se_test_errors[subsample_idx, k_idx] <- k_data$se_test_error
        mean_train_errors[subsample_idx, k_idx] <- k_data$mean_train_error
        se_train_errors[subsample_idx, k_idx] <- k_data$se_train_error
      }
    }

    # Average across subsamples
    avg_mean_test_error <- colMeans(mean_test_errors)
    avg_se_test_error <- apply(mean_test_errors, 2, sd) / sqrt(n_subsamples)
    avg_mean_train_error <- colMeans(mean_train_errors)
    avg_se_train_error <- apply(mean_train_errors, 2, sd) / sqrt(n_subsamples)

    # Find best k (based on test error)
    best_k_idx <- which.min(avg_mean_test_error)
    best_k <- k_range[best_k_idx]

    aggregated[[method_name]] <- list(
      clustering = method_combinations$clustering[i],
      prediction = method_combinations$prediction[i],
      mean_test_error = avg_mean_test_error,
      se_test_error = avg_se_test_error,
      mean_train_error = avg_mean_train_error,
      se_train_error = avg_se_train_error,
      test_error_matrix = mean_test_errors,  # n_subsamples × n_k
      train_error_matrix = mean_train_errors,  # n_subsamples × n_k
      best_k = best_k,
      best_k_error = avg_mean_test_error[best_k_idx]
    )
  }

  return(aggregated)
}
