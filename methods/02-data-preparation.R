# ==============================================================================
# Data Preparation Functions
# ==============================================================================
# This file contains functions for loading, preparing, and transforming data
# for the clustering and prediction pipeline.

ensure_p_id_yr <- function(outcome_obj) {
  if (!is.null(outcome_obj$p_id_yr)) {
    return(outcome_obj)
  }

  if (!is.null(outcome_obj$id_yr)) {
    outcome_obj$p_id_yr <- sub("_(\\d)$", "_yr\\1", outcome_obj$id_yr)
    return(outcome_obj)
  }

  if (!is.null(outcome_obj$pid) && !is.null(outcome_obj$yr)) {
    outcome_obj$p_id_yr <- paste0(outcome_obj$pid, "_yr", outcome_obj$yr)
    return(outcome_obj)
  }

  stop("Cannot construct p_id_yr: need one of {p_id_yr}, {id_yr}, or {pid, yr} columns.")
}

#' Prepare outcome variables with log transformation
#'
#' Extracts baseline and target outcomes, log-transforms them, and filters for
#' complete cases. The response variable is calculated as log(v3/v1).
#'
#' @param outcome_obj Data frame containing outcome measurements
#' @param outcome_col_1 Column name for baseline measurement (v1)
#' @param outcome_col_2 Column name for target measurement (v3)
#' @param p_id Vector of patient IDs from the cell data
#' @return List with y (response: log(v3/v1)), x (baseline: log(v1)), and filtered patient IDs
prepare_outcome <- function(outcome_obj, outcome_col_1, outcome_col_2, p_id) {
  outcome_obj <- ensure_p_id_yr(outcome_obj)

  # Get both outcomes
  outcome_v1 <- outcome_obj[[outcome_col_1]]  # Baseline (visit 1)
  outcome_v3 <- outcome_obj[[outcome_col_2]]  # Target (visit 3)

  # Calculate response as log(v3/v1) = log(v3) - log(v1)
  # This represents fold-change on the log scale
  outcome_y <- log(outcome_v3) - log(outcome_v1)
  outcome_x <- log(outcome_v1)  # Baseline covariate (log-transformed)

  # Check for completeness - both must be non-NA and finite
  has_both <- !is.na(outcome_y) & !is.na(outcome_x) &
              is.finite(outcome_x) & is.finite(outcome_y)
  has_cell <- outcome_obj$p_id_yr %in% p_id
  complete_index <- has_both & has_cell

  # Return both outcomes for complete cases
  outcome_y <- outcome_y[complete_index]
  outcome_x <- outcome_x[complete_index]
  patient_ids <- outcome_obj$p_id_yr[complete_index]

  patient_outcome_y <- setNames(outcome_y, patient_ids)
  patient_outcome_x <- setNames(outcome_x, patient_ids)

  return(list(y = patient_outcome_y, x = patient_outcome_x))
}


#' Prepare outcome variables - ALTERNATIVE: predict log(v3) using log(v1) as feature
#'
#' This is an alternative formulation where:
#' - Response (y) = log(v3)
#' - Baseline feature (x) = log(v1) [used as a regular feature, not subtracted]
#'
#' This allows the model to learn the relationship between baseline and outcome
#' rather than forcing it to be a fold-change.
#'
#' @param outcome_obj Data frame containing outcome measurements
#' @param outcome_col_1 Column name for baseline measurement (v1)
#' @param outcome_col_2 Column name for target measurement (v3)
#' @param p_id Vector of patient IDs from the cell data
#' @return List with y (response: log(v3)), x (feature: log(v1)), and filtered patient IDs
prepare_outcome_absolute <- function(outcome_obj, outcome_col_1, outcome_col_2, p_id) {
  outcome_obj <- ensure_p_id_yr(outcome_obj)

  # Get both outcomes
  outcome_v1 <- outcome_obj[[outcome_col_1]]  # Baseline (visit 1)
  outcome_v3 <- outcome_obj[[outcome_col_2]]  # Target (visit 3)

  # Response is log(v3) directly
  outcome_y <- log(outcome_v3)
  # Baseline is log(v1) - will be used as a feature
  outcome_x <- log(outcome_v1)

  # Check for completeness
  has_both <- !is.na(outcome_y) & !is.na(outcome_x) &
              is.finite(outcome_x) & is.finite(outcome_y)
  has_cell <- outcome_obj$p_id_yr %in% p_id
  complete_index <- has_both & has_cell

  # Return outcomes for complete cases
  outcome_y <- outcome_y[complete_index]
  outcome_x <- outcome_x[complete_index]
  patient_ids <- outcome_obj$p_id_yr[complete_index]

  patient_outcome_y <- setNames(outcome_y, patient_ids)
  patient_outcome_x <- setNames(outcome_x, patient_ids)

  return(list(y = patient_outcome_y, x = patient_outcome_x))
}


#' Split data into training and test sets for cross-validation
#'
#' Splits patients (not cells) into train/test based on CV folds to avoid
#' data leakage from the nested structure.
#'
#' @param fold_id Current fold index
#' @param folds List of fold indices from createFolds
#' @param unique_p_ids Unique patient IDs
#' @param p_id Vector mapping cells to patients
#' @param mini_data Cell-level feature matrix
#' @param patient_outcomes List containing y and x outcomes
#' @return List with train/test splits for data, responses, and patient IDs
split_data <- function(fold_id, folds, unique_p_ids, p_id, mini_data, patient_outcomes) {
  test_unique_ids <- unique_p_ids[folds[[fold_id]]]
  train_unique_ids <- setdiff(unique_p_ids, test_unique_ids)

  train_indices <- which(p_id %in% train_unique_ids)
  test_indices <- which(p_id %in% test_unique_ids)

  train_data <- mini_data[train_indices, ]
  test_data <- mini_data[test_indices, ]

  # Split both y and x
  train_response_y <- patient_outcomes$y[train_unique_ids]
  test_response_y <- patient_outcomes$y[test_unique_ids]
  train_response_x <- patient_outcomes$x[train_unique_ids]
  test_response_x <- patient_outcomes$x[test_unique_ids]

  train_p_id <- p_id[train_indices]
  test_p_id <- p_id[test_indices]

  list(
    train_data = train_data, test_data = test_data,
    train_response_y = train_response_y, test_response_y = test_response_y,
    train_response_x = train_response_x, test_response_x = test_response_x,
    train_p_id = train_p_id, test_p_id = test_p_id
  )
}


#' Build patient-level feature matrix from cell clusters
#'
#' Aggregates cell-level cluster assignments into patient-level features by
#' counting cells per cluster for each patient, then log-transforming for
#' variance stabilization.
#'
#' @param p_id Vector of patient IDs for each cell
#' @param clusters Vector of cluster assignments for each cell
#' @param baseline_outcome Named vector of baseline outcomes (log(v1)) per patient
#' @param cluster_names Optional vector of cluster names/labels
#' @return Matrix with patients as rows, clusters as columns, plus baseline
build_patient_features <- function(p_id, clusters, baseline_outcome, cluster_names = NULL) {
  # Build patient × cluster count matrix
  feature_df <- as.data.frame.matrix(table(p_id, clusters))

  # Log-transform for variance stabilization: log(count + 1)
  feature_df <- log(feature_df + 1)

  # Rename columns if cluster names provided
  if (!is.null(cluster_names)) {
    colnames(feature_df) <- cluster_names
  } else {
    colnames(feature_df) <- paste0("cluster_", 1:ncol(feature_df))
  }

  # Convert to matrix and add baseline covariate
  feature_matrix <- as.matrix(feature_df)
  patient_ids <- rownames(feature_matrix)
  feature_matrix_with_baseline <- cbind(feature_matrix, baseline = baseline_outcome[patient_ids])

  return(feature_matrix_with_baseline)
}


#' Align test features with training features
#'
#' Ensures test feature matrix has same columns as training, adding zeros for
#' missing clusters.
#'
#' @param test_features Test feature data frame
#' @param train_colnames Column names from training features
#' @return Aligned test feature data frame
align_features <- function(test_features, train_colnames) {
  # Add missing columns with zeros
  missing_columns <- setdiff(train_colnames, colnames(test_features))
  for (col in missing_columns) {
    test_features[[col]] <- 0
  }

  # Reorder to match training
  test_features <- test_features[, train_colnames, drop = FALSE]

  return(test_features)
}
