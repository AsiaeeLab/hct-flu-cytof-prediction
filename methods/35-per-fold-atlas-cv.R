# ==============================================================================
# Per-Fold Atlas Building for Nested CV
# ==============================================================================
# This module implements atlas building INSIDE each outer CV fold to avoid
# data leakage. The atlas (micro-cluster centroids + scaling) is learned only
# on training patients, then used to transform both train and test cells.
#
# Key difference from 11-atlas-build-cache-multivirus.R:
#   - Atlas built per-fold, not once globally
#   - No test patient cells ever seen during atlas construction
#   - Proper generalization estimation
#
# Functions:
#   - build_atlas_for_fold(): Build atlas from training patient cells only
#   - run_nested_cv_perfold_atlas(): Full nested CV with per-fold atlas
#
# Design for computational efficiency:
#   - Sample cells for k-means fitting (default: 200k cells)
#   - Assign ALL cells to centroids in chunks
#   - Parallelization ready (outer folds can run in parallel)
# ==============================================================================

suppressPackageStartupMessages({
  library(RANN)
  library(caret)
  library(glmnet)
})

# Source two-stage clustering helpers if not already loaded
if (!exists("two_stage_build_centroid_knn_graph")) {
  source(here::here("methods", "10-two-stage-clustering.R"))
}

base_pid_from_p_id_yr <- function(p_id_yr) {
  sub("_yr[0-9]+$", "", as.character(p_id_yr))
}

infer_cohort_from_p_id <- function(p_id_yr) {
  base <- base_pid_from_p_id_yr(p_id_yr)
  prefix <- sub("_.*$", "", base)
  out <- ifelse(prefix %in% c("ASCT", "ASOT"), "adult",
                ifelse(prefix == "PSCT", "peds", "unknown"))
  names(out) <- as.character(p_id_yr)
  out
}

create_group_folds <- function(p_ids,
                               group_ids,
                               k,
                               seed,
                               strata = NULL) {
  stopifnot(length(p_ids) > 0)
  if (is.null(names(group_ids))) {
    stop("group_ids must be a named vector keyed by p_id.")
  }
  group_ids <- as.character(group_ids[p_ids])
  if (any(is.na(group_ids))) stop("group_ids missing for some p_ids.")

  groups <- unique(group_ids)
  if (length(groups) < k) {
    stop("Too few groups (", length(groups), ") for k=", k, " folds.")
  }

  if (!is.null(strata)) {
    if (is.null(names(strata))) names(strata) <- p_ids
    strata <- strata[p_ids]
    group_strata <- tapply(strata, group_ids, function(x) x[which(!is.na(x))[1]])
    group_strata <- factor(group_strata[groups])
    set.seed(seed)
    fold_groups <- caret::createFolds(group_strata, k = k, list = TRUE)
  } else {
    set.seed(seed)
    fold_groups <- caret::createFolds(seq_along(groups), k = k, list = TRUE)
  }

  lapply(fold_groups, function(g_idx) {
    test_groups <- groups[g_idx]
    which(group_ids %in% test_groups)
  })
}

as_pred_matrix <- function(pred, colnames_out) {
  if (length(dim(pred)) == 3L) {
    pred <- pred[, , 1, drop = FALSE]
    dim(pred) <- dim(pred)[1:2]
  }
  pred <- as.matrix(pred)
  colnames(pred) <- colnames_out
  pred
}

fit_predict_multitask_en <- function(X_train,
                                     Y_train,
                                     X_test,
                                     seed,
                                     alpha = 0.5,
                                     n_cores = 1,
                                     foldid = NULL,
                                     lambda_rule = c("lambda.min", "lambda.1se")) {
  lambda_rule <- match.arg(lambda_rule)
  set.seed(seed)
  cv_args <- list(
    x = X_train,
    y = Y_train,
    family = "mgaussian",
    alpha = alpha,
    standardize = TRUE,
    parallel = n_cores > 1
  )
  if (!is.null(foldid)) {
    cv_args$foldid <- foldid
  } else {
    cv_args$nfolds <- 5L
  }
  cvfit <- do.call(cv.glmnet, cv_args)
  s_val <- if (lambda_rule == "lambda.min") cvfit$lambda.min else cvfit$lambda.1se
  pred_train <- as_pred_matrix(
    predict(cvfit, newx = X_train, s = s_val),
    colnames(Y_train)
  )
  pred_test <- as_pred_matrix(
    predict(cvfit, newx = X_test, s = s_val),
    colnames(Y_train)
  )
  list(pred_train = pred_train, pred_test = pred_test, model = cvfit)
}

# ==============================================================================
# Helper functions
# ==============================================================================

#' Transform meta-cluster counts to features
#'
#' C4-LOGCOUNTS fix: Supports log(count+1) (legacy), proportions, and CLR.
#'
#' @param meta_counts Patient x meta-cluster count matrix (raw integer counts)
#' @param method One of "log_count", "proportion", "clr"
#' @return Transformed feature matrix (same dimensions)
transform_counts <- function(meta_counts, method = c("log_count", "proportion", "clr")) {
  method <- match.arg(method)
  if (method == "log_count") {
    return(as.matrix(log(meta_counts + 1)))
  }
  if (method == "proportion") {
    row_sums <- rowSums(meta_counts)
    row_sums[row_sums == 0] <- 1  # avoid division by zero
    return(as.matrix(meta_counts / row_sums))
  }
  if (method == "clr") {
    # Centered log-ratio: log(x_i / geometric_mean(x))
    props <- meta_counts + 0.5  # pseudocount for zeros
    log_props <- log(props)
    geo_means <- rowMeans(log_props)
    return(as.matrix(log_props - geo_means))
  }
}

#' OLS baseline model (single response)
#'
#' C2-BASELINE fix: Deterministic baseline — no seed, no lambda, no cv.glmnet.
#' Appropriate for 4-5 covariates (dose, age, sex, time, baseline titer).
#'
#' @param X_train Training covariate matrix
#' @param y_train Training response vector
#' @param X_test Test covariate matrix
#' @return List with pred_train, pred_test, model
fit_ols_baseline <- function(X_train, y_train, X_test) {
  df_train <- data.frame(X_train, y = y_train, check.names = FALSE)
  model <- lm(y ~ ., data = df_train)
  pred_train <- predict(model, newdata = data.frame(X_train, check.names = FALSE))
  pred_test <- predict(model, newdata = data.frame(X_test, check.names = FALSE))
  list(pred_train = as.numeric(pred_train),
       pred_test = as.numeric(pred_test),
       model = model)
}

#' OLS baseline model (multi-response, separate OLS per virus)
#'
#' C2-BASELINE fix for multitask: Fits independent OLS per virus column.
#'
#' @param X_train Training covariate matrix
#' @param Y_train Training response matrix (n x n_viruses)
#' @param X_test Test covariate matrix
#' @return List with pred_train, pred_test (matrices with virus columns)
fit_ols_baseline_multi <- function(X_train, Y_train, X_test) {
  virus_names <- colnames(Y_train)
  pred_train <- matrix(NA_real_, nrow = nrow(X_train), ncol = length(virus_names))
  pred_test <- matrix(NA_real_, nrow = nrow(X_test), ncol = length(virus_names))
  colnames(pred_train) <- colnames(pred_test) <- virus_names
  rownames(pred_train) <- rownames(X_train)
  rownames(pred_test) <- rownames(X_test)

  for (v in virus_names) {
    ols <- fit_ols_baseline(X_train, Y_train[, v], X_test)
    pred_train[, v] <- ols$pred_train
    pred_test[, v] <- ols$pred_test
  }
  list(pred_train = pred_train, pred_test = pred_test)
}

#' Compute grouped foldid for cv.glmnet lambda CV
#'
#' M3-LAMBDA-GROUP fix: Ensures yr1 and yr2 of the same patient stay
#' in the same lambda-CV fold.
#'
#' @param patient_ids Character vector of p_id_yr identifiers
#' @param n_folds Number of CV folds (default 5)
#' @param seed Random seed for reproducible fold assignment
#' @return Integer vector of fold IDs (same length as patient_ids)
compute_grouped_foldid <- function(patient_ids, n_folds = 5L, seed = NULL) {
  base_pids <- base_pid_from_p_id_yr(patient_ids)
  unique_bases <- unique(base_pids)
  if (length(unique_bases) < n_folds) {
    n_folds <- length(unique_bases)
  }
  if (!is.null(seed)) set.seed(seed)
  group_folds <- sample(rep(seq_len(n_folds), length.out = length(unique_bases)))
  names(group_folds) <- unique_bases
  as.integer(group_folds[base_pids])
}

# ==============================================================================

build_ratio_feature_matrices <- function(train_meta_counts,
                                         test_meta_counts = NULL,
                                         max_ratio_features = 561,
                                         selector = c("variance", "all")) {
  selector <- match.arg(selector)

  train_log <- log(train_meta_counts + 1)
  test_log <- if (is.null(test_meta_counts)) NULL else log(test_meta_counts + 1)
  k <- ncol(train_log)
  if (k < 2) {
    empty_train <- matrix(numeric(0), nrow = nrow(train_log), ncol = 0)
    rownames(empty_train) <- rownames(train_log)
    empty_test <- NULL
    if (!is.null(test_log)) {
      empty_test <- matrix(numeric(0), nrow = nrow(test_log), ncol = 0)
      rownames(empty_test) <- rownames(test_log)
    }
    return(list(train = empty_train, test = empty_test, pairs = matrix(numeric(0), ncol = 2)))
  }

  pairs <- t(utils::combn(k, 2))
  n_pairs <- nrow(pairs)
  keep_idx <- seq_len(n_pairs)

  if (selector == "variance" && is.finite(max_ratio_features) && max_ratio_features < n_pairs) {
    var_scores <- apply(pairs, 1, function(idx) {
      stats::var(train_log[, idx[1]] - train_log[, idx[2]], na.rm = TRUE)
    })
    keep_idx <- order(var_scores, decreasing = TRUE)[seq_len(max_ratio_features)]
  }

  if (selector == "all" && is.finite(max_ratio_features) && max_ratio_features < n_pairs) {
    keep_idx <- seq_len(max_ratio_features)
  }

  keep_pairs <- pairs[sort(keep_idx), , drop = FALSE]
  if (nrow(keep_pairs) == 0) {
    empty_train <- matrix(numeric(0), nrow = nrow(train_log), ncol = 0)
    rownames(empty_train) <- rownames(train_log)
    empty_test <- NULL
    if (!is.null(test_log)) {
      empty_test <- matrix(numeric(0), nrow = nrow(test_log), ncol = 0)
      rownames(empty_test) <- rownames(test_log)
    }
    return(list(train = empty_train, test = empty_test, pairs = keep_pairs))
  }

  ratio_names <- paste0("ratio_c", keep_pairs[, 1], "_c", keep_pairs[, 2])
  ratio_train <- sapply(seq_len(nrow(keep_pairs)), function(i) {
    idx <- keep_pairs[i, ]
    train_log[, idx[1]] - train_log[, idx[2]]
  })
  ratio_train <- as.matrix(ratio_train)
  colnames(ratio_train) <- ratio_names
  rownames(ratio_train) <- rownames(train_log)

  ratio_test <- NULL
  if (!is.null(test_log)) {
    ratio_test <- sapply(seq_len(nrow(keep_pairs)), function(i) {
      idx <- keep_pairs[i, ]
      test_log[, idx[1]] - test_log[, idx[2]]
    })
    ratio_test <- as.matrix(ratio_test)
    colnames(ratio_test) <- ratio_names
    rownames(ratio_test) <- rownames(test_log)
  }

  list(train = ratio_train, test = ratio_test, pairs = keep_pairs)
}

compute_patient_micro_counts_from_assignments <- function(p_ids,
                                                          micro_assignments,
                                                          patient_ids,
                                                          n_micro) {
  stopifnot(length(p_ids) == length(micro_assignments))
  patient_ids <- as.character(patient_ids)
  micro_assignments <- as.integer(micro_assignments)
  p_ids <- as.character(p_ids)

  out <- matrix(0L, nrow = length(patient_ids), ncol = n_micro)
  rownames(out) <- patient_ids

  ord <- order(p_ids)
  p_sorted <- p_ids[ord]
  micro_sorted <- micro_assignments[ord]
  r <- rle(p_sorted)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  pid_to_row <- setNames(seq_len(length(patient_ids)), patient_ids)

  for (i in seq_along(r$values)) {
    pid <- r$values[i]
    row <- pid_to_row[[pid]]
    if (is.null(row)) next
    seg <- micro_sorted[starts[i]:ends[i]]
    out[row, ] <- tabulate(seg, nbins = n_micro)
  }

  out
}


#' Build atlas for a single fold (training patients only)
#'
#' Constructs micro-cluster centroids and scaling parameters from cells
#' belonging to training patients only. Test patient cells are NEVER seen.
#'
#' @param full_data Full cell-level matrix (cells x markers)
#' @param full_p_id Patient ID vector (length = nrow(full_data))
#' @param train_patient_ids Character vector of training patient IDs
#' @param n_micro Number of micro-clusters (k-means centers)
#' @param n_neighbors Number of neighbors for centroid kNN graph
#' @param kmeans_fit_n_cells Max cells to sample for k-means fitting
#' @param kmeans_nstart Number of random starts for k-means
#' @param kmeans_iter_max Max iterations for k-means
#' @param seed Random seed for reproducibility
#' @param verbose Print progress
#' @return List with:
#'   - centroids: micro-cluster centroids (n_micro x p)
#'   - scale_center: scaling center (length p)
#'   - scale_scale: scaling scale (length p)
#'   - centroid_graph: igraph kNN graph on centroids
#'   - n_train_cells: number of training cells used
#'   - timing_s: time taken in seconds
build_atlas_for_fold <- function(full_data,
                                  full_p_id,
                                  train_patient_ids,
                                  n_micro = 5000,
                                  n_neighbors = 30,
                                  kmeans_fit_n_cells = 200000,
                                  kmeans_nstart = 1,
                                  kmeans_iter_max = 50,
                                  seed = 123,
                                  verbose = FALSE) {

  tic <- Sys.time()

  # Identify training cells
  train_cell_idx <- which(full_p_id %in% train_patient_ids)
  n_train_cells <- length(train_cell_idx)

  if (n_train_cells == 0) {
    stop("No cells found for training patients.")
  }

  if (verbose) {
    cat(sprintf("    [build_atlas_for_fold] Training cells: %d from %d patients\n",
                n_train_cells, length(train_patient_ids)))
  }

  # Sample cells for k-means fitting (patient-stratified)
  set.seed(seed)
  fit_n <- min(kmeans_fit_n_cells, n_train_cells)

  if (fit_n < n_train_cells) {
    # Patient-stratified sampling
    train_p_id_subset <- full_p_id[train_cell_idx]
    by_patient <- split(train_cell_idx, train_p_id_subset)
    target_per_patient <- ceiling(fit_n / length(by_patient))

    fit_idx <- unlist(lapply(by_patient, function(ix) {
      sample(ix, size = min(length(ix), target_per_patient))
    }))

    if (length(fit_idx) > fit_n) {
      fit_idx <- sample(fit_idx, size = fit_n)
    }
  } else {
    fit_idx <- train_cell_idx
  }

  if (verbose) {
    cat(sprintf("    [build_atlas_for_fold] Fitting k-means on %d cells\n", length(fit_idx)))
  }

  # Compute scaling parameters from fitting sample
  X_fit <- full_data[fit_idx, , drop = FALSE]
  scale_center <- colMeans(X_fit)
  scale_scale <- apply(X_fit, 2, sd)
  scale_scale[scale_scale == 0] <- 1  # Avoid division by zero

  X_fit_scaled <- scale(X_fit, center = scale_center, scale = scale_scale)

  # K-means clustering
  set.seed(seed + 1)
  km <- kmeans(
    X_fit_scaled,
    centers = n_micro,
    nstart = kmeans_nstart,
    iter.max = kmeans_iter_max
  )
  centroids <- km$centers
  fit_micro_assignments <- km$cluster
  fit_p_id <- full_p_id[fit_idx]

  if (verbose) {
    cat(sprintf("    [build_atlas_for_fold] K-means complete: %d centroids\n", nrow(centroids)))
  }

  # Build centroid kNN graph for Leiden clustering
  centroid_graph <- two_stage_build_centroid_knn_graph(
    centroids = centroids,
    n_neighbors = n_neighbors
  )

  timing_s <- as.numeric(difftime(Sys.time(), tic, units = "secs"))

  if (verbose) {
    cat(sprintf("    [build_atlas_for_fold] Complete in %.1f seconds\n", timing_s))
  }

  list(
    centroids = centroids,
    scale_center = scale_center,
    scale_scale = scale_scale,
    centroid_graph = centroid_graph,
    n_train_cells = n_train_cells,
    n_fit_cells = length(fit_idx),
    fit_p_id = fit_p_id,
    fit_micro_assignments = fit_micro_assignments,
    timing_s = timing_s
  )
}


#' Compute patient x micro-cluster counts for a set of patients
#'
#' Assigns all cells from specified patients to pre-built centroids and
#' returns aggregated counts per patient.
#'
#' @param full_data Full cell-level matrix
#' @param full_p_id Patient ID vector
#' @param patient_ids Character vector of patient IDs to process
#' @param centroids Micro-cluster centroids (n_micro x p)
#' @param scale_center Scaling center
#' @param scale_scale Scaling scale
#' @param chunk_size Cells per chunk for memory efficiency
#' @param verbose Print progress
#' @return Integer matrix (patients x n_micro) with micro-cluster counts
compute_patient_micro_counts <- function(full_data,
                                          full_p_id,
                                          patient_ids,
                                          centroids,
                                          scale_center,
                                          scale_scale,
                                          chunk_size = 300000,
                                          max_cells_per_patient = NULL,
                                          seed = NULL,
                                          verbose = FALSE) {

  # Use the existing helper from 10-two-stage-clustering.R
  two_stage_patient_micro_counts_full(
    full_data = full_data,
    full_p_id = full_p_id,
    patient_ids = patient_ids,
    centroids = centroids,
    scale_center = scale_center,
    scale_scale = scale_scale,
    chunk_size = chunk_size,
    max_cells_per_patient = max_cells_per_patient,
    seed = seed,
    verbose = verbose
  )
}


#' Compute targeted marker summary features (per-patient per-meta-cluster means)
#'
#' For each meta-cluster, selects the top N most variable markers (by cross-patient
#' variance on training data) and computes per-patient mean marker expression.
#' Uses chunked cell processing for memory efficiency.
#'
#' @param full_data Full cell-level marker matrix (cells x markers)
#' @param full_p_id Patient ID per cell
#' @param patient_ids Character vector of patients to compute features for
#' @param centroids Micro-cluster centroids (n_micro x p)
#' @param scale_center Atlas scaling center
#' @param scale_scale Atlas scaling scale
#' @param micro_to_meta Named integer vector mapping micro-cluster -> meta-cluster
#' @param n_markers_per_cluster Top N markers to keep per meta-cluster
#' @param selection_patient_ids Patient IDs used for variance-based marker selection
#'        (should be training patients only to avoid leakage)
#' @param selected_markers Optional pre-computed marker selection (list of character
#'        vectors per meta-cluster). If provided, skips variance-based selection.
#' @param chunk_size Cells per chunk for RANN assignment
#' @param verbose Print progress
#' @return List with:
#'   - features: Matrix (patients x selected_features), columns named marker_<name>_c<cluster>
#'   - selected_markers: Named list (meta-cluster -> character vector of marker names)
#'   - n_features: Total number of marker features
compute_targeted_marker_features <- function(full_data,
                                              full_p_id,
                                              patient_ids,
                                              centroids,
                                              scale_center,
                                              scale_scale,
                                              micro_to_meta,
                                              n_markers_per_cluster = 5L,
                                              selection_patient_ids = NULL,
                                              selected_markers = NULL,
                                              chunk_size = 300000,
                                              verbose = FALSE) {

  patient_ids <- as.character(patient_ids)
  n_patients <- length(patient_ids)
  marker_names <- colnames(full_data)
  n_markers <- ncol(full_data)
  n_micro <- nrow(centroids)

  meta_ids <- sort(unique(micro_to_meta))
  n_meta <- length(meta_ids)

  # Pre-allocate accumulator arrays: patients x meta-clusters x markers
  marker_sums <- array(0, dim = c(n_patients, n_meta, n_markers))
  cell_counts <- matrix(0L, nrow = n_patients, ncol = n_meta)

  pid_to_row <- setNames(seq_along(patient_ids), patient_ids)
  meta_to_col <- setNames(seq_along(meta_ids), as.character(meta_ids))

  # Find cells belonging to requested patients
  idx <- which(full_p_id %in% patient_ids)
  if (length(idx) == 0) stop("No cells found for requested patient_ids in marker features.")

  row_idx_all <- pid_to_row[full_p_id[idx]]

  # Process cells in chunks (same pattern as two_stage_patient_micro_counts_full)
  n_chunks <- ceiling(length(idx) / chunk_size)
  for (chunk in seq_len(n_chunks)) {
    start <- (chunk - 1) * chunk_size + 1
    end <- min(chunk * chunk_size, length(idx))
    sel <- start:end

    idx_chunk <- idx[sel]
    row_idx <- row_idx_all[sel]

    X_raw <- full_data[idx_chunk, , drop = FALSE]
    X_scaled <- scale(X_raw, center = scale_center, scale = scale_scale)

    # Assign to nearest micro-centroid, then map to meta-cluster
    micro_assign <- RANN::nn2(data = centroids, query = X_scaled, k = 1)$nn.idx[, 1]
    meta_assign <- micro_to_meta[micro_assign]
    meta_col <- meta_to_col[as.character(meta_assign)]

    # Accumulate marker sums and cell counts per patient x meta-cluster
    for (i in seq_along(sel)) {
      r <- row_idx[i]
      m <- meta_col[i]
      marker_sums[r, m, ] <- marker_sums[r, m, ] + X_raw[i, ]
      cell_counts[r, m] <- cell_counts[r, m] + 1L
    }

    if (verbose && (chunk %% max(1, floor(n_chunks / 5)) == 0)) {
      cat(sprintf("    [marker_features] Chunk %d/%d\n", chunk, n_chunks))
    }
  }

  # Compute marker means: sums / counts (NA where count == 0)
  marker_means <- array(NA_real_, dim = c(n_patients, n_meta, n_markers))
  for (j in seq_len(n_meta)) {
    for (m in seq_len(n_markers)) {
      nonzero <- cell_counts[, j] > 0
      marker_means[nonzero, j, m] <- marker_sums[nonzero, j, m] / cell_counts[nonzero, j]
    }
  }

  # Select top N markers per meta-cluster by cross-patient variance
  if (is.null(selected_markers)) {
    if (is.null(selection_patient_ids)) {
      selection_patient_ids <- patient_ids
    }
    sel_rows <- match(as.character(selection_patient_ids), patient_ids)
    sel_rows <- sel_rows[!is.na(sel_rows)]

    selected_markers <- list()
    for (j in seq_len(n_meta)) {
      variances <- numeric(n_markers)
      for (m in seq_len(n_markers)) {
        vals <- marker_means[sel_rows, j, m]
        vals <- vals[!is.na(vals)]
        variances[m] <- if (length(vals) >= 3) var(vals) else 0
      }
      names(variances) <- marker_names
      n_sel <- min(n_markers_per_cluster, sum(variances > 0))
      if (n_sel > 0) {
        top_idx <- order(variances, decreasing = TRUE)[seq_len(n_sel)]
        selected_markers[[as.character(meta_ids[j])]] <- marker_names[top_idx]
      } else {
        selected_markers[[as.character(meta_ids[j])]] <- character(0)
      }
    }
  }

  # Build feature matrix from selected markers
  feat_cols <- character(0)
  feat_list <- list()
  for (j in seq_len(n_meta)) {
    cluster_id <- as.character(meta_ids[j])
    markers_sel <- selected_markers[[cluster_id]]
    if (length(markers_sel) == 0) next
    for (mk in markers_sel) {
      mk_idx <- match(mk, marker_names)
      col_name <- paste0("marker_", mk, "_c", cluster_id)
      vals <- marker_means[, j, mk_idx]
      # Replace NA (no cells in that cluster for that patient) with 0
      vals[is.na(vals)] <- 0
      feat_list[[col_name]] <- vals
      feat_cols <- c(feat_cols, col_name)
    }
  }

  if (length(feat_list) == 0) {
    features <- matrix(0, nrow = n_patients, ncol = 0)
    rownames(features) <- patient_ids
  } else {
    features <- do.call(cbind, feat_list)
    rownames(features) <- patient_ids
    colnames(features) <- feat_cols
  }

  if (verbose) {
    cat(sprintf("    [marker_features] %d features (%d clusters x up to %d markers)\n",
                ncol(features), n_meta, n_markers_per_cluster))
  }

  list(
    features = features,
    selected_markers = selected_markers,
    n_features = ncol(features)
  )
}


#' Run nested CV with per-fold atlas building
#'
#' Proper nested CV where the atlas (micro-cluster centroids) is built
#' INSIDE each outer fold using only training patients. This avoids
#' the data leakage present in the global atlas approach.
#'
#' Architecture:
#'   Outer fold (performance estimation):
#'     1. Split patients into train/test
#'     2. Build atlas on TRAIN patients only
#'     3. Compute micro-counts for train and test patients using this atlas
#'     4. Inner CV on train patients to select Leiden resolution
#'     5. Final evaluation on test patients
#'
#' @param full_data Full cell-level matrix (cells x markers)
#' @param full_p_id Patient ID vector
#' @param patient_outcomes List with y (response), x (baseline), covariates (optional)
#' @param resolution_grid Vector of Leiden resolution values to try
#' @param n_outer_folds Number of outer CV folds
#' @param n_inner_folds Number of inner CV folds for resolution selection
#' @param n_micro Number of micro-clusters
#' @param n_neighbors Number of neighbors for centroid graph
#' @param kmeans_fit_n_cells Max cells for k-means fitting
#' @param assignment_chunk_size Chunk size for cell assignment
#' @param prediction_method "elasticnet" or "random_forest"
#' @param lambda_grid Lambda grid for elasticnet
#' @param alpha Elasticnet mixing parameter
#' @param seed_base Base random seed
#' @param n_cores Number of cores for parallel cv.glmnet
#' @param verbose Print progress
#' @return List with CV results
run_nested_cv_perfold_atlas <- function(full_data,
                                         full_p_id,
                                         patient_outcomes,
                                         resolution_grid = seq(0.5, 4.0, by = 0.5),
                                         n_outer_folds = 5,
                                         n_inner_folds = 3,
                                         n_micro = 5000,
                                         n_neighbors = 30,
                                         kmeans_fit_n_cells = 200000,
                                         assignment_chunk_size = 300000,
                                         assignment_max_cells_per_patient = NULL,
                                         train_micro_counts_source = c("full", "fit"),
                                         prediction_method = "elasticnet",
                                         lambda_grid = 10^seq(-3, 1, length.out = 20),
                                         alpha = 0.5,
                                         seed_base = 2027,
                                         group_ids = NULL,
                                         fold_strata = NULL,
                                         n_cores = 1,
                                         verbose = TRUE) {

  tic_total <- Sys.time()

  unique_p_ids <- names(patient_outcomes$y)
  n_patients <- length(unique_p_ids)

  if (verbose) {
    cat("\n============================================================\n")
    cat("PER-FOLD ATLAS NESTED CV\n")
    cat("============================================================\n")
    cat(sprintf("Patients: %d\n", n_patients))
    cat(sprintf("Outer folds: %d\n", n_outer_folds))
    cat(sprintf("Inner folds: %d\n", n_inner_folds))
    cat(sprintf("Micro-clusters: %d\n", n_micro))
    cat(sprintf("Resolution grid: %s\n", paste(resolution_grid, collapse = ", ")))
    cat(sprintf("Prediction method: %s\n", prediction_method))
    cat(sprintf("Cores: %d\n", n_cores))
    cat("\n")
  }

  train_micro_counts_source <- match.arg(train_micro_counts_source)

  # Create outer folds (optionally grouped, e.g. by base patient id)
  if (!is.null(group_ids)) {
    outer_folds <- create_group_folds(
      p_ids = unique_p_ids,
      group_ids = group_ids,
      k = n_outer_folds,
      seed = seed_base,
      strata = fold_strata
    )
  } else {
    set.seed(seed_base)
    outer_folds <- createFolds(seq_len(n_patients), k = n_outer_folds, list = TRUE)
  }

  # Storage for results
  results <- list(
    test_errors = numeric(n_outer_folds),
    train_errors = numeric(n_outer_folds),
    resolution_selected = numeric(n_outer_folds),
    k_found = integer(n_outer_folds),
    atlas_timing_s = numeric(n_outer_folds),
    test_predictions = list(),
    test_true = list(),
    inner_cv_validation = list(),
    inner_k_found = list()
  )

  # Optionally register parallel backend for cv.glmnet
  cl <- NULL
  if (n_cores > 1) {
    suppressPackageStartupMessages(library(doParallel))
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
  }

  # Process each outer fold
  for (fold_id in seq_len(n_outer_folds)) {

    if (verbose) {
      cat(sprintf("\n--- OUTER FOLD %d/%d ---\n", fold_id, n_outer_folds))
    }

    # Split patients
    test_p_ids <- unique_p_ids[outer_folds[[fold_id]]]
    train_p_ids <- setdiff(unique_p_ids, test_p_ids)

    if (verbose) {
      cat(sprintf("  Train patients: %d, Test patients: %d\n",
                  length(train_p_ids), length(test_p_ids)))
    }

    # =========================================================================
    # STEP 1: Build atlas on TRAINING patients only
    # =========================================================================
    atlas <- build_atlas_for_fold(
      full_data = full_data,
      full_p_id = full_p_id,
      train_patient_ids = train_p_ids,
      n_micro = n_micro,
      n_neighbors = n_neighbors,
      kmeans_fit_n_cells = kmeans_fit_n_cells,
      seed = seed_base + fold_id * 1000,
      verbose = verbose
    )
    results$atlas_timing_s[fold_id] <- atlas$timing_s

    # =========================================================================
    # STEP 2: Compute micro-counts for TRAIN patients
    # =========================================================================
    if (train_micro_counts_source == "fit") {
      if (verbose) cat("  Computing train patient micro-counts (from k-means fit sample)...\n")
      train_micro_counts <- compute_patient_micro_counts_from_assignments(
        p_ids = atlas$fit_p_id,
        micro_assignments = atlas$fit_micro_assignments,
        patient_ids = train_p_ids,
        n_micro = n_micro
      )
    } else {
      if (verbose) cat("  Computing train patient micro-counts (full assignment)...\n")
      train_micro_counts <- compute_patient_micro_counts(
        full_data = full_data,
        full_p_id = full_p_id,
        patient_ids = train_p_ids,
        centroids = atlas$centroids,
        scale_center = atlas$scale_center,
        scale_scale = atlas$scale_scale,
        chunk_size = assignment_chunk_size,
        max_cells_per_patient = assignment_max_cells_per_patient,
        seed = seed_base + fold_id * 1000 + 11,
        verbose = FALSE
      )
    }

    # =========================================================================
    # STEP 3: Inner CV to select Leiden resolution
    # =========================================================================
    if (verbose) cat("  Inner CV for resolution selection...\n")

    # Create inner folds on training patients
    set.seed(seed_base + fold_id * 100)
    if (!is.null(group_ids)) {
      inner_folds <- create_group_folds(
        p_ids = train_p_ids,
        group_ids = group_ids,
        k = n_inner_folds,
        seed = seed_base + fold_id * 100,
        strata = if (is.null(fold_strata)) NULL else fold_strata[train_p_ids]
      )
    } else {
      inner_folds <- createFolds(seq_len(length(train_p_ids)),
                                 k = n_inner_folds, list = TRUE)
    }

    # Pre-compute Leiden mappings for each resolution
    leiden_mappings <- lapply(resolution_grid, function(res) {
      meta <- two_stage_leiden_micro_to_meta(
        centroid_graph = atlas$centroid_graph,
        resolution = res
      )
      meta
    })

    k_found_vec <- sapply(leiden_mappings, function(m) m$k_found)

    # Evaluate each resolution
    resolution_scores <- numeric(length(resolution_grid))

    for (r_idx in seq_along(resolution_grid)) {
      micro_to_meta <- leiden_mappings[[r_idx]]$micro_to_meta

      # Collapse to meta-counts
      train_meta_counts <- two_stage_collapse_micro_to_meta_counts(
        micro_counts = train_micro_counts,
        micro_to_meta = micro_to_meta
      )

      # Build features
      X_train <- log(train_meta_counts + 1)
      X_train <- cbind(X_train, baseline = patient_outcomes$x[train_p_ids])
      if (!is.null(patient_outcomes$covariates)) {
        X_train <- cbind(X_train, patient_outcomes$covariates[train_p_ids, , drop = FALSE])
      }
      X_train <- as.matrix(X_train)
      y_train <- patient_outcomes$y[train_p_ids]

      # Inner CV error
      fold_errors <- numeric(n_inner_folds)
      for (ii in seq_along(inner_folds)) {
        val_idx <- inner_folds[[ii]]
        tr_idx <- setdiff(seq_len(nrow(X_train)), val_idx)

        set.seed(seed_base + fold_id * 100 + r_idx * 10 + ii)
        cvfit <- cv.glmnet(
          x = X_train[tr_idx, , drop = FALSE],
          y = y_train[tr_idx],
          alpha = alpha,
          nfolds = 5,
          standardize = TRUE,
          parallel = n_cores > 1
        )

        pred_val <- predict(cvfit, newx = X_train[val_idx, , drop = FALSE],
                            s = cvfit$lambda.min)
        fold_errors[ii] <- mean((pred_val - y_train[val_idx])^2)
      }

      resolution_scores[r_idx] <- mean(fold_errors)

      if (verbose) {
        cat(sprintf("    res=%.2f k=%d inner_mse=%.4f\n",
                    resolution_grid[r_idx], k_found_vec[r_idx], resolution_scores[r_idx]))
      }
    }

    # Select best resolution
    best_r_idx <- which.min(resolution_scores)
    best_resolution <- resolution_grid[best_r_idx]
    best_k_found <- k_found_vec[best_r_idx]
    best_micro_to_meta <- leiden_mappings[[best_r_idx]]$micro_to_meta

    results$resolution_selected[fold_id] <- best_resolution
    results$k_found[fold_id] <- best_k_found
    results$inner_cv_validation[[fold_id]] <- setNames(resolution_scores, resolution_grid)
    results$inner_k_found[[fold_id]] <- setNames(k_found_vec, resolution_grid)

    if (verbose) {
      cat(sprintf("  Selected: resolution=%.2f, k=%d\n", best_resolution, best_k_found))
    }

    # =========================================================================
    # STEP 4: Compute TEST patient micro-counts using SAME atlas
    # =========================================================================
    if (verbose) cat("  Computing test patient micro-counts...\n")
    test_micro_counts <- compute_patient_micro_counts(
      full_data = full_data,
      full_p_id = full_p_id,
      patient_ids = test_p_ids,
      centroids = atlas$centroids,
      scale_center = atlas$scale_center,
      scale_scale = atlas$scale_scale,
      chunk_size = assignment_chunk_size,
      max_cells_per_patient = assignment_max_cells_per_patient,
      seed = seed_base + fold_id * 1000 + 17,
      verbose = FALSE
    )

    # =========================================================================
    # STEP 5: Final model fitting and evaluation
    # =========================================================================
    # Collapse to meta-counts with best resolution
    train_meta_final <- two_stage_collapse_micro_to_meta_counts(
      train_micro_counts, best_micro_to_meta
    )
    test_meta_final <- two_stage_collapse_micro_to_meta_counts(
      test_micro_counts, best_micro_to_meta
    )

    # Build final features
    X_train_final <- log(train_meta_final + 1)
    X_train_final <- cbind(X_train_final, baseline = patient_outcomes$x[train_p_ids])
    if (!is.null(patient_outcomes$covariates)) {
      X_train_final <- cbind(X_train_final,
                             patient_outcomes$covariates[train_p_ids, , drop = FALSE])
    }

    X_test_final <- log(test_meta_final + 1)
    X_test_final <- cbind(X_test_final, baseline = patient_outcomes$x[test_p_ids])
    if (!is.null(patient_outcomes$covariates)) {
      X_test_final <- cbind(X_test_final,
                            patient_outcomes$covariates[test_p_ids, , drop = FALSE])
    }

    X_train_final <- as.matrix(X_train_final)
    X_test_final <- as.matrix(X_test_final)
    y_train <- patient_outcomes$y[train_p_ids]
    y_test <- patient_outcomes$y[test_p_ids]

    # Fit final model
    set.seed(seed_base + fold_id)
    if (prediction_method == "elasticnet") {
      cvfit <- cv.glmnet(
        x = X_train_final,
        y = y_train,
        alpha = alpha,
        nfolds = 5,
        standardize = TRUE,
        parallel = n_cores > 1
      )

      pred_train <- as.numeric(predict(cvfit, newx = X_train_final, s = cvfit$lambda.min))
      pred_test <- as.numeric(predict(cvfit, newx = X_test_final, s = cvfit$lambda.min))

    } else if (prediction_method == "random_forest") {
      library(randomForest)
      train_df <- as.data.frame(X_train_final)
      train_df$y <- y_train

      rf_model <- randomForest(y ~ ., data = train_df, ntree = 500)
      pred_train <- predict(rf_model, newdata = train_df)
      pred_test <- predict(rf_model, newdata = as.data.frame(X_test_final))
    }

    # Store results
    results$train_errors[fold_id] <- mean((pred_train - y_train)^2)
    results$test_errors[fold_id] <- mean((pred_test - y_test)^2)
    results$test_predictions[[fold_id]] <- pred_test
    results$test_true[[fold_id]] <- y_test

    if (verbose) {
      cat(sprintf("  Fold %d: train_MSE=%.4f, test_MSE=%.4f\n",
                  fold_id, results$train_errors[fold_id], results$test_errors[fold_id]))
    }
  }

  # Cleanup parallel backend
  if (!is.null(cl)) {
    parallel::stopCluster(cl)
  }

  # Compute summary statistics
  total_time <- as.numeric(difftime(Sys.time(), tic_total, units = "secs"))

  summary_results <- list(
    mean_test_error = mean(results$test_errors),
    se_test_error = sd(results$test_errors) / sqrt(n_outer_folds),
    mean_train_error = mean(results$train_errors),
    se_train_error = sd(results$train_errors) / sqrt(n_outer_folds),
    resolution_selected = results$resolution_selected,
    k_found = results$k_found,
    test_errors_per_fold = results$test_errors,
    train_errors_per_fold = results$train_errors,
    test_predictions = results$test_predictions,
    test_true = results$test_true,
    inner_cv_validation = results$inner_cv_validation,
    inner_k_found = results$inner_k_found,
    atlas_timing_s = results$atlas_timing_s,
    total_time_s = total_time,
    parameters = list(
      n_outer_folds = n_outer_folds,
      n_inner_folds = n_inner_folds,
      n_micro = n_micro,
      n_neighbors = n_neighbors,
      resolution_grid = resolution_grid,
      prediction_method = prediction_method,
      assignment_max_cells_per_patient = assignment_max_cells_per_patient
    )
  )

  if (verbose) {
    cat("\n============================================================\n")
    cat("SUMMARY\n")
    cat("============================================================\n")
    cat(sprintf("Mean test MSE:  %.4f (+/- %.4f)\n",
                summary_results$mean_test_error, summary_results$se_test_error))
    cat(sprintf("Mean train MSE: %.4f (+/- %.4f)\n",
                summary_results$mean_train_error, summary_results$se_train_error))
    cat(sprintf("Mean test RMSE: %.4f\n", sqrt(summary_results$mean_test_error)))
    cat(sprintf("Resolutions selected: %s\n",
                paste(summary_results$resolution_selected, collapse = ", ")))
    cat(sprintf("k_found: %s\n", paste(summary_results$k_found, collapse = ", ")))
    cat(sprintf("Total time: %.1f seconds\n", total_time))
    cat("============================================================\n")
  }

  summary_results
}

fit_predict_model <- function(method,
                              X_train,
                              y_train,
                              X_test,
                              seed,
                              alpha = 0.5,
                              n_cores = 1,
                              rf_ntree = 500,
                              foldid = NULL,
                              lambda_rule = c("lambda.min", "lambda.1se")) {
  method <- match.arg(method, c("elasticnet", "random_forest"))
  lambda_rule <- match.arg(lambda_rule)
  set.seed(seed)

  if (method == "elasticnet") {
    cv_args <- list(
      x = X_train,
      y = y_train,
      alpha = alpha,
      standardize = TRUE,
      parallel = n_cores > 1
    )
    if (!is.null(foldid)) {
      cv_args$foldid <- foldid
    } else {
      cv_args$nfolds <- 5L
    }
    cvfit <- do.call(cv.glmnet, cv_args)
    s_val <- if (lambda_rule == "lambda.min") cvfit$lambda.min else cvfit$lambda.1se
    pred_train <- as.numeric(predict(cvfit, newx = X_train, s = s_val))
    pred_test <- as.numeric(predict(cvfit, newx = X_test, s = s_val))
    return(list(pred_train = pred_train, pred_test = pred_test, model = cvfit))
  }

  use_ranger <- requireNamespace("ranger", quietly = TRUE)
  use_rf_pkg <- requireNamespace("randomForest", quietly = TRUE)
  if (!use_ranger && !use_rf_pkg) {
    stop("Need either 'ranger' or 'randomForest' installed for random_forest.")
  }

  train_df <- as.data.frame(X_train)
  train_df$y <- y_train
  test_df <- as.data.frame(X_test)

  if (use_ranger) {
    rf_model <- ranger::ranger(
      y ~ .,
      data = train_df,
      num.trees = rf_ntree,
      mtry = max(1L, floor(sqrt(ncol(X_train)))),
      min.node.size = 5,
      respect.unordered.factors = "order"
    )
    pred_train <- rf_model$predictions
    pred_test <- predict(rf_model, data = test_df)$predictions
    return(list(pred_train = pred_train, pred_test = pred_test, model = rf_model))
  }

  rf_model <- randomForest::randomForest(y ~ ., data = train_df, ntree = rf_ntree)
  pred_train <- predict(rf_model, newdata = train_df)
  pred_test <- predict(rf_model, newdata = test_df)
  list(pred_train = pred_train, pred_test = pred_test, model = rf_model)
}

#' Run per-fold atlas nested CV for multiple viruses (shared atlas per fold)
#'
#' This is the efficient multi-virus version: build atlas + micro-counts once per
#' outer fold, select a shared Leiden resolution using inner CV aggregated across
#' viruses, then evaluate per virus on the outer test fold.
#'
#' By default, folds are grouped by base patient ID (pid) to avoid leakage across
#' repeated years (e.g., *_yr1 and *_yr2 from the same pid landing in different folds).
#'
#' @param full_data Full cell-level matrix
#' @param full_p_id Patient ID per cell (p_id_yr)
#' @param outcome_obj Outcome data frame
#' @param meta_obj Metadata data frame
#' @param viruses Vector of virus names (e.g. c("h1n1","h3n2","bvic","byam"))
#' @param resolution_grid Leiden resolutions to try
#' @param n_outer_folds Outer CV folds
#' @param n_inner_folds Inner CV folds
#' @param n_micro Number of micro-clusters
#' @param n_neighbors Centroid kNN neighbors
#' @param kmeans_fit_n_cells Max training cells sampled for k-means
#' @param assignment_chunk_size Chunk size for test assignment
#' @param train_micro_counts_source "fit" uses only the k-means fit sample for train counts (fast)
#' @param prediction_method "elasticnet" or "random_forest"
#' @param model_mode "single_task" (one model per virus) or "multitask" (joint mgaussian EN)
#' @param alpha Elastic net alpha
#' @param n_cores Cores for cv.glmnet parallel folds
#' @param group_by_pid If TRUE, group folds by base pid
#' @param cohort_filter "all", "adult", or "peds"
#' @param add_cluster_ratios If TRUE, append log-ratio features from meta-cluster counts
#' @param max_ratio_features Max number of ratio features retained (training-variance ranked)
#' @param rf_ntree_inner RF trees for inner CV
#' @param rf_ntree_final RF trees for final fit
#' @param verbose Print progress
#' @return List with per-virus results and aggregated selection details
run_nested_cv_perfold_atlas_multivirus <- function(full_data,
                                                   full_p_id,
                                                   outcome_obj,
                                                   meta_obj,
                                                   viruses = c("h1n1", "h3n2", "bvic", "byam"),
                                                   resolution_grid = seq(0.5, 5.0, by = 0.25),
                                                   n_outer_folds = 5,
                                                   n_inner_folds = 3,
                                                   n_micro = 5000,
                                                   n_neighbors = 30,
                                                   kmeans_fit_n_cells = 1000000,
                                                   assignment_chunk_size = 300000,
                                                   assignment_max_cells_per_patient = NULL,
                                                   train_micro_counts_source = c("fit", "full"),
                                                   prediction_method = "elasticnet",
                                                   model_mode = c("single_task", "multitask"),
                                                   alpha = 0.5,
                                                   seed_base = 2027,
                                                   n_cores = 1,
                                                   group_by_pid = TRUE,
                                                   cohort_filter = c("all", "adult", "peds"),
                                                   add_cluster_ratios = FALSE,
                                                   max_ratio_features = 561,
                                                   add_marker_summaries = FALSE,
                                                   n_markers_per_cluster = 5L,
                                                   feature_transform = c("log_count", "proportion", "clr"),
                                                   lambda_rule = c("lambda.min", "lambda.1se"),
                                                   fold_seed = NULL,
                                                   permute_seed = NULL,
                                                   rf_ntree_inner = 300,
                                                   rf_ntree_final = 1000,
                                                   verbose = TRUE) {

  if (!exists("prepare_outcome_with_covariates")) {
    source(here::here("methods", "08-multi-outcome-covariates.R"))
  }

  train_micro_counts_source <- match.arg(train_micro_counts_source)
  model_mode <- match.arg(model_mode)
  cohort_filter <- match.arg(cohort_filter)
  feature_transform <- match.arg(feature_transform)
  lambda_rule <- match.arg(lambda_rule)

  # C1-SEED fix: use fold_seed for fold creation, seed_base for atlas/model seeds
  fold_seed_val <- if (!is.null(fold_seed)) fold_seed else seed_base

  if (model_mode == "multitask" && prediction_method != "elasticnet") {
    stop("model_mode='multitask' currently requires prediction_method='elasticnet'.")
  }

  # Prepare outcomes per virus
  outcomes <- list()
  for (v in viruses) {
    baseline_col <- paste0("geo_mean_", v, "_v1")
    target_col <- paste0("geo_mean_", v, "_v3")
    outcomes[[v]] <- prepare_outcome_with_covariates(
      outcome_obj, meta_obj,
      outcome_col_1 = baseline_col,
      outcome_col_2 = target_col,
      p_id = full_p_id,
      include_dose = TRUE,
      include_age = TRUE,
      include_sex = TRUE,
      include_time = TRUE
    )
  }

  common_p_ids <- Reduce(intersect, lapply(outcomes, function(o) names(o$y)))
  if (length(common_p_ids) < 50) {
    stop("Too few common observations across viruses: ", length(common_p_ids))
  }

  # S4-PERMUTATION: shuffle y across patients to break signal (null test)
  if (!is.null(permute_seed)) {
    set.seed(permute_seed)
    perm_order <- sample(length(common_p_ids))
    for (v in names(outcomes)) {
      y_vals <- outcomes[[v]]$y[common_p_ids]
      outcomes[[v]]$y[common_p_ids] <- y_vals[perm_order]
    }
    if (verbose) cat("  PERMUTATION: shuffled y with seed=", permute_seed, "\n")
  }

  cohort_all <- infer_cohort_from_p_id(common_p_ids)
  if (cohort_filter != "all") {
    common_p_ids <- common_p_ids[cohort_all[common_p_ids] == cohort_filter]
  }
  if (length(common_p_ids) < max(50, n_outer_folds * 8)) {
    stop("Too few observations after cohort_filter='", cohort_filter, "': ", length(common_p_ids))
  }

  # Shared covariates matrix (should be identical across viruses; take first)
  covariates <- outcomes[[viruses[1]]]$covariates[common_p_ids, , drop = FALSE]

  # Cohort stratification at observation level
  fold_strata <- infer_cohort_from_p_id(common_p_ids)
  if (length(unique(fold_strata)) <= 1) fold_strata <- NULL

  group_ids <- NULL
  if (group_by_pid) {
    group_ids <- setNames(base_pid_from_p_id_yr(common_p_ids), common_p_ids)
  }

  if (verbose) {
    cat("\n============================================================\n")
    cat("PER-FOLD ATLAS NESTED CV (MULTI-VIRUS)\n")
    cat("============================================================\n")
    cat("Viruses:            ", paste(viruses, collapse = ", "), "\n", sep = "")
    cat("Common observations:", length(common_p_ids), "\n", sep = "")
    cat("Cohort filter:      ", cohort_filter, "\n", sep = "")
    cat("Grouped by pid:     ", group_by_pid, "\n", sep = "")
    cat("Outer folds:        ", n_outer_folds, "\n", sep = "")
    cat("Inner folds:        ", n_inner_folds, "\n", sep = "")
    cat("Micro-clusters:     ", n_micro, "\n", sep = "")
    cat("K-means fit cells:  ", kmeans_fit_n_cells, "\n", sep = "")
    cat("Res grid:           ", paste(resolution_grid, collapse = ", "), "\n", sep = "")
    cat("Train counts:       ", train_micro_counts_source, "\n", sep = "")
    cat("Assign max cells/p: ",
        ifelse(is.null(assignment_max_cells_per_patient), "all", assignment_max_cells_per_patient),
        "\n", sep = "")
    cat("Method:             ", prediction_method, "\n", sep = "")
    cat("Model mode:         ", model_mode, "\n", sep = "")
    cat("Add ratios:         ", add_cluster_ratios, "\n", sep = "")
    if (isTRUE(add_cluster_ratios)) {
      cat("Max ratio features: ", max_ratio_features, "\n", sep = "")
    }
    cat("Feature transform:  ", feature_transform, "\n", sep = "")
    cat("Lambda rule:        ", lambda_rule, "\n", sep = "")
    cat("Fold seed:          ", ifelse(is.null(fold_seed), "same as seed_base", fold_seed), "\n", sep = "")
    cat("Baseline model:     OLS (deterministic)\n", sep = "")
    cat("Lambda CV grouping: patient-grouped foldid\n", sep = "")
    cat("Cores (glmnet):     ", n_cores, "\n\n", sep = "")
  }

  # doParallel for glmnet
  cl <- NULL
  if (prediction_method == "elasticnet" && n_cores > 1) {
    suppressPackageStartupMessages(library(doParallel))
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    on.exit({
      try(parallel::stopCluster(cl), silent = TRUE)
    }, add = TRUE)
  }

  # Outer folds (C1-SEED fix: use fold_seed_val so all experiments share folds)
  if (!is.null(group_ids)) {
    outer_folds <- create_group_folds(
      p_ids = common_p_ids,
      group_ids = group_ids,
      k = n_outer_folds,
      seed = fold_seed_val,
      strata = fold_strata
    )
  } else {
    set.seed(fold_seed_val)
    if (is.null(fold_strata)) {
      outer_folds <- caret::createFolds(seq_along(common_p_ids), k = n_outer_folds, list = TRUE)
    } else {
      outer_folds <- caret::createFolds(fold_strata, k = n_outer_folds, list = TRUE)
    }
  }

  per_fold <- vector("list", n_outer_folds)

  for (fold_id in seq_len(n_outer_folds)) {
    if (verbose) cat(sprintf("\n--- OUTER FOLD %d/%d ---\n", fold_id, n_outer_folds))

    test_idx <- outer_folds[[fold_id]]
    test_p_ids <- common_p_ids[test_idx]
    train_p_ids <- setdiff(common_p_ids, test_p_ids)

    if (verbose) {
      cat(sprintf("  Train obs: %d, Test obs: %d\n", length(train_p_ids), length(test_p_ids)))
      if (!is.null(group_ids)) {
        cat(sprintf("  Train pid: %d, Test pid: %d\n",
                    length(unique(group_ids[train_p_ids])),
                    length(unique(group_ids[test_p_ids]))))
      }
    }

    atlas <- build_atlas_for_fold(
      full_data = full_data,
      full_p_id = full_p_id,
      train_patient_ids = train_p_ids,
      n_micro = n_micro,
      n_neighbors = n_neighbors,
      kmeans_fit_n_cells = kmeans_fit_n_cells,
      seed = seed_base + fold_id * 1000,
      verbose = verbose
    )

    if (train_micro_counts_source == "fit") {
      train_micro_counts <- compute_patient_micro_counts_from_assignments(
        p_ids = atlas$fit_p_id,
        micro_assignments = atlas$fit_micro_assignments,
        patient_ids = train_p_ids,
        n_micro = n_micro
      )
    } else {
      train_micro_counts <- compute_patient_micro_counts(
        full_data = full_data,
        full_p_id = full_p_id,
        patient_ids = train_p_ids,
        centroids = atlas$centroids,
        scale_center = atlas$scale_center,
        scale_scale = atlas$scale_scale,
        chunk_size = assignment_chunk_size,
        max_cells_per_patient = assignment_max_cells_per_patient,
        seed = seed_base + fold_id * 1000 + 11,
        verbose = FALSE
      )
    }

    if (verbose) cat("  Computing test micro-counts (full assignment)...\n")
    test_micro_counts <- compute_patient_micro_counts(
      full_data = full_data,
      full_p_id = full_p_id,
      patient_ids = test_p_ids,
      centroids = atlas$centroids,
      scale_center = atlas$scale_center,
      scale_scale = atlas$scale_scale,
      chunk_size = assignment_chunk_size,
      max_cells_per_patient = assignment_max_cells_per_patient,
      seed = seed_base + fold_id * 1000 + 17,
      verbose = FALSE
    )

    # Precompute Leiden mappings for each resolution (outer fold)
    leiden_mappings <- lapply(resolution_grid, function(res) {
      two_stage_leiden_micro_to_meta(atlas$centroid_graph, resolution = res)
    })
    k_found_vec <- sapply(leiden_mappings, function(m) m$k_found)

    # Inner folds on training observations
    if (!is.null(group_ids)) {
      inner_folds <- create_group_folds(
        p_ids = train_p_ids,
        group_ids = group_ids,
        k = n_inner_folds,
        seed = fold_seed_val + fold_id * 100,
        strata = if (is.null(fold_strata)) NULL else fold_strata[train_p_ids]
      )
    } else {
      set.seed(fold_seed_val + fold_id * 100)
      if (is.null(fold_strata)) {
        inner_folds <- caret::createFolds(seq_along(train_p_ids), k = n_inner_folds, list = TRUE)
      } else {
        inner_folds <- caret::createFolds(fold_strata[train_p_ids], k = n_inner_folds, list = TRUE)
      }
    }

    # C2-BASELINE fix: Pre-compute OLS baseline MSE for each inner fold
    # (resolution-independent, computed ONCE outside the r_idx loop)
    inner_base_mse <- matrix(NA_real_, nrow = length(inner_folds), ncol = length(viruses),
                             dimnames = list(NULL, viruses))
    for (ii in seq_along(inner_folds)) {
      val_idx <- inner_folds[[ii]]
      val_p <- train_p_ids[val_idx]
      tr_p <- setdiff(train_p_ids, val_p)

      if (model_mode == "multitask") {
        base_train_mat <- sapply(viruses, function(v) outcomes[[v]]$x[tr_p])
        base_val_mat <- sapply(viruses, function(v) outcomes[[v]]$x[val_p])
        colnames(base_train_mat) <- paste0("baseline_", viruses)
        colnames(base_val_mat) <- paste0("baseline_", viruses)
        cov_tr <- covariates[tr_p, , drop = FALSE]
        cov_val <- covariates[val_p, , drop = FALSE]
        X_base_tr <- as.matrix(cbind(base_train_mat, cov_tr))
        X_base_val <- as.matrix(cbind(base_val_mat, cov_val))
        Y_tr <- sapply(viruses, function(v) outcomes[[v]]$y[tr_p])
        Y_val <- sapply(viruses, function(v) outcomes[[v]]$y[val_p])
        colnames(Y_tr) <- viruses
        colnames(Y_val) <- viruses
        base_fit <- fit_ols_baseline_multi(X_base_tr, Y_tr, X_base_val)
        for (v in viruses) {
          inner_base_mse[ii, v] <- mean((base_fit$pred_test[, v] - Y_val[, v])^2)
        }
      } else {
        for (vv in seq_along(viruses)) {
          v <- viruses[vv]
          y_tr <- outcomes[[v]]$y[tr_p]
          y_val <- outcomes[[v]]$y[val_p]
          base_tr <- outcomes[[v]]$x[tr_p]
          base_val <- outcomes[[v]]$x[val_p]
          cov_tr <- covariates[tr_p, , drop = FALSE]
          cov_val <- covariates[val_p, , drop = FALSE]
          X_base_tr <- as.matrix(cbind(baseline = base_tr, cov_tr))
          X_base_val <- as.matrix(cbind(baseline = base_val, cov_val))
          base_fit <- fit_ols_baseline(X_base_tr, y_tr, X_base_val)
          inner_base_mse[ii, v] <- mean((base_fit$pred_test - y_val)^2)
        }
      }
    }

    # Score each resolution using median ratio across viruses
    resolution_scores <- numeric(length(resolution_grid))
    for (r_idx in seq_along(resolution_grid)) {
      micro_to_meta <- leiden_mappings[[r_idx]]$micro_to_meta
      train_meta_counts <- two_stage_collapse_micro_to_meta_counts(
        micro_counts = train_micro_counts,
        micro_to_meta = micro_to_meta
      )

      # C4-LOGCOUNTS fix: use transform_counts instead of log(count+1)
      cluster_train_all <- transform_counts(train_meta_counts, method = feature_transform)
      if (isTRUE(add_cluster_ratios)) {
        ratio_feats <- build_ratio_feature_matrices(
          train_meta_counts = train_meta_counts,
          max_ratio_features = max_ratio_features,
          selector = "variance"
        )
        cluster_train_all <- cbind(cluster_train_all, ratio_feats$train)
      }

      fold_scores <- numeric(length(inner_folds))
      for (ii in seq_along(inner_folds)) {
        val_idx <- inner_folds[[ii]]
        val_p <- train_p_ids[val_idx]
        tr_p <- setdiff(train_p_ids, val_p)

        if (model_mode == "multitask") {
          base_train_mat <- sapply(viruses, function(v) outcomes[[v]]$x[tr_p])
          base_val_mat <- sapply(viruses, function(v) outcomes[[v]]$x[val_p])
          colnames(base_train_mat) <- paste0("baseline_", viruses)
          colnames(base_val_mat) <- paste0("baseline_", viruses)
          cov_tr <- covariates[tr_p, , drop = FALSE]
          cov_val <- covariates[val_p, , drop = FALSE]
          X_base_tr <- as.matrix(cbind(base_train_mat, cov_tr))
          X_base_val <- as.matrix(cbind(base_val_mat, cov_val))
          X_clust_tr <- as.matrix(cbind(cluster_train_all[tr_p, , drop = FALSE], X_base_tr))
          X_clust_val <- as.matrix(cbind(cluster_train_all[val_p, , drop = FALSE], X_base_val))
          Y_tr <- sapply(viruses, function(v) outcomes[[v]]$y[tr_p])
          Y_val <- sapply(viruses, function(v) outcomes[[v]]$y[val_p])
          colnames(Y_tr) <- viruses
          colnames(Y_val) <- viruses

          # M3 fix: grouped foldid for lambda CV
          inner_foldid <- compute_grouped_foldid(tr_p, n_folds = 5L,
                                                 seed = seed_base + fold_id * 10000 + ii)
          seed_fit <- seed_base + fold_id * 10000 + r_idx * 100 + ii * 10
          clust_fit <- fit_predict_multitask_en(
            X_train = X_clust_tr, Y_train = Y_tr, X_test = X_clust_val,
            seed = seed_fit + 1, alpha = alpha, n_cores = n_cores,
            foldid = inner_foldid, lambda_rule = lambda_rule
          )

          ratios <- sapply(viruses, function(v) {
            mse_clust <- mean((clust_fit$pred_test[, v] - Y_val[, v])^2)
            mse_clust / inner_base_mse[ii, v]
          })
          fold_scores[ii] <- median(ratios)
        } else {
          ratios <- numeric(length(viruses))
          for (vv in seq_along(viruses)) {
            v <- viruses[vv]
            y_tr <- outcomes[[v]]$y[tr_p]
            y_val <- outcomes[[v]]$y[val_p]
            base_tr <- outcomes[[v]]$x[tr_p]
            base_val <- outcomes[[v]]$x[val_p]
            cov_tr <- covariates[tr_p, , drop = FALSE]
            cov_val <- covariates[val_p, , drop = FALSE]
            X_clust_tr <- as.matrix(cbind(cluster_train_all[tr_p, , drop = FALSE],
                                          baseline = base_tr, cov_tr))
            X_clust_val <- as.matrix(cbind(cluster_train_all[val_p, , drop = FALSE],
                                           baseline = base_val, cov_val))

            inner_foldid <- compute_grouped_foldid(tr_p, n_folds = 5L,
                                                   seed = seed_base + fold_id * 10000 + ii + vv)
            seed_fit <- seed_base + fold_id * 10000 + r_idx * 100 + ii * 10 + vv
            clust_fit <- fit_predict_model(
              method = prediction_method,
              X_train = X_clust_tr, y_train = y_tr, X_test = X_clust_val,
              seed = seed_fit + 1, alpha = alpha, n_cores = n_cores,
              rf_ntree = rf_ntree_inner,
              foldid = inner_foldid, lambda_rule = lambda_rule
            )

            mse_clust <- mean((clust_fit$pred_test - y_val)^2)
            ratios[vv] <- mse_clust / inner_base_mse[ii, v]
          }
          fold_scores[ii] <- median(ratios)
        }
      }

      resolution_scores[r_idx] <- mean(fold_scores)
      if (verbose) {
        cat(sprintf("    res=%.2f k=%d inner_score=%.4f\n",
                    resolution_grid[r_idx], k_found_vec[r_idx], resolution_scores[r_idx]))
      }
    }

    best_r_idx <- which.min(resolution_scores)
    best_res <- resolution_grid[best_r_idx]
    best_k_found <- k_found_vec[best_r_idx]
    best_micro_to_meta <- leiden_mappings[[best_r_idx]]$micro_to_meta

    if (verbose) cat(sprintf("  Selected: res=%.2f (k=%d)\n", best_res, best_k_found))

    # Final evaluation (baseline vs clusters) per virus on outer test fold
    train_meta_final <- two_stage_collapse_micro_to_meta_counts(train_micro_counts, best_micro_to_meta)
    test_meta_final <- two_stage_collapse_micro_to_meta_counts(test_micro_counts, best_micro_to_meta)

    # C4-LOGCOUNTS fix: use transform_counts instead of log(count+1)
    cluster_train_final <- transform_counts(train_meta_final, method = feature_transform)
    cluster_test_final <- transform_counts(test_meta_final, method = feature_transform)
    n_ratio_features_used <- 0L
    if (isTRUE(add_cluster_ratios)) {
      ratio_feats <- build_ratio_feature_matrices(
        train_meta_counts = train_meta_final,
        test_meta_counts = test_meta_final,
        max_ratio_features = max_ratio_features,
        selector = "variance"
      )
      cluster_train_final <- cbind(cluster_train_final, ratio_feats$train)
      cluster_test_final <- cbind(cluster_test_final, ratio_feats$test)
      n_ratio_features_used <- ncol(ratio_feats$train)
    }
    n_marker_features_used <- 0L
    if (isTRUE(add_marker_summaries)) {
      if (verbose) cat("  Computing targeted marker summary features (train)...\n")
      train_marker_res <- compute_targeted_marker_features(
        full_data = full_data,
        full_p_id = full_p_id,
        patient_ids = train_p_ids,
        centroids = atlas$centroids,
        scale_center = atlas$scale_center,
        scale_scale = atlas$scale_scale,
        micro_to_meta = best_micro_to_meta,
        n_markers_per_cluster = n_markers_per_cluster,
        selection_patient_ids = train_p_ids,
        chunk_size = assignment_chunk_size,
        verbose = verbose
      )
      if (verbose) cat("  Computing targeted marker summary features (test)...\n")
      test_marker_res <- compute_targeted_marker_features(
        full_data = full_data,
        full_p_id = full_p_id,
        patient_ids = test_p_ids,
        centroids = atlas$centroids,
        scale_center = atlas$scale_center,
        scale_scale = atlas$scale_scale,
        micro_to_meta = best_micro_to_meta,
        n_markers_per_cluster = n_markers_per_cluster,
        selection_patient_ids = NULL,
        selected_markers = train_marker_res$selected_markers,
        chunk_size = assignment_chunk_size,
        verbose = verbose
      )
      cluster_train_final <- cbind(cluster_train_final, train_marker_res$features)
      cluster_test_final <- cbind(cluster_test_final, test_marker_res$features)
      n_marker_features_used <- train_marker_res$n_features
    }

    # M3 fix: grouped foldid for outer lambda CV
    outer_foldid <- compute_grouped_foldid(train_p_ids, n_folds = 5L,
                                           seed = seed_base + fold_id * 1000 + 99)

    virus_results <- list()
    if (model_mode == "multitask") {
      cov_train <- covariates[train_p_ids, , drop = FALSE]
      cov_test <- covariates[test_p_ids, , drop = FALSE]

      base_train_mat <- sapply(viruses, function(v) outcomes[[v]]$x[train_p_ids])
      base_test_mat <- sapply(viruses, function(v) outcomes[[v]]$x[test_p_ids])
      colnames(base_train_mat) <- paste0("baseline_", viruses)
      colnames(base_test_mat) <- paste0("baseline_", viruses)

      X_base_train <- as.matrix(cbind(base_train_mat, cov_train))
      X_base_test <- as.matrix(cbind(base_test_mat, cov_test))
      X_clust_train <- as.matrix(cbind(cluster_train_final[train_p_ids, , drop = FALSE], X_base_train))
      X_clust_test <- as.matrix(cbind(cluster_test_final[test_p_ids, , drop = FALSE], X_base_test))

      Y_train <- sapply(viruses, function(v) outcomes[[v]]$y[train_p_ids])
      Y_test <- sapply(viruses, function(v) outcomes[[v]]$y[test_p_ids])
      colnames(Y_train) <- viruses
      colnames(Y_test) <- viruses

      # C2 fix: OLS baseline (deterministic, no seed needed)
      base_fit <- fit_ols_baseline_multi(X_base_train, Y_train, X_base_test)

      seed_fit <- seed_base + fold_id * 1000 + best_r_idx * 10
      clust_fit <- fit_predict_multitask_en(
        X_train = X_clust_train, Y_train = Y_train, X_test = X_clust_test,
        seed = seed_fit + 1, alpha = alpha, n_cores = n_cores,
        foldid = outer_foldid, lambda_rule = lambda_rule
      )

      for (v in viruses) {
        y_train <- Y_train[, v]
        y_test <- Y_test[, v]
        virus_results[[v]] <- list(
          baseline_mse = mean((base_fit$pred_test[, v] - y_test)^2),
          clusters_mse = mean((clust_fit$pred_test[, v] - y_test)^2),
          baseline_train_mse = mean((base_fit$pred_train[, v] - y_train)^2),
          clusters_train_mse = mean((clust_fit$pred_train[, v] - y_train)^2),
          n_train = length(y_train),
          n_test = length(y_test)
        )
      }
    } else {
      for (v in viruses) {
        y_train <- outcomes[[v]]$y[train_p_ids]
        y_test <- outcomes[[v]]$y[test_p_ids]
        base_train <- outcomes[[v]]$x[train_p_ids]
        base_test <- outcomes[[v]]$x[test_p_ids]

        cov_train <- covariates[train_p_ids, , drop = FALSE]
        cov_test <- covariates[test_p_ids, , drop = FALSE]

        X_base_train <- as.matrix(cbind(baseline = base_train, cov_train))
        X_base_test <- as.matrix(cbind(baseline = base_test, cov_test))
        X_clust_train <- as.matrix(cbind(cluster_train_final[train_p_ids, , drop = FALSE],
                                         baseline = base_train, cov_train))
        X_clust_test <- as.matrix(cbind(cluster_test_final[test_p_ids, , drop = FALSE],
                                        baseline = base_test, cov_test))

        # C2 fix: OLS baseline (deterministic)
        base_fit <- fit_ols_baseline(X_base_train, y_train, X_base_test)

        seed_fit <- seed_base + fold_id * 1000 + best_r_idx * 10
        clust_fit <- fit_predict_model(
          method = prediction_method,
          X_train = X_clust_train, y_train = y_train, X_test = X_clust_test,
          seed = seed_fit + 1, alpha = alpha, n_cores = n_cores,
          rf_ntree = rf_ntree_final,
          foldid = outer_foldid, lambda_rule = lambda_rule
        )

        virus_results[[v]] <- list(
          baseline_mse = mean((base_fit$pred_test - y_test)^2),
          clusters_mse = mean((clust_fit$pred_test - y_test)^2),
          baseline_train_mse = mean((base_fit$pred_train - y_train)^2),
          clusters_train_mse = mean((clust_fit$pred_train - y_train)^2),
          n_train = length(y_train),
          n_test = length(y_test)
        )
      }
    }

    per_fold[[fold_id]] <- list(
      fold_id = fold_id,
      train_obs = length(train_p_ids),
      test_obs = length(test_p_ids),
      selected_resolution = best_res,
      k_found = best_k_found,
      n_ratio_features = n_ratio_features_used,
      n_marker_features = n_marker_features_used,
      resolution_scores = setNames(resolution_scores, resolution_grid),
      k_found_by_res = setNames(k_found_vec, resolution_grid),
      viruses = virus_results
    )
  }

  # Aggregate results across folds
  agg <- data.frame(
    virus = viruses,
    baseline_test_mse = NA_real_,
    clusters_test_mse = NA_real_,
    baseline_test_rmse = NA_real_,
    clusters_test_rmse = NA_real_,
    stringsAsFactors = FALSE
  )

  # M2-FOLD-WEIGHT fix: weight MSE by fold test size (pool squared errors)
  for (v in viruses) {
    base_mse <- sapply(per_fold, function(f) f$viruses[[v]]$baseline_mse)
    clust_mse <- sapply(per_fold, function(f) f$viruses[[v]]$clusters_mse)
    n_test <- sapply(per_fold, function(f) f$viruses[[v]]$n_test)
    w <- n_test / sum(n_test)
    agg[agg$virus == v, "baseline_test_mse"] <- sum(w * base_mse)
    agg[agg$virus == v, "clusters_test_mse"] <- sum(w * clust_mse)
    agg[agg$virus == v, "baseline_test_rmse"] <- sqrt(sum(w * base_mse))
    agg[agg$virus == v, "clusters_test_rmse"] <- sqrt(sum(w * clust_mse))
  }

  selected_k <- sapply(per_fold, function(f) f$k_found)
  selected_res <- sapply(per_fold, function(f) f$selected_resolution)

  list(
    params = list(
      viruses = viruses,
      resolution_grid = resolution_grid,
      n_outer_folds = n_outer_folds,
      n_inner_folds = n_inner_folds,
      n_micro = n_micro,
      n_neighbors = n_neighbors,
      kmeans_fit_n_cells = kmeans_fit_n_cells,
      assignment_chunk_size = assignment_chunk_size,
      assignment_max_cells_per_patient = assignment_max_cells_per_patient,
      train_micro_counts_source = train_micro_counts_source,
      prediction_method = prediction_method,
      model_mode = model_mode,
      alpha = alpha,
      seed_base = seed_base,
      n_cores = n_cores,
      group_by_pid = group_by_pid,
      cohort_filter = cohort_filter,
      add_cluster_ratios = add_cluster_ratios,
      max_ratio_features = max_ratio_features,
      feature_transform = feature_transform,
      lambda_rule = lambda_rule,
      fold_seed = fold_seed,
      rf_ntree_inner = rf_ntree_inner,
      rf_ntree_final = rf_ntree_final
    ),
    common_observations = common_p_ids,
    n_groups = if (is.null(group_ids)) NA_integer_ else length(unique(group_ids[common_p_ids])),
    selected_resolutions = selected_res,
    selected_k_found = selected_k,
    per_fold = per_fold,
    aggregated = agg
  )
}


#' Run per-fold atlas CV for multiple viruses
#'
#' Runs the per-fold atlas nested CV for each virus independently and
#' aggregates results.
#'
#' @param full_data Full cell-level matrix
#' @param full_p_id Patient ID vector
#' @param outcome_obj Outcome data frame
#' @param meta_obj Metadata data frame
#' @param viruses Vector of virus names (e.g., c("h1n1", "h3n2"))
#' @param ... Additional arguments passed to run_nested_cv_perfold_atlas
#' @return List with results per virus and aggregated summary
run_perfold_atlas_multivirus <- function(full_data,
                                          full_p_id,
                                          outcome_obj,
                                          meta_obj,
                                          viruses = c("h1n1", "h3n2", "bvic", "byam"),
                                          n_outer_folds = 5,
                                          n_inner_folds = 3,
                                          resolution_grid = seq(0.5, 4.0, by = 0.5),
                                          n_micro = 5000,
                                          n_neighbors = 30,
                                          kmeans_fit_n_cells = 200000,
                                          prediction_method = "elasticnet",
                                          n_cores = 1,
                                          verbose = TRUE) {

  # Source outcome preparation if needed
  if (!exists("prepare_outcome_with_covariates")) {
    source(here::here("methods", "08-multi-outcome-covariates.R"))
  }

  results_by_virus <- list()

  for (v in viruses) {
    if (verbose) {
      cat("\n################################################################\n")
      cat(sprintf("VIRUS: %s\n", toupper(v)))
      cat("################################################################\n")
    }

    # Prepare outcomes for this virus
    baseline_col <- paste0("geo_mean_", v, "_v1")
    target_col <- paste0("geo_mean_", v, "_v3")

    patient_outcomes <- prepare_outcome_with_covariates(
      outcome_obj, meta_obj,
      outcome_col_1 = baseline_col,
      outcome_col_2 = target_col,
      p_id = full_p_id,
      include_dose = TRUE,
      include_age = TRUE,
      include_sex = TRUE,
      include_time = TRUE
    )

    # Run per-fold atlas CV
    results_by_virus[[v]] <- run_nested_cv_perfold_atlas(
      full_data = full_data,
      full_p_id = full_p_id,
      patient_outcomes = patient_outcomes,
      resolution_grid = resolution_grid,
      n_outer_folds = n_outer_folds,
      n_inner_folds = n_inner_folds,
      n_micro = n_micro,
      n_neighbors = n_neighbors,
      kmeans_fit_n_cells = kmeans_fit_n_cells,
      prediction_method = prediction_method,
      n_cores = n_cores,
      verbose = verbose
    )
  }

  # Aggregate across viruses
  aggregated <- data.frame(
    virus = viruses,
    mean_test_mse = sapply(results_by_virus, function(r) r$mean_test_error),
    se_test_mse = sapply(results_by_virus, function(r) r$se_test_error),
    mean_train_mse = sapply(results_by_virus, function(r) r$mean_train_error),
    test_rmse = sapply(results_by_virus, function(r) sqrt(r$mean_test_error)),
    train_rmse = sapply(results_by_virus, function(r) sqrt(r$mean_train_error)),
    median_k_found = sapply(results_by_virus, function(r) median(r$k_found)),
    total_time_s = sapply(results_by_virus, function(r) r$total_time_s)
  )

  list(
    by_virus = results_by_virus,
    aggregated = aggregated
  )
}


# ==============================================================================
# Quick test function for development
# ==============================================================================

#' Quick test of per-fold atlas with reduced parameters
#'
#' For testing the implementation with:
#'   - 2 outer folds (not 5)
#'   - 50-100 patients (not all)
#'   - 1000 micro-clusters (not 5000)
#'   - Reduced resolution grid
#'
#' @param full_data Full cell-level matrix
#' @param full_p_id Patient ID vector
#' @param patient_outcomes Prepared patient outcomes
#' @param n_patients_test Number of patients to use (NULL = all)
#' @param n_micro Number of micro-clusters (default: 1000)
#' @param n_outer_folds Number of outer folds (default: 2)
#' @param verbose Print progress
#' @return CV results
quick_test_perfold_atlas <- function(full_data,
                                      full_p_id,
                                      patient_outcomes,
                                      n_patients_test = 100,
                                      n_micro = 1000,
                                      n_outer_folds = 2,
                                      verbose = TRUE) {

  # Optionally subset patients
  all_patient_ids <- names(patient_outcomes$y)

  if (!is.null(n_patients_test) && n_patients_test < length(all_patient_ids)) {
    set.seed(12345)
    use_patient_ids <- sample(all_patient_ids, n_patients_test)

    # Subset patient_outcomes
    patient_outcomes$y <- patient_outcomes$y[use_patient_ids]
    patient_outcomes$x <- patient_outcomes$x[use_patient_ids]
    if (!is.null(patient_outcomes$covariates)) {
      patient_outcomes$covariates <- patient_outcomes$covariates[use_patient_ids, , drop = FALSE]
    }

    if (verbose) {
      cat(sprintf("Testing with %d patients (subset from %d)\n",
                  n_patients_test, length(all_patient_ids)))
    }
  }

  # Run with reduced parameters
  run_nested_cv_perfold_atlas(
    full_data = full_data,
    full_p_id = full_p_id,
    patient_outcomes = patient_outcomes,
    resolution_grid = c(0.5, 1.0, 2.0),  # Reduced grid
    n_outer_folds = n_outer_folds,
    n_inner_folds = 2,
    n_micro = n_micro,
    n_neighbors = 30,
    kmeans_fit_n_cells = 100000,
    assignment_chunk_size = 200000,
    prediction_method = "elasticnet",
    n_cores = 1,
    verbose = verbose
  )
}
