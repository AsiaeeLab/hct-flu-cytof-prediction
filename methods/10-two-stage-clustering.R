# ==============================================================================
# Two-Stage Clustering for Large-Scale Nested CV
# ==============================================================================
# This module implements a scalable two-stage clustering scheme:
#   1) Micro-clustering of cells (e.g., k-means or FlowSOM) with many clusters
#   2) Hierarchical clustering on micro-cluster centroids
#      → Meta-clusters obtained by cutting the dendrogram at different k
#
# The goal is to:
#   - Avoid O(n^2) hierarchical clustering on millions of cells
#   - Still have a nested, interpretable k path via the dendrogram
#   - Plug into a proper nested CV pipeline without touching legacy code
#
# This file defines:
#   - two_stage_microclustering()
#   - tune_k_on_microstructure()
#   - proper_nested_cv_two_stage()
#   - proper_nested_cv_two_stage_with_subsampling()
#
# It reuses:
#   - Outcome preparation / feature builders from methods/02-data-preparation.R
#   - Covariate helpers from methods/08-multi-outcome-covariates.R
#   - Baseline + aggregation helpers from methods/09-proper-nested-cv.R
# ==============================================================================


#' Two-stage microclustering + hierarchical meta-clustering
#'
#' Stage 1: micro-clustering using a scalable method (k-means or FlowSOM). \cr
#' Stage 2: hierarchical clustering of micro-cluster centroids (Ward.D2) to
#'          define a nested tree that can be cut at arbitrary meta k.
#'
#' @param data Numeric matrix of scaled cell-level features (cells × markers)
#' @param micro_method Micro-clustering method: "kmeans" or "flowsom"
#' @param n_microclusters Number of micro-clusters to form
#' @param kmeans_nstart Number of random starts for k-means
#' @param kmeans_iter_max Maximum iterations for k-means
#' @param verbose Print progress messages
#' @return List with micro_clusters (per cell), centroids (n_micro × p),
#'         hc (hclust object), and metadata
two_stage_microclustering <- function(data,
                                      micro_method = "kmeans",
                                      n_microclusters = 500,
                                      kmeans_nstart = 10,
                                      kmeans_iter_max = 100,
                                      build_hc = TRUE,
                                      hc_method = "ward.D2",
                                      verbose = FALSE) {

  if (verbose) {
    cat("      [two-stage] Microclustering with",
        micro_method, "(", n_microclusters, "clusters )...\n")
  }

  n_cells <- nrow(data)
  if (n_microclusters >= n_cells) {
    stop("n_microclusters (", n_microclusters,
         ") must be < number of cells (", n_cells, ").")
  }

  if (micro_method == "kmeans") {
    set.seed(123)
    km <- kmeans(data,
                 centers = n_microclusters,
                 nstart = kmeans_nstart,
                 iter.max = kmeans_iter_max)
    micro_clusters <- km$cluster
    centroids <- km$centers

  } else if (micro_method == "flowsom") {
    if (!exists("flowsom_clustering")) {
      stop("flowsom_clustering() not found. ",
           "Source methods/03-clustering-methods.R before calling.")
    }

    # Reuse existing FlowSOM wrapper; treat its k clusters as micro-clusters
    fs <- flowsom_clustering(data, k = n_microclusters)
    micro_clusters <- fs$clusters

    # Compute centroids per micro-cluster
    centroids <- sapply(1:n_microclusters, function(cl_id) {
      colMeans(data[micro_clusters == cl_id, , drop = FALSE])
    })
    centroids <- t(centroids) # k × p

  } else {
    stop("Unsupported micro_method: ", micro_method,
         ". Use 'kmeans' or 'flowsom'.")
  }

  # Optional hierarchical clustering on micro-centroids (can be very expensive
  # for large n_microclusters)
  hc <- NULL
  if (build_hc) {
    if (verbose) {
      cat("      [two-stage] Building hierarchical tree on centroids...\n")
    }
    hc <- hclust(dist(centroids), method = hc_method)
  }

  list(
    micro_clusters = micro_clusters,
    centroids = centroids,
    hc = hc,
    micro_method = micro_method,
    n_microclusters = n_microclusters
  )
}


#' Build meta-cluster labels for all training cells
#'
#' @param micro_clusters Integer vector of micro-cluster IDs per cell
#' @param hc hclust object built on micro-centroids
#' @param k_meta Number of meta-clusters
#' @return Integer vector of meta-cluster IDs (1..k_meta) per cell
two_stage_meta_clusters <- function(micro_clusters, hc, k_meta) {
  micro_to_meta <- cutree(hc, k = k_meta)
  meta_clusters <- micro_to_meta[micro_clusters]
  as.integer(meta_clusters)
}


#' Assign test cells to meta-clusters via micro-centroids
#'
#' 1) Assign each test cell to nearest micro-centroid \cr
#' 2) Map micro-cluster to meta-cluster via hclust cut
#'
#' @param test_data Scaled test cell matrix
#' @param centroids Micro-cluster centroids (n_micro × p)
#' @param hc hclust object on centroids
#' @param k_meta Number of meta-clusters
#' @return Integer vector of meta-cluster IDs for each test cell
two_stage_assign_test_meta <- function(test_data, centroids, hc, k_meta) {
  # Fast nearest-centroid assignment if available, otherwise fall back to the
  # (slow) euclidean_assignment() helper.
  if (requireNamespace("RANN", quietly = TRUE)) {
    micro_assign <- RANN::nn2(data = centroids, query = test_data, k = 1)$nn.idx[, 1]
  } else if (requireNamespace("FNN", quietly = TRUE)) {
    micro_assign <- FNN::get.knnx(centroids, test_data, k = 1)$nn.index[, 1]
  } else if (exists("euclidean_assignment")) {
    micro_assign <- euclidean_assignment(test_data, centroids)
  } else {
    micro_assign <- apply(test_data, 1, function(row) {
      dists <- apply(centroids, 1, function(centroid) sum((row - centroid)^2))
      which.min(dists)
    })
  }

  micro_to_meta <- cutree(hc, k = k_meta)
  meta_assign <- micro_to_meta[micro_assign]
  as.integer(meta_assign)
}

#
# Centroid graph meta-clustering (Leiden)
#

#' Build a kNN graph on centroids (for Leiden meta-clustering)
#'
#' @param centroids Numeric matrix (n_micro × p)
#' @param n_neighbors Number of neighbors for kNN graph
#' @return igraph graph object
two_stage_build_centroid_knn_graph <- function(centroids, n_neighbors = 30) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph package is required for Leiden meta-clustering.")
  }
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("RANN package is required to build a kNN graph efficiently.")
  }

  if (n_neighbors >= nrow(centroids)) {
    n_neighbors <- max(1, nrow(centroids) - 1)
  }

  # RANN includes the query point itself as the nearest neighbor.
  knn <- RANN::nn2(data = centroids, query = centroids, k = n_neighbors + 1)
  nn_index <- knn$nn.idx[, -1, drop = FALSE]
  edges <- cbind(
    rep(seq_len(nrow(centroids)), each = n_neighbors),
    as.vector(t(nn_index))
  )

  g <- igraph::graph_from_edgelist(edges, directed = FALSE)
  igraph::simplify(g)
}


#' Leiden meta-clustering on the centroid kNN graph
#'
#' @param centroid_graph igraph graph object (built on micro-centroids)
#' @param resolution Leiden resolution parameter (higher = more clusters)
#' @param seed Random seed for reproducibility
#' @return List with micro_to_meta mapping (length n_micro) and k_found
two_stage_leiden_micro_to_meta <- function(centroid_graph,
                                           resolution,
                                           seed = 123) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph package is required for Leiden meta-clustering.")
  }

  set.seed(seed)
  res <- igraph::cluster_leiden(
    centroid_graph,
    resolution = resolution,
    objective_function = "modularity"
  )

  membership_raw <- igraph::membership(res)
  micro_to_meta <- as.integer(factor(membership_raw))

  list(
    micro_to_meta = micro_to_meta,
    k_found = length(unique(micro_to_meta)),
    resolution = resolution
  )
}


#' Map per-cell micro cluster IDs to meta cluster IDs
#'
#' @param micro_clusters Integer vector of micro cluster IDs per cell (1..n_micro)
#' @param micro_to_meta Integer vector mapping micro ID → meta ID (length n_micro)
#' @return Integer vector of meta cluster IDs per cell
two_stage_meta_clusters_from_map <- function(micro_clusters, micro_to_meta) {
  as.integer(micro_to_meta[micro_clusters])
}


#' Assign test cells to nearest micro-centroids
#'
#' @param test_data Scaled test cell matrix
#' @param centroids Micro-centroid matrix (n_micro × p)
#' @return Integer vector of micro cluster IDs (1..n_micro) per test cell
two_stage_assign_test_micro <- function(test_data, centroids) {
  if (requireNamespace("RANN", quietly = TRUE)) {
    return(RANN::nn2(data = centroids, query = test_data, k = 1)$nn.idx[, 1])
  } else if (requireNamespace("FNN", quietly = TRUE)) {
    return(FNN::get.knnx(centroids, test_data, k = 1)$nn.index[, 1])
  }

  if (exists("euclidean_assignment")) {
    return(euclidean_assignment(test_data, centroids))
  }

  apply(test_data, 1, function(row) {
    dists <- apply(centroids, 1, function(centroid) sum((row - centroid)^2))
    which.min(dists)
  })
}


#' Build patient × micro-count matrix (dense)
#'
#' @param p_id Patient IDs per cell
#' @param micro_clusters Micro cluster IDs per cell (1..n_micro)
#' @param n_microclusters Total number of micro clusters
#' @return Integer matrix (n_patients × n_micro) with rownames = patient IDs
two_stage_patient_micro_counts <- function(p_id, micro_clusters, n_microclusters) {
  micro_factor <- factor(micro_clusters, levels = seq_len(n_microclusters))
  counts <- as.matrix(table(p_id, micro_factor))
  storage.mode(counts) <- "integer"
  colnames(counts) <- paste0("micro_", seq_len(n_microclusters))
  counts
}

#
# Full-data helpers (assign many cells to fixed centroids)
#

#' Compute patient × micro-counts from full data by nearest-centroid assignment
#'
#' This is used for strict outer-fold evaluation:
#' centroids are learned on outer-train only, and then all cells from a set of
#' patients are assigned to those centroids and aggregated.
#'
#' @param full_data Cell-level matrix (cells × markers)
#' @param full_p_id Patient ID vector (length = nrow(full_data))
#' @param patient_ids Character vector of patient IDs to include
#' @param centroids Micro-centroid matrix (n_micro × p) in scaled space
#' @param scale_center Numeric vector (length p) for scaling
#' @param scale_scale Numeric vector (length p) for scaling
#' @param chunk_size Number of cells per chunk
#' @param max_cells_per_patient Optional cap on cells per patient before assignment
#' @param seed Optional seed for patient-level subsampling
#' @param verbose Print progress
#' @return Integer matrix (patients × n_micro) with rownames = patient_ids
two_stage_patient_micro_counts_full <- function(full_data, full_p_id,
                                                patient_ids,
                                                centroids,
                                                scale_center,
                                                scale_scale,
                                                chunk_size = 500000,
                                                max_cells_per_patient = NULL,
                                                seed = NULL,
                                                verbose = FALSE) {
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("RANN package is required for nearest-centroid assignment.")
  }

  patient_ids <- as.character(patient_ids)
  n_patients <- length(patient_ids)
  n_micro <- nrow(centroids)

  if (n_patients == 0) {
    stop("patient_ids is empty.")
  }

  idx <- which(full_p_id %in% patient_ids)
  if (length(idx) == 0) {
    stop("No cells found for requested patient_ids.")
  }

  if (!is.null(max_cells_per_patient) && is.finite(max_cells_per_patient)) {
    max_cells_per_patient <- as.integer(max_cells_per_patient)
    if (max_cells_per_patient > 0) {
      if (!is.null(seed)) set.seed(seed)
      by_patient <- split(idx, full_p_id[idx])
      idx <- unlist(lapply(by_patient, function(ix) {
        if (length(ix) <= max_cells_per_patient) {
          ix
        } else {
          sample(ix, size = max_cells_per_patient)
        }
      }), use.names = FALSE)
    }
  }

  # Precompute mapping patient_id -> row index
  row_idx_all <- match(full_p_id[idx], patient_ids)
  if (any(is.na(row_idx_all))) {
    stop("Internal error: failed to match some patient IDs.")
  }

  nbins <- n_patients * n_micro
  counts_vec <- numeric(nbins)

  n_chunks <- ceiling(length(idx) / chunk_size)
  for (chunk in seq_len(n_chunks)) {
    start <- (chunk - 1) * chunk_size + 1
    end <- min(chunk * chunk_size, length(idx))
    sel <- start:end

    idx_chunk <- idx[sel]
    row_idx <- row_idx_all[sel]

    X_raw <- full_data[idx_chunk, , drop = FALSE]
    X_scaled <- scale(X_raw, center = scale_center, scale = scale_scale)

    micro_assign <- RANN::nn2(data = centroids, query = X_scaled, k = 1)$nn.idx[, 1]

    lin <- (row_idx - 1) * n_micro + micro_assign
    counts_vec <- counts_vec + tabulate(lin, nbins = nbins)

    if (verbose && (chunk %% max(1, floor(n_chunks / 10)) == 0)) {
      cat(sprintf("        Assigned chunk %d/%d (%d cells)\n",
                  chunk, n_chunks, end - start + 1))
    }
  }

  counts <- matrix(counts_vec, nrow = n_patients, ncol = n_micro, byrow = TRUE)
  storage.mode(counts) <- "integer"
  rownames(counts) <- patient_ids
  colnames(counts) <- paste0("micro_", seq_len(n_micro))
  counts
}


#' Collapse patient micro-counts into meta-counts via a micro→meta mapping
#'
#' @param micro_counts Integer matrix (patients × n_micro)
#' @param micro_to_meta Integer vector mapping micro→meta (length n_micro)
#' @return Integer matrix (patients × k_meta)
two_stage_collapse_micro_to_meta_counts <- function(micro_counts, micro_to_meta) {
  if (length(micro_to_meta) != ncol(micro_counts)) {
    stop("micro_to_meta length (", length(micro_to_meta),
         ") must equal ncol(micro_counts) (", ncol(micro_counts), ").")
  }

  meta_counts <- t(rowsum(t(micro_counts), group = micro_to_meta, reorder = FALSE))
  storage.mode(meta_counts) <- "integer"
  colnames(meta_counts) <- paste0("cluster_", seq_len(ncol(meta_counts)))
  meta_counts
}


#' Tune k on a fixed microclustering + dendrogram
#'
#' This mirrors tune_k_inside_outer() but:
#'   - Uses meta-clusters from two-stage micro + hclust
#'   - Reuses the same micro structure for all k in k_range
#'
#' @param train_data Scaled training cell data
#' @param train_p_id Training patient IDs per cell
#' @param patient_outcomes List with y, x, and optionally covariates
#' @param micro_clusters Micro-cluster assignments per training cell
#' @param hc hclust object on micro-centroids
#' @param k_range Vector of meta k values
#' @param lambda_grid Lambda grid for elasticnet
#' @param mtry_grid mtry grid for RF
#' @param nodesize_grid nodesize grid for RF
#' @param alpha Elastic net alpha
#' @param n_inner_folds Number of inner CV folds
#' @return List with best k, best hyperparameters, and inner CV curves
tune_k_on_microstructure <- function(train_data, train_p_id, patient_outcomes,
                                     micro_clusters, hc, k_range,
                                     prediction_method,
                                     lambda_grid = 10^seq(-3, 1, length.out = 20),
                                     mtry_grid = NULL,
                                     nodesize_grid = c(1, 5, 10),
                                     alpha = 0.5,
                                     n_inner_folds = 5) {

  best <- list(k = NA, inner_cv_error = Inf, hyper = NULL)

  inner_cv_validation_errors <- numeric(length(k_range))
  inner_cv_training_errors <- numeric(length(k_range))
  names(inner_cv_validation_errors) <- k_range
  names(inner_cv_training_errors) <- k_range

  # Patients with available outcomes
  train_p_ids_with_outcome <- unique(train_p_id)
  train_p_ids_with_outcome <- train_p_ids_with_outcome[
    train_p_ids_with_outcome %in% names(patient_outcomes$y)
  ]

  for (k_idx in seq_along(k_range)) {
    k_meta <- k_range[k_idx]

    # Meta-cluster labels for all training cells
    tr_clusters <- two_stage_meta_clusters(micro_clusters, hc, k_meta)

    # Build patient-level features (with covariates if present)
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
      gs <- grid_search_elasticnet(
        Xtr, ytr,
        alpha = alpha,
        lambda_grid = lambda_grid,
        n_inner_folds = n_inner_folds
      )
      inner_validation_err <- min(gs$cv_errors, na.rm = TRUE)
      hyper <- list(lambda = gs$best_lambda, alpha = alpha)

      model_train <- glmnet(Xtr, ytr,
                            alpha = hyper$alpha,
                            lambda = hyper$lambda)
      train_pred <- predict(model_train, newx = Xtr, s = hyper$lambda)
      inner_training_err <- mean((train_pred - ytr)^2)

    } else if (prediction_method == "random_forest") {
      gs <- grid_search_rf(
        Xtr, ytr,
        mtry_grid = mtry_grid,
        nodesize_grid = nodesize_grid,
        ntree = 500,
        n_inner_folds = max(3, floor(n_inner_folds / 2))
      )
      inner_validation_err <- min(gs$grid_results$cv_error, na.rm = TRUE)
      hyper <- list(
        mtry = gs$best_mtry,
        nodesize = gs$best_nodesize,
        ntree = 500
      )

      train_df <- as.data.frame(Xtr)
      train_df$y <- ytr
      model_train <- randomForest(
        y ~ ., data = train_df,
        ntree = hyper$ntree,
        mtry = hyper$mtry,
        nodesize = hyper$nodesize
      )
      train_pred <- predict(model_train, newdata = train_df)
      inner_training_err <- mean((train_pred - ytr)^2)

    } else {
      stop("Unsupported prediction_method in tune_k_on_microstructure(): ",
           prediction_method)
    }

    inner_cv_validation_errors[k_idx] <- inner_validation_err
    inner_cv_training_errors[k_idx] <- inner_training_err

    if (inner_validation_err < best$inner_cv_error) {
      best <- list(k = k_meta, inner_cv_error = inner_validation_err, hyper = hyper)
    }
  }

  best$inner_cv_validation_errors <- inner_cv_validation_errors
  best$inner_cv_training_errors <- inner_cv_training_errors
  best
}


#' Tune Leiden resolution on a fixed microclustering + centroid graph
#'
#' Analogous to tune_k_on_microstructure(), but:
#'   - Builds patient × micro-counts once
#'   - Clusters the centroids (micro clusters) via Leiden at each resolution
#'   - Collapses micro-counts into meta-counts using the Leiden membership
#'
#' @param train_p_id Training patient IDs per cell
#' @param patient_outcomes List with y, x, and optionally covariates
#' @param micro_clusters Micro cluster ID per training cell (1..n_microclusters)
#' @param n_microclusters Total number of micro clusters
#' @param centroid_graph igraph graph built on micro-centroids
#' @param resolution_grid Vector of Leiden resolution values to try
#' @param prediction_method "elasticnet" or "random_forest"
#' @param lambda_grid Lambda grid for elasticnet
#' @param mtry_grid mtry grid for RF
#' @param nodesize_grid nodesize grid for RF
#' @param alpha Elastic net alpha
#' @param n_inner_folds Number of inner CV folds
#' @return List with best resolution and inner CV curves
tune_resolution_on_microstructure <- function(train_p_id,
                                             patient_outcomes,
                                             micro_clusters,
                                             n_microclusters,
                                             centroid_graph,
                                             resolution_grid,
                                             prediction_method,
                                             lambda_grid = 10^seq(-3, 1, length.out = 20),
                                             mtry_grid = NULL,
                                             nodesize_grid = c(1, 5, 10),
                                             alpha = 0.5,
                                             n_inner_folds = 5) {

  best <- list(
    resolution = NA_real_,
    k_found_best = NA_integer_,
    inner_cv_error = Inf,
    hyper = NULL,
    micro_to_meta = NULL
  )

  # Restrict to patients with outcomes
  train_p_ids_with_outcome <- unique(train_p_id)
  train_p_ids_with_outcome <- train_p_ids_with_outcome[
    train_p_ids_with_outcome %in% names(patient_outcomes$y)
  ]

  # Precompute patient × micro-counts once (dense)
  micro_counts_all <- two_stage_patient_micro_counts(
    p_id = train_p_id,
    micro_clusters = micro_clusters,
    n_microclusters = n_microclusters
  )
  micro_counts_all <- micro_counts_all[train_p_ids_with_outcome, , drop = FALSE]

  inner_cv_validation_errors <- numeric(length(resolution_grid))
  inner_cv_training_errors <- numeric(length(resolution_grid))
  k_found_vec <- integer(length(resolution_grid))
  names(inner_cv_validation_errors) <- resolution_grid
  names(inner_cv_training_errors) <- resolution_grid
  names(k_found_vec) <- resolution_grid

  for (r_idx in seq_along(resolution_grid)) {
    resolution <- resolution_grid[r_idx]

    meta <- two_stage_leiden_micro_to_meta(
      centroid_graph = centroid_graph,
      resolution = resolution
    )
    micro_to_meta <- meta$micro_to_meta
    k_found_vec[r_idx] <- meta$k_found

    meta_counts <- two_stage_collapse_micro_to_meta_counts(
      micro_counts = micro_counts_all,
      micro_to_meta = micro_to_meta
    )

    # Features: log(meta_counts+1) + baseline (+ covariates)
    Xtr <- log(meta_counts + 1)
    Xtr <- cbind(
      Xtr,
      baseline = patient_outcomes$x[train_p_ids_with_outcome]
    )
    if (!is.null(patient_outcomes$covariates)) {
      covars <- patient_outcomes$covariates[train_p_ids_with_outcome, , drop = FALSE]
      Xtr <- cbind(Xtr, covars)
    }

    ytr <- patient_outcomes$y[rownames(Xtr)]

    if (prediction_method == "elasticnet") {
      gs <- grid_search_elasticnet(
        Xtr, ytr,
        alpha = alpha,
        lambda_grid = lambda_grid,
        n_inner_folds = n_inner_folds
      )
      inner_validation_err <- min(gs$cv_errors, na.rm = TRUE)
      hyper <- list(lambda = gs$best_lambda, alpha = alpha)

      model_train <- glmnet(Xtr, ytr,
                            alpha = hyper$alpha,
                            lambda = hyper$lambda)
      train_pred <- predict(model_train, newx = Xtr, s = hyper$lambda)
      inner_training_err <- mean((train_pred - ytr)^2)

    } else if (prediction_method == "random_forest") {
      gs <- grid_search_rf(
        Xtr, ytr,
        mtry_grid = mtry_grid,
        nodesize_grid = nodesize_grid,
        ntree = 500,
        n_inner_folds = max(3, floor(n_inner_folds / 2))
      )
      inner_validation_err <- min(gs$grid_results$cv_error, na.rm = TRUE)
      hyper <- list(
        mtry = gs$best_mtry,
        nodesize = gs$best_nodesize,
        ntree = 500
      )

      train_df <- as.data.frame(Xtr)
      train_df$y <- ytr
      model_train <- randomForest(
        y ~ ., data = train_df,
        ntree = hyper$ntree,
        mtry = hyper$mtry,
        nodesize = hyper$nodesize
      )
      train_pred <- predict(model_train, newdata = train_df)
      inner_training_err <- mean((train_pred - ytr)^2)

    } else {
      stop("Unsupported prediction_method in tune_resolution_on_microstructure(): ",
           prediction_method)
    }

    inner_cv_validation_errors[r_idx] <- inner_validation_err
    inner_cv_training_errors[r_idx] <- inner_training_err

    if (inner_validation_err < best$inner_cv_error) {
      best <- list(
        resolution = resolution,
        k_found_best = meta$k_found,
        inner_cv_error = inner_validation_err,
        hyper = hyper,
        micro_to_meta = micro_to_meta
      )
    }
  }

  best$inner_cv_validation_errors <- inner_cv_validation_errors
  best$inner_cv_training_errors <- inner_cv_training_errors
  best$k_found_by_resolution <- k_found_vec
  best
}


#' Proper nested CV using centroid microclustering + Leiden meta-clustering
#'
#' Outer-fold safe: micro-centroids and Leiden meta-clusters are fit on outer
#' training cells, and test cells are assigned by nearest centroid.
#'
#' @param mini_data Cell-level feature matrix (subsampled, scaled later)
#' @param p_id Vector mapping cells to patient IDs
#' @param patient_outcomes Outcome + baseline (+ covariates) list
#' @param resolution_grid Vector of Leiden resolution values to try
#' @param outer_folds List of outer CV folds (patient indices, as from createFolds)
#' @param prediction_method Prediction method ("elasticnet" or "random_forest")
#' @param micro_method Micro-clustering method ("kmeans" or "flowsom")
#' @param n_microclusters Number of micro-centroids
#' @param n_neighbors Number of neighbors for centroid kNN graph
#' @param lambda_grid Lambda grid for elasticnet
#' @param mtry_grid mtry grid for RF
#' @param nodesize_grid nodesize grid for RF
#' @param alpha Elastic net alpha
#' @param scale_features Scale features before clustering
#' @param verbose Print progress messages
#' @return List with mean_test_error, resolution_selected per fold, etc.
proper_nested_cv_two_stage_leiden <- function(mini_data, p_id, patient_outcomes,
                                              resolution_grid,
                                              outer_folds,
                                              prediction_method,
                                              micro_method = "kmeans",
                                              n_microclusters = 500,
                                              n_neighbors = 30,
                                              lambda_grid = 10^seq(-3, 1, length.out = 20),
                                              mtry_grid = NULL,
                                              nodesize_grid = c(1, 5, 10),
                                              alpha = 0.5,
                                              scale_features = TRUE,
                                              verbose = FALSE) {

  n_outer_folds <- length(outer_folds)
  unique_p_ids <- names(patient_outcomes$y)

  results <- list(
    test_errors = numeric(n_outer_folds),
    train_errors = numeric(n_outer_folds),
    resolution_selected = numeric(n_outer_folds),
    k_found = integer(n_outer_folds),
    test_predictions = list(),
    test_true = list(),
    inner_cv_training = list(),
    inner_cv_validation = list(),
    inner_k_found = list()
  )

  for (fold_id in 1:n_outer_folds) {
    if (verbose) {
      cat(sprintf(
        "    [Two-stage+Leiden] Fold %d/%d: tuning resolution in {%s}\n",
        fold_id, n_outer_folds, paste(resolution_grid, collapse = ",")
      ))
    }

    test_p_ids <- unique_p_ids[outer_folds[[fold_id]]]
    train_p_ids <- setdiff(unique_p_ids, test_p_ids)

    train_idx <- which(p_id %in% train_p_ids)
    test_idx <- which(p_id %in% test_p_ids)

    train_data_raw <- mini_data[train_idx, , drop = FALSE]
    test_data_raw <- mini_data[test_idx, , drop = FALSE]
    train_p_id <- p_id[train_idx]
    test_p_id <- p_id[test_idx]

    # Scale features inside outer fold
    if (scale_features) {
      sc <- scale(train_data_raw)
      train_data <- sc
      test_data <- scale(
        test_data_raw,
        center = attr(sc, "scaled:center"),
        scale = attr(sc, "scaled:scale")
      )
    } else {
      train_data <- train_data_raw
      test_data <- test_data_raw
    }

    # Microclustering on outer training data (skip hierarchical build)
    micro_res <- two_stage_microclustering(
      data = train_data,
      micro_method = micro_method,
      n_microclusters = n_microclusters,
      build_hc = FALSE,
      verbose = verbose
    )

    centroid_graph <- two_stage_build_centroid_knn_graph(
      centroids = micro_res$centroids,
      n_neighbors = n_neighbors
    )

    # Tune resolution using fixed microstructure + centroid graph
    best <- tune_resolution_on_microstructure(
      train_p_id = train_p_id,
      patient_outcomes = patient_outcomes,
      micro_clusters = micro_res$micro_clusters,
      n_microclusters = n_microclusters,
      centroid_graph = centroid_graph,
      resolution_grid = resolution_grid,
      prediction_method = prediction_method,
      lambda_grid = lambda_grid,
      mtry_grid = mtry_grid,
      nodesize_grid = nodesize_grid,
      alpha = alpha,
      n_inner_folds = 5
    )

    best_resolution <- best$resolution
    best_hyper <- best$hyper
    micro_to_meta <- best$micro_to_meta

    results$resolution_selected[fold_id] <- best_resolution
    results$k_found[fold_id] <- best$k_found_best
    results$inner_cv_training[[fold_id]] <- best$inner_cv_training_errors
    results$inner_cv_validation[[fold_id]] <- best$inner_cv_validation_errors
    results$inner_k_found[[fold_id]] <- best$k_found_by_resolution

    # --- Refit on full outer training set with best resolution ---
    train_p_ids_with_outcome <- train_p_ids[train_p_ids %in% names(patient_outcomes$y)]
    test_p_ids_with_outcome <- test_p_ids[test_p_ids %in% names(patient_outcomes$y)]

    # Train patient × micro counts
    train_micro_counts <- two_stage_patient_micro_counts(
      p_id = train_p_id,
      micro_clusters = micro_res$micro_clusters,
      n_microclusters = n_microclusters
    )
    train_micro_counts <- train_micro_counts[train_p_ids_with_outcome, , drop = FALSE]

    train_meta_counts <- two_stage_collapse_micro_to_meta_counts(
      micro_counts = train_micro_counts,
      micro_to_meta = micro_to_meta
    )

    train_features <- log(train_meta_counts + 1)
    train_features <- cbind(
      train_features,
      baseline = patient_outcomes$x[train_p_ids_with_outcome]
    )
    if (!is.null(patient_outcomes$covariates)) {
      train_covariates <- patient_outcomes$covariates[train_p_ids_with_outcome, , drop = FALSE]
      train_features <- cbind(train_features, train_covariates)
    }

    # Test patient × micro counts (assign to nearest centroids)
    test_micro_assign <- two_stage_assign_test_micro(
      test_data = test_data,
      centroids = micro_res$centroids
    )
    test_micro_counts <- two_stage_patient_micro_counts(
      p_id = test_p_id,
      micro_clusters = test_micro_assign,
      n_microclusters = n_microclusters
    )
    test_micro_counts <- test_micro_counts[test_p_ids_with_outcome, , drop = FALSE]

    test_meta_counts <- two_stage_collapse_micro_to_meta_counts(
      micro_counts = test_micro_counts,
      micro_to_meta = micro_to_meta
    )

    test_features <- log(test_meta_counts + 1)
    test_features <- cbind(
      test_features,
      baseline = patient_outcomes$x[test_p_ids_with_outcome]
    )
    if (!is.null(patient_outcomes$covariates)) {
      test_covariates <- patient_outcomes$covariates[test_p_ids_with_outcome, , drop = FALSE]
      test_features <- cbind(test_features, test_covariates)
    }

    train_features <- as.matrix(train_features)
    test_features <- as.matrix(test_features)

    train_y <- patient_outcomes$y[rownames(train_features)]
    test_y <- patient_outcomes$y[rownames(test_features)]

    # Fit final model with best hyperparameters
    if (prediction_method == "elasticnet") {
      model <- glmnet(
        x = train_features,
        y = train_y,
        alpha = best_hyper$alpha,
        lambda = best_hyper$lambda
      )
      test_pred <- predict(model, newx = test_features, s = best_hyper$lambda)
      train_pred <- predict(model, newx = train_features, s = best_hyper$lambda)

    } else if (prediction_method == "random_forest") {
      train_df <- as.data.frame(train_features)
      train_df$y <- train_y
      model <- randomForest(
        y ~ ., data = train_df,
        ntree = best_hyper$ntree,
        mtry = best_hyper$mtry,
        nodesize = best_hyper$nodesize
      )
      test_pred <- predict(model, newdata = as.data.frame(test_features))
      train_pred <- predict(model, newdata = train_df)
    }

    results$test_errors[fold_id] <- mean((as.numeric(test_pred) - test_y)^2)
    results$train_errors[fold_id] <- mean((as.numeric(train_pred) - train_y)^2)
    results$test_predictions[[fold_id]] <- as.numeric(test_pred)
    results$test_true[[fold_id]] <- as.numeric(test_y)

    if (verbose) {
      cat("      [two-stage+Leiden] Fold", fold_id,
          "res=", best_resolution,
          "k=", best$k_found_best,
          " Test:", round(results$test_errors[fold_id], 4),
          " Train:", round(results$train_errors[fold_id], 4), "\n")
    }
  }

  list(
    mean_test_error = mean(results$test_errors),
    se_test_error = sd(results$test_errors) / sqrt(n_outer_folds),
    mean_train_error = mean(results$train_errors),
    se_train_error = sd(results$train_errors) / sqrt(n_outer_folds),
    # For compatibility with aggregate_proper_nested_cv()
    k_selected = results$k_found,
    resolution_selected = results$resolution_selected,
    k_found = results$k_found,
    test_errors_per_fold = results$test_errors,
    train_errors_per_fold = results$train_errors,
    test_predictions = results$test_predictions,
    test_true = results$test_true,
    inner_cv_training = results$inner_cv_training,
    inner_cv_validation = results$inner_cv_validation,
    inner_k_found = results$inner_k_found
  )
}


#' Proper nested CV using centroid microclustering + Leiden with full-cell aggregation
#'
#' Strict outer-fold pipeline:
#'   - Fit scaling + k-means centroids on a SAMPLE of outer-train cells
#'   - Build centroid kNN graph and tune Leiden resolution inside outer-train
#'   - Assign ALL cells (train/test patients) to centroids and aggregate features
#'
#' This avoids the "empty clusters due to subsampling" issue by using all cells
#' for patient-level counts while keeping the representation leakage-safe.
proper_nested_cv_two_stage_leiden_fullcells <- function(full_data, full_p_id,
                                                        patient_outcomes,
                                                        resolution_grid,
                                                        outer_folds,
                                                        prediction_method,
                                                        micro_method = "kmeans",
                                                        n_microclusters = 2000,
                                                        n_neighbors = 30,
                                                        kmeans_fit_n_cells = 200000,
                                                        kmeans_nstart = 1,
                                                        kmeans_iter_max = 50,
                                                        assignment_chunk_size = 500000,
                                                        lambda_grid = 10^seq(-3, 1, length.out = 20),
                                                        mtry_grid = NULL,
                                                        nodesize_grid = c(1, 5, 10),
                                                        alpha = 0.5,
                                                        verbose = FALSE) {

  if (micro_method != "kmeans") {
    stop("Currently only micro_method='kmeans' is supported in fullcells mode.")
  }

  n_outer_folds <- length(outer_folds)
  unique_p_ids <- names(patient_outcomes$y)

  results <- list(
    test_errors = numeric(n_outer_folds),
    train_errors = numeric(n_outer_folds),
    resolution_selected = numeric(n_outer_folds),
    k_found = integer(n_outer_folds),
    test_predictions = list(),
    test_true = list(),
    inner_cv_training = list(),
    inner_cv_validation = list(),
    inner_k_found = list()
  )

  for (fold_id in 1:n_outer_folds) {
    if (verbose) {
      cat(sprintf(
        "    [Strict fullcells] Fold %d/%d\n",
        fold_id, n_outer_folds
      ))
    }

    test_p_ids <- unique_p_ids[outer_folds[[fold_id]]]
    train_p_ids <- setdiff(unique_p_ids, test_p_ids)

    # Fit scaling + kmeans on a training-only sample
    set.seed(1000 + fold_id)
    train_cell_idx <- which(full_p_id %in% train_p_ids)
    fit_n <- min(kmeans_fit_n_cells, length(train_cell_idx))
    fit_idx <- sample(train_cell_idx, size = fit_n)

    X_fit <- full_data[fit_idx, , drop = FALSE]
    scale_center <- colMeans(X_fit)
    scale_scale <- apply(X_fit, 2, sd)
    scale_scale[scale_scale == 0] <- 1
    X_fit_scaled <- scale(X_fit, center = scale_center, scale = scale_scale)

    if (verbose) {
      cat(sprintf(
        "      Fitting k-means: n_fit=%d, micro=%d, nstart=%d\n",
        nrow(X_fit_scaled), n_microclusters, kmeans_nstart
      ))
    }

    set.seed(123 + fold_id)
    km <- kmeans(
      X_fit_scaled,
      centers = n_microclusters,
      nstart = kmeans_nstart,
      iter.max = kmeans_iter_max
    )
    centroids <- km$centers

    centroid_graph <- two_stage_build_centroid_knn_graph(
      centroids = centroids,
      n_neighbors = n_neighbors
    )

    # Build patient × micro-counts for TRAIN patients (all their cells)
    if (verbose) cat("      Assigning all training cells to centroids...\n")
    train_micro_counts <- two_stage_patient_micro_counts_full(
      full_data = full_data,
      full_p_id = full_p_id,
      patient_ids = train_p_ids,
      centroids = centroids,
      scale_center = scale_center,
      scale_scale = scale_scale,
      chunk_size = assignment_chunk_size,
      verbose = verbose
    )

    # NOTE: tune_resolution_on_microstructure() expects per-cell micro_clusters.
    # In fullcells mode we instead precompute train_micro_counts, so we tune
    # resolution by reusing these counts (no test touched).
    best <- list(
      resolution = NA_real_,
      k_found_best = NA_integer_,
      inner_cv_error = Inf,
      hyper = NULL,
      micro_to_meta = NULL,
      inner_cv_validation_errors = numeric(length(resolution_grid)),
      inner_cv_training_errors = numeric(length(resolution_grid)),
      k_found_by_resolution = integer(length(resolution_grid))
    )
    names(best$inner_cv_validation_errors) <- resolution_grid
    names(best$inner_cv_training_errors) <- resolution_grid
    names(best$k_found_by_resolution) <- resolution_grid

    train_p_ids_with_outcome <- train_p_ids[train_p_ids %in% names(patient_outcomes$y)]
    train_micro_counts <- train_micro_counts[train_p_ids_with_outcome, , drop = FALSE]

    for (r_idx in seq_along(resolution_grid)) {
      resolution <- resolution_grid[r_idx]
      meta <- two_stage_leiden_micro_to_meta(centroid_graph, resolution = resolution)
      micro_to_meta <- meta$micro_to_meta
      best$k_found_by_resolution[r_idx] <- meta$k_found

      train_meta_counts <- two_stage_collapse_micro_to_meta_counts(
        micro_counts = train_micro_counts,
        micro_to_meta = micro_to_meta
      )

      Xtr <- log(train_meta_counts + 1)
      Xtr <- cbind(Xtr, baseline = patient_outcomes$x[train_p_ids_with_outcome])
      if (!is.null(patient_outcomes$covariates)) {
        covars <- patient_outcomes$covariates[train_p_ids_with_outcome, , drop = FALSE]
        Xtr <- cbind(Xtr, covars)
      }
      ytr <- patient_outcomes$y[rownames(Xtr)]

      if (prediction_method == "elasticnet") {
        gs <- grid_search_elasticnet(
          Xtr, ytr,
          alpha = alpha,
          lambda_grid = lambda_grid,
          n_inner_folds = 5
        )
        inner_validation_err <- min(gs$cv_errors, na.rm = TRUE)
        hyper <- list(lambda = gs$best_lambda, alpha = alpha)

        model_train <- glmnet(Xtr, ytr, alpha = hyper$alpha, lambda = hyper$lambda)
        train_pred <- predict(model_train, newx = Xtr, s = hyper$lambda)
        inner_training_err <- mean((train_pred - ytr)^2)
      } else if (prediction_method == "random_forest") {
        gs <- grid_search_rf(
          Xtr, ytr,
          mtry_grid = mtry_grid,
          nodesize_grid = nodesize_grid,
          ntree = 500,
          n_inner_folds = 3
        )
        inner_validation_err <- min(gs$grid_results$cv_error, na.rm = TRUE)
        hyper <- list(mtry = gs$best_mtry, nodesize = gs$best_nodesize, ntree = 500)

        train_df <- as.data.frame(Xtr)
        train_df$y <- ytr
        model_train <- randomForest(y ~ ., data = train_df,
                                    ntree = hyper$ntree,
                                    mtry = hyper$mtry,
                                    nodesize = hyper$nodesize)
        train_pred <- predict(model_train, newdata = train_df)
        inner_training_err <- mean((train_pred - ytr)^2)
      } else {
        stop("Unsupported prediction_method in fullcells mode: ", prediction_method)
      }

      best$inner_cv_validation_errors[r_idx] <- inner_validation_err
      best$inner_cv_training_errors[r_idx] <- inner_training_err

      if (inner_validation_err < best$inner_cv_error) {
        best$resolution <- resolution
        best$k_found_best <- meta$k_found
        best$inner_cv_error <- inner_validation_err
        best$hyper <- hyper
        best$micro_to_meta <- micro_to_meta
      }
    }

    # Store tuning curves
    results$resolution_selected[fold_id] <- best$resolution
    results$k_found[fold_id] <- best$k_found_best
    results$inner_cv_training[[fold_id]] <- best$inner_cv_training_errors
    results$inner_cv_validation[[fold_id]] <- best$inner_cv_validation_errors
    results$inner_k_found[[fold_id]] <- best$k_found_by_resolution

    # Build final TRAIN features
    train_meta_counts <- two_stage_collapse_micro_to_meta_counts(
      micro_counts = train_micro_counts,
      micro_to_meta = best$micro_to_meta
    )
    X_train <- log(train_meta_counts + 1)
    X_train <- cbind(X_train, baseline = patient_outcomes$x[train_p_ids_with_outcome])
    if (!is.null(patient_outcomes$covariates)) {
      covars <- patient_outcomes$covariates[train_p_ids_with_outcome, , drop = FALSE]
      X_train <- cbind(X_train, covars)
    }

    # Build TEST features from ALL test cells
    if (verbose) cat("      Assigning all test cells to centroids...\n")
    test_p_ids_with_outcome <- test_p_ids[test_p_ids %in% names(patient_outcomes$y)]
    test_micro_counts <- two_stage_patient_micro_counts_full(
      full_data = full_data,
      full_p_id = full_p_id,
      patient_ids = test_p_ids_with_outcome,
      centroids = centroids,
      scale_center = scale_center,
      scale_scale = scale_scale,
      chunk_size = assignment_chunk_size,
      verbose = verbose
    )
    test_meta_counts <- two_stage_collapse_micro_to_meta_counts(
      micro_counts = test_micro_counts,
      micro_to_meta = best$micro_to_meta
    )
    X_test <- log(test_meta_counts + 1)
    X_test <- cbind(X_test, baseline = patient_outcomes$x[test_p_ids_with_outcome])
    if (!is.null(patient_outcomes$covariates)) {
      covars <- patient_outcomes$covariates[test_p_ids_with_outcome, , drop = FALSE]
      X_test <- cbind(X_test, covars)
    }

    X_train <- as.matrix(X_train)
    X_test <- as.matrix(X_test)
    y_train <- patient_outcomes$y[rownames(X_train)]
    y_test <- patient_outcomes$y[rownames(X_test)]

    # Fit final model with tuned hyperparameters
    if (prediction_method == "elasticnet") {
      model <- glmnet(X_train, y_train,
                      alpha = best$hyper$alpha,
                      lambda = best$hyper$lambda)
      pred_test <- predict(model, newx = X_test, s = best$hyper$lambda)
      pred_train <- predict(model, newx = X_train, s = best$hyper$lambda)
    } else if (prediction_method == "random_forest") {
      train_df <- as.data.frame(X_train)
      train_df$y <- y_train
      model <- randomForest(
        y ~ ., data = train_df,
        ntree = best$hyper$ntree,
        mtry = best$hyper$mtry,
        nodesize = best$hyper$nodesize
      )
      pred_test <- predict(model, newdata = as.data.frame(X_test))
      pred_train <- predict(model, newdata = train_df)
    }

    results$test_errors[fold_id] <- mean((as.numeric(pred_test) - y_test)^2)
    results$train_errors[fold_id] <- mean((as.numeric(pred_train) - y_train)^2)
    results$test_predictions[[fold_id]] <- as.numeric(pred_test)
    results$test_true[[fold_id]] <- as.numeric(y_test)

    if (verbose) {
      cat(sprintf(
        "      Fold %d: res=%.3f k=%d test=%.4f train=%.4f\n",
        fold_id, best$resolution, best$k_found_best,
        results$test_errors[fold_id], results$train_errors[fold_id]
      ))
    }
  }

  list(
    mean_test_error = mean(results$test_errors),
    se_test_error = sd(results$test_errors) / sqrt(n_outer_folds),
    mean_train_error = mean(results$train_errors),
    se_train_error = sd(results$train_errors) / sqrt(n_outer_folds),
    # For compatibility with aggregation helpers
    k_selected = results$k_found,
    resolution_selected = results$resolution_selected,
    k_found = results$k_found,
    test_errors_per_fold = results$test_errors,
    train_errors_per_fold = results$train_errors,
    test_predictions = results$test_predictions,
    test_true = results$test_true,
    inner_cv_training = results$inner_cv_training,
    inner_cv_validation = results$inner_cv_validation,
    inner_k_found = results$inner_k_found
  )
}


#' Proper nested CV using two-stage clustering
#'
#' This mirrors proper_nested_cv() but replaces direct clustering with:
#'   - Micro-clustering on training data
#'   - Hierarchical clustering on micro-centroids
#'   - Meta-clusters obtained by cutting the tree at k
#'
#' @param mini_data Cell-level feature matrix (subsampled, scaled later)
#' @param p_id Vector mapping cells to patient IDs
#' @param patient_outcomes Outcome + baseline (+ covariates) list
#' @param k_range Vector of meta-cluster counts to try
#' @param outer_folds List of outer CV folds (patient indices, as from createFolds)
#' @param prediction_method Prediction method ("elasticnet" or "random_forest")
#' @param micro_method Micro-clustering method ("kmeans" or "flowsom")
#' @param n_microclusters Number of micro-clusters
#' @param lambda_grid Lambda grid for elasticnet
#' @param mtry_grid mtry grid for RF
#' @param nodesize_grid nodesize grid for RF
#' @param alpha Elastic net alpha
#' @param scale_features Scale features before clustering
#' @param verbose Print progress messages
#' @return List with mean_test_error, k_selected per fold, etc.
proper_nested_cv_two_stage <- function(mini_data, p_id, patient_outcomes,
                                       k_range,
                                       outer_folds,
                                       prediction_method,
                                       micro_method = "kmeans",
                                       n_microclusters = 500,
                                       lambda_grid = 10^seq(-3, 1, length.out = 20),
                                       mtry_grid = NULL,
                                       nodesize_grid = c(1, 5, 10),
                                       alpha = 0.5,
                                       scale_features = TRUE,
                                       verbose = FALSE) {

  n_outer_folds <- length(outer_folds)
  unique_p_ids <- names(patient_outcomes$y)

  results <- list(
    test_errors = numeric(n_outer_folds),
    train_errors = numeric(n_outer_folds),
    k_selected = integer(n_outer_folds),
    test_predictions = list(),
    test_true = list(),
    inner_cv_training = list(),
    inner_cv_validation = list()
  )

  for (fold_id in 1:n_outer_folds) {
    if (verbose) {
      cat(sprintf(
        "    [Two-stage] Fold %d/%d: tuning k in {%s}\n",
        fold_id, n_outer_folds, paste(k_range, collapse = ",")
      ))
    }

    test_p_ids <- unique_p_ids[outer_folds[[fold_id]]]
    train_p_ids <- setdiff(unique_p_ids, test_p_ids)

    train_idx <- which(p_id %in% train_p_ids)
    test_idx <- which(p_id %in% test_p_ids)

    train_data_raw <- mini_data[train_idx, , drop = FALSE]
    test_data_raw <- mini_data[test_idx, , drop = FALSE]
    train_p_id <- p_id[train_idx]
    test_p_id <- p_id[test_idx]

    # Scale features inside outer fold
    if (scale_features) {
      sc <- scale(train_data_raw)
      train_data <- sc
      test_data <- scale(
        test_data_raw,
        center = attr(sc, "scaled:center"),
        scale = attr(sc, "scaled:scale")
      )
    } else {
      train_data <- train_data_raw
      test_data <- test_data_raw
    }

    # Microclustering on outer training data
    micro_res <- two_stage_microclustering(
      data = train_data,
      micro_method = micro_method,
      n_microclusters = n_microclusters,
      verbose = verbose
    )

    # Tune k using this microstructure
    best <- tune_k_on_microstructure(
      train_data = train_data,
      train_p_id = train_p_id,
      patient_outcomes = patient_outcomes,
      micro_clusters = micro_res$micro_clusters,
      hc = micro_res$hc,
      k_range = k_range,
      prediction_method = prediction_method,
      lambda_grid = lambda_grid,
      mtry_grid = mtry_grid,
      nodesize_grid = nodesize_grid,
      alpha = alpha,
      n_inner_folds = 5
    )

    best_k <- best$k
    best_hyper <- best$hyper
    results$k_selected[fold_id] <- best_k
    results$inner_cv_training[[fold_id]] <- best$inner_cv_training_errors
    results$inner_cv_validation[[fold_id]] <- best$inner_cv_validation_errors

    # --- Refit on full outer training set with best_k ---
    # Meta-cluster labels for all training cells
    train_clusters <- two_stage_meta_clusters(
      micro_clusters = micro_res$micro_clusters,
      hc = micro_res$hc,
      k_meta = best_k
    )

    # Build training features (with covariates if available)
    train_covariates <- if (!is.null(patient_outcomes$covariates)) {
      patient_outcomes$covariates[train_p_ids, , drop = FALSE]
    } else {
      NULL
    }

    if (!is.null(train_covariates)) {
      train_features <- build_patient_features_with_covariates(
        p_id = train_p_id,
        clusters = train_clusters,
        baseline_outcome = patient_outcomes$x[train_p_ids],
        covariates = train_covariates
      )
    } else {
      train_features <- build_patient_features(
        p_id = train_p_id,
        clusters = train_clusters,
        baseline_outcome = patient_outcomes$x[train_p_ids]
      )
    }

    # Assign test cells via micro-centroids then map to meta-clusters
    test_clusters <- two_stage_assign_test_meta(
      test_data = test_data,
      centroids = micro_res$centroids,
      hc = micro_res$hc,
      k_meta = best_k
    )

    # Build test features
    test_counts <- as.data.frame.matrix(table(test_p_id, test_clusters))
    test_counts <- log(test_counts + 1)

    non_cluster_cols <- c("baseline")
    if (!is.null(train_covariates)) {
      non_cluster_cols <- c(non_cluster_cols, colnames(train_covariates))
    }
    train_cluster_cols <- setdiff(colnames(train_features), non_cluster_cols)

    test_features_aligned <- align_features(
      test_features = test_counts,
      train_colnames = train_cluster_cols
    )

    test_patient_ids <- rownames(test_features_aligned)
    test_baseline <- patient_outcomes$x[test_patient_ids]
    test_features <- cbind(
      as.matrix(test_features_aligned),
      baseline = test_baseline
    )

    if (!is.null(train_covariates)) {
      test_covariates <- patient_outcomes$covariates[test_patient_ids, , drop = FALSE]
      test_features <- cbind(test_features, test_covariates)
    }

    train_features <- as.matrix(train_features)
    test_features <- as.matrix(test_features)

    train_y <- patient_outcomes$y[rownames(train_features)]
    test_y <- patient_outcomes$y[rownames(test_features)]

    # Fit final model with best hyperparameters
    if (prediction_method == "elasticnet") {
      model <- glmnet(
        x = train_features,
        y = train_y,
        alpha = best_hyper$alpha,
        lambda = best_hyper$lambda
      )
      test_pred <- predict(model, newx = test_features, s = best_hyper$lambda)
      train_pred <- predict(model, newx = train_features, s = best_hyper$lambda)

    } else if (prediction_method == "random_forest") {
      train_df <- as.data.frame(train_features)
      train_df$y <- train_y
      model <- randomForest(
        y ~ ., data = train_df,
        ntree = best_hyper$ntree,
        mtry = best_hyper$mtry,
        nodesize = best_hyper$nodesize
      )
      test_pred <- predict(model, newdata = as.data.frame(test_features))
      train_pred <- predict(model, newdata = train_df)
    }

    results$test_errors[fold_id] <- mean((as.numeric(test_pred) - test_y)^2)
    results$train_errors[fold_id] <- mean((as.numeric(train_pred) - train_y)^2)
    results$test_predictions[[fold_id]] <- as.numeric(test_pred)
    results$test_true[[fold_id]] <- as.numeric(test_y)

    if (verbose) {
      cat("      [two-stage] Fold", fold_id,
          "k=", best_k,
          " Test:", round(results$test_errors[fold_id], 4),
          " Train:", round(results$train_errors[fold_id], 4), "\n")
    }
  }

  list(
    mean_test_error = mean(results$test_errors),
    se_test_error = sd(results$test_errors) / sqrt(n_outer_folds),
    mean_train_error = mean(results$train_errors),
    se_train_error = sd(results$train_errors) / sqrt(n_outer_folds),
    k_selected = results$k_selected,
    test_errors_per_fold = results$test_errors,
    train_errors_per_fold = results$train_errors,
    test_predictions = results$test_predictions,
    test_true = results$test_true,
    inner_cv_training = results$inner_cv_training,
    inner_cv_validation = results$inner_cv_validation
  )
}


# NOTE:
# The full two-stage nested CV, including prediction_method handling and
# aggregation across subsamples, is implemented in
# proper_nested_cv_two_stage_with_subsampling() below. That function manages:
#   - prediction_method selection
#   - micro-method / n_microclusters configuration
#   - parallel execution across subsamples
#   - aggregation consistent with existing proper_nested_cv_with_subsampling()


#' Proper nested CV with subsampling (two-stage clustering)
#'
#' Wrapper around proper_nested_cv_two_stage to:
#'   - Perform patient-stratified subsampling of cells
#'   - Run nested CV for each prediction method
#'   - Optionally run in parallel over subsamples
#'   - Aggregate results across subsamples
#'
#' @param full_data Full cell-level feature matrix
#' @param full_p_id Patient ID for each cell
#' @param outcome_obj Outcome data frame
#' @param meta_obj Metadata with covariates (optional)
#' @param outcome_col_1 Baseline column (v1)
#' @param outcome_col_2 Target column (v3)
#' @param prediction_methods Vector of prediction method names
#' @param k_range Meta-cluster sizes to test
#' @param n_cells Target cells per subsample
#' @param n_subsamples Number of subsamples
#' @param n_outer_folds Number of outer CV folds
#' @param micro_method Micro-clustering method ("kmeans" or "flowsom")
#' @param n_microclusters Number of micro-clusters
#' @param lambda_grid Lambda grid (elasticnet)
#' @param alpha Elastic net alpha
#' @param include_dose Include dose covariate
#' @param include_age Include age covariate
#' @param include_sex Include sex covariate
#' @param include_time Include time covariate
#' @param scale_features Scale features before clustering
#' @param patient_stratified_sampling Patient-stratified subsampling
#' @param outcome_type "fold_change" or "absolute"
#' @param n_cores Number of parallel cores (1 = sequential)
#' @param verbose Print progress
#' @return List with raw_results and aggregated (same structure as
#'         proper_nested_cv_with_subsampling())
proper_nested_cv_two_stage_with_subsampling <- function(
    full_data, full_p_id, outcome_obj,
    meta_obj = NULL,
    outcome_col_1, outcome_col_2,
    prediction_methods = c("elasticnet", "random_forest"),
    k_range = seq(20, 40, by = 10),
    n_cells = 1000000,
    n_subsamples = 3,
    n_outer_folds = 5,
    micro_method = "kmeans",
    n_microclusters = 500,
    lambda_grid = 10^seq(-3, 1, length.out = 20),
    alpha = 0.5,
    include_dose = TRUE,
    include_age = TRUE,
    include_sex = TRUE,
    include_time = TRUE,
    scale_features = TRUE,
    patient_stratified_sampling = TRUE,
    outcome_type = "fold_change",
    n_cores = 1,
    verbose = TRUE) {

  if (verbose) {
    cat("\n", rep("=", 70), "\n", sep = "")
    cat("TWO-STAGE PROPER NESTED CV WITH SUBSAMPLING\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("Parameters:\n")
    cat("  - n_cells per subsample:", n_cells, "\n")
    cat("  - n_subsamples:", n_subsamples, "\n")
    cat("  - k_range (meta clusters):", paste(k_range, collapse = ", "), "\n")
    cat("  - micro_method:", micro_method, "\n")
    cat("  - n_microclusters:", n_microclusters, "\n")
    cat("  - prediction_methods:", paste(prediction_methods, collapse = ", "), "\n")
    cat("  - n_outer_folds:", n_outer_folds, "\n")
    cat("  - n_cores:", n_cores, "\n")
    cat("  - outcome_type:", outcome_type, "\n")
    covars <- c()
    if (include_dose) covars <- c(covars, "dose")
    if (include_age) covars <- c(covars, "age")
    if (include_sex) covars <- c(covars, "sex")
    if (include_time) covars <- c(covars, "time")
    cat("  - Covariates:", paste(covars, collapse = ", "), "\n\n")
  }

  # Prepare outcomes once
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
      patient_outcomes <- prepare_outcome_absolute(
        outcome_obj, outcome_col_1, outcome_col_2, full_p_id
      )
    } else {
      patient_outcomes <- prepare_outcome(
        outcome_obj, outcome_col_1, outcome_col_2, full_p_id
      )
    }
  }

  # Label for clustering method (for aggregation)
  clustering_label <- paste0("two_stage_", micro_method, "_micro", n_microclusters)

  # Storage
  results <- list()

  # Helper to process one subsample (sequential, called directly or via foreach)
  process_subsample <- function(subsample_idx) {
    if (verbose && n_cores == 1) {
      cat("\n", rep("-", 70), "\n", sep = "")
      cat("SUBSAMPLE ", subsample_idx, " / ", n_subsamples, "\n", sep = "")
      cat(rep("-", 70), "\n", sep = "")
    }

    set.seed(2025 + subsample_idx)

    # Patient-stratified subsampling
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

    # Outer CV folds (patient-level) for this subsample
    set.seed(2025 + subsample_idx)
    unique_p_ids <- names(patient_outcomes$y)
    outer_folds <- createFolds(seq_along(unique_p_ids),
                               k = n_outer_folds,
                               list = TRUE)

    subsample_results <- list()

    # Two-stage methods (one clustering label, multiple prediction methods)
    for (pred_method in prediction_methods) {
      method_name <- paste(clustering_label, pred_method, sep = "_")

      # Progress logging
      cat(sprintf(
        "[Two-stage] Subsample %d/%d, Method: %s\n",
        subsample_idx, n_subsamples, pred_method
      ))

      cv_res <- proper_nested_cv_two_stage(
        mini_data = mini_data,
        p_id = p_id,
        patient_outcomes = patient_outcomes,
        k_range = k_range,
        outer_folds = outer_folds,
        prediction_method = pred_method,
        micro_method = micro_method,
        n_microclusters = n_microclusters,
        lambda_grid = lambda_grid,
        mtry_grid = NULL,
        nodesize_grid = c(1, 5, 10),
        alpha = alpha,
        scale_features = scale_features,
        verbose = verbose  # Pass through verbose from parent
      )

      subsample_results[[method_name]] <- cv_res
    }

    # Baseline-only controls
    for (pred_method in prediction_methods) {
      method_name <- paste0("baseline_only_", pred_method)

      base_res <- baseline_only_cv(
        patient_outcomes = patient_outcomes,
        outer_folds = outer_folds,
        prediction_method = pred_method,
        lambda_grid = lambda_grid,
        mtry_grid = NULL,
        nodesize_grid = c(1, 5, 10),
        alpha = alpha,
        verbose = FALSE
      )

      subsample_results[[method_name]] <- base_res
    }

    subsample_results
  }

  # Run subsamples (optionally in parallel)
  if (n_cores > 1) {
    if (verbose) {
      cat("Running", n_subsamples, "subsamples in parallel using",
          n_cores, "cores...\n")
    }

    library(parallel)
    library(doParallel)
    library(foreach)

    cl <- makeCluster(n_cores)
    registerDoParallel(cl)

    clusterEvalQ(cl, {
      library(here)
      source(here::here("methods", "00-paths.R"))
      source(here::here("methods", "02-data-preparation.R"))
      source(here::here("methods", "03-clustering-methods.R"))
      source(here::here("methods", "04-prediction-models.R"))
      source(here::here("methods", "07-nested-cv-pipeline.R"))
      source(here::here("methods", "08-multi-outcome-covariates.R"))
      source(here::here("methods", "09-proper-nested-cv.R"))
      source(here::here("methods", "10-two-stage-clustering.R"))
    })

    clusterExport(
      cl,
      c("full_data", "full_p_id", "patient_outcomes", "n_cells",
        "patient_stratified_sampling", "n_outer_folds",
        "k_range", "lambda_grid", "alpha", "scale_features",
        "prediction_methods", "micro_method", "n_microclusters",
        "baseline_only_cv", "proper_nested_cv_two_stage",
        "two_stage_microclustering", "two_stage_meta_clusters",
        "two_stage_assign_test_meta", "tune_k_on_microstructure",
        "clustering_label"),
      envir = environment()
    )

    subsample_results_list <- foreach(subsample_idx = 1:n_subsamples,
                                      .packages = c("glmnet", "randomForest", "caret", "FlowSOM")) %dopar% {
      process_subsample(subsample_idx)
    }

    stopCluster(cl)
    if (verbose) cat("✓ Parallel processing complete\n")

  } else {
    subsample_results_list <- vector("list", n_subsamples)
    for (s in seq_len(n_subsamples)) {
      subsample_results_list[[s]] <- process_subsample(s)
      if (verbose) cat("✓ Subsample", s, "complete\n")
    }
  }

  # Reorganize into method-centric list
  for (s in seq_len(n_subsamples)) {
    sub_res <- subsample_results_list[[s]]
    for (method_name in names(sub_res)) {
      if (is.null(results[[method_name]])) {
        results[[method_name]] <- list()
      }
      results[[method_name]][[paste0("subsample", s)]] <- sub_res[[method_name]]
    }
  }

  # Aggregate using existing helper (requires 09-proper-nested-cv.R)
  aggregated <- aggregate_proper_nested_cv(
    results = results,
    clustering_methods = c(clustering_label),
    prediction_methods = prediction_methods,
    n_subsamples = n_subsamples,
    include_baseline = TRUE
  )

  list(
    raw_results = results,
    aggregated = aggregated
  )
}


#' Proper nested CV with subsampling (two-stage clustering + Leiden on centroids)
#'
#' Wrapper around proper_nested_cv_two_stage_leiden() to:
#'   - subsample cells (optionally patient-stratified)
#'   - generate patient-level outer folds
#'   - aggregate results across subsamples
#'
#' Note: For meta-clustering we tune the Leiden resolution parameter.
proper_nested_cv_two_stage_leiden_with_subsampling <- function(
    full_data, full_p_id,
    outcome_obj, meta_obj = NULL,
    outcome_col_1, outcome_col_2,
    prediction_methods = c("elasticnet"),
    resolution_grid = seq(0.5, 4.0, by = 0.5),
    n_cells = 1000000,
    n_subsamples = 3,
    n_outer_folds = 3,
    micro_method = "kmeans",
    n_microclusters = 5000,
    n_neighbors = 30,
    lambda_grid = 10^seq(-3, 1, length.out = 20),
    alpha = 0.5,
    scale_features = TRUE,
    patient_stratified_sampling = TRUE,
    outcome_type = c("fold_change", "absolute"),
    include_dose = TRUE,
    include_age = TRUE,
    include_sex = TRUE,
    include_time = FALSE,
    n_cores = 1,
    verbose = FALSE) {

  outcome_type <- match.arg(outcome_type)

  # Prepare outcomes (+ covariates if provided)
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
      patient_outcomes <- prepare_outcome_absolute(
        outcome_obj, outcome_col_1, outcome_col_2, full_p_id
      )
    } else {
      patient_outcomes <- prepare_outcome(
        outcome_obj, outcome_col_1, outcome_col_2, full_p_id
      )
    }
  }

  # Label for aggregation
  clustering_label <- paste0(
    "two_stage_leiden_", micro_method,
    "_micro", n_microclusters,
    "_nn", n_neighbors
  )

  results <- list()

  process_subsample <- function(subsample_idx) {
    if (verbose && n_cores == 1) {
      cat("\n", rep("-", 70), "\n", sep = "")
      cat("SUBSAMPLE ", subsample_idx, " / ", n_subsamples, "\n", sep = "")
      cat(rep("-", 70), "\n", sep = "")
    }

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

    # Outer folds on patients with outcomes (consistent with other pipelines)
    set.seed(2025 + subsample_idx)
    unique_p_ids <- names(patient_outcomes$y)
    outer_folds <- createFolds(seq_along(unique_p_ids),
                              k = n_outer_folds,
                              list = TRUE)

    subsample_results <- list()

    for (pred_method in prediction_methods) {
      method_name <- paste(clustering_label, pred_method, sep = "_")

      cat(sprintf(
        "[Two-stage+Leiden] Subsample %d/%d, Method: %s\n",
        subsample_idx, n_subsamples, pred_method
      ))

      cv_res <- proper_nested_cv_two_stage_leiden(
        mini_data = mini_data,
        p_id = p_id,
        patient_outcomes = patient_outcomes,
        resolution_grid = resolution_grid,
        outer_folds = outer_folds,
        prediction_method = pred_method,
        micro_method = micro_method,
        n_microclusters = n_microclusters,
        n_neighbors = n_neighbors,
        lambda_grid = lambda_grid,
        mtry_grid = NULL,
        nodesize_grid = c(1, 5, 10),
        alpha = alpha,
        scale_features = scale_features,
        verbose = verbose
      )

      subsample_results[[method_name]] <- cv_res
    }

    for (pred_method in prediction_methods) {
      method_name <- paste0("baseline_only_", pred_method)
      base_res <- baseline_only_cv(
        patient_outcomes = patient_outcomes,
        outer_folds = outer_folds,
        prediction_method = pred_method,
        lambda_grid = lambda_grid,
        mtry_grid = NULL,
        nodesize_grid = c(1, 5, 10),
        alpha = alpha,
        verbose = FALSE
      )
      subsample_results[[method_name]] <- base_res
    }

    subsample_results
  }

  # Run subsamples (optionally in parallel)
  if (n_cores > 1) {
    if (verbose) {
      cat("Running", n_subsamples, "subsamples in parallel using",
          n_cores, "cores...\n")
    }

    library(parallel)
    library(doParallel)
    library(foreach)

    cl <- makeCluster(n_cores)
    registerDoParallel(cl)

    clusterEvalQ(cl, {
      library(here)
      source(here::here("methods", "00-paths.R"))
      source(here::here("methods", "02-data-preparation.R"))
      source(here::here("methods", "03-clustering-methods.R"))
      source(here::here("methods", "04-prediction-models.R"))
      source(here::here("methods", "07-nested-cv-pipeline.R"))
      source(here::here("methods", "08-multi-outcome-covariates.R"))
      source(here::here("methods", "09-proper-nested-cv.R"))
      source(here::here("methods", "10-two-stage-clustering.R"))
    })

    clusterExport(
      cl,
      c("full_data", "full_p_id", "patient_outcomes", "n_cells",
        "patient_stratified_sampling", "n_outer_folds",
        "resolution_grid", "lambda_grid", "alpha", "scale_features",
        "prediction_methods", "micro_method", "n_microclusters",
        "n_neighbors",
        "baseline_only_cv", "proper_nested_cv_two_stage_leiden",
        "two_stage_microclustering", "two_stage_build_centroid_knn_graph",
        "two_stage_leiden_micro_to_meta", "two_stage_assign_test_micro",
        "two_stage_patient_micro_counts", "two_stage_collapse_micro_to_meta_counts",
        "tune_resolution_on_microstructure",
        "clustering_label"),
      envir = environment()
    )

    subsample_results_list <- foreach(subsample_idx = 1:n_subsamples,
                                      .packages = c("glmnet", "randomForest", "caret")) %dopar% {
      process_subsample(subsample_idx)
    }

    stopCluster(cl)
    if (verbose) cat("✓ Parallel processing complete\n")

  } else {
    subsample_results_list <- vector("list", n_subsamples)
    for (s in seq_len(n_subsamples)) {
      subsample_results_list[[s]] <- process_subsample(s)
      if (verbose) cat("✓ Subsample", s, "complete\n")
    }
  }

  # Reorganize into method-centric list
  for (s in seq_len(n_subsamples)) {
    sub_res <- subsample_results_list[[s]]
    for (method_name in names(sub_res)) {
      if (is.null(results[[method_name]])) {
        results[[method_name]] <- list()
      }
      results[[method_name]][[paste0("subsample", s)]] <- sub_res[[method_name]]
    }
  }

  aggregated <- aggregate_proper_nested_cv(
    results = results,
    clustering_methods = c(clustering_label),
    prediction_methods = prediction_methods,
    n_subsamples = n_subsamples,
    include_baseline = TRUE
  )

  list(
    raw_results = results,
    aggregated = aggregated
  )
}
