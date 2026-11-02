#!/usr/bin/env Rscript
# E14-01: Per-fold atlas nested CV on FIRST vaccination episodes only.
#
# Inclusion: each recipient's first-season episode (visit_type "first" in the
# episode metadata; 262 episodes). Second-season episodes (repeaters) are not
# used anywhere in model building, including atlas construction: k-means
# micro-clusters and Leiden candidates are rebuilt within each outer training
# fold from first-episode cells only. Outer folds keep each recipient's fold
# from the primary all-episode analysis, so the two analyses are paired by fold.
#
# Repeater validation: each second-season episode is assigned to the fold model
# that held out its recipient, mapped to that fold's atlas, and scored with the
# fitted M1-M4 models (multitask and composite). Nothing about these episodes
# enters any fit.
#
# Everything else (resolution grid, inner CV rule, elastic-net settings, seeds)
# follows E2-perfold-4model.R. Run one fold per process with FOLD_ID=f, then
# aggregate with AGGREGATE_ONLY=1.
suppressPackageStartupMessages({
  library(here)
})

log_msg <- function(fmt, ...) {
  msg <- sprintf(fmt, ...)
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
  flush.console()
}

build_baseline_mat <- function(p_ids, outcomes, viruses) {
  mat <- vapply(
    viruses,
    function(v) as.numeric(outcomes[[v]]$x[p_ids]),
    numeric(length(p_ids))
  )
  dim(mat) <- c(length(p_ids), length(viruses))
  rownames(mat) <- p_ids
  colnames(mat) <- paste0("baseline_", viruses)
  mat
}

build_y_mat <- function(p_ids, outcomes, viruses) {
  mat <- vapply(
    viruses,
    function(v) as.numeric(outcomes[[v]]$y[p_ids]),
    numeric(length(p_ids))
  )
  dim(mat) <- c(length(p_ids), length(viruses))
  rownames(mat) <- p_ids
  colnames(mat) <- viruses
  mat
}

extract_mgaussian_coefs <- function(fit_obj, viruses, lambda_rule) {
  s_val <- if (lambda_rule == "lambda.min") fit_obj$model$lambda.min else fit_obj$model$lambda.1se
  coef_list <- coef(fit_obj$model, s = s_val)
  if (is.null(names(coef_list)) || !all(viruses %in% names(coef_list))) {
    names(coef_list) <- viruses[seq_along(coef_list)]
  }

  out <- setNames(vector("list", length(viruses)), viruses)
  for (v in viruses) {
    cf <- coef_list[[v]]
    vec <- as.numeric(cf)
    names(vec) <- rownames(cf)
    out[[v]] <- vec
  }
  out
}

# ============================================================================
# Parameters (must match prior experiments)
# ============================================================================
VIRUSES <- c("h1n1", "h3n2", "bvic")
N_MICRO <- as.integer(Sys.getenv("N_MICRO", "5000"))
N_NEIGHBORS <- 30L
KMEANS_FIT_CELLS <- 1000000L
N_OUTER_FOLDS <- 5L
N_INNER_FOLDS <- 3L
ALPHA <- 0.9
FOLD_SEED <- 7000L
SEED_BASE <- 7100L
N_CORES <- as.integer(Sys.getenv("N_CORES", "3"))
FEATURE_TRANSFORM <- "clr"
LAMBDA_RULE <- "lambda.min"
TRAIN_COUNTS_SOURCE <- "full"
RESOLUTION_GRID <- seq(0.5, 5.0, by = 0.25)
ASSIGNMENT_CHUNK_SIZE <- 300000L

OUTPUT_ROOT <- Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/first_episode_2026-10-14")
PRIMARY_ROOT <- "/fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06"
FOLD_ONLY <- as.integer(Sys.getenv("FOLD_ID", "0"))
AGGREGATE_ONLY <- identical(Sys.getenv("AGGREGATE_ONLY", "0"), "1")
CACHE_ROOT <- Sys.getenv("FLU_ATLAS_CACHE_ROOT", "/fast/datasets/flu_cellular_prediction")
PER_FOLD_ROOT <- file.path(OUTPUT_ROOT, paste0("per_fold_atlas_", N_MICRO))
AGGREGATED_DIR <- file.path(PER_FOLD_ROOT, "aggregated")

for (d in c(OUTPUT_ROOT, PER_FOLD_ROOT, AGGREGATED_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ============================================================================
# Load code and data
# ============================================================================
log_msg("Sourcing methods files...")
source(here("methods", "00-paths.R"))
source(here("methods", "08-multi-outcome-covariates.R"))
source(here("methods", "35-per-fold-atlas-cv.R"))

# Episode universe of the primary analysis (343 episodes) and its outer folds.
cache_agg <- new.env()
load(file.path(PRIMARY_ROOT, paste0("per_fold_atlas_", N_MICRO), "aggregated", "E2_aggregated.RData"), envir = cache_agg)
full_p_id <- cache_agg$common_p_ids
primary_fold_of <- setNames(integer(0), character(0))
for (f in seq_len(N_OUTER_FOLDS)) {
  pe <- new.env(); load(file.path(PRIMARY_ROOT, paste0("per_fold_atlas_", N_MICRO), sprintf("fold_%d", f),
                                  sprintf("fold_%d_all.RData", f)), envir = pe)
  primary_fold_of <- c(primary_fold_of, setNames(rep(f, length(pe$test_p_ids)), pe$test_p_ids))
  rm(pe)
}

if (!AGGREGATE_ONLY) {
  log_msg("Loading cell-level object (this may take a while)...")
  load(file.path(paths$raw, "hsct_transplant_object_harmony_final_small.RData"))
  full_data <- obj$batch_correction$data
  cell_p_id <- obj$source
  rm(obj)
  invisible(gc())
}

log_msg("Loading outcomes and metadata CSV files...")
outcome_obj <- read.csv(
  file.path(paths$raw, "transplant_cytof_counts_hai_12MAY2024.csv"),
  stringsAsFactors = FALSE
)
meta_obj <- read.csv(
  file.path(paths$raw, "hct_limited-dataset_2025-07-09.csv"),
  stringsAsFactors = FALSE
)

timing_source_md5 <- unname(tools::md5sum(Sys.getenv("FLU_EPISODE_METADATA",
  file.path(paths$raw, "hsct_metadata_15JUL2025.csv"))))

log_msg("Preparing per-virus outcomes and shared covariates...")
outcomes <- setNames(vector("list", length(VIRUSES)), VIRUSES)
for (v in VIRUSES) {
  baseline_col <- paste0("geo_mean_", v, "_v1")
  target_col <- paste0("geo_mean_", v, "_v3")
  outcomes[[v]] <- prepare_outcome_with_covariates(
    outcome_obj = outcome_obj,
    meta_obj = meta_obj,
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

covariates <- outcomes[[VIRUSES[1]]]$covariates[common_p_ids, , drop = FALSE]
covariates_all <- covariates
all_p_ids <- common_p_ids                     # 343 episodes (for repeater scoring)

episode_meta <- read.csv(Sys.getenv("FLU_EPISODE_METADATA", file.path(paths$raw, "hsct_metadata_15JUL2025.csv")),
                         stringsAsFactors = FALSE)
episode_meta$p_id_yr <- sub("_([0-9]+)$", "_yr\\1", episode_meta$id_yr)
visit_type <- setNames(episode_meta$visit_type, episode_meta$p_id_yr)[all_p_ids]
stopifnot(!anyNA(visit_type), all(visit_type %in% c("first", "repeat")))
common_p_ids <- all_p_ids[visit_type == "first"]
repeat_p_ids_all <- all_p_ids[visit_type == "repeat"]
log_msg("First episodes: %d; second-season episodes held out: %d", length(common_p_ids), length(repeat_p_ids_all))

fold_strata <- infer_cohort_from_p_id(common_p_ids)
if (length(unique(fold_strata)) <= 1L) {
  fold_strata <- NULL
}

group_ids <- setNames(base_pid_from_p_id_yr(all_p_ids), all_p_ids)

# Outer folds: each recipient keeps its fold from the primary analysis.
fold_of_recipient <- tapply(primary_fold_of[all_p_ids], group_ids[all_p_ids], function(x) unique(x))
stopifnot(all(lengths(fold_of_recipient) == 1L))
first_fold <- unlist(fold_of_recipient)[group_ids[common_p_ids]]
outer_folds <- lapply(seq_len(N_OUTER_FOLDS), function(f) which(first_fold == f))
repeat_fold <- unlist(fold_of_recipient)[group_ids[repeat_p_ids_all]]
names(repeat_fold) <- repeat_p_ids_all

log_msg("Outer fold test sizes: %s", paste(vapply(outer_folds, length, integer(1)), collapse = ", "))

# cv.glmnet parallel backend
if (N_CORES > 1L) {
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop("Package 'doParallel' is required when N_CORES > 1.")
  }
  log_msg("Registering doParallel backend with %d workers...", N_CORES)
  cl <- parallel::makeCluster(N_CORES)
  doParallel::registerDoParallel(cl)
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
  }, add = TRUE)
}

# ============================================================================
# Per-fold run (resumable)
# ============================================================================
for (fold_id in if (AGGREGATE_ONLY) integer(0) else if (FOLD_ONLY > 0) FOLD_ONLY else seq_len(N_OUTER_FOLDS)) {
  fold_dir <- file.path(PER_FOLD_ROOT, sprintf("fold_%d", fold_id))
  dir.create(fold_dir, recursive = TRUE, showWarnings = FALSE)
  fold_out_path <- file.path(fold_dir, sprintf("fold_%d_all.RData", fold_id))

  if (file.exists(fold_out_path)) {
    log_msg("Fold %d already exists at %s. Skipping (resumable).", fold_id, fold_out_path)
    next
  }

  log_msg("================ Fold %d/%d ================", fold_id, N_OUTER_FOLDS)

  test_idx <- outer_folds[[fold_id]]
  test_p_ids <- common_p_ids[test_idx]
  train_p_ids <- setdiff(common_p_ids, test_p_ids)
  log_msg("Fold %d: train=%d, test=%d", fold_id, length(train_p_ids), length(test_p_ids))

  repeat_p_ids <- repeat_p_ids_all[repeat_fold == fold_id]
  log_msg("Fold %d: second-season episodes to score: %d", fold_id, length(repeat_p_ids))

  # 1) Build atlas on training first-episode cells only
  log_msg("Fold %d: building atlas...", fold_id)
  atlas <- build_atlas_for_fold(
    full_data = full_data,
    full_p_id = cell_p_id,
    train_patient_ids = train_p_ids,
    n_micro = N_MICRO,
    n_neighbors = N_NEIGHBORS,
    kmeans_fit_n_cells = KMEANS_FIT_CELLS,
    seed = SEED_BASE + fold_id * 1000,
    verbose = TRUE
  )
  atlas_centroids <- atlas$centroids
  atlas_scale_center <- atlas$scale_center
  atlas_scale_scale <- atlas$scale_scale
  atlas_centroid_graph <- atlas$centroid_graph
  atlas_timing_s <- atlas$timing_s
  rm(atlas); invisible(gc())

  # 2) Micro-counts for train, test and held-out second-season episodes
  micro_counts_for <- function(ids, seed_offset) compute_patient_micro_counts(
    full_data = full_data, full_p_id = cell_p_id, patient_ids = ids,
    centroids = atlas_centroids, scale_center = atlas_scale_center, scale_scale = atlas_scale_scale,
    chunk_size = ASSIGNMENT_CHUNK_SIZE, seed = SEED_BASE + fold_id * 1000 + seed_offset, verbose = FALSE)
  log_msg("Fold %d: computing train micro-counts...", fold_id)
  train_micro_counts <- micro_counts_for(train_p_ids, 11)
  log_msg("Fold %d: computing test micro-counts...", fold_id)
  test_micro_counts <- micro_counts_for(test_p_ids, 17)
  log_msg("Fold %d: computing second-season micro-counts...", fold_id)
  repeat_micro_counts <- micro_counts_for(repeat_p_ids, 23)

  # 3) Inner CV resolution selection (score by M4 vs OLS baseline ratio)
  log_msg("Fold %d: inner CV resolution selection...", fold_id)
  inner_folds <- create_group_folds(
    p_ids = train_p_ids,
    group_ids = group_ids,
    k = N_INNER_FOLDS,
    seed = FOLD_SEED + fold_id * 100,
    strata = if (is.null(fold_strata)) NULL else fold_strata[train_p_ids]
  )

  leiden_mappings <- lapply(RESOLUTION_GRID, function(res) {
    two_stage_leiden_micro_to_meta(atlas_centroid_graph, resolution = res)
  })
  saveRDS(leiden_mappings, file.path(fold_dir, "resolution_candidates.rds"))
  k_found_vec <- vapply(leiden_mappings, function(m) as.integer(m$k_found), integer(1))

  inner_base_mse <- matrix(
    NA_real_,
    nrow = length(inner_folds),
    ncol = length(VIRUSES),
    dimnames = list(NULL, VIRUSES)
  )

  for (ii in seq_along(inner_folds)) {
    val_idx <- inner_folds[[ii]]
    val_p <- train_p_ids[val_idx]
    tr_p <- setdiff(train_p_ids, val_p)

    base_train_mat <- build_baseline_mat(tr_p, outcomes, VIRUSES)
    base_val_mat <- build_baseline_mat(val_p, outcomes, VIRUSES)
    cov_tr <- covariates[tr_p, , drop = FALSE]
    cov_val <- covariates[val_p, , drop = FALSE]

    X_base_tr <- as.matrix(cbind(base_train_mat, cov_tr))
    X_base_val <- as.matrix(cbind(base_val_mat, cov_val))
    Y_tr <- build_y_mat(tr_p, outcomes, VIRUSES)
    Y_val <- build_y_mat(val_p, outcomes, VIRUSES)

    base_fit <- fit_ols_baseline_multi(X_base_tr, Y_tr, X_base_val)
    for (v in VIRUSES) {
      inner_base_mse[ii, v] <- mean((base_fit$pred_test[, v] - Y_val[, v])^2)
    }
  }

  resolution_scores <- numeric(length(RESOLUTION_GRID))
  names(resolution_scores) <- as.character(RESOLUTION_GRID)

  for (r_idx in seq_along(RESOLUTION_GRID)) {
    micro_to_meta_candidate <- leiden_mappings[[r_idx]]$micro_to_meta
    train_meta_counts_candidate <- two_stage_collapse_micro_to_meta_counts(
      micro_counts = train_micro_counts,
      micro_to_meta = micro_to_meta_candidate
    )
    cluster_train_all <- transform_counts(train_meta_counts_candidate, method = FEATURE_TRANSFORM)

    fold_scores <- numeric(length(inner_folds))

    for (ii in seq_along(inner_folds)) {
      val_idx <- inner_folds[[ii]]
      val_p <- train_p_ids[val_idx]
      tr_p <- setdiff(train_p_ids, val_p)

      base_train_mat <- build_baseline_mat(tr_p, outcomes, VIRUSES)
      base_val_mat <- build_baseline_mat(val_p, outcomes, VIRUSES)
      cov_tr <- covariates[tr_p, , drop = FALSE]
      cov_val <- covariates[val_p, , drop = FALSE]

      X_base_tr <- as.matrix(cbind(base_train_mat, cov_tr))
      X_base_val <- as.matrix(cbind(base_val_mat, cov_val))

      X_clust_tr <- as.matrix(cbind(cluster_train_all[tr_p, , drop = FALSE], X_base_tr))
      X_clust_val <- as.matrix(cbind(cluster_train_all[val_p, , drop = FALSE], X_base_val))

      Y_tr <- build_y_mat(tr_p, outcomes, VIRUSES)
      Y_val <- build_y_mat(val_p, outcomes, VIRUSES)

      inner_foldid <- compute_grouped_foldid(
        tr_p,
        n_folds = 5L,
        seed = SEED_BASE + fold_id * 10000 + ii
      )

      seed_fit <- SEED_BASE + fold_id * 10000 + r_idx * 100 + ii * 10
      clust_fit <- fit_predict_multitask_en(
        X_train = X_clust_tr,
        Y_train = Y_tr,
        X_test = X_clust_val,
        seed = seed_fit + 1,
        alpha = ALPHA,
        n_cores = N_CORES,
        foldid = inner_foldid,
        lambda_rule = LAMBDA_RULE
      )

      ratios <- vapply(VIRUSES, function(v) {
        mse_clust <- mean((clust_fit$pred_test[, v] - Y_val[, v])^2)
        denom <- inner_base_mse[ii, v]
        if (!is.finite(denom) || denom <= 0) {
          return(NA_real_)
        }
        mse_clust / denom
      }, numeric(1))

      fold_scores[ii] <- median(ratios, na.rm = TRUE)
    }

    resolution_scores[r_idx] <- mean(fold_scores, na.rm = TRUE)
    if (!is.finite(resolution_scores[r_idx])) {
      resolution_scores[r_idx] <- Inf
    }

    log_msg(
      "Fold %d: res=%.2f k=%d inner_score=%.4f",
      fold_id,
      RESOLUTION_GRID[r_idx],
      k_found_vec[r_idx],
      resolution_scores[r_idx]
    )
  }

  best_r_idx <- which.min(resolution_scores)[1]
  selected_resolution <- RESOLUTION_GRID[best_r_idx]
  k_found <- as.integer(k_found_vec[best_r_idx])
  k_found_by_res <- setNames(as.integer(k_found_vec), as.character(RESOLUTION_GRID))
  micro_to_meta <- leiden_mappings[[best_r_idx]]$micro_to_meta

  log_msg("Fold %d: selected resolution=%.2f (k=%d)", fold_id, selected_resolution, k_found)

  # 4) Final features at selected resolution
  train_meta_counts <- two_stage_collapse_micro_to_meta_counts(train_micro_counts, micro_to_meta)
  test_meta_counts <- two_stage_collapse_micro_to_meta_counts(test_micro_counts, micro_to_meta)

  clr_train <- transform_counts(train_meta_counts, method = FEATURE_TRANSFORM)
  clr_test <- transform_counts(test_meta_counts, method = FEATURE_TRANSFORM)
  repeat_meta_counts <- two_stage_collapse_micro_to_meta_counts(repeat_micro_counts, micro_to_meta)
  clr_repeat <- transform_counts(repeat_meta_counts, method = FEATURE_TRANSFORM)

  base_train_mat <- build_baseline_mat(train_p_ids, outcomes, VIRUSES)
  base_test_mat <- build_baseline_mat(test_p_ids, outcomes, VIRUSES)
  cov_train <- covariates[train_p_ids, , drop = FALSE]
  cov_test <- covariates[test_p_ids, , drop = FALSE]

  X_M1_train <- as.matrix(base_train_mat)
  X_M1_test <- as.matrix(base_test_mat)

  X_M2_train <- as.matrix(cbind(base_train_mat, cov_train))
  X_M2_test <- as.matrix(cbind(base_test_mat, cov_test))

  X_M3_train <- as.matrix(cbind(clr_train[train_p_ids, , drop = FALSE], base_train_mat))
  X_M3_test <- as.matrix(cbind(clr_test[test_p_ids, , drop = FALSE], base_test_mat))

  X_M4_train <- as.matrix(cbind(clr_train[train_p_ids, , drop = FALSE], base_train_mat, cov_train))
  X_M4_test <- as.matrix(cbind(clr_test[test_p_ids, , drop = FALSE], base_test_mat, cov_test))

  Y_train <- build_y_mat(train_p_ids, outcomes, VIRUSES)
  Y_test <- build_y_mat(test_p_ids, outcomes, VIRUSES)

  base_repeat_mat <- build_baseline_mat(repeat_p_ids, outcomes, VIRUSES)
  cov_repeat <- covariates_all[repeat_p_ids, , drop = FALSE]
  X_M1_repeat <- as.matrix(base_repeat_mat)
  X_M2_repeat <- as.matrix(cbind(base_repeat_mat, cov_repeat))
  X_M3_repeat <- as.matrix(cbind(clr_repeat[repeat_p_ids, , drop = FALSE], base_repeat_mat))
  X_M4_repeat <- as.matrix(cbind(clr_repeat[repeat_p_ids, , drop = FALSE], base_repeat_mat, cov_repeat))
  Y_repeat <- build_y_mat(repeat_p_ids, outcomes, VIRUSES)

  # 5) Fit all models
  outer_foldid <- compute_grouped_foldid(
    train_p_ids,
    n_folds = 5L,
    seed = SEED_BASE + fold_id * 1000 + 99
  )

  seed_fit_base <- SEED_BASE + fold_id * 1000 + best_r_idx * 50
  log_msg("Fold %d: fitting EN ladder models...", fold_id)

  fit_M1_en <- fit_predict_multitask_en(
    X_train = X_M1_train,
    Y_train = Y_train,
    X_test = X_M1_test,
    seed = seed_fit_base + 1,
    alpha = ALPHA,
    n_cores = N_CORES,
    foldid = outer_foldid,
    lambda_rule = LAMBDA_RULE
  )

  fit_M2_en <- fit_predict_multitask_en(
    X_train = X_M2_train,
    Y_train = Y_train,
    X_test = X_M2_test,
    seed = seed_fit_base + 2,
    alpha = ALPHA,
    n_cores = N_CORES,
    foldid = outer_foldid,
    lambda_rule = LAMBDA_RULE
  )

  fit_M3_en <- fit_predict_multitask_en(
    X_train = X_M3_train,
    Y_train = Y_train,
    X_test = X_M3_test,
    seed = seed_fit_base + 3,
    alpha = ALPHA,
    n_cores = N_CORES,
    foldid = outer_foldid,
    lambda_rule = LAMBDA_RULE
  )

  fit_M4_en <- fit_predict_multitask_en(
    X_train = X_M4_train,
    Y_train = Y_train,
    X_test = X_M4_test,
    seed = seed_fit_base + 4,
    alpha = ALPHA,
    n_cores = N_CORES,
    foldid = outer_foldid,
    lambda_rule = LAMBDA_RULE
  )

  log_msg("Fold %d: fitting OLS sensitivity models...", fold_id)
  fit_M1_ols <- fit_ols_baseline_multi(X_M1_train, Y_train, X_M1_test)
  fit_M2_ols <- fit_ols_baseline_multi(X_M2_train, Y_train, X_M2_test)

  pred_M1_en_train <- fit_M1_en$pred_train
  pred_M1_en_test <- fit_M1_en$pred_test
  pred_M2_en_train <- fit_M2_en$pred_train
  pred_M2_en_test <- fit_M2_en$pred_test
  pred_M3_en_train <- fit_M3_en$pred_train
  pred_M3_en_test <- fit_M3_en$pred_test
  pred_M4_en_train <- fit_M4_en$pred_train
  pred_M4_en_test <- fit_M4_en$pred_test

  pred_M1_ols_train <- fit_M1_ols$pred_train
  pred_M1_ols_test <- fit_M1_ols$pred_test
  pred_M2_ols_train <- fit_M2_ols$pred_train
  pred_M2_ols_test <- fit_M2_ols$pred_test

  for (obj_name in c(
    "pred_M1_en_train", "pred_M1_en_test", "pred_M2_en_train", "pred_M2_en_test",
    "pred_M3_en_train", "pred_M3_en_test", "pred_M4_en_train", "pred_M4_en_test",
    "pred_M1_ols_train", "pred_M1_ols_test", "pred_M2_ols_train", "pred_M2_ols_test"
  )) {
    obj <- get(obj_name)
    if (grepl("_train$", obj_name)) {
      rownames(obj) <- train_p_ids
    } else {
      rownames(obj) <- test_p_ids
    }
    assign(obj_name, obj)
  }

  coef_M1_en <- extract_mgaussian_coefs(fit_M1_en, VIRUSES, LAMBDA_RULE)
  coef_M2_en <- extract_mgaussian_coefs(fit_M2_en, VIRUSES, LAMBDA_RULE)
  coef_M3_en <- extract_mgaussian_coefs(fit_M3_en, VIRUSES, LAMBDA_RULE)
  coef_M4_en <- extract_mgaussian_coefs(fit_M4_en, VIRUSES, LAMBDA_RULE)

  lambda_M1 <- if (LAMBDA_RULE == "lambda.min") fit_M1_en$model$lambda.min else fit_M1_en$model$lambda.1se
  lambda_M2 <- if (LAMBDA_RULE == "lambda.min") fit_M2_en$model$lambda.min else fit_M2_en$model$lambda.1se
  lambda_M3 <- if (LAMBDA_RULE == "lambda.min") fit_M3_en$model$lambda.min else fit_M3_en$model$lambda.1se
  lambda_M4 <- if (LAMBDA_RULE == "lambda.min") fit_M4_en$model$lambda.min else fit_M4_en$model$lambda.1se

  # Second-season episodes scored by this fold's models (never used in fitting)
  pred_rep <- function(fit, X) {
    s <- if (LAMBDA_RULE == "lambda.min") fit$model$lambda.min else fit$model$lambda.1se
    m <- as_pred_matrix(predict(fit$model, newx = X, s = s), VIRUSES); rownames(m) <- rownames(X); m }
  pred_M1_en_repeat <- pred_rep(fit_M1_en, X_M1_repeat)
  pred_M2_en_repeat <- pred_rep(fit_M2_en, X_M2_repeat)
  pred_M3_en_repeat <- pred_rep(fit_M3_en, X_M3_repeat)
  pred_M4_en_repeat <- pred_rep(fit_M4_en, X_M4_repeat)

  # Composite outcome (mean of the three antigen log fold-rises): single-outcome
  # elastic nets with the same predictors, penalty folds and selected resolution.
  comp_fit <- function(Xtr, Xte, Xrep, seed) {
    set.seed(seed)
    fit <- glmnet::cv.glmnet(Xtr, rowMeans(Y_train), family = "gaussian", alpha = ALPHA,
                             foldid = outer_foldid, standardize = TRUE)
    list(test = as.vector(predict(fit, newx = Xte, s = LAMBDA_RULE)),
         repeat_ = as.vector(predict(fit, newx = Xrep, s = LAMBDA_RULE)),
         coef = as.matrix(coef(fit, s = LAMBDA_RULE)), lambda = fit[[LAMBDA_RULE]])
  }
  comp_M1 <- comp_fit(X_M1_train, X_M1_test, X_M1_repeat, seed_fit_base + 11)
  comp_M2 <- comp_fit(X_M2_train, X_M2_test, X_M2_repeat, seed_fit_base + 12)
  comp_M3 <- comp_fit(X_M3_train, X_M3_test, X_M3_repeat, seed_fit_base + 13)
  comp_M4 <- comp_fit(X_M4_train, X_M4_test, X_M4_repeat, seed_fit_base + 14)

  # Ensure required names and types for resolution objects
  resolution_scores <- setNames(as.numeric(resolution_scores), as.character(RESOLUTION_GRID))

  log_msg("Fold %d: saving full artifact bundle to %s", fold_id, fold_out_path)
  save(
    timing_source_md5,
    # Fold metadata
    fold_id,
    train_p_ids,
    test_p_ids,

    # Atlas artifacts
    atlas_centroids,
    atlas_scale_center,
    atlas_scale_scale,
    atlas_centroid_graph,
    atlas_timing_s,

    # Resolution selection
    selected_resolution,
    k_found,
    resolution_scores,
    k_found_by_res,
    micro_to_meta,

    # Count matrices
    train_micro_counts,
    test_micro_counts,
    train_meta_counts,
    test_meta_counts,

    # CLR features
    clr_train,
    clr_test,

    # Feature matrices (all models)
    X_M1_train, X_M1_test,
    X_M2_train, X_M2_test,
    X_M3_train, X_M3_test,
    X_M4_train, X_M4_test,

    # Y matrices
    Y_train,
    Y_test,

    # Model objects
    fit_M1_en,
    fit_M2_en,
    fit_M3_en,
    fit_M4_en,

    # Predictions
    pred_M1_en_train,
    pred_M1_en_test,
    pred_M2_en_train,
    pred_M2_en_test,
    pred_M3_en_train,
    pred_M3_en_test,
    pred_M4_en_train,
    pred_M4_en_test,
    pred_M1_ols_train,
    pred_M1_ols_test,
    pred_M2_ols_train,
    pred_M2_ols_test,

    # Coefficients
    coef_M1_en,
    coef_M2_en,
    coef_M3_en,
    coef_M4_en,

    # Lambda selected
    lambda_M1,
    lambda_M2,
    lambda_M3,
    lambda_M4,

    # Second-season episodes (held out) and composite models
    repeat_p_ids, repeat_micro_counts, repeat_meta_counts, clr_repeat,
    X_M1_repeat, X_M2_repeat, X_M3_repeat, X_M4_repeat, Y_repeat,
    pred_M1_en_repeat, pred_M2_en_repeat, pred_M3_en_repeat, pred_M4_en_repeat,
    comp_M1, comp_M2, comp_M3, comp_M4,

    file = fold_out_path
  )

  log_msg("Fold %d complete.", fold_id)
  invisible(gc())
}

# ============================================================================
# Aggregated out-of-fold summary
# ============================================================================
if (FOLD_ONLY > 0 && !AGGREGATE_ONLY) { log_msg("Fold %d done; aggregate separately.", FOLD_ONLY); quit(save = "no") }
log_msg("Aggregating out-of-fold predictions across folds...")

n_obs <- length(common_p_ids)
n_vir <- length(VIRUSES)

new_oof_mat <- function() {
  m <- matrix(NA_real_, nrow = n_obs, ncol = n_vir)
  rownames(m) <- common_p_ids
  colnames(m) <- VIRUSES
  m
}

oof_pred_M1_en <- new_oof_mat()
oof_pred_M2_en <- new_oof_mat()
oof_pred_M3_en <- new_oof_mat()
oof_pred_M4_en <- new_oof_mat()
oof_pred_M1_ols <- new_oof_mat()
oof_pred_M2_ols <- new_oof_mat()
oof_Y <- new_oof_mat()

selected_k_per_fold <- integer(N_OUTER_FOLDS)
selected_res_per_fold <- numeric(N_OUTER_FOLDS)
fold_assignments <- outer_folds

for (fold_id in seq_len(N_OUTER_FOLDS)) {
  fold_out_path <- file.path(PER_FOLD_ROOT, sprintf("fold_%d", fold_id), sprintf("fold_%d_all.RData", fold_id))
  if (!file.exists(fold_out_path)) {
    stop("Missing expected fold artifact: ", fold_out_path)
  }

  e <- new.env(parent = emptyenv())
  load(fold_out_path, envir = e)

  test_p_ids <- get("test_p_ids", envir = e)
  idx <- match(test_p_ids, common_p_ids)
  if (any(is.na(idx))) {
    stop("Failed to map fold ", fold_id, " test_p_ids to common_p_ids.")
  }

  Y_test <- get("Y_test", envir = e)
  if (is.null(rownames(Y_test))) {
    rownames(Y_test) <- test_p_ids
  }

  oof_pred_M1_en[idx, ] <- as.matrix(get("pred_M1_en_test", envir = e))
  oof_pred_M2_en[idx, ] <- as.matrix(get("pred_M2_en_test", envir = e))
  oof_pred_M3_en[idx, ] <- as.matrix(get("pred_M3_en_test", envir = e))
  oof_pred_M4_en[idx, ] <- as.matrix(get("pred_M4_en_test", envir = e))
  oof_pred_M1_ols[idx, ] <- as.matrix(get("pred_M1_ols_test", envir = e))
  oof_pred_M2_ols[idx, ] <- as.matrix(get("pred_M2_ols_test", envir = e))
  oof_Y[idx, ] <- as.matrix(Y_test)

  selected_k_per_fold[fold_id] <- as.integer(get("k_found", envir = e))
  selected_res_per_fold[fold_id] <- as.numeric(get("selected_resolution", envir = e))
}

if (anyNA(oof_Y)) {
  stop("oof_Y contains NA values after aggregation; fold outputs may be incomplete.")
}

model_preds <- list(
  M1_en = oof_pred_M1_en,
  M2_en = oof_pred_M2_en,
  M3_en = oof_pred_M3_en,
  M4_en = oof_pred_M4_en,
  M1_ols = oof_pred_M1_ols,
  M2_ols = oof_pred_M2_ols
)

mse_rows <- list()
for (v in VIRUSES) {
  for (m in names(model_preds)) {
    mse_rows[[length(mse_rows) + 1L]] <- data.frame(
      virus = v,
      model = m,
      mse = mean((model_preds[[m]][, v] - oof_Y[, v])^2),
      stringsAsFactors = FALSE
    )
  }
}
for (m in names(model_preds)) {
  mse_rows[[length(mse_rows) + 1L]] <- data.frame(
    virus = "overall",
    model = m,
    mse = mean((as.vector(model_preds[[m]]) - as.vector(oof_Y))^2),
    stringsAsFactors = FALSE
  )
}
mse_table <- do.call(rbind, mse_rows)

overall_mse <- setNames(mse_table$mse[mse_table$virus == "overall"], mse_table$model[mse_table$virus == "overall"])

pairwise_improvement_pct <- do.call(
  rbind,
  lapply(combn(names(model_preds), 2L, simplify = FALSE), function(pair) {
    m_a <- pair[1]
    m_b <- pair[2]
    mse_a <- overall_mse[[m_a]]
    mse_b <- overall_mse[[m_b]]
    data.frame(
      model_a = m_a,
      model_b = m_b,
      mse_a = mse_a,
      mse_b = mse_b,
      delta_mse = mse_a - mse_b,
      improvement_pct = 100 * (mse_a - mse_b) / mse_a,
      stringsAsFactors = FALSE
    )
  })
)

aggregated_rdata <- file.path(AGGREGATED_DIR, "E2_aggregated.RData")
aggregated_csv <- file.path(AGGREGATED_DIR, "E2_mse_summary.csv")

save(
  oof_pred_M1_en,
  oof_pred_M2_en,
  oof_pred_M3_en,
  oof_pred_M4_en,
  oof_pred_M1_ols,
  oof_pred_M2_ols,
  oof_Y,
  mse_table,
  overall_mse,
  pairwise_improvement_pct,
  fold_assignments,
  common_p_ids,
  selected_k_per_fold,
  selected_res_per_fold,
  file = aggregated_rdata
)

write.csv(mse_table, aggregated_csv, row.names = FALSE)

log_msg("Saved aggregated RData: %s", aggregated_rdata)
log_msg("Saved MSE summary CSV: %s", aggregated_csv)
log_msg("Done.")
