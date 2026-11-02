#!/usr/bin/env Rscript
# =============================================================================
# E3-interaction-baseline.R
#
# Post-hoc analysis addressing Andrew's review comment that a linear-only
# clinical baseline (M2) may underestimate clinical predictive power. We fit
# an interaction-rich clinical baseline M2_prime (all pairwise interactions of
# the 7 clinical features) and compare its 5-fold test MSE to the existing M3
# (titers + clusters). If clusters still beat the interaction-rich baseline,
# the cluster contribution cannot be explained by un-modeled pairwise
# nonlinearity in the clinical block.
#
# Data: 5 outer-fold artifacts at
#   /fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06/per_fold_atlas_5000/fold_{i}/fold_{i}_all.RData
# Each holds X_M2_train, X_M2_test, Y_train, Y_test, train_p_ids, test_p_ids,
# and pred_M3_en_test (the existing M3 elastic-net predictions).
#
# Model: cv.glmnet(family = "mgaussian", alpha = 0.9), patient-grouped inner
# foldids (3 inner folds), lambda.min. Deterministic seed = 7000.
#
# Outputs (no existing artifacts modified):
#   /fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06/post_hoc_5000/interaction_baseline/
#     E3_interaction_M2prime.RData
#     E3_interaction_M2prime.csv
# =============================================================================

FLU_ROOT <- Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06")
suppressPackageStartupMessages({
  library(glmnet)
})

# -----------------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------------
SEED        <- 7000L
ALPHA       <- 0.9
N_INNER     <- 3L          # inner CV folds for lambda selection on training
N_CORES     <- 5L          # parallel across outer folds (5 folds total)
VIRUSES     <- c("h1n1", "h3n2", "bvic")

FOLD_ROOT   <- file.path(FLU_ROOT, "per_fold_atlas_5000")
OUT_DIR     <- file.path(FLU_ROOT, "post_hoc_5000/interaction_baseline")
OUT_RDATA   <- file.path(OUT_DIR, "E3_interaction_M2prime.RData")
OUT_CSV     <- file.path(OUT_DIR, "E3_interaction_M2prime.csv")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(fmt, ...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")),
      sprintf(fmt, ...), "\n", sep = "")
}

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

# Reuses the project's convention: p_id_yr -> base patient id (strip _yrN).
base_pid_from_p_id_yr <- function(p_id_yr) {
  sub("_yr[0-9]+$", "", as.character(p_id_yr))
}

# Patient-grouped fold ids for cv.glmnet lambda selection.
compute_grouped_foldid <- function(patient_ids, n_folds = 3L, seed = NULL) {
  base_pids <- base_pid_from_p_id_yr(patient_ids)
  unique_bases <- unique(base_pids)
  if (length(unique_bases) < n_folds) n_folds <- length(unique_bases)
  if (!is.null(seed)) set.seed(seed)
  group_folds <- sample(rep(seq_len(n_folds), length.out = length(unique_bases)))
  names(group_folds) <- unique_bases
  as.integer(group_folds[base_pids])
}

# Build pairwise-interaction design matrix from a 7-column clinical matrix.
# Uses model.matrix(~ .^2 - 1) on a data.frame. Returns a numeric matrix.
build_interaction_mat <- function(X) {
  df <- as.data.frame(X, stringsAsFactors = FALSE)
  M  <- stats::model.matrix(~ .^2 - 1, data = df)
  # Strip the "(Intercept)" attribute carry-over so glmnet sees a plain matrix.
  attr(M, "assign")    <- NULL
  attr(M, "contrasts") <- NULL
  M
}

# Drop constant (zero-variance) columns based on a training matrix; apply the
# same column subset to test. Returns a list with aligned train/test matrices.
drop_constant_columns <- function(X_train, X_test) {
  sds <- apply(X_train, 2, stats::sd)
  keep <- which(is.finite(sds) & sds > 0)
  if (length(keep) == 0L) stop("All interaction columns are constant on training set.")
  list(train = X_train[, keep, drop = FALSE],
       test  = X_test[, keep, drop = FALSE],
       kept  = colnames(X_train)[keep])
}

# Fit a single outer fold. Returns per-fold predictions and MSE by virus.
fit_one_fold <- function(fold_id) {
  fold_path <- file.path(FOLD_ROOT, sprintf("fold_%d", fold_id),
                        sprintf("fold_%d_all.RData", fold_id))
  if (!file.exists(fold_path)) {
    stop(sprintf("Missing fold artifact: %s", fold_path))
  }
  e <- new.env(parent = emptyenv())
  load(fold_path, envir = e)

  X_M2_train      <- e$X_M2_train
  X_M2_test       <- e$X_M2_test
  Y_train         <- e$Y_train
  Y_test          <- e$Y_test
  train_p_ids     <- e$train_p_ids
  test_p_ids      <- e$test_p_ids
  pred_M3_en_test <- e$pred_M3_en_test

  stopifnot(identical(colnames(Y_train), VIRUSES))
  stopifnot(identical(colnames(Y_test),  VIRUSES))
  stopifnot(nrow(X_M2_train) == nrow(Y_train))
  stopifnot(nrow(X_M2_test)  == nrow(Y_test))

  # Build M2_prime = clinical + all pairwise interactions
  Xp_train <- build_interaction_mat(X_M2_train)
  Xp_test  <- build_interaction_mat(X_M2_test)

  # Align test columns to training columns (defensive: model.matrix should yield
  # identical column orders here, but guard against ordering drift).
  common_cols <- intersect(colnames(Xp_train), colnames(Xp_test))
  Xp_train <- Xp_train[, common_cols, drop = FALSE]
  Xp_test  <- Xp_test[,  common_cols, drop = FALSE]

  # Drop training-constant columns and mirror on test
  cleaned <- drop_constant_columns(Xp_train, Xp_test)
  Xp_train <- cleaned$train
  Xp_test  <- cleaned$test
  kept_features <- cleaned$kept

  # Patient-grouped inner foldids for lambda CV (3 folds)
  inner_foldid <- compute_grouped_foldid(
    train_p_ids,
    n_folds = N_INNER,
    seed    = SEED + fold_id * 13L
  )

  # Fit multi-task elastic net
  set.seed(SEED + fold_id * 17L)
  cvfit <- cv.glmnet(
    x        = Xp_train,
    y        = Y_train,
    family   = "mgaussian",
    alpha    = ALPHA,
    foldid   = inner_foldid,
    standardize = TRUE,
    parallel = FALSE
  )

  lambda_used <- cvfit$lambda.min

  pred_array <- predict(cvfit, newx = Xp_test, s = lambda_used)
  # mgaussian returns array [n, n_virus, 1]; reduce to matrix [n, n_virus]
  pred_test <- as.matrix(pred_array[, , 1L])
  colnames(pred_test) <- VIRUSES
  rownames(pred_test) <- rownames(Xp_test)

  # Per-virus test MSE for M2_prime and M3
  mse_M2p <- vapply(VIRUSES, function(v) {
    mean((pred_test[, v] - Y_test[, v])^2)
  }, numeric(1))
  mse_M3 <- vapply(VIRUSES, function(v) {
    mean((pred_M3_en_test[, v] - Y_test[, v])^2)
  }, numeric(1))

  list(
    fold_id        = fold_id,
    n_train        = nrow(Xp_train),
    n_test         = nrow(Xp_test),
    n_features     = ncol(Xp_train),
    lambda_used    = lambda_used,
    kept_features  = kept_features,
    pred_M2p_test  = pred_test,
    pred_M3_test   = pred_M3_en_test,
    Y_test         = Y_test,
    test_p_ids     = test_p_ids,
    mse_M2p_virus  = mse_M2p,
    mse_M3_virus   = mse_M3
  )
}

# -----------------------------------------------------------------------------
# Run across folds (parallel via mclapply when available)
# -----------------------------------------------------------------------------
log_msg("E3 interaction baseline starting (seed=%d, alpha=%.2f, inner_folds=%d).",
        SEED, ALPHA, N_INNER)

use_parallel <- requireNamespace("parallel", quietly = TRUE) &&
                .Platform$OS.type == "unix"

fold_ids <- 1:5

if (use_parallel) {
  log_msg("Running 5 outer folds in parallel via mclapply (%d cores).", N_CORES)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(SEED)
  fold_results <- parallel::mclapply(
    fold_ids,
    fit_one_fold,
    mc.cores = min(N_CORES, length(fold_ids)),
    mc.preschedule = FALSE
  )
} else {
  log_msg("Running 5 outer folds sequentially.")
  fold_results <- lapply(fold_ids, fit_one_fold)
}

# Sanity-check no fold errored out
bad <- vapply(fold_results, function(r) inherits(r, "try-error") || is.null(r), logical(1))
if (any(bad)) {
  msgs <- vapply(fold_results[bad], function(x) as.character(x), character(1))
  stop("Some folds failed:\n", paste(msgs, collapse = "\n"))
}

for (r in fold_results) {
  log_msg("Fold %d: n_train=%d n_test=%d n_features=%d lambda=%.4g",
          r$fold_id, r$n_train, r$n_test, r$n_features, r$lambda_used)
  log_msg("  M2': mse h1n1=%.4f h3n2=%.4f bvic=%.4f",
          r$mse_M2p_virus[["h1n1"]], r$mse_M2p_virus[["h3n2"]], r$mse_M2p_virus[["bvic"]])
  log_msg("  M3 : mse h1n1=%.4f h3n2=%.4f bvic=%.4f",
          r$mse_M3_virus[["h1n1"]], r$mse_M3_virus[["h3n2"]], r$mse_M3_virus[["bvic"]])
}

# -----------------------------------------------------------------------------
# Aggregate global MSE (n_test-weighted across folds; equivalent to pooling
# squared errors across all test rows from all folds).
# -----------------------------------------------------------------------------
n_test_vec <- vapply(fold_results, `[[`, integer(1), "n_test")
total_n    <- sum(n_test_vec)

agg_virus <- function(field) {
  num <- 0
  for (r in fold_results) {
    num <- num + r$n_test * r[[field]]
  }
  num / total_n
}

mse_M2p_virus_global <- agg_virus("mse_M2p_virus")
mse_M3_virus_global  <- agg_virus("mse_M3_virus")

# Overall MSE = mean across viruses of the per-virus n-weighted MSE
mse_M2p_overall <- mean(mse_M2p_virus_global)
mse_M3_overall  <- mean(mse_M3_virus_global)

# Per-virus and overall improvement (positive = M3 beats M2_prime)
improvement_virus_pct <- 100 * (mse_M2p_virus_global - mse_M3_virus_global) /
                                mse_M2p_virus_global
improvement_overall_pct <- 100 * (mse_M2p_overall - mse_M3_overall) /
                                  mse_M2p_overall

log_msg("---- Global aggregate (n_test-weighted across folds) ----")
log_msg("M2' (interaction baseline) overall MSE = %.5f", mse_M2p_overall)
log_msg("M3  (titers + clusters)    overall MSE = %.5f", mse_M3_overall)
log_msg("Overall improvement of M3 over M2' = %.2f%%", improvement_overall_pct)
for (v in VIRUSES) {
  log_msg("  %-4s: M2'=%.4f  M3=%.4f  improvement=%.2f%%",
          v,
          mse_M2p_virus_global[[v]],
          mse_M3_virus_global[[v]],
          improvement_virus_pct[[v]])
}

# -----------------------------------------------------------------------------
# Persist results
# -----------------------------------------------------------------------------
summary_tbl <- data.frame(
  virus               = c(VIRUSES, "overall"),
  mse_M2_prime        = c(unname(mse_M2p_virus_global), mse_M2p_overall),
  mse_M3              = c(unname(mse_M3_virus_global),  mse_M3_overall),
  improvement_pct     = c(unname(improvement_virus_pct), improvement_overall_pct),
  stringsAsFactors    = FALSE
)

per_fold_tbl <- do.call(rbind, lapply(fold_results, function(r) {
  data.frame(
    fold_id      = rep(r$fold_id, length(VIRUSES)),
    n_train      = rep(r$n_train, length(VIRUSES)),
    n_test       = rep(r$n_test,  length(VIRUSES)),
    n_features   = rep(r$n_features, length(VIRUSES)),
    lambda_used  = rep(r$lambda_used, length(VIRUSES)),
    virus        = VIRUSES,
    mse_M2_prime = unname(r$mse_M2p_virus[VIRUSES]),
    mse_M3       = unname(r$mse_M3_virus[VIRUSES]),
    stringsAsFactors = FALSE
  )
}))

utils::write.csv(summary_tbl, OUT_CSV, row.names = FALSE)
log_msg("Wrote summary CSV: %s", OUT_CSV)

save(
  fold_results,
  summary_tbl,
  per_fold_tbl,
  mse_M2p_virus_global,
  mse_M3_virus_global,
  mse_M2p_overall,
  mse_M3_overall,
  improvement_virus_pct,
  improvement_overall_pct,
  SEED,
  ALPHA,
  N_INNER,
  VIRUSES,
  file = OUT_RDATA
)
log_msg("Wrote RData: %s", OUT_RDATA)

log_msg("Done.")
