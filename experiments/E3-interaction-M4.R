# =============================================================================
# E3-interaction-M4: M4_int = clinical-with-interactions + cluster main effects
# Compares against M2_int (clinical-with-interactions only) to ask:
#   Do clusters add value beyond a generously-specified clinical baseline?
# Pairs with E3-interaction-baseline.R (M2_int vs M3 comparison).
# =============================================================================

FLU_ROOT <- Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06")
FLU_ANN_FILE <- Sys.getenv("FLU_ANN_FILE", "/home/amir/repo/flu-cellular-prediction/results/reports/cluster_annotations_corrected.csv")
suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(parallel)
})

PER_FOLD_ROOT <- file.path(FLU_ROOT, "per_fold_atlas_5000")
OUT_DIR <- file.path(FLU_ROOT, "post_hoc_5000/interaction_baseline")
ANN_FILE <- FLU_ANN_FILE

set.seed(7000)
n_folds <- 5
inner_folds <- 3
VIRUSES <- c("h1n1", "h3n2", "bvic")

build_grouped_foldid <- function(p_ids, k = 3, seed = 7000) {
  p_ids <- sub("_yr[0-9]+$", "", p_ids)
  set.seed(seed)
  uniq <- unique(p_ids)
  fold_assign <- sample(rep(seq_len(k), length.out = length(uniq)))
  names(fold_assign) <- uniq
  unname(fold_assign[p_ids])
}

fit_one_fold <- function(i) {
  f <- file.path(PER_FOLD_ROOT, sprintf("fold_%d", i), sprintf("fold_%d_all.RData", i))
  e <- new.env(parent = emptyenv())
  load(f, envir = e)

  # Build M2_int: clinical + pairwise interactions among M2 features
  Xm2_tr <- as.data.frame(e$X_M2_train)
  Xm2_te <- as.data.frame(e$X_M2_test)
  Xint_tr <- model.matrix(~ .^2 - 1, data = Xm2_tr)
  Xint_te <- model.matrix(~ .^2 - 1, data = Xm2_te)

  # Cluster features from M3 train/test (M3 = baseline_titers + cluster CLR features)
  # The cluster columns are in X_M3 minus the baseline titer columns
  Xm3_tr <- e$X_M3_train
  Xm3_te <- e$X_M3_test
  titer_cols <- grep("^baseline_", colnames(Xm3_tr), value = TRUE)
  cluster_cols <- setdiff(colnames(Xm3_tr), titer_cols)
  Xclu_tr <- Xm3_tr[, cluster_cols, drop = FALSE]
  Xclu_te <- Xm3_te[, cluster_cols, drop = FALSE]

  # M4_int = interaction-rich clinical + cluster main effects
  X_tr <- cbind(Xint_tr, Xclu_tr)
  X_te <- cbind(Xint_te, Xclu_te)

  Y_tr <- e$Y_train
  Y_te <- e$Y_test

  foldid <- build_grouped_foldid(e$train_p_ids, k = inner_folds, seed = 7000 + i * 13L)

  fit <- cv.glmnet(X_tr, Y_tr, family = "mgaussian", alpha = 0.9,
                   foldid = foldid, standardize = TRUE)
  lambda_min <- fit$lambda.min

  pred_te <- predict(fit, newx = X_te, s = "lambda.min")
  pred_te <- drop(pred_te)  # n_test x 3

  # Per-virus and overall SE
  err <- (pred_te - Y_te)^2
  list(
    fold = i,
    n_test = nrow(Y_te),
    mse_overall = mean(err),
    mse_per_v   = colMeans(err),
    n_features  = ncol(X_tr)
  )
}

cat("Fitting M4_int across 5 outer folds (parallel, 5 cores)...\n")
results <- mclapply(1:n_folds, fit_one_fold, mc.cores = 5, mc.set.seed = TRUE)

# Aggregate weighted MSE
n_total <- sum(sapply(results, function(r) r$n_test))
mse_overall <- sum(sapply(results, function(r) r$mse_overall * r$n_test)) / n_total

mse_per_v <- sapply(VIRUSES, function(v) {
  sum(sapply(results, function(r) r$mse_per_v[v] * r$n_test)) / n_total
})

# Also load the M2_int MSE from the prior run
m2_int_path <- file.path(OUT_DIR, "E3_interaction_M2prime.RData")
e_prior <- new.env(parent = emptyenv())
load(m2_int_path, envir = e_prior)
m2_int_overall <- e_prior$mse_M2p_overall
m2_int_per_v <- e_prior$mse_M2p_virus_global
stopifnot(length(m2_int_overall)==1L, length(m2_int_per_v)==3L)

improvement_overall <- (m2_int_overall - mse_overall) / m2_int_overall * 100
improvement_per_v <- (m2_int_per_v - mse_per_v) / m2_int_per_v * 100

cat("\n=== M4_int results (titers + clinical + interactions + cluster main effects) ===\n")
cat(sprintf("M4_int overall MSE: %.4f\n", mse_overall))
cat(sprintf("M2_int overall MSE: %.4f\n", m2_int_overall))
cat(sprintf("M2_int -> M4_int improvement: %.2f%%\n", improvement_overall))
cat("\nPer-virus:\n")
for (v in VIRUSES) {
  cat(sprintf("  %s: M2_int=%.3f, M4_int=%.3f, improvement=%.2f%%\n",
              v, m2_int_per_v[v], mse_per_v[v], improvement_per_v[v]))
}
cat(sprintf("\nN features in M4_int: %d (28 interactions + %d clusters)\n",
            results[[1]]$n_features, results[[1]]$n_features - 28))

# Save
m4_int_overall_mse <- mse_overall
m4_int_per_virus_mse <- mse_per_v
m4_int_per_fold <- results
m2int_to_m4int_improvement_pct <- improvement_overall
m2int_to_m4int_per_virus_pct <- improvement_per_v

save(m4_int_overall_mse, m4_int_per_virus_mse, m4_int_per_fold,
     m2int_to_m4int_improvement_pct, m2int_to_m4int_per_virus_pct,
     file = file.path(OUT_DIR, "E3_interaction_M4int.RData"))

out_csv <- data.table(
  model = c("M2_int", "M4_int"),
  overall_mse = c(m2_int_overall, mse_overall),
  h1n1_mse = c(m2_int_per_v["h1n1"], mse_per_v["h1n1"]),
  h3n2_mse = c(m2_int_per_v["h3n2"], mse_per_v["h3n2"]),
  bvic_mse = c(m2_int_per_v["bvic"], mse_per_v["bvic"])
)
fwrite(out_csv, file.path(OUT_DIR, "E3_interaction_M4int.csv"))

cat("\nSaved:", file.path(OUT_DIR, "E3_interaction_M4int.RData"), "\n")
cat("Saved CSV.\n")
