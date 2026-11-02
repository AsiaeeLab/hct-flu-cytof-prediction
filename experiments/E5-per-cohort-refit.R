# =============================================================================
# E5: Per-cohort refit of M1 to M4 with recipient-clustered bootstrap inference.
#
# Question: does the pooling of pediatric and adult HCT cohorts hold up under
# a within-cohort analysis? Refit M1 to M4 separately on each cohort using the
# same production-atlas CLR features, then compare per-cohort percent
# improvement and selection frequency to the pooled result.
#
# Frozen atlas: production atlas (33 meta-clusters), shared across cohorts.
# This is the apples-to-apples comparison; it isolates the prediction-model
# refit per cohort while holding the feature construction identical.
# =============================================================================

FLU_ROOT <- Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06")
FLU_BIO_FILE <- Sys.getenv("FLU_BIO_FILE", "/home/amir/repo/flu-cellular-prediction/results/final_atlas_production/biology_production_90203bf3.RData")
FLU_ANN_FILE <- Sys.getenv("FLU_ANN_FILE", "/home/amir/repo/flu-cellular-prediction/results/reports/cluster_annotations_corrected.csv")
suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(parallel)
})

set.seed(7000)

OUT_DIR <- file.path(FLU_ROOT, "post_hoc_5000/per_cohort")
FIG_DIR <- "/home/amir/repo/flu-cellular-prediction/paper/images"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

BIO_FILE  <- FLU_BIO_FILE
FOLD_FILE <- file.path(FLU_ROOT, "per_fold_atlas_5000/fold_1/fold_1_all.RData")
ANN_FILE  <- FLU_ANN_FILE

# ---------------------------------------------------------------------------
# Load shared feature matrix, Y, and clinical covariates
# ---------------------------------------------------------------------------
e_bio  <- new.env(parent = emptyenv()); load(BIO_FILE,  envir = e_bio)
e_fold <- new.env(parent = emptyenv()); load(FOLD_FILE, envir = e_fold)

ids       <- e_bio$common_p_ids                  # 343
cohort    <- e_bio$cohort_vec                    # "peds" / "adult"
clr_feat  <- as.matrix(e_bio$clr_features)       # 343 x 33

all_Y     <- rbind(e_fold$Y_train, e_fold$Y_test)
all_X_M1  <- rbind(e_fold$X_M1_train, e_fold$X_M1_test)
all_X_M2  <- rbind(e_fold$X_M2_train, e_fold$X_M2_test)
all_X_M3  <- cbind(clr_feat[rownames(all_Y), , drop = FALSE], all_X_M1)
all_X_M4  <- cbind(clr_feat[rownames(all_Y), , drop = FALSE], all_X_M2)

ord <- match(ids, rownames(all_Y))
all_Y    <- all_Y[ord, , drop = FALSE]
all_X_M1 <- all_X_M1[ord, , drop = FALSE]
all_X_M2 <- all_X_M2[ord, , drop = FALSE]
all_X_M3 <- all_X_M3[ord, , drop = FALSE]
all_X_M4 <- all_X_M4[ord, , drop = FALSE]
stopifnot(identical(rownames(all_Y), ids))

base_id <- sub("_yr[0-9]+$", "", ids)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
make_recipient_folds <- function(base_recipients, k = 5, seed = 7000) {
  set.seed(seed)
  uniq <- unique(base_recipients)
  fold <- sample(rep(seq_len(k), length.out = length(uniq)))
  names(fold) <- uniq
  unname(fold[base_recipients])
}

inner_grouped_foldid <- function(p_ids, k = 3, seed) {
  set.seed(seed)
  uniq <- unique(p_ids)
  fold_assign <- sample(rep(seq_len(k), length.out = length(uniq)))
  names(fold_assign) <- uniq
  unname(fold_assign[p_ids])
}

refit_one_outer_fold <- function(fold_id, fold_assign, X, Y, p_ids, seed_inner) {
  test_idx  <- which(fold_assign == fold_id)
  train_idx <- which(fold_assign != fold_id)
  X_tr <- X[train_idx, , drop = FALSE]
  X_te <- X[test_idx,  , drop = FALSE]
  Y_tr <- Y[train_idx, , drop = FALSE]
  Y_te <- Y[test_idx,  , drop = FALSE]
  fid  <- inner_grouped_foldid(p_ids[train_idx], k = 3, seed = seed_inner)

  fit  <- cv.glmnet(X_tr, Y_tr, family = "mgaussian", alpha = 0.9,
                    foldid = fid, standardize = TRUE)
  pred_te <- drop(predict(fit, newx = X_te, s = "lambda.min"))
  list(test_ids = rownames(X_te),
       test_p_ids = p_ids[test_idx],
       Y_test = Y_te,
       pred = pred_te,
       n_test = nrow(X_te))
}

run_cohort <- function(cohort_label, mask, n_outer = 5) {
  cat(sprintf("\n=== Cohort: %s (n_episodes = %d) ===\n", cohort_label, sum(mask)))
  ids_c    <- ids[mask]
  base_c   <- base_id[mask]
  Y_c      <- all_Y[mask, , drop = FALSE]
  X1_c     <- all_X_M1[mask, , drop = FALSE]
  X2_c     <- all_X_M2[mask, , drop = FALSE]
  X3_c     <- all_X_M3[mask, , drop = FALSE]
  X4_c     <- all_X_M4[mask, , drop = FALSE]
  fold     <- make_recipient_folds(base_c, k = n_outer, seed = 7000)

  models <- list(M1 = X1_c, M2 = X2_c, M3 = X3_c, M4 = X4_c)
  per_model_pred <- list()
  for (mname in names(models)) {
    X_c <- models[[mname]]
    results <- lapply(seq_len(n_outer), function(f)
      refit_one_outer_fold(f, fold, X_c, Y_c, base_c,
                            seed_inner = 7000 + f * 17))
    pred_all <- do.call(rbind, lapply(results, function(r) r$pred))
    rownames(pred_all) <- unlist(lapply(results, function(r) r$test_ids))
    pred_all <- pred_all[ids_c, , drop = FALSE]
    per_model_pred[[mname]] <- pred_all
  }
  list(ids_c = ids_c, base_c = base_c, Y_c = Y_c, pred = per_model_pred)
}

# Run for adult and peds
res_peds  <- run_cohort("peds",  cohort == "peds")
res_adult <- run_cohort("adult", cohort == "adult")

# ---------------------------------------------------------------------------
# Recipient-clustered bootstrap for M1->M3, M2->M4 per cohort
# ---------------------------------------------------------------------------
delta_pct_recipient_boot <- function(pred_A, pred_B, Y, base_ids,
                                      B = 2000L, seed = 7000) {
  err_A <- (pred_A - Y)^2
  err_B <- (pred_B - Y)^2
  n_rows <- nrow(Y)
  observed <- 100 * (mean(err_A) - mean(err_B)) / mean(err_A)
  uniq <- unique(base_ids)
  rows_by_recip <- split(seq_len(n_rows), base_ids)
  set.seed(seed)
  boot <- replicate(B, {
    smpl <- sample(uniq, length(uniq), replace = TRUE)
    idx <- unlist(rows_by_recip[smpl], use.names = FALSE)
    ea <- err_A[idx, , drop = FALSE]
    eb <- err_B[idx, , drop = FALSE]
    100 * (mean(ea) - mean(eb)) / mean(ea)
  })
  ci <- quantile(boot, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
  p_two_sided <- 2 * min(mean(boot <= 0, na.rm = TRUE),
                          mean(boot >= 0, na.rm = TRUE))
  list(delta_pct = observed, ci_low = ci[1], ci_high = ci[2],
       p_two_sided = p_two_sided, boot = boot)
}

run_inference <- function(res, label) {
  cat(sprintf("\n[%s] M1 -> M3 overall:\n", label))
  d1 <- delta_pct_recipient_boot(res$pred$M1, res$pred$M3, res$Y_c, res$base_c)
  cat(sprintf("  +%.2f%% [%.2f, %.2f], p = %.4f\n",
              d1$delta_pct, d1$ci_low, d1$ci_high, d1$p_two_sided))
  cat(sprintf("[%s] M2 -> M4 overall:\n", label))
  d2 <- delta_pct_recipient_boot(res$pred$M2, res$pred$M4, res$Y_c, res$base_c)
  cat(sprintf("  +%.2f%% [%.2f, %.2f], p = %.4f\n",
              d2$delta_pct, d2$ci_low, d2$ci_high, d2$p_two_sided))
  # Per-antigen for M1 -> M3
  per_antigen <- list()
  for (v in c("h1n1", "h3n2", "bvic")) {
    pa <- res$pred$M1[, v, drop = FALSE]
    pb <- res$pred$M3[, v, drop = FALSE]
    Y  <- res$Y_c[,   v, drop = FALSE]
    d  <- delta_pct_recipient_boot(pa, pb, Y, res$base_c)
    per_antigen[[v]] <- d
    cat(sprintf("[%s] M1 -> M3 %s: +%.2f%% [%.2f, %.2f], p = %.4f\n",
                label, v, d$delta_pct, d$ci_low, d$ci_high, d$p_two_sided))
  }
  list(M1_M3 = d1, M2_M4 = d2, per_antigen = per_antigen)
}

inf_peds  <- run_inference(res_peds,  "peds")
inf_adult <- run_inference(res_adult, "adult")

# ---------------------------------------------------------------------------
# Per-cohort baseline-titer quartile stratification (using per-cohort
# refit predictions). Quartiles computed on the within-cohort baseline.
# ---------------------------------------------------------------------------
stratify_cohort <- function(res, label) {
  Y      <- res$Y_c
  pred_1 <- res$pred$M1
  pred_3 <- res$pred$M3
  # Within-cohort mean baseline log titer per recipient (proxy: just use
  # the M1 features average, since the M1 design includes only baseline titers)
  # We use the X1 baseline log titer matrix saved earlier
  # The recipients' average baseline log titer:
  base_avg <- rowMeans(all_X_M1[match(res$ids_c, ids), , drop = FALSE])
  q <- cut(base_avg, breaks = quantile(base_avg, probs = c(0, 0.25, 0.5, 0.75, 1.0),
                                        na.rm = TRUE),
            include.lowest = TRUE,
            labels = c("Q1", "Q2", "Q3", "Q4"))
  out <- data.table()
  for (qq in levels(q)) {
    idx <- which(q == qq)
    if (length(idx) < 5) next
    mse1 <- mean((pred_1[idx, ] - Y[idx, ])^2)
    mse3 <- mean((pred_3[idx, ] - Y[idx, ])^2)
    out <- rbind(out, data.table(cohort = label, quartile = qq,
                                  n = length(idx),
                                  M1 = mse1, M3 = mse3,
                                  improvement_pct = 100 * (mse1 - mse3) / mse1))
  }
  out
}
strat_peds  <- stratify_cohort(res_peds,  "peds")
strat_adult <- stratify_cohort(res_adult, "adult")
strat_all   <- rbind(strat_peds, strat_adult)
print(strat_all)

# ---------------------------------------------------------------------------
# Per-cohort feature selection frequency (bootstrap)
# ---------------------------------------------------------------------------
boot_selection_freq <- function(X, Y, p_ids, lambda_min,
                                  B = 200L, seed = 7000) {
  set.seed(seed)
  uniq <- unique(p_ids)
  rows_by_recip <- split(seq_len(nrow(X)), p_ids)
  feat_names <- colnames(X)
  counts <- matrix(0, nrow = length(feat_names), ncol = 3,
                    dimnames = list(feat_names, c("h1n1", "h3n2", "bvic")))
  for (b in seq_len(B)) {
    smpl <- sample(uniq, length(uniq), replace = TRUE)
    idx  <- unlist(rows_by_recip[smpl], use.names = FALSE)
    Xb <- X[idx, , drop = FALSE]
    Yb <- Y[idx, , drop = FALSE]
    fit <- glmnet(Xb, Yb, family = "mgaussian", alpha = 0.9,
                   lambda = lambda_min, standardize = TRUE)
    cf <- coef(fit, s = lambda_min)
    for (v in c("h1n1", "h3n2", "bvic")) {
      cv <- as.matrix(cf[[v]])
      cv <- cv[rownames(cv) != "(Intercept)", , drop = FALSE]
      counts[, v] <- counts[, v] + (abs(cv[, 1]) > 0)
    }
  }
  counts / B
}

# Get a reasonable lambda by running cv.glmnet on the full cohort feature set
# for a quick estimate; use the same lambda for both cohorts' selection-freq
# bootstraps so the regularization strength is fixed.
fit_cv_M3 <- cv.glmnet(all_X_M3, all_Y, family = "mgaussian", alpha = 0.9,
                        foldid = inner_grouped_foldid(base_id, 3, 7000),
                        standardize = TRUE)
lambda_ref <- fit_cv_M3$lambda.min
cat(sprintf("\nlambda_ref (from full cohort, used for both subset selection-freq): %.4f\n",
            lambda_ref))

sel_peds  <- boot_selection_freq(all_X_M3[cohort == "peds",  , drop = FALSE],
                                  all_Y[cohort == "peds",   , drop = FALSE],
                                  base_id[cohort == "peds"],  lambda_ref, B = 200)
sel_adult <- boot_selection_freq(all_X_M3[cohort == "adult", , drop = FALSE],
                                  all_Y[cohort == "adult",  , drop = FALSE],
                                  base_id[cohort == "adult"], lambda_ref, B = 200)

# Mean over the three antigens (multitask penalty produces identical sparsity
# across antigens in the pooled fit; in per-antigen single-task we report the
# average for compactness)
sel_summary <- data.table(
  feature   = rownames(sel_peds),
  peds_freq  = rowMeans(sel_peds),
  adult_freq = rowMeans(sel_adult)
)
sel_summary <- sel_summary[order(-pmax(peds_freq, adult_freq))]
cat("\nTop 15 features by max selection frequency across cohorts:\n")
print(head(sel_summary, 15))

# ---------------------------------------------------------------------------
# Save artifacts
# ---------------------------------------------------------------------------
save(res_peds, res_adult, inf_peds, inf_adult, strat_all, sel_summary,
     sel_peds, sel_adult, lambda_ref,
     file = file.path(OUT_DIR, "E5_per_cohort.RData"))
fwrite(strat_all,   file.path(OUT_DIR, "E5_strat.csv"))
fwrite(sel_summary, file.path(OUT_DIR, "E5_selection.csv"))

summary_tbl <- data.table(
  cohort = c("peds", "peds", "adult", "adult"),
  comparison = c("M1_to_M3", "M2_to_M4", "M1_to_M3", "M2_to_M4"),
  delta_pct = c(inf_peds$M1_M3$delta_pct,  inf_peds$M2_M4$delta_pct,
                inf_adult$M1_M3$delta_pct, inf_adult$M2_M4$delta_pct),
  ci_low    = c(inf_peds$M1_M3$ci_low,     inf_peds$M2_M4$ci_low,
                inf_adult$M1_M3$ci_low,    inf_adult$M2_M4$ci_low),
  ci_high   = c(inf_peds$M1_M3$ci_high,    inf_peds$M2_M4$ci_high,
                inf_adult$M1_M3$ci_high,   inf_adult$M2_M4$ci_high),
  p_two_sided = c(inf_peds$M1_M3$p_two_sided,  inf_peds$M2_M4$p_two_sided,
                   inf_adult$M1_M3$p_two_sided, inf_adult$M2_M4$p_two_sided)
)
fwrite(summary_tbl, file.path(OUT_DIR, "E5_summary.csv"))
cat("\n=== Summary saved ===\n"); print(summary_tbl)
