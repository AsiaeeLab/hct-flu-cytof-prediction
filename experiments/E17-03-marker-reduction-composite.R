#!/usr/bin/env Rscript
# E17-03: Marker reduction (Figure 6) refitted on the composite outcome.
# E2-marker-reduction-perfold.R reassigned the sampled cells (up to 2,000 per
# recipient) to each fold's atlas with each reduced marker set and saved the
# per-fold, per-budget population counts; that step does not depend on the
# outcome. Here the M4 elastic net is refitted to the composite response (mean of
# the three antigen log fold-rises) from those saved counts, with the same
# seeds, grouped penalty-tuning folds and alpha as E2. The reference line is the
# composite M1 error of the primary analysis (E12-01).
# Output: post_hoc_5000/marker_reduction/E2_marker_perf_composite.csv, read by
# E2-marker-reduction-perfold.R when OUTCOME=composite.
suppressPackageStartupMessages({ library(here); library(data.table); library(glmnet); library(parallel) })
source(here("methods", "00-paths.R")); source(here("methods", "35-per-fold-atlas-cv.R"))
ROOT <- Sys.getenv("FLU_RESULTS_ROOT"); stopifnot(nzchar(ROOT))
MR <- file.path(ROOT, "post_hoc_5000/marker_reduction")
SEED <- 7000L; ALPHA <- 0.9; BUDGETS <- c(30, 25, 20, 15, 12, 10, 8, 6)
N_CORES <- as.integer(Sys.getenv("N_CORES", "40"))
agg <- new.env(); load(file.path(ROOT, "per_fold_atlas_5000/aggregated/E2_aggregated.RData"), agg)
ids <- agg$common_p_ids
y <- setNames(rowMeans(agg$oof_Y), ids)
jobs <- CJ(fold = 1:5, budget = BUDGETS)
res <- mclapply(seq_len(nrow(jobs)), function(j) {
  f <- jobs$fold[j]; b <- jobs$budget[j]
  env <- new.env(); load(file.path(ROOT, "per_fold_atlas_5000", sprintf("fold_%d/fold_%d_all.RData", f, f)), env)
  Xb <- rbind(env$X_M2_train, env$X_M2_test)[ids, , drop = FALSE]
  cnt <- readRDS(file.path(MR, sprintf("counts_fold%d_budget%d.rds", f, b)))[ids, , drop = FALSE]
  X <- cbind(Xb, transform_counts(cnt, method = "clr")); storage.mode(X) <- "double"
  tr <- match(env$train_p_ids, ids); te <- match(env$test_p_ids, ids)
  set.seed(SEED + b * 10 + f * 100)
  fit <- cv.glmnet(X[tr, ], y[tr], family = "gaussian", alpha = ALPHA, standardize = TRUE,
                   foldid = compute_grouped_foldid(env$train_p_ids, 5L, SEED + f * 100))
  p <- as.vector(predict(fit, newx = X[te, ], s = "lambda.min"))
  data.table(fold = f, budget = b, n_test = length(te), test_mse = mean((p - y[te])^2))
}, mc.cores = min(N_CORES, nrow(jobs)))
perf <- rbindlist(res)
ct <- as.data.table(readRDS(file.path(ROOT, "post_hoc_5000/composite/composite_episode_table.rds")))
m1 <- ct[, mean((y_comp - pred_M1)^2)]
agg_b <- perf[, .(mse = weighted.mean(test_mse, n_test)), by = budget][, `:=`(rmse = sqrt(mse), gain_vs_M1 = 100 * (1 - mse / m1))]
fwrite(perf, file.path(MR, "E2_marker_perf_composite.csv"))
fwrite(cbind(agg_b, m1_mse = m1), file.path(MR, "E2_marker_summary_composite.csv"))
print(agg_b[order(-budget)]); cat(sprintf("Composite M1 MSE %.3f (RMSE %.3f)\n", m1, sqrt(m1)))
