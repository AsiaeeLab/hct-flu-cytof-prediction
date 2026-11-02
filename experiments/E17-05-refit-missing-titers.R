#!/usr/bin/env Rscript
# E17-05: Per-fold refit after removing episodes whose HAI titers are missing.
#
# The counts file codes missing titers as 1 (real values below detection are 5);
# methods/08 now treats any titer below 5 as missing. That removes 2 first-season
# episodes (post-vaccination titers missing) and 6 second-season episodes
# (3 missing pre-vaccination, 3 missing post-vaccination titers).
#
# The per-fold atlases of E14-01 (k-means micro-clusters and Leiden candidates)
# are built from cells only and never see an outcome, so they are kept as they
# are. Everything that uses an outcome is refitted on the remaining episodes,
# with the E14-01 seeds and settings:
#   1. inner 3-fold cross-validation at every resolution of the grid, recording
#      both the composite criterion (single-outcome M4 elastic net on the mean log
#      fold-rise versus OLS on titers and clinical factors) and the multitask
#      criterion (as in E15-01);
#   2. held-out M3 and M4 predictions at every resolution (diagnostic, E15-03);
#   3. selection of the resolution minimizing the composite criterion (as E16-01);
#   4. M1 to M4 (multitask and composite), OLS baselines and second-season scores
#      at the selected resolution.
# Output: <FLU_RESULTS_ROOT>/per_fold_atlas_5000/fold_f/fold_f_all.RData with the
# E16-01 object names, and post_hoc_5000/resolution_sweep/ as in E15-01.
suppressPackageStartupMessages({ library(here); library(data.table); library(glmnet) })
OLD_ROOT <- "/fast/datasets/flu_cellular_prediction/first_episode_2026-10-14"
NEW_ROOT <- Sys.getenv("FLU_RESULTS_ROOT"); stopifnot(nzchar(NEW_ROOT), NEW_ROOT != OLD_ROOT)
PRIMARY_ROOT <- "/fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06"
VIRUSES <- c("h1n1", "h3n2", "bvic"); N_OUTER_FOLDS <- 5L; N_INNER_FOLDS <- 3L; ALPHA <- 0.9
FOLD_SEED <- 7000L; SEED_BASE <- 7100L; LAMBDA_RULE <- "lambda.min"; FEATURE_TRANSFORM <- "clr"
RESOLUTION_GRID <- seq(0.5, 5.0, by = 0.25)
N_WORKERS <- as.integer(Sys.getenv("N_WORKERS", "5"))
SWEEP <- file.path(NEW_ROOT, "post_hoc_5000/resolution_sweep"); dir.create(SWEEP, recursive = TRUE, showWarnings = FALSE)
source(here("methods", "00-paths.R")); source(here("methods", "08-multi-outcome-covariates.R"))
source(here("methods", "35-per-fold-atlas-cv.R"))

build_baseline_mat <- function(p_ids, outcomes, viruses) {
  mat <- vapply(viruses, function(v) as.numeric(outcomes[[v]]$x[p_ids]), numeric(length(p_ids)))
  dim(mat) <- c(length(p_ids), length(viruses)); rownames(mat) <- p_ids; colnames(mat) <- paste0("baseline_", viruses); mat
}
build_y_mat <- function(p_ids, outcomes, viruses) {
  mat <- vapply(viruses, function(v) as.numeric(outcomes[[v]]$y[p_ids]), numeric(length(p_ids)))
  dim(mat) <- c(length(p_ids), length(viruses)); rownames(mat) <- p_ids; colnames(mat) <- viruses; mat
}
extract_mgaussian_coefs <- function(fit_obj, viruses) {
  cl <- coef(fit_obj$model, s = fit_obj$model$lambda.min); names(cl) <- viruses
  lapply(cl, function(m) setNames(as.numeric(m), rownames(m)))
}
en_comp <- function(Xtr, ytr, Xte, foldid, seed, Xrep = NULL) {
  set.seed(seed)
  fit <- glmnet::cv.glmnet(Xtr, ytr, family = "gaussian", alpha = ALPHA, foldid = foldid, standardize = TRUE)
  list(test = as.vector(predict(fit, newx = Xte, s = LAMBDA_RULE)),
       repeat_ = if (is.null(Xrep)) NULL else as.vector(predict(fit, newx = Xrep, s = LAMBDA_RULE)),
       coef = as.matrix(coef(fit, s = LAMBDA_RULE)), lambda = fit[[LAMBDA_RULE]])
}

# ---- Episodes with complete titers (methods/08 drops titers coded as 1) -------
cache_agg <- new.env(); load(file.path(PRIMARY_ROOT, "per_fold_atlas_5000/aggregated/E2_aggregated.RData"), cache_agg)
outcome_obj <- read.csv(file.path(paths$raw, "transplant_cytof_counts_hai_12MAY2024.csv"), stringsAsFactors = FALSE)
meta_obj <- read.csv(file.path(paths$raw, "hct_limited-dataset_2025-07-09.csv"), stringsAsFactors = FALSE)
outcomes <- setNames(lapply(VIRUSES, function(v) prepare_outcome_with_covariates(
  outcome_obj = outcome_obj, meta_obj = meta_obj, outcome_col_1 = paste0("geo_mean_", v, "_v1"),
  outcome_col_2 = paste0("geo_mean_", v, "_v3"), p_id = cache_agg$common_p_ids,
  include_dose = TRUE, include_age = TRUE, include_sex = TRUE, include_time = TRUE)), VIRUSES)
complete_ids <- Reduce(intersect, lapply(outcomes, function(o) names(o$y)))
covariates <- outcomes[[VIRUSES[1]]]$covariates[complete_ids, , drop = FALSE]
group_ids <- setNames(base_pid_from_p_id_yr(complete_ids), complete_ids)
first_all <- unlist(lapply(1:N_OUTER_FOLDS, function(f) { e <- new.env()
  load(file.path(OLD_ROOT, sprintf("per_fold_atlas_5000/fold_%d/fold_%d_all.RData", f, f)), e); e$test_p_ids }))
fold_strata <- infer_cohort_from_p_id(intersect(first_all, complete_ids))
cat(sprintf("Complete episodes: %d; first-season kept %d of %d; dropped: %s\n", length(complete_ids),
            length(intersect(first_all, complete_ids)), length(first_all),
            paste(setdiff(first_all, complete_ids), collapse = ", ")))

run_fold <- function(fold_id) {
  e <- new.env(); load(file.path(OLD_ROOT, sprintf("per_fold_atlas_5000/fold_%d/fold_%d_all.RData", fold_id, fold_id)), e)
  cands <- readRDS(file.path(OLD_ROOT, sprintf("per_fold_atlas_5000/fold_%d/resolution_candidates.rds", fold_id)))
  train_p_ids <- intersect(e$train_p_ids, complete_ids); test_p_ids <- intersect(e$test_p_ids, complete_ids)
  repeat_p_ids <- intersect(e$repeat_p_ids, complete_ids)
  dropped <- c(setdiff(e$train_p_ids, train_p_ids), setdiff(e$test_p_ids, test_p_ids), setdiff(e$repeat_p_ids, repeat_p_ids))
  train_micro_counts <- e$train_micro_counts[train_p_ids, , drop = FALSE]
  test_micro_counts <- e$test_micro_counts[test_p_ids, , drop = FALSE]
  repeat_micro_counts <- e$repeat_micro_counts[repeat_p_ids, , drop = FALSE]
  inner_folds <- create_group_folds(p_ids = train_p_ids, group_ids = group_ids, k = N_INNER_FOLDS,
    seed = FOLD_SEED + fold_id * 100, strata = fold_strata[train_p_ids])
  Ytr_all <- build_y_mat(train_p_ids, outcomes, VIRUSES)
  Xbase_all <- as.matrix(cbind(build_baseline_mat(train_p_ids, outcomes, VIRUSES), covariates[train_p_ids, , drop = FALSE]))
  base_inner <- lapply(seq_along(inner_folds), function(ii) {
    val_p <- train_p_ids[inner_folds[[ii]]]; tr_p <- setdiff(train_p_ids, val_p)
    bf <- fit_ols_baseline_multi(Xbase_all[tr_p, , drop = FALSE], Ytr_all[tr_p, , drop = FALSE], Xbase_all[val_p, , drop = FALSE])
    yv <- Ytr_all[val_p, , drop = FALSE]
    oc <- lm.fit(cbind(1, Xbase_all[tr_p, , drop = FALSE]), rowMeans(Ytr_all[tr_p, , drop = FALSE]))
    cf <- oc$coefficients; cf[is.na(cf)] <- 0
    pc <- drop(cbind(1, Xbase_all[val_p, , drop = FALSE]) %*% cf)
    list(mse = colMeans((bf$pred_test - yv)^2), mse_comp = mean((pc - rowMeans(yv))^2))
  })
  outer_foldid <- compute_grouped_foldid(train_p_ids, n_folds = 5L, seed = SEED_BASE + fold_id * 1000 + 99)
  inner <- list(); outer <- list()
  for (r_idx in seq_along(RESOLUTION_GRID)) {
    m2m <- cands[[r_idx]]$micro_to_meta; k <- as.integer(cands[[r_idx]]$k_found)
    clr_tr <- transform_counts(two_stage_collapse_micro_to_meta_counts(train_micro_counts, m2m), method = FEATURE_TRANSFORM)
    clr_te <- transform_counts(two_stage_collapse_micro_to_meta_counts(test_micro_counts, m2m), method = FEATURE_TRANSFORM)
    for (ii in seq_along(inner_folds)) {
      val_p <- train_p_ids[inner_folds[[ii]]]; tr_p <- setdiff(train_p_ids, val_p)
      Xtr <- as.matrix(cbind(clr_tr[tr_p, , drop = FALSE], Xbase_all[tr_p, , drop = FALSE]))
      Xva <- as.matrix(cbind(clr_tr[val_p, , drop = FALSE], Xbase_all[val_p, , drop = FALSE]))
      fid <- compute_grouped_foldid(tr_p, n_folds = 5L, seed = SEED_BASE + fold_id * 10000 + ii)
      seed_fit <- SEED_BASE + fold_id * 10000 + r_idx * 100 + ii * 10
      mt <- fit_predict_multitask_en(X_train = Xtr, Y_train = Ytr_all[tr_p, , drop = FALSE], X_test = Xva,
        seed = seed_fit + 1, alpha = ALPHA, n_cores = 1, foldid = fid, lambda_rule = LAMBDA_RULE)
      yv <- Ytr_all[val_p, , drop = FALSE]
      pc <- en_comp(Xtr, rowMeans(Ytr_all[tr_p, , drop = FALSE]), Xva, fid, seed_fit + 2)$test
      inner[[length(inner) + 1]] <- data.table(fold = fold_id, res = RESOLUTION_GRID[r_idx], k = k, inner = ii,
        crit_multi = median(colMeans((mt$pred_test - yv)^2) / base_inner[[ii]]$mse),
        crit_comp = mean((pc - rowMeans(yv))^2) / base_inner[[ii]]$mse_comp)
    }
    X3tr <- as.matrix(cbind(clr_tr[train_p_ids, , drop = FALSE], build_baseline_mat(train_p_ids, outcomes, VIRUSES)))
    X3te <- as.matrix(cbind(clr_te[test_p_ids, , drop = FALSE], build_baseline_mat(test_p_ids, outcomes, VIRUSES)))
    X4tr <- as.matrix(cbind(X3tr, covariates[train_p_ids, , drop = FALSE])); X4te <- as.matrix(cbind(X3te, covariates[test_p_ids, , drop = FALSE]))
    sb <- SEED_BASE + fold_id * 1000 + r_idx * 50
    f3 <- fit_predict_multitask_en(X3tr, Ytr_all, X3te, seed = sb + 3, alpha = ALPHA, n_cores = 1, foldid = outer_foldid, lambda_rule = LAMBDA_RULE)
    f4 <- fit_predict_multitask_en(X4tr, Ytr_all, X4te, seed = sb + 4, alpha = ALPHA, n_cores = 1, foldid = outer_foldid, lambda_rule = LAMBDA_RULE)
    c3 <- en_comp(X3tr, rowMeans(Ytr_all), X3te, outer_foldid, sb + 13)$test
    c4 <- en_comp(X4tr, rowMeans(Ytr_all), X4te, outer_foldid, sb + 14)$test
    outer[[r_idx]] <- data.table(fold = fold_id, res = RESOLUTION_GRID[r_idx], k = k, episode_id = test_p_ids,
      M3_h1n1 = f3$pred_test[, 1], M3_h3n2 = f3$pred_test[, 2], M3_bvic = f3$pred_test[, 3],
      M4_h1n1 = f4$pred_test[, 1], M4_h3n2 = f4$pred_test[, 2], M4_bvic = f4$pred_test[, 3], M3_comp = c3, M4_comp = c4)
  }
  inner <- rbindlist(inner)
  crit <- inner[, .(comp = mean(crit_comp), multi = mean(crit_multi)), by = res][order(res)]
  best_r_idx <- which.min(crit$comp); selected_resolution <- RESOLUTION_GRID[best_r_idx]
  k_found_by_res <- e$k_found_by_res; k_found <- as.integer(cands[[best_r_idx]]$k_found)
  micro_to_meta <- cands[[best_r_idx]]$micro_to_meta
  resolution_scores <- setNames(crit$comp, as.character(RESOLUTION_GRID))
  resolution_scores_multitask <- setNames(crit$multi, as.character(RESOLUTION_GRID))
  selection_criterion <- "composite"
  # ---- Final models at the selected resolution (E14-01 steps 4 and 5) ----------
  train_meta_counts <- two_stage_collapse_micro_to_meta_counts(train_micro_counts, micro_to_meta)
  test_meta_counts <- two_stage_collapse_micro_to_meta_counts(test_micro_counts, micro_to_meta)
  repeat_meta_counts <- two_stage_collapse_micro_to_meta_counts(repeat_micro_counts, micro_to_meta)
  clr_train <- transform_counts(train_meta_counts, method = FEATURE_TRANSFORM)
  clr_test <- transform_counts(test_meta_counts, method = FEATURE_TRANSFORM)
  clr_repeat <- transform_counts(repeat_meta_counts, method = FEATURE_TRANSFORM)
  bt <- build_baseline_mat(train_p_ids, outcomes, VIRUSES); bs <- build_baseline_mat(test_p_ids, outcomes, VIRUSES)
  br <- build_baseline_mat(repeat_p_ids, outcomes, VIRUSES)
  ct <- covariates[train_p_ids, , drop = FALSE]; cs <- covariates[test_p_ids, , drop = FALSE]; cr <- covariates[repeat_p_ids, , drop = FALSE]
  X_M1_train <- bt; X_M1_test <- bs; X_M1_repeat <- br
  X_M2_train <- as.matrix(cbind(bt, ct)); X_M2_test <- as.matrix(cbind(bs, cs)); X_M2_repeat <- as.matrix(cbind(br, cr))
  X_M3_train <- as.matrix(cbind(clr_train[train_p_ids, , drop = FALSE], bt)); X_M3_test <- as.matrix(cbind(clr_test[test_p_ids, , drop = FALSE], bs))
  X_M3_repeat <- as.matrix(cbind(clr_repeat[repeat_p_ids, , drop = FALSE], br))
  X_M4_train <- as.matrix(cbind(X_M3_train, ct)); X_M4_test <- as.matrix(cbind(X_M3_test, cs)); X_M4_repeat <- as.matrix(cbind(X_M3_repeat, cr))
  Y_train <- build_y_mat(train_p_ids, outcomes, VIRUSES); Y_test <- build_y_mat(test_p_ids, outcomes, VIRUSES)
  Y_repeat <- build_y_mat(repeat_p_ids, outcomes, VIRUSES)
  seed_fit_base <- SEED_BASE + fold_id * 1000 + best_r_idx * 50
  mt <- function(Xtr, Xte, s) fit_predict_multitask_en(X_train = Xtr, Y_train = Y_train, X_test = Xte, seed = s,
    alpha = ALPHA, n_cores = 1, foldid = outer_foldid, lambda_rule = LAMBDA_RULE)
  fit_M1_en <- mt(X_M1_train, X_M1_test, seed_fit_base + 1); fit_M2_en <- mt(X_M2_train, X_M2_test, seed_fit_base + 2)
  fit_M3_en <- mt(X_M3_train, X_M3_test, seed_fit_base + 3); fit_M4_en <- mt(X_M4_train, X_M4_test, seed_fit_base + 4)
  fit_M1_ols <- fit_ols_baseline_multi(X_M1_train, Y_train, X_M1_test); fit_M2_ols <- fit_ols_baseline_multi(X_M2_train, Y_train, X_M2_test)
  nm <- function(m, ids) { rownames(m) <- ids; m }
  pred_M1_en_train <- nm(fit_M1_en$pred_train, train_p_ids); pred_M1_en_test <- nm(fit_M1_en$pred_test, test_p_ids)
  pred_M2_en_train <- nm(fit_M2_en$pred_train, train_p_ids); pred_M2_en_test <- nm(fit_M2_en$pred_test, test_p_ids)
  pred_M3_en_train <- nm(fit_M3_en$pred_train, train_p_ids); pred_M3_en_test <- nm(fit_M3_en$pred_test, test_p_ids)
  pred_M4_en_train <- nm(fit_M4_en$pred_train, train_p_ids); pred_M4_en_test <- nm(fit_M4_en$pred_test, test_p_ids)
  pred_M1_ols_train <- nm(fit_M1_ols$pred_train, train_p_ids); pred_M1_ols_test <- nm(fit_M1_ols$pred_test, test_p_ids)
  pred_M2_ols_train <- nm(fit_M2_ols$pred_train, train_p_ids); pred_M2_ols_test <- nm(fit_M2_ols$pred_test, test_p_ids)
  coef_M1_en <- extract_mgaussian_coefs(fit_M1_en, VIRUSES); coef_M2_en <- extract_mgaussian_coefs(fit_M2_en, VIRUSES)
  coef_M3_en <- extract_mgaussian_coefs(fit_M3_en, VIRUSES); coef_M4_en <- extract_mgaussian_coefs(fit_M4_en, VIRUSES)
  lambda_M1 <- fit_M1_en$model$lambda.min; lambda_M2 <- fit_M2_en$model$lambda.min
  lambda_M3 <- fit_M3_en$model$lambda.min; lambda_M4 <- fit_M4_en$model$lambda.min
  pr <- function(fit, X) { m <- as_pred_matrix(predict(fit$model, newx = X, s = fit$model$lambda.min), VIRUSES); rownames(m) <- rownames(X); m }
  pred_M1_en_repeat <- pr(fit_M1_en, X_M1_repeat); pred_M2_en_repeat <- pr(fit_M2_en, X_M2_repeat)
  pred_M3_en_repeat <- pr(fit_M3_en, X_M3_repeat); pred_M4_en_repeat <- pr(fit_M4_en, X_M4_repeat)
  cf <- function(Xtr, Xte, Xrep, s) en_comp(Xtr, rowMeans(Y_train), Xte, outer_foldid, s, Xrep)
  comp_M1 <- cf(X_M1_train, X_M1_test, X_M1_repeat, seed_fit_base + 11); comp_M2 <- cf(X_M2_train, X_M2_test, X_M2_repeat, seed_fit_base + 12)
  comp_M3 <- cf(X_M3_train, X_M3_test, X_M3_repeat, seed_fit_base + 13); comp_M4 <- cf(X_M4_train, X_M4_test, X_M4_repeat, seed_fit_base + 14)
  # Held-out predictions at the selected resolution must match the sweep
  o <- outer[[best_r_idx]]; stopifnot(max(abs(o$M3_comp - comp_M3$test)) < 1e-8, max(abs(o$M4_comp - comp_M4$test)) < 1e-8,
                                       max(abs(o$M3_h1n1 - pred_M3_en_test[, 1])) < 1e-8)
  fold_dir <- file.path(NEW_ROOT, "per_fold_atlas_5000", sprintf("fold_%d", fold_id)); dir.create(fold_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(OLD_ROOT, sprintf("per_fold_atlas_5000/fold_%d/resolution_candidates.rds", fold_id)), fold_dir, overwrite = TRUE)
  timing_source_md5 <- e$timing_source_md5
  atlas_centroids <- e$atlas_centroids; atlas_scale_center <- e$atlas_scale_center; atlas_scale_scale <- e$atlas_scale_scale
  atlas_centroid_graph <- e$atlas_centroid_graph; atlas_timing_s <- e$atlas_timing_s
  save(timing_source_md5, fold_id, train_p_ids, test_p_ids, dropped,
       atlas_centroids, atlas_scale_center, atlas_scale_scale, atlas_centroid_graph, atlas_timing_s,
       selected_resolution, k_found, resolution_scores, resolution_scores_multitask, selection_criterion, k_found_by_res, micro_to_meta,
       train_micro_counts, test_micro_counts, train_meta_counts, test_meta_counts, clr_train, clr_test,
       X_M1_train, X_M1_test, X_M2_train, X_M2_test, X_M3_train, X_M3_test, X_M4_train, X_M4_test, Y_train, Y_test,
       fit_M1_en, fit_M2_en, fit_M3_en, fit_M4_en,
       pred_M1_en_train, pred_M1_en_test, pred_M2_en_train, pred_M2_en_test, pred_M3_en_train, pred_M3_en_test,
       pred_M4_en_train, pred_M4_en_test, pred_M1_ols_train, pred_M1_ols_test, pred_M2_ols_train, pred_M2_ols_test,
       coef_M1_en, coef_M2_en, coef_M3_en, coef_M4_en, lambda_M1, lambda_M2, lambda_M3, lambda_M4,
       repeat_p_ids, repeat_micro_counts, repeat_meta_counts, clr_repeat,
       X_M1_repeat, X_M2_repeat, X_M3_repeat, X_M4_repeat, Y_repeat,
       pred_M1_en_repeat, pred_M2_en_repeat, pred_M3_en_repeat, pred_M4_en_repeat,
       comp_M1, comp_M2, comp_M3, comp_M4,
       file = file.path(fold_dir, sprintf("fold_%d_all.RData", fold_id)))
  ref <- data.table(fold = fold_id, episode_id = test_p_ids,
    M1_h1n1 = pred_M1_en_test[, 1], M1_h3n2 = pred_M1_en_test[, 2], M1_bvic = pred_M1_en_test[, 3],
    M2_h1n1 = pred_M2_en_test[, 1], M2_h3n2 = pred_M2_en_test[, 2], M2_bvic = pred_M2_en_test[, 3],
    M1_comp = comp_M1$test, M2_comp = comp_M2$test)
  cat(sprintf("Fold %d: dropped %s; train %d, test %d, second season %d; resolution %.2f (k = %d)\n", fold_id,
              paste(dropped, collapse = " "), length(train_p_ids), length(test_p_ids), length(repeat_p_ids), selected_resolution, k_found))
  list(inner = inner, outer = rbindlist(outer), ref = ref)
}
t0 <- Sys.time()
res <- parallel::mclapply(seq_len(N_OUTER_FOLDS), function(f) tryCatch(run_fold(f), error = function(err) err), mc.cores = N_WORKERS)
bad <- vapply(res, inherits, logical(1), "error")
if (any(bad)) stop("Fold(s) failed: ", paste(which(bad), collapse = ", "), "; ", conditionMessage(res[[which(bad)[1]]]))
first_ids <- intersect(first_all, complete_ids); Y <- build_y_mat(first_ids, outcomes, VIRUSES)
fwrite(rbindlist(lapply(res, `[[`, "inner")), file.path(SWEEP, "inner_criteria_by_fold.csv"))
saveRDS(list(outer = rbindlist(lapply(res, `[[`, "outer")), ref = rbindlist(lapply(res, `[[`, "ref")),
             y = data.table(episode_id = first_ids, y_h1n1 = Y[, 1], y_h3n2 = Y[, 2], y_bvic = Y[, 3], y_comp = rowMeans(Y))),
        file.path(SWEEP, "outer_predictions_by_resolution.rds"))
cat(sprintf("Done in %.1f min.\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
