#!/usr/bin/env Rscript
# E14-11: Descriptive numbers quoted in the first-season manuscript that no
# other script prints. Read-only; prints to the console (logs/E14-numbers.log).
suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("FLU_RESULTS_ROOT"); PH <- file.path(ROOT, "post_hoc_5000")
a <- new.env(); load(file.path(ROOT, "per_fold_atlas_5000/aggregated/E2_aggregated.RData"), a)
cat("== aggregated objects:", paste(ls(a), collapse = ", "), "\n")
Y <- a$oof_Y; X0 <- if (!is.null(a$oof_baseline)) a$oof_baseline else NULL
cat("\n== Correlations of antigen log fold-rises (n =", nrow(Y), ")\n"); print(round(cor(Y), 2))
cs <- readRDS(file.path(PH, "clinical_score/episode_table.rds"))
cat("\n== episode table columns:", paste(names(cs), collapse = ", "), "\n")
bcols <- grep("^base_(h1n1|h3n2|bvic)$", names(cs), value = TRUE)
if (length(bcols)) { cat("Geometric mean baseline titers:\n"); print(round(exp(colMeans(cs[, ..bcols])), 1)) }
cat("Months since transplant: "); print(summary(cs$months)); cat("missing months:", sum(is.na(cs$months)), "\n")
# Per-episode prediction benefit (M1 minus M3 squared error, mean over antigens) vs months
P1 <- a$oof_pred_M1_en; P3 <- a$oof_pred_M3_en
if (!is.null(P1) && !is.null(P3)) {
  ben <- rowMeans((Y - P1)^2) - rowMeans((Y - P3)^2)
  ids <- rownames(Y); mo <- cs$months[match(ids, cs$episode_id)]
  sl <- coef(lm(ben ~ mo))[2]; set.seed(20260728)
  bs <- replicate(2000, { i <- sample(length(ben), replace = TRUE); coef(lm(ben[i] ~ mo[i]))[2] })
  cat(sprintf("\n== Benefit slope per month: %.4f (95%% CI %.4f to %.4f)\n", sl, quantile(bs, .025), quantile(bs, .975)))
} else cat("\n(no oof_pred M1_en/M3_en in aggregated file)\n")
# Raw composite score by fold
cmp <- readRDS(file.path(PH, "composite/composite_episode_table.rds"))
cat("\n== composite table columns:", paste(names(cmp), collapse = ", "), "\n")
sc <- intersect(c("score_comp", "s_comp", "comp_score", "score"), names(cmp))[1]
fc <- intersect(c("fold", "outer_fold"), names(cmp))[1]
if (!is.na(sc) && !is.na(fc)) print(cmp[, .(mean = round(mean(get(sc)), 2), sd = round(sd(get(sc)), 2), n = .N), by = fc][order(get(fc))])
# Shared-atlas bootstrap: stability and baseline coefficients
e <- new.env(); load(Sys.getenv("FLU_BOOT_FILE"), e)
cat("\n== bootstrap objects:", paste(ls(e), collapse = ", "), "\n")
st <- as.data.table(e$stability)
print(st[!grepl("^cluster_", feature), .(sel = round(max(selection_freq), 3)), by = feature])
cl <- st[grepl("^cluster_", feature), .(sel = max(selection_freq)), by = feature][order(-sel)]
cat("Populations selected in >= 50% of refits:", sum(cl$sel >= .5), "of", nrow(cl), "\n"); print(head(cl, 6))
if (!is.null(e$lambda_min)) cat("lambda_min:", e$lambda_min, "\n")
cd <- as.data.table(e$coef_dt)
print(cd[feature == paste0("baseline_", virus), .(mean_same_antigen_coef = round(mean(coef), 3)), by = virus])
# Pooled and per-cohort signs
ps <- fread(file.path(PH, "predictor_signs/pooled_signs.csv"))
ps[, avg := (mean_coef_h1n1 + mean_coef_h3n2 + mean_coef_bvic) / 3]
cat("\n== Pooled signs (stably selected, any antigen >= 50%)\n")
print(ps[pmax(selection_freq_h1n1, selection_freq_h3n2, selection_freq_bvic) >= .5][order(-avg),
         .(cluster_id, avg = round(avg, 3), sel = pmax(selection_freq_h1n1, selection_freq_h3n2, selection_freq_bvic), direction, annotation)])
ann <- fread(Sys.getenv("FLU_ANN_FILE"))
bc <- ann[broad_lineage == "B cell", cluster_id]
cat("B-cell populations and pooled average coefficient:\n"); print(ps[cluster_id %in% bc, .(cluster_id, avg = round(avg, 3), annotation)][order(-avg)])
cat("Plasmablast and cTfh:\n"); print(ps[cluster_id %in% c(28, 42, 25), .(cluster_id, round(mean_coef_h1n1, 3), round(mean_coef_h3n2, 3), round(mean_coef_bvic, 3), direction)])
pc <- fread(file.path(PH, "predictor_signs/per_cohort_signs.csv"))
cat("\n== Per-cohort signs, populations selected >= 70% in both cohorts\n")
print(pc[freq_adult >= .7 & freq_peds >= .7, .(cluster_id, freq_peds, freq_adult, peds = round(imp_mean_peds, 3), adult = round(imp_mean_adult, 3), same_sign, annotation)])
# Top-k ranks by fold
rk <- fread(file.path(PH, "top_k_populations/population_ranks_by_fold.csv"))
cat("\n== Populations in the top 6 by fold\n"); print(rk[rank <= 6, .(folds = .N, folds_list = paste(fold, collapse = ",")), by = population][order(-folds)])
# Example episodes
exf <- file.path(PH, "example_episodes/example_titers.csv")
if (file.exists(exf)) { ex <- fread(exf); cat("\n== Example episodes\n"); print(ex) }
# Per-fold selected k
cat("\n== selected k by fold:", paste(unlist(a$selected_k_per_fold), collapse = ", "), "; res:", paste(unlist(a$selected_res_per_fold), collapse = ", "), "\n")
for (f in 1:5) { fe <- new.env(); ff <- file.path(ROOT, "per_fold_atlas_5000", paste0("fold_", f))
  rf <- list.files(ff, pattern = "\\.RData$", full.names = TRUE)[1]; if (is.na(rf)) next; load(rf, fe)
  cat(sprintf("fold %d: k = %s, resolution = %s\n", f, paste(fe$selected_k, collapse = "/"), paste(fe$selected_res, collapse = "/"))) }
# Repeat-episode season interval
rv <- readRDS(file.path(PH, "repeater_validation/repeater_pairs.rds"))
cat("\n== repeater pair columns:", paste(names(rv), collapse = ", "), "\n")
mc <- grep("months", names(rv), value = TRUE)
if (length(mc) >= 2) { d <- rv[[mc[2]]] - rv[[mc[1]]]; cat("Months between seasons:"); print(summary(d)) }
cat("Cohorts among pairs:"); if ("cohort" %in% names(rv)) print(table(rv$cohort))
