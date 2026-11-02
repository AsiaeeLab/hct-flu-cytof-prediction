#!/usr/bin/env Rscript
# Clinical reading of the cellular composition score (October 2026 round).
#
# Uses the saved out-of-fold scores and predictions of the primary pipeline; no
# model is refitted. Each primary fold has its own atlas, so raw scores carry a
# fold-specific offset. Every analysis here uses the score standardized within
# its fold (z, in standard deviations), which puts all episodes on one scale.
#
# Outputs (episode table and summary tables) feed E13-03 (figures) and E13-04
# (training-aware p-values from the E13-01 null replicates).
#
#   1. Episode table: cohort, first episode, time since transplant, dose, age,
#      routine lymphocyte counts, titers, responses, standardized scores.
#   2. Rank correlation of score with composite fold-rise (Spearman), overall,
#      by cohort and in first episodes; partial rank correlation given baseline.
#   3. Fourfold-response rate by fifth of the score; bottom fifth versus rest.
#   4. Logistic models: response ~ baseline titer + score, then + cohort
#      interaction (does the association differ between children and adults?),
#      then + time since transplant, dose and age.
#   5. ROC: titers-only (M1) versus titers-plus-composition (M3) predictions.
#   6. Pooled-model error reduction computed within each cohort.
#   7. Routine counts (ALC, CD4, CD19) compared with the score.
# Recipient-level bootstrap (B = 2000) for intervals; cluster-robust Wald tests
# by recipient for model terms.
FLU_ROOT <- Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06")
suppressPackageStartupMessages({library(data.table); library(sandwich); library(lmtest); library(pROC)})
ROOT <- FLU_ROOT
OUT <- file.path(ROOT, "post_hoc_5000/clinical_score"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
V <- c("h1n1","h3n2","bvic"); B <- 2000L; L4 <- log(4) - 1e-12; L40 <- log(40) - 1e-12

# ---- 1. Episode table ------------------------------------------------------
d <- as.data.table(readRDS(file.path(ROOT, "post_hoc_5000/composite/composite_episode_table.rds")))
a <- new.env(); load(file.path(ROOT, "per_fold_atlas_5000/aggregated/E2_aggregated.RData"), a)
stopifnot(setequal(a$common_p_ids, d$episode_id))
for (m in c("M1","M2","M3","M4")) { p <- a[[paste0("oof_pred_", m, "_en")]]; rownames(p) <- a$common_p_ids
  for (i in 1:3) d[, (paste0("p", m, "_", V[i])) := p[episode_id, i]] }
meta <- fread("~/datasets/flu/raw/hsct_metadata_15JUL2025.csv")
meta[, episode_id := sub("_([0-9]+)$", "_yr\\1", id_yr)]
meta <- meta[match(d$episode_id, episode_id)]; stopifnot(identical(meta$episode_id, d$episode_id))
num <- function(x) suppressWarnings(as.numeric(x))
d[, `:=`(cohort = factor(ifelse(grepl("^PSCT", episode_id), "Children", "Adults"), c("Children","Adults")),
         months = meta$mo_post, dose = meta$intervention_recode, age = num(meta$age),
         alc = num(meta$alc_visit), cd4 = num(meta$cd4_visit), cd19 = num(meta$cd19_visit),
         yr = as.integer(sub(".*_yr", "", episode_id)))]
d[, first_episode := yr == min(yr), by = base_id]
for (s in c("comp", V)) d[, (paste0("z_", s)) := (get(paste0("score_", s)) - mean(get(paste0("score_", s)))) /
                                               sd(get(paste0("score_", s))), by = fold]
d[, resp_comp := y_comp >= L4]
for (v in V) d[, c(paste0("resp_", v), paste0("prot_", v), paste0("unprot0_", v)) :=
                 list(get(paste0("y_", v)) >= L4, get(paste0("base_", v)) + get(paste0("y_", v)) >= L40,
                      get(paste0("base_", v)) < L40)]
d[, n_resp := resp_h1n1 + resp_h3n2 + resp_bvic]
d[, n_prot := prot_h1n1 + prot_h3n2 + prot_bvic]
d[, fifth := cut(z_comp, quantile(z_comp, 0:5/5), include.lowest = TRUE, labels = paste0("F", 1:5))]
d[, bottom_fifth := fifth == "F1"]
# Within-fold SD of the raw composite score, to translate one SD into fold-rise units
sd_raw <- d[, sd(score_comp), by = fold][, mean(V1)]
saveRDS(d, file.path(OUT, "episode_table.rds"))
cat(sprintf("Episodes %d (children %d, adults %d); first episodes %d; within-fold SD of raw composite score %.3f (x%.2f predicted fold-rise per SD)\n",
            nrow(d), sum(d$cohort=="Children"), sum(d$cohort=="Adults"), sum(d$first_episode), sd_raw, exp(sd_raw)))

# ---- helpers -----------------------------------------------------------------
recips <- unique(d$base_id); rows_by <- split(seq_len(nrow(d)), d$base_id)[recips]
set.seed(20261013L)
boot_idx <- replicate(B, unlist(rows_by[sample.int(length(recips), replace = TRUE)], use.names = FALSE), simplify = FALSE)
boot_ci <- function(fun, idx_subset = seq_len(nrow(d))) {
  keep <- logical(nrow(d)); keep[idx_subset] <- TRUE
  est <- fun(d[keep])
  bs <- vapply(boot_idx, function(ix) { ix <- ix[keep[ix]]; fun(d[ix]) }, numeric(length(est)))
  bs <- matrix(bs, nrow = length(est))
  data.table(estimate = est, ci_low = apply(bs, 1, quantile, .025, na.rm = TRUE),
             ci_high = apply(bs, 1, quantile, .975, na.rm = TRUE),
             p_boot = apply(bs, 1, function(b) min(1, 2 * min(mean(b <= 0, na.rm=TRUE), mean(b >= 0, na.rm=TRUE)))))
}
prank <- function(x) (rank(x) - .5) / length(x)
partial_spearman <- function(x, y, z) cor(resid(lm(prank(x) ~ prank(z))), resid(lm(prank(y) ~ prank(z))))
subsets <- list(All = which(TRUE | d$cohort == "Children"), Children = which(d$cohort == "Children"),
                Adults = which(d$cohort == "Adults"), `First episodes` = which(d$first_episode))

# ---- 2. Rank correlations ----------------------------------------------------
rho <- rbindlist(lapply(names(subsets), function(s) {
  rbind(cbind(subset = s, statistic = "Spearman", boot_ci(function(x) cor(x$z_comp, x$y_comp, method = "spearman"), subsets[[s]])),
        cbind(subset = s, statistic = "Partial Spearman given baseline", boot_ci(function(x) partial_spearman(x$z_comp, x$y_comp, x$base_comp), subsets[[s]])))
}))
rho[, n := lengths(subsets)[subset]]
fwrite(rho, file.path(OUT, "rank_correlations.csv")); print(rho)

# ---- 3. Response by fifth ----------------------------------------------------
fifths <- d[, .(n = .N, responders = sum(resp_comp), rate = mean(resp_comp)), by = .(cohort, fifth)][order(cohort, fifth)]
fifths_all <- d[, .(cohort = "All", n = .N, responders = sum(resp_comp), rate = mean(resp_comp)), by = fifth][order(fifth)]
fifths <- rbind(fifths_all, fifths)
# Wilson intervals (descriptive)
wilson <- function(k, n) { p <- k/n; z <- 1.96; c <- (p + z^2/(2*n)) / (1 + z^2/n); h <- z*sqrt(p*(1-p)/n + z^2/(4*n^2)) / (1 + z^2/n); cbind(c - h, c + h) }
fifths[, c("lo","hi") := as.data.table(wilson(responders, n))]
fwrite(fifths, file.path(OUT, "response_by_fifth.csv")); print(fifths)
bottom <- rbindlist(lapply(names(subsets), function(s) cbind(subset = s,
  boot_ci(function(x) c(mean(x$resp_comp[x$bottom_fifth]), mean(x$resp_comp[!x$bottom_fifth]),
                        mean(x$resp_comp[!x$bottom_fifth]) - mean(x$resp_comp[x$bottom_fifth])), subsets[[s]]),
  quantity = c("rate_bottom_fifth","rate_rest","difference"))))
fwrite(bottom, file.path(OUT, "bottom_fifth.csv")); print(bottom)

# ---- 4. Logistic models --------------------------------------------------------
rob <- function(fit, term) { ct <- coeftest(fit, vcov = vcovCL(fit, cluster = d$base_id[as.integer(rownames(model.frame(fit)))]))
  data.table(term = term, estimate = ct[term, 1], se = ct[term, 2], p = ct[term, 4]) }
d[, row := .I]
lg1 <- glm(resp_comp ~ base_comp + z_comp, binomial, d)
lg2 <- glm(resp_comp ~ base_comp + z_comp * cohort, binomial, d)
lg3 <- glm(resp_comp ~ base_comp + z_comp + log(months) + I(dose == "HD") + cohort + age, binomial, d)
lg4 <- glm(resp_comp ~ base_comp * z_comp, binomial, d)
logit <- rbind(cbind(model = "baseline + score", rob(lg1, "z_comp")),
               cbind(model = "baseline + score x cohort", rob(lg2, "z_comp:cohortAdults")),
               cbind(model = "baseline + score + months + dose + cohort + age", rob(lg3, "z_comp")),
               cbind(model = "baseline + score + months + dose + cohort + age", rob(lg3, "log(months)")),
               cbind(model = "baseline x score", rob(lg4, "base_comp:z_comp")))
for (ch in c("Children","Adults")) { f <- glm(resp_comp ~ base_comp + z_comp, binomial, d[cohort == ch])
  ct <- coeftest(f, vcov = vcovCL(f, cluster = d[cohort == ch]$base_id))
  logit <- rbind(logit, data.table(model = paste("baseline + score,", ch), term = "z_comp", estimate = ct["z_comp",1], se = ct["z_comp",2], p = ct["z_comp",4])) }
logit[, `:=`(or = exp(estimate), or_lo = exp(estimate - 1.96*se), or_hi = exp(estimate + 1.96*se))]
fwrite(logit, file.path(OUT, "logistic_models.csv")); print(logit)
fe <- glm(resp_comp ~ base_comp + z_comp, binomial, d[first_episode == TRUE])
cte <- coeftest(fe, vcov = vcovCL(fe, cluster = d[first_episode == TRUE]$base_id))
fwrite(data.table(model = "baseline + score, first episodes", term = "z_comp", estimate = cte["z_comp",1], se = cte["z_comp",2],
                  p = cte["z_comp",4], or = exp(cte["z_comp",1]), or_lo = exp(cte["z_comp",1] - 1.96*cte["z_comp",2]),
                  or_hi = exp(cte["z_comp",1] + 1.96*cte["z_comp",2])), file.path(OUT, "logistic_first_episode.csv"))
# Linear: observed log fold-rise per SD of score at fixed baseline (fold-rise ratio per SD)
ln <- lm(y_comp ~ base_comp + z_comp, d); ctl <- coeftest(ln, vcov = vcovCL(ln, cluster = d$base_id))
lin <- data.table(term = "z_comp", ratio_per_sd = exp(ctl["z_comp",1]), lo = exp(ctl["z_comp",1]-1.96*ctl["z_comp",2]),
                  hi = exp(ctl["z_comp",1]+1.96*ctl["z_comp",2]), p = ctl["z_comp",4], predicted_ratio_per_sd = exp(sd_raw))
fwrite(lin, file.path(OUT, "foldrise_per_sd.csv")); print(lin)

# ---- 5. ROC -------------------------------------------------------------------
auc_fun <- function(x, pred) as.numeric(pROC::auc(pROC::roc(x$resp_comp, x[[pred]], quiet = TRUE, direction = "<")))
roc_tab <- rbindlist(lapply(names(subsets), function(s) {
  r <- boot_ci(function(x) { a1 <- auc_fun(x, "pred_M1"); a3 <- auc_fun(x, "pred_M3"); a2 <- auc_fun(x, "pred_M2"); a4 <- auc_fun(x, "pred_M4")
                             c(a1, a3, a3 - a1, a2, a4, a4 - a2) }, subsets[[s]])
  x <- d[subsets[[s]]]
  dl <- pROC::roc.test(pROC::roc(x$resp_comp, x$pred_M1, quiet=TRUE, direction="<"), pROC::roc(x$resp_comp, x$pred_M3, quiet=TRUE, direction="<"), method = "delong", paired = TRUE)$p.value
  dl2 <- pROC::roc.test(pROC::roc(x$resp_comp, x$pred_M2, quiet=TRUE, direction="<"), pROC::roc(x$resp_comp, x$pred_M4, quiet=TRUE, direction="<"), method = "delong", paired = TRUE)$p.value
  cbind(subset = s, quantity = c("AUC titers (M1)","AUC titers + composition (M3)","Delta AUC M1->M3",
                                 "AUC titers + clinical (M2)","AUC titers + clinical + composition (M4)","Delta AUC M2->M4"), r,
        p_delong = c(NA, NA, dl, NA, NA, dl2))
}))
fwrite(roc_tab, file.path(OUT, "roc.csv")); print(roc_tab)
# Per-antigen ROC from the multitask predictions (supplement)
roc_ag <- rbindlist(lapply(V, function(v) { x <- d
  r1 <- pROC::roc(x[[paste0("resp_", v)]], x[[paste0("pM1_", v)]], quiet=TRUE, direction="<")
  r3 <- pROC::roc(x[[paste0("resp_", v)]], x[[paste0("pM3_", v)]], quiet=TRUE, direction="<")
  data.table(antigen = v, auc_m1 = as.numeric(r1$auc), auc_m3 = as.numeric(r3$auc),
             p_delong = pROC::roc.test(r1, r3, method="delong", paired=TRUE)$p.value) }))
fwrite(roc_ag, file.path(OUT, "roc_by_antigen.csv")); print(roc_ag)

# ---- 6. Pooled model within each cohort -----------------------------------------
mse_gain <- function(x, a, b) { ea <- rowSums(sapply(V, function(v) (x[[paste0("p", a, "_", v)]] - x[[paste0("y_", v)]])^2))
  eb <- rowSums(sapply(V, function(v) (x[[paste0("p", b, "_", v)]] - x[[paste0("y_", v)]])^2)); 100 * (sum(ea) - sum(eb)) / sum(ea) }
cohort_gain <- rbindlist(lapply(names(subsets), function(s) cbind(subset = s, comparison = c("M1->M3","M2->M4"),
  boot_ci(function(x) c(mse_gain(x, "M1", "M3"), mse_gain(x, "M2", "M4")), subsets[[s]]))))
fwrite(cohort_gain, file.path(OUT, "pooled_model_gain_by_cohort.csv")); print(cohort_gain)

# ---- 7. Routine counts -----------------------------------------------------------
labs <- rbindlist(lapply(c("alc","cd4","cd19"), function(l) { x <- d[is.finite(get(l)) & get(l) > 0]
  f0 <- glm(resp_comp ~ base_comp, binomial, x); fl <- glm(as.formula(paste("resp_comp ~ base_comp + log(", l, ")")), binomial, x)
  fz <- glm(resp_comp ~ base_comp + z_comp, binomial, x); fb <- glm(as.formula(paste("resp_comp ~ base_comp + z_comp + log(", l, ")")), binomial, x)
  ctb <- coeftest(fb, vcov = vcovCL(fb, cluster = x$base_id))
  data.table(lab = l, n = nrow(x), rho_score = cor(x$z_comp, x[[l]], method = "spearman"),
             rho_outcome = cor(x[[l]], x$y_comp, method = "spearman"),
             partial_rho_outcome = partial_spearman(x[[l]], x$y_comp, x$base_comp),
             auc_baseline = as.numeric(pROC::auc(x$resp_comp, fitted(f0), quiet=TRUE)),
             auc_baseline_lab = as.numeric(pROC::auc(x$resp_comp, fitted(fl), quiet=TRUE)),
             auc_baseline_score = as.numeric(pROC::auc(x$resp_comp, fitted(fz), quiet=TRUE)),
             p_score_adjusted_for_lab = ctb["z_comp", 4], p_lab_adjusted_for_score = ctb[paste0("log(", l, ")"), 4]) }))
fwrite(labs, file.path(OUT, "routine_counts.csv")); print(labs)
cat("Done.\n")
