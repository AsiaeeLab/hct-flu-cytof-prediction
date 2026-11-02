#!/usr/bin/env Rscript
# E15-02: Number-of-populations selection rules, from the E15-01 sweep.
#
# Rules compared, each applied inside every outer training fold (inner CV only,
# so the outer estimate stays nested):
#   min      resolution with the lowest mean inner criterion (E14-01 rule);
#   one_se   coarsest resolution whose mean inner criterion is within one
#            standard error (SD over the three inner folds / sqrt(3)) of the
#            minimum (Hastie, Tibshirani and Friedman's one-standard-error rule,
#            with fewer populations as the simpler model).
# Each rule is applied to the multitask criterion (primary) and to the
# composite-outcome criterion. For each rule the outer held-out gains
# (M1->M3, M2->M4; three antigens and composite) come with recipient-bootstrap
# 95% intervals and two-sided p (B = 2,000, seed 20260728).
# Also written: the shared-atlas version of each rule (criterion averaged over
# the five folds, SE over folds) and the held-out gain at every fixed
# resolution (diagnostic only; not used to choose anything).
suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/first_episode_2026-10-14")
OUT <- file.path(ROOT, "post_hoc_5000/resolution_sweep")
inner <- fread(file.path(OUT, "inner_criteria_by_fold.csv"))
o <- readRDS(file.path(OUT, "outer_predictions_by_resolution.rds"))
outer <- o$outer; ref <- o$ref; y <- o$y
V <- c("h1n1", "h3n2", "bvic")

# ---- Per-fold selection -------------------------------------------------------
pick <- function(dt, crit) {
  s <- dt[, .(m = mean(get(crit)), se = sd(get(crit)) / sqrt(.N), k = k[1]), by = .(fold, res)]
  s[, by = fold, `:=`(best = m == min(m))]
  thr <- s[best == TRUE, .(fold, thr = m + se)]
  s <- merge(s, thr, by = "fold")
  rbind(s[best == TRUE, .(fold, res, k, rule = "min")][, .SD[1], by = fold],
        s[m <= thr][order(res)][, .SD[1], by = fold][, .(fold, res, k, rule = "one_se")])
}
sel <- rbind(cbind(pick(inner, "crit_multi"), criterion = "multitask"),
             cbind(pick(inner, "crit_comp"), criterion = "composite"))
setorder(sel, criterion, rule, fold)
print(dcast(sel, criterion + rule ~ fold, value.var = "k"))
fwrite(sel, file.path(OUT, "selected_k_by_rule.csv"))

# ---- Outer gains for a set of per-fold resolutions ----------------------------
mse_table <- function(pred) {
  d <- merge(merge(pred, ref, by = c("fold", "episode_id")), y, by = "episode_id")
  se <- function(m, a) if (a == "comp") (d[[paste0(m, "_comp")]] - d$y_comp)^2 else
    rowMeans(sapply(V, function(v) (d[[paste0(m, "_", v)]] - d[[paste0("y_", v)]])^2))
  data.table(episode_id = d$episode_id, M1 = se("M1", "multi"), M2 = se("M2", "multi"),
             M3 = se("M3", "multi"), M4 = se("M4", "multi"),
             cM1 = se("M1", "comp"), cM2 = se("M2", "comp"), cM3 = se("M3", "comp"), cM4 = se("M4", "comp"))
}
B <- 2000L; set.seed(20260728)
N_EP <- length(unique(outer$episode_id))
boot_idx <- replicate(B, sample.int(N_EP, replace = TRUE))   # one episode per recipient
gain <- function(a, b, idx) 100 * (1 - sum(b[idx]) / sum(a[idx]))
summ <- function(e) {
  pairs <- list(c("M1", "M3"), c("M2", "M4"), c("cM1", "cM3"), c("cM2", "cM4"))
  rbindlist(lapply(pairs, function(p) {
    a <- e[[p[1]]]; b <- e[[p[2]]]; est <- gain(a, b, seq_along(a))
    bs <- apply(boot_idx, 2, function(i) gain(a, b, i))
    data.table(comparison = paste0(p[1], "->", p[2]), gain = est, lo = quantile(bs, .025), hi = quantile(bs, .975),
               p = min(1, 2 * min(mean(bs <= 0), mean(bs >= 0))))
  }))
}
rule_results <- rbindlist(lapply(split(sel, by = c("criterion", "rule")), function(s) {
  pred <- merge(outer, s[, .(fold, res)], by = c("fold", "res"))
  stopifnot(nrow(pred) == N_EP)
  cbind(criterion = s$criterion[1], rule = s$rule[1], k_by_fold = paste(s$k, collapse = "/"), summ(mse_table(pred)))
}))
print(rule_results[, .(criterion, rule, k_by_fold, comparison, gain = round(gain, 1), lo = round(lo, 1), hi = round(hi, 1), p = round(p, 3))])
fwrite(rule_results, file.path(OUT, "gains_by_selection_rule.csv"))

# ---- Shared-atlas version of each rule (criterion averaged over folds) --------
shared <- rbindlist(lapply(c("crit_multi", "crit_comp"), function(cr) {
  pf <- inner[, .(m = mean(get(cr)), k = k[1]), by = .(fold, res)]
  a <- pf[, .(mean = mean(m), se = sd(m) / sqrt(.N), median_k = median(k)), by = res][order(res)]
  b <- a[which.min(mean)]
  data.table(criterion = cr, rule = c("min", "one_se"),
             res = c(b$res, a[mean <= b$mean + b$se][order(res)][1, res]),
             median_k = c(b$median_k, a[mean <= b$mean + b$se][order(res)][1, median_k]))
}))
print(shared); fwrite(shared, file.path(OUT, "shared_atlas_k_by_rule.csv"))

# ---- Held-out gain at every fixed resolution (diagnostic) ---------------------
fixed <- rbindlist(lapply(sort(unique(outer$res)), function(r) {
  pred <- outer[res == r]
  e <- mse_table(pred)
  data.table(res = r, median_k = median(pred[, k[1], by = fold]$V1),
             M1_M3 = gain(e$M1, e$M3, seq_len(nrow(e))), M2_M4 = gain(e$M2, e$M4, seq_len(nrow(e))),
             cM1_cM3 = gain(e$cM1, e$cM3, seq_len(nrow(e))), cM2_cM4 = gain(e$cM2, e$cM4, seq_len(nrow(e))))
}))
print(fixed[, lapply(.SD, function(x) round(x, 1))])
fwrite(fixed, file.path(OUT, "heldout_gain_by_fixed_resolution.csv"))
cat("Done.\n")
