#!/usr/bin/env Rscript
# E6-01-cellular-score-clinical.R
#
# Translates the multitask elastic-net model (M3) into a single per-episode
# "cellular response score" and reports what that score buys in clinical units:
# observed geometric-mean fold-rise and observed seroconversion rate, stratified
# by pre-vaccination titer quartile.
#
# Score definition
#   score_i = sum_j  CLR_ij * beta_j
# where beta_j are the M3 elastic-net coefficients on the cluster features only.
# The intercept and the three baseline-titer coefficients are deliberately left
# out, so the score carries the cellular signal and nothing else.
#
# Everything is out-of-fold: episodes in fold f are scored with the model that
# was fit without fold f. Higher score = model predicts a stronger response;
# the clinically interesting group is therefore the LOW-score group.
#
# Baseline-titer quartiles are taken verbatim from the stratified analysis that
# backs Table 3 of the manuscript, so the strata match exactly.
#
# Confidence intervals come from a recipient-level bootstrap: unique base
# recipients (episode id with the trailing _yr<n> stripped) are resampled with
# replacement, so repeat vaccinations from the same person move together.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

cat("=== E6-01: cellular response score -> clinical readout ===\n")

PER_FOLD_ROOT <- file.path(Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06"), "per_fold_atlas_5000")
STRAT_FILE <- file.path(Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06"), "post_hoc_5000/stratified/E2_stratified.RData")
OUT_DIR <- Sys.getenv("FLU_MEMO_DIR_E6_01", "/home/amir/repo/flu-cellular-prediction/paper/memo_figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

VIRUSES <- c("h1n1", "h3n2", "bvic")
VIRUS_LABEL <- c(h1n1 = "H1N1", h3n2 = "H3N2", bvic = "B/Victoria")
PRIMARY <- "h3n2"
N_FOLDS <- 5L
B_BOOT <- 2000L
SEED <- 20260728L
SEROCONV_LOG <- log(4) - 1e-12  # numerical tolerance at the fourfold boundary      # outcomes are natural-log fold-rise

QUART_LABEL <- c(Q1_low = "Q1\n(lowest baseline)", Q2 = "Q2", Q3 = "Q3",
                 Q4_high = "Q4\n(highest baseline)")
GROUP_COLORS <- c(High = "#2166ac", Low = "#b2182b")
GROUP_LABELS <- c(High = "High cellular score",
                  Low  = "Low cellular score")

# ---------------------------------------------------------------------------
# 1. Out-of-fold cellular response score
# ---------------------------------------------------------------------------
cat("\n[1] Building out-of-fold cellular response scores...\n")

score_rows <- list()
for (f in seq_len(N_FOLDS)) {
  env <- new.env()
  load(file.path(PER_FOLD_ROOT, sprintf("fold_%d", f), sprintf("fold_%d_all.RData", f)),
       envir = env)

  clr <- env$clr_test
  stopifnot(identical(rownames(clr), env$test_p_ids))
  cluster_cols <- colnames(clr)

  fold_scores <- matrix(NA_real_, nrow = nrow(clr), ncol = length(VIRUSES),
                        dimnames = list(rownames(clr), VIRUSES))
  for (v in VIRUSES) {
    cf <- env$coef_M3_en[[v]]
    beta <- cf[cluster_cols]
    if (any(is.na(beta))) {
      stop("Fold ", f, " / ", v, ": cluster coefficients missing for ",
           paste(cluster_cols[is.na(beta)], collapse = ", "))
    }
    fold_scores[, v] <- as.vector(clr %*% beta)
  }

  score_rows[[f]] <- data.table(
    episode_id = rownames(clr),
    fold = f,
    n_clusters = ncol(clr),
    n_nonzero_h3n2 = sum(env$coef_M3_en[["h3n2"]][cluster_cols] != 0),
    score_h1n1 = fold_scores[, "h1n1"],
    score_h3n2 = fold_scores[, "h3n2"],
    score_bvic = fold_scores[, "bvic"],
    y_h1n1 = env$Y_test[, "h1n1"],
    y_h3n2 = env$Y_test[, "h3n2"],
    y_bvic = env$Y_test[, "bvic"]
  )
  cat(sprintf("  fold %d: n=%d, k=%d, nonzero cluster coefs (H3N2)=%d\n",
              f, nrow(clr), ncol(clr), score_rows[[f]]$n_nonzero_h3n2[1]))
  rm(env)
}
dt <- rbindlist(score_rows)
stopifnot(!anyDuplicated(dt$episode_id))
cat(sprintf("  total episodes scored: %d\n", nrow(dt)))

# Per-fold score distributions. Each fold has its own atlas and its own number of
# clusters, so the score is not guaranteed to sit on a common scale across folds.
# Report the spread so the reader can judge how much fold-to-fold offset there is.
cat("\n  Per-fold H3N2 score distribution (raw, as defined):\n")
print(dt[, .(n = .N, mean = round(mean(score_h3n2), 3), sd = round(sd(score_h3n2), 3),
             min = round(min(score_h3n2), 3), max = round(max(score_h3n2), 3)), by = fold][order(fold)])

# ---------------------------------------------------------------------------
# 2. Baseline-titer quartiles (reused from the Table 3 stratification)
# ---------------------------------------------------------------------------
cat("\n[2] Attaching baseline-titer quartiles...\n")
strat <- new.env()
load(STRAT_FILE, envir = strat)
baseline_vec <- strat$avg_baseline               # named: log geometric-mean HAI across strains
quartile_vec <- as.character(strat$quartile)     # same order as baseline_vec

stopifnot(setequal(names(baseline_vec), dt$episode_id))
match_idx <- match(dt$episode_id, names(baseline_vec))
stopifnot(!anyNA(match_idx))
set(dt, j = "avg_baseline", value = unname(baseline_vec[match_idx]))
set(dt, j = "quartile",
    value = factor(quartile_vec[match_idx], levels = c("Q1_low", "Q2", "Q3", "Q4_high")))
stopifnot(!anyNA(dt$quartile))
print(dt[, .N, by = quartile][order(quartile)])

# Base recipient id: strip the _yr<n> episode suffix.
dt[, base_id := sub("_yr[0-9]+$", "", episode_id)]
cat(sprintf("  unique base recipients: %d (episodes: %d)\n",
            uniqueN(dt$base_id), nrow(dt)))

# ---------------------------------------------------------------------------
# 3. Within-quartile median split on the cellular score
# ---------------------------------------------------------------------------
cat("\n[3] Splitting each quartile at its median cellular score...\n")
for (v in VIRUSES) {
  sc <- paste0("score_", v)
  gc_col <- paste0("group_", v)
  dt[, (gc_col) := {
    med <- median(get(sc))
    factor(ifelse(get(sc) > med, "High", "Low"), levels = c("High", "Low"))
  }, by = quartile]
}
print(dcast(dt[, .N, by = .(quartile, group_h3n2)], quartile ~ group_h3n2, value.var = "N"))

# Pooled Q1-Q3 label: keep the within-quartile split (so baseline is still
# controlled for), and additionally record a single unstratified median split
# across Q1-Q3 as a sensitivity check.
idx_q13 <- which(dt$quartile %in% c("Q1_low", "Q2", "Q3"))
for (v in VIRUSES) {
  sc <- dt[[paste0("score_", v)]]
  grp <- rep(NA_character_, nrow(dt))
  med <- median(sc[idx_q13])
  grp[idx_q13] <- ifelse(sc[idx_q13] > med, "High", "Low")
  set(dt, j = paste0("group_pooled_", v), value = grp)
}

# ---------------------------------------------------------------------------
# 4. Point estimates and recipient-level bootstrap CIs
# ---------------------------------------------------------------------------
cat("\n[4] Recipient-level bootstrap (B =", B_BOOT, ", seed =", SEED, ")...\n")

base_levels <- sort(unique(dt$base_id))
rows_by_base <- split(seq_len(nrow(dt)), factor(dt$base_id, levels = base_levels))

set.seed(SEED)
boot_index <- vector("list", B_BOOT)
for (b in seq_len(B_BOOT)) {
  picked <- sample.int(length(base_levels), size = length(base_levels), replace = TRUE)
  boot_index[[b]] <- unlist(rows_by_base[picked], use.names = FALSE)
}

geo_mean_fr <- function(y) exp(mean(y))
seroconv <- function(y) mean(y >= SEROCONV_LOG)

ci_of <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 100) return(c(NA_real_, NA_real_))
  as.numeric(stats::quantile(x, c(0.025, 0.975), names = FALSE))
}

# One stratum = one (analysis set, quartile label) pair; returns a tidy row per
# score group plus the between-group contrast, all with bootstrap CIs.
summarize_stratum <- function(rows_high, rows_low, y, stratum, virus) {
  gm_hi <- geo_mean_fr(y[rows_high]); gm_lo <- geo_mean_fr(y[rows_low])
  sc_hi <- seroconv(y[rows_high]);    sc_lo <- seroconv(y[rows_low])

  bh <- bl <- integer(0)
  b_gm_hi <- b_gm_lo <- b_sc_hi <- b_sc_lo <- numeric(B_BOOT)
  b_ratio <- b_scdiff <- numeric(B_BOOT)
  in_hi <- logical(nrow(dt)); in_hi[rows_high] <- TRUE
  in_lo <- logical(nrow(dt)); in_lo[rows_low] <- TRUE

  for (b in seq_len(B_BOOT)) {
    ib <- boot_index[[b]]
    bh <- ib[in_hi[ib]]
    bl <- ib[in_lo[ib]]
    if (length(bh) < 2 || length(bl) < 2) {
      b_gm_hi[b] <- b_gm_lo[b] <- b_sc_hi[b] <- b_sc_lo[b] <- NA_real_
      b_ratio[b] <- b_scdiff[b] <- NA_real_
      next
    }
    b_gm_hi[b] <- geo_mean_fr(y[bh]); b_gm_lo[b] <- geo_mean_fr(y[bl])
    b_sc_hi[b] <- seroconv(y[bh]);    b_sc_lo[b] <- seroconv(y[bl])
    b_ratio[b] <- b_gm_hi[b] / b_gm_lo[b]
    b_scdiff[b] <- b_sc_hi[b] - b_sc_lo[b]
  }

  ci_gm_hi <- ci_of(b_gm_hi); ci_gm_lo <- ci_of(b_gm_lo)
  ci_sc_hi <- ci_of(b_sc_hi); ci_sc_lo <- ci_of(b_sc_lo)
  ci_ratio <- ci_of(b_ratio); ci_scdiff <- ci_of(b_scdiff)

  out <- rbind(
    data.table(virus = virus, stratum = stratum, group = "High", n = length(rows_high),
               gmfr = gm_hi, gmfr_lo = ci_gm_hi[1], gmfr_hi = ci_gm_hi[2],
               seroconv = sc_hi, seroconv_lo = ci_sc_hi[1], seroconv_hi = ci_sc_hi[2]),
    data.table(virus = virus, stratum = stratum, group = "Low", n = length(rows_low),
               gmfr = gm_lo, gmfr_lo = ci_gm_lo[1], gmfr_hi = ci_gm_lo[2],
               seroconv = sc_lo, seroconv_lo = ci_sc_lo[1], seroconv_hi = ci_sc_lo[2])
  )
  contrast <- data.table(
    virus = virus, stratum = stratum,
    gmfr_ratio = gm_hi / gm_lo, gmfr_ratio_lo = ci_ratio[1], gmfr_ratio_hi = ci_ratio[2],
    seroconv_diff_pp = 100 * (sc_hi - sc_lo),
    seroconv_diff_lo = 100 * ci_scdiff[1], seroconv_diff_hi = 100 * ci_scdiff[2],
    fraction_boot_gmfr_ratio_gt1 = mean(b_ratio > 1, na.rm = TRUE),
    fraction_boot_seroconv_diff_gt0 = mean(b_scdiff > 0, na.rm = TRUE)
  )
  list(groups = out, contrast = contrast)
}

group_tab <- list(); contrast_tab <- list()
for (v in VIRUSES) {
  y <- dt[[paste0("y_", v)]]
  gcol <- dt[[paste0("group_", v)]]
  pcol <- dt[[paste0("group_pooled_", v)]]

  for (q in levels(dt$quartile)) {
    res <- summarize_stratum(
      rows_high = which(dt$quartile == q & gcol == "High"),
      rows_low  = which(dt$quartile == q & gcol == "Low"),
      y = y, stratum = q, virus = v
    )
    group_tab[[length(group_tab) + 1L]] <- res$groups
    contrast_tab[[length(contrast_tab) + 1L]] <- res$contrast
  }

  # Pooled Q1-Q3, keeping the within-quartile median split.
  in_q13 <- dt$quartile %in% c("Q1_low", "Q2", "Q3")
  res <- summarize_stratum(
    rows_high = which(in_q13 & gcol == "High"),
    rows_low  = which(in_q13 & gcol == "Low"),
    y = y, stratum = "Q1-Q3 pooled (within-quartile split)", virus = v
  )
  group_tab[[length(group_tab) + 1L]] <- res$groups
  contrast_tab[[length(contrast_tab) + 1L]] <- res$contrast

  # Pooled Q1-Q3, single median split across the pooled set (sensitivity).
  res <- summarize_stratum(
    rows_high = which(in_q13 & pcol == "High"),
    rows_low  = which(in_q13 & pcol == "Low"),
    y = y, stratum = "Q1-Q3 pooled (single median split)", virus = v
  )
  group_tab[[length(group_tab) + 1L]] <- res$groups
  contrast_tab[[length(contrast_tab) + 1L]] <- res$contrast

  # All four quartiles pooled, within-quartile split.
  res <- summarize_stratum(
    rows_high = which(gcol == "High"),
    rows_low  = which(gcol == "Low"),
    y = y, stratum = "All episodes (within-quartile split)", virus = v
  )
  group_tab[[length(group_tab) + 1L]] <- res$groups
  contrast_tab[[length(contrast_tab) + 1L]] <- res$contrast
}

group_tab <- rbindlist(group_tab)
contrast_tab <- rbindlist(contrast_tab)

fwrite(group_tab, file.path(OUT_DIR, "memo_tab1_score_groups.csv"))
fwrite(contrast_tab, file.path(OUT_DIR, "memo_tab1_score_contrasts.csv"))

cat("\n--- H3N2 group summaries ---\n")
print(group_tab[virus == "h3n2"], nrows = 100)
cat("\n--- H3N2 contrasts ---\n")
print(contrast_tab[virus == "h3n2"], nrows = 100)

# ---------------------------------------------------------------------------
# 5. Continuous association between score and outcome
# ---------------------------------------------------------------------------
cat("\n[5] Continuous score-outcome association...\n")
corr_rows <- list()
for (v in VIRUSES) {
  sc <- dt[[paste0("score_", v)]]; y <- dt[[paste0("y_", v)]]
  corr_rows[[length(corr_rows) + 1L]] <- data.table(
    virus = v, stratum = "All", n = length(y),
    pearson = cor(sc, y), spearman = cor(sc, y, method = "spearman"))
  for (q in levels(dt$quartile)) {
    i <- which(dt$quartile == q)
    corr_rows[[length(corr_rows) + 1L]] <- data.table(
      virus = v, stratum = q, n = length(i),
      pearson = cor(sc[i], y[i]), spearman = cor(sc[i], y[i], method = "spearman"))
  }
}
score_corr <- rbindlist(corr_rows)
fwrite(score_corr, file.path(OUT_DIR, "memo_tab1_score_correlation.csv"))
print(score_corr[virus == "h3n2"])

# ---------------------------------------------------------------------------
# 6. Figures
# ---------------------------------------------------------------------------
cat("\n[6] Drawing figures...\n")

plot_dt <- group_tab[stratum %in% levels(dt$quartile)]
plot_dt[, quartile := factor(stratum, levels = levels(dt$quartile))]
plot_dt[, group := factor(group, levels = c("High", "Low"))]

base_theme <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())

dodge <- position_dodge(width = 0.72)

# Fold-rise spans two orders of magnitude across quartiles, so the bars are drawn
# on a log2 axis with the bar base at a fold-rise of 1 (no change) and the axis
# relabeled in fold-rise units.
gm_breaks <- log2(c(0.5, 1, 2, 4, 8, 16, 32, 64, 128))
gm_labels <- c("0.5", "1", "2", "4", "8", "16", "32", "64", "128")
for (cc in c("gmfr", "gmfr_lo", "gmfr_hi")) {
  plot_dt[[paste0("l2_", cc)]] <- log2(plot_dt[[cc]])
}
gm_lim <- function(d) c(min(d$l2_gmfr_lo, 0) - 0.25, max(d$l2_gmfr_hi) + 0.25)

p1a <- ggplot(plot_dt[virus == PRIMARY],
              aes(x = quartile, y = l2_gmfr, fill = group)) +
  geom_col(position = dodge, width = 0.66) +
  geom_errorbar(aes(ymin = l2_gmfr_lo, ymax = l2_gmfr_hi), position = dodge,
                width = 0.16, linewidth = 0.5) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.4) +
  geom_hline(yintercept = log2(4), linetype = "dashed", color = "grey35", linewidth = 0.45) +
  annotate("text", x = 4.35, y = log2(4) + 0.22, label = "4-fold rise", hjust = 1,
           size = 3.1, color = "grey35") +
  scale_fill_manual(values = GROUP_COLORS, labels = GROUP_LABELS) +
  scale_x_discrete(labels = QUART_LABEL) +
  scale_y_continuous(breaks = gm_breaks, labels = gm_labels,
                     limits = gm_lim(plot_dt[virus == PRIMARY]),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Pre-vaccination titer quartile",
       y = expression("Geometric mean fold-rise (" * log[2] * " scale)")) +
  base_theme

p1b <- ggplot(plot_dt[virus == PRIMARY],
              aes(x = quartile, y = 100 * seroconv, fill = group)) +
  geom_col(position = dodge, width = 0.66) +
  geom_errorbar(aes(ymin = 100 * seroconv_lo, ymax = 100 * seroconv_hi),
                position = dodge, width = 0.16, linewidth = 0.5) +
  scale_fill_manual(values = GROUP_COLORS, labels = GROUP_LABELS) +
  scale_x_discrete(labels = QUART_LABEL) +
  labs(x = "Pre-vaccination titer quartile",
       y = "Seroconversion rate (%)") +
  base_theme

fig1 <- (p1a | p1b) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = list(c("A.", "B."))) &
  theme(plot.tag = element_text(face = "bold", size = 13), legend.position = "bottom")

ggsave(file.path(OUT_DIR, "memo_fig1_score_by_quartile.png"), fig1,
       width = 10.5, height = 5.2, dpi = 300)
cat("  wrote memo_fig1_score_by_quartile.png\n")

# --- Figure 2: continuous score vs observed fold-rise -----------------------
fr_breaks <- log(c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256))
fr_labels <- c("0.125", "0.25", "0.5", "1", "2", "4", "8", "16", "32", "64", "128", "256")
QCOL <- c(Q1_low = "#2166ac", Q2 = "#67a9cf", Q3 = "#ef8a62", Q4_high = "#b2182b")
QLAB <- c(Q1_low = "Q1 (lowest baseline)", Q2 = "Q2", Q3 = "Q3",
          Q4_high = "Q4 (highest baseline)")

sc2 <- data.table(score = dt[[paste0("score_", PRIMARY)]],
                  y = dt[[paste0("y_", PRIMARY)]],
                  quartile = dt$quartile)

fig2 <- ggplot(sc2, aes(x = score, y = y)) +
  geom_hline(yintercept = SEROCONV_LOG, linetype = "dashed",
             color = "grey35", linewidth = 0.4) +
  geom_point(aes(color = quartile), alpha = 0.55, size = 1.9) +
  geom_smooth(aes(color = quartile), method = "lm", formula = y ~ x,
              se = FALSE, linewidth = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              color = "black", linetype = "solid", linewidth = 1.1, fill = "grey70") +
  annotate("text", x = min(sc2$score), y = SEROCONV_LOG + 0.22,
           label = "Seroconversion threshold", hjust = 0, size = 3.1, color = "grey35") +
  scale_color_manual(values = QCOL, labels = QLAB, name = "Baseline quartile") +
  scale_y_continuous(breaks = fr_breaks, labels = fr_labels) +
  labs(x = "Cellular response score (out-of-fold)",
       y = "Observed fold-rise") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

ggsave(file.path(OUT_DIR, "memo_fig2_score_continuous.png"), fig2,
       width = 8.2, height = 6.0, dpi = 300)
cat("  wrote memo_fig2_score_continuous.png\n")

# --- Figure 3: all three antigens ------------------------------------------
plot_dt[, virus_lab := factor(VIRUS_LABEL[virus], levels = VIRUS_LABEL)]

fig3 <- ggplot(plot_dt, aes(x = quartile, y = l2_gmfr, fill = group)) +
  geom_col(position = dodge, width = 0.66) +
  geom_errorbar(aes(ymin = l2_gmfr_lo, ymax = l2_gmfr_hi), position = dodge,
                width = 0.16, linewidth = 0.45) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.4) +
  geom_hline(yintercept = log2(4), linetype = "dashed", color = "grey35", linewidth = 0.45) +
  facet_wrap(~ virus_lab, nrow = 1) +
  scale_fill_manual(values = GROUP_COLORS, labels = GROUP_LABELS) +
  scale_x_discrete(labels = c(Q1_low = "Q1", Q2 = "Q2", Q3 = "Q3", Q4_high = "Q4")) +
  scale_y_continuous(breaks = gm_breaks, labels = gm_labels,
                     limits = gm_lim(plot_dt), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Pre-vaccination titer quartile",
       y = expression("Geometric mean fold-rise (" * log[2] * " scale)")) +
  base_theme +
  theme(strip.text = element_text(face = "bold", size = 11))

ggsave(file.path(OUT_DIR, "memo_fig3_all_antigens.png"), fig3,
       width = 11.5, height = 4.6, dpi = 300)
cat("  wrote memo_fig3_all_antigens.png\n")

# ---------------------------------------------------------------------------
# 7. Sensitivity: fold-centered score
# ---------------------------------------------------------------------------
# The per-fold atlases differ in size, so the raw score can carry a fold-level
# offset. Repeat the H3N2 median split after centering the score within fold to
# check that the headline contrast is not an artifact of that offset.
cat("\n[7] Sensitivity: fold-centered score (H3N2)...\n")
dt[, score_h3n2_c := score_h3n2 - mean(score_h3n2), by = fold]
dt[, group_h3n2_c := {
  med <- median(score_h3n2_c)
  factor(ifelse(score_h3n2_c > med, "High", "Low"), levels = c("High", "Low"))
}, by = quartile]

sens_rows <- list()
for (q in levels(dt$quartile)) {
  res <- summarize_stratum(
    rows_high = which(dt$quartile == q & dt$group_h3n2_c == "High"),
    rows_low  = which(dt$quartile == q & dt$group_h3n2_c == "Low"),
    y = dt$y_h3n2, stratum = q, virus = "h3n2_fold_centered")
  sens_rows[[length(sens_rows) + 1L]] <- res$groups
}
res <- summarize_stratum(
  rows_high = which(dt$quartile %in% c("Q1_low", "Q2", "Q3") & dt$group_h3n2_c == "High"),
  rows_low  = which(dt$quartile %in% c("Q1_low", "Q2", "Q3") & dt$group_h3n2_c == "Low"),
  y = dt$y_h3n2, stratum = "Q1-Q3 pooled", virus = "h3n2_fold_centered")
sens_rows[[length(sens_rows) + 1L]] <- res$groups
sens_tab <- rbindlist(sens_rows)
fwrite(sens_tab, file.path(OUT_DIR, "memo_tab1_sensitivity_fold_centered.csv"))
print(sens_tab)

saveRDS(dt, file.path(OUT_DIR, "memo_episode_scores.rds"))
cat("\nDone. Outputs in", OUT_DIR, "\n")
