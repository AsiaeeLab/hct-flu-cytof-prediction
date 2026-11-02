#!/usr/bin/env Rscript

FLU_ROOT <- Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06")
FLU_PROD_ATLAS <- Sys.getenv("FLU_PROD_ATLAS", "/home/amir/repo/flu-cellular-prediction/results/timing_corrected_2026-09-06/final_atlas_timing_corrected.RData")
suppressPackageStartupMessages({
  library(here)
  library(ggplot2)
})

log_msg <- function(fmt, ...) {
  msg <- sprintf(fmt, ...)
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
  flush.console()
}

reorder_within <- function(x, by, within, fun = mean, sep = "___") {
  stats::reorder(paste(x, within, sep = sep), by, FUN = fun)
}

scale_x_reordered <- function(..., sep = "___") {
  scale_x_discrete(labels = function(x) gsub(paste0(sep, ".*$"), "", x), ...)
}

compute_bootstrap_delta <- function(pred_a,
                                    pred_b,
                                    y_true,
                                    B = 2000L,
                                    seed = 42L,
                                    virus_idx = NULL) {
  err_a <- (pred_a - y_true)^2
  err_b <- (pred_b - y_true)^2

  if (!is.null(virus_idx)) {
    err_a <- err_a[, virus_idx, drop = FALSE]
    err_b <- err_b[, virus_idx, drop = FALSE]
  }

  delta <- err_a - err_b
  observed_delta <- mean(delta)
  mse_a <- mean(err_a)

  n <- nrow(delta)
  set.seed(seed)
  boot <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    boot[b] <- mean(delta[idx, , drop = FALSE])
  }

  ci <- as.numeric(stats::quantile(boot, probs = c(0.025, 0.975), names = FALSE, na.rm = TRUE))

  list(
    delta_mse = observed_delta,
    delta_pct = 100 * observed_delta / mse_a,
    ci_low = ci[1],
    ci_high = ci[2],
    prob_delta_gt0 = mean(boot > 0),
    p_boot = mean(boot <= 0),
    boot = boot,
    mse_a = mse_a
  )
}

safe_quartile <- function(x) {
  q <- stats::quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  if (any(diff(q) == 0)) {
    r <- rank(x, ties.method = "average")
    q_r <- stats::quantile(r, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
    return(cut(r,
               breaks = q_r,
               include.lowest = TRUE,
               labels = c("Q1_low", "Q2", "Q3", "Q4_high")))
  }
  cut(x,
      breaks = q,
      include.lowest = TRUE,
      labels = c("Q1_low", "Q2", "Q3", "Q4_high"))
}

# ============================================================================
# Paths and inputs
# ============================================================================
OUTPUT_ROOT <- FLU_ROOT
AGG_PATH <- file.path(OUTPUT_ROOT, "per_fold_atlas_5000", "aggregated", "E2_aggregated.RData")
POSTHOC_ROOT <- file.path(OUTPUT_ROOT, "post_hoc_5000")

BOOT_DIR <- file.path(POSTHOC_ROOT, "bootstrap_ci")
STRAT_DIR <- file.path(POSTHOC_ROOT, "stratified")
CTIME_DIR <- file.path(POSTHOC_ROOT, "cluster_time")
FIG_DIR <- file.path(POSTHOC_ROOT, "figures")

for (d in c(POSTHOC_ROOT, BOOT_DIR, STRAT_DIR, CTIME_DIR, FIG_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(AGG_PATH)) {
  stop("Missing aggregated input file: ", AGG_PATH)
}

log_msg("Loading aggregated OOF results from %s", AGG_PATH)
load(AGG_PATH)

if (!exists("oof_pred_M1_en") || !exists("oof_pred_M4_en") || !exists("oof_Y")) {
  stop("Aggregated RData is missing required objects.")
}

VIRUSES <- colnames(oof_Y)
if (is.null(VIRUSES)) {
  VIRUSES <- c("h1n1", "h3n2", "bvic")
  colnames(oof_Y) <- VIRUSES
}

model_preds <- list(
  M1 = oof_pred_M1_en,
  M2 = oof_pred_M2_en,
  M3 = oof_pred_M3_en,
  M4 = oof_pred_M4_en
)

# ============================================================================
# 2A: Paired bootstrap CIs
# ============================================================================
log_msg("Task 2A: Paired bootstrap CIs...")

comparisons <- data.frame(
  comparison = c("M1_vs_M3", "M1_vs_M2", "M2_vs_M4", "M1_vs_M4", "M3_vs_M4", "M2_vs_M3"),
  model_a = c("M1", "M1", "M2", "M1", "M3", "M2"),
  model_b = c("M3", "M2", "M4", "M4", "M4", "M3"),
  stringsAsFactors = FALSE
)

B <- 2000L
bootstrap_results_overall <- list()
bootstrap_results_per_virus <- list()
bootstrap_draws <- list()

for (i in seq_len(nrow(comparisons))) {
  cmp <- comparisons[i, ]
  pred_a <- model_preds[[cmp$model_a]]
  pred_b <- model_preds[[cmp$model_b]]

  res_overall <- compute_bootstrap_delta(
    pred_a = pred_a,
    pred_b = pred_b,
    y_true = oof_Y,
    B = B,
    seed = 42L + i
  )

  bootstrap_results_overall[[i]] <- data.frame(
    comparison = cmp$comparison,
    model_a = cmp$model_a,
    model_b = cmp$model_b,
    virus = "overall",
    delta_mse = res_overall$delta_mse,
    delta_pct = res_overall$delta_pct,
    ci_low = res_overall$ci_low,
    ci_high = res_overall$ci_high,
    prob_delta_gt0 = res_overall$prob_delta_gt0,
    p_boot = res_overall$p_boot,
    stringsAsFactors = FALSE
  )

  bootstrap_draws[[cmp$comparison]] <- res_overall$boot

  for (v_idx in seq_along(VIRUSES)) {
    res_v <- compute_bootstrap_delta(
      pred_a = pred_a,
      pred_b = pred_b,
      y_true = oof_Y,
      B = B,
      seed = 4200L + i * 100L + v_idx,
      virus_idx = v_idx
    )

    bootstrap_results_per_virus[[length(bootstrap_results_per_virus) + 1L]] <- data.frame(
      comparison = cmp$comparison,
      model_a = cmp$model_a,
      model_b = cmp$model_b,
      virus = VIRUSES[v_idx],
      delta_mse = res_v$delta_mse,
      delta_pct = res_v$delta_pct,
      ci_low = res_v$ci_low,
      ci_high = res_v$ci_high,
      prob_delta_gt0 = res_v$prob_delta_gt0,
      p_boot = res_v$p_boot,
      stringsAsFactors = FALSE
    )
  }
}

bootstrap_results_overall <- do.call(rbind, bootstrap_results_overall)
bootstrap_results_per_virus <- do.call(rbind, bootstrap_results_per_virus)
bootstrap_results_all <- rbind(bootstrap_results_overall, bootstrap_results_per_virus)

boot_rdata <- file.path(BOOT_DIR, "E2_bootstrap_ci.RData")
boot_csv <- file.path(BOOT_DIR, "E2_bootstrap_ci.csv")

save(
  B,
  comparisons,
  bootstrap_results_overall,
  bootstrap_results_per_virus,
  bootstrap_results_all,
  bootstrap_draws,
  file = boot_rdata
)
write.csv(bootstrap_results_all, boot_csv, row.names = FALSE)

log_msg("Saved: %s", boot_rdata)
log_msg("Saved: %s", boot_csv)

# RUN_FULL_POSTHOC gates the raw-data / production-atlas sections (2B stratified,
# 2C cluster-time, 2D feature importance). The model-ladder figure (2F) needs
# only the OOF predictions and the 2A bootstrap results, so a figure-only refresh
# can skip the heavy sections. Set RUN_FULL_POSTHOC=TRUE to reproduce every output.
RUN_FULL_POSTHOC <- as.logical(Sys.getenv("RUN_FULL_POSTHOC", "FALSE"))

if (RUN_FULL_POSTHOC) {

# ============================================================================
# 2B: Stratified analysis by baseline titer quartile
# ============================================================================
log_msg("Task 2B: Stratified analysis by baseline titer quartile...")

source(here("methods", "00-paths.R"))
source(here("methods", "08-multi-outcome-covariates.R"))

outcome_obj <- read.csv(
  file.path(paths$raw, "transplant_cytof_counts_hai_12MAY2024.csv"),
  stringsAsFactors = FALSE
)
meta_obj <- read.csv(
  file.path(paths$raw, "hct_limited-dataset_2025-07-09.csv"),
  stringsAsFactors = FALSE
)

outcomes <- setNames(vector("list", length(VIRUSES)), VIRUSES)
for (v in VIRUSES) {
  outcomes[[v]] <- prepare_outcome_with_covariates(
    outcome_obj = outcome_obj,
    meta_obj = meta_obj,
    outcome_col_1 = paste0("geo_mean_", v, "_v1"),
    outcome_col_2 = paste0("geo_mean_", v, "_v3"),
    p_id = common_p_ids,
    include_dose = TRUE,
    include_age = TRUE,
    include_sex = TRUE,
    include_time = TRUE
  )
}

avg_baseline <- rowMeans(cbind(
  outcomes[["h1n1"]]$x[common_p_ids],
  outcomes[["h3n2"]]$x[common_p_ids],
  outcomes[["bvic"]]$x[common_p_ids]
))

quartile <- safe_quartile(avg_baseline)
quartile <- factor(quartile, levels = c("Q1_low", "Q2", "Q3", "Q4_high"))

mse_in_subset <- function(pred, truth, idx) {
  mean((as.vector(pred[idx, , drop = FALSE]) - as.vector(truth[idx, , drop = FALSE]))^2)
}

strat_rows <- list()
for (q in levels(quartile)) {
  idx <- which(quartile == q)
  for (m in names(model_preds)) {
    strat_rows[[length(strat_rows) + 1L]] <- data.frame(
      quartile = q,
      model = m,
      mse = mse_in_subset(model_preds[[m]], oof_Y, idx),
      n_patients = length(idx),
      stringsAsFactors = FALSE
    )
  }
}
stratified_mse <- do.call(rbind, strat_rows)

stratified_improvement <- do.call(
  rbind,
  lapply(levels(quartile), function(q) {
    mse_m1 <- stratified_mse$mse[stratified_mse$quartile == q & stratified_mse$model == "M1"]
    mse_m3 <- stratified_mse$mse[stratified_mse$quartile == q & stratified_mse$model == "M3"]
    mse_m4 <- stratified_mse$mse[stratified_mse$quartile == q & stratified_mse$model == "M4"]
    data.frame(
      quartile = q,
      improvement_M3_vs_M1_pct = 100 * (mse_m1 - mse_m3) / mse_m1,
      improvement_M4_vs_M1_pct = 100 * (mse_m1 - mse_m4) / mse_m1,
      stringsAsFactors = FALSE
    )
  })
)

strat_rdata <- file.path(STRAT_DIR, "E2_stratified.RData")
strat_csv <- file.path(STRAT_DIR, "E2_stratified.csv")

save(
  avg_baseline,
  quartile,
  stratified_mse,
  stratified_improvement,
  file = strat_rdata
)
write.csv(
  merge(stratified_mse, stratified_improvement, by = "quartile", all.x = TRUE),
  strat_csv,
  row.names = FALSE
)

log_msg("Saved: %s", strat_rdata)
log_msg("Saved: %s", strat_csv)

# ============================================================================
# 2C: Cluster -> time prediction (production atlas)
# ============================================================================
log_msg("Task 2C: Cluster -> time prediction...")

prod_atlas_path <- FLU_PROD_ATLAS
if (!file.exists(prod_atlas_path)) {
  stop("Missing production atlas file: ", prod_atlas_path)
}
load(prod_atlas_path)

if (!exists("clr_features") || !exists("X_base")) {
  stop("Production atlas RData must contain clr_features and X_base.")
}

if (!("time_post_transplant" %in% colnames(X_base))) {
  stop("X_base is missing 'time_post_transplant' column.")
}

time_col <- as.numeric(X_base[, "time_post_transplant"])
fit_time <- lm(time_col ~ ., data = as.data.frame(clr_features))
fit_time_summary <- summary(fit_time)

coef_tab <- as.data.frame(coef(fit_time_summary))
coef_tab$feature <- rownames(coef_tab)
coef_tab <- coef_tab[coef_tab$feature != "(Intercept)", , drop = FALSE]
coef_tab$abs_t <- abs(coef_tab[, "t value"])
coef_tab <- coef_tab[order(coef_tab$abs_t, decreasing = TRUE), , drop = FALSE]
top5_time_predictors <- head(coef_tab, 5)

spearman_rho <- vapply(seq_len(ncol(clr_features)), function(j) {
  cor(clr_features[, j], time_col, method = "spearman")
}, numeric(1))

spearman_p <- vapply(seq_len(ncol(clr_features)), function(j) {
  suppressWarnings(cor.test(clr_features[, j], time_col, method = "spearman", exact = FALSE)$p.value)
}, numeric(1))

cluster_time_df <- data.frame(
  cluster = colnames(clr_features),
  spearman_rho = spearman_rho,
  spearman_p = spearman_p,
  stringsAsFactors = FALSE
)

coef_lookup <- setNames(coef_tab[, "Estimate"], coef_tab$feature)
t_lookup <- setNames(coef_tab[, "t value"], coef_tab$feature)
cluster_time_df$lm_coef <- unname(coef_lookup[cluster_time_df$cluster])
cluster_time_df$lm_t <- unname(t_lookup[cluster_time_df$cluster])
cluster_time_df$abs_lm_t <- abs(cluster_time_df$lm_t)
cluster_time_df$is_top5_abs_t <- cluster_time_df$cluster %in% top5_time_predictors$feature
cluster_time_df$model_r_squared <- fit_time_summary$r.squared
cluster_time_df$model_adj_r_squared <- fit_time_summary$adj.r.squared
cluster_time_df <- cluster_time_df[order(abs(cluster_time_df$spearman_rho), decreasing = TRUE), , drop = FALSE]

ctime_rdata <- file.path(CTIME_DIR, "E2_cluster_time.RData")
ctime_csv <- file.path(CTIME_DIR, "E2_cluster_time.csv")

save(
  fit_time,
  fit_time_summary,
  top5_time_predictors,
  cluster_time_df,
  file = ctime_rdata
)
write.csv(cluster_time_df, ctime_csv, row.names = FALSE)

log_msg("Saved: %s", ctime_rdata)
log_msg("Saved: %s", ctime_csv)

# ============================================================================
# 2D: All-feature importance plot from M4 coefficients
# ============================================================================
log_msg("Task 2D: All-feature importance plot...")

fold_paths <- file.path(
  OUTPUT_ROOT,
  "per_fold_atlas_5000",
  paste0("fold_", seq_len(5)),
  paste0("fold_", seq_len(5), "_all.RData")
)
missing_folds <- fold_paths[!file.exists(fold_paths)]
if (length(missing_folds) > 0) {
  stop("Missing fold artifacts needed for importance plot: ", paste(missing_folds, collapse = ", "))
}

importance_rows <- list()
for (i in seq_along(fold_paths)) {
  e <- new.env(parent = emptyenv())
  load(fold_paths[i], envir = e)

  X_M4_train <- get("X_M4_train", envir = e)
  coef_M4_en <- get("coef_M4_en", envir = e)

  feature_sds <- apply(X_M4_train, 2, sd)
  feature_sds[is.na(feature_sds)] <- 0

  if (is.null(names(coef_M4_en)) || !all(VIRUSES %in% names(coef_M4_en))) {
    names(coef_M4_en) <- VIRUSES[seq_along(coef_M4_en)]
  }

  for (v in VIRUSES) {
    coef_v <- coef_M4_en[[v]]
    coef_v <- coef_v[names(coef_v) != "(Intercept)"]
    common_feat <- intersect(names(feature_sds), names(coef_v))
    if (length(common_feat) == 0) {
      next
    }

    imp <- abs(coef_v[common_feat] * feature_sds[common_feat])
    importance_rows[[length(importance_rows) + 1L]] <- data.frame(
      fold = i,
      virus = v,
      feature = common_feat,
      importance = as.numeric(imp),
      stringsAsFactors = FALSE
    )
  }
}

importance_df <- do.call(rbind, importance_rows)
importance_avg <- aggregate(importance ~ virus + feature, data = importance_df, FUN = mean)
importance_avg <- importance_avg[order(importance_avg$virus, importance_avg$importance, decreasing = TRUE), , drop = FALSE]

p_importance <- ggplot(
  importance_avg,
  aes(x = reorder_within(feature, importance, virus), y = importance, fill = virus)
) +
  geom_col(show.legend = FALSE, width = 0.8) +
  coord_flip() +
  facet_wrap(~virus, scales = "free_y", ncol = 1) +
  scale_x_reordered() +
  labs(
    title = "M4 Standardized Feature Importance (|coef x SD|)",
    x = "Feature",
    y = "Average Importance Across Folds"
  ) +
  theme_bw(base_size = 11)

importance_pdf <- file.path(FIG_DIR, "E2_all_feature_importance.pdf")
importance_png <- file.path(FIG_DIR, "E2_all_feature_importance.png")

ggsave(importance_pdf, p_importance, width = 11, height = 12)
ggsave(importance_png, p_importance, width = 11, height = 12, dpi = 300)

log_msg("Saved: %s", importance_pdf)
log_msg("Saved: %s", importance_png)

}  # end RUN_FULL_POSTHOC guard (2B stratified / 2C cluster-time / 2D importance)

# ============================================================================
# 2E: Noise floor recalculation
# ============================================================================
log_msg("Task 2E: Noise floor recalculation...")

sigma2_noise <- log(2)^2

mean_mse_M1 <- mean((as.vector(oof_pred_M1_en) - as.vector(oof_Y))^2)
mean_mse_M4 <- mean((as.vector(oof_pred_M4_en) - as.vector(oof_Y))^2)
var_Y <- stats::var(as.vector(oof_Y))

noise_floor_overall <- data.frame(
  virus = "overall",
  var_Y = var_Y,
  sigma2_noise = sigma2_noise,
  R2_ceiling = 1 - sigma2_noise / var_Y,
  mean_mse_M1 = mean_mse_M1,
  mean_mse_M4 = mean_mse_M4,
  R2_M1 = 1 - mean_mse_M1 / var_Y,
  R2_M4 = 1 - mean_mse_M4 / var_Y,
  stringsAsFactors = FALSE
)
noise_floor_overall$gap_closed <- (noise_floor_overall$R2_M4 - noise_floor_overall$R2_M1) /
  (noise_floor_overall$R2_ceiling - noise_floor_overall$R2_M1)

noise_floor_per_virus <- do.call(
  rbind,
  lapply(VIRUSES, function(v) {
    var_v <- stats::var(oof_Y[, v])
    mse_m1_v <- mean((oof_pred_M1_en[, v] - oof_Y[, v])^2)
    mse_m4_v <- mean((oof_pred_M4_en[, v] - oof_Y[, v])^2)
    r2_ceiling_v <- 1 - sigma2_noise / var_v
    r2_m1_v <- 1 - mse_m1_v / var_v
    r2_m4_v <- 1 - mse_m4_v / var_v

    data.frame(
      virus = v,
      var_Y = var_v,
      sigma2_noise = sigma2_noise,
      R2_ceiling = r2_ceiling_v,
      mean_mse_M1 = mse_m1_v,
      mean_mse_M4 = mse_m4_v,
      R2_M1 = r2_m1_v,
      R2_M4 = r2_m4_v,
      gap_closed = (r2_m4_v - r2_m1_v) / (r2_ceiling_v - r2_m1_v),
      stringsAsFactors = FALSE
    )
  })
)

noise_floor_df <- rbind(noise_floor_overall, noise_floor_per_virus)
noise_floor_csv <- file.path(POSTHOC_ROOT, "E2_noise_floor.csv")
write.csv(noise_floor_df, noise_floor_csv, row.names = FALSE)

log_msg("Saved: %s", noise_floor_csv)

# ============================================================================
# 2F: Model ladder figure (+ bootstrap CI error bars)
# ============================================================================
log_msg("Task 2F: Model ladder figure...")

mse_model_virus <- do.call(
  rbind,
  lapply(c(VIRUSES, "overall"), function(v) {
    data.frame(
      virus = v,
      model = c("M1", "M2", "M3", "M4"),
      mse = c(
        if (v == "overall") mean((as.vector(oof_pred_M1_en) - as.vector(oof_Y))^2) else mean((oof_pred_M1_en[, v] - oof_Y[, v])^2),
        if (v == "overall") mean((as.vector(oof_pred_M2_en) - as.vector(oof_Y))^2) else mean((oof_pred_M2_en[, v] - oof_Y[, v])^2),
        if (v == "overall") mean((as.vector(oof_pred_M3_en) - as.vector(oof_Y))^2) else mean((oof_pred_M3_en[, v] - oof_Y[, v])^2),
        if (v == "overall") mean((as.vector(oof_pred_M4_en) - as.vector(oof_Y))^2) else mean((oof_pred_M4_en[, v] - oof_Y[, v])^2)
      ),
      stringsAsFactors = FALSE
    )
  })
)

boot_ci_lookup <- rbind(bootstrap_results_overall, bootstrap_results_per_virus)
pair_map <- c(M2 = "M1_vs_M2", M3 = "M1_vs_M3", M4 = "M1_vs_M4")

ladder_rows <- list()
for (v in c(VIRUSES, "overall")) {
  mse_m1 <- mse_model_virus$mse[mse_model_virus$virus == v & mse_model_virus$model == "M1"]

  ladder_rows[[length(ladder_rows) + 1L]] <- data.frame(
    virus = v,
    model = "M1",
    mse = mse_m1,
    mse_ci_low = NA_real_,
    mse_ci_high = NA_real_,
    improvement_pct = 0,
    ci_low_pct = NA_real_,
    ci_high_pct = NA_real_,
    stringsAsFactors = FALSE
  )

  for (m in c("M2", "M3", "M4")) {
    row_ci <- boot_ci_lookup[
      boot_ci_lookup$comparison == pair_map[[m]] & boot_ci_lookup$virus == v,
      ,
      drop = FALSE
    ]

    if (nrow(row_ci) != 1L) {
      stop("Could not find unique bootstrap CI row for model ", m, " and virus ", v)
    }

    delta <- row_ci$delta_mse
    ci_low <- row_ci$ci_low
    ci_high <- row_ci$ci_high

    mse_m <- mse_m1 - delta

    ladder_rows[[length(ladder_rows) + 1L]] <- data.frame(
      virus = v,
      model = m,
      mse = mse_m,
      mse_ci_low = mse_m1 - ci_high,
      mse_ci_high = mse_m1 - ci_low,
      improvement_pct = 100 * delta / mse_m1,
      ci_low_pct = 100 * ci_low / mse_m1,
      ci_high_pct = 100 * ci_high / mse_m1,
      stringsAsFactors = FALSE
    )
  }
}

ladder_df <- do.call(rbind, ladder_rows)
ladder_df$model <- factor(ladder_df$model, levels = c("M1", "M2", "M3", "M4"))
ladder_df$virus <- factor(ladder_df$virus, levels = c(VIRUSES, "overall"))

virus_labels <- c(h1n1 = "H1N1", h3n2 = "H3N2", bvic = "B/Victoria", overall = "Overall")
p_ladder <- ggplot(ladder_df, aes(x = model, y = improvement_pct, fill = model)) +
  geom_col(width = 0.75, color = "grey30") +
  geom_errorbar(
    data = ladder_df[ladder_df$model != "M1", , drop = FALSE],
    aes(ymin = ci_low_pct, ymax = ci_high_pct),
    width = 0.20,
    linewidth = 0.5
  ) +
  facet_wrap(~virus, ncol = 2, labeller = labeller(virus = virus_labels)) +
  labs(
    x = "Model",
    y = "MSE improvement over M1 (%)"
  ) +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "none")

ladder_pdf <- file.path(FIG_DIR, "E2_model_ladder.pdf")
ladder_png <- file.path(FIG_DIR, "E2_model_ladder.png")

ggsave(ladder_pdf, p_ladder, width = 10, height = 7)
ggsave(ladder_png, p_ladder, width = 10, height = 7, dpi = 300)
file.copy(ladder_png,
          "/home/amir/repo/flu-cellular-prediction/paper/images/E2_model_ladder.png",
          overwrite = TRUE)

log_msg("Saved: %s", ladder_pdf)
log_msg("Saved: %s", ladder_png)
log_msg("Done.")
