#!/usr/bin/env Rscript
# E2-marker-reduction-perfold.R
# Marker reduction analysis using per-fold atlas (5000 micro-clusters).
# For each fold, uses that fold's atlas to rank markers, reassign cells with
# reduced panels, build new CLR features, and evaluate prediction.

FLU_ROOT <- Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06")
suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(glmnet)
  library(RANN)
  library(ggplot2)
  library(patchwork)
})

cat("=== Marker Reduction (Per-Fold Atlas, 5000 micro) ===\n")

SEED <- 7000L
ALPHA <- 0.9
VIRUSES <- c("h1n1", "h3n2", "bvic")
CELL_PER_PATIENT <- 2000L
CELL_SAMPLE_PRES <- 100000L
BUDGETS <- c(30, 25, 20, 15, 12, 10, 8, 6)
n_folds <- 5

PER_FOLD_ROOT <- file.path(FLU_ROOT, "per_fold_atlas_5000")
OUT_DIR <- file.path(FLU_ROOT, "post_hoc_5000")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(file.path(OUT_DIR, "marker_reduction"), recursive = TRUE, showWarnings = FALSE)

SAVED_RDATA <- file.path(OUT_DIR, "marker_reduction", "E2_marker_reduction_perfold.RData")

if (file.exists(SAVED_RDATA)) {
  cat("Cached marker-reduction results found; loading and skipping recompute.\n")
  load(SAVED_RDATA)   # restores pres_dt, perf_dt, rank_dt, BUDGETS, timings
} else {

# Methods needed only for the raw-data recompute path (some pull in caret).
source(here("methods", "00-paths.R"))
source(here("methods", "08-multi-outcome-covariates.R"))
source(here("methods", "35-per-fold-atlas-cv.R"))

tic_total <- Sys.time()

# --- Load raw cell data ---
cat("Loading raw cell data...\n")
load(file.path(paths$raw, "hsct_transplant_object_harmony_final_small.RData"))
full_data <- obj$batch_correction$data
full_p_id <- obj$source
marker_names <- colnames(full_data)
rm(obj); gc()
cat(sprintf("  Cells: %d x %d markers\n", nrow(full_data), ncol(full_data)))

# --- Load aggregated data for fold assignments and outcomes ---
load(file.path(PER_FOLD_ROOT, "aggregated", "E2_aggregated.RData"))

# --- Prepare outcomes ---
cat("Preparing outcomes...\n")
outcome_obj <- read.csv(
  file.path(paths$raw, "transplant_cytof_counts_hai_12MAY2024.csv"),
  stringsAsFactors = FALSE
)
meta_obj <- read.csv(
  file.path(paths$raw, "hct_limited-dataset_2025-07-09.csv"),
  stringsAsFactors = FALSE
)

outcomes_list <- setNames(vector("list", length(VIRUSES)), VIRUSES)
for (v in VIRUSES) {
  outcomes_list[[v]] <- prepare_outcome_with_covariates(
    outcome_obj = outcome_obj,
    meta_obj = meta_obj,
    outcome_col_1 = paste0("geo_mean_", v, "_v1"),
    outcome_col_2 = paste0("geo_mean_", v, "_v3"),
    p_id = full_p_id,
    include_dose = TRUE,
    include_age = TRUE,
    include_sex = TRUE,
    include_time = TRUE
  )
}

eval_patients <- common_p_ids
eval_patients <- eval_patients[eval_patients %in% unique(full_p_id)]
cat("  Eval patients:", length(eval_patients), "\n")

Y_mat <- oof_Y  # 343 x 3

# Build baseline covariate matrix
base_train_mat <- sapply(VIRUSES, function(v) outcomes_list[[v]]$x[eval_patients])
colnames(base_train_mat) <- paste0("baseline_", VIRUSES)
covariates <- outcomes_list[[VIRUSES[1]]]$covariates[eval_patients, , drop = FALSE]
X_base <- as.matrix(cbind(base_train_mat, covariates))
storage.mode(X_base) <- "double"

# --- Cell index by patient ---
cat("Indexing cells per patient...\n")
idx_by_patient <- split(
  which(full_p_id %in% eval_patients),
  full_p_id[full_p_id %in% eval_patients]
)

# Sample cells per patient for prediction features
set.seed(SEED + 1)
sample_idx <- unlist(lapply(eval_patients, function(pid) {
  ix <- idx_by_patient[[pid]]
  if (is.null(ix) || length(ix) == 0) return(integer())
  sample(ix, size = min(CELL_PER_PATIENT, length(ix)))
}), use.names = FALSE)
sample_pid <- full_p_id[sample_idx]
cat("  Sampled cells for features:", length(sample_idx), "\n")

# Sample cells for preservation evaluation
set.seed(SEED + 2)
idx_pres_all <- which(full_p_id %in% eval_patients)
idx_pres <- sample(idx_pres_all, size = min(CELL_SAMPLE_PRES, length(idx_pres_all)))
cat("  Sampled cells for preservation:", length(idx_pres), "\n")

# --- Iterate over folds ---
pres_rows <- list()
perf_rows <- list()
rank_list <- list()

for (f in 1:n_folds) {
  cat(sprintf("\n=== Fold %d ===\n", f))
  env <- new.env()
  load(file.path(PER_FOLD_ROOT, paste0("fold_", f), paste0("fold_", f, "_all.RData")), envir = env)

  centroids <- env$atlas_centroids   # 5000 x 30
  micro_to_meta <- env$micro_to_meta # length 5000
  scale_center <- env$atlas_scale_center
  scale_scale <- env$atlas_scale_scale
  train_ids <- env$train_p_ids
  test_ids <- env$test_p_ids
  n_meta <- length(unique(micro_to_meta))
  meta_levels <- sort(unique(micro_to_meta))
  cat("  k =", n_meta, "\n")

  # Fold indices in eval_patients
  tr_idx <- match(train_ids, eval_patients)
  te_idx <- match(test_ids, eval_patients)

  # Rank markers by between-meta-cluster variance of this fold's meta-centroids
  meta_centroids <- rowsum(centroids, group = micro_to_meta) /
    as.numeric(table(micro_to_meta))
  btw_var <- apply(meta_centroids, 2, var)
  fold_rank <- data.table(marker = marker_names, between_meta_var = btw_var)[order(-between_meta_var)]
  fold_rank[, rank := .I]
  rank_list[[f]] <- copy(fold_rank)[, fold := f]

  # Full-marker baseline for preservation
  X_pres_raw <- full_data[idx_pres, , drop = FALSE]
  X_pres_scaled <- scale(X_pres_raw, center = scale_center, scale = scale_scale)
  micro_pres_full <- RANN::nn2(data = centroids, query = X_pres_scaled, k = 1)$nn.idx[, 1]
  meta_pres_full <- micro_to_meta[micro_pres_full]

  for (b in BUDGETS) {
    cat(sprintf("  Budget %d...", b))
    markers_keep <- fold_rank[rank <= b, marker]
    col_idx <- match(markers_keep, marker_names)

    cent_sub <- centroids[, col_idx, drop = FALSE]
    center_sub <- scale_center[col_idx]
    scale_sub <- scale_scale[col_idx]

    # --- Cluster preservation ---
    X_pres_sub <- X_pres_raw[, col_idx, drop = FALSE]
    X_pres_scaled_sub <- scale(X_pres_sub, center = center_sub, scale = scale_sub)
    micro_pres_sub <- RANN::nn2(data = cent_sub, query = X_pres_scaled_sub, k = 1)$nn.idx[, 1]
    meta_pres_sub <- micro_to_meta[micro_pres_sub]
    agree <- mean(meta_pres_sub == meta_pres_full)
    pres_rows[[length(pres_rows) + 1]] <- data.table(fold = f, budget = b, meta_match = agree)

    # --- Build patient meta counts from cell sample ---
    X_raw_samp <- full_data[sample_idx, col_idx, drop = FALSE]
    X_scaled_samp <- scale(X_raw_samp, center = center_sub, scale = scale_sub)
    micro_samp <- RANN::nn2(data = cent_sub, query = X_scaled_samp, k = 1)$nn.idx[, 1]
    meta_samp <- micro_to_meta[micro_samp]

    meta_counts <- as.matrix(table(
      sample_pid,
      factor(meta_samp, levels = meta_levels)
    ))
    meta_counts <- meta_counts[eval_patients, , drop = FALSE]
    storage.mode(meta_counts) <- "integer"
    colnames(meta_counts) <- paste0("cluster_", meta_levels)

    saveRDS(meta_counts, file.path(OUT_DIR, "marker_reduction", sprintf("counts_fold%d_budget%d.rds", f, b)))

    # CLR transform
    clr_feat <- transform_counts(meta_counts, method = "clr")

    # Full feature matrix (clusters + all covariates = M4)
    X_full <- as.matrix(cbind(X_base, clr_feat))
    storage.mode(X_full) <- "double"

    # --- Fit M4 elastic net on this fold ---
    set.seed(SEED + b * 10 + f * 100)
    cvfit <- cv.glmnet(
      x = X_full[tr_idx, , drop = FALSE],
      y = Y_mat[tr_idx, , drop = FALSE],
      family = "mgaussian",
      alpha = ALPHA,
      standardize = TRUE,
      foldid = compute_grouped_foldid(train_ids, 5L, SEED + f * 100)
    )
    pred_raw <- predict(cvfit, newx = X_full[te_idx, , drop = FALSE], s = cvfit$lambda.min)
    pred_te <- do.call(cbind, lapply(1:length(VIRUSES), function(i) pred_raw[, i, 1]))
    colnames(pred_te) <- VIRUSES

    for (v in VIRUSES) {
      clust_mse <- mean((pred_te[, v] - Y_mat[te_idx, v])^2)
      perf_rows[[length(perf_rows) + 1]] <- data.table(
        fold = f, budget = b, virus = v, cluster_test_mse = clust_mse
      )
    }
    cat(" agree=", round(agree, 3), "\n")
  }
  rm(env); gc()
}

pres_dt <- rbindlist(pres_rows)
perf_dt <- rbindlist(perf_rows)
rank_dt <- rbindlist(rank_list)

timings <- list(total_s = as.numeric(difftime(Sys.time(), tic_total, units = "secs")))
cat(sprintf("\nTotal time: %.1f minutes\n", timings$total_s / 60))

# --- Save ---
save(pres_dt, perf_dt, rank_dt, BUDGETS, timings,
     file = file.path(OUT_DIR, "marker_reduction", "E2_marker_reduction_perfold.RData"))
fwrite(pres_dt, file.path(OUT_DIR, "marker_reduction", "E2_preservation.csv"))
fwrite(perf_dt, file.path(OUT_DIR, "marker_reduction", "E2_marker_perf.csv"))
fwrite(rank_dt, file.path(OUT_DIR, "marker_reduction", "E2_marker_ranking.csv"))

}  # end recompute-or-load guard

# --- Figure ---
cat("Generating figure...\n")

# Aggregate across folds
pres_agg <- pres_dt[, .(mean_match = mean(meta_match), sd_match = sd(meta_match)), by = budget]

# M1 baseline MSE (from aggregated results)
load(file.path(PER_FOLD_ROOT, "aggregated", "E2_aggregated.RData"))
m1_overall <- overall_mse["M1_en"]

perf_dt[, overall_mse := mean(cluster_test_mse), by = .(fold, budget)]
fold_n <- setNames(lengths(fold_assignments), seq_along(fold_assignments))
perf_dt[, n_test := as.numeric(fold_n[as.character(fold)])]
stopifnot(!anyNA(perf_dt$n_test), sum(fold_n) %in% c(343, 335, 262, 260))  # all episodes or first season
perf_agg <- perf_dt[virus == VIRUSES[1], .(
  mean_mse = weighted.mean(overall_mse, n_test),
  sd_mse = sd(overall_mse)
), by = budget]
perf_agg[, improvement := 100 * (1 - mean_mse / m1_overall)]
perf_agg[, rmse := sqrt(mean_mse)]
RMSE_LAB <- "Test RMSE (M4, overall)"
if (Sys.getenv("OUTCOME") == "composite") {
  # Composite-outcome refit from the saved counts (E17-03-marker-reduction-composite.R)
  cs <- fread(file.path(OUT_DIR, "marker_reduction", "E2_marker_summary_composite.csv"))
  perf_agg <- cs[, .(budget, mean_mse = mse, improvement = gain_vs_M1, rmse)]
  m1_overall <- cs$m1_mse[1]
  RMSE_LAB <- "Test RMSE (M4, composite)"
}
fwrite(perf_agg, file.path(OUT_DIR, "marker_reduction", "E2_marker_summary.csv"))
fwrite(pres_agg, file.path(OUT_DIR, "marker_reduction", "E2_preservation_summary.csv"))

# Panel A: Preservation
p_a <- ggplot(pres_agg, aes(x = budget, y = 100 * mean_match)) +
  geom_line(color = "#1B4F72", linewidth = 0.7) +
  geom_point(color = "#1B4F72", size = 2) +
  geom_errorbar(aes(ymin = 100 * (mean_match - sd_match),
                     ymax = 100 * (mean_match + sd_match)),
                width = 0.5, color = "#1B4F72") +
  geom_vline(xintercept = 15, linetype = "dashed", color = "#7F8C8D") +
  scale_x_reverse(breaks = BUDGETS) +
  labs(x = "Number of markers", y = "Cluster preservation (%)") +
  theme_minimal(base_size = 9)

# Panel B: RMSE vs budget
rmse_full <- sqrt(m1_overall)
p_b <- ggplot(perf_agg, aes(x = budget, y = rmse)) +
  geom_line(color = "#E67E22", linewidth = 0.7) +
  geom_point(color = "#E67E22", size = 2) +
  geom_hline(yintercept = rmse_full, linetype = "dotted", color = "gray50") +
  geom_vline(xintercept = 15, linetype = "dashed", color = "#7F8C8D") +
  scale_x_reverse(breaks = BUDGETS) +
  labs(x = "Number of markers", y = RMSE_LAB) +
  theme_minimal(base_size = 9)

# Panel C: consensus top-15 marker panel, colored by marker type
consensus <- rank_dt[, .(mean_var = mean(between_meta_var)), by = marker][order(-mean_var)]
consensus[, rank := .I]
top15 <- consensus[rank <= 15]
fwrite(consensus, file.path(OUT_DIR, "marker_reduction", "E2_marker_consensus.csv"))

marker_type <- c(
  CD19 = "Lineage", CD20 = "Lineage", CD3 = "Lineage", CD4 = "Lineage",
  CD8 = "Lineage", CD21 = "Lineage",
  CXCR5 = "Homing", CCR6 = "Homing", CCR2 = "Homing",
  `HLA-DR` = "Activation", CD86 = "Activation", CD25 = "Activation", CD38 = "Activation",
  CD27 = "Differentiation", CD57 = "Differentiation"
)
top15[, mtype := marker_type[as.character(marker)]]
top15[is.na(mtype), mtype := "Other"]
top15[, mtype := factor(mtype,
  levels = c("Lineage", "Homing", "Activation", "Differentiation", "Other"))]
type_cols <- c(Lineage = "#4393c3", Homing = "#5aae61",
               Activation = "#d6604d", Differentiation = "#f4a582", Other = "gray60")

p_c <- ggplot(top15, aes(x = reorder(marker, mean_var), y = mean_var, fill = mtype)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = type_cols, name = "Marker type", drop = TRUE) +
  labs(x = NULL, y = "Mean between-meta-cluster variance") +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom")

fig <- (p_a | p_b) / p_c +
  plot_layout(heights = c(1, 0.8)) +
  plot_annotation(tag_levels = "A", tag_suffix = ".")
fig <- fig & theme(plot.tag = element_text(face = "bold", size = 10),
                   panel.grid.minor = element_blank(),
                   axis.text = element_text(size = 7),
                   axis.title = element_text(size = 8))
ggsave(file.path(FIG_DIR, "E2_marker_reduction_perfold.pdf"), fig, width = 7.2, height = 5.7)
ggsave(file.path(FIG_DIR, "E2_marker_reduction_perfold.png"), fig, width = 7.2, height = 5.7, dpi = 600)
file.copy(file.path(FIG_DIR, "E2_marker_reduction_perfold.png"),
          "/home/amir/repo/flu-cellular-prediction/paper/images/E2_marker_reduction_perfold.png",
          overwrite = TRUE)

cat("\nFigure saved to:", file.path(FIG_DIR, "E2_marker_reduction_perfold.png"), "\n")
cat("Results saved to:", file.path(OUT_DIR, "marker_reduction"), "\n")
cat("Done.\n")
