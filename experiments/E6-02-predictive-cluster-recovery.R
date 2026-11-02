#!/usr/bin/env Rscript
# E6-02-predictive-cluster-recovery.R
#
# The published marker-reduction curve asks how well a reduced antibody panel
# reproduces the FULL-panel cluster labels, averaged over every cell. That is a
# stricter question than the one a reader actually cares about: does the reduced
# panel still recover the populations the model relies on?
#
# This script recomputes cluster preservation on the production atlas (k = 33)
# and splits the curve two ways:
#   - all clusters      : every cell, i.e. the published quantity
#   - predictive only   : cells whose full-panel label is one of the clusters
#                         selected in at least 50% of bootstrap refits
#
# Cells come from sample_dt inside the production atlas object (200,000 cells
# with their stored micro and meta assignments). Reassigning those cells with
# the full 30-marker panel reproduces the stored labels exactly, which is what
# licenses using them here instead of re-reading the 5.9 GB raw cell matrix.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(RANN)
})

cat("=== E6-02: marker-budget recovery of the predictive clusters ===\n")

ATLAS <- Sys.getenv("FLU_PROD_ATLAS", "/home/amir/repo/flu-cellular-prediction/results/final_atlas_production/final_atlas_e862228d.RData")
BOOT  <- Sys.getenv("FLU_BOOT_FILE", "/home/amir/repo/flu-cellular-prediction/results/timing_corrected_2026-09-06/bootstrap_importance_3f0f3ae8.RData")
CACHED_PRES <- "/home/amir/repo/flu-cellular-prediction/results/marker_reduction_production/tables/marker_preservation_b4e8f6ed.csv"
CACHED_RANK <- "/home/amir/repo/flu-cellular-prediction/results/marker_reduction_production/tables/marker_ranking_b4e8f6ed.csv"
OUT_DIR <- Sys.getenv("FLU_MEMO_DIR", "/home/amir/repo/flu-cellular-prediction/paper/memo_figures")
# The cached production checks apply only to the all-episode shared atlas.
CHECK_CACHED <- !nzchar(Sys.getenv("FLU_PROD_ATLAS"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

BUDGETS <- c(30, 25, 20, 15, 12, 10, 8, 6)
SELECTION_THRESHOLD <- 0.5

# ---------------------------------------------------------------------------
# 1. Predictive cluster set
# ---------------------------------------------------------------------------
cat("\n[1] Identifying the stably selected (predictive) clusters...\n")
be <- new.env(); load(BOOT, envir = be)
stab <- as.data.table(be$stability)
cluster_stab <- stab[grepl("^cluster_", feature)]
max_freq <- cluster_stab[, .(max_selection_freq = max(selection_freq)), by = feature]
predictive_features <- max_freq[max_selection_freq >= SELECTION_THRESHOLD, feature]
predictive_ids <- sort(as.integer(sub("cluster_", "", predictive_features)))
cat(sprintf("  clusters with selection frequency >= %.0f%% in at least one antigen: %d\n",
            100 * SELECTION_THRESHOLD, length(predictive_ids)))
cat("  ", paste(predictive_ids, collapse = ", "), "\n")

# ---------------------------------------------------------------------------
# 2. Production atlas and marker ranking
# ---------------------------------------------------------------------------
cat("\n[2] Loading the production atlas...\n")
ae <- new.env(); load(ATLAS, envir = ae)
centroids <- ae$centroids                # 5000 micro-centroids in scaled space
micro_to_meta <- ae$micro_to_meta
marker_names <- ae$marker_names
scale_center <- ae$scale_center
scale_scale <- ae$scale_scale
sample_dt <- ae$sample_dt
n_meta <- length(unique(micro_to_meta))
cat(sprintf("  k = %d meta-clusters, %d micro-clusters, %d sampled cells\n",
            n_meta, nrow(centroids), nrow(sample_dt)))
stopifnot(!CHECK_CACHED || n_meta == 33, all(predictive_ids <= n_meta))

# Marker ranking: between-meta-cluster variance of the meta-centroid means,
# exactly as in the production marker-reduction run.
meta_centroids <- rowsum(centroids, group = micro_to_meta) / as.numeric(table(micro_to_meta))
btw_var <- apply(meta_centroids, 2, var)
rank_dt <- data.table(marker = marker_names, between_meta_var = btw_var)[order(-between_meta_var)]
rank_dt[, rank := .I]

if (CHECK_CACHED) {
cached_rank <- fread(CACHED_RANK)
rank_agree <- identical(rank_dt$marker, cached_rank$marker)
cat("  marker ranking reproduces the cached production ranking:", rank_agree, "\n")
if (!rank_agree) {
  cat("  cached top 15:", paste(cached_rank$marker[1:15], collapse = ", "), "\n")
  cat("  recomputed   :", paste(rank_dt$marker[1:15], collapse = ", "), "\n")
}
}

# ---------------------------------------------------------------------------
# 3. Full-panel reference labels for the sampled cells
# ---------------------------------------------------------------------------
cat("\n[3] Verifying full-panel reassignment reproduces the stored labels...\n")
X_raw <- as.matrix(sample_dt[, marker_names, with = FALSE])
X_scaled <- scale(X_raw, center = scale_center, scale = scale_scale)
micro_full <- RANN::nn2(data = centroids, query = X_scaled, k = 1)$nn.idx[, 1]
meta_full <- micro_to_meta[micro_full]
cat(sprintf("  micro agreement %.4f, meta agreement %.4f\n",
            mean(micro_full == sample_dt$micro), mean(meta_full == sample_dt$meta)))
stopifnot(all(meta_full == sample_dt$meta))

is_predictive_cell <- meta_full %in% predictive_ids
cat(sprintf("  cells in predictive clusters: %d of %d (%.1f%%)\n",
            sum(is_predictive_cell), length(meta_full), 100 * mean(is_predictive_cell)))

# ---------------------------------------------------------------------------
# 4. Preservation vs marker budget
# ---------------------------------------------------------------------------
cat("\n[4] Sweeping marker budgets...\n")
curve_rows <- list(); percluster_rows <- list()

for (b in BUDGETS) {
  keep <- rank_dt[rank <= b, marker]
  col_idx <- match(keep, marker_names)

  X_sub <- scale(X_raw[, col_idx, drop = FALSE],
                 center = scale_center[col_idx], scale = scale_scale[col_idx])
  micro_sub <- RANN::nn2(data = centroids[, col_idx, drop = FALSE],
                         query = X_sub, k = 1)$nn.idx[, 1]
  meta_sub <- micro_to_meta[micro_sub]

  match_vec <- meta_sub == meta_full
  curve_rows[[length(curve_rows) + 1L]] <- data.table(
    budget = b,
    set = c("All clusters", "Predictive clusters"),
    preservation = c(mean(match_vec), mean(match_vec[is_predictive_cell])),
    n_cells = c(length(match_vec), sum(is_predictive_cell))
  )

  pc <- data.table(cluster = meta_full, matched = match_vec)[
    , .(n_cells = .N, preservation = mean(matched)), by = cluster][order(cluster)]
  pc[, budget := b]
  pc[, predictive := cluster %in% predictive_ids]
  percluster_rows[[length(percluster_rows) + 1L]] <- pc

  cat(sprintf("  budget %2d: all %.4f | predictive %.4f\n",
              b, mean(match_vec), mean(match_vec[is_predictive_cell])))
}

curve_dt <- rbindlist(curve_rows)
percluster_dt <- rbindlist(percluster_rows)

# Cross-check the "all clusters" curve against the cached production run, which
# used a different random 100,000-cell sample drawn from the full cell matrix.
chk <- data.table()
if (CHECK_CACHED) {
cached <- fread(CACHED_PRES)
chk <- merge(curve_dt[set == "All clusters", .(budget, recomputed = preservation)],
             cached[, .(budget, cached = meta_match)], by = "budget")
chk[, abs_diff := abs(recomputed - cached)]
cat("\n  Cross-check against the cached production preservation curve:\n")
print(chk)
}

fwrite(curve_dt, file.path(OUT_DIR, "memo_tab4_recovery_curve.csv"))
fwrite(percluster_dt, file.path(OUT_DIR, "memo_tab4_recovery_per_cluster.csv"))
fwrite(chk, file.path(OUT_DIR, "memo_tab4_recovery_crosscheck.csv"))

# ---------------------------------------------------------------------------
# 5. Figure
# ---------------------------------------------------------------------------
cat("\n[5] Drawing the figure...\n")

SET_COLORS <- c("All clusters" = "#7f7f7f", "Predictive clusters" = "#2166ac")

fig4 <- ggplot(curve_dt, aes(x = budget, y = 100 * preservation,
                             color = set, shape = set)) +
  geom_vline(xintercept = 15, linetype = "dashed", color = "#b2182b", linewidth = 0.45) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.8) +
  annotate("text", x = 15.4, y = 30, label = "15-marker panel", hjust = 0,
           size = 3.2, color = "#b2182b") +
  scale_x_reverse(breaks = BUDGETS) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  scale_color_manual(values = SET_COLORS, name = NULL) +
  scale_shape_manual(values = c("All clusters" = 16, "Predictive clusters" = 17), name = NULL) +
  labs(x = "Number of markers in the panel",
       y = "Cluster label preservation (%)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

ggsave(file.path(OUT_DIR, "memo_fig4_predictive_recovery.png"), fig4,
       width = 7.6, height = 5.2, dpi = 300)
cat("  wrote memo_fig4_predictive_recovery.png\n")

cat("\nPer-cluster preservation at the 15-marker budget (predictive clusters):\n")
print(percluster_dt[budget == 15 & predictive == TRUE][order(preservation)])

cat("\nDone. Outputs in", OUT_DIR, "\n")
