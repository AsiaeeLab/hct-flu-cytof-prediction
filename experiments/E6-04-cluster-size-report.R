#!/usr/bin/env Rscript
# E6-04-cluster-size-report.R
#
# Diagnostic only: how large is each production-atlas meta-cluster in a typical
# recipient, and which clusters are small enough that their CLR feature is
# mostly measurement noise?
#
# Frequencies are the per-episode fraction of that recipient's profiled cells
# assigned to the cluster (meta_props_all, 343 episodes x 33 clusters). Nothing
# is dropped and no model is refit here; the output is a watch list.

FLU_BIO_FILE <- Sys.getenv("FLU_BIO_FILE", "/home/amir/repo/flu-cellular-prediction/results/final_atlas_production/biology_production_90203bf3.RData")
FLU_BOOT_FILE <- Sys.getenv("FLU_BOOT_FILE", "/home/amir/repo/flu-cellular-prediction/results/timing_corrected_2026-09-06/bootstrap_importance_3f0f3ae8.RData")
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

cat("=== E6-04: cluster size / noise report ===\n")

BIO <- FLU_BIO_FILE
BOOT <- FLU_BOOT_FILE
OUT_DIR <- Sys.getenv("FLU_MEMO_DIR", "/home/amir/repo/flu-cellular-prediction/paper/memo_figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

WATCH_THRESHOLD <- 0.005   # median frequency below 0.5% = watch list
NOISE_THRESHOLD <- 0.001   # median frequency below 0.1% = strong candidate

be <- new.env(); load(BIO, envir = be)
props <- be$meta_props_all
counts <- be$meta_counts_all
cat(sprintf("Episodes: %d, clusters: %d\n", nrow(props), ncol(props)))
stopifnot(nrow(props) %in% c(343, 335, 262, 260))  # all episodes or first season
cat(sprintf("Row sums of the proportion matrix: min %.4f, max %.4f\n",
            min(rowSums(props)), max(rowSums(props))))
cells_per_episode <- rowSums(counts)
cat(sprintf("Cells profiled per episode: median %.0f, range %.0f-%.0f\n",
            median(cells_per_episode), min(cells_per_episode), max(cells_per_episode)))
cat(sprintf("Episodes with < 10,000 cells: %d; < 5,000 cells: %d; < 1,000 cells: %d\n",
            sum(cells_per_episode < 10000), sum(cells_per_episode < 5000),
            sum(cells_per_episode < 1000)))
depth_dt <- data.table(episode_id = rownames(counts), n_cells = cells_per_episode)[
  order(n_cells)]
fwrite(depth_dt, file.path(OUT_DIR, "memo_tab6_cells_per_episode.csv"))

se <- new.env(); load(BOOT, envir = se)
stab <- as.data.table(se$stability)
cl_stab <- stab[grepl("^cluster_", feature),
                .(max_selection_freq = max(selection_freq)), by = feature]
predictive <- cl_stab[max_selection_freq >= 0.5, feature]

cluster_names <- colnames(props)
if (is.null(cluster_names)) cluster_names <- paste0("cluster_", seq_len(ncol(props)))

size_dt <- data.table(
  cluster = cluster_names,
  cluster_id = seq_len(ncol(props)),
  mean_freq = colMeans(props),
  median_freq = apply(props, 2, median),
  q25_freq = apply(props, 2, quantile, 0.25),
  q75_freq = apply(props, 2, quantile, 0.75),
  min_freq = apply(props, 2, min),
  max_freq = apply(props, 2, max),
  pct_episodes_zero = 100 * colMeans(counts == 0),
  median_cells_per_episode = apply(counts, 2, median)
)
size_dt[, predictive := cluster %in% predictive]
size_dt[, max_selection_freq := cl_stab$max_selection_freq[match(cluster, cl_stab$feature)]]
size_dt[, flag := fifelse(median_freq < NOISE_THRESHOLD, "Noise candidate (<0.1%)",
                  fifelse(median_freq < WATCH_THRESHOLD, "Watch list (<0.5%)", "Adequate"))]
setorder(size_dt, median_freq)

cat("\nPer-cluster frequency, sorted ascending (percent of a recipient's cells):\n")
print(size_dt[, .(cluster, mean_pct = round(100 * mean_freq, 3),
                  median_pct = round(100 * median_freq, 3),
                  iqr_pct = paste0(round(100 * q25_freq, 3), "-", round(100 * q75_freq, 3)),
                  median_cells = median_cells_per_episode,
                  pct_zero = round(pct_episodes_zero, 1),
                  predictive, flag)], nrows = 40)

cat("\nFlag counts:\n"); print(size_dt[, .N, by = flag])
cat("\nFlagged clusters that are also in the predictive set:\n")
print(size_dt[flag != "Adequate" & predictive == TRUE,
              .(cluster, median_pct = round(100 * median_freq, 3),
                max_selection_freq, flag)])

fwrite(size_dt, file.path(OUT_DIR, "memo_tab6_cluster_sizes.csv"))

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
plot_dt <- copy(size_dt)
plot_dt[, cluster_f := factor(cluster, levels = cluster)]
plot_dt[, role := fifelse(predictive, "Stably selected", "Not selected")]

fig6 <- ggplot(plot_dt, aes(x = cluster_f)) +
  geom_hline(yintercept = 100 * WATCH_THRESHOLD, linetype = "dashed",
             color = "#b2182b", linewidth = 0.45) +
  geom_hline(yintercept = 100 * NOISE_THRESHOLD, linetype = "dotted",
             color = "#b2182b", linewidth = 0.45) +
  geom_linerange(aes(ymin = 100 * q25_freq, ymax = 100 * q75_freq),
                 color = "grey60", linewidth = 0.7) +
  geom_point(aes(y = 100 * median_freq, color = role, shape = role), size = 2.4) +
  annotate("text", x = ncol(props) - 0.2, y = 100 * WATCH_THRESHOLD * 1.25,
           label = "0.5%", hjust = 1, size = 3, color = "#b2182b") +
  annotate("text", x = ncol(props) - 0.2, y = 100 * NOISE_THRESHOLD * 1.25,
           label = "0.1%", hjust = 1, size = 3, color = "#b2182b") +
  scale_y_log10(breaks = c(0.01, 0.1, 0.5, 1, 5, 10, 25),
                labels = c("0.01", "0.1", "0.5", "1", "5", "10", "25")) +
  scale_color_manual(values = c("Stably selected" = "#2166ac",
                                "Not selected" = "grey45"), name = NULL) +
  scale_shape_manual(values = c("Stably selected" = 17, "Not selected" = 16),
                     name = NULL) +
  labs(x = "Meta-cluster (sorted by median frequency)",
       y = "Frequency within a recipient (%), median and IQR") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

ggsave(file.path(OUT_DIR, "memo_fig6_cluster_sizes.png"), fig6,
       width = 10.5, height = 5.4, dpi = 300)
cat("\nWrote memo_fig6_cluster_sizes.png\n")
cat("Done. Outputs in", OUT_DIR, "\n")
