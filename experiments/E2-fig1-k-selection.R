#!/usr/bin/env Rscript
# E2-fig1-k-selection.R
# Regenerate Figure 1 (k-selection performance) from 5000-micro per-fold data.
# 3-panel figure: (A) Inner CV heatmap, (B) k per fold, (C) Per-fold MSE improvement

library(ggplot2)
library(patchwork)

cat("=== Generating Figure 1 from per-fold atlas (5000 micro) ===\n")

PER_FOLD_ROOT <- file.path(Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06"), "per_fold_atlas_5000")
OUT_DIR <- file.path(Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06"), "post_hoc_5000/figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Load all fold data ---
n_folds <- 5
fold_data <- list()
for (f in 1:n_folds) {
  rdata_file <- file.path(PER_FOLD_ROOT, paste0("fold_", f), paste0("fold_", f, "_all.RData"))
  cat("Loading", rdata_file, "\n")
  env <- new.env()
  load(rdata_file, envir = env)
  fold_data[[f]] <- env
}

# --- Panel A: Inner CV heatmap (resolution_scores across folds) ---
# Extract resolution_scores from each fold
# resolution_scores is a named numeric vector: names = resolutions, values = inner CV scores
heatmap_df <- do.call(rbind, lapply(1:n_folds, function(f) {
  scores <- fold_data[[f]]$resolution_scores
  k_by_res <- fold_data[[f]]$k_found_by_res
  sel_res <- fold_data[[f]]$selected_resolution
  data.frame(
    fold = paste0("Fold ", f),
    resolution = as.numeric(names(scores)),
    score = as.numeric(scores),
    k = as.integer(k_by_res),
    is_selected = abs(as.numeric(names(scores)) - sel_res) < 0.01,
    stringsAsFactors = FALSE
  )
}))

# The score is the inner CV median MSE ratio (lower = better, < 1 means clusters beat baseline)
p_a <- ggplot(heatmap_df, aes(x = factor(resolution), y = fold, fill = score)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_point(data = subset(heatmap_df, is_selected),
             aes(x = factor(resolution), y = fold),
             shape = 8, size = 3, color = "red", inherit.aes = FALSE) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 1.0, name = "Median\nMSE ratio",
                       limits = c(min(heatmap_df$score, 0.85), max(heatmap_df$score, 1.15))) +
  labs(x = "Leiden resolution", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "right",
        panel.grid = element_blank())

# --- Panel B: k selected per fold ---
k_df <- data.frame(
  fold = paste0("Fold ", 1:n_folds),
  k = sapply(1:n_folds, function(f) fold_data[[f]]$k_found),
  resolution = sapply(1:n_folds, function(f) fold_data[[f]]$selected_resolution)
)

p_b <- ggplot(k_df, aes(x = fold, y = k)) +
  geom_col(fill = "#4393c3", width = 0.6) +
  geom_text(aes(label = k), vjust = -0.3, size = 4, fontface = "bold") +
  geom_hline(yintercept = median(k_df$k), linetype = "dashed", color = "red", linewidth = 0.5) +
  annotate("text", x = 5.4, y = median(k_df$k) + 1.5,
           label = paste0("median = ", median(k_df$k)), color = "red", hjust = 1, size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = expression("Meta-clusters (" * italic(k) * ")")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

# --- Panel C: Per-fold per-virus MSE improvement (M3 vs M1) ---
# Compute MSE per fold per virus
perf_list <- list()
for (f in 1:n_folds) {
  env <- fold_data[[f]]
  Y_test <- env$Y_test
  pred_M1 <- env$pred_M1_en_test
  pred_M3 <- env$pred_M3_en_test

  viruses <- colnames(Y_test)
  for (v in viruses) {
    mse_m1 <- mean((Y_test[, v] - pred_M1[, v])^2)
    mse_m3 <- mean((Y_test[, v] - pred_M3[, v])^2)
    improvement <- 100 * (1 - mse_m3 / mse_m1)
    perf_list[[length(perf_list) + 1]] <- data.frame(
      fold = paste0("Fold ", f),
      virus = toupper(v),
      mse_m1 = mse_m1,
      mse_m3 = mse_m3,
      improvement = improvement,
      stringsAsFactors = FALSE
    )
  }
}
perf_df <- do.call(rbind, perf_list)
perf_df$virus <- gsub("BVIC", "B/Vic", perf_df$virus)

# Compute overall improvement per fold
overall_perf <- do.call(rbind, lapply(1:n_folds, function(f) {
  env <- fold_data[[f]]
  Y_test <- env$Y_test
  n_test <- nrow(Y_test)
  pred_M1 <- env$pred_M1_en_test
  pred_M3 <- env$pred_M3_en_test
  mse_m1 <- mean((Y_test - pred_M1)^2)
  mse_m3 <- mean((Y_test - pred_M3)^2)
  data.frame(fold = paste0("Fold ", f), virus = "Overall",
             mse_m1 = mse_m1, mse_m3 = mse_m3,
             improvement = 100 * (1 - mse_m3 / mse_m1))
}))
perf_df <- rbind(perf_df, overall_perf)
perf_df$virus <- factor(perf_df$virus, levels = c("H1N1", "H3N2", "B/Vic", "Overall"))

p_c <- ggplot(perf_df, aes(x = fold, y = improvement, fill = virus)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c("H1N1" = "#66c2a5", "H3N2" = "#fc8d62",
                                "B/Vic" = "#8da0cb", "Overall" = "#e78ac3"),
                    labels = c("H1N1" = "H1N1", "H3N2" = "H3N2",
                               "B/Vic" = "B/Victoria", "Overall" = "Overall"),
                    name = "Antigen") +
  labs(x = NULL, y = "MSE improvement (%), M3 vs M1") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")

# --- Combine panels ---
fig1 <- p_a / (p_b | p_c) +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(tag_levels = "A", tag_suffix = ".")
fig1 <- fig1 & theme(plot.tag = element_text(face = "bold", size = 13))

ggsave(file.path(OUT_DIR, "E2_fig1_k_selection.pdf"), fig1, width = 12, height = 8)
ggsave(file.path(OUT_DIR, "E2_fig1_k_selection.png"), fig1, width = 12, height = 8, dpi = 300)
# Copy into the paper image directory
file.copy(file.path(OUT_DIR, "E2_fig1_k_selection.png"),
          "/home/amir/repo/flu-cellular-prediction/paper/images/E2_fig1_k_selection.png",
          overwrite = TRUE)

cat("\nFigure 1 saved to:\n")
cat("  PDF:", file.path(OUT_DIR, "E2_fig1_k_selection.pdf"), "\n")
cat("  PNG:", file.path(OUT_DIR, "E2_fig1_k_selection.png"), "\n")

# Also save as summary table
write.csv(perf_df, file.path(OUT_DIR, "E2_perfold_improvement.csv"), row.names = FALSE)
cat("Per-fold improvements saved to:", file.path(OUT_DIR, "E2_perfold_improvement.csv"), "\n")
