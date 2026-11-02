#!/usr/bin/env Rscript
# E14-10: More supplementary figures for the first-season analysis, adapted from
# E8-review-figures.R and E9-refresh-supporting.R (plotting code unchanged):
#   E8_predictive_recovery.png  cell-label agreement on the shared atlas (E6-02)
#   E8_cluster_sizes.png        population abundance and selection (E6-04)
#   figS9_counterfactual_scatter.png  fitted responses under the two dose models
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(randomForest)})
src <- Sys.getenv("FLU_MEMO_DIR"); out <- "paper/images"; ROOT <- Sys.getenv("FLU_RESULTS_ROOT")
stopifnot(nzchar(src), nzchar(ROOT))
theme_pub <- theme_minimal(base_size = 9, base_family = "sans") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = "bottom", legend.title = element_blank(),
        axis.text = element_text(size = 7, color = "grey20"), axis.title = element_text(size = 8),
        strip.text = element_text(face = "bold", size = 9), plot.tag = element_text(face = "bold", size = 11))
save_plot <- function(plot, name, height) ggsave(file.path(out, name), plot, width = 7.2, height = height, dpi = 600, bg = "white")
recovery <- fread(file.path(src, "memo_tab4_recovery_curve.csv"))
recovery_plot <- ggplot(recovery, aes(budget, 100 * preservation, color = set, shape = set)) +
  geom_vline(xintercept = 15, linetype = "dashed", color = "grey55", linewidth = 0.35) +
  geom_line(linewidth = 0.7) + geom_point(size = 2) +
  scale_x_reverse(breaks = c(30, 25, 20, 15, 12, 10, 8, 6)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  scale_color_manual(values = c("All clusters" = "#7F8C8D", "Predictive clusters" = "#1B4F72"),
                     labels = c("All cells", "Cells in selected populations")) +
  scale_shape_manual(values = c("All clusters" = 16, "Predictive clusters" = 17),
                     labels = c("All cells", "Cells in selected populations")) +
  labs(x = "Number of retained markers", y = "Cells retaining their population label (%)") + theme_pub
save_plot(recovery_plot, "E8_predictive_recovery.png", 3.8)
sizes <- fread(file.path(src, "memo_tab6_cluster_sizes.csv"))
sizes[, cluster_f := factor(paste0("C", cluster_id), levels = paste0("C", cluster_id))]
sizes[, role := factor(ifelse(predictive, "Stably selected", "Not selected"), levels = c("Stably selected", "Not selected"))]
size_plot <- ggplot(sizes, aes(cluster_f, 100 * median_freq, color = role, shape = role)) +
  geom_hline(yintercept = c(0.1, 0.5), linetype = c("dotted", "dashed"), color = "grey55", linewidth = 0.35) +
  geom_linerange(aes(ymin = 100 * q25_freq, ymax = 100 * q75_freq), color = "grey65", linewidth = 0.5) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Stably selected" = "#1B4F72", "Not selected" = "#7F8C8D")) +
  scale_shape_manual(values = c("Stably selected" = 17, "Not selected" = 16)) +
  scale_y_continuous(trans = "log10", breaks = c(0.1, 0.5, 1, 5, 10)) +
  labs(x = "Population (ordered by median frequency)", y = "Cells per episode (%, median and IQR)") + theme_pub +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
save_plot(size_plot, "E8_cluster_sizes.png", 3.8)
# Fitted, partly in-sample predictions of the composite response from random
# forests fitted separately in each dose arm (E16-02, composite outcome).
tl <- fread(file.path(ROOT, "post_hoc_5000/composite_main/dose_models.csv"))
rho <- cor(tl$y_hat_hd, tl$y_hat_sd, method = "spearman"); cat(sprintf("Spearman rho between dose models: %.2f\n", rho))
lim <- range(c(tl$y_hat_hd, tl$y_hat_sd))
p <- ggplot(tl, aes(y_hat_sd, y_hat_hd, color = subgroup)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(alpha = .7, size = 1.3) + coord_equal(xlim = lim, ylim = lim) +
  scale_color_manual(values = c("Higher predicted response" = "#1B4F72", "Intermediate" = "grey65", "Lower predicted response" = "#E67E22")) +
  labs(x = "Fitted composite log fold-rise, standard-dose model", y = "Fitted composite log fold-rise, high-dose model", color = NULL) +
  theme_minimal(base_size = 8, base_family = "sans") + theme(panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave(file.path(out, "figS9_counterfactual_scatter.png"), p, width = 4.6, height = 4.6, dpi = 600, bg = "white")
cat("Done.\n")
