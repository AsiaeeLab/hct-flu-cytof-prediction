#!/usr/bin/env Rscript
# Publication views of the September 6 timing-corrected analyses.
# Plotting only: E2/E6/E9 create the fitted models and numerical summaries.
# Prior memo figures are preserved in paper/archive/2026-09-06_pre-timing-correction/.
FLU_ROOT <- Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06")
FLU_ANN_FILE <- Sys.getenv("FLU_ANN_FILE", "/home/amir/repo/flu-cellular-prediction/results/reports/cluster_annotations_corrected.csv")
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
src <- "paper/memo_figures"
out <- "paper/images"
groups <- fread(file.path(src, "memo_tab1_score_groups.csv"))
episodes <- as.data.table(readRDS(file.path(src, "memo_episode_scores.rds")))
stopifnot(nrow(episodes) > 0)  # 343 episodes (all) or 262 (first season)
qlevels <- c("Q1_low", "Q2", "Q3", "Q4_high")
pooled <- "Q1-Q3 pooled (within-quartile split)"
palette <- c(High = "#1B4F72", Low = "#E67E22")
qpalette <- c(Q1_low = "#1B4F72", Q2 = "#85C1E9", Q3 = "#E67E22", Q4_high = "#C0392B")
theme_pub <- theme_minimal(base_size = 9, base_family = "sans") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = "bottom", legend.title = element_blank(),
        axis.text = element_text(size = 7, color = "grey20"),
        axis.title = element_text(size = 8),
        strip.text = element_text(face = "bold", size = 9),
        plot.tag = element_text(face = "bold", size = 11))
dodge <- position_dodge(width = 0.5)
save_plot <- function(plot, name, height) {
  ggsave(file.path(out, name), plot, width = 7.2, height = height,
         dpi = 600, bg = "white")
}

# Main: raw clinical outcomes, with the pooled lower-three-quartile summary.
d <- groups[virus == "h3n2" & stratum %in% c(qlevels, pooled)]
d[, stratum := factor(stratum, levels = c(qlevels, pooled))]
d[, group := factor(group, levels = c("High", "Low"))]
qlabels <- c("Q1\nLowest", "Q2\n", "Q3\n", "Q4\nHighest", "Q1-Q3\nPooled")
tests_path <- file.path(FLU_ROOT, "score_null_bootstrap/score_tests.csv")
tests <- if (file.exists(tests_path)) fread(tests_path) else NULL
if (!is.null(tests)) stopifnot(nrow(tests) == 30L, all(tests$n_null == 999L))
p_labels <- function(virus_name, endpoint_name, scopes, label) {
  if (is.null(tests)) return(label)
  tt <- tests[virus == virus_name & endpoint == endpoint_name]
  pp <- tt$p_holm[match(scopes, tt$scope)]
  stopifnot(!anyNA(pp))
  paste0(label, "\n", "pH=", sub("^0", "", sprintf("%.3f", pp)))
}
base <- ggplot(d, aes(stratum, color = group, shape = group, group = group)) +
  scale_color_manual(values = palette, labels = c("Higher cellular score", "Lower cellular score")) +
  scale_shape_manual(values = c(16, 17), labels = c("Higher cellular score", "Lower cellular score")) +
  scale_x_discrete(labels = qlabels) + theme_pub +
  labs(x = "Pre-vaccination titer quartile")
p1 <- base +
  geom_hline(yintercept = 4, linetype = "dashed", color = "grey55", linewidth = 0.35) +
  geom_linerange(aes(ymin = gmfr_lo, ymax = gmfr_hi), position = dodge, linewidth = 0.6) +
  geom_point(aes(y = gmfr), position = dodge, size = 2.3) +
  scale_y_continuous(trans = "log2", breaks = c(1, 2, 4, 8, 16, 32, 64, 128)) +
  labs(y = "Geometric mean H3N2 fold-rise")
if (!is.null(tests)) p1 <- p1 + scale_x_discrete(labels = p_labels(
  "h3n2", "log_fold_ratio", c(paste0("Q", 1:4), "Q1-Q3"), qlabels))
p2 <- base +
  geom_linerange(aes(ymin = 100 * seroconv_lo, ymax = 100 * seroconv_hi),
                 position = dodge, linewidth = 0.6) +
  geom_point(aes(y = 100 * seroconv), position = dodge, size = 2.3) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(y = "Fourfold H3N2 response (%)")
if (!is.null(tests)) p2 <- p2 + scale_x_discrete(labels = p_labels(
  "h3n2", "response_difference", c(paste0("Q", 1:4), "Q1-Q3"), qlabels))
main <- (p1 | p2) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") & theme(legend.position = "bottom")
save_plot(main, "E8_cellular_score_clinical.png", 4.1)

# Supplement: continuous score, showing the substantial individual overlap.
# The descriptive regression lines have no row-independent confidence band.
episodes[, quartile := factor(quartile, levels = qlevels)]
continuous <- ggplot(episodes, aes(score_h3n2, exp(y_h3n2), color = quartile)) +
  geom_hline(yintercept = 4, linetype = "dashed", color = "grey55", linewidth = 0.35) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.6) +
  scale_y_continuous(trans = "log2", breaks = c(0.125, 0.5, 2, 8, 32, 128, 512, 2048),
                     labels = c("0.125", "0.5", "2", "8", "32", "128", "512", "2048")) +
  scale_color_manual(values = qpalette, labels = c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)")) +
  labs(x = "Cellular response score (out-of-fold)", y = "Observed H3N2 fold-rise") +
  theme_pub
save_plot(continuous, "E8_score_continuous.png", 4.2)

# Supplement: each antigen uses its own cellular score and median split.
all_d <- groups[stratum %in% qlevels]
all_d[, stratum := factor(stratum, levels = qlevels)]
all_d[, antigen := factor(virus, levels = c("h1n1", "h3n2", "bvic"),
                          labels = c("H1N1", "H3N2", "B/Victoria"))]
all_d[, group := factor(group, levels = c("High", "Low"))]
all_d[, x_label := p_labels(virus[1], "log_fold_ratio", paste0("Q", 1:4), paste0("Q", 1:4))[match(stratum,qlevels)], by = virus]
all_d[, x_key := factor(paste(virus,stratum,sep="_"),
  levels=unlist(lapply(c("h1n1","h3n2","bvic"),function(v) paste(v,qlevels,sep="_"))))]
label_map <- setNames(all_d$x_label,as.character(all_d$x_key))
all_antigens <- ggplot(all_d, aes(x_key, gmfr, color = group, shape = group)) +
  geom_hline(yintercept = 4, linetype = "dashed", color = "grey55", linewidth = 0.35) +
  geom_linerange(aes(ymin = gmfr_lo, ymax = gmfr_hi), position = dodge, linewidth = 0.5) +
  geom_point(position = dodge, size = 2) +
  facet_wrap(~antigen, nrow = 1, scales = "free_x") +
  scale_x_discrete(labels = function(x) unname(label_map[x])) +
  scale_y_continuous(trans = "log2", breaks = c(1, 2, 4, 8, 16, 32, 64, 128)) +
  scale_color_manual(values = palette, labels = c("Higher cellular score", "Lower cellular score")) +
  scale_shape_manual(values = c(16, 17), labels = c("Higher cellular score", "Lower cellular score")) +
  labs(x = "Pre-vaccination titer quartile", y = "Geometric mean fold-rise") + theme_pub
save_plot(all_antigens, "E8_score_all_antigens.png", 3.3)
combined <- (continuous / all_antigens) + plot_layout(heights = c(1.2, 1)) +
  plot_annotation(tag_levels = "A")
save_plot(combined, "E8_score_supporting.png", 6.7)

# Supplement: label agreement is cell-weighted, not a count of populations.
recovery <- fread(file.path(src, "memo_tab4_recovery_curve.csv"))
recovery_plot <- ggplot(recovery, aes(budget, 100 * preservation, color = set, shape = set)) +
  geom_vline(xintercept = 15, linetype = "dashed", color = "grey55", linewidth = 0.35) +
  geom_line(linewidth = 0.7) + geom_point(size = 2) +
  scale_x_reverse(breaks = c(30, 25, 20, 15, 12, 10, 8, 6)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  scale_color_manual(values = c("All clusters" = "#7F8C8D", "Predictive clusters" = "#1B4F72"),
                     labels = c("All cells", "Cells in selected clusters")) +
  scale_shape_manual(values = c("All clusters" = 16, "Predictive clusters" = 17),
                     labels = c("All cells", "Cells in selected clusters")) +
  labs(x = "Number of retained markers", y = "Cells retaining their cluster label (%)") + theme_pub
save_plot(recovery_plot, "E8_predictive_recovery.png", 3.8)

# Supplement: size is descriptive QC, not evidence that a cluster is noise-free.
sizes <- fread(file.path(src, "memo_tab6_cluster_sizes.csv"))
sizes[, cluster_f := factor(paste0("C", cluster_id), levels = paste0("C", cluster_id))]
sizes[, role := factor(ifelse(predictive, "Stably selected", "Not selected"),
                       levels = c("Stably selected", "Not selected"))]
size_plot <- ggplot(sizes, aes(cluster_f, 100 * median_freq, color = role, shape = role)) +
  geom_hline(yintercept = c(0.1, 0.5), linetype = c("dotted", "dashed"),
             color = "grey55", linewidth = 0.35) +
  geom_linerange(aes(ymin = 100 * q25_freq, ymax = 100 * q75_freq), color = "grey65", linewidth = 0.5) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Stably selected" = "#1B4F72", "Not selected" = "#7F8C8D")) +
  scale_shape_manual(values = c("Stably selected" = 17, "Not selected" = 16)) +
  scale_y_continuous(trans = "log10", breaks = c(0.1, 0.5, 1, 5, 10)) +
  labs(x = "Meta-cluster (ordered by median frequency)",
       y = "Cells per episode (%, median and IQR)") + theme_pub +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
save_plot(size_plot, "E8_cluster_sizes.png", 3.8)

# Refresh the longitudinal figure's annotations from its saved summaries.
# Numerical results and original dimensions are retained.
repeater_env <- new.env(parent = emptyenv())
load(file.path(FLU_ROOT, "post_hoc_5000/repeater_analysis/E4_repeater_workspace.RData"),
     envir = repeater_env)
delta_summary <- copy(as.data.table(repeater_env$delta_summary))
ann <- fread(FLU_ANN_FILE)
delta_summary[, label := sprintf("C%d: %s", cluster_id,
                                ann$annotation[match(cluster_id, ann$cluster_id)])]
top_movers <- delta_summary[!is.na(trans_h3n2),
  .(rng = max(abs(mean_delta))), by = cluster_id][rng > 0.15, cluster_id]
delta_top <- delta_summary[cluster_id %in% top_movers & !is.na(trans_h3n2)]
trans_colors <- c("NR -> R" = "#2ca02c", "R -> R" = "#1f77b4",
                  "NR -> NR" = "#d6604d", "R -> NR" = "#762a83")
delta_plot <- ggplot(delta_top, aes(mean_delta, reorder(label, mean_delta), color = trans_h3n2)) +
  geom_point(size = 2.6, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(xmin = mean_delta - se_delta, xmax = mean_delta + se_delta),
                orientation = "y", position = position_dodge(width = 0.5), width = 0.25) +
  geom_vline(xintercept = 0, color = "gray40", linetype = "dashed") +
  scale_color_manual(values = trans_colors, name = "H3N2 transition",
                     labels = function(x) gsub(" -> ", " → ", x)) +
  labs(x = "Mean CLR change (year 2 − year 1)", y = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 7))
ggsave(file.path(out, "E4_repeater_cluster_delta_h3n2.png"), delta_plot,
       width = 9, height = 8, dpi = 300)
