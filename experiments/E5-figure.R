# Per-cohort comparison figure: side-by-side percent improvement (M1->M3
# and M2->M4) with recipient-clustered bootstrap CIs, and side-by-side
# selection frequency for the top shared features.

FLU_ROOT <- Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06")
FLU_ANN_FILE <- Sys.getenv("FLU_ANN_FILE", "/home/amir/repo/flu-cellular-prediction/results/reports/cluster_annotations_corrected.csv")
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

OUT_DIR <- file.path(FLU_ROOT, "post_hoc_5000/per_cohort")
FIG_DIR <- "/home/amir/repo/flu-cellular-prediction/paper/images"
load(file.path(OUT_DIR, "E5_per_cohort.RData"))

ANN_FILE <- FLU_ANN_FILE
ann <- fread(ANN_FILE); ann[, cluster_id := as.integer(cluster_id)]
short_label <- function(annotation) {
  out <- gsub("\\s*\\(.*?\\)", "", annotation)
  trimws(gsub("\\s+", " ", out))
}
ann[, short := short_label(annotation)]

# ---------------------------------------------------------------------------
# Pooled-cohort reference values
# (from the recipient-clustered bootstrap run in earlier E3 work)
# ---------------------------------------------------------------------------
pooled <- fread(file.path(FLU_ROOT, "post_hoc_5000/recipient_bootstrap/E3_recipient_boot.csv"))[scope == "overall" & comparison %in% c("M1_to_M3", "M2_to_M4"), .(cohort="pooled",comparison,delta_pct,ci_low=ci_low_pct,ci_high=ci_high_pct,p_two_sided)]
per <- fread(file.path(OUT_DIR, "E5_summary.csv"))
summary_all <- rbind(per, pooled, fill = TRUE)
summary_all[, cohort := factor(cohort,
  levels = c("peds", "adult", "pooled"),
  labels = c("Pediatric", "Adult", "Pooled"))]
summary_all[, comparison := factor(comparison,
  levels = c("M1_to_M3", "M2_to_M4"),
  labels = c("M1 → M3", "M2 → M4"))]

cohort_colors <- c(
  "Pediatric" = "#2166ac",
  "Adult"     = "#d6604d",
  "Pooled"    = "#4d4d4d"
)

p_imp <- ggplot(summary_all, aes(x = cohort, y = delta_pct, fill = cohort)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.18, linewidth = 0.5) +
  geom_hline(yintercept = 0, color = "gray30", linewidth = 0.3) +
  geom_text(aes(label = sprintf("italic(p) == %.3f", p_two_sided),
                 y = ci_high + 1.5),
            parse = TRUE, size = 2.8, color = "gray20") +
  facet_wrap(~ comparison, nrow = 1) +
  scale_fill_manual(values = cohort_colors, guide = "none") +
  labs(x = NULL, y = "MSE improvement (%)") +
  theme_bw(base_size = 9.5) +
  theme(axis.text.x = element_text(angle = 0, size = 8),
        strip.text = element_text(face = "bold"),
        strip.background = element_rect(fill = "gray92"))

# ---------------------------------------------------------------------------
# Selection-frequency panel
# Show all cluster features (drop the 3 baseline titers since they are 100% in
# both cohorts and not informative for the cohort-comparison story).
# ---------------------------------------------------------------------------
sel <- sel_summary[!grepl("^baseline_", feature)]
sel[, cluster_id := as.integer(sub("cluster_", "", feature))]
sel <- merge(sel, ann[, .(cluster_id, short)], by = "cluster_id", all.x = TRUE)
sel[, label := sprintf("C%d: %s", cluster_id,
                       ifelse(is.na(short), as.character(cluster_id), short))]
# Order by max selection frequency across the two cohorts
sel[, max_freq := pmax(peds_freq, adult_freq)]
sel <- sel[order(-max_freq)]

# Long format for plotting
sel_long <- melt(sel[1:18, .(label, peds_freq, adult_freq)],
                  id.vars = "label",
                  variable.name = "cohort_var",
                  value.name = "freq")
sel_long[, cohort := factor(cohort_var,
  levels = c("peds_freq", "adult_freq"),
  labels = c("Pediatric", "Adult"))]
sel_long[, label := factor(label, levels = rev(sel$label[1:18]))]

p_sel <- ggplot(sel_long, aes(x = freq, y = label, fill = cohort)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65,
            color = "gray30", linewidth = 0.2) +
  geom_vline(xintercept = 0.5, color = "red", linetype = "dashed",
              linewidth = 0.3) +
  scale_fill_manual(values = c("Pediatric" = "#2166ac", "Adult" = "#d6604d"),
                     name = NULL) +
  scale_x_continuous(limits = c(0, 1.02), expand = expansion(0.01)) +
  labs(x = "Selection frequency", y = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 7))

# ---------------------------------------------------------------------------
# Stratified quartile panel
# ---------------------------------------------------------------------------
strat_all[, cohort_lbl := factor(cohort, levels = c("peds", "adult"),
                                   labels = c("Pediatric", "Adult"))]
p_strat <- ggplot(strat_all, aes(x = quartile, y = improvement_pct,
                                   fill = cohort_lbl)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
            color = "gray30", linewidth = 0.3) +
  geom_hline(yintercept = 0, color = "gray30", linewidth = 0.3) +
  geom_text(aes(label = sprintf("n=%d", n),
                 y = ifelse(improvement_pct >= 0,
                              improvement_pct + 1.5, improvement_pct - 1.5)),
            position = position_dodge(width = 0.7),
            size = 2.6, color = "gray25") +
  scale_fill_manual(values = c("Pediatric" = "#2166ac", "Adult" = "#d6604d"),
                     name = NULL) +
  labs(x = "Baseline titer quartile (within cohort)",
       y = "MSE improvement (%)") +
  theme_bw(base_size = 9.5) +
  theme(legend.position = "bottom")

# Combine: 3-panel figure
fig <- (p_imp / p_sel / p_strat) +
  plot_layout(heights = c(0.95, 1.4, 0.95)) +
  plot_annotation(tag_levels = "A", tag_suffix = ".")
fig <- fig & theme(plot.tag = element_text(face = "bold", size = 13))

ggsave(file.path(FIG_DIR, "E5_per_cohort_comparison.png"), fig,
       width = 9.5, height = 10, dpi = 300)
ggsave(file.path(FIG_DIR, "E5_per_cohort_comparison.pdf"), fig,
       width = 9.5, height = 10, device = cairo_pdf)
cat("Saved E5_per_cohort_comparison figure.\n")
