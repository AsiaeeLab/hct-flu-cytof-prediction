#!/usr/bin/env Rscript
# E16-04: Figures for the composite-outcome primary analysis.
#   E16_ladder.png            Figure 2: added predictive value, composite outcome
#   E16_predictors.png        Figure 3: standardized coefficients of the composite
#                             model with titers, clinical factors and composition
#   E16_cohort_selection.png  supplement: selection frequency by cohort
#   E16_cohort_quartiles.png  supplement: gains by baseline quartile in each cohort
# Plotting only; estimates come from E12-01, E14-02c (OUTCOME=composite) and E16-02.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
ROOT <- Sys.getenv("FLU_RESULTS_ROOT"); PH <- file.path(ROOT, "post_hoc_5000")
fig_dir <- "paper/images"; blue <- "#1B4F72"; orange <- "#E67E22"; grey <- "#777777"
coh_col <- c(Children = "#2166ac", Adults = "#d6604d")
th <- theme_minimal(base_size = 8, base_family = "sans") + theme(
  panel.grid.minor = element_blank(), axis.text = element_text(size = 6.5, color = "grey20"),
  axis.title = element_text(size = 7.5), strip.text = element_text(face = "bold", size = 8),
  plot.title = element_text(face = "bold", size = 8.5), plot.tag = element_text(face = "bold", size = 11),
  legend.position = "bottom", legend.title = element_blank())
savefig <- function(p, name, h) ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = 7.2, height = h, dpi = 600, bg = "white")
pfmt <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %s", formatC(p, digits = 2, format = "fg")))
ann <- fread(Sys.getenv("FLU_ANN_FILE"))
short <- function(x) trimws(gsub("\\s+", " ", gsub("\\s*\\(.*?\\)", "", x)))
ann[, short := gsub("CD38hi", "CD38-high", short(annotation))]

# ---- Figure 2: composite four-model comparisons ------------------------------------
ci <- fread(file.path(PH, "composite/composite_ladder_bootstrap.csv"))
key <- data.table(comparison = c("M1_to_M3", "M2_to_M4", "M1_to_M2", "M3_to_M4", "M1_to_M4", "M2_to_M3"),
  label = c("Add immune composition to titers (M1 → M3)", "Add immune composition to titers + clinical factors (M2 → M4)",
            "Add clinical factors to titers (M1 → M2)", "Add clinical factors to titers + immune composition (M3 → M4)",
            "Add clinical factors and immune composition to titers (M1 → M4)", "Titers + immune composition vs titers + clinical factors (M2 → M3)"))
cc <- merge(key, ci, by = "comparison", sort = FALSE); stopifnot(nrow(cc) == 6)
cc[, label := factor(label, levels = rev(key$label))]; cc[, primary := comparison %in% c("M1_to_M3", "M2_to_M4")]
p2 <- ggplot(cc, aes(delta_pct, label, color = primary)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = .35) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), linewidth = .7) + geom_point(size = 2.5) +
  geom_text(aes(x = 27, label = sprintf("%.1f%% (%.1f to %.1f), %s", delta_pct, ci_low, ci_high, pfmt(p_two_sided))),
            hjust = 0, size = 2.2, color = "grey20") +
  scale_color_manual(values = c(`TRUE` = blue, `FALSE` = grey), guide = "none") +
  scale_x_continuous(breaks = seq(-5, 25, 5)) + coord_cartesian(xlim = c(-6, 26), clip = "off") +
  labs(x = "Reduction in mean squared prediction error (%)", y = NULL) + th +
  theme(panel.grid.major.y = element_blank(), plot.margin = margin(5, 120, 5, 5))
savefig(p2, "E16_ladder", 3.0)

# ---- Figure 3: standardized coefficients of the composite M4 -----------------------
a <- new.env(); load(Sys.getenv("FLU_PROD_ATLAS"), a)
bt <- new.env(); load(Sys.getenv("FLU_BOOT_FILE_COMP"), bt)
cd <- as.data.table(bt$coef_dt)[feature != "(Intercept)"]
sdv <- apply(a$X_clust, 2, sd)
agg <- cd[, .(mean_coef = mean(coef), sel = mean(coef != 0)), by = feature]
agg[, imp := mean_coef * sdv[feature]]; stopifnot(!anyNA(agg$imp))
agg[, is_cluster := grepl("^cluster_", feature)]
agg[, cid := ifelse(is_cluster, as.integer(sub("cluster_", "", feature)), NA_integer_)]
cov_pretty <- c(baseline_h1n1 = "Baseline H1N1 titer", baseline_h3n2 = "Baseline H3N2 titer", baseline_bvic = "Baseline B/Victoria titer",
                dose_HD = "High-dose vaccine", age = "Age", sex_Male = "Male sex", time_post_transplant = "Months since transplantation")
agg[, label := ifelse(is_cluster, sprintf("C%d: %s", cid, ann$short[match(cid, ann$cluster_id)]), cov_pretty[feature])]
agg[, lineage := ifelse(is_cluster, ann$broad_lineage[match(cid, ann$cluster_id)], "Clinical")]
fwrite(agg[order(-abs(imp))], file.path(PH, "composite_main/composite_importance.csv"))
cl <- agg[is_cluster == TRUE & sel >= 0.5][order(imp)]
cl[, label := factor(label, levels = label)]
cv <- agg[is_cluster == FALSE][order(imp)]; cv[, label := factor(label, levels = label)]
bar <- function(d, title) ggplot(d, aes(imp, label, fill = imp > 0)) + geom_col(width = .72) +
  geom_vline(xintercept = 0, linewidth = .3) +
  scale_fill_manual(values = c(`TRUE` = orange, `FALSE` = blue), breaks = c("TRUE", "FALSE"),
                    labels = c("Higher predicted fold-rise", "Lower predicted fold-rise")) +
  labs(x = "Standardized coefficient (log fold-rise per SD)", y = NULL, title = title) + th +
  theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 6.5))
p3 <- (bar(cv, "Baseline Titers and Clinical Factors") / bar(cl, "Immune-Cell Populations")) +
  plot_layout(heights = c(nrow(cv), nrow(cl)) + 1.5, guides = "collect") +
  plot_annotation(tag_levels = "A", tag_suffix = ".") & theme(legend.position = "bottom")
savefig(p3, "E16_predictors", 1.2 + 0.17 * (nrow(cv) + nrow(cl)) + 0.8)
cat("Populations selected in >= 50% of refits:", nrow(cl), "of", sum(agg$is_cluster), "\n")
print(agg[is_cluster == TRUE & sel >= 0.5][order(-imp), .(cid, label, sel = round(sel, 3), imp = round(imp, 3))])
print(agg[is_cluster == FALSE, .(feature, sel = round(sel, 3), imp = round(imp, 3))])

# ---- Supplement: selection frequency by cohort --------------------------------------
sw <- fread(file.path(PH, "composite_main/per_cohort_selection.csv"))[grepl("^cluster_", feature)]
sw[, cid := as.integer(sub("cluster_", "", feature))]
sw[, label := sprintf("C%d: %s", cid, ann$short[match(cid, ann$cluster_id)])]
sw <- sw[pmax(freq_Children, freq_Adults) >= 0.5][order(freq_Children + freq_Adults)]
sw[, label := factor(label, levels = label)]
sl <- melt(sw, id.vars = "label", measure.vars = c("freq_Children", "freq_Adults"), variable.name = "cohort", value.name = "freq")
sl[, cohort := factor(sub("freq_", "", cohort), c("Children", "Adults"))]
pS <- ggplot(sl, aes(freq, label, color = cohort)) +
  geom_vline(xintercept = .7, linetype = "dotted", color = "grey55") +
  geom_line(aes(group = label), color = "grey80", linewidth = .5) + geom_point(size = 2) +
  scale_color_manual(values = coh_col) + scale_x_continuous(limits = c(0, 1), labels = function(x) paste0(100 * x, "%")) +
  labs(x = "Selected in bootstrap refits", y = NULL) + th + theme(panel.grid.major.y = element_blank())
savefig(pS, "E16_cohort_selection", 0.6 + 0.16 * nrow(sw) + 0.6)
fwrite(sw, file.path(PH, "composite_main/per_cohort_selection_shown.csv"))

# ---- Supplement: gains by baseline quartile within each cohort ----------------------
pq <- fread(file.path(PH, "composite_main/per_cohort_quartiles.csv"))
pq[, cohort := factor(cohort, c("Children", "Adults"))]
pQ <- ggplot(pq, aes(quartile, improvement_pct, fill = cohort)) + geom_hline(yintercept = 0, color = "grey40", linewidth = .3) +
  geom_col(position = position_dodge(.7), width = .62) +
  geom_text(aes(y = ifelse(improvement_pct >= 0, improvement_pct + 1.5, improvement_pct - 1.5), label = paste0("n=", n)),
            position = position_dodge(.7), size = 2.1) +
  scale_fill_manual(values = coh_col) + labs(x = "Baseline-titer quartile within cohort", y = "Reduction in prediction error (%)") + th +
  theme(panel.grid.major.x = element_blank())
savefig(pQ, "E16_cohort_quartiles", 3.2)
cat("Done.\n")
