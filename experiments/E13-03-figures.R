#!/usr/bin/env Rscript
# Figures for the clinical reading of the cellular composition score
# (October 2026 round). Presentation only: reads the outputs of E13-02, E5 and
# E12-05; nothing is refitted. Run from the repository root.
#
#   E13_score_clinical.png       main: what the score measures and what it buys
#   E13_cohort_similarity.png    main: children and adults share the profile
#   E13_score_continuous.png     supplement: score versus fold-rise, continuous
#   E13_time_and_protection.png  supplement: response by score within time bands;
#                                protection among recipients unprotected at baseline
#   E13_score_populations.png    supplement: every population, lowest vs highest fifth
FLU_ROOT <- Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06")
FLU_BIO_FILE <- Sys.getenv("FLU_BIO_FILE", "/home/amir/repo/flu-cellular-prediction/results/final_atlas_production/biology_production_90203bf3.RData")
FLU_ANN_FILE <- Sys.getenv("FLU_ANN_FILE", "/home/amir/repo/flu-cellular-prediction/results/reports/cluster_annotations_corrected.csv")
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(pROC); library(scales)})
fig_dir <- "paper/images"
ROOT <- file.path(FLU_ROOT, "post_hoc_5000")
CS <- file.path(ROOT, "clinical_score")
navy <- "#1B4F72"; orange <- "#E67E22"; grey <- "#8C8C8C"
coh_col <- c(Children = "#2166ac", Adults = "#d6604d")
th <- theme_minimal(base_size = 8, base_family = "sans") + theme(
  panel.grid.minor = element_blank(), axis.text = element_text(size = 6.5, color = "grey20"),
  axis.title = element_text(size = 7.5), strip.text = element_text(face = "bold", size = 8),
  plot.title = element_text(face = "bold", size = 8.5), plot.tag = element_text(face = "bold", size = 11),
  legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 7))
savefig <- function(p, name, h) ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = 7.2, height = h, dpi = 600, bg = "white")
pfmt <- function(p) ifelse(p < 0.001, "italic(p) < 0.001", sprintf("italic(p) == %s", formatC(p, format = "fg", digits = 2)))
pfmt2 <- function(p) ifelse(p < 0.001, 'italic(p)*" < 0.001"', ifelse(p <= 0.0010001, 'italic(p)<=0.001', sprintf('italic(p)*" = %s"', formatC(p, format = "fg", digits = 2))))
pfmt_txt <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %s", formatC(p, format = "fg", digits = 2)))

d <- readRDS(file.path(CS, "episode_table.rds"))
pnull <- if (file.exists(file.path(CS, "null_pvalues.csv"))) fread(file.path(CS, "null_pvalues.csv")) else NULL
getp <- function(stat, sub = "All", fallback) {
  if (is.null(pnull)) return(fallback)
  hit <- pnull[pnull$statistic == stat & pnull$subset == sub]
  if (nrow(hit)) hit$p_null[1] else fallback }
fifth_lab <- c(F1 = "1\n(lowest)", F2 = "2", F3 = "3", F4 = "4", F5 = "5\n(highest)")

# ---- Figure 4: the score in clinical terms -----------------------------------
# A. What a low score looks like: population frequencies, lowest vs highest fifth
bio <- new.env(); load(FLU_BIO_FILE, bio)
props <- bio$meta_props_all; rownames(props) <- bio$common_p_ids
ann <- fread(FLU_ANN_FILE)
short <- function(x) gsub("CD38hi", "CD38-high", trimws(sub("\\s*\\(.*\\)$", "", x)), fixed = TRUE)
lin_lab <- c(`B cell` = "B cells", `CD4+ T` = "CD4+ T cells", `CD8+ T` = "CD8+ T cells",
             `Other T` = "Double-negative T cells", `NK cell` = "NK-like cells",
             `APC/Monocyte` = "Monocytes and other APCs", Other = "Lineage-negative populations")
lin <- ann$broad_lineage[match(seq_len(ncol(props)), ann$cluster_id)]
L <- sapply(names(lin_lab), function(l) rowSums(props[d$episode_id, paste0("cluster_", which(lin == l)), drop = FALSE])) * 100
an <- rbindlist(lapply(names(lin_lab), function(l) data.table(lineage = lin_lab[[l]],
  low = median(L[d$fifth == "F1", l]), high = median(L[d$fifth == "F5", l]))))
an[, ratio := high / low]
an <- an[order(-ratio)]; an[, label := factor(lineage, rev(lineage))]
anl <- melt(an, id.vars = c("lineage","label","ratio"), measure.vars = c("low","high"), variable.name = "grp", value.name = "pct")
anl[, grp := factor(ifelse(grp == "low", "Lowest-fifth scores", "Highest-fifth scores"), c("Lowest-fifth scores","Highest-fifth scores"))]
an[, rlab := sprintf("%.1f×", ifelse(ratio >= 1, ratio, 1 / ratio))]
pA <- ggplot(an) +
  geom_segment(aes(x = low, xend = high, y = label, yend = label), color = "grey70", linewidth = .7) +
  geom_point(data = anl, aes(pct, label, color = grp), size = 2.4) +
  geom_text(aes(x = 75, y = label, label = rlab, color = ifelse(ratio >= 1, "Highest-fifth scores", "Lowest-fifth scores")),
            size = 2.4, hjust = 0, show.legend = FALSE, fontface = "bold") +
  scale_x_log10(labels = function(x) paste0(format(x, drop0trailing = TRUE, scientific = FALSE), "%"),
                breaks = c(.1, .3, 1, 3, 10, 30)) +
  coord_cartesian(xlim = c(.08, 60), clip = "off") +
  scale_color_manual(values = c("Lowest-fifth scores" = orange, "Highest-fifth scores" = navy)) +
  labs(x = "Median share of PBMCs (log scale)", y = NULL) + th +
  theme(panel.grid.major.y = element_blank(), legend.position = "top", legend.box.margin = margin(0, 0, -6, 0),
        plot.margin = margin(5, 24, 5, 5))

# Supplement: every population, lowest versus highest fifth
NPOP <- ncol(props)
pop <- rbindlist(lapply(seq_len(NPOP), function(j) { x <- props[d$episode_id, paste0("cluster_", j)] * 100
  data.table(label = sprintf("C%d: %s", j, short(ann$annotation[ann$cluster_id == j])), lineage = lin_lab[[lin[j]]],
             low = median(x[d$fifth == "F1"]), high = median(x[d$fifth == "F5"])) }))
pop[, ratio := high / low]; pop <- pop[order(-ratio)]; pop[, label := factor(label, rev(label))]
popl <- melt(pop, id.vars = c("label","lineage","ratio"), measure.vars = c("low","high"), variable.name = "grp", value.name = "pct")
popl[, grp := factor(ifelse(grp == "low", "Lowest-fifth scores", "Highest-fifth scores"), c("Lowest-fifth scores","Highest-fifth scores"))]
pPop <- ggplot(pop) + geom_segment(aes(x = low, xend = high, y = label, yend = label), color = "grey70", linewidth = .5) +
  geom_point(data = popl, aes(pct, label, color = grp), size = 1.7) +
  geom_text(aes(x = 60, y = label, label = sprintf("%.1f×", ifelse(ratio >= 1, ratio, 1 / ratio)),
                color = ifelse(ratio >= 1, "Highest-fifth scores", "Lowest-fifth scores")), size = 2, hjust = 1, show.legend = FALSE) +
  scale_x_log10(labels = function(x) paste0(format(x, drop0trailing = TRUE, scientific = FALSE), "%"), breaks = c(.03, .1, .3, 1, 3, 10)) +
  coord_cartesian(xlim = c(.03, 60), clip = "off") +
  scale_color_manual(values = c("Lowest-fifth scores" = orange, "Highest-fifth scores" = navy)) +
  labs(x = "Median share of PBMCs (log scale)", y = NULL) + th + theme(panel.grid.major.y = element_blank(), legend.position = "top")
savefig(pPop, "E13_score_populations", 7.6)
fwrite(rbind(cbind(level = "lineage", an[, .(name = lineage, low, high, ratio)]), cbind(level = "population", pop[, .(name = as.character(label), low, high, ratio)])),
       file.path(CS, "composition_by_score_fifth.csv"))

# B. Fourfold response by fifth of the score (all episodes)
ff <- fread(file.path(CS, "response_by_fifth.csv"))
fa <- ff[cohort == "All"]; fa[, lab := fifth_lab[fifth]]
bt <- fread(file.path(CS, "bottom_fifth.csv"))[subset == "All"]
overall <- mean(d$resp_comp)
p_or <- getp("logistic_or_per_sd", "All", fread(file.path(CS, "logistic_models.csv"))[model == "baseline + score", p])
pB <- ggplot(fa, aes(fifth, 100 * rate)) +
  geom_col(fill = ifelse(fa$fifth == "F1", orange, navy), width = .62) +
  geom_errorbar(aes(ymin = 100 * lo, ymax = 100 * hi), width = .18, linewidth = .4, color = "grey25") +
  geom_text(aes(y = 100 * hi + 4, label = sprintf("%.0f%%", 100 * rate)), size = 2.4) +
  geom_text(aes(y = 3, label = sprintf("%d/%d", responders, n)), size = 1.9, color = "white") +
  geom_hline(yintercept = 100 * overall, linetype = "dashed", color = "grey50", linewidth = .35) +
  scale_x_discrete(labels = fifth_lab) + scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, .02))) +
  labs(x = "Cellular composition score (fifths)", y = "Recipients with fourfold composite rise (%)") + th +
  theme(panel.grid.major.x = element_blank())

# C. Response map: baseline titer by score, with modeled probability
lg <- glm(resp_comp ~ base_comp + z_comp, binomial, d)
grid <- CJ(base_comp = seq(0, log(6500), length.out = 150), z_comp = seq(-3.4, 2.5, length.out = 150))
grid[, prob := predict(lg, grid, type = "response")]
d[, status := factor(ifelse(resp_comp, "Fourfold rise", "Less than fourfold rise"), c("Fourfold rise","Less than fourfold rise"))]
tick <- c(1, 5, 20, 80, 320, 1280, 5120)
pC <- ggplot() +
  geom_raster(data = grid, aes(base_comp, z_comp, fill = prob), interpolate = TRUE, alpha = .9) +
  geom_contour(data = grid, aes(base_comp, z_comp, z = prob), breaks = c(.2, .4, .6, .8), color = "white", linewidth = .35) +
  # fourfold rise: filled, area scaled by the composite fold-rise; below fourfold: open, fixed size
  geom_point(data = d[resp_comp == FALSE], aes(base_comp, z_comp, shape = "Less than fourfold rise"), color = orange, size = 1.1, stroke = .45) +
  geom_point(data = d[resp_comp == TRUE], aes(base_comp, z_comp, size = exp(y_comp)), shape = 21,
             fill = scales::alpha(navy, .3), color = scales::alpha(navy, .75), stroke = .3) +
  scale_shape_manual(values = c(`Less than fourfold rise` = 1), name = NULL) +
  scale_size(range = c(1.1, 3.4), trans = "log2", breaks = c(4, 16, 64, 256), name = "Fold-rise,\nif fourfold or more") +
  scale_fill_gradient(low = "#F3F6F9", high = "#9DB9D3", limits = c(0, 1), breaks = c(0, .5, 1), labels = percent, name = "Modeled chance\nof fourfold rise") +
  scale_x_continuous(breaks = log(tick), labels = tick, expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Baseline HAI titer (geometric mean of three antigens)", y = "Cellular composition score (SD)") + th +
  guides(fill = guide_colorbar(barheight = unit(2.6, "mm"), barwidth = unit(22, "mm"), title.position = "left",
                               title.theme = element_text(size = 6.2)),
         size = guide_legend(title.position = "left", title.theme = element_text(size = 6.2), override.aes = list(color = navy)),
         shape = guide_legend(override.aes = list(color = orange, size = 1.6))) +
  theme(legend.box = "vertical", legend.spacing.y = unit(0, "mm"), legend.margin = margin(0, 0, 0, 0))
cl <- rbindlist(lapply(c(.2, .4, .6, .8), function(p) {
  zs <- seq(-3, 2.2, by = .05); bs <- (qlogis(p) - coef(lg)[1] - coef(lg)[3] * zs) / coef(lg)[2]
  ok <- which(bs > log(2) & bs < log(4000)); if (!length(ok)) return(NULL)
  i <- ok[which.min(abs(bs[ok] - log(60)))]; data.table(lbl = sprintf("%d%%", 100 * p), base_comp = bs[i], z_comp = zs[i]) }))
pC <- pC + geom_label(data = cl, aes(base_comp, z_comp, label = lbl), size = 1.9, linewidth = 0, fill = "white", alpha = .85,
                      label.padding = unit(.06, "lines"))

# D. ROC: titers alone versus titers plus composition
mk_roc <- function(pred, nm) { r <- pROC::roc(d$resp_comp, d[[pred]], quiet = TRUE, direction = "<")
  data.table(fpr = 1 - r$specificities, tpr = r$sensitivities, model = nm, auc = as.numeric(r$auc)) }
rr <- rbind(mk_roc("pred_M1", "Titers"), mk_roc("pred_M3", "Titers + composition"))
rt <- fread(file.path(CS, "roc.csv"))[subset == "All"]
p_auc <- getp("delta_auc_M1_M3", "All", rt[quantity == "Delta AUC M1->M3", p_delong])
rr <- rr[order(model, fpr, tpr)]
aucs <- rr[, .(auc = auc[1]), by = model]
aucs[, lab := sprintf("%s: AUC %.2f", model, auc)]
pD <- ggplot(rr, aes(fpr, tpr, color = model)) + geom_abline(linetype = "dotted", color = "grey60") +
  geom_step(linewidth = .6, direction = "hv") +
  annotate("text", x = .98, y = .16, hjust = 1, size = 2.3, label = aucs$lab[1], color = grey) +
  annotate("text", x = .98, y = .09, hjust = 1, size = 2.3, label = aucs$lab[2], color = navy) +
  annotate("text", x = .98, y = .02, hjust = 1, size = 2.3, parse = TRUE,
           label = sprintf('Delta*"AUC = %.2f, "*%s', rt[quantity == "Delta AUC M1->M3", estimate], pfmt2(p_auc))) +
  scale_color_manual(values = c(Titers = grey, `Titers + composition` = navy), guide = "none") +
  coord_equal(expand = FALSE) + labs(x = "False-positive rate", y = "True-positive rate") + th

tg <- function(p, t) p + labs(tag = t)
fig4 <- wrap_elements(full = (tg(pA, "A.") | tg(pB, "B.")) + plot_layout(widths = c(1.2, 1))) /
        wrap_elements(full = (tg(pC, "C.") | tg(pD, "D.")) + plot_layout(widths = c(1.25, 1))) + plot_layout(heights = c(1, 1.08))
savefig(fig4, "E13_score_clinical", 7.0)

# ---- Figure 5: children and adults ----------------------------------------------
cg <- fread(file.path(CS, "pooled_model_gain_by_cohort.csv"))[subset %in% c("Children","Adults")]
cg[, model := "Pooled model"]
e5 <- fread(file.path(ROOT, "per_cohort/E5_summary.csv"))
e5[, `:=`(subset = ifelse(cohort == "peds", "Children", "Adults"), comparison = sub("_to_", "->", comparison),
          estimate = delta_pct, p_boot = p_two_sided, model = "Model fitted within cohort")]
g <- rbind(cg[, .(subset, comparison, estimate, ci_low, ci_high, p_boot, model)], e5[, .(subset, comparison, estimate, ci_low, ci_high, p_boot, model)])
# Composite-outcome primary analysis (E16-02): pooled models within each cohort
# and models refitted within each cohort, both on the composite outcome.
CM <- file.path(ROOT, "composite_main")
if (file.exists(file.path(CM, "pooled_within_cohort.csv"))) {
  pw <- fread(file.path(CM, "pooled_within_cohort.csv"))[cohort %in% c("Children", "Adults")]
  rf <- fread(file.path(CM, "per_cohort_refit.csv"))
  g <- rbind(pw[, .(subset = cohort, comparison, estimate = gain, ci_low = lo, ci_high = hi, p_boot = p, model = "Pooled model")],
             rf[, .(subset = cohort, comparison, estimate = gain, ci_low = lo, ci_high = hi, p_boot = p, model = "Model fitted within cohort")])
}
g[, subset := factor(subset, c("Children","Adults"))]; g[, model := factor(model, c("Pooled model","Model fitted within cohort"))]
g[, comparison := factor(ifelse(comparison == "M1->M3", "Added to baseline titers", "Added to titers and clinical factors"),
                        c("Added to baseline titers", "Added to titers and clinical factors"))]
pd <- position_dodge(.7)
p5A <- ggplot(g, aes(subset, estimate, fill = subset, alpha = model, group = model)) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = .3) +
  geom_col(position = pd, width = .62, color = NA) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), position = pd, width = .18, linewidth = .35, alpha = 1, color = "grey25") +
  geom_text(aes(y = ci_high + 1.8, label = pfmt_txt(p_boot)), position = pd, size = 2.0, alpha = 1, fontface = "italic") +
  facet_wrap(~comparison) + scale_fill_manual(values = coh_col, guide = "none") +
  scale_alpha_manual(values = c(`Pooled model` = 1, `Model fitted within cohort` = .45)) +
  labs(x = NULL, y = "Reduction in prediction error (%)") + th + theme(panel.grid.major.x = element_blank())
fc <- ff[cohort %in% c("Children","Adults")]; fc[, cohort := factor(cohort, c("Children","Adults"))]
pd2 <- position_dodge(.45)
lgm <- fread(file.path(CS, "logistic_models.csv"))
p_int <- lgm[model == "baseline + score x cohort", p]
p5B <- ggplot(fc, aes(fifth, 100 * rate, color = cohort, group = cohort)) +
  geom_line(position = pd2, linewidth = .45) +
  geom_errorbar(aes(ymin = 100 * lo, ymax = 100 * hi), position = pd2, width = .2, linewidth = .4) +
  geom_point(position = pd2, size = 1.8) +
  annotate("text", x = 5.4, y = 4, hjust = 1, size = 2.2, label = sprintf("Score-by-cohort\ninteraction, %s", pfmt_txt(p_int))) +
  scale_color_manual(values = coh_col) + scale_x_discrete(labels = fifth_lab) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Cellular composition score (fifths)", y = "Fourfold composite rise (%)") + th
fo <- rbind(
  lgm[model == "baseline + score", .(lab = "All recipients", or, or_lo, or_hi, p)],
  lgm[model == "baseline + score, Children", .(lab = "Children", or, or_lo, or_hi, p)],
  lgm[model == "baseline + score, Adults", .(lab = "Adults", or, or_lo, or_hi, p)],
  lgm[model == "baseline + score + months + dose + cohort + age" & term == "z_comp", .(lab = "All, adjusted for months since\ntransplant, dose, age and cohort", or, or_lo, or_hi, p)])
fe <- fread(file.path(CS, "logistic_first_episode.csv"))
# In the first-season analysis this subset repeats the all-episode row.
if (abs(fe$or[1] - fo$or[1]) > 1e-8)
  fo <- rbind(fo[1:3], fe[, .(lab = "First episodes only", or, or_lo, or_hi, p)], fo[4])
fo[, lab := factor(lab, rev(lab))]
fo[, col := ifelse(lab == "Children", coh_col["Children"], ifelse(lab == "Adults", coh_col["Adults"], navy))]
XMAX5C <- max(4.2, max(fo$or_hi) * 1.05)  # whiskers end before the labels
p5C <- ggplot(fo, aes(or, lab)) + geom_vline(xintercept = 1, color = "grey50", linewidth = .3) +
  geom_errorbar(aes(xmin = or_lo, xmax = or_hi), orientation = "y", width = .18, linewidth = .45, color = fo$col) +
  geom_point(size = 2, color = fo$col) +
  geom_text(aes(x = XMAX5C * 1.04, label = sprintf("%.2f (%.2f–%.2f), %s", or, or_lo, or_hi, pfmt_txt(p))), hjust = 0, size = 2.1) +
  scale_x_log10(breaks = c(1, 1.5, 2, 3, 4, 6)) + coord_cartesian(xlim = c(.9, XMAX5C), clip = "off") +
  labs(x = "Odds ratio for fourfold rise per SD of score", y = NULL) + th +
  theme(panel.grid.major.y = element_blank(), plot.margin = margin(5, 92, 5, 5))
fig5 <- wrap_elements(full = tg(p5A, "A.")) / wrap_elements(full = (tg(p5B, "B.") | tg(p5C, "C.")) + plot_layout(widths = c(1, 1))) +
  plot_layout(heights = c(1, 1.15))
savefig(fig5, "E13_cohort_similarity", 5.6)

# ---- Supplement: continuous score versus fold-rise -------------------------------
rho <- fread(file.path(CS, "rank_correlations.csv"))
# every first-season recipient contributes one episode, so the third panel pools both cohorts
panels <- list(Children = d[cohort == "Children"], Adults = d[cohort == "Adults"], `Both cohorts` = d[first_episode == TRUE])
sc <- rbindlist(lapply(names(panels), function(n) cbind(panels[[n]][, .(z_comp, y_comp, cohort)], panel = n)))
sc[, panel := factor(panel, names(panels))]
rl <- rho[statistic == "Spearman" & subset %in% c("Children", "Adults", "First episodes")]
rl[, p_show := mapply(function(s, fb) getp("spearman", s, fb), subset, p_boot)]
rl[, lab := sprintf('rho*" = %.2f, "*%s', estimate, pfmt2(p_show))]; rl[, panel := factor(ifelse(subset == "First episodes", "Both cohorts", subset), names(panels))]
sc[, col := ifelse(panel == "Both cohorts", "First", as.character(cohort))]
pS1 <- ggplot(sc, aes(z_comp, exp(y_comp))) +
  geom_hline(yintercept = 4, linetype = "dashed", color = "grey50", linewidth = .35) +
  geom_point(aes(color = col), size = .9, alpha = .6) +
  geom_smooth(method = "loess", formula = y ~ x, span = .9, color = "black", fill = "grey80", linewidth = .6) +
  geom_text(data = rl, aes(x = -3, y = 900, label = lab), parse = TRUE, hjust = 0, size = 2.4) +
  scale_y_log10(breaks = c(1/16, .25, 1, 4, 16, 64, 256, 1024), labels = c("1/16","1/4","1","4","16","64","256","1024")) +
  scale_color_manual(values = c(coh_col, First = navy), guide = "none") +
  facet_wrap(~panel) + labs(x = "Cellular composition score (SD)", y = "Composite fold-rise") + th
# per-antigen response by fifth
pa <- rbindlist(lapply(c("h1n1","h3n2","bvic"), function(v) { zz <- d[[paste0("z_", v)]]
  q <- cut(zz, quantile(zz, 0:5/5), include.lowest = TRUE, labels = paste0("F", 1:5))
  d[, .(antigen = v, rate = mean(get(paste0("resp_", v))), n = .N, k = sum(get(paste0("resp_", v)))), by = .(fifth = q)] }))
pa[, antigen := factor(c(h1n1 = "H1N1", h3n2 = "H3N2", bvic = "B/Victoria")[antigen], c("H1N1","H3N2","B/Victoria"))]
w <- function(k, n) { p <- k/n; z <- 1.96; c <- (p + z^2/(2*n)) / (1 + z^2/n); h <- z*sqrt(p*(1-p)/n + z^2/(4*n^2)) / (1 + z^2/n); list(c - h, c + h) }
pa[, c("lo","hi") := w(k, n)]
pS2 <- ggplot(pa, aes(fifth, 100 * rate)) + geom_col(fill = navy, width = .6) +
  geom_errorbar(aes(ymin = 100 * lo, ymax = 100 * hi), width = .18, linewidth = .35, color = "grey25") +
  geom_text(aes(y = 100 * hi + 5, label = sprintf("%.0f%%", 100 * rate)), size = 2.1) +
  facet_wrap(~antigen) + scale_x_discrete(labels = c(F1 = "1", F2 = "2", F3 = "3", F4 = "4", F5 = "5")) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Antigen-specific cellular composition score (fifths, lowest to highest)", y = "Recipients with fourfold rise (%)") + th +
  theme(panel.grid.major.x = element_blank())
savefig(pS1 / pS2 + plot_annotation(tag_levels = "A", tag_suffix = "."), "E13_score_continuous", 5.6)

# ---- Supplement: time since transplant ---------------------------------------------
d[, tgrp := factor(cut(months, c(0, 6, 12, 100), labels = c("Less than 6 months", "6 to 12 months", "More than 12 months")))]
d[, sgrp := factor(ifelse(fifth == "F1", "Lowest fifth", ifelse(fifth == "F5", "Highest fifth", "Middle three fifths")),
                   c("Lowest fifth", "Middle three fifths", "Highest fifth"))]
tt <- d[, .(n = .N, k = sum(resp_comp)), by = .(tgrp, sgrp)]
tt[, rate := k / n]; tt[, c("lo","hi") := w(k, n)]
p_adj <- lgm[model == "baseline + score + months + dose + cohort + age" & term == "z_comp", p]
pT <- ggplot(tt, aes(tgrp, 100 * rate, color = sgrp, group = sgrp)) +
  geom_errorbar(aes(ymin = 100 * lo, ymax = 100 * hi), position = position_dodge(.55), width = .2, linewidth = .4) +
  geom_point(position = position_dodge(.55), size = 2) +
  geom_text(aes(y = 100 * hi + 5, label = sprintf("%d/%d", k, n)), position = position_dodge(.55), size = 2, show.legend = FALSE) +
  scale_color_manual(values = c(`Lowest fifth` = orange, `Middle three fifths` = grey, `Highest fifth` = navy)) +
  scale_y_continuous(limits = c(0, 105)) +
  labs(x = "Time from transplantation to vaccination", y = "Recipients with fourfold composite rise (%)") + th


# ---- Supplement: seroprotection among the unprotected --------------------------------
sp <- rbindlist(lapply(c("h1n1","h3n2","bvic"), function(v) { x <- d[get(paste0("unprot0_", v)) == TRUE]
  x[, .(antigen = v, n = .N, k = sum(get(paste0("prot_", v)))), by = .(grp = ifelse(bottom_fifth, "Lowest fifth", "Other four fifths"))] }))
sp[, rate := k / n]; sp[, c("lo","hi") := w(k, n)]
sp[, antigen := factor(c(h1n1 = "H1N1", h3n2 = "H3N2", bvic = "B/Victoria")[antigen], c("H1N1","H3N2","B/Victoria"))]
sp[, grp := factor(grp, c("Lowest fifth","Other four fifths"))]
pP <- ggplot(sp, aes(antigen, 100 * rate, fill = grp)) +
  geom_col(position = position_dodge(.7), width = .62) +
  geom_errorbar(aes(ymin = 100 * lo, ymax = 100 * hi), position = position_dodge(.7), width = .18, linewidth = .35, color = "grey25") +
  geom_text(aes(y = 100 * hi + 5, label = sprintf("%.0f%% (%d/%d)", 100 * rate, k, n)), position = position_dodge(.7), size = 2.1) +
  scale_fill_manual(values = c(`Lowest fifth` = orange, `Other four fifths` = navy)) + scale_y_continuous(limits = c(0, 105)) +
  labs(x = NULL, y = "Reached HAI titer of at least 1:40 (%)") + th + theme(panel.grid.major.x = element_blank())
savefig(wrap_elements(full = tg(pT, "A.")) / wrap_elements(full = tg(pP, "B.")), "E13_time_and_protection", 5.8)
cat("Figures written.\n")

# ---- Supplement: how many populations carry the prediction ----------------------
tk <- fread(file.path(ROOT, "top_k_populations/top_k_mse.csv"))
tk[, scope := factor(ifelse(scope == "overall", "Three antigens (multitask)", "Composite outcome"),
                     c("Three antigens (multitask)", "Composite outcome"))]
KMAX <- max(tk$k)
ref <- tk[k == KMAX]
pK <- ggplot(tk, aes(k, delta_pct, color = scope, fill = scope)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = .3) +
  geom_hline(data = ref, aes(yintercept = delta_pct, color = scope), linetype = "dashed", linewidth = .35, show.legend = FALSE) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = .12, color = NA) +
  geom_line(linewidth = .6) + geom_point(size = 1.5) +
  scale_x_continuous(trans = "log2", breaks = c(1, 2, 3, 4, 6, 8, 12, 20, KMAX)) +
  scale_color_manual(values = c(`Three antigens (multitask)` = navy, `Composite outcome` = grey)) +
  scale_fill_manual(values = c(`Three antigens (multitask)` = navy, `Composite outcome` = grey)) +
  labs(x = "Number of top-ranked populations (ranked within each training fold)", y = "Reduction in prediction error\nbeyond baseline titers (%)") + th
savefig(pK, "E13_top_k", 3.0)
cat("Top-k figure written.\n")
