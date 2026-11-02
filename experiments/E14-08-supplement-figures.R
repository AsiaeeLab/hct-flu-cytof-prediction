#!/usr/bin/env Rscript
# E14-08: Supplementary figures for the first-season analysis that the older
# refresh scripts cannot rebuild directly (they assume 33 populations or read
# all-episode-only inputs). Adapted from E9-refresh-supporting.R.
#   figS4_stability_landscape.png, figS5_selection_frequency.png,
#   figS6_coef_distributions.png  (bootstrap stability on the shared atlas)
#   figS7_transfer_results.png    (training in one cohort, testing in the other)
#   E14_same_score.png            (responders vs non-responders at the same score)
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
ROOT <- Sys.getenv("FLU_RESULTS_ROOT")
out <- "paper/images"
e <- new.env(); load(Sys.getenv("FLU_BOOT_FILE"), e)
K <- max(as.integer(sub("cluster_", "", grep("^cluster_", e$stability$feature, value = TRUE))))
s <- as.data.table(e$stability)[grepl("^cluster_", feature)]
s[, cluster_id := as.integer(sub("cluster_", "", feature))]; s[, label := paste0("C", cluster_id)]
labs_v <- c(h1n1 = "H1N1", h3n2 = "H3N2", bvic = "B/Victoria", comp = "Composite outcome")
VORD <- intersect(c("h1n1", "h3n2", "bvic", "comp"), unique(e$stability$virus))
theme_set(theme_minimal(base_size = 8, base_family = "sans") +
  theme(panel.grid.minor = element_blank(), axis.text = element_text(size = 6.5),
        strip.text = element_text(face = "bold"), legend.position = "bottom"))
savefig <- function(p, name, h = 4.4) ggsave(file.path(out, name), p, width = 7.2, height = h, dpi = 600, bg = "white")
s[, virus := factor(virus, levels = VORD)]
p <- ggplot(s, aes(selection_freq, sign_consistency, size = abs(mean_coef))) +
  geom_vline(xintercept = .5, linetype = "dashed", color = "grey55") +
  geom_hline(yintercept = .8, linetype = "dotted", color = "grey55") +
  geom_point(color = "#1B4F72", alpha = .65) + facet_wrap(~virus, labeller = labeller(virus = labs_v)) +
  scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Selection frequency", y = "Sign consistency", size = "Absolute mean coefficient")
savefig(p, "figS4_stability_landscape.png", 3.4)
s[, cluster_f := factor(label, levels = paste0("C", 1:K))]
p <- ggplot(s[virus == if ("comp" %in% VORD) "comp" else "h3n2"], aes(cluster_f, selection_freq)) + geom_col(fill = "#1B4F72") +
  geom_hline(yintercept = .5, linetype = "dashed", color = "grey55") + scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Population", y = "Selection frequency") + theme(axis.text.x = element_text(angle = 60, hjust = 1))
savefig(p, "figS5_selection_frequency.png", 3.5)
selected <- unique(s[selection_freq >= .5, feature])
d <- as.data.table(e$coef_dt)[feature %in% selected]
d[, cluster_f := factor(sub("cluster_", "C", feature), levels = paste0("C", sort(as.integer(sub("cluster_", "", selected)))))]
d[, virus := factor(virus, levels = VORD)]
p <- ggplot(d, aes(cluster_f, coef)) + geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
  geom_boxplot(fill = "#85C1E9", outlier.size = .4, linewidth = .3) + coord_flip() +
  facet_wrap(~virus, labeller = labeller(virus = labs_v)) + labs(x = "Population", y = "M4 coefficient")
savefig(p, "figS6_coef_distributions.png", 5.0)
a <- new.env(); load(Sys.getenv("FLU_PROD_ATLAS"), a)
t <- as.data.table(a$transfer_dt)
CT <- file.path(ROOT, "post_hoc_5000/composite_main/transfer.csv")
if (file.exists(CT)) t <- fread(CT)[, .(direction, virus = "comp", improve_pct = -gain, lo = -hi, hi = -lo, n_train, n_test)]
fwrite(t, file.path(ROOT, "post_hoc_5000/transfer_first_episode.csv"))
t[, direction := factor(direction, levels = c("peds_to_adult", "adult_to_peds"), labels = c("Pediatric → adult", "Adult → pediatric"))]
t[, virus := factor(virus, levels = intersect(c("h1n1", "h3n2", "bvic", "comp"), unique(virus)))]
p <- ggplot(t, aes(virus, -improve_pct, fill = direction)) + geom_hline(yintercept = 0, color = "grey50") +
  geom_col(position = position_dodge(.7), width = .62) + scale_x_discrete(labels = labs_v) +
  (if ("lo" %in% names(t)) geom_errorbar(aes(ymin = -hi, ymax = -lo), position = position_dodge(.7), width = .15, linewidth = .35, color = "grey25") else NULL) +
  scale_fill_manual(values = c("#1B4F72", "#E67E22")) +
  labs(x = NULL, y = "Test MSE improvement over clinical baseline (%)", fill = NULL)
savefig(p, "figS7_transfer_results.png", 3.6)
print(t[, .(direction, virus, improvement = round(-improve_pct, 1), n_train, n_test)])
# Same-score comparison (E14-06)
sm <- fread(file.path(ROOT, "post_hoc_5000/same_score/smd_pooled_with_adjusted.csv"))
sm <- sm[order(-abs(smd))][1:20]
sm[, feature := gsub("CD38hi", "CD38-high", feature)]
sm[, feature := factor(feature, rev(feature))]
lg <- melt(sm, id.vars = "feature", measure.vars = c("smd", "smd_baseline_adjusted"), variable.name = "type", value.name = "v")
lg[, type := factor(type, c("smd", "smd_baseline_adjusted"), c("Unadjusted", "Adjusted for baseline titer"))]
p <- ggplot(lg[is.finite(v)], aes(v, feature, color = type, shape = type)) + geom_vline(xintercept = 0, color = "grey50", linewidth = .3) +
  geom_point(size = 1.8) + scale_color_manual(values = c("#1B4F72", "#E67E22")) +
  labs(x = "Standardized difference, responders minus non-responders, within score fifths", y = NULL, color = NULL, shape = NULL) +
  theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 6))
savefig(p, "E14_same_score.png", 4.6)
cat("Done.\n")
