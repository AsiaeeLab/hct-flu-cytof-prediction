#!/usr/bin/env Rscript
# E14-07: Figure for the second-season validation (reads E14-05 outputs).
#   A. ROC for a fourfold composite rise in second-season episodes: first-season
#      models with titers alone (M1) or with immune-cell composition (M3).
#   B. Fourfold composite rise four weeks after the first dose: all first-season
#      recipients, and the returning recipients in each of their two seasons.
#   C. Each returning recipient's move between seasons on the probability map of
#      Figure 4C (first-season logistic model of a fourfold rise on baseline
#      titer and score), by response transition, with the median predicted
#      probability in each season.
#   D. Median HAI titer before and after vaccination in both seasons, by
#      response transition.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(pROC)})
ROOT <- Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/first_episode_2026-10-14")
RV <- file.path(ROOT, "post_hoc_5000/repeater_validation")
navy <- "#1B4F72"; orange <- "#E67E22"; grey <- "#8C8C8C"
th <- theme_minimal(base_size = 8, base_family = "sans") + theme(
  panel.grid.minor = element_blank(), axis.text = element_text(size = 6.5, color = "grey20"),
  axis.title = element_text(size = 7.5), strip.text = element_text(face = "bold", size = 7.5),
  plot.tag = element_text(face = "bold", size = 11), legend.position = "bottom", legend.title = element_blank(),
  legend.text = element_text(size = 7))
tg <- function(p, t) p + labs(tag = t)
d <- readRDS(file.path(RV, "repeater_episode_table.rds")); w <- readRDS(file.path(RV, "repeater_pairs.rds"))
val <- fread(file.path(RV, "year2_validation.csv"))
y2 <- d[season == "second"]
mk <- function(pred, nm) { r <- pROC::roc(y2$resp, y2[[pred]], quiet = TRUE, direction = "<")
  data.table(fpr = 1 - r$specificities, tpr = r$sensitivities, model = nm) }
rr <- rbind(mk("pred_M1", "Titers"), mk("pred_M3", "Titers + composition"))[order(model, fpr, tpr)]
pfmt <- function(p) if (p < 0.001) "p < 0.001" else sprintf("p = %.2g", p)
pA <- ggplot(rr, aes(fpr, tpr, color = model)) + geom_abline(linetype = "dotted", color = "grey60") +
  geom_step(linewidth = .6) +
  annotate("text", x = .98, y = .23, hjust = 1, size = 2.2, color = grey, label = sprintf("Titers: AUC %.2f", val$auc_M1)) +
  annotate("text", x = .98, y = .15, hjust = 1, size = 2.2, color = navy, label = sprintf("Titers + composition: AUC %.2f", val$auc_M3)) +
  annotate("text", x = .98, y = .07, hjust = 1, size = 2.2, label = sprintf("Error reduction %.0f%% (95%% CI, %.0f-%.0f%%)", val$mse_gain_pct, val$mse_gain_lo, val$mse_gain_hi)) +
  scale_color_manual(values = c(Titers = grey, `Titers + composition` = navy), guide = "none") +
  coord_equal(expand = FALSE) + labs(x = "False-positive rate", y = "True-positive rate") + th
lab <- c(`NR to NR` = "Non-responder, both seasons", `NR to R` = "Non-responder to responder",
         `R to NR` = "Responder to non-responder", `R to R` = "Responder, both seasons")
w[, grp := factor(lab[as.character(transition)], lab)]; w[, n := .N, by = grp]
w[, grpn := sprintf("%s (n = %d)", grp, n)]; lv <- unique(w[order(grp)]$grpn); w[, grpn := factor(grpn, lv)]
# C. Same model and color scale as Figure 4C, fitted to the first-season recipients
s1 <- d[season == "first"]
lg <- glm(resp ~ base_comp + z, binomial, s1)
w[, `:=`(p1 = predict(lg, data.frame(base_comp = base_comp_first, z = z_first), type = "response"),
         p2 = predict(lg, data.frame(base_comp = base_comp_second, z = z_second), type = "response"))]
XL <- log(c(4, 8000)); YL <- c(-2.8, 2.6)
stopifnot(all(range(c(w$base_comp_first, w$base_comp_second)) > XL[1]), all(range(c(w$base_comp_first, w$base_comp_second)) < XL[2]))
grid <- CJ(base_comp = seq(XL[1], XL[2], length.out = 150), z = seq(YL[1], YL[2], length.out = 150))
grid[, prob := predict(lg, grid, type = "response")]
gg <- rbindlist(lapply(lv, function(k) cbind(grid, grpn = k))); gg[, grpn := factor(grpn, lv)]
med <- w[, .(x1 = median(base_comp_first), y1 = median(z_first), x2 = median(base_comp_second), y2 = median(z_second),
             mp1 = median(p1), mp2 = median(p2)), by = grpn]
med[, lab := sprintf("Modeled chance %.2f → %.2f", mp1, mp2)]
fwrite(med, file.path(RV, "transition_probability_moves.csv")); print(med)
crimson <- "#B2182B"; tick <- c(5, 20, 80, 320, 1280, 5120)
pB <- ggplot() +
  geom_raster(data = gg, aes(base_comp, z, fill = prob), interpolate = TRUE, alpha = .9) +
  geom_contour(data = gg, aes(base_comp, z, z = prob), breaks = c(.2, .4, .6, .8), color = "white", linewidth = .35) +
  geom_segment(data = w, aes(base_comp_first, z_first, xend = base_comp_second, yend = z_second), color = "grey35",
               alpha = .55, linewidth = .3, arrow = arrow(length = unit(.04, "in"), type = "closed")) +
  geom_segment(data = med, aes(x1, y1, xend = x2, yend = y2), color = crimson, linewidth = 1.1,
               arrow = arrow(length = unit(.09, "in"), type = "closed")) +
  geom_label(data = med, aes(x = log(5.5), y = YL[2] - .35, label = lab), hjust = 0, size = 2.3, color = crimson,
             fill = "white", linewidth = 0, label.padding = unit(.12, "lines")) +
  facet_wrap(~grpn, nrow = 2) +
  scale_fill_gradient(low = "#F3F6F9", high = "#9DB9D3", limits = c(0, 1), breaks = c(0, .5, 1),
                      labels = function(x) paste0(100 * x, "%"), name = "Modeled chance\nof fourfold rise") +
  scale_x_continuous(breaks = log(tick), labels = tick, expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = XL, ylim = YL) +
  labs(x = "Baseline HAI titer (geometric mean of three antigens)", y = "Cellular composition score (SD)") + th +
  theme(panel.spacing = unit(.7, "lines"), legend.position = "right", legend.title = element_text(size = 6.5))
# D. Titer before and after vaccination in both seasons
st <- w[, .(v1 = median(exp(base_comp_first)), v2 = median(exp(base_comp_first + y_comp_first)),
            v3 = median(exp(base_comp_second)), v4 = median(exp(base_comp_second + y_comp_second))), by = grpn]
sl <- melt(st, id.vars = "grpn", variable.name = "visit", value.name = "titer"); sl[, x := as.integer(visit)]
fwrite(st, file.path(RV, "transition_titer_medians.csv")); print(st)
pC <- ggplot(sl, aes(x, titer, color = grpn)) + geom_hline(yintercept = 40, linetype = "dashed", color = "grey50", linewidth = .3) +
  geom_line(linewidth = .8) + geom_point(size = 1.8) +
  scale_y_continuous(trans = "log2", breaks = c(10, 20, 40, 80, 160, 320, 640, 1280)) +
  scale_x_continuous(breaks = 1:4, labels = c("First season,\nbefore", "First season,\nafter", "Second season,\nbefore", "Second season,\nafter")) +
  scale_color_manual(values = c("#7F7F7F", navy, orange, "#5DADE2")) +
  labs(x = NULL, y = "Median HAI titer") + th + theme(legend.position = "right", legend.text = element_text(size = 6.5))
# B. One-dose response: all first-season recipients, and the returning recipients
# with a one-dose titer in both seasons (the same recipients in both bars)
rid <- intersect(d[season == "first" & is.finite(y_v2), base_id], d[season == "second" & is.finite(y_v2), base_id])
od <- rbind(d[season == "first" & is.finite(y_v2), .(grp = "All recipients,\nfirst season", n = .N, k = sum(resp_v2))],
            d[season == "first" & base_id %in% rid, .(grp = "Returning recipients,\nfirst season", n = .N, k = sum(resp_v2))],
            d[season == "second" & base_id %in% rid, .(grp = "Returning recipients,\nsecond season", n = .N, k = sum(resp_v2))])
w2 <- function(k, n) { p <- k/n; z <- 1.96; c <- (p + z^2/(2*n)) / (1 + z^2/n); h <- z*sqrt(p*(1-p)/n + z^2/(4*n^2)) / (1 + z^2/n); list(c - h, c + h) }
od[, rate := k / n]; od[, c("lo", "hi") := w2(k, n)]; od[, grp := factor(grp, grp)]
od[, kind := factor(c("all", "ret", "ret"))]
fwrite(od, file.path(RV, "one_dose_three_groups.csv")); print(od)
pD <- ggplot(od, aes(grp, 100 * rate, fill = kind)) + geom_col(width = .55) +
  geom_errorbar(aes(ymin = 100 * lo, ymax = 100 * hi), width = .15, linewidth = .4, color = "grey25") +
  geom_text(aes(y = 100 * hi + 6, label = sprintf("%.0f%% (%d/%d)", 100 * rate, k, n)), size = 2.2) +
  annotate("segment", x = 1.75, xend = 3.25, y = 88, yend = 88, linewidth = .35, color = "grey30") +
  annotate("text", x = 2.5, y = 94, label = "Same recipients", size = 2.3, color = "grey20") +
  scale_fill_manual(values = c(all = "#9DB9D3", ret = navy), guide = "none") +
  scale_y_continuous(limits = c(0, 100)) + labs(x = NULL, y = "Fourfold rise after one dose (%)") + th +
  theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(size = 6))
fig <- wrap_elements(full = (tg(pA, "A.") | tg(pD, "B.")) + plot_layout(widths = c(1, 1.15))) /
       wrap_elements(full = tg(pB, "C.")) / wrap_elements(full = tg(pC, "D.")) + plot_layout(heights = c(1, 1.55, .8))
ggsave("paper/images/E14_repeater_validation.png", fig, width = 7.2, height = 9.2, dpi = 600, bg = "white")
cat("Figure written.\n")
