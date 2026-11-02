#!/usr/bin/env Rscript
# E15-03: Held-out prediction gain at each fixed clustering resolution
# (E15-01 sweep), with recipient-bootstrap 95% intervals. Descriptive: shows
# how the gain depends on the number of populations; it is not used to choose
# the number. Writes paper/images/E15_gain_by_k.png.
#   A. Composite outcome.  B. Three antigens (multitask).
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
ROOT <- Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/first_episode_2026-10-14")
OUT <- file.path(ROOT, "post_hoc_5000/resolution_sweep")
o <- readRDS(file.path(OUT, "outer_predictions_by_resolution.rds")); V <- c("h1n1", "h3n2", "bvic")
B <- 2000L; set.seed(20260728); bi <- replicate(B, sample.int(nrow(o$y), replace = TRUE))
g <- function(a, b, i) 100 * (1 - sum(b[i]) / sum(a[i]))
rows <- rbindlist(lapply(sort(unique(o$outer$res)), function(r) {
  d <- merge(merge(o$outer[res == r], o$ref, by = c("fold", "episode_id")), o$y, by = "episode_id")
  sm <- function(m) rowMeans(sapply(V, function(v) (d[[paste0(m, "_", v)]] - d[[paste0("y_", v)]])^2))
  sc <- function(m) (d[[paste0(m, "_comp")]] - d$y_comp)^2
  k <- median(d[, k[1], by = fold]$V1)
  rbindlist(lapply(list(list("Composite outcome", "Added to titers", sc("M1"), sc("M3")),
                        list("Composite outcome", "Added to titers and clinical factors", sc("M2"), sc("M4")),
                        list("Three antigens", "Added to titers", sm("M1"), sm("M3")),
                        list("Three antigens", "Added to titers and clinical factors", sm("M2"), sm("M4"))), function(x) {
    bs <- apply(bi, 2, function(i) g(x[[3]], x[[4]], i))
    data.table(res = r, k = k, outcome = x[[1]], comparison = x[[2]], gain = g(x[[3]], x[[4]], seq_along(x[[3]])),
               lo = quantile(bs, .025), hi = quantile(bs, .975)) }))
}))
fwrite(rows, file.path(OUT, "heldout_gain_by_fixed_resolution_ci.csv"))
cols <- c(`Added to titers` = "#1B4F72", `Added to titers and clinical factors` = "#E67E22")
th <- theme_minimal(base_size = 8, base_family = "sans") + theme(panel.grid.minor = element_blank(),
  axis.text = element_text(size = 6.5, color = "grey20"), axis.title = element_text(size = 7.5),
  plot.title = element_text(face = "bold", size = 8), plot.tag = element_text(face = "bold", size = 11),
  legend.position = "bottom", legend.title = element_blank())
pl <- function(oc, tag) ggplot(rows[outcome == oc], aes(k, gain, color = comparison, fill = comparison)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = .3) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .12, color = NA) +
  geom_line(linewidth = .6) + geom_point(size = 1.3) +
  scale_x_continuous(trans = "log2", breaks = c(5, 7, 10, 14, 20, 28, 35, 42)) +
  scale_color_manual(values = cols) + scale_fill_manual(values = cols) +
  labs(x = "Number of populations (median across folds)", y = "Held-out reduction in prediction error (%)",
       title = oc, tag = tag) + th
fig <- (pl("Composite outcome", "A.") | pl("Three antigens", "B.")) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave("paper/images/E15_gain_by_k.png", fig, width = 7.2, height = 3.4, dpi = 600, bg = "white")
print(rows[, .(k, outcome, comparison, gain = round(gain, 1), lo = round(lo, 1), hi = round(hi, 1))][order(outcome, comparison, k)])
