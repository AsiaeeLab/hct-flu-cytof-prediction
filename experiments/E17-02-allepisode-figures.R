#!/usr/bin/env Rscript
# E17-02: Replot the two all-episode supplementary displays (supplemental Section on
# all episodes) from their saved September 6, 2026 results, with the paper's label
# conventions (bold "A." panel tags, "B/Victoria"). Plotting only; nothing refitted.
#   E2_multitask_vs_single.png   (from E2-single-task-comparison.R output)
#   E2_sensitivity_2000v5000.png (from E9-refresh-supporting.R inputs)
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
root <- Sys.getenv("FLU_ALLEP_ROOT", "/fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06")
out <- "paper/images"
th <- theme_minimal(base_size = 8, base_family = "sans") +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom", legend.title = element_blank(),
        plot.tag = element_text(face = "bold", size = 11))
savefig <- function(p, name, h) ggsave(file.path(out, name), p, width = 7.2, height = h, dpi = 600, bg = "white")
vlab <- c(h1n1 = "H1N1", h3n2 = "H3N2", bvic = "B/Victoria")

# ---- Joint versus separate antigen models ------------------------------------------
e <- new.env(); load(file.path(root, "per_fold_atlas_5000/aggregated/E2_single_task.RData"), e)
mse <- as.data.table(e$mse_df)
m1 <- mse[model == "M1_multi"]
ov <- mse[virus == "overall" & model != "M1_multi"]
ov[, `:=`(task = ifelse(grepl("multi", model), "Multitask", "Single-task"), model_base = sub("_.*", "", model),
          improvement = 100 * (1 - mse / m1[virus == "overall", mse]))]
pv <- merge(mse[model %in% c("M3_multi", "M3_single") & virus != "overall"], m1[virus != "overall", .(virus, mse_m1 = mse)], by = "virus")
pv[, `:=`(task = ifelse(grepl("multi", model), "Multitask", "Single-task"), improvement = 100 * (1 - mse / mse_m1),
          virus = factor(vlab[virus], vlab))]
print(ov[, .(model, improvement = round(improvement, 1))]); print(pv[, .(virus, task, improvement = round(improvement, 1))])
bars <- function(d, x) ggplot(d, aes(.data[[x]], improvement, fill = task)) +
  geom_col(position = position_dodge(.6), width = .5) +
  geom_text(aes(label = sprintf("%.1f%%", improvement)), position = position_dodge(.6), vjust = -.3, size = 2.3) +
  scale_fill_manual(values = c(Multitask = "#1B4F72", `Single-task` = "#85C1E9")) +
  scale_y_continuous(expand = expansion(mult = c(0, .15))) +
  labs(x = NULL, y = "Reduction in prediction error vs M1 (%)") + th
savefig((bars(ov, "model_base") | bars(pv, "virus")) + plot_layout(guides = "collect") +
          plot_annotation(tag_levels = "A", tag_suffix = ".") & theme(legend.position = "bottom"),
        "E2_multitask_vs_single.png", 3.2)

# ---- 2,000 versus 5,000 micro-clusters ----------------------------------------------
ci <- rbind(fread(file.path(root, "post_hoc_5000/recipient_bootstrap/E3_recipient_boot.csv"))[comparison == "M1_to_M3"][, n_micro := "5,000"],
            fread(file.path(root, "post_hoc_2000/recipient_bootstrap/E3_recipient_boot.csv"))[comparison == "M1_to_M3"][, n_micro := "2,000"])
ci[, scope := factor(scope, c("h1n1", "h3n2", "bvic", "overall"), c(vlab, "Overall"))]
print(ci[scope == "Overall", .(n_micro, delta_pct = round(delta_pct, 1))])
cols <- c(`2,000` = "#E67E22", `5,000` = "#1B4F72")
p1 <- ggplot(ci, aes(scope, delta_pct, color = n_micro)) + geom_hline(yintercept = 0, color = "grey60") +
  geom_linerange(aes(ymin = ci_low_pct, ymax = ci_high_pct), position = position_dodge(.5)) +
  geom_point(position = position_dodge(.5), size = 2) + scale_color_manual(values = cols) +
  labs(x = NULL, y = expression("M1" %->% "M3 reduction in prediction error (%)")) + th
ks <- rbindlist(lapply(c(2000, 5000), function(n) { a <- new.env(); load(file.path(root, paste0("per_fold_atlas_", n), "aggregated/E2_aggregated.RData"), a)
  data.table(fold = factor(1:5), k = a$selected_k_per_fold, n_micro = format(n, big.mark = ",")) }))
p2 <- ggplot(ks, aes(fold, k, fill = n_micro)) + geom_col(position = "dodge", width = .7) +
  scale_fill_manual(values = cols, guide = "none") + labs(x = "Outer fold", y = "Selected meta-clusters") + th
savefig((p1 | p2) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A", tag_suffix = ".") &
          theme(legend.position = "bottom"), "E2_sensitivity_2000v5000.png", 3.2)
cat("Done.\n")
