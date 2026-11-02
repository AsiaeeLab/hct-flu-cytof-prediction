#!/usr/bin/env Rscript
# E14-05: Second-season episodes as a within-person validation set.
#
# The primary models (E14-01) were fitted on first-season episodes only. Each
# second-season episode was mapped to the atlas of the fold that held out its
# recipient and scored with that fold's models, so neither the episode nor its
# recipient's first episode was used in fitting the model that scores it. Both
# seasons of a recipient are therefore scored by the same model.
#
#   1. Validation in year 2: composite fourfold response against the
#      year-2 score at a given year-2 baseline titer (odds ratio per SD), AUC of
#      titers-only versus titers-plus-composition predictions, response by
#      score group with cut points taken from the first-episode distribution.
#   2. Within-person change: score, baseline titer and lineage composition in
#      year 1 and year 2 by response transition (non-responder or responder
#      in each season).
#   3. One dose: response four weeks after the first dose (V2/V1) in each
#      season, against the score.
# Scores are standardized with the mean and SD of the fold's first-episode
# out-of-fold scores, so both seasons share the first-episode scale.
suppressPackageStartupMessages({library(data.table); library(sandwich); library(lmtest); library(pROC)})
ROOT <- Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/first_episode_2026-10-14")
FOLD_ROOT <- file.path(ROOT, "per_fold_atlas_5000")
OUT <- file.path(ROOT, "post_hoc_5000/repeater_validation"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
BIO_FILE <- Sys.getenv("FLU_BIO_FILE")
ANN_FILE <- Sys.getenv("FLU_ANN_FILE")
V <- c("h1n1", "h3n2", "bvic"); L4 <- log(4) - 1e-12; B <- 2000L

cl_score <- function(clr, coef) { cl <- intersect(colnames(clr), rownames(coef)); drop(clr[, cl, drop = FALSE] %*% coef[cl, 1]) }
rows <- list()
for (f in 1:5) {
  e <- new.env(); load(file.path(FOLD_ROOT, sprintf("fold_%d", f), sprintf("fold_%d_all.RData", f)), e)
  s_first <- cl_score(e$clr_test[e$test_p_ids, , drop = FALSE], e$comp_M3$coef)
  s_rep <- cl_score(e$clr_repeat[e$repeat_p_ids, , drop = FALSE], e$comp_M3$coef)
  mu <- mean(s_first); sdv <- sd(s_first)
  mk <- function(ids, s, X1, Y, p1, p3, season) data.table(episode_id = ids, fold = f, season = season,
    score_raw = s, z = (s - mu) / sdv, base_comp = rowMeans(X1), min_base = apply(X1, 1, min), y_comp = rowMeans(Y), pred_M1 = p1, pred_M3 = p3)
  rows[[length(rows) + 1]] <- mk(e$test_p_ids, s_first, e$X_M1_test, e$Y_test, e$comp_M1$test, e$comp_M3$test, "first")
  rows[[length(rows) + 1]] <- mk(e$repeat_p_ids, s_rep, e$X_M1_repeat, e$Y_repeat, e$comp_M1$repeat_, e$comp_M3$repeat_, "second")
}
d <- rbindlist(rows)
# Baseline (V1) titers recorded as 1 lie below the assay's reporting floor (values
# below 1:10 are set to 1:5), so they are missing values, not measurements. Three
# second-season adult episodes have V1 = 1 for every antigen; their fold-rises are
# undefined and they are excluded. No first-season episode is affected.
bad <- d[min_base < log(5) - 1e-9]
cat("Excluded for missing baseline titer:", paste(bad$episode_id, collapse = ", "), "\n")
stopifnot(all(bad$season == "second"))
d <- d[min_base >= log(5) - 1e-9]
d[, base_id := sub("_yr[0-9]+$", "", episode_id)]
d[, resp := y_comp >= L4]
meta <- fread("~/datasets/flu/raw/hsct_metadata_15JUL2025.csv"); meta[, episode_id := sub("_([0-9]+)$", "_yr\\1", id_yr)]
meta <- meta[match(d$episode_id, episode_id)]
d[, cohort := factor(ifelse(grepl("^PSCT", episode_id), "Children", "Adults"), c("Children", "Adults"))]
d[, months := meta$mo_post]
# One-dose response: V2/V1 composite (geometric mean over the three antigens),
# with titers from the REDCap exports (experiments/E17-07-redcap-titers.R)
hai <- fread("~/datasets/flu/raw/hct_hai_redcap_2026-10-29.csv"); hai[, episode_id := sub("_([0-9]+)$", "_yr\\1", id_yr)]
hai <- hai[match(d$episode_id, episode_id)]; stopifnot(!anyNA(hai$episode_id))
v2 <- sapply(V, function(v) log(hai[[paste0("geo_mean_", v, "_v2")]]) - log(hai[[paste0("geo_mean_", v, "_v1")]]))
d[, y_v2 := rowMeans(v2)]; d[, resp_v2 := y_v2 >= L4]
# Score groups with cut points from the first-episode distribution
cuts <- quantile(d[season == "first", z], 0:5/5); cuts[1] <- -Inf; cuts[6] <- Inf
d[, fifth := cut(z, cuts, labels = paste0("F", 1:5))]
saveRDS(d, file.path(OUT, "repeater_episode_table.rds"))
cat(sprintf("First-season episodes %d; second-season episodes %d (recipients with both: %d)\n",
            d[season == "first", .N], d[season == "second", .N],
            length(intersect(d[season == "first", base_id], d[season == "second", base_id]))))

# 1. Validation in year 2
y2 <- d[season == "second"]
rob <- function(fit, cl) { ct <- coeftest(fit, vcov = vcovCL(fit, cluster = cl)); ct }
lg <- glm(resp ~ base_comp + z, binomial, y2); ct <- rob(lg, y2$base_id)
r1 <- pROC::roc(y2$resp, y2$pred_M1, quiet = TRUE, direction = "<"); r3 <- pROC::roc(y2$resp, y2$pred_M3, quiet = TRUE, direction = "<")
val <- data.table(
  n = nrow(y2), responders = sum(y2$resp),
  or_per_sd = exp(ct["z", 1]), or_lo = exp(ct["z", 1] - 1.96 * ct["z", 2]), or_hi = exp(ct["z", 1] + 1.96 * ct["z", 2]), p_or = ct["z", 4],
  auc_M1 = as.numeric(r1$auc), auc_M3 = as.numeric(r3$auc),
  p_delong = pROC::roc.test(r1, r3, method = "delong", paired = TRUE)$p.value,
  spearman = cor(y2$z, y2$y_comp, method = "spearman"),
  p_spearman = cor.test(y2$z, y2$y_comp, method = "spearman", exact = FALSE)$p.value,
  mse_gain_pct = 100 * (mean((y2$pred_M1 - y2$y_comp)^2) - mean((y2$pred_M3 - y2$y_comp)^2)) / mean((y2$pred_M1 - y2$y_comp)^2))
# recipient bootstrap for the MSE gain and AUC difference in year 2
set.seed(20261014L)
bt <- replicate(B, { ix <- sample.int(nrow(y2), replace = TRUE); x <- y2[ix]
  c(100 * (mean((x$pred_M1 - x$y_comp)^2) - mean((x$pred_M3 - x$y_comp)^2)) / mean((x$pred_M1 - x$y_comp)^2),
    tryCatch(as.numeric(pROC::auc(x$resp, x$pred_M3, quiet = TRUE, direction = "<")) -
             as.numeric(pROC::auc(x$resp, x$pred_M1, quiet = TRUE, direction = "<")), error = function(e) NA)) })
val[, `:=`(mse_gain_lo = quantile(bt[1, ], .025), mse_gain_hi = quantile(bt[1, ], .975),
           p_mse_gain = min(1, 2 * min(mean(bt[1, ] <= 0), mean(bt[1, ] >= 0))),
           dauc_lo = quantile(bt[2, ], .025, na.rm = TRUE), dauc_hi = quantile(bt[2, ], .975, na.rm = TRUE))]
fwrite(val, file.path(OUT, "year2_validation.csv")); print(val)
fifths <- d[, .(n = .N, k = sum(resp), rate = mean(resp)), by = .(season, fifth)][order(season, fifth)]
fwrite(fifths, file.path(OUT, "response_by_fifth_by_season.csv")); print(fifths)

# 2. Within-person change
w <- dcast(d[base_id %in% d[season == "second", base_id] & base_id %in% d[season == "first", base_id]],
           base_id + cohort ~ season, value.var = c("z", "resp", "base_comp", "y_comp", "resp_v2", "y_v2", "months"))
w[, transition := factor(paste0(ifelse(resp_first, "R", "NR"), " to ", ifelse(resp_second, "R", "NR")),
                         c("NR to NR", "NR to R", "R to NR", "R to R"))]
w[, `:=`(dz = z_second - z_first, dbase = exp(base_comp_second - base_comp_first))]
tr <- w[, .(n = .N, z_first = median(z_first), z_second = median(z_second), median_dz = median(dz),
            base_first = median(exp(base_comp_first)), base_second = median(exp(base_comp_second)),
            median_base_fold_change = median(dbase)), by = transition][order(transition)]
fwrite(tr, file.path(OUT, "transitions.csv")); print(tr)
# Does the change in score predict becoming a responder, given the change in baseline?
nr1 <- w[resp_first == FALSE]
lg2 <- glm(resp_second ~ I(base_comp_second - base_comp_first) + dz, binomial, nr1)
print(summary(lg2)$coefficients)
# year-2 response by year-2 score among year-1 non-responders
saveRDS(w, file.path(OUT, "repeater_pairs.rds"))

# 3. One dose (V2)
od <- d[is.finite(y_v2), .(n = .N, one_dose_fourfold = mean(resp_v2), spearman_z_v2 = cor(z, y_v2, method = "spearman")), by = season]
print(od); fwrite(od, file.path(OUT, "one_dose.csv"))
lg3 <- glm(resp_v2 ~ base_comp + z, binomial, d[season == "second" & is.finite(y_v2)])
ct3 <- rob(lg3, d[season == "second" & is.finite(y_v2)]$base_id); print(ct3)
lg4 <- glm(resp_v2 ~ base_comp + z, binomial, d[season == "first" & is.finite(y_v2)])
ct4 <- rob(lg4, d[season == "first" & is.finite(y_v2)]$base_id); print(ct4)
fwrite(rbind(data.table(season = "second", or = exp(ct3["z", 1]), lo = exp(ct3["z", 1] - 1.96 * ct3["z", 2]), hi = exp(ct3["z", 1] + 1.96 * ct3["z", 2]), p = ct3["z", 4]),
             data.table(season = "first", or = exp(ct4["z", 1]), lo = exp(ct4["z", 1] - 1.96 * ct4["z", 2]), hi = exp(ct4["z", 1] + 1.96 * ct4["z", 2]), p = ct4["z", 4])),
       file.path(OUT, "one_dose_logistic.csv"))

# Composition change (shared first-episode atlas, lineage shares)
if (nzchar(BIO_FILE) && file.exists(BIO_FILE) && nzchar(ANN_FILE) && file.exists(ANN_FILE)) {
  bio <- new.env(); load(BIO_FILE, bio)
  P <- bio$meta_props_every; rownames(P) <- bio$all_p_ids
  ann <- fread(ANN_FILE)
  lin <- ann$broad_lineage[match(seq_len(ncol(P)), ann$cluster_id)]
  L <- sapply(unique(lin), function(l) rowSums(P[, which(lin == l), drop = FALSE])) * 100
  ep <- d[base_id %in% w$base_id]
  comp <- cbind(ep[, .(episode_id, base_id, season)], L[ep$episode_id, , drop = FALSE])
  long <- melt(comp, id.vars = c("episode_id", "base_id", "season"), variable.name = "lineage", value.name = "pct")
  cw <- dcast(long, base_id + lineage ~ season, value.var = "pct")
  cw <- merge(cw, w[, .(base_id, transition)], by = "base_id")
  lc <- cw[, .(n = .N, first = median(first), second = median(second), median_change = median(second - first)), by = .(transition, lineage)]
  fwrite(lc, file.path(OUT, "lineage_change.csv")); print(lc[order(lineage, transition)])
}
cat("Done.\n")
