#!/usr/bin/env Rscript
# E14-06: Responders and non-responders at the same cellular score.
#
# Question raised by Jim Chappell (October 13 meeting): episodes with similar
# cellular composition scores can differ sharply in response, for example a
# low-baseline responder and a high-baseline non-responder at the same score.
# What else differs between them?
#
# Within each fifth of the composite score (first episodes), responders
# (composite fourfold rise) are compared with non-responders on baseline titer,
# age, months since transplantation, high-dose vaccine, cohort, lineage shares
# and population shares (shared first-episode atlas). Differences are
# standardized mean differences (responders minus non-responders, pooled SD).
# Population shares use log(share + 0.01%). A pooled estimate across fifths
# is the average of the within-fifth differences weighted by fifth size.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/first_episode_2026-10-14")
CS <- file.path(ROOT, "post_hoc_5000/clinical_score")
OUT <- file.path(ROOT, "post_hoc_5000/same_score"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
d <- readRDS(file.path(CS, "episode_table.rds"))
bio <- new.env(); load(Sys.getenv("FLU_BIO_FILE"), bio)
P <- bio$meta_props_all; rownames(P) <- bio$common_p_ids
ann <- fread(Sys.getenv("FLU_ANN_FILE"))
lin <- ann$broad_lineage[match(seq_len(ncol(P)), ann$cluster_id)]
L <- sapply(unique(lin), function(l) rowSums(P[, which(lin == l), drop = FALSE]))
feat <- data.table(episode_id = d$episode_id,
  `Baseline titer (log)` = d$base_comp, `Age` = d$age, `Months since transplant (log)` = log(d$months),
  `High-dose vaccine` = as.numeric(d$dose == "HD"), `Adult cohort` = as.numeric(d$cohort == "Adults"))
for (l in colnames(L)) feat[, (paste0("Lineage: ", l)) := log(L[d$episode_id, l] + 1e-4)]
for (j in seq_len(ncol(P))) feat[, (sprintf("C%d: %s", j, ann$annotation[ann$cluster_id == j])) := log(P[d$episode_id, j] + 1e-4)]
fcols <- setdiff(names(feat), "episode_id")
smd <- function(x, g) { a <- x[g]; b <- x[!g]; s <- sqrt(((length(a) - 1) * var(a) + (length(b) - 1) * var(b)) / (length(x) - 2))
  if (!is.finite(s) || s == 0) return(NA_real_); (mean(a) - mean(b)) / s }
res <- rbindlist(lapply(levels(d$fifth), function(q) { ix <- which(d$fifth == q); g <- d$resp_comp[ix]
  data.table(fifth = q, n_resp = sum(g), n_nonresp = sum(!g), feature = fcols,
             smd = vapply(fcols, function(f) smd(feat[[f]][ix], g), numeric(1))) }))
pooled <- res[, .(smd = weighted.mean(smd, n_resp + n_nonresp, na.rm = TRUE)), by = feature][order(-abs(smd))]
fwrite(res, file.path(OUT, "smd_by_fifth.csv")); fwrite(pooled, file.path(OUT, "smd_pooled.csv"))
print(head(pooled, 25))
# Same comparison after accounting for baseline titer: residualize each
# feature on log baseline titer within the fifth.
res_adj <- rbindlist(lapply(levels(d$fifth), function(q) { ix <- which(d$fifth == q); g <- d$resp_comp[ix]; b <- d$base_comp[ix]
  data.table(fifth = q, feature = fcols, smd = vapply(fcols, function(f) {
    x <- feat[[f]][ix]; if (f == "Baseline titer (log)" || sd(x) == 0) return(NA_real_); smd(resid(lm(x ~ b)), g) }, numeric(1))) }))
pooled_adj <- res_adj[, .(smd_baseline_adjusted = mean(smd, na.rm = TRUE)), by = feature]
fwrite(merge(pooled, pooled_adj, by = "feature")[order(-abs(smd))], file.path(OUT, "smd_pooled_with_adjusted.csv"))
# Jim's pair: closest scores, one responder at low baseline, one non-responder at high baseline
cand <- CJ(i = which(d$resp_comp & d$base_comp <= quantile(d$base_comp, .25)),
           j = which(!d$resp_comp & d$base_comp >= quantile(d$base_comp, .75)))
cand[, dz := abs(d$z_comp[i] - d$z_comp[j])]
pair <- cand[order(dz)][1:5]
pair[, `:=`(resp_ep = d$episode_id[i], nonresp_ep = d$episode_id[j], z_resp = d$z_comp[i], z_non = d$z_comp[j],
            base_resp = exp(d$base_comp[i]), base_non = exp(d$base_comp[j]))]
fwrite(pair, file.path(OUT, "matched_pairs.csv")); print(pair)
cat("Done.\n")
