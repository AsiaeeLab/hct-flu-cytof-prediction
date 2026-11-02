#!/usr/bin/env Rscript
# E14-03b: Carry population labels from the all-episode shared atlas (33
# populations, annotated by Josh Simmons and Spyros Kalams) to the
# first-episode shared atlas.
#
# Each new population is compared with each old population by its 30-marker
# mean expression profile (original marker scale, each marker standardized
# across the union of old and new populations). The best match by Pearson
# correlation donates its label and lineage when the correlation is at least
# 0.90 and the match is unambiguous (margin over the second-best match of at
# least 0.05); otherwise the label is marked for review. One old label may be
# carried by more than one new population (a split) and some old labels may
# have no match (a merge); both are listed. Josh rechecks every label.
suppressPackageStartupMessages(library(data.table))
OLD_BIO <- "/home/amir/repo/flu-cellular-prediction/results/final_atlas_production/biology_production_90203bf3.RData"
OLD_ANN <- "/home/amir/repo/flu-cellular-prediction/results/reports/cluster_annotations_corrected.csv"
NEW_BIO <- Sys.getenv("FLU_BIO_FILE")
OUT <- Sys.getenv("FLU_ANN_FILE")
stopifnot(file.exists(NEW_BIO), nzchar(OUT))
prof <- function(file) { e <- new.env(); load(file, e)
  m <- as.data.frame(e$means_original); rownames(m) <- m$meta; as.matrix(m[order(m$meta), e$marker_names]) }
old <- prof(OLD_BIO); new <- prof(NEW_BIO)
stopifnot(identical(colnames(old), colnames(new)))
both <- rbind(old, new); mu <- colMeans(both); sdv <- apply(both, 2, sd)
zo <- sweep(sweep(old, 2, mu), 2, sdv, "/"); zn <- sweep(sweep(new, 2, mu), 2, sdv, "/")
R <- cor(t(zn), t(zo))
ann <- fread(OLD_ANN)
res <- rbindlist(lapply(seq_len(nrow(R)), function(i) { r <- sort(R[i, ], decreasing = TRUE)
  best <- as.integer(names(r)[1]); a <- ann[cluster_id == best]
  ok <- r[1] >= 0.90 && (r[1] - r[2]) >= 0.05
  data.table(cluster_id = as.integer(rownames(R)[i]), matched_old_cluster = best, match_r = r[1],
             second_old_cluster = as.integer(names(r)[2]), second_r = r[2],
             annotation = a$annotation, broad_lineage = a$broad_lineage,
             confidence = if (ok) "carried" else "review", reasoning = a$reasoning) }))
res[, carried_by_n := .N, by = matched_old_cluster]
# Populations that share a carried label (splits of one old population) are
# distinguished by a subset number in the order of their IDs.
res <- res[order(cluster_id)]
res[, subset_no := seq_len(.N), by = annotation]
res[, n_same := .N, by = annotation]
res[n_same > 1 & subset_no > 1, annotation := sprintf("%s, subset %d", annotation, subset_no)]
res[, c("subset_no", "n_same") := NULL]
fwrite(res, OUT)
cat(sprintf("New populations: %d; carried: %d; to review: %d\n", nrow(res), sum(res$confidence == "carried"), sum(res$confidence == "review")))
print(res[, .(cluster_id, matched_old_cluster, match_r = round(match_r, 3), second_r = round(second_r, 3), confidence, carried_by_n, annotation)])
unmatched <- setdiff(ann$cluster_id, res$matched_old_cluster)
cat("Old populations with no new match:", if (length(unmatched)) paste(unmatched, collapse = ", ") else "none", "\n")
if (length(unmatched)) print(ann[cluster_id %in% unmatched, .(cluster_id, annotation)])
