#!/usr/bin/env Rscript
# Training-aware p-values for the clinical score statistics of E13-02.
#
# Each E13-01 null replicate carries out-of-fold scores and predictions from a
# full refit of the pipeline (resolution selection, penalty, coefficients) on
# outcomes generated under the null of no cellular effect beyond baseline
# titers. Every statistic below is recomputed on each replicate exactly as on
# the observed data, and p = (1 + #{|null| >= |observed|}) / (B + 1).
# Scores are standardized within fold, as in E13-02.
FLU_ROOT <- Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06")
suppressPackageStartupMessages({library(data.table); library(pROC)})
ROOT <- FLU_ROOT
NULL_DIR <- file.path(ROOT, "score_null_bootstrap_scores")
CS <- file.path(ROOT, "post_hoc_5000/clinical_score")
d <- readRDS(file.path(CS, "episode_table.rds"))
obs <- readRDS(file.path(NULL_DIR, "observed_scores.rds"))
ord <- match(d$episode_id, obs$ids); stopifnot(!anyNA(ord))
stopifnot(max(abs(obs$cscore[ord] - d$score_comp)) < 1e-9)
fold <- d$fold; base <- d$base_comp; cohort <- d$cohort; first <- d$first_episode
L4 <- log(4) - 1e-12
# A null fit can select no population in a fold (all scores zero); its
# standardized score is then zero (no information), not undefined.
zf <- function(x) ave(x, fold, FUN = function(v) if (sd(v) > 0) (v - mean(v)) / sd(v) else 0 * v)
safe_cor <- function(x, y) if (sd(x) > 0) cor(x, y, method = "spearman") else 0
prank <- function(x) (rank(x) - .5) / length(x)
pspear <- function(x, y, z) if (sd(x) == 0) 0 else cor(resid(lm(prank(x) ~ prank(z))), resid(lm(prank(y) ~ prank(z))))
lor <- function(resp, z, ix) { b <- coef(suppressWarnings(glm(resp[ix] ~ base[ix] + z[ix], family = binomial)))[3]; if (is.na(b)) 0 else unname(b) }
auc <- function(resp, p, ix) as.numeric(pROC::auc(pROC::roc(resp[ix], p[ix], quiet = TRUE, direction = "<")))
sets <- list(All = rep(TRUE, nrow(d)), Children = cohort == "Children", Adults = cohort == "Adults", `First episodes` = first)
stats <- function(y, cscore, cp1, cp3, score) {
  ybar <- rowMeans(y); z <- zf(cscore); resp <- ybar >= L4
  bottom <- rank(z, ties.method = "first") <= round(.2 * length(z))
  out <- c()
  for (s in names(sets)) { ix <- sets[[s]]
    out[paste0("spearman|", s)] <- safe_cor(z[ix], ybar[ix])
    out[paste0("logistic_or_per_sd|", s)] <- lor(resp, z, ix)
    out[paste0("delta_auc_M1_M3|", s)] <- auc(resp, cp3, ix) - auc(resp, cp1, ix) }
  out["partial_spearman|All"] <- pspear(z, ybar, base)
  out["bottom_fifth_difference|All"] <- mean(resp[!bottom]) - mean(resp[bottom])
  for (v in 1:3) { zv <- zf(score[, v])
    out[paste0("spearman_", c("h1n1","h3n2","bvic")[v], "|All")] <- safe_cor(zv, y[, v]) }
  out
}
o <- stats(obs$y[ord, ], obs$cscore[ord], obs$cp1[ord], obs$cp3[ord], obs$score[ord, ])
files <- sort(list.files(NULL_DIR, pattern = "^replicate_\\d+\\.rds$", full.names = TRUE))
cat("Null replicates found:", length(files), "\n")
nul <- parallel::mclapply(files, function(f) { r <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(r)) return(NULL)
  stats(r$y[ord, ], r$cscore[ord], r$cp1[ord], r$cp3[ord], r$score[ord, ]) }, mc.cores = 12)
nul <- do.call(rbind, Filter(function(x) is.numeric(x), nul))
B <- nrow(nul)
res <- data.table(stat_key = names(o), observed = o,
                  p_null = (1 + colSums(abs(nul) >= rep(abs(o), each = B))) / (B + 1),
                  null_mean = colMeans(nul), null_q975 = apply(nul, 2, quantile, .975), n_null = B)
res[, c("statistic", "subset") := tstrsplit(stat_key, "|", fixed = TRUE)]
res[, stat_key := NULL]
fwrite(res, file.path(CS, "null_pvalues.csv")); print(res)
