#!/usr/bin/env Rscript
# E14-04: Score-group tests (Supplementary Tables of per-antigen and composite
# score groups) recomputed from the E13-01 null replicates.
#
# E13-01 saves every replicate's out-of-fold scores and null outcomes, so the
# E9 statistics (high-minus-low mean log fold-rise and fourfold-response
# difference, within-quartile median split, three antigens x five strata) and
# the E12-02 composite statistics (same, composite score and outcome) can be
# computed from them directly. Outputs use the E9 and E12-02 file formats so
# the table generators read them unchanged.
suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/first_episode_2026-10-14")
NULL_DIR <- file.path(ROOT, "score_null_bootstrap_scores")
V <- c("h1n1", "h3n2", "bvic"); L4 <- log(4) - 1e-12
obs <- readRDS(file.path(NULL_DIR, "observed_scores.rds"))
quartile <- obs$quartile
scopes <- c(paste0("Q", 1:4), "Q1-Q3")
keys <- CJ(virus = V, scope = scopes, endpoint = c("log_fold_ratio", "response_difference"), sorted = FALSE)
ckeys <- CJ(outcome = "composite", scope = scopes, endpoint = c("log_fold_ratio", "response_difference"), sorted = FALSE)
grp_stat <- function(y, score, scope, endpoint) {
  g <- ave(score, quartile, FUN = function(x) as.numeric(x > median(x)))
  ix <- if (scope == "Q1-Q3") which(quartile != "Q4") else which(quartile == scope)
  hi <- ix[g[ix] == 1]; lo <- ix[g[ix] == 0]
  if (!length(hi) || !length(lo)) return(0)
  yy <- if (endpoint == "log_fold_ratio") y else as.numeric(y >= L4)
  mean(yy[hi]) - mean(yy[lo])
}
stats <- function(y, score, cscore) {
  a <- vapply(seq_len(nrow(keys)), function(j) { v <- match(keys$virus[j], V)
    grp_stat(y[, v], score[, v], keys$scope[j], keys$endpoint[j]) }, numeric(1))
  ybar <- rowMeans(y)
  b <- vapply(seq_len(nrow(ckeys)), function(j) grp_stat(ybar, cscore, ckeys$scope[j], ckeys$endpoint[j]), numeric(1))
  list(a = a, b = b)
}
o <- stats(obs$y, obs$score, obs$cscore)
files <- sort(list.files(NULL_DIR, pattern = "^replicate_\\d+\\.rds$", full.names = TRUE))
res <- parallel::mclapply(files, function(f) { r <- readRDS(f); stats(r$y, r$score, r$cscore) }, mc.cores = 10)
B <- length(res)
na <- do.call(rbind, lapply(res, `[[`, "a")); nb <- do.call(rbind, lapply(res, `[[`, "b"))
fin <- function(k, obsv, nul) {
  k[, estimate := obsv]
  k[, p_null := (1 + colSums(abs(nul) >= rep(abs(obsv), each = B))) / (B + 1)]
  k[, p_holm := p.adjust(p_null, "holm")]
  k[, monte_carlo_se := sqrt(p_null * (1 - p_null) / (B + 1))]
  k[, n_null := B]; k }
dir.create(file.path(ROOT, "score_null_bootstrap"), showWarnings = FALSE)
dir.create(file.path(ROOT, "score_null_bootstrap_composite"), showWarnings = FALSE)
fwrite(fin(keys, o$a, na), file.path(ROOT, "score_null_bootstrap/score_tests.csv"))
fwrite(fin(ckeys, o$b, nb), file.path(ROOT, "score_null_bootstrap_composite/composite_score_tests.csv"))
saveRDS(list(ids = obs$ids, cscore = obs$cscore, ybar = rowMeans(obs$y), quartile = quartile, k = obs$k),
        file.path(ROOT, "score_null_bootstrap_composite/observed_composite_scores.rds"))
print(keys); print(ckeys)
