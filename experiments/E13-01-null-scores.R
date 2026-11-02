#!/usr/bin/env Rscript
# Training-aware null replicates that keep the fitted scores and predictions.
#
# Same null model, seeds, folds, resolution selection, penalties and score
# construction as E9-score-null-bootstrap.R and E12-02-composite-null-bootstrap.R.
# The difference: every replicate saves its out-of-fold cellular scores (three
# antigens and composite), its out-of-fold M1 (titers only) and M3 (titers plus
# cellular composition) predictions, and the null outcomes. Any statistic of the
# scores (rank correlations, response rates by score level, logistic odds
# ratios, areas under ROC curves, cohort interactions) can then be given a
# p-value that accounts for the whole model-fitting pipeline, without refitting.
# Each replicate also reproduces the 30 per-antigen E9 statistics, checked
# against the saved E9 replicate.
FLU_ROOT <- Sys.getenv("FLU_RESULTS_ROOT", "/fast/datasets/flu_cellular_prediction/timing_corrected_2026-09-06")
suppressPackageStartupMessages({library(data.table); library(glmnet); library(here)})
source(here("methods", "35-per-fold-atlas-cv.R"))
ROOT <- FLU_ROOT
E9_OUT <- file.path(ROOT, "score_null_bootstrap")
OUT <- file.path(ROOT, "score_null_bootstrap_scores")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
B <- as.integer(Sys.getenv("B_NULL", "999"))
N_CORES <- as.integer(Sys.getenv("N_CORES", "30"))
V <- c("h1n1", "h3n2", "bvic")
# CRITERION=composite reselects the resolution in each replicate with the
# composite-outcome inner criterion (single-outcome M4 elastic net on the mean
# log fold-rise versus OLS on titers and clinical factors), matching E16-01.
CRITERION <- Sys.getenv("CRITERION", "multitask")
a <- new.env(); load(file.path(ROOT, "per_fold_atlas_5000/aggregated/E2_aggregated.RData"), a)
ids <- a$common_p_ids; Y <- a$oof_Y
gid <- setNames(base_pid_from_p_id_yr(ids), ids)
cohort <- infer_cohort_from_p_id(ids)

# --- Fold structures (identical to E9) --------------------------------------
folds <- vector("list", 5L)
for (f in 1:5) {
  e <- new.env(); fd <- file.path(ROOT, "per_fold_atlas_5000", paste0("fold_", f))
  load(file.path(fd, sprintf("fold_%d_all.RData", f)), e)
  tr <- e$train_p_ids; te <- e$test_p_ids
  maps <- readRDS(file.path(fd, "resolution_candidates.rds"))
  inner <- create_group_folds(tr, gid, 3L, 7000 + f * 100, cohort[tr])
  splits <- lapply(seq_along(inner), function(ii) {
    va <- tr[inner[[ii]]]; tt <- setdiff(tr, va)
    Xtr <- cbind(1, e$X_M2_train[tt, , drop=FALSE])
    Xva <- cbind(1, e$X_M2_train[va, , drop=FALSE])
    list(tr=match(tt, ids), va=match(va, ids), ti=match(tt, tr), vi=match(va, tr),
         baseline_map=Xva %*% qr.solve(Xtr, diag(nrow(Xtr))),
         fid=compute_grouped_foldid(tt, 5L, 7100 + f*10000 + ii))
  })
  candidates <- lapply(maps, function(m) {
    ct <- transform_counts(two_stage_collapse_micro_to_meta_counts(e$train_micro_counts, m$micro_to_meta), "clr")
    ce <- transform_counts(two_stage_collapse_micro_to_meta_counts(e$test_micro_counts, m$micro_to_meta), "clr")
    list(x4=cbind(ct[tr,,drop=FALSE], e$X_M2_train),
         x3=cbind(ct[tr,,drop=FALSE], e$X_M1_train),
         ctest=ce[te,,drop=FALSE])
  })
  folds[[f]] <- list(tr=match(tr,ids), te=match(te,ids), splits=splits, candidates=candidates,
    fid=compute_grouped_foldid(tr, 5L, 7100 + f*1000 + 99),
    m1_train=e$X_M1_train[tr,,drop=FALSE], m1_test=e$X_M1_test[te,,drop=FALSE],
    observed_k=e$k_found, observed_score=sapply(V,function(v)
      as.vector(e$clr_test %*% e$coef_M3_en[[v]][colnames(e$clr_test)])))
}
base <- rbind(e$X_M1_train,e$X_M1_test)[ids,,drop=FALSE]
avg_base <- rowMeans(base)
quartile <- cut(avg_base, quantile(avg_base, 0:4/4), include.lowest=TRUE, labels=paste0("Q",1:4))
scopes <- c(paste0("Q",1:4),"Q1-Q3")
keys <- CJ(virus=V, scope=scopes, endpoint=c("log_fold_ratio","response_difference"), sorted=FALSE)
ckeys <- CJ(outcome="composite", scope=scopes, endpoint=c("log_fold_ratio","response_difference"), sorted=FALSE)

# Per-antigen statistics, exactly as in E9.
statistics <- function(y, score) {
  ans <- numeric(nrow(keys))
  groups <- sapply(1:3, function(v) ave(score[,v], quartile, FUN=function(x) as.numeric(x>median(x))))
  for (j in seq_len(nrow(keys))) {
    v <- match(keys$virus[j], V)
    ix <- if(keys$scope[j]=="Q1-Q3") which(quartile!="Q4") else which(quartile==keys$scope[j])
    hi <- ix[groups[ix,v]==1]; lo <- ix[groups[ix,v]==0]
    if (!length(hi) || !length(lo)) {ans[j]<-0; next}
    yy <- if(keys$endpoint[j]=="log_fold_ratio") y[,v] else as.numeric(y[,v]>=log(4)-1e-12)
    ans[j] <- mean(yy[hi])-mean(yy[lo])
  }
  ans
}

# Composite statistics: same within-quartile median split, applied to the
# composite score and the composite outcome.
composite_statistics <- function(ybar, cscore) {
  ans <- numeric(nrow(ckeys))
  groups <- ave(cscore, quartile, FUN=function(x) as.numeric(x>median(x)))
  for (j in seq_len(nrow(ckeys))) {
    ix <- if(ckeys$scope[j]=="Q1-Q3") which(quartile!="Q4") else which(quartile==ckeys$scope[j])
    hi <- ix[groups[ix]==1]; lo <- ix[groups[ix]==0]
    if (!length(hi) || !length(lo)) {ans[j]<-0; next}
    yy <- if(ckeys$endpoint[j]=="log_fold_ratio") ybar else as.numeric(ybar>=log(4)-1e-12)
    ans[j] <- mean(yy[hi])-mean(yy[lo])
  }
  ans
}

fit_scores <- function(y) {
  score <- matrix(NA_real_,nrow(Y),3,dimnames=dimnames(Y))
  p1 <- p3 <- score
  cscore <- cp1 <- cp3 <- setNames(rep(NA_real_, nrow(Y)), rownames(Y))
  ybar <- rowMeans(y)
  selected <- integer(5)
  for (f in 1:5) {
    ff <- folds[[f]]
    bmse <- lapply(ff$splits, function(s) colMeans((s$baseline_map %*% y[s$tr,,drop=FALSE]-y[s$va,,drop=FALSE])^2))
    if (CRITERION == "composite") {
      cbmse <- lapply(ff$splits, function(s) mean((s$baseline_map %*% ybar[s$tr] - ybar[s$va])^2))
      losses <- vapply(ff$candidates, function(ca) {
        mean(vapply(seq_along(ff$splits), function(ii) {
          s <- ff$splits[[ii]]
          fit <- cv.glmnet(ca$x4[s$ti,,drop=FALSE], ybar[s$tr], family="gaussian", alpha=.9,
            foldid=s$fid, standardize=TRUE, parallel=FALSE)
          mean((as.vector(predict(fit, ca$x4[s$vi,,drop=FALSE], s="lambda.min")) - ybar[s$va])^2) / cbmse[[ii]]
        }, numeric(1)))
      }, numeric(1))
    } else {
    losses <- vapply(ff$candidates, function(ca) {
      mean(vapply(seq_along(ff$splits),function(ii) {
        s<-ff$splits[[ii]]
        fit<-cv.glmnet(ca$x4[s$ti,,drop=FALSE], y[s$tr,,drop=FALSE],
          family="mgaussian",alpha=.9,foldid=s$fid,standardize=TRUE,parallel=FALSE)
        pred<-as_pred_matrix(predict(fit,ca$x4[s$vi,,drop=FALSE],s="lambda.min"),V)
        median(colMeans((pred-y[s$va,,drop=FALSE])^2)/bmse[[ii]])
      },numeric(1)))
    },numeric(1))
    }
    ca <- ff$candidates[[which.min(losses)]]
    selected[f] <- ncol(ca$ctest)
    x3te <- cbind(ca$ctest, ff$m1_test)
    # Multitask M3 (per-antigen scores, as in E9)
    fit <- cv.glmnet(ca$x3, y[ff$tr,,drop=FALSE],family="mgaussian",alpha=.9,
      foldid=ff$fid,standardize=TRUE,parallel=FALSE)
    cf <- coef(fit,s="lambda.min")
    for(v in 1:3) score[ff$te,v]<-ca$ctest %*% as.matrix(cf[[v]])[colnames(ca$ctest),1]
    p3[ff$te,] <- as_pred_matrix(predict(fit, x3te, s="lambda.min"), V)
    # Multitask M1 (titers only)
    fit1 <- cv.glmnet(ff$m1_train, y[ff$tr,,drop=FALSE],family="mgaussian",alpha=.9,
      foldid=ff$fid,standardize=TRUE,parallel=FALSE)
    p1[ff$te,] <- as_pred_matrix(predict(fit1, ff$m1_test, s="lambda.min"), V)
    # Composite M3 and M1 (single outcome)
    cfit <- cv.glmnet(ca$x3, ybar[ff$tr], family="gaussian", alpha=.9,
      foldid=ff$fid, standardize=TRUE, parallel=FALSE)
    ccf <- as.matrix(coef(cfit, s="lambda.min"))
    cscore[ff$te] <- ca$ctest %*% ccf[colnames(ca$ctest),1]
    cp3[ff$te] <- as.vector(predict(cfit, x3te, s="lambda.min"))
    cfit1 <- cv.glmnet(ff$m1_train, ybar[ff$tr], family="gaussian", alpha=.9,
      foldid=ff$fid, standardize=TRUE, parallel=FALSE)
    cp1[ff$te] <- as.vector(predict(cfit1, ff$m1_test, s="lambda.min"))
  }
  list(score=score, cscore=cscore, p1=p1, p3=p3, cp1=cp1, cp3=cp3, k=selected)
}

cat("Verifying exact observed pipeline reproduction...\n")
tic <- proc.time()[3]
obs <- fit_scores(Y)
saved_score <- matrix(NA_real_, nrow(Y),3)
for(f in 1:5) saved_score[folds[[f]]$te,] <- folds[[f]]$observed_score
stopifnot(max(abs(obs$score-saved_score))<1e-9,
          identical(obs$k, as.integer(a$selected_k_per_fold)))
cat("Verified. Seconds per serial replicate:",proc.time()[3]-tic,"\n")
ybar_obs <- rowMeans(Y)
observed_c <- composite_statistics(ybar_obs, obs$cscore)
ref_c <- file.path(ROOT, "score_null_bootstrap_composite/observed_composite_scores.rds")
if (file.exists(ref_c)) stopifnot(max(abs(obs$cscore - readRDS(ref_c)$cscore)) < 1e-9)
saveRDS(c(list(ids=ids, y=Y, quartile=quartile, cohort=cohort, gid=gid), obs),
        file.path(OUT, "observed_scores.rds"))

Xnull <- cbind(1,base)
nullfit <- lm.fit(Xnull,Y)
mu <- nullfit$fitted.values; residual <- nullfit$residuals
leverage <- rowSums(qr.Q(qr(Xnull))^2)
residual <- residual/sqrt(1-leverage)
stopifnot(all(abs(mu+residual-log(4))>1e-12), all(abs(mu-residual-log(4))>1e-12))
unique_ids <- unique(gid)
one <- function(b) {
  path <- file.path(OUT,sprintf("replicate_%04d.rds",b))
  if(file.exists(path)) return(readRDS(path))
  set.seed(20260906+b)
  w <- sample(c(-1,1),length(unique_ids),replace=TRUE)[match(gid,unique_ids)]
  yb <- mu+residual*w
  fit <- fit_scores(yb)
  stat <- statistics(yb, fit$score)
  # Exact agreement with the saved E9 replicate confirms identical null data.
  e9 <- file.path(E9_OUT, sprintf("replicate_%04d.rds", b))
  e9_diff <- if (file.exists(e9)) max(abs(readRDS(e9)$stat - stat)) else NA_real_
  z <- c(list(y=yb, e9_max_abs_diff=e9_diff), fit)
  saveRDS(z,path)
  z
}
cat("Running",B,"null replicates (scores saved) with",N_CORES,"workers...\n");flush.console()
res <- parallel::mclapply(seq_len(B),one,mc.cores=N_CORES,mc.preschedule=FALSE,mc.set.seed=FALSE)
stopifnot(!any(vapply(res,inherits,logical(1),"try-error")))
e9_check <- vapply(res, `[[`, numeric(1), "e9_max_abs_diff")
cat("Max abs difference from saved per-antigen null statistics:", max(e9_check, na.rm=TRUE),
    " (replicates compared:", sum(!is.na(e9_check)), ")\n")
stopifnot(max(e9_check, na.rm=TRUE) < 1e-9)
cat("Done:", B, "replicates saved in", OUT, "\n")
