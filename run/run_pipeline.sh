#!/bin/bash
# Full analysis in run order. Each step writes to the folders set in run/env.sh
# and logs to logs/<step>.log. Steps 1-3 are the expensive ones (per-fold atlases
# take about 40 minutes per fold); later steps reuse their outputs.
cd "$(dirname "$0")/.."
source run/env.sh
mkdir -p logs paper/images paper/tables
run() { local name=$1; shift; echo "[$(date +%H:%M:%S)] START $name"; env "$@" > logs/$name.log 2>&1; local rc=$?; echo "[$(date +%H:%M:%S)] END $name exit=$rc"; }

# 1. HAI titers from the trial REDCap exports
run titers Rscript experiments/E17-07-redcap-titers.R

# 2. Per-fold atlases: in each outer training fold, 5,000 k-means micro-clusters
#    and Leiden partitions over the resolution grid (no outcomes used)
for f in 1 2 3 4 5; do
  run atlas-fold$f FLU_RESULTS_ROOT=$FLU_ATLAS_ROOT FOLD_ID=$f N_CORES=3 Rscript experiments/E14-01-first-episode-perfold.R
done

# 3. Resolution selection, model fits and held-out predictions in each fold
run refit Rscript experiments/E17-05-refit-missing-titers.R
run aggregate AGGREGATE_ONLY=1 Rscript experiments/E14-01-first-episode-perfold.R

# 4. Shared interpretation atlas (geometry, outcome objects, biology, labels)
run atlas-stage1 STAGE=1 Rscript experiments/E14-02-shared-atlas-first-episode.R
run atlas-stage2 STAGE=2 TARGET_K=35 FORCE=1 N_CORES=8 Rscript experiments/E14-02-shared-atlas-first-episode.R
run biology Rscript experiments/E14-02b-shared-atlas-biology.R
run labels Rscript experiments/E14-03b-annotation-transfer.R

# 5. Training-aware null bootstrap for the score tests (999 replicates)
run null CRITERION=composite N_CORES=36 Rscript experiments/E13-01-null-scores.R

# 6. Selection stability on the shared atlas
run bootmulti N_CORES=6 FORCE=1 Rscript experiments/E14-02c-bootstrap-importance.R
run bootcomp OUTCOME=composite N_CORES=6 FORCE=1 Rscript experiments/E14-02c-bootstrap-importance.R

# 7. Out-of-fold analyses, second-season validation and per-cohort analyses
run posthoc RUN_FULL_POSTHOC=TRUE Rscript experiments/E2-posthoc.R
run e3boot N_CORES=6 Rscript experiments/E3-recipient-bootstrap.R
run e12-01 Rscript experiments/E12-01-composite-endpoint.R
run e13-02 Rscript experiments/E13-02-clinical-score.R
run e14-05 Rscript experiments/E14-05-repeater-validation.R
run e14-06 Rscript experiments/E14-06-same-score-comparison.R
run e5 Rscript experiments/E5-per-cohort-refit.R
run e12-05 Rscript experiments/E12-05-predictor-signs.R
run e12-04 Rscript experiments/E12-04-shared-atlas-signature.R
run e13-06 Rscript experiments/E13-06-top-k-populations.R
run e6-01 Rscript experiments/E6-01-cellular-score-clinical.R
run e6-02 Rscript experiments/E6-02-predictive-cluster-recovery.R
run e6-04 Rscript experiments/E6-04-cluster-size-report.R
run e3int N_CORES=6 Rscript experiments/E3-interaction-baseline.R
run e3intM4 Rscript experiments/E3-interaction-M4.R
run e16-02 Rscript experiments/E16-02-composite-analyses.R
run e15-02 Rscript experiments/E15-02-selection-rules.R
run e14-04 Rscript experiments/E14-04-score-group-tests.R
run e13-04 Rscript experiments/E13-04-null-pvalues.R

# 8. Marker reduction
run e17-03 N_CORES=40 Rscript experiments/E17-03-marker-reduction-composite.R
run marker OUTCOME=composite Rscript experiments/E2-marker-reduction-perfold.R
run e12-03 Rscript experiments/E12-03-score-15-vs-30-markers.R
run e13-05 N_CORES=8 Rscript experiments/E13-05-marker-topdown.R

# 9. All-episode sensitivity analysis (reuses the cached per-fold geometry)
for n in 5000 2000; do
  run allep-$n FLU_RESULTS_ROOT=$FLU_ALLEP_ROOT N_MICRO=$n N_CORES=8 Rscript experiments/E2-perfold-4model.R
  run allep-e3-$n FLU_RESULTS_ROOT=$FLU_ALLEP_ROOT N_MICRO=$n N_CORES=8 Rscript experiments/E3-recipient-bootstrap.R
done
run allep-posthoc FLU_RESULTS_ROOT=$FLU_ALLEP_ROOT RUN_FULL_POSTHOC=TRUE Rscript experiments/E2-posthoc.R
run allep-e12-01 FLU_RESULTS_ROOT=$FLU_ALLEP_ROOT Rscript experiments/E12-01-composite-endpoint.R
run allep-e13-02 FLU_RESULTS_ROOT=$FLU_ALLEP_ROOT Rscript experiments/E13-02-clinical-score.R
run allep-single FLU_RESULTS_ROOT=$FLU_ALLEP_ROOT Rscript experiments/E2-single-task-comparison.R

# 10. Figures and tables
run e14-09 Rscript experiments/E14-09-annotation-table.R
run e15-03 Rscript experiments/E15-03-gain-by-k-figure.R
run e16-04 Rscript experiments/E16-04-composite-figures.R
run e16-06 python3 experiments/E16-06-multitask-table.py
run e10 Rscript experiments/E10-main-figure-presentation.R
run e11 Rscript experiments/E11-supplement-presentation.R
run e7 Rscript experiments/E7-heatmap-catalyst-scaling.R
run e8u Rscript experiments/E8-refresh-umap-labels.R
run e14-07 Rscript experiments/E14-07-repeater-figure.R
run e6-02c FLU_BOOT_FILE=$FLU_BOOT_FILE_COMP Rscript experiments/E6-02-predictive-cluster-recovery.R
run e6-04c FLU_BOOT_FILE=$FLU_BOOT_FILE_COMP Rscript experiments/E6-04-cluster-size-report.R
run e14-08 FLU_BOOT_FILE=$FLU_BOOT_FILE_COMP Rscript experiments/E14-08-supplement-figures.R
run e14-10 Rscript experiments/E14-10-supplement-figures-2.R
run e17-02 Rscript experiments/E17-02-allepisode-figures.R
run e13-07 python3 experiments/E13-07-clinical-table.py
run e12-07 python3 experiments/E12-07-supplement-tables.py
run e16-05 python3 experiments/E16-05-cohort-signs-table.py
run e9rep python3 experiments/E9-score-test-report.py
run e13-03 Rscript experiments/E13-03-figures.R
run e14-11 Rscript experiments/E14-11-numbers.R
echo "[$(date +%H:%M:%S)] DONE"
