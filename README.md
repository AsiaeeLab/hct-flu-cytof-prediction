# Deep immune cell phenotyping predicts influenza vaccine immunogenicity in HCT recipients

Analysis code for the manuscript of the same title (Asiaee, Pal, Simmons,
Spieker, Chappell, Halasa, Kalams and the Influenza Vaccine in HCT Study Group;
preprint on medRxiv, 2026).

The study asks whether pre-vaccination immune-cell composition, measured by
30-marker mass cytometry, predicts the influenza antibody response of allogeneic
hematopoietic cell transplant (HCT) recipients beyond their baseline titers and
clinical factors. It uses the first vaccination season of 260 recipients (147
children, 113 adults) from two multicenter phase 2 trials of high- versus
standard-dose vaccine (NCT02860039, NCT03179761), and holds out 74 second-season
episodes for validation. The outcome is the composite log fold-rise in HAI titer
from before vaccination to 4 weeks after the second dose, averaged over H1N1,
H3N2 and B/Victoria.

## Pipeline

1. **Titers.** HAI titers are taken from the trial REDCap exports
   (`experiments/E17-07-redcap-titers.R`): geometric mean of duplicate wells,
   below-detection results set to 5.
2. **Per-fold atlases.** Recipients are split into five outer folds. In each
   outer training set, about a million cells are clustered into 5,000 k-means
   micro-clusters, and Leiden clustering of the micro-cluster graph gives
   candidate partitions over a grid of resolutions
   (`experiments/E14-01-first-episode-perfold.R`, `methods/35-per-fold-atlas-cv.R`).
   No outcome is used in this step.
3. **Model fitting.** Inside each training set, three inner folds choose the
   resolution, and elastic nets are fitted for four nested models: titers (M1),
   titers and clinical factors (M2), titers and composition (M3), and all three
   (M4). Held-out recipients' cells are assigned to the frozen atlas and their
   responses predicted (`experiments/E17-05-refit-missing-titers.R`).
4. **Inference and interpretation.** Recipient-level bootstrap intervals, a
   training-aware null bootstrap for the score tests, a shared 35-population
   atlas used only to name and compare populations, selection stability,
   per-cohort and transfer analyses, the second-season validation, marker
   reduction and the all-episode sensitivity analysis.

`run/run_pipeline.sh` runs every step in order; `run/env.sh` sets the data and
results folders. Figures are written to `paper/images/` and tables to
`paper/tables/`.

## Layout

- `methods/` shared functions: data preparation, two-stage clustering, per-fold
  atlas and cross-validation helpers, outcome and covariate handling.
- `experiments/` one script per analysis step, figure or table, named by the
  analysis round in which it was introduced.
- `run/` environment settings and the full run order.
- `notebooks_2026-07/` earlier R Markdown notebooks for the July 2026 analysis of
  343 episodes. They do not reproduce the current results and are kept for
  reference.

## Requirements

R (4.4) with data.table, here, glmnet, igraph, RANN, FNN, pROC, ggplot2,
patchwork, ggrepel, scales, uwot, randomForest, ranger, caret, sandwich,
lmtest, digest, doParallel and rjson; Python 3 with pandas and numpy for the
LaTeX table scripts.

## Data availability

The de-identified CyTOF data, trial data and intermediate analysis objects are
available from the corresponding author upon reasonable request. No
patient-derived data and no analysis outputs are included in this repository.

## License

Released under the MIT License (see `LICENSE`).
