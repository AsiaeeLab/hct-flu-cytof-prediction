# Deep immune cell phenotyping predicts influenza vaccine immunogenicity in HCT recipients

Analysis code accompanying the manuscript
*Deep Immune Cell Phenotyping Predicts Influenza Vaccine Immunogenicity in HCT
Recipients.*

This repository holds the R Markdown notebooks that regenerate every figure and
table in the manuscript and its supplement, together with the knitted HTML
reports. The patient-derived cytometry data and the precomputed intermediate
objects are not distributed here; see [Data availability](#data-availability).

## Overview

We analyzed 343 influenza-vaccination episodes from 263 allogeneic
hematopoietic-cell-transplant (HCT) recipients enrolled in two multi-center,
NIH-funded Phase II randomized trials of high- versus standard-dose inactivated
influenza vaccine (one pediatric, one adult). Pre-vaccination blood was profiled
by 30-marker mass cytometry (CyTOF). The analysis converts each single-cell
profile into a centered-log-ratio cellular-composition vector using a
recipient-grouped, per-fold atlas pipeline (feature construction is kept strictly
separate from prediction), and quantifies the contribution of cellular
composition beyond the pre-vaccination HAI titer and clinical covariates with a
four-model elastic-net ladder and recipient-clustered bootstrap inference.

The notebooks are organized by theme. Each one is commented for a general
scientific reader: it explains, in plain language, what each step does and what
the output means, so that a hematologist or immunologist can follow the logic
without a machine-learning background.

## Notebook-to-figure map

Every notebook is designed to knit by **loading** the precomputed intermediate
objects rather than re-running the expensive clustering, cross-validation, or
bootstrap. Compute-intensive steps are guarded so that, when the cached objects
are present, the notebook simply loads them.

| Notebook | Manuscript outputs regenerated |
|----------|--------------------------------|
| `01-atlas-construction-and-k-selection.Rmd` | **Figure 1** (atlas pipeline and resolution/k-selection across folds); production-atlas summary |
| `02-four-model-ladder-and-inference.Rmd` | **Figure 2** (four-model ladder with bootstrap CIs); **Table 1** (ladder MSE); **Table 2** (recipient-clustered bootstrap); **Table 3** (baseline-titer stratification); **Supplementary Table** (HAI noise floor) |
| `03-cluster-biology-and-annotation.Rmd` | **Figure 3** (annotated UMAP); **Supplementary Figures** (faceted UMAP; biclustered marker heatmap); **Supplementary Table** (meta-cluster annotations); cluster-vs-time-post-transplant regression |
| `04-feature-importance-and-stability.Rmd` | **Figure 4** (signed feature importance, M4); **Supplementary Figures** (bootstrap stability landscape; selection frequency; coefficient distributions) |
| `05-per-cohort-and-transfer.Rmd` | **Figure 5** (per-cohort gain, feature selection, and stratification); **Supplementary Figure** (cross-cohort transfer) |
| `06-repeater-longitudinal.Rmd` | **Figure 6** (repeater antigen-count transitions); **Supplementary Figures** (reconstitution schematic; PCA trajectories; per-cluster change; cross-sectional predicted-response scatter) |
| `07-marker-panel-reduction.Rmd` | **Figure 7** (15-marker panel reduction) |
| `08-sensitivity-and-dose-supplement.Rmd` | **Supplementary Figures** (multitask vs single-task; 2,000-vs-5,000 micro-cluster sensitivity; counterfactual dose subgroups); one-episode-per-recipient and interaction-rich-baseline sensitivity analyses |

Together the eight notebooks regenerate all 7 main figures, all 13 supplementary
figures, the 3 main tables, and both supplementary tables (the meta-cluster
annotation table and the noise-floor table).

## Repository contents

- `*.Rmd` -- the eight analysis notebooks (source).
- `*.html` -- the knitted reports. These are the rendered outputs; because the
  data is not distributed, the HTML files are the way to see the figures and
  tables without access to the intermediate objects.
- `README.md`, `LICENSE`, `.gitignore`.

No data files are included. The `.gitignore` explicitly excludes `*.RData`,
`*.rds`, `*.csv`, `*.fcs`, the `/fast` data drive, and any local results
directories, so that no patient-derived or intermediate data can be committed by
accident.

## Environment and dependencies

The analyses were performed in R 4.5.2; the HTML reports in this repository were
last knitted under R 4.6.1. Rendering the notebooks with the cached intermediate
objects present requires:

- `rmarkdown`, `knitr` -- notebook rendering
- `ggplot2`, `patchwork`, `scales` -- figures
- `data.table` -- data manipulation
- `randomForest` -- the counterfactual dose T-learner (notebook 08)

Reproducing the intermediate objects from raw data (the guarded "heavy compute"
branches) additionally requires `glmnet` (elastic net), `igraph` (Leiden
community detection), and `RANN` (nearest-neighbor cell assignment), along with
the atlas helper functions from the project code and the raw CyTOF object.

Approximate versions used for the knitted reports: R 4.6.1, ggplot2 4.0.3,
patchwork 1.3.2, data.table 1.18.4, scales 1.4.0, randomForest 4.7.1.2,
knitr 1.51, rmarkdown 2.31 (glmnet 5.0 and igraph 2.3.2 for full recomputation).

## How to run

1. Obtain the intermediate objects (see Data availability) and place them on a
   local drive.
2. At the top of each notebook, edit the two path variables (`FAST` and `REPO`)
   to point at your copy of the intermediate data. The notebooks reference these
   paths but do not copy any data into this repository.
3. Knit a notebook with, for example:

   ```r
   rmarkdown::render("02-four-model-ladder-and-inference.Rmd")
   ```

   With the cached objects present, each notebook loads them and produces the
   figures, tables, and a complete HTML report. It does not re-run the
   clustering, cross-validation, or bootstrap.

## Data availability

The de-identified CyTOF data and the precomputed intermediate objects that
support the findings of this study are available from the corresponding author
upon reasonable request. Consistent with the manuscript's Data Sharing
Statement, no patient-derived data and no intermediate analysis objects are
included in this repository; only the source notebooks and their knitted HTML
reports are provided.

## License

Released under the MIT License (see `LICENSE`).
