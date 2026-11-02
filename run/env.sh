# Locations used by every script. Edit these for your system; the cell-level
# object and trial files they point to are not distributed (see README).
# methods/00-paths.R also reads ~/Paths/flu_cellular_prediction.json for the
# raw-data folder.
export FLU_ATLAS_ROOT=/fast/datasets/flu_cellular_prediction/first_episode_2026-10-14
export FLU_RESULTS_ROOT=/fast/datasets/flu_cellular_prediction/composite_2026-10-29
export FLU_SHARED_ATLAS_DIR=$FLU_RESULTS_ROOT/shared_atlas
export FLU_MEMO_DIR=$FLU_RESULTS_ROOT/post_hoc_5000/memo
export FLU_PROD_ATLAS=$FLU_SHARED_ATLAS_DIR/final_atlas_d4b296c0.RData
export FLU_BIO_FILE=$FLU_SHARED_ATLAS_DIR/biology_production_cbc9036d.RData
export FLU_ANN_FILE=$FLU_SHARED_ATLAS_DIR/cluster_annotations.csv
export FLU_BOOT_FILE=$FLU_SHARED_ATLAS_DIR/bootstrap_importance_c3c56f77.RData
export FLU_BOOT_FILE_COMP=$FLU_SHARED_ATLAS_DIR/bootstrap_importance_7644c0b1.RData
export FLU_ALLEP_ROOT=/fast/datasets/flu_cellular_prediction/timing_corrected_2026-10-29
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
