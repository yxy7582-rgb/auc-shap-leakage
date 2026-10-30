# AUC–SHAP dissociation under preprocessing leakage

R scripts for the simulation and clinical analyses in the accompanying manuscript.

## Scripts
- `simulate_compare.R` — Experiment 1 (redundancy × leakage)
- `simulate_selectors_bin=30.R` — Experiments 2–3 (leakage-strength sweep, n=30/bin; selectors)
- `real_data_test.R` — Real-world dermatomyositis/CADM cohort (ILD endpoint)

## Data
Simulation data are generated in-script. Clinical data are not shared publicly; contact the corresponding author with ethics approval. Place approved data at `data/clinical.xlsx` or set `CLINICAL_DATA_PATH`.

## Requirements
R packages: tidyverse, ranger, fastshap, pROC, glmnet, caret, readxl, ggplot2, patchwork
