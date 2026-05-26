# PAH regional threshold derivation: R code package

This folder contains cleaned, manuscript-aligned R scripts for the PAH machine-learning, model interpretation, regional prediction, probabilistic SSD/HC5 derivation and uncertainty-propagation workflow. The scripts were organized from the working analysis files without changing the main computational parameters.

## Included files

- `01_full_model_training_comparison.R` — trains Elastic Net, SVM, RF and XGBoost models and exports model-performance summaries, hyperparameters, selected model and prediction files.
- `02_model_interpretation_shap_pdp.R` — performs SHAP, PDP and feeding-group interpretation. Run after Script 01 in the same R session if the trained model object is needed.
- `03_streamlining_m4_m8_m12.R` — compares M4_Core, M8_Extended and M12_Full models and exports the final four-variable prediction bundle.
- `04_regional_prediction_bootstrap_uncertainty.R` — applies the final M4_Core XGBoost model to regional assemblages and propagates prediction uncertainty by bootstrap refitting.
- `05_probabilistic_ssd_hc5.R` — derives regional acute HC5 values using probabilistic SSDs and Monte Carlo simulation.
- `06_ssd_figures_and_hc5_tables.R` — generates SSD figure panels and HC5/ACR/solubility summary tables used for manuscript and Supplementary Information preparation.

## Data availability

This code-only package does not include the processed toxicity dataset, regional assemblage prediction inputs, ACR tables, HC5 summaries, PAH property tables, or monitoring/risk data. These data are reported or summarized in the manuscript and Supplementary Information, and expected input file names are listed at the top of each script.

Users who wish to rerun the workflow should place the required processed input files in the working directory using the file names specified in each script.

## Main computational settings retained

- Full-model training: repeated 10-fold cross-validation with 10 repeats.
- Full-model algorithms: Elastic Net, SVM, RF and XGBoost.
- Bayesian optimization settings and search bounds were retained from the working scripts.
- Streamlining feature sets were retained:
  - `M4_Core`: Solubility, Size, log.Kow, Feeding.
  - `M8_Extended`: Solubility, log.Kow, Complexity, Respiration, Locomotion, Feeding, Size, HLC.
  - `M12_Full`: all 12 predictors.
- Regional prediction uncertainty: 1,000 bootstrap refits (`B = 1000`).
- pSSD/HC5 simulation: 10,000 Monte Carlo iterations (`n_iter = 10000`).
- Candidate SSD distributions: log-normal, log-logistic, Burr type III, gamma and Weibull.

## Suggested run order

1. Run `01_full_model_training_comparison.R` using the processed model-training dataset (`datasets.csv`).
2. Run `02_model_interpretation_shap_pdp.R` in the same R session if SHAP/PDP outputs are needed.
3. Run `03_streamlining_m4_m8_m12.R` to generate the final M4_Core prediction bundle.
4. Prepare the regional prediction input files in the working directory:
   - `England_PAH_input_4vars.csv`
   - `Jiangsu_PAH_input_4vars.csv`
5. Run `04_regional_prediction_bootstrap_uncertainty.R`.
6. Run `05_probabilistic_ssd_hc5.R`.
7. Run `06_ssd_figures_and_hc5_tables.R` for SSD figures and HC5-related summary tables.

## Notes

- The cleaned scripts remove hard-coded local `setwd()` calls. Set the working directory to this folder before running.
- File names were kept close to the working analysis to avoid changing the computational workflow.
- This repository is intended to document the main analytical workflow. Processed input data and result tables are provided or summarized separately in the manuscript and Supplementary Information.
