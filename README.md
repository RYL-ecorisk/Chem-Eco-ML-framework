# Chemical–Ecological ML Framework

Public R code supporting a machine-learning framework for deriving assemblage-specific freshwater protection thresholds under toxicity-data scarcity.

## Background

Conventional ecological protection thresholds are often derived from a limited set of surrogate species and may not fully represent the sensitivity of regional freshwater assemblages. This framework combines chemical descriptors and ecological traits to extend toxicity prediction across chemicals and taxa, then links those predictions with probabilistic species sensitivity distributions (pSSDs) to derive regional HC5 values.

Polycyclic aromatic hydrocarbons (PAHs) and freshwater macroinvertebrates are used as a case study. The workflow compares multiple machine-learning algorithms, streamlines the best-performing model to a small set of interpretable predictors, applies it to regional species pools, and propagates uncertainty through toxicity prediction, SSD fitting, and acute-to-chronic extrapolation. The overall aim is to provide a transparent and transferable approach for region-specific ecological threshold derivation when direct toxicity data are sparse.

## Modules

1. `01_full_model_training_comparison.R`  
   Train and compare Elastic Net, SVM, Random Forest, and XGBoost models.

2. `02_model_interpretation_shap.R`  
   Interpret the final XGBoost model using native TreeSHAP.

3. `03_streamlined_model_selection.R`  
   Streamline the full model directly to the final four-predictor XGBoost model.

4. `04_regional_prediction_bootstrap_uncertainty.R`  
   Generate regional toxicity predictions and quantify prediction uncertainty by bootstrap resampling.

5. `05_probabilistic_ssd_hc5.R`  
   Construct probabilistic species sensitivity distributions and derive acute HC5 values using 10,000 Monte Carlo simulations.

6. `06_ssd_figures.R`  
   Produce the main-text and supplementary acute/chronic SSD figures.

7. `07_acr_model.R`  
   Fit the Gamma GLM for acute-to-chronic ratios and obtain final ACR values.

8. `08_uncertainty_propagation.R`  
   Convert acute HC5 values to chronic HC5 values and propagate uncertainty associated with predicted ACRs.

## Notes

- Toxicity is modelled on the log10 concentration scale.
- The streamlined XGBoost model uses water solubility, maximum body size, log Kow, and feeding group.
- Regional prediction uncertainty is estimated using bootstrap resampling together with repeated-cross-validation error.
- Acute HC5 values are derived from 10,000 Monte Carlo probabilistic SSD simulations.
- ACRs are modelled with a Gamma GLM using molecular complexity as the predictor.
- Chronic HC5 is calculated as acute HC5 divided by the final ACR.
- Predicted-ACR uncertainty is propagated when estimating chronic HC5 uncertainty.
