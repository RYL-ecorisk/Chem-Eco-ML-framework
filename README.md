# Chemical–Ecological ML Framework

Public R code for the PAH toxicity-prediction and HC5 workflow.

## Modules

1. `01_full_model_training_comparison.R`  
   Train and compare Elastic Net, SVM, Random Forest and XGBoost.

2. `02_model_interpretation_shap.R`  
   Native TreeSHAP interpretation of the final XGBoost model.

3. `03_streamlined_model_selection.R`  
   Streamline Full12 directly to Four4: Solubility, Size, log.Kow and Feeding.

4. `04_regional_prediction_bootstrap_uncertainty.R`  
   Predict toxicity for England and Jiangsu with bootstrap uncertainty.

5. `05_probabilistic_ssd_hc5.R`  
   Run 10,000 probabilistic SSD simulations and derive acute HC5.

6. `06_ssd_figures.R`  
   Produce the main-text and supplementary Jiangsu acute/chronic SSD figures using 10,000 simulated curves.

7. `07_acr_model.R`  
   Fit the Gamma GLM for ACR using molecular Complexity and generate final ACR values.

8. `08_uncertainty_propagation.R`  
   Convert acute HC5 to chronic HC5 and propagate predicted-ACR uncertainty.

## Core inputs

- `datasets.csv`
- `England_PAH_input_4vars.csv`
- `Jiangsu_PAH_input_4vars.csv`
- `Observed_ACR.csv`
- `PAHs4.20.csv`

## Notes

- Toxicity is modelled on the log10 concentration scale.
- Four4 uses Solubility, Size, log.Kow and Feeding.
- Acute HC5 is based on 10,000 Monte Carlo pSSD simulations.
- ACR is modelled with a Gamma GLM using molecular Complexity.
- Chronic HC5 = Acute HC5 / Final ACR.
- Module 06 uses `Final_ACR_all16.csv`; therefore run Module 07 before Module 06 when reproducing the figures.
