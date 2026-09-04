# =============================================================================
# 08_uncertainty_propagation.R
# Acute-to-chronic HC5 conversion and total uncertainty propagation
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

OUT_DIR <- "outputs"
Z95 <- 1.96

HC5_FILE <- file.path(
  OUT_DIR,
  "England_Jiangsu_PAH_HC5_Summary_10k_Median.csv"
)

ACR_FILE <- file.path(
  OUT_DIR,
  "Final_ACR_all16.csv"
)

ACR_MODEL_FILE <- file.path(
  OUT_DIR,
  "Gamma_GLM_ACR_model.rds"
)

hc5 <- read_csv(
  HC5_FILE,
  show_col_types = FALSE
)

acr <- read_csv(
  ACR_FILE,
  show_col_types = FALSE
)

acr_model <- readRDS(
  ACR_MODEL_FILE
)

gamma_phi <- summary(acr_model)$dispersion

sigma_ln_acr_model <- sqrt(
  log(1 + gamma_phi)
)

# -----------------------------------------------------------------------------
# 1. Direct acute-to-chronic conversion
# -----------------------------------------------------------------------------

hc5_table <- hc5 |>
  left_join(
    acr |>
      transmute(
        PAH = Abb,
        Name,
        ACR_source,
        Final_ACR
      ),
    by = "PAH"
  ) |>
  mutate(
    Acute_HC5 = HC5_Median,
    Acute_HC5_LCL = HC5_LCL,
    Acute_HC5_UCL = HC5_UCL,

    Chronic_HC5 =
      Acute_HC5 / Final_ACR,

    Chronic_HC5_LCL =
      Acute_HC5_LCL / Final_ACR,

    Chronic_HC5_UCL =
      Acute_HC5_UCL / Final_ACR
  )

write_csv(
  hc5_table,
  file.path(
    OUT_DIR,
    "HC5_Acute_Chronic_Table.csv"
  )
)

# -----------------------------------------------------------------------------
# 2. Total chronic HC5 uncertainty
# -----------------------------------------------------------------------------
# Acute uncertainty is inherited from the 10,000 pSSD Monte Carlo simulations.
#
# For predicted ACRs:
#   sigma_ln_ACR = sqrt(log(1 + phi))
#
# For observed ACRs:
#   no additional ACR-model uncertainty is added.
#
# Independent components are combined on the natural-log scale.
# -----------------------------------------------------------------------------

uncertainty <- hc5_table |>
  mutate(
    Acute_fold =
      Acute_HC5_UCL / Acute_HC5_LCL,

    sigma_ln_acute =
      log(Acute_fold) / (2 * Z95),

    sigma_ln_acr = if_else(
      ACR_source == "Predicted",
      sigma_ln_acr_model,
      0
    ),

    ACR_fold =
      exp(
        2 * Z95 * sigma_ln_acr
      ),

    sigma_ln_chronic =
      sqrt(
        sigma_ln_acute^2 +
          sigma_ln_acr^2
      ),

    Chronic_factor =
      exp(
        Z95 * sigma_ln_chronic
      ),

    Chronic_HC5_total_LCL =
      Chronic_HC5 / Chronic_factor,

    Chronic_HC5_total_UCL =
      Chronic_HC5 * Chronic_factor,

    Chronic_fold =
      exp(
        2 * Z95 * sigma_ln_chronic
      )
  ) |>
  select(
    Region,
    PAH,
    Name,
    ACR_source,
    Final_ACR,

    Acute_HC5,
    Acute_HC5_LCL,
    Acute_HC5_UCL,
    Acute_fold,

    ACR_fold,

    Chronic_HC5,
    Chronic_HC5_total_LCL,
    Chronic_HC5_total_UCL,
    Chronic_fold
  )

write_csv(
  uncertainty,
  file.path(
    OUT_DIR,
    "HC5_ACR_Uncertainty.csv"
  )
)

message("08 complete.")
