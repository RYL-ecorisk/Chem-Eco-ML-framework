# =============================================================================
# 07_acr_model.R
# Gamma GLM for PAH acute-to-chronic ratios (ACRs)
# Predictor: molecular Complexity
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

OBSERVED_FILE <- "Observed_ACR.csv"
PAH_FILE <- "PAHs4.20.csv"
OUT_DIR <- "outputs"

observed <- read_csv(
  OBSERVED_FILE,
  show_col_types = FALSE
)

pahs <- read_csv(
  PAH_FILE,
  show_col_types = FALSE
)

stopifnot(
  all(c("PAH_code", "Observed_ACR") %in% names(observed))
)

stopifnot(
  all(c("CAS number", "Name", "Abb", "Complexity") %in% names(pahs))
)

observed <- observed |>
  transmute(
    PAH_code = trimws(PAH_code),
    Observed_ACR = as.numeric(Observed_ACR)
  ) |>
  filter(
    is.finite(Observed_ACR),
    Observed_ACR > 0
  )

pahs <- pahs |>
  mutate(
    Abb = trimws(Abb),
    Complexity = as.numeric(Complexity)
  )

# Multiple observed ACR records for one PAH:
# use the arithmetic mean, matching the final analysis.
observed_summary <- observed |>
  group_by(PAH_code) |>
  summarise(
    Observed_ACR = mean(Observed_ACR),
    n_records = n(),
    .groups = "drop"
  )

acr_data <- pahs |>
  left_join(
    observed_summary,
    by = c("Abb" = "PAH_code")
  )

model_data <- acr_data |>
  filter(
    is.finite(Observed_ACR),
    Observed_ACR > 0,
    is.finite(Complexity)
  )

gamma_model <- glm(
  Observed_ACR ~ Complexity,
  data = model_data,
  family = Gamma(link = "log")
)

gamma_phi <- summary(gamma_model)$dispersion

acr_data <- acr_data |>
  mutate(
    Predicted_ACR = predict(
      gamma_model,
      newdata = acr_data,
      type = "response"
    ),
    ACR_source = if_else(
      is.finite(Observed_ACR),
      "Observed",
      "Predicted"
    ),
    Final_ACR = if_else(
      ACR_source == "Observed",
      Observed_ACR,
      Predicted_ACR
    )
  )

final_acr <- acr_data |>
  select(
    `CAS number`,
    Name,
    Abb,
    Complexity,
    Observed_ACR,
    n_records,
    Predicted_ACR,
    Final_ACR,
    ACR_source
  )

write_csv(
  final_acr,
  file.path(
    OUT_DIR,
    "Final_ACR_all16.csv"
  )
)

saveRDS(
  gamma_model,
  file.path(
    OUT_DIR,
    "Gamma_GLM_ACR_model.rds"
  )
)

model_info <- tibble(
  Metric = c(
    "Intercept",
    "Complexity_slope",
    "Gamma_dispersion_phi"
  ),
  Value = c(
    coef(gamma_model)[1],
    coef(gamma_model)[2],
    gamma_phi
  )
)

write_csv(
  model_info,
  file.path(
    OUT_DIR,
    "Gamma_GLM_ACR_model_info.csv"
  )
)

p <- ggplot(
  model_data,
  aes(Complexity, Observed_ACR)
) +
  geom_point(size = 2.4) +
  stat_function(
    fun = function(x) {
      exp(
        coef(gamma_model)[1] +
          coef(gamma_model)[2] * x
      )
    },
    linewidth = 0.9
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = "Molecular complexity",
    y = "Acute-to-chronic ratio (ACR)"
  )

ggsave(
  file.path(
    OUT_DIR,
    "Figure_ACR_Gamma_GLM.png"
  ),
  p,
  width = 6.0,
  height = 4.6,
  dpi = 600
)

message(
  "07 complete. Gamma dispersion phi = ",
  signif(gamma_phi, 6)
)
