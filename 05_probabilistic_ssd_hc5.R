# =============================================================================
# 05_probabilistic_ssd_hc5.R
# Probabilistic SSD and acute HC5 derivation
# 10,000 Monte Carlo iterations per PAH and region
# =============================================================================

suppressPackageStartupMessages({
  library(ssdtools)
  library(dplyr)
  library(readr)
})

OUT_DIR <- "outputs"

ENG_FILE <- file.path(OUT_DIR, "England_PAH_pred_with_uncertainty_Four4.csv")
JS_FILE  <- file.path(OUT_DIR, "Jiangsu_PAH_pred_with_uncertainty_Four4.csv")

N_MC <- 10000L
SEED <- 2026L

DISTS <- c("lnorm", "llogis", "burrIII3", "gamma", "weibull")

set.seed(SEED)

eng <- read_csv(ENG_FILE, show_col_types = FALSE)
js  <- read_csv(JS_FILE, show_col_types = FALSE)

required <- c("PAHs", "Pred_LogTox_Mean", "Sigma_Total")
stopifnot(all(required %in% names(eng)))
stopifnot(all(required %in% names(js)))

run_region <- function(dat, region) {
  summary_list <- list()
  draw_list <- list()

  for (pah in unique(dat$PAHs)) {
    d <- dat |>
      filter(
        PAHs == pah,
        is.finite(Pred_LogTox_Mean),
        is.finite(Sigma_Total)
      )

    if (nrow(d) < 5) next

    hc5 <- rep(NA_real_, N_MC)
    best_dist <- rep(NA_character_, N_MC)

    for (i in seq_len(N_MC)) {
      conc <- 10^rnorm(
        nrow(d),
        mean = d$Pred_LogTox_Mean,
        sd = d$Sigma_Total
      )

      sim <- data.frame(Conc = conc[is.finite(conc) & conc > 0])
      if (nrow(sim) < 5) next

      fits <- try(
        ssd_fit_dists(sim, dists = DISTS, silent = TRUE),
        silent = TRUE
      )
      if (inherits(fits, "try-error")) next

      gof <- try(ssd_gof(fits, wt = TRUE), silent = TRUE)
      if (inherits(gof, "try-error") || nrow(gof) == 0) next

      bd <- gof$dist[which.min(gof$aicc)]
      best_dist[i] <- bd

      h <- try(
        ssd_hc(fits, proportion = 0.05, average = FALSE) |>
          filter(dist == bd) |>
          pull(est),
        silent = TRUE
      )

      if (!inherits(h, "try-error") && length(h) > 0 && is.finite(h[1])) {
        hc5[i] <- h[1]
      }
    }

    valid <- is.finite(hc5) & hc5 > 0
    hv <- hc5[valid]
    dv <- best_dist[valid]

    freq <- sort(table(dv), decreasing = TRUE)

    summary_list[[pah]] <- tibble(
      Region = region,
      PAH = pah,
      N_species = nrow(d),
      HC5_Median = median(hv),
      HC5_LCL = quantile(hv, 0.025),
      HC5_UCL = quantile(hv, 0.975),
      Best_Dist = names(freq)[1],
      Win_Rate_Pct = 100 * as.numeric(freq[1]) / sum(freq),
      Success_Rate = 100 * length(hv) / N_MC
    )

    draw_list[[pah]] <- tibble(
      Region = region,
      PAH = pah,
      Iteration = which(valid),
      Acute_HC5 = hv,
      Best_Dist = dv
    )
  }

  list(
    summary = bind_rows(summary_list),
    draws = bind_rows(draw_list)
  )
}

res_eng <- run_region(eng, "England")
res_js  <- run_region(js, "Jiangsu")

summary_all <- bind_rows(res_eng$summary, res_js$summary)
draws_all <- bind_rows(res_eng$draws, res_js$draws)

write_csv(
  summary_all,
  file.path(OUT_DIR, "England_Jiangsu_PAH_HC5_Summary_10k_Median.csv")
)

# Saved because script 06 propagates ACR uncertainty using the actual HC5 draws.
write_csv(
  draws_all,
  file.path(OUT_DIR, "England_Jiangsu_PAH_HC5_MC_Draws_10k.csv")
)

message("05 complete.")
