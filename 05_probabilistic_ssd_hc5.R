# ==============================================================================
# 05_probabilistic_ssd_hc5.R
# ------------------------------------------------------------------------------
# Probabilistic SSD and HC5 derivation for England and Jiangsu.
# Computational settings are kept unchanged:
#   - dist_set = lnorm, llogis, burrIII3, gamma, weibull
#   - n_iter = 10000 Monte Carlo simulations
#   - set.seed(2026)
#   - ssdtools syntax updated as in the working script
# ==============================================================================

# ==============================================================================
# Probabilistic SSD / HC5 for England and Jiangsu
#
# Based on XGB 4-variable predictions with uncertainty
#
# Logic:
#   - Monte Carlo sample from Pred_LogTox_Mean and Sigma_Total
#   - Convert back to original concentration scale
#   - Fit 5 candidate SSD distributions
#   - Select best distribution by minimum AICc
#   - Extract HC5 for each iteration
#   - Summarise HC5 median, CI, best-dist win rate
#
# Updated ssdtools syntax:
#   - ssd_gof(fits, wt = TRUE)
#   - ssd_hc(fits, proportion = 0.05, average = FALSE)
# ==============================================================================

library(ssdtools)
library(dplyr)
library(tidyr)
library(readr)

# 1. Settings ------------------------------------------------------------------
# Set the working directory to the project root before running this script.
dist_set <- c("lnorm", "llogis", "burrIII3", "gamma", "weibull")
n_iter <- 10000

eng_file <- "England_PAH_pred_with_uncertainty_XGB4.csv"
js_file  <- "Jiangsu_PAH_pred_with_uncertainty_XGB4.csv"

save_aicc_logs <- FALSE
set.seed(2026)

# 2. Load prediction files ------------------------------------------------------

eng_raw <- read_csv(eng_file, show_col_types = FALSE)
js_raw  <- read_csv(js_file,  show_col_types = FALSE)

required_cols <- c("PAHs", "Latin", "Pred_LogTox_Mean", "Sigma_Total")

if (!all(required_cols %in% names(eng_raw))) {
  stop(
    "England file missing required columns: ",
    paste(setdiff(required_cols, names(eng_raw)), collapse = ", ")
  )
}

if (!all(required_cols %in% names(js_raw))) {
  stop(
    "Jiangsu file missing required columns: ",
    paste(setdiff(required_cols, names(js_raw)), collapse = ", ")
  )
}

# 3. Core function --------------------------------------------------------------

run_prob_ssd_hc5 <- function(raw_data,
                             region_name,
                             dist_set,
                             n_iter,
                             save_aicc_logs = FALSE) {
  
  pah_list <- unique(raw_data$PAHs)
  
  all_pah_hc5_stats <- list()
  all_aicc_logs <- list()
  
  for (p in pah_list) {
    
    message("\n>>> [", region_name, "] analysing PAH: ", p)
    
    sub_df <- raw_data %>%
      filter(PAHs == p) %>%
      filter(is.finite(Pred_LogTox_Mean), is.finite(Sigma_Total))
    
    n_spp <- nrow(sub_df)
    
    if (n_spp < 5) {
      warning(region_name, " / ", p, ": fewer than 5 records; skipped.")
      next
    }
    
    iter_hc5 <- numeric(n_iter)
    iter_best_dist <- character(n_iter)
    iter_aicc_list <- vector("list", n_iter)
    
    pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
    
    for (i in seq_len(n_iter)) {
      
      # Monte Carlo sampling on log10 scale, then back-transform
      sim_data <- data.frame(
        Conc = 10^rnorm(
          n    = n_spp,
          mean = sub_df$Pred_LogTox_Mean,
          sd   = sub_df$Sigma_Total
        )
      )
      
      # Defensive filter: keep only finite positive concentrations
      sim_data <- sim_data %>%
        filter(is.finite(Conc), Conc > 0)
      
      round_log <- data.frame(
        Region = region_name,
        PAH = p,
        Iteration = i,
        dist = dist_set,
        aicc = NA_real_,
        Fit_Status = "Failed",
        stringsAsFactors = FALSE
      )
      
      if (nrow(sim_data) >= 5) {
        
        try({
          
          fits <- ssd_fit_dists(
            sim_data,
            dists = dist_set,
            silent = TRUE
          )
          
          # Updated syntax for ssdtools >= 2.3.1
          gof <- ssd_gof(
            fits,
            wt = TRUE
          )
          
          # Update AICc log
          for (d in gof$dist) {
            round_log$aicc[round_log$dist == d] <- gof$aicc[gof$dist == d]
            round_log$Fit_Status[round_log$dist == d] <- "Success"
          }
          
          # Choose best distribution by minimum AICc
          if (nrow(gof) > 0) {
            
            best_d <- gof$dist[which.min(gof$aicc)]
            iter_best_dist[i] <- best_d
            
            # Updated syntax for ssdtools >= 2.0.0
            hc5_val <- ssd_hc(
              fits,
              proportion = 0.05,
              average = FALSE
            ) %>%
              filter(dist == best_d) %>%
              pull(est)
            
            if (length(hc5_val) > 0 &&
                is.finite(hc5_val[1]) &&
                hc5_val[1] > 0) {
              iter_hc5[i] <- hc5_val[1]
            }
          }
          
        }, silent = TRUE)
      }
      
      iter_aicc_list[[i]] <- round_log
      setTxtProgressBar(pb, i)
    }
    
    close(pb)
    
    # AICc logs
    pah_aicc_full <- bind_rows(iter_aicc_list)
    all_aicc_logs[[p]] <- pah_aicc_full
    
    # Distribution win rates
    dist_freq <- as.data.frame(table(iter_best_dist[iter_best_dist != ""]))
    colnames(dist_freq) <- c("dist", "Count")
    
    if (nrow(dist_freq) > 0) {
      dist_freq <- dist_freq %>%
        mutate(Prop = Count / sum(Count) * 100) %>%
        arrange(desc(Prop))
    } else {
      dist_freq <- data.frame(
        dist = NA_character_,
        Count = NA_integer_,
        Prop = NA_real_
      )
    }
    
    # Valid HC5 values
    v_hc5 <- iter_hc5[is.finite(iter_hc5) & iter_hc5 > 0]
    
    if (length(v_hc5) == 0) {
      
      all_pah_hc5_stats[[p]] <- data.frame(
        Region = region_name,
        PAH = p,
        N_species = n_spp,
        HC5_Median = NA_real_,
        HC5_LCL = NA_real_,
        HC5_UCL = NA_real_,
        Best_Dist = NA_character_,
        Win_Rate_Pct = NA_real_,
        Second_Best = NA_character_,
        Second_Rate_Pct = NA_real_,
        Success_Rate = 0,
        stringsAsFactors = FALSE
      )
      
    } else {
      
      all_pah_hc5_stats[[p]] <- data.frame(
        Region = region_name,
        PAH = p,
        N_species = n_spp,
        HC5_Median = median(v_hc5),
        HC5_LCL = quantile(v_hc5, 0.025, na.rm = TRUE),
        HC5_UCL = quantile(v_hc5, 0.975, na.rm = TRUE),
        Best_Dist = dist_freq$dist[1],
        Win_Rate_Pct = dist_freq$Prop[1],
        Second_Best = ifelse(
          nrow(dist_freq) > 1,
          as.character(dist_freq$dist[2]),
          NA
        ),
        Second_Rate_Pct = ifelse(
          nrow(dist_freq) > 1,
          dist_freq$Prop[2],
          NA
        ),
        Success_Rate = length(v_hc5) / n_iter * 100,
        stringsAsFactors = FALSE
      )
    }
  }
  
  final_summary <- bind_rows(all_pah_hc5_stats)
  
  out <- list(
    summary = final_summary,
    aicc_logs = bind_rows(all_aicc_logs)
  )
  
  if (save_aicc_logs) {
    write_csv(
      out$aicc_logs,
      paste0(region_name, "_PAH_AICc_Full_10k.csv")
    )
  }
  
  write_csv(
    out$summary,
    paste0(region_name, "_PAH_HC5_Summary_10k_Median.csv")
  )
  
  return(out)
}

# 4. Run England ---------------------------------------------------------------

eng_res <- run_prob_ssd_hc5(
  raw_data = eng_raw,
  region_name = "England",
  dist_set = dist_set,
  n_iter = n_iter,
  save_aicc_logs = save_aicc_logs
)

# 5. Run Jiangsu ---------------------------------------------------------------

js_res <- run_prob_ssd_hc5(
  raw_data = js_raw,
  region_name = "Jiangsu",
  dist_set = dist_set,
  n_iter = n_iter,
  save_aicc_logs = save_aicc_logs
)

# 6. Combined summary ----------------------------------------------------------

final_summary_all <- bind_rows(
  eng_res$summary,
  js_res$summary
)

write_csv(
  final_summary_all,
  "England_Jiangsu_PAH_HC5_Summary_10k_Median.csv"
)

cat("\n>>> Done.\n")
cat("Saved files:\n")
cat(" - England_PAH_HC5_Summary_10k_Median.csv\n")
cat(" - Jiangsu_PAH_HC5_Summary_10k_Median.csv\n")
cat(" - England_Jiangsu_PAH_HC5_Summary_10k_Median.csv\n")

print(final_summary_all)