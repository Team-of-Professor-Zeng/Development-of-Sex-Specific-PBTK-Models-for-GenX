# =====================================================================
# Compute fitting error metrics by sex / dose / tissue
# Metrics: RMSE, MAE, AFE (Average Fold Error), AAFE (Absolute Average Fold Error)
#
# AFE / AAFE use log-scale definitions (common in PBPK):
#   AFE  = 10^mean(log10(PRE/OBS))      -> <1 under-predicted, >1 over-predicted, =1 perfect
#   AAFE = 10^mean(|log10(PRE/OBS)|)    -> >=1, =1 perfect, larger means more deviation
# =====================================================================

library(dplyr)

setwd("c:/Users/15960/Desktop/Modfit/Mice")

## ---- 1. Read data ----
Male_raw   <- read.csv("Male/Male.csv")            # raw observation data (incl. Dose)
Female_raw <- read.csv("Female/Female.csv")

FDataA_M <- read.csv("Male/FDataA_M.csv",   row.names = 1)   # fitting residuals (incl. OBS/PRE/log.OPR)
FDataA_F <- read.csv("Female/FDataA_F.csv", row.names = 1)

FDataA_M$Gender <- "Male"
FDataA_F$Gender <- "Female"

## ---- 2. Match back to raw data via OBS value to get Dose / Study / original Sample ----
match_meta <- function(fdata, raw) {
  Dose      <- rep(NA_real_, nrow(fdata))
  Study     <- rep(NA_real_, nrow(fdata))
  RawSample <- rep(NA_character_, nrow(fdata))
  for (i in seq_len(nrow(fdata))) {
    d <- abs(raw$Conc - fdata$OBS[i])
    j <- which.min(d)
    if (d[j] < 1e-6) {
      Dose[i]      <- raw$Dose[j]
      Study[i]     <- raw$Study[j]
      RawSample[i] <- as.character(raw$Sample[j])
    }
  }
  fdata$Dose      <- Dose
  fdata$Study     <- Study
  fdata$RawSample <- RawSample
  fdata
}

FDataA_M <- match_meta(FDataA_M, Male_raw)
FDataA_F <- match_meta(FDataA_F, Female_raw)

# Verify matching completeness
n_miss_M <- sum(is.na(FDataA_M$Dose))
n_miss_F <- sum(is.na(FDataA_F$Dose))
if (n_miss_M) warning("Male: ", n_miss_M, " points not matched to dose info")
if (n_miss_F) warning("Female: ", n_miss_F, " points not matched to dose info")

FDataA_all <- rbind(FDataA_M, FDataA_F)

## ---- 3. Metric calculation function ----
calc_metrics <- function(df) {
  obs <- df$OBS
  pre <- df$PRE
  n   <- length(obs)

  # RMSE / MAE on log10 scale
  log_pre <- log10(pre)
  log_obs <- log10(obs)
  valid_log <- is.finite(log_pre) & is.finite(log_obs)
  if (sum(valid_log) > 0) {
    RMSE <- sqrt(mean((log_pre[valid_log] - log_obs[valid_log])^2))
    MAE  <- mean(abs(log_pre[valid_log] - log_obs[valid_log]))
  } else {
    RMSE <- MAE <- NA
  }

  # Use existing log.OPR column ( = log10(PRE/OBS) )
  valid <- is.finite(df$log.OPR)
  if (sum(valid) > 0) {
    AFE  <- 10 ^ mean(df$log.OPR[valid])
    AAFE <- 10 ^ mean(abs(df$log.OPR[valid]))
  } else {
    AFE <- AAFE <- NA
  }

  # R² (adjusted) based on log10-log10 linear regression, consistent with Figure 2
  if (n >= 3 && sum(is.finite(df$Log.OBS) & is.finite(df$Log.PRE)) >= 3) {
    fit <- lm(Log.OBS ~ Log.PRE, data = df)
    R2  <- summary(fit)$adj.r.squared
  } else {
    R2 <- NA
  }

  # Percentage of predictions within 2-fold of observations (0.5 <= PRE/OBS <= 2)
  if (sum(valid) > 0) {
    pct_2fold <- mean(df$OPR[valid] >= 1/2 & df$OPR[valid] <= 2) * 100
  } else {
    pct_2fold <- NA
  }

  data.frame(N = n, RMSE = RMSE, MAE = MAE, AFE = AFE, AAFE = AAFE, R2 = R2,
             pct_2fold = pct_2fold)
}

## ---- 4. Grouped calculations (single dimension) ----
# 4a. By sex x dose (pooled across tissues)
results_dose <- FDataA_all %>%
  group_by(Gender, Dose) %>%
  group_modify(~calc_metrics(.x)) %>%
  ungroup() %>%
  arrange(Gender, Dose)

# 4b. By sex x tissue (pooled across doses)
results_tissue <- FDataA_all %>%
  group_by(Gender, Tissue = name) %>%
  group_modify(~calc_metrics(.x)) %>%
  ungroup() %>%
  arrange(Gender, Tissue)

# 4c. By sex overall
results_gender <- FDataA_all %>%
  group_by(Gender) %>%
  group_modify(~calc_metrics(.x)) %>%
  ungroup()

## ---- 5. Output ----
write.csv(results_dose,   "Error_metrics_dose.csv",   row.names = FALSE)
write.csv(results_tissue, "Error_metrics_tissue.csv", row.names = FALSE)
write.csv(results_gender, "Error_metrics_gender.csv", row.names = FALSE)

## =====================================================================
## 6. Error metrics for model evaluation data (evadata.csv)
##    Trailing columns correspond to: Study, Time, Dose, Conc(=OBS), SD, Sample
##    Key columns: V2=PRE, V3=OBS, V1=Tissue (model output name), Dose=V[ncol-3]
## =====================================================================

eva_M_raw <- read.csv("Male/evadata.csv",   skip = 1, header = FALSE,
                       stringsAsFactors = FALSE, fill = TRUE)
eva_F_raw <- read.csv("Female/evadata.csv", skip = 1, header = FALSE,
                       stringsAsFactors = FALSE, fill = TRUE)

process_eva <- function(eva, gender) {
  nc <- ncol(eva)
  data.frame(
    Gender = gender,
    name   = eva[[1]],                          # model output channel (Plasma/Liver)
    OBS    = as.numeric(eva[[3]]),
    PRE    = as.numeric(eva[[2]]),
    Dose   = as.numeric(eva[[nc - 3]]),
    Study  = as.numeric(eva[[nc - 5]]),
    stringsAsFactors = FALSE
  )
}

eva_M <- process_eva(eva_M_raw, "Male")
eva_F <- process_eva(eva_F_raw, "Female")

# Compute log-scale columns (for use in calc_metrics)
add_log_cols <- function(df) {
  df$Log.OBS <- log10(df$OBS)
  df$Log.PRE <- log10(df$PRE)
  df$OPR     <- df$PRE / df$OBS              # predicted/observed ratio
  df$log.OPR <- df$Log.PRE - df$Log.OBS      # = log10(PRE/OBS)
  df
}
eva_M <- add_log_cols(eva_M)
eva_F <- add_log_cols(eva_F)

eva_all <- rbind(eva_M, eva_F)

# By sex overall
eva_gender <- eva_all %>%
  group_by(Gender) %>%
  group_modify(~calc_metrics(.x)) %>%
  ungroup()

write.csv(eva_gender, "Error_metrics_eva_gender.csv", row.names = FALSE)
