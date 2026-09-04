## ========================================================================
## Single oral dose: human male model prediction vs experimental observations
## Reference: Abraham et al. (2024) Environ Int. 193:109047
##   A healthy male volunteer received a single oral dose of an HFPO-DA
##   containing 15-PFAS mixture; plasma concentrations tracked for 450 days
## Body weight: 82 kg | single oral dose: 19.58 ug (total)
## ========================================================================

## ---- 1. Load libraries and model ----
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(showtext)
library(scales)

setwd("C:/Users/15960/Desktop/Modfit/Human/Male")
source(file = "GenX Hmod_M.R")
mod <- mcode("HumanPBPK.code", HumanPBPK.code)

## ---- 2. Parameter settings (override BW=82) ----
## Model parameters (from Fitcode final values), only BW modified
Human.theta <- log(c(
  BW      = 82,        # body weight 82 kg (overrides default 73)
  QCC     = 15.62,
  QLC     = 0.25,
  QLuC    = 0.025,
  QKC     = 0.19,
  Htc     = 0.434,
  VPlasC  = 0.0411,
  VLC     = 0.0247,
  VLuC    = 0.0068,
  VKC     = 0.0042,
  VFilC   = 0.00042,
  GEC     = 3.51,
  GFRC    = 24.19,
  FVBK    = 0.160,
  Free    = 0.0106,
  PL      = 0.753,
  PK      = 0.854,
  PLu     = 2.0399,
  PRest   = 0.14,
  K0C     = 0.25,
  Kabsc   = 0.142,
  KunabsC = 0.00019,
  KurineC = 0.00854
))

## ---- 3. Single oral dose simulation ----
## Dose: 19.58 ug = 0.01958 mg (total dose)
## DOSE (mg/kg) = 0.01958 / 82 = 2.3878e-4 mg/kg
BW_user   <- 82                       # kg
DOSE_ug   <- 19.58                    # micrograms (total dose)
DOSE_mg   <- DOSE_ug / 1000           # 0.01958 mg
DOSE_mgkg <- DOSE_mg / BW_user        # 2.3878e-4 mg/kg

## Single-dose regimen
tinterval <- 24                       # time interval (h), irrelevant for single dose
GBW       <- BW_user
GDOSE     <- DOSE_mgkg                # mg/kg
GDOSEoral <- GDOSE * GBW              # 0.01958 mg (total oral dose)

## Simulate 30 days = 720 hours (matches observed data range)
sim_days  <- 30
sim_hours <- sim_days * 24

Gex.oral <- ev(ID = 1,
               time = 0,             # dosing time
               amt  = GDOSEoral,     # total dose (mg)
               ii   = tinterval,
               addl = 0,             # single dose: no additional doses
               cmt  = "AST",
               replicate = FALSE)

## Sampling: every 6 hours (fine resolution in first 7 days, sufficient later)
tsamp <- tgrid(0, sim_hours, 6)

Gout <- mod %>%
  param(lapply(Human.theta, exp)) %>%
  update(atol = 1E-8, maxsteps = 500000) %>%
  mrgsim_d(data = Gex.oral, tgrid = tsamp)

## Extract concentrations (mg/L -> ng/mL = ug/L, multiply by 1000)
pred_df <- data.frame(
  Time_day = Gout$time / 24,
  Time_hr  = Gout$time,
  Plasma   = Gout$Plasma * 1000,     # ng/mL
  Liver    = Gout$Liver  * 1000      # ng/mL
)

## ---- 4. Experimental observation data (human measured values) ----
## Source: 人体GENX实验.csv
## time: day, concentration: ng/mL
obs_raw <- read.csv("人体GENX实验.csv")
obs_df <- data.frame(
  Time_day   = obs_raw$time,
  Plasma_obs = obs_raw$concentration
)

## ---- 5. Error metrics + half-life + Cmax ----
calc_metrics_single <- function(pred, obs) {
  ## pred: predicted value vector, obs: observed value vector (same length or interpolated by time)
  valid <- is.finite(pred) & is.finite(obs) & obs > 0 & pred > 0
  n <- sum(valid)
  if (n == 0) return(data.frame(N = 0, RMSE = NA, MAE = NA, AFE = NA, AAFE = NA, R2 = NA, pct_2fold = NA))
  p <- pred[valid]; o <- obs[valid]
  log_p <- log10(p); log_o <- log10(o)
  OPR <- p / o
  RMSE <- sqrt(mean((log_p - log_o)^2))
  MAE  <- mean(abs(log_p - log_o))
  AFE  <- 10^mean(log10(OPR))
  AAFE <- 10^mean(abs(log10(OPR)))
  R2 <- if (n >= 3) summary(lm(log_o ~ log_p))$adj.r.squared else NA
  pct_2fold <- mean(OPR >= 0.5 & OPR <= 2) * 100
  data.frame(N = n, RMSE = RMSE, MAE = MAE, AFE = AFE, AAFE = AAFE, R2 = R2, pct_2fold = pct_2fold)
}

## Interpolate predictions at observed time points
obs_df$Plasma_pred <- approx(pred_df$Time_day, pred_df$Plasma,
                              xout = obs_df$Time_day, rule = 2)$y

## Compute error metrics (log10 scale)
metrics <- calc_metrics_single(obs_df$Plasma_pred, obs_df$Plasma_obs)

## ---- 5a. Cmax and Tmax ----
Cmax_pred <- max(pred_df$Plasma)
Tmax_pred <- pred_df$Time_day[which.max(pred_df$Plasma)]
Cmax_obs  <- max(obs_df$Plasma_obs)
Tmax_obs  <- obs_df$Time_day[which.max(obs_df$Plasma_obs)]

metrics$Cmax_pred  <- Cmax_pred
metrics$Tmax_pred  <- Tmax_pred
metrics$Cmax_obs   <- Cmax_obs
metrics$Tmax_obs   <- Tmax_obs
metrics$Cmax_ratio <- Cmax_pred / Cmax_obs

## ---- 5b. Half-life (t1/2) ----
## Based on elimination phase (terminal log-linear regression): log10(C) = a + b*t
## Slope b = -Kel/ln(10), so t1/2 = log10(2) / |b| = 0.301/|b|
calc_halflife <- function(time, conc, t_range) {
  idx <- time >= t_range[1] & time <= t_range[2] & conc > 0 & is.finite(conc)
  if (sum(idx) < 3) return(NA)
  fit <- lm(log10(conc[idx]) ~ time[idx])
  b <- coef(fit)[2]           # log10 regression slope (negative)
  if (b >= 0) return(NA)
  t_half <- -log10(2) / b     # = 0.301 / |b|
  return(unname(t_half))
}

# Predicted half-life: 5-30 day elimination phase
t1_2_pred <- calc_halflife(pred_df$Time_day, pred_df$Plasma, c(5, 30))
# Observed half-life: 5-28 day elimination phase
t1_2_obs  <- calc_halflife(obs_df$Time_day, obs_df$Plasma_obs, c(5, 28))

metrics$t1_2_pred <- t1_2_pred
metrics$t1_2_obs  <- t1_2_obs

write.csv(metrics, "Single_dose_error_metrics.csv", row.names = FALSE)

## ---- 6. Figure S5: prediction vs observation ----
showtext_auto()

p_plot <- ggplot() +
  ## PBPK prediction curve
  geom_line(data = pred_df, aes(x = Time_day, y = Plasma, color = "Predicted"),
            linewidth = 1.2) +
  ## Observed points (scatter)
  geom_point(data = obs_df, aes(x = Time_day, y = Plasma_obs, color = "Observed"),
             size = 3, shape = 16, alpha = 0.85) +
  scale_color_manual(name = NULL,
                     values = c("Predicted" = "#2166AC", "Observed" = "#B2182B"),
                     breaks = c("Observed", "Predicted")) +
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 30),
                labels = c("0.01", "0.1", "1", "10", "30"),
                expand = c(0.02, 0.02)) +
  scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1, 2),
                labels = c("0.001", "0.01", "0.1", "1", "2"),
                expand = c(0.02, 0.02)) +
  coord_cartesian(xlim = c(0.008, 35), ylim = c(0.0008, 2.5)) +
  labs(x = "Time (days)",
       y = "Plasma concentration (ng/mL)") +
  theme_classic(base_size = 36) +
  theme(
    text = element_text(family = "sans", face = "bold"),
    axis.line   = element_blank(),
    axis.ticks  = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(-0.12, "cm"),
    axis.text.x = element_text(margin = margin(t = 0.18, unit = "cm"), color = "black",
                               size = 28, face = "bold"),
    axis.text.y = element_text(margin = margin(r = 0.18, unit = "cm"), color = "black",
                               size = 28, face = "bold"),
    axis.title  = element_text(face = "bold", size = 30),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.margin = margin(8, 8, 8, 8),
    legend.position = c(0.95, 0.95),
    legend.justification = c(1, 1),
    legend.background = element_rect(color = "black", fill = "white",
                                     linewidth = 0.4, linetype = "solid"),
    legend.text = element_text(size = 22, face = "bold"),
    legend.title = element_text(size = 22, face = "bold"),
    legend.key.size = unit(0.5, "cm"),
    legend.margin = margin(2, 4, 2, 4),
    legend.box.margin = margin(1, 1, 1, 1)
  )

ggsave("Figure_S5.png", p_plot, width = 7.48, height = 5.5, dpi = 300)
