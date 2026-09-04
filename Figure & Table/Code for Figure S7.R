#!/usr/bin/env Rscript

# ==============================================================================
# Generate Figure S7: single-parameter swapping (GenX Human PBPK)
# ==============================================================================

suppressPackageStartupMessages({
  library(mrgsolve)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
})

# ------------------------------------------------------------------------------
# 1. Path setup
# ------------------------------------------------------------------------------

script_dir <- tryCatch({
  args_full <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("--file=", "", grep("--file=", args_full, value = TRUE))
  if (length(file_arg) > 0 && nzchar(file_arg)) {
    dirname(normalizePath(file_arg))
  } else {
    NA_character_
  }
}, error = function(e) NA_character_)

if (!is.na(script_dir)) {
  setwd(script_dir)
}

male_model_file   <- "GenX Hmod_M.R"
female_model_file <- "../Female/GenX Hmod_F.R"
outdir            <- "./GenX_parameter_swapping"
DOSE_A            <- 3.074e-06

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 2. Compile male and female models
# ------------------------------------------------------------------------------

env_m <- new.env(parent = globalenv())
sys.source(male_model_file, envir = env_m)
male_code <- get("HumanPBPK.code", envir = env_m)
mod_male <- mcode("HumanPBPK_Swap_M", male_code)

env_f <- new.env(parent = globalenv())
sys.source(female_model_file, envir = env_f)
female_code <- get("HumanPBPK.code", envir = env_f)
mod_female <- mcode("HumanPBPK_Swap_F", female_code)

# ------------------------------------------------------------------------------
# 3. Full nominal human parameter sets
# ------------------------------------------------------------------------------

male_par <- c(
  BW      = 73,
  QCC     = 15.62,
  QLC     = 0.25,
  QLuC    = 0.025,
  QKC     = 0.19,
  Htc     = 0.434,
  VPlasC  = 0.0411,
  VLC     = 0.0247,
  VLuC    = 0.0068,
  VKC     = 0.0042,
  VfilC   = 0.00042,
  GFRC    = 24.19,
  GEC     = 3.51,
  FVBK    = 0.160,
  Free    = 0.0106,
  PL      = 0.753,
  PLu     = 2.0399,
  PK      = 0.854,
  PRest   = 0.14,
  Kabsc   = 0.142,
  K0C     = 0.25,
  KunabsC = 0.00019,
  KurineC = 0.00854
)

female_par <- c(
  BW      = 60,
  QCC     = 16.42,
  QLC     = 0.27,
  QLuC    = 0.025,
  QKC     = 0.17,
  Htc     = 0.385,
  VPlasC  = 0.04,
  VLC     = 0.0233,
  VLuC    = 0.007,
  VKC     = 0.0046,
  VfilC   = 0.00046,
  GFRC    = 27.28,
  GEC     = 3.51,
  FVBK    = 0.160,
  Free    = 0.0704,
  PL      = 1.0649,
  PLu     = 0.331,
  PK      = 0.854,
  PRest   = 0.595,
  Kabsc   = 0.286,
  K0C     = 0.135,
  KunabsC = 0.0036,
  KurineC = 0.0165
)

# ------------------------------------------------------------------------------
# 4. Child phase
# ------------------------------------------------------------------------------

pred_child <- function(mod, pars, dose_mgkgday = DOSE_A) {

  tinterval <- 24
  TDoses <- 365 * 17

  ex <- tibble(
    ID = rep(1, TDoses - 365 + 1),
    time = seq(
      from = 24 * 365,
      to = tinterval * TDoses,
      by = tinterval
    )
  ) %>%
    mutate(
      DAY = time / 24,
      YEAR = DAY / 365,
      BW = if_else(
        YEAR <= 17,
        (9.86 + 0.370 * YEAR) /
          (1 - 0.0789 * YEAR + 0.00205 * YEAR^2),
        63
      ),
      amt = dose_mgkgday * BW,
      cmt = "AST",
      ii = tinterval,
      evid = 1
    )

  tsamp <- tgrid(
    start = 24 * 365,
    end = tinterval * TDoses,
    delta = tinterval
  )

  out <- mod %>%
    param(as.list(pars)) %>%
    update(atol = 1e-8, maxsteps = 5000) %>%
    mrgsim_d(data = ex, tgrid = tsamp) %>%
    as.data.frame() %>%
    filter(time > 0)

  out
}

# ------------------------------------------------------------------------------
# 5. Extract child final state for adult initialization
# ------------------------------------------------------------------------------

get_child_init <- function(child_out) {

  child_out %>%
    slice_tail(n = 1) %>%
    select(
      APlas_free,
      AFil,
      Aurine,
      ARest,
      AST,
      ASI,
      Afeces,
      AL,
      ALu,
      AKb,
      AT
    )
}

# ------------------------------------------------------------------------------
# 6. Adult phase: 33 years
# ------------------------------------------------------------------------------

pred_adult <- function(
    mod,
    pars,
    dose_mgkgday,
    Init,
    return_timecourse = FALSE
) {

  adult_BW <- unname(pars["BW"])

  tinterval <- 24
  GTDOSE <- 365 * 33
  GDOSEoral <- dose_mgkgday * adult_BW

  Gex.oral <- ev(
    ID = 1,
    time = 0,
    amt = GDOSEoral,
    ii = tinterval,
    addl = GTDOSE - 1,
    cmt = "AST",
    replicate = FALSE
  )

  Gtsamp <- tgrid(
    0,
    tinterval * (GTDOSE - 1) + 24 * 1,
    24
  )

  Gout <- mod %>%
    init(
      APlas_free = Init$APlas_free,
      AFil       = Init$AFil,
      Aurine     = Init$Aurine,
      ARest      = Init$ARest,
      AST        = Init$AST,
      ASI        = Init$ASI,
      Afeces     = Init$Afeces,
      AL         = Init$AL,
      ALu        = Init$ALu,
      AKb        = Init$AKb,
      ADOSE      = Init$AT + GDOSEoral
    ) %>%
    param(as.list(pars)) %>%
    update(atol = 1e-3, maxsteps = 500000) %>%
    mrgsim_d(data = Gex.oral, tgrid = Gtsamp) %>%
    as.data.frame()

  Goutdf <- Gout %>%
    transmute(
      Time = time / (24 * 365),
      CPlas = Plasma * 1000,
      CL = Liver * 1000,
      CK = Kidney * 1000,
      CLu = Lung * 1000
    )

  if (return_timecourse) {
    return(Goutdf)
  }

  Goutdf %>%
    filter(Time == 33)
}

# ------------------------------------------------------------------------------
# 7. Complete 17+33 scenario
# ------------------------------------------------------------------------------

run_full_scenario <- function(
    sex = c("Male", "Female"),
    pars,
    dose_mgkgday = DOSE_A,
    return_timecourse = FALSE
) {

  sex <- match.arg(sex)
  mod <- if (sex == "Male") mod_male else mod_female

  child_out <- pred_child(
    mod = mod,
    pars = pars,
    dose_mgkgday = dose_mgkgday
  )

  Init <- get_child_init(child_out)

  adult_out <- pred_adult(
    mod = mod,
    pars = pars,
    dose_mgkgday = dose_mgkgday,
    Init = Init,
    return_timecourse = return_timecourse
  )

  list(
    child_final_state = Init,
    adult = adult_out
  )
}

# ------------------------------------------------------------------------------
# 8. Baseline male/female simulations
# ------------------------------------------------------------------------------

baseline_male <- run_full_scenario(sex = "Male", pars = male_par)
baseline_female <- run_full_scenario(sex = "Female", pars = female_par)

base_m <- baseline_male$adult
base_f <- baseline_female$adult

baseline_ratios <- tibble(
  endpoint = c("Plasma", "Liver", "Kidney", "Lung"),
  male = c(base_m$CPlas, base_m$CL, base_m$CK, base_m$CLu),
  female = c(base_f$CPlas, base_f$CL, base_f$CK, base_f$CLu)
) %>%
  mutate(Male_Female_ratio = male / female)

baseline_plasma_ratio <- baseline_ratios$Male_Female_ratio[
  baseline_ratios$endpoint == "Plasma"
]

# ------------------------------------------------------------------------------
# 9. Helper for a single counterfactual swap
# ------------------------------------------------------------------------------

run_swap <- function(
    direction = c("Male_to_Female", "Female_to_Male"),
    changed_parameters,
    scenario_name
) {

  direction <- match.arg(direction)

  if (direction == "Male_to_Female") {
    scenario_sex <- "Male"
    base_par <- male_par
    donor_par <- female_par
  } else {
    scenario_sex <- "Female"
    base_par <- female_par
    donor_par <- male_par
  }

  swapped_par <- base_par
  swapped_par[changed_parameters] <- donor_par[changed_parameters]

  res <- run_full_scenario(
    sex = scenario_sex,
    pars = swapped_par
  )

  end <- res$adult

  tibble(
    direction = direction,
    scenario = scenario_name,
    changed_parameters = paste(changed_parameters, collapse = ";"),
    CPlas = end$CPlas,
    CL = end$CL,
    CK = end$CK,
    CLu = end$CLu
  )
}

# ------------------------------------------------------------------------------
# 10. Single-parameter swaps
# ------------------------------------------------------------------------------

single_parameters <- names(male_par)[male_par != female_par]

single_list <- list()
idx <- 1

for (p in single_parameters) {
  single_list[[idx]] <- run_swap(
    direction = "Male_to_Female",
    changed_parameters = p,
    scenario_name = p
  )
  idx <- idx + 1

  single_list[[idx]] <- run_swap(
    direction = "Female_to_Male",
    changed_parameters = p,
    scenario_name = p
  )
  idx <- idx + 1
}

single_summary <- bind_rows(single_list)

# Always orient ratios as Male/Female:
# M->F : counterfactual male / baseline female
# F->M : baseline male / counterfactual female
single_summary <- single_summary %>%
  mutate(
    Plasma_MF_ratio = case_when(
      direction == "Male_to_Female" ~ CPlas / base_f$CPlas,
      direction == "Female_to_Male" ~ base_m$CPlas / CPlas
    ),
    Liver_MF_ratio = case_when(
      direction == "Male_to_Female" ~ CL / base_f$CL,
      direction == "Female_to_Male" ~ base_m$CL / CL
    ),
    Delta_Plasma_ratio =
      Plasma_MF_ratio - baseline_ratios$Male_Female_ratio[
        baseline_ratios$endpoint == "Plasma"
      ],
    Abs_change_Plasma_ratio = abs(Delta_Plasma_ratio)
  )

# ------------------------------------------------------------------------------
# 11. Figure S7: single-parameter swapping plot
# ------------------------------------------------------------------------------

single_plot <- single_summary %>%
  select(direction, scenario, Plasma_MF_ratio, Liver_MF_ratio) %>%
  pivot_longer(
    cols = c(Plasma_MF_ratio, Liver_MF_ratio),
    names_to = "endpoint",
    values_to = "MF_ratio"
  ) %>%
  mutate(
    endpoint = recode(
      endpoint,
      Plasma_MF_ratio = "Plasma",
      Liver_MF_ratio = "Liver"
    )
  )

single_plot_plasma <- single_plot %>%
  filter(endpoint == "Plasma") %>%
  mutate(
    direction_label = case_when(
      direction == "Male_to_Female" ~ "M -> F",
      direction == "Female_to_Male" ~ "F -> M"
    )
  )

# Order by influence (weak -> strong, strongest at top)
single_order_abs <- single_summary %>%
  group_by(scenario) %>%
  summarise(max_abs = max(Abs_change_Plasma_ratio, na.rm = TRUE)) %>%
  arrange(max_abs) %>%
  pull(scenario)

single_plot_plasma$scenario <- factor(
  single_plot_plasma$scenario,
  levels = single_order_abs
)

# Annotate Top 3 influential parameters (top 3 in each direction)
top_labels_mf <- single_plot_plasma %>%
  filter(direction == "Male_to_Female") %>%
  mutate(abs_delta = abs(MF_ratio - baseline_plasma_ratio)) %>%
  arrange(desc(abs_delta)) %>%
  head(3)

top_labels_fm <- single_plot_plasma %>%
  filter(direction == "Female_to_Male") %>%
  mutate(abs_delta = abs(MF_ratio - baseline_plasma_ratio)) %>%
  arrange(desc(abs_delta)) %>%
  head(3)

# Top labels: union of top 3 from each direction, annotate both directions
top_scenarios <- union(top_labels_mf$scenario, top_labels_fm$scenario)
top_labels_all <- single_plot_plasma %>%
  filter(scenario %in% top_scenarios) %>%
  mutate(label_side = ifelse(MF_ratio < baseline_plasma_ratio, "left", "right"))

top_labels_left  <- top_labels_all %>% filter(label_side == "left")
top_labels_right <- top_labels_all %>% filter(label_side == "right")

p_single <- ggplot(
  single_plot_plasma,
  aes(
    x = MF_ratio,
    y = scenario,
    color = direction_label,
    shape = direction_label
  )
) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_vline(
    xintercept = baseline_plasma_ratio,
    linetype = "longdash", color = "black", linewidth = 0.6
  ) +
  geom_point(size = 2.6, position = position_dodge(width = 0.6)) +
  # Points left of the dashed line: label on the left side
  geom_text(
    data = top_labels_left,
    aes(
      label = sprintf("%.2f", MF_ratio),
      color = direction_label
    ),
    size = 2.6,
    hjust = 1.4,
    vjust = 0.45,
    family = "sans",
    position = position_dodge(width = 0.6),
    show.legend = FALSE
  ) +
  # Points right of the dashed line: label on the right side
  geom_text(
    data = top_labels_right,
    aes(
      label = sprintf("%.2f", MF_ratio),
      color = direction_label
    ),
    size = 2.6,
    hjust = -0.4,
    vjust = 0.45,
    family = "sans",
    position = position_dodge(width = 0.6),
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c("M -> F" = "#D55E00", "F -> M" = "#0072B2")
  ) +
  scale_shape_manual(
    values = c("M -> F" = 16, "F -> M" = 17)
  ) +
  # X-axis directly displays baseline value
  scale_x_log10(
    breaks = c(1, 2, 3, 5, baseline_plasma_ratio),
    labels = c("1", "2", "3", "5", sprintf("%.2f", baseline_plasma_ratio))
  ) +
  labs(
    x = "Male/Female plasma concentration ratio after swap",
    y = NULL,
    color = "Swap direction",
    shape = "Swap direction"
  ) +
  theme_classic(base_size = 9) +
  theme(
    text = element_text(family = "sans"),
    legend.position = c(0.02, 0.03),
    legend.justification = c(0, 0),
    legend.background = element_rect(
      color = "black", fill = "white", linewidth = 0.3
    ),
    legend.margin = margin(t = 4, r = 6, b = 4, l = 6, unit = "pt"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.32, "cm"),
    axis.ticks.length = unit(-0.15, "cm"),
    axis.text.x = element_text(margin = margin(t = 0.15, unit = "cm")),
    axis.text.y = element_text(margin = margin(r = 0.15, unit = "cm")),
    plot.background = element_rect(
      color = "black", fill = "white", linewidth = 0.8
    ),
    plot.margin = margin(t = 14, r = 36, b = 14, l = 36, unit = "pt")
  )

ggsave(
  file.path(outdir, "Figure_S7.jpg"),
  p_single,
  width = 7.48,
  height = 6.0,
  dpi = 300,
  quality = 95
)
