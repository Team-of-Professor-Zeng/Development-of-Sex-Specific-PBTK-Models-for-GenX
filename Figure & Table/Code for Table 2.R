##############################  HED Derivation  ##############################################
### Male
## Human
library(mrgsolve)    # Needed for Loading mrgsolve code into r via mcode from the 'mrgsolve' pckage
library(magrittr)    # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)       # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(ggplot2)     # Needed for plot
library(FME)         # Package for MCMC simulation and model fitting
library(minpack.lm)  # Package for model fitting
library(reshape)     # Package for melt function to reshape the table
library(truncnorm)   # Package for the truncated normal distribution function   
library(EnvStats)    # Package for Environmental Statistics, Including US EPA Guidance
library(invgamma)    # Package for inverse gamma distribution function
library(foreach)     # Package for parallel computing
library(doParallel)  # Package for parallel computing
library(bayesplot)   # Package for MCMC traceplot
library(gridExtra)   # Package for combine ggplot
library(ggExtra)     # R-package for extend the ggplot2
library(reshape2)    # melt function to reshape the table
library(tidyverse)   # Needed for the pipe %>% operator
library(coda)        # for convergence diagnosis
library(DescTools)   # descriptive statistics and exploratory data analysis 
library(scales)      # for plotting the figure
library(grid)        # for plotting the figure
library(lattice)     # for plotting the figure

setwd("C:/Users/15960/Desktop/2025GenX/Modfit/Human/Male")
source (file = "GenX Hmod_M.R")
mod <- mcode ("HumanPBPK.code", HumanPBPK.code)
theta.int.M <- log(c(
  KurineC                        = 0.00854,                        
  Free                           = 0.0106,                        
  PL                             = 0.753,                        
  PK                             = 0.854,                         
  PLu                            = 2.0399,
  PRest                          = 0.14,                         
  K0C                            = 0.136,                           
  Kabsc                          = 0.288,                      
  KunabsC                        = 0.00019                     
))

Pred.child <- function(pars,DOSE_A) {
  
  ## Get out of log domain
  pars <- exp(pars)
  
  tinterval <- 24  
  TDoses    <- 365 * 17 
  ex <- tibble(
    ID   = rep(1, TDoses - 365 + 1), 
    time = seq(from = 24 * 365, to = tinterval * TDoses, by = tinterval)  
  ) %>%
    mutate(
      DAY   = time / 24,  
      YEAR  = DAY / 365, 
      BW    = if_else(YEAR <= 17, 
                      true  = (9.86 + 0.370*YEAR)/(1 - 0.0789*YEAR + 0.00205*YEAR^2),
                      ifelse(YEAR > 17, 63, 63)),           ##Equation from Chang 2022
      amt   = DOSE_A * BW,  
      cmt   = "AST",        
      ii    = tinterval,    
      evid  = 1)
  
  tsamp <- tgrid(start = 24 * 365, end = tinterval * TDoses, delta = tinterval)
  
  out <- mod %>%
    param (pars) %>%
    update(atol = 1E-8, maxsteps = 5000) %>%
    mrgsim_d(data = ex, tgrid = tsamp)%>%
    filter (time > 0)
  
  return(out)
}

HED_M <- function(pars, DOSE_A, N) {
  
  Init <- Pred.child (theta.int.M,DOSE_A) %>% filter (row_number()== n()) %>% select(-c("Plasma", "Liver"))
  
  ## Get out of log domain
  pars <- exp(pars)
  
  ## Time parameters
  tinterval <- 24  
  TDoses    <- 365 * 30 
  
  ## Generate individual-specific parameters
  idata <- tibble(ID = 1:N) %>% 
    mutate(
      BW        = rnormTrunc  (N, min = 30.0768, max = 115.923, mean = 73,      sd = 21.9),      #CV=0.3
      VKC       = rnormTrunc  (N, min = 0.00173, max = 0.00667, mean = 0.0042,  sd = 0.00126),   #CV=0.3
      GFRC      = rnormTrunc  (N, min = 9.967  , max =  38.413, mean = 24.19 ,  sd = 7.257) ,    #CV=0.3
      DOSEoral.A = DOSE_A * BW
    )
  
  ## Define time grid for simulation
  tsamp <- tgrid(start = 24 * 365, end = tinterval * TDoses, delta = 24 * 365)
  
  ## Define dosing events using idata
  Gex.oral_1 <- ev(
    ID   = 1:N,             # One individual
    time = 0,               # Dose start time
    amt  = idata$DOSEoral.A, # Individual-specific dose
    ii   = tinterval,       # Time interval
    addl = TDoses - 1,      # Additional dosing
    cmt  = "AST",           # The dosing compartment: AST Stomach  
    replicate = FALSE       # No replicate
  )
  
  ## Run the simulation
  out <- mod %>%
    init(APlas_free = Init$APlas_free, AFil = Init$AFil, AKb = Init$AKb, ARest = Init$ARest,AL = Init$AL, ALu = Init$ALu, 
         Aurine = Init$Aurine, AST = Init$AST, ASI= Init$ASI, Afeces= Init$Afeces) %>% 
    param(pars) %>%
    data_set(Gex.oral_1) %>%
    idata_set(idata) %>%
    update(atol = 1E-8, maxsteps = 5000) %>%
    mrgsim(obsonly = TRUE, tgrid = tsamp)
  
  ## Format output
  out <- cbind.data.frame(
    ID        = out$ID,
    Time      = out$time / (24 * 365),  # Convert time to years
    CPlas     = out$Plasma ,          
    CL        = out$Liver,
    AUC_CPlas = out$AUC_CA,
    AUC_CL    = out$AUC_CL)   
  
  out <- out %>% filter(Time == 30) }

N = 1000
DOSE_B = 0.06 #mg/kg/day， Thompson.et al.2026
HAUC.M    <- HED_M (theta.int.M,DOSE_B,N = N) %>% select(c("ID","Time","CPlas","CL","AUC_CPlas","AUC_CL")) #CPlas, CL	mg/L; AUC_CPlas, AUC_CL	mg·h/L

HAUC.M$Avg.CPlas = HAUC.M$AUC_CPlas/(30*365*24) #mg/L
HAUC.M$Avg.CL    = HAUC.M$AUC_CL   /(30*365*24) #mg/L
write.csv(HAUC.M, file = "HAUC.M3.csv", row.names = FALSE)

## Mice
setwd("C:/Users/15960/Desktop/2025GenX/Modfit/Mice/Male")
source (file = "GenX MMod_M.R")
mod <- mcode ("micepbpk", MicePBPK_M.code)

mice.theta.int.M <- log(c(
  KurineC                        = 0.0628,                        
  Free                           = 0.0106,                 ### fitting parameters
  PL                             = 0.753,                  ### fitting parameters
  PK                             = 0.854,
  PLu                            = 2.0399,                 ### fitting parameters
  PRest                          = 0.14 ,                  ### fitting parameters
  K0C                            = 1.838,                  ### fitting parameters               
  Kabsc                          = 1.0461,                 ### fitting parameters                       
  KunabsC                        = 0.00140                   
))
mice.HED_M <- function(pars, DOSE_A, N) {
  
  ## Get out of log domain
  pars <- exp(pars)
  
  ## Time parameters
  tinterval <- 24  
  TDoses    <- 540 
  
  ## Generate individual-specific parameters
  idata <- tibble(ID = 1:N) %>% 
    mutate(
      BW        = rnormTrunc  (N, min = 0.0103, max = 0.0397, mean = 0.025,      sd = 0.0075),   #CV=0.3
      VKC       = rnormTrunc  (N, min = 0.007 , max = 0.0270, mean = 0.017,      sd = 0.0051),   #CV=0.3
      GFRC      = rnormTrunc  (N, min = 25.586, max = 98.614, mean = 62.1 ,      sd = 18.63) ,   #CV=0.3
      DOSEoral.A = DOSE_A * BW
    )
  
  ## Define time grid for simulation
  tsamp <- tgrid(start = 0, end = tinterval * TDoses, delta = 1)
  
  ## Define dosing events using idata
  Gex.oral_1 <- ev(
    ID   = 1:N,             # One individual
    time = 0,               # Dose start time
    amt  = idata$DOSEoral.A, # Individual-specific dose
    ii   = tinterval,       # Time interval
    addl = TDoses - 1,      # Additional dosing
    cmt  = "AST",           # The dosing compartment: AST Stomach  
    replicate = FALSE       # No replicate
  )
  
  ## Run the simulation
  out <- mod %>%
    param(pars) %>%
    data_set(Gex.oral_1) %>%
    idata_set(idata) %>%
    update(atol = 1E-8, maxsteps = 5000) %>%
    mrgsim(obsonly = TRUE, tgrid = tsamp)
  
  ## Format output
  out <- cbind.data.frame(
    ID        = out$ID,
    Time      = out$time,
    CPlas     = out$Plasma ,          
    CL        = out$Liver,
    AUC_CPlas = out$AUC_CA,
    AUC_CL    = out$AUC_CL)   
  
  out <- out %>% filter(Time == 540) }

N = 1000
DOSE_B = 0.06 #mg/kg/day，Thompson.et al.2026
MAUC.M    <- mice.HED_M (mice.theta.int.M,DOSE_B,N = N) %>% select(c("ID","Time","CPlas","CL","AUC_CPlas","AUC_CL"))

MAUC.M$Avg.CPlas = MAUC.M$AUC_CPlas/(540*24) 
MAUC.M$Avg.CL    = MAUC.M$AUC_CL   /(540*24)
write.csv(MAUC.M, file = "MAUC.M.csv", row.names = FALSE)

##########################################################################################################################
### Female
## Human
setwd("C:/Users/15960/Desktop/2025GenX/Modfit/Human/Female")
source (file = "GenX Hmod_F.R")
mod <- mcode ("HumanPBPK.code", HumanPBPK.code)
theta.int.F <- log(c(
  KurineC                        = 0.0166,                        
  Free                           = 0.0704,                        
  PL                             = 1.0649,                        
  PK                             = 0.854,                         
  PLu                            = 0.331,
  PRest                          = 0.595,                         
  K0C                            = 0.136,                           
  Kabsc                          = 0.288,                      
  KunabsC                        = 0.0036                     
))

Pred.child <- function(pars,DOSE_A) {
  
  ## Get out of log domain
  pars <- exp(pars)
  
  tinterval <- 24  
  TDoses    <- 365 * 17 
  ex <- tibble(
    ID   = rep(1, TDoses - 365 + 1), 
    time = seq(from = 24 * 365, to = tinterval * TDoses, by = tinterval)  
  ) %>%
    mutate(
      DAY   = time / 24,  
      YEAR  = DAY / 365, 
      BW    = if_else(YEAR <= 17, 
                      true  = (9.86 + 0.370*YEAR)/(1 - 0.0789*YEAR + 0.00205*YEAR^2),
                      ifelse(YEAR > 17, 63, 63)),           ##Equation from Chang 2022
      amt   = DOSE_A * BW,  
      cmt   = "AST",        
      ii    = tinterval,    
      evid  = 1)
  
  tsamp <- tgrid(start = 24 * 365, end = tinterval * TDoses, delta = tinterval)
  
  out <- mod %>%
    param (pars) %>%
    update(atol = 1E-8, maxsteps = 5000) %>%
    mrgsim_d(data = ex, tgrid = tsamp)%>%
    filter (time > 0)
  
  return(out)
}

HED_F <- function(pars, DOSE_A, N) {
  
  Init <- Pred.child (theta.int.F,DOSE_A) %>% filter (row_number()== n()) %>% select(-c("Plasma", "Liver"))
  
  ## Get out of log domain
  pars <- exp(pars)
  
  ## Time parameters
  tinterval <- 24  
  TDoses    <- 365 * 30 
  
  ## Generate individual-specific parameters
  idata <- tibble(ID = 1:N) %>% 
    mutate(
      BW        = rnormTrunc  (N, min = 24.721, max = 95.279, mean = 60,      sd = 18),       #CV=0.3
      VKC       = rnormTrunc  (N, min = 0.0019, max = 0.0073, mean = 0.0046,  sd = 0.00138),  #CV=0.3
      GFRC      = rnormTrunc  (N, min = 11.24 , max = 43.32 , mean = 27.28 ,  sd = 8.184) ,   #CV=0.3
      DOSEoral.A = DOSE_A * BW
    )
  
  ## Define time grid for simulation
  tsamp <- tgrid(start = 24 * 365, end = tinterval * TDoses, delta = 24 * 365)
  
  ## Define dosing events using idata
  Gex.oral_1 <- ev(
    ID   = 1:N,             # One individual
    time = 0,               # Dose start time
    amt  = idata$DOSEoral.A, # Individual-specific dose
    ii   = tinterval,       # Time interval
    addl = TDoses - 1,      # Additional dosing
    cmt  = "AST",           # The dosing compartment: AST Stomach  
    replicate = FALSE       # No replicate
  )
  
  ## Run the simulation
  out <- mod %>%
    init(APlas_free = Init$APlas_free, AFil = Init$AFil, AKb = Init$AKb, ARest = Init$ARest,AL = Init$AL, ALu = Init$ALu, 
         Aurine = Init$Aurine, AST = Init$AST, ASI= Init$ASI, Afeces= Init$Afeces) %>% 
    param(pars) %>%
    data_set(Gex.oral_1) %>%
    idata_set(idata) %>%
    update(atol = 1E-8, maxsteps = 5000) %>%
    mrgsim(obsonly = TRUE, tgrid = tsamp)
  
  ## Format output
  out <- cbind.data.frame(
    ID        = out$ID,
    Time      = out$time / (24 * 365),  # Convert time to years
    CPlas     = out$Plasma ,          
    CL        = out$Liver,
    AUC_CPlas = out$AUC_CA,
    AUC_CL    = out$AUC_CL)   
  
  out <- out %>% filter(Time == 30) }

N = 1000
DOSE_B = 0.19 #mg/kg/day，Thompson.et al.2026
HAUC.F    <- HED_F (theta.int.F,DOSE_B,N = N) %>% select(c("ID","Time","CPlas","CL","AUC_CPlas","AUC_CL"))

HAUC.F$Avg.CPlas = HAUC.F$AUC_CPlas/(30*365*24) 
HAUC.F$Avg.CL    = HAUC.F$AUC_CL   /(30*365*24)
write.csv(HAUC.F, file = "HAUC.F.csv", row.names = FALSE)

## Mice
setwd("C:/Users/15960/Desktop/2025GenX/Modfit/Mice/Female")
source (file = "GenX MMod_F.R")
mod <- mcode ("micepbpk", MicePBPK_F.code)

mice.theta.int.F <- log(c(
  KurineC                        = 0.122,                        
  Free                           = 0.0704,                ### fitting parameters
  PL                             = 1.0649,                 ### fitting parameters
  PK                             = 0.854,
  PLu                            = 0.331,                 ### fitting parameters
  PRest                          = 0.595,                  
  K0C                            = 1,                
  Kabsc                          = 2.12,                         
  KunabsC                        = 0.0265                   
))
mice.HED_F <- function(pars, DOSE_A, N) {
  
  ## Get out of log domain
  pars <- exp(pars)
  
  ## Time parameters
  tinterval <- 24  
  TDoses    <- 540 
  
  ## Generate individual-specific parameters
  idata <- tibble(ID = 1:N) %>% 
    mutate(
      BW        = rnormTrunc  (N, min = 0.0082, max = 0.0318, mean = 0.02 ,      sd = 0.006 ),   #CV=0.3
      VKC       = rnormTrunc  (N, min = 0.007 , max = 0.0270, mean = 0.017,      sd = 0.0051),   #CV=0.3
      GFRC      = rnormTrunc  (N, min = 16.913, max = 65.167, mean = 41.04,      sd = 12.31 ) ,  #CV=0.3
      DOSEoral.A = DOSE_A * BW
    )
  
  ## Define time grid for simulation
  tsamp <- tgrid(start = 0, end = tinterval * TDoses, delta = 1)
  
  ## Define dosing events using idata
  Gex.oral_1 <- ev(
    ID   = 1:N,             # One individual
    time = 0,               # Dose start time
    amt  = idata$DOSEoral.A, # Individual-specific dose
    ii   = tinterval,       # Time interval
    addl = TDoses - 1,      # Additional dosing
    cmt  = "AST",           # The dosing compartment: AST Stomach  
    replicate = FALSE       # No replicate
  )
  
  ## Run the simulation
  out <- mod %>%
    param(pars) %>%
    data_set(Gex.oral_1) %>%
    idata_set(idata) %>%
    update(atol = 1E-8, maxsteps = 5000) %>%
    mrgsim(obsonly = TRUE, tgrid = tsamp)
  
  ## Format output
  out <- cbind.data.frame(
    ID        = out$ID,
    Time      = out$time,
    CPlas     = out$Plasma ,          
    CL        = out$Liver,
    AUC_CPlas = out$AUC_CA,
    AUC_CL    = out$AUC_CL)   
  
  out <- out %>% filter(Time == 540) }

N = 1000
DOSE_B = 0.19 #mg/kg/day，Thompson.et al.2026
MAUC.F    <- mice.HED_F (mice.theta.int.F,DOSE_B,N = N) %>% select(c("ID","Time","CPlas","CL","AUC_CPlas","AUC_CL"))
MAUC.F$Avg.CPlas = MAUC.F$AUC_CPlas/(540*24) 
MAUC.F$Avg.CL    = MAUC.F$AUC_CL   /(540*24)
write.csv(MAUC.F, file = "MAUC.F.csv", row.names = FALSE)

################################################################################################################################
## Calculated the HED 
## This study method
## Male
HAUC.M <- read.csv("C:/Users/15960/Desktop/2025GenX/Modfit/Human/Male/HAUC.M.csv", check.names = FALSE)
MAUC.M <- read.csv("C:/Users/15960/Desktop/2025GenX/Modfit/Mice/Male/MAUC.M.csv", check.names = FALSE)

AUC.M <- data.frame(
  HED.CA.M = ((MAUC.M$Avg.CPlas)/(HAUC.M$Avg.CPlas)) * 0.06,
  HED.CL.M = ((MAUC.M$Avg.CL   )/(HAUC.M$Avg.CL   ))* 0.06,
  ASC = (MAUC.M$Avg.CPlas)/(HAUC.M$Avg.CPlas),
  ALC = (MAUC.M$Avg.CL   )/(HAUC.M$Avg.CL   ))

## Estimated the median (95% CI) of ASC in mice
M.ASC.M.range <- cbind.data.frame (
  Median.CA    = quantile(MAUC.M$Avg.CPlas, probs = 0.5, names = FALSE,na.rm=T),
  lower.CA     = quantile(MAUC.M$Avg.CPlas, probs = 0.025, names = FALSE,na.rm=T),
  upper.CA     = quantile(MAUC.M$Avg.CPlas, probs = 0.975, names = FALSE,na.rm=T),
  Median.CL    = quantile(MAUC.M$Avg.CL   , probs = 0.5, names = FALSE,na.rm=T),
  lower.CL     = quantile(MAUC.M$Avg.CL   , probs = 0.025, names = FALSE,na.rm=T),
  upper.CL     = quantile(MAUC.M$Avg.CL   , probs = 0.975, names = FALSE,na.rm=T)
)
M.ASC.M.range <- signif(M.ASC.M.range, 3)
write.csv(M.ASC.M.range, "C:/Users/15960/Desktop/2025GenX/Modfit/Human/M.ASC.M.range.csv", row.names = FALSE)
CA_string <- sprintf("%.2e (%.2e, %.2e)", M.ASC.M.range$Median.CA, M.ASC.M.range$lower.CA, M.ASC.M.range$upper.CA)
CL_string <- sprintf("%.2e (%.2e, %.2e)", M.ASC.M.range$Median.CL, M.ASC.M.range$lower.CL, M.ASC.M.range$upper.CL)
cat("CA:", CA_string, "\n")
cat("CL:", CL_string, "\n")

## Estimated the median (95% CI) of ASC mice/ASC human
ASC_ASC.M.range <- cbind.data.frame (
  Median.CA    = quantile(AUC.M$ASC, probs = 0.5, names = FALSE,na.rm=T),
  lower.CA     = quantile(AUC.M$ASC, probs = 0.025, names = FALSE,na.rm=T),
  upper.CA     = quantile(AUC.M$ASC, probs = 0.975, names = FALSE,na.rm=T),
  Median.CL    = quantile(AUC.M$ALC, probs = 0.5, names = FALSE,na.rm=T),
  lower.CL     = quantile(AUC.M$ALC, probs = 0.025, names = FALSE,na.rm=T),
  upper.CL     = quantile(AUC.M$ALC, probs = 0.975, names = FALSE,na.rm=T)
)
ASC_ASC.M.range <- signif(ASC_ASC.M.range, 3)
CA_string <- sprintf("%.3g (%.3g, %.3g)", ASC_ASC.M.range$Median.CA, ASC_ASC.M.range$lower.CA, ASC_ASC.M.range$upper.CA)
CL_string <- sprintf("%.3g (%.3g, %.3g)", ASC_ASC.M.range$Median.CL, ASC_ASC.M.range$lower.CL, ASC_ASC.M.range$upper.CL)
cat("CA:", CA_string, "\n")
cat("CL:", CL_string, "\n")

## Estimated HED based on the AUC Method
HED.human.M.range <- cbind.data.frame (
  Median.CA  = quantile(AUC.M$HED.CA.M , probs = 0.5, names = FALSE,na.rm=T),
  lower.CA   = quantile(AUC.M$HED.CA.M , probs = 0.025, names = FALSE,na.rm=T),
  upper.CA   = quantile(AUC.M$HED.CA.M , probs = 0.975, names = FALSE,na.rm=T),
  Median.CL  = quantile(AUC.M$HED.CL.M , probs = 0.5, names = FALSE,na.rm=T),
  lower.CL   = quantile(AUC.M$HED.CL.M , probs = 0.025, names = FALSE,na.rm=T),
  upper.CL   = quantile(AUC.M$HED.CL.M , probs = 0.975, names = FALSE,na.rm=T)
)
HED.human.M.range <- signif(HED.human.M.range, 3)
write.csv(HED.human.M.range, "C:/Users/15960/Desktop/2025GenX/Modfit/Human/HED.human.M.range.csv", row.names = FALSE)
CA_string <- sprintf("%.3g (%.3g, %.3g)", HED.human.M.range$Median.CA, HED.human.M.range$lower.CA, HED.human.M.range$upper.CA)
CL_string <- sprintf("%.3g (%.3g, %.3g)", HED.human.M.range$Median.CL, HED.human.M.range$lower.CL, HED.human.M.range$upper.CL)
cat("CA:", CA_string, "\n")
cat("CL:", CL_string, "\n")


## Female
HAUC.F <- read.csv("C:/Users/15960/Desktop/2025GenX/Modfit/Human/Female/HAUC.F.csv", check.names = FALSE)
MAUC.F <- read.csv("C:/Users/15960/Desktop/2025GenX/Modfit/Mice/Female/MAUC.F.csv", check.names = FALSE)
AUC.F <- data.frame(
  HED.CA.F = ((MAUC.F$Avg.CPlas)/(HAUC.F$Avg.CPlas))*0.19,
  HED.CL.F = ((MAUC.F$Avg.CL   )/(HAUC.F$Avg.CL   ))*0.19,
  ASC      = (MAUC.F$Avg.CPlas)/(HAUC.F$Avg.CPlas),
  ALC      = (MAUC.F$Avg.CL   )/(HAUC.F$Avg.CL   ))

## Estimated the median (95% CI) of ASC
M.ASC.F.range <- cbind.data.frame (
  Median.CA    = quantile(MAUC.F$Avg.CPlas, probs = 0.5, names = FALSE,na.rm=T),
  lower.CA     = quantile(MAUC.F$Avg.CPlas, probs = 0.025, names = FALSE,na.rm=T),
  upper.CA     = quantile(MAUC.F$Avg.CPlas, probs = 0.975, names = FALSE,na.rm=T),
  Median.CL    = quantile(MAUC.F$Avg.CL   , probs = 0.5, names = FALSE,na.rm=T),
  lower.CL     = quantile(MAUC.F$Avg.CL   , probs = 0.025, names = FALSE,na.rm=T),
  upper.CL     = quantile(MAUC.F$Avg.CL   , probs = 0.975, names = FALSE,na.rm=T)
)
M.ASC.F.range <- signif(M.ASC.F.range, 3)
write.csv(M.ASC.F.range, "C:/Users/15960/Desktop/2025GenX/Modfit/Human/M.ASC.F.range.csv", row.names = FALSE)
CA_string <- sprintf("%.3g (%.3g, %.3g)", M.ASC.F.range$Median.CA, M.ASC.F.range$lower.CA, M.ASC.F.range$upper.CA)
CL_string <- sprintf("%.3g (%.3g, %.3g)", M.ASC.F.range$Median.CL, M.ASC.F.range$lower.CL, M.ASC.F.range$upper.CL)
cat("CA:", CA_string, "\n")
cat("CL:", CL_string, "\n")

# Estimated the median (95% CI) of ASC mice/ASC human
ASC_ASC.F.range <- cbind.data.frame (
  Median.CA    = quantile(AUC.F$ASC, probs = 0.5, names = FALSE,na.rm=T),
  lower.CA     = quantile(AUC.F$ASC, probs = 0.025, names = FALSE,na.rm=T),
  upper.CA     = quantile(AUC.F$ASC, probs = 0.975, names = FALSE,na.rm=T),
  Median.CL    = quantile(AUC.F$ALC, probs = 0.5, names = FALSE,na.rm=T),
  lower.CL     = quantile(AUC.F$ALC, probs = 0.025, names = FALSE,na.rm=T),
  upper.CL     = quantile(AUC.F$ALC, probs = 0.975, names = FALSE,na.rm=T)
)
ASC_ASC.F.range <- signif(ASC_ASC.F.range, 3)
CA_string <- sprintf("%.3g (%.3g, %.3g)", ASC_ASC.F.range$Median.CA, ASC_ASC.F.range$lower.CA, ASC_ASC.F.range$upper.CA)
CL_string <- sprintf("%.3g (%.3g, %.3g)", ASC_ASC.F.range$Median.CL, ASC_ASC.F.range$lower.CL, ASC_ASC.F.range$upper.CL)
cat("CA:", CA_string, "\n")
cat("CL:", CL_string, "\n")

## Estimated HED based on the AUC Method
HED.human.F.range <- cbind.data.frame (
  Median.CA  = quantile(AUC.F$HED.CA.F , probs = 0.5, names = FALSE,na.rm=T),
  lower.CA   = quantile(AUC.F$HED.CA.F , probs = 0.025, names = FALSE,na.rm=T),
  upper.CA   = quantile(AUC.F$HED.CA.F , probs = 0.975, names = FALSE,na.rm=T),
  Median.CL  = quantile(AUC.F$HED.CL.F , probs = 0.5, names = FALSE,na.rm=T),
  lower.CL   = quantile(AUC.F$HED.CL.F , probs = 0.025, names = FALSE,na.rm=T),
  upper.CL   = quantile(AUC.F$HED.CL.F , probs = 0.975, names = FALSE,na.rm=T)
)
HED.human.F.range <- signif(HED.human.F.range, 3)
write.csv(HED.human.F.range, "C:/Users/15960/Desktop/2025GenX/Modfit/Human/HED.human.F.range.csv", row.names = FALSE)
CA_string <- sprintf("%.3g (%.3g, %.3g)", HED.human.F.range$Median.CA, HED.human.F.range$lower.CA, HED.human.F.range$upper.CA)
CL_string <- sprintf("%.3g (%.3g, %.3g)", HED.human.F.range$Median.CL, HED.human.F.range$lower.CL, HED.human.F.range$upper.CL)
cat("CA:", CA_string, "\n")
cat("CL:", CL_string, "\n")





