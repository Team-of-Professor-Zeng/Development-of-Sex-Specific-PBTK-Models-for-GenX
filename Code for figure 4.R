#Predicted value

### Male ####################################################################################
setwd("D:/zs/PBPK/1-方案/GenX/小鼠/code/Human/Male")
source (file = "GenX Hmod_M.R")

## Build mrgsolve-based PBPK Model
mod <- mcode ("HumanPBPK.code", HumanPBPK.code)

theta.int.H <- log(c(
  KurineC                        = 0.00854,                        
  Free                           = 0.0106,                        
  PL                             = 0.753,                        
  PK                             = 0.854,                         
  PLu                            = 2.0399,
  PRest                          = 0.14,                         
  K0C                            = 0.25,                           
  Kabsc                          = 0.142,                      
  KunabsC                        = 0.00019                     
))

#saveRDS(theta.int, file = "int_H.rds") 

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

## The assumed exposure dose were estiamted from previous literatures
# Init value: Chinese population     ;  Dose: 3.074 ng/kg/day estimated from Dong et al. (2024);
DOSE_A <- 3.074E-06
Init <- Pred.child(theta.int.H,DOSE_A) %>% filter (row_number()== n()) %>% select(-c("Plasma", "Liver", "Kidney", "Lung"))

## Exposure sceinario 
pred.adult <- function(pars, DOSE,Init) {
  
  ## Get out of log domain
  pars <- lapply(pars, exp)## Return a list of exp (parametrs for gestational model) from log scale
  
  ## Exposure scenario for gestational exposure
  GBW          = 73                 ## Body weight during gestation 
  tinterval    = 24                 ## Time interval; 
  GTDOSE       = 365 * 33           ## Total dosing/Dose times; 
  GDOSE        = DOSE               ## Input oral dose  
  GDOSEoral    = GDOSE*GBW          ## Amount of oral dose
  
  # To create exposure scenario
  Gex.oral <- ev (ID   = 1, 
                  time = 0,             ## Dosing start time
                  amt  = GDOSEoral,     ## Amount of dose 
                  ii   = tinterval,     ## Time interval
                  addl = GTDOSE - 1,    ## Additional doseing 
                  cmt  = "AST",         ## dosing: AST Stomach  
                  replicate = FALSE)    ## No replicate
  
  Gtsamp  = tgrid(0, tinterval*(GTDOSE - 1) + 24*1, 24) ## Simulation time 
  
  ## Simulation of exposure scenario
  Gout <- 
    mod %>%
    init(APlas_free = Init$APlas_free,AFil = Init$AFil, Aurine = Init$Aurine, ARest= Init$ARest, AST = Init$AST, ASI= Init$ASI, Afeces= Init$Afeces, 
         AL = Init$AL, ALu = Init$ALu, AKb = Init$AKb, ADOSE = Init$AT+GDOSEoral) %>% 
    param (pars) %>%
    update(atol = 1E-3, maxsteps = 500000) %>%          
    mrgsim_d (data = Gex.oral, tgrid = Gtsamp)
  
  Goutdf = cbind.data.frame(Time   = Gout$time/(24 * 365), 
                            CPlas  = Gout$Plasma*1000,
                            CL     = Gout$Liver*1000,
                            CK     = Gout$Kidney*1000,
                            CLu    = Gout$Lung*1000)
  
  Goutdf <- Goutdf %>% filter (Time == 33) 
  return (Goutdf)
}

## Estimate the model residual with experimental data by modCost function (from FME package)
outdf.A <- pred.adult (pars = theta.int.H, DOSE = DOSE_A, Init = Init)

## MCMC
PBPK_H_pop <- function(pars, DOSE_A, N) {
  
  Init <- Pred.child (theta.int.H,DOSE_A) %>% filter (row_number()== n()) %>% select(-c("Plasma", "Liver", "Kidney", "Lung"))
  
  ## Get out of log domain
  pars <- exp(pars)
  
  ## Time parameters
  tinterval <- 24  
  TDoses    <- 365 * 33 
  
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
    init(APlas_free = Init$APlas_free,AFil = Init$AFil, Aurine = Init$Aurine, ARest= Init$ARest, AST = Init$AST, ASI= Init$ASI, Afeces= Init$Afeces, 
         AL = Init$AL, ALu = Init$ALu, AKb = Init$AKb) %>% 
    param(pars) %>%
    data_set(Gex.oral_1) %>%
    idata_set(idata) %>%
    update(atol = 1E-8, maxsteps = 5000) %>%
    mrgsim(obsonly = TRUE, tgrid = tsamp)
  
  ## Format output
  out <- cbind.data.frame(
    ID    = out$ID,
    Time  = out$time / (24 * 365),  # Convert time to years
    CPlas = out$Plasma * 1000,          # Convert concentration to ng/mL
    CL    = out$Liver* 1000,
    CK    = out$Kidney* 1000,
    CLu   = out$Lung* 1000)         
  
  out <- out %>% filter(Time == 33) 
  
}

N = 3000
PlotDat_A1   <- PBPK_H_pop (theta.int.H,DOSE_A,N = N)
write.csv(PlotDat_A1, file = "PlotDat_A1_M3000.csv", row.names = FALSE)
H_CPlas <- PlotDat_A1  %>% summarize (min    = min(CPlas, na.rm = TRUE),
                                      P2.5   = quantile (CPlas, probs = 0.025),
                                      P25    = quantile (CPlas, probs = 0.25),
                                      median = quantile (CPlas, probs = 0.5), 
                                      P75    = quantile (CPlas, probs = 0.75),
                                      P97.5  = quantile (CPlas, probs = 0.975),
                                      max    = max(CPlas, na.rm = TRUE),
                                      mean   = mean(CPlas, na.rm = TRUE),
                                      sd     = sd(CPlas, na.rm = TRUE))
H_CPlas

### Female ##################################################################################################################################################
setwd("D:/zs/PBPK/1-方案/GenX/小鼠/code/Human/Female")
source (file = "GenX Hmod_F.R")

## Build mrgsolve-based PBPK Model
mod <- mcode ("HumanPBPK.code", HumanPBPK.code)

theta.int.H <- log(c(
  KurineC                        = 0.0165,                        
  Free                           = 0.0704,                        
  PL                             = 1.0649,                        
  PK                             = 0.854,                         
  PLu                            = 0.331,
  PRest                          = 0.595,                         
  K0C                            = 0.135,                           
  Kabsc                          = 0.286,                      
  KunabsC                        = 0.0036                     
))

#saveRDS(theta.int, file = "int_H.rds") 

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

## The assumed exposure dose were estiamted from previous literatures
# Init value: Chinese population     ;  Dose: 3.074 ng/kg/day estimated from Dong et al. (2024);
DOSE_A <- 3.074E-06
Init <- Pred.child(theta.int.H,DOSE_A) %>% filter (row_number()== n()) %>% select(-c("Plasma", "Liver", "Kidney", "Lung"))

## Exposure sceinario 
pred.adult <- function(pars, DOSE,Init) {
  
  ## Get out of log domain
  pars <- lapply(pars, exp)## Return a list of exp (parametrs for gestational model) from log scale
  
  ## Exposure scenario for gestational exposure
  GBW          = 60                 ## Body weight during gestation 
  tinterval    = 24                 ## Time interval; 
  GTDOSE       = 365 * 33           ## Total dosing/Dose times; 
  GDOSE        = DOSE               ## Input oral dose  
  GDOSEoral    = GDOSE*GBW          ## Amount of oral dose
  
  # To create exposure scenario
  Gex.oral <- ev (ID   = 1, 
                  time = 0,             ## Dosing start time
                  amt  = GDOSEoral,     ## Amount of dose 
                  ii   = tinterval,     ## Time interval
                  addl = GTDOSE - 1,    ## Additional doseing 
                  cmt  = "AST",         ## dosing: AST Stomach  
                  replicate = FALSE)    ## No replicate
  
  Gtsamp  = tgrid(0, tinterval*(GTDOSE - 1) + 24*1, 24) ## Simulation time 
  
  ## Simulation of exposure scenario
  Gout <- 
    mod %>%
    init(APlas_free = Init$APlas_free,AFil = Init$AFil, Aurine = Init$Aurine, ARest= Init$ARest, AST = Init$AST, ASI= Init$ASI, Afeces= Init$Afeces, 
         AL = Init$AL, ALu = Init$ALu, AKb = Init$AKb, ADOSE = Init$AT+GDOSEoral) %>% 
    param (pars) %>%
    update(atol = 1E-3, maxsteps = 500000) %>%          
    mrgsim_d (data = Gex.oral, tgrid = Gtsamp)
  
  Goutdf = cbind.data.frame(Time   = Gout$time/(24 * 365), 
                            CPlas  = Gout$Plasma*1000,
                            CL     = Gout$Liver*1000,
                            CK     = Gout$Kidney*1000,
                            CLu    = Gout$Lung*1000)
  
  Goutdf <- Goutdf %>% filter (Time == 33) 
  return (Goutdf)
}

## Estimate the model residual with experimental data by modCost function (from FME package)
outdf.A <- pred.adult (pars = theta.int.H, DOSE = DOSE_A, Init = Init)

## MCMC
PBPK_H_pop <- function(pars, DOSE_A, N) {
  
  Init <- Pred.child (theta.int.H,DOSE_A) %>% filter (row_number()== n()) %>% select(-c("Plasma", "Liver", "Kidney", "Lung"))
  
  ## Get out of log domain
  pars <- exp(pars)
  
  ## Time parameters
  tinterval <- 24  
  TDoses    <- 365 * 33 
  
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
    init(APlas_free = Init$APlas_free,AFil = Init$AFil, Aurine = Init$Aurine, ARest= Init$ARest, AST = Init$AST, ASI= Init$ASI, Afeces= Init$Afeces, 
         AL = Init$AL, ALu = Init$ALu, AKb = Init$AKb) %>% 
    param(pars) %>%
    data_set(Gex.oral_1) %>%
    idata_set(idata) %>%
    update(atol = 1E-8, maxsteps = 5000) %>%
    mrgsim(obsonly = TRUE, tgrid = tsamp)
  
  ## Format output
  out <- cbind.data.frame(
    ID    = out$ID,
    Time  = out$time / (24 * 365),  # Convert time to years
    CPlas = out$Plasma * 1000,          # Convert concentration to ng/mL
    CL    = out$Liver* 1000,
    CK    = out$Kidney* 1000,
    CLu   = out$Lung* 1000)         
  
  out <- out %>% filter(Time == 33) 
  
}

N = 4920
PlotDat_A2   <- PBPK_H_pop (theta.int.H,DOSE_A,N = N)
write.csv(PlotDat_A2, file = "PlotDat_A2_F4920.csv", row.names = FALSE)
H_CPlas <- PlotDat_A2  %>% summarize (min    = min(CPlas, na.rm = TRUE),
                                      P2.5   = quantile (CPlas, probs = 0.025),
                                      P25    = quantile (CPlas, probs = 0.25),
                                      median = quantile (CPlas, probs = 0.5), 
                                      P75    = quantile (CPlas, probs = 0.75),
                                      P97.5  = quantile (CPlas, probs = 0.975),
                                      max    = max(CPlas, na.rm = TRUE),
                                      mean   = mean(CPlas, na.rm = TRUE),
                                      sd     = sd(CPlas, na.rm = TRUE))
H_CPlas


## Combined ##########################################################################################################################
setwd("D:/zs/PBPK/1-方案/GenX/小鼠/code/Human/Male")
M508 <- read.csv(file ="PlotDat_A1_M508.csv")
setwd("D:/zs/PBPK/1-方案/GenX/小鼠/code/Human/Female")
F492 <- read.csv(file ="PlotDat_A2_F492.csv")
PlotDat1   <- rbind(M508, F492)
Combined1 <- PlotDat1  %>% summarize (min    = min(CPlas, na.rm = TRUE),
                                      P2.5   = quantile (CPlas, probs = 0.025),
                                      P25    = quantile (CPlas, probs = 0.25),
                                      median = quantile (CPlas, probs = 0.5), 
                                      P75    = quantile (CPlas, probs = 0.75),
                                      P97.5  = quantile (CPlas, probs = 0.975),
                                      max    = max(CPlas, na.rm = TRUE),
                                      mean   = mean(CPlas, na.rm = TRUE),
                                      sd     = sd(CPlas, na.rm = TRUE))
Combined1

setwd("D:/zs/PBPK/1-方案/GenX/小鼠/code/Human/Male")
M300 <- read.csv(file ="PlotDat_A1_M300.csv")
setwd("D:/zs/PBPK/1-方案/GenX/小鼠/code/Human/Female")
F700 <- read.csv(file ="PlotDat_A2_F700.csv")
#write.csv(F700, file = "PlotDat_A2_F700.csv", row.names = FALSE)
PlotDat2   <- rbind(M300, F700)
Combined2 <- PlotDat2  %>% summarize (min    = min(CPlas, na.rm = TRUE),
                                      P2.5   = quantile (CPlas, probs = 0.025),
                                      P25    = quantile (CPlas, probs = 0.25),
                                      median = quantile (CPlas, probs = 0.5), 
                                      P75    = quantile (CPlas, probs = 0.75),
                                      P97.5  = quantile (CPlas, probs = 0.975),
                                      max    = max(CPlas, na.rm = TRUE),
                                      mean   = mean(CPlas, na.rm = TRUE),
                                      sd     = sd(CPlas, na.rm = TRUE))
Combined2

##############################################################################################################
##Male
setwd("D:/zs/PBPK/1-方案/GenX/小鼠/code/Human/Male")
source (file = "GenX Hmod_M.R")

## Build mrgsolve-based PBPK Model
mod <- mcode ("HumanPBPK.code", HumanPBPK.code)

theta.int.H <- log(c(
  KurineC                        = 0.00854,                        
  Free                           = 0.0106,                        
  PL                             = 0.753,                        
  PK                             = 0.854,                         
  PLu                            = 2.0399,
  PRest                          = 0.14,                         
  K0C                            = 0.25,                           
  Kabsc                          = 0.142,                      
  KunabsC                        = 0.00019                     
))

#saveRDS(theta.int, file = "int_H.rds") 

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

## The assumed exposure dose were estiamted from previous literatures
# Init value: Chinese population     ;  Dose: 3.074 ng/kg/day estimated from Dong et al. (2024);
DOSE_A <- 3.074E-06
Init <- Pred.child(theta.int.H,DOSE_A) %>% filter (row_number()== n()) %>% select(-c("Plasma", "Liver", "Kidney", "Lung"))

## Exposure sceinario 
pred.adult <- function(pars, DOSE,Init) {
  
  ## Get out of log domain
  pars <- lapply(pars, exp)## Return a list of exp (parametrs for gestational model) from log scale
  
  ## Exposure scenario for gestational exposure
  GBW          = 73                 ## Body weight during gestation 
  tinterval    = 24                 ## Time interval; 
  GTDOSE       = 365 * 33           ## Total dosing/Dose times; 
  GDOSE        = DOSE               ## Input oral dose  
  GDOSEoral    = GDOSE*GBW          ## Amount of oral dose
  
  # To create exposure scenario
  Gex.oral <- ev (ID   = 1, 
                  time = 0,             ## Dosing start time
                  amt  = GDOSEoral,     ## Amount of dose 
                  ii   = tinterval,     ## Time interval
                  addl = GTDOSE - 1,    ## Additional doseing 
                  cmt  = "AST",         ## dosing: AST Stomach  
                  replicate = FALSE)    ## No replicate
  
  Gtsamp  = tgrid(0, tinterval*(GTDOSE - 1) + 24*1, 24) ## Simulation time 
  
  ## Simulation of exposure scenario
  Gout <- 
    mod %>%
    init(APlas_free = Init$APlas_free,AFil = Init$AFil, Aurine = Init$Aurine, ARest= Init$ARest, AST = Init$AST, ASI= Init$ASI, Afeces= Init$Afeces, 
         AL = Init$AL, ALu = Init$ALu, AKb = Init$AKb, ADOSE = Init$AT+GDOSEoral) %>% 
    param (pars) %>%
    update(atol = 1E-3, maxsteps = 500000) %>%          
    mrgsim_d (data = Gex.oral, tgrid = Gtsamp)
  
  Goutdf = cbind.data.frame(Time   = Gout$time/(24 * 365), 
                            CPlas  = Gout$Plasma*1000,
                            CL     = Gout$Liver*1000,
                            CK     = Gout$Kidney*1000,
                            CLu    = Gout$Lung*1000)
  
  Goutdf <- Goutdf %>% filter (Time == 33) 
  return (Goutdf)
}

## Estimate the model residual with experimental data by modCost function (from FME package)
outdf.A <- pred.adult (pars = theta.int.H, DOSE = DOSE_A, Init = Init)

## MCMC
PBPK_H_annual <- function(pars, DOSE_A, N, years = 55) {
  
  Init <- Pred.child (theta.int.H,DOSE_A) %>% filter (row_number()== n()) %>% select(-c("Plasma"))
  
  ## Get out of log domain
  pars <- exp(pars)
  
  ## Time parameters
  tinterval <- 24  
  TDoses    <- 365 * years 
  
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
    init(APlas_free = Init$APlas_free,AFil = Init$AFil, Aurine = Init$Aurine, ARest= Init$ARest, AST = Init$AST, ASI= Init$ASI, Afeces= Init$Afeces, 
         AL = Init$AL, ALu = Init$ALu, AKb = Init$AKb) %>% 
    param(pars) %>%
    data_set(Gex.oral_1) %>%
    idata_set(idata) %>%
    update(atol = 1E-8, maxsteps = 5000) %>%
    mrgsim(obsonly = TRUE, tgrid = tsamp)
  
  ## Format output
  out_df <- data.frame(
    ID = out$ID,
    Time = out$time / (24 * 365),  # Convert time to years
    CPlas = out$Plasma * 1000      # Convert concentration to ng/mL
  )        
  return(out_df)
}

N = 1000
Annual_M <- PBPK_H_annual(theta.int.H, DOSE_A, N = N, years = 55)
Median_M <- Annual_M %>% group_by(Time) %>% summarize (Median = quantile (CPlas, probs = 0.5))
Mean_M <- Annual_M %>% group_by(Time) %>% summarize (Mean   = mean(CPlas, na.rm = TRUE))

write.csv(Annual_M, "Annual_M.csv", row.names = FALSE)
write.csv(Median_M, "Median_M.csv", row.names = FALSE)
write.csv(Mean_M  , "Mean_M.csv", row.names = FALSE)
cat("Male Annual Median Concentration Summary:\n")
print(summary(Median_M$Median))

##Female
setwd("D:/zs/PBPK/1-方案/GenX/小鼠/code/Human/Female")
source (file = "GenX Hmod_F.R")
mod <- mcode ("HumanPBPK.code", HumanPBPK.code)

theta.int.H <- log(c(
  KurineC                        = 0.0165,                        
  Free                           = 0.0704,                        
  PL                             = 1.0649,                        
  PK                             = 0.854,                         
  PLu                            = 0.331,
  PRest                          = 0.595,                         
  K0C                            = 0.135,                           
  Kabsc                          = 0.286,                      
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

## The assumed exposure dose were estiamted from previous literatures
# Init value: Chinese population     ;  Dose: 3.074 ng/kg/day estimated from Dong et al. (2024);
DOSE_A <- 3.074E-06
Init <- Pred.child(theta.int.H,DOSE_A) %>% filter (row_number()== n()) %>% select(-c("Plasma"))

## Exposure sceinario 
pred.adult <- function(pars, DOSE,Init) {
  
  ## Get out of log domain
  pars <- lapply(pars, exp)## Return a list of exp (parametrs for gestational model) from log scale
  
  ## Exposure scenario for gestational exposure
  GBW          = 60                 ## Body weight during gestation 
  tinterval    = 24                 ## Time interval; 
  GTDOSE       = 365 * 33           ## Total dosing/Dose times; 
  GDOSE        = DOSE               ## Input oral dose  
  GDOSEoral    = GDOSE*GBW          ## Amount of oral dose
  
  # To create exposure scenario
  Gex.oral <- ev (ID   = 1, 
                  time = 0,             ## Dosing start time
                  amt  = GDOSEoral,     ## Amount of dose 
                  ii   = tinterval,     ## Time interval
                  addl = GTDOSE - 1,    ## Additional doseing 
                  cmt  = "AST",         ## dosing: AST Stomach  
                  replicate = FALSE)    ## No replicate
  
  Gtsamp  = tgrid(0, tinterval*(GTDOSE - 1) + 24*1, 24) ## Simulation time 
  
  ## Simulation of exposure scenario
  Gout <- 
    mod %>%
    init(APlas_free = Init$APlas_free,AFil = Init$AFil, Aurine = Init$Aurine, ARest= Init$ARest, AST = Init$AST, ASI= Init$ASI, Afeces= Init$Afeces, 
         AL = Init$AL, ALu = Init$ALu, AKb = Init$AKb, ADOSE = Init$AT+GDOSEoral) %>% 
    param (pars) %>%
    update(atol = 1E-3, maxsteps = 500000) %>%          
    mrgsim_d (data = Gex.oral, tgrid = Gtsamp)
  
  Goutdf = cbind.data.frame(Time   = Gout$time/(24 * 365), 
                            CPlas  = Gout$Plasma*1000,
                            CL     = Gout$Liver*1000,
                            CK     = Gout$Kidney*1000,
                            CLu    = Gout$Lung*1000)
  
  Goutdf <- Goutdf %>% filter (Time == 33) 
  return (Goutdf)
}

## Estimate the model residual with experimental data by modCost function (from FME package)
outdf.A <- pred.adult (pars = theta.int.H, DOSE = DOSE_A, Init = Init)

PBPK_F_annual <- function(pars, DOSE_A, N, years = 55) {
  Init <- Pred.child (theta.int.H,DOSE_A) %>% filter (row_number()== n()) %>% select(-c("Plasma"))
  ## Get out of log domain
  pars <- exp(pars)
  
  ## Time parameters
  tinterval   = 24  
  TDoses      = 365 * years 
  
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
    init(APlas_free = Init$APlas_free,AFil = Init$AFil, Aurine = Init$Aurine, ARest= Init$ARest, AST = Init$AST, ASI= Init$ASI, Afeces= Init$Afeces, 
         AL = Init$AL, ALu = Init$ALu, AKb = Init$AKb) %>% 
    param(pars) %>%
    data_set(Gex.oral_1) %>%
    idata_set(idata) %>%
    update(atol = 1E-8, maxsteps = 5000) %>%
    mrgsim(obsonly = TRUE, tgrid = tsamp)
  
  ## Format output
  out_df <- data.frame(
    ID = out$ID,
    Time = out$time / (24 * 365),  # Convert time to years
    CPlas = out$Plasma * 1000      # Convert concentration to ng/mL
  )        
  return(out_df)
}

N <- 1000
Annual_F <- PBPK_F_annual(theta.int.H, DOSE_A, N = N, years = 55)
Median_F <- Annual_F %>% group_by(Time) %>% summarize (Median = quantile (CPlas, probs = 0.5))
Mean_F   <- Annual_F %>% group_by(Time) %>% summarize (Mean   = mean(CPlas, na.rm = TRUE))

write.csv(Annual_F, "Annual_F.csv", row.names = FALSE)
write.csv(Median_F, "Median_F.csv", row.names = FALSE)
write.csv(Mean_F, "Mean_F.csv", row.names = FALSE)
cat("Female Annual Median Concentration Summary:\n")
print(summary(Median_F$Median))
Mean_M   <- read.csv("D:/zs/PBPK/1-方案/GenX/小鼠/code/Human/Male/Mean_M.csv")
Mean_F   <- read.csv("D:/zs/PBPK/1-方案/GenX/小鼠/code/Human/Female/Mean_F.csv")

### Plot figure 
#################################################################################################################################
### Figure S5
Median_M <- read.csv("D:/zs/PBPK/1-方案/GenX/小鼠/code/Human/Male/Median_M.csv")%>% filter(Time>=18)
Median_F <- read.csv("D:/zs/PBPK/1-方案/GenX/小鼠/code/Human/Female/Median_F.csv")%>% filter(Time>=18)

Median_M$Gender <- "Male"
Median_F$Gender <- "Female"
combined_data1 <- rbind(Median_M, Median_F)
combined_data1$Gender <- factor(combined_data1$Gender, levels = c("Male", "Female"))

p_male <- ggplot(combined_data1 %>% filter(Gender == "Male"), 
                 aes(x = Time, y = Median)) +
  geom_smooth(method = "loess", span = 0.3, se = TRUE, alpha = 0.3, 
              linewidth = 1.8, aes(color = "Male", fill = "Male")) +  
  scale_color_manual(values = c("Male" = "steelblue"), name = "Sex") +  
  scale_fill_manual(values = c("Male" = "lightblue"), name = "Sex") +
  labs(tag = "(A)", x = "Time (years)", y = "Predicted GenX Concentration in Human Plasma (ng/mL)") +
  theme_bw(base_family = "Times New Roman") +
  theme(
    panel.background = element_rect(fill = "White"),
    plot.background = element_rect(fill = "White"),
    text = element_text(family = "Times New Roman", size = 12),
    plot.tag = element_text(size = 30, face = "bold", family = "Times New Roman"),
    plot.tag.position = c(0.03, 1),
    plot.margin = margin(t = 20, r = 30, b = 10, l = 30, unit = "pt"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
    axis.text.y = element_text(size = 20, face = "bold", colour = "black", margin = margin(r = 8, unit = "pt")),
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    legend.position = c(0.8, 0.1),
    legend.title = element_text(size = 28, face = "bold", family = "Times New Roman"),  
    legend.text = element_text(size = 25, face = "bold", family = "Times New Roman"),
    axis.ticks.length.y = unit(0.2, "cm") 
  ) +
  scale_x_continuous(breaks = seq(0, 50, by = 10), expand = c(0.05, 0)) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE, nsmall = 11), expand = c(0, 0))+
  annotate("segment", x = -Inf, xend = Inf, y = Inf, yend = Inf, color = "black", linewidth = 1.2) +  
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf, color = "black", linewidth = 1.2) + 
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = Inf, color = "black", linewidth = 1.2)
print(p_male)

p_female <- ggplot(combined_data1 %>% filter(Gender == "Female"), 
                   aes(x = Time, y = Median)) +
  geom_smooth(method = "loess", span = 0.3, se = TRUE, alpha = 0.3, 
              linewidth = 1.8, aes(color = "Female", fill = "Female")) +
  scale_color_manual(values = c("Female" = "#E9687A"), name = NULL) +  
  scale_fill_manual(values = c("Female" = "#F6B3AC"), name = NULL) + 
  labs(x = "Time (years)", y = "Predicted GenX Concentration in Human Plasma (ng/mL)") +
  theme_bw(base_family = "Times New Roman") +
  theme(
    panel.background = element_rect(fill = "White"),
    plot.background = element_rect(fill = "White"),
    text = element_text(family = "Times New Roman", size = 12),
    plot.tag = element_text(size = 30, face = "bold", family = "Times New Roman"),
    plot.tag.position = c(0.03, 0.98),
    plot.margin = margin(t = 10, r = 30, b = 10, l = 30, unit = "pt"),
    axis.text.x = element_text(size = 20, face = "bold", colour = "black", margin = margin(t = 8, unit = "pt")),
    axis.text.y = element_text(size = 20, face = "bold", colour = "black", margin = margin(r = 8, unit = "pt")),
    axis.title.x = element_text(size = 22, face = "bold", margin = margin(t = 15, unit = "pt")),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    legend.position = c(0.81, 0.99),
    legend.title = element_blank(),
    legend.text = element_text(size = 25, face = "bold", family = "Times New Roman"), 
    axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
    axis.ticks.length.y = unit(0.2, "cm") 
  ) +
  scale_x_continuous(breaks = seq(0, 50, by = 10), expand = c(0.05, 0)) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE, nsmall = 11), expand = c(0.05, 0), breaks = c(0.01350326, 0.01350327, 0.01350328))+
  coord_cartesian(ylim = c(0.01350326, 0.01350328)) +
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf, color = "black", linewidth = 1.2) +  
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = 0.01350328, color = "black", linewidth = 1.2) + 
  annotate("segment", x = Inf, xend = Inf, y = -Inf, yend = 0.01350328, color = "black", linewidth = 1.2)
print(p_female)
p_combined <-grid.arrange(p_male, p_female, nrow = 2)

### ORININ ###############################################################################################################################
p1 <- ggplot   (combined_data1, aes(x = Time, y = Median, color = Gender, fill = Gender)) +
  geom_smooth(method = "loess", span = 0.3, se = TRUE, alpha = 0.3, linewidth = 1.8) +
  scale_color_manual(values = c("Male" = "steelblue", "Female" = "#E9687A")) +
  scale_fill_manual (values = c("Male" = "lightblue", "Female" = "#F6B3AC")) +
  labs(tag = "(A)", x = "Time (years)", y ="Predicted GenX Concentration in Human Plasma (ng/mL)", color = "Gender", fill = "Gender") +
  theme_bw(base_family = "Times New Roman")+
  theme(
    panel.background = element_rect(fill = "White"),
    plot.background = element_rect(fill = "White"),
    text = element_text(family = "Times New Roman", size = 12),
    plot.tag = element_text(size = 30, face = "bold", family = "Times New Roman"),
    plot.tag.position = c(0.03, 0.98),
    plot.margin = margin(t = 20, r = 5, b = 5, l = 30, unit = "pt"),
    axis.text.x = element_text(size = 17, face = "bold", colour = "black", margin = margin(t = 8, unit = "pt")),
    axis.text.y = element_text(size = 17, face = "bold", colour = "black", margin = margin(r = 8, unit = "pt")),
    axis.title.x = element_text(size = 22, face = "bold", margin = margin(t = 15, unit = "pt")),
    axis.title.y = element_text(size = 22, face = "bold", margin = margin(r = 15, unit = "pt")),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "none",  
    strip.text = element_text(size = 22, face = "bold", family = "Times New Roman"), 
    strip.background = element_blank(),
    panel.spacing = unit(0.6, "cm")) +
  scale_x_continuous(breaks = seq(0, 50, by = 10), expand = c(0.05, 0)) +
  scale_y_continuous(labels =  scales::scientific_format()) +
  facet_wrap(~ Gender, scales = "free_y")
print(p1)

########################################################################################################################################
### Figure 4
install.packages("ggtext")  
library(ggtext)
Human_Evaluation <- read.csv("D:/zs/PBPK/1-方案/GenX/小鼠/code/Human/Male/Human Evaluation.csv")
scatter <- Human_Evaluation[1:7,c(1, 3, 4, 5, 8, 9, 10)]
rows_a <- c(1, 2, 5, 7)
scatter[rows_a, 1] <- paste0(scatter[rows_a, 1], "<sup>a</sup>")
rows_b <- c(3, 4, 6)
scatter[rows_b, 1] <- paste0(scatter[rows_b, 1], "<sup>b</sup>")
print(scatter[, 1])
scatter$References <- factor(scatter$References, levels = scatter$References)
scatter$legend_group <- paste(scatter$Type, scatter$color, scatter$shape, sep = "|")

custom_labels <- as.character(levels(scatter$References))
custom_labels <- gsub("Prediction1<sup>a</sup>", "Prediction<sup>a</sup>", custom_labels)
custom_labels <- gsub("Prediction2<sup>a</sup>", "Prediction<sup>a</sup>", custom_labels)
custom_labels <- gsub("Prediction3<sup>a</sup>", "Prediction<sup>a</sup>", custom_labels)

p_scatter <- ggplot(scatter, aes(x = References, y = Median)) + 
  geom_point(aes(color = legend_group, shape = legend_group), size = 10) +
  geom_errorbar(aes(ymin = Median - P25, ymax = Median + P75, color = legend_group), width = 0.15, size = 1.6) +
  labs(tag = "(B)", x = "References", y = "GenX Concentration in Human Plasma (ng/mL)") +
  ylim(-0.01, 0.2) +
  scale_color_manual(name = "Type",
    values = c(
      "Female|a|e" = "#E9687A",
      "Female-Prediction|b|f" = "#F8B2A2", 
      "General Population|c|e" = "#8768A6",
      "General Population|c|e" = "#8768A6", 
      "General-Prediction|d|f" = "gray50",#5A5A5A
      "General Population|c|e" = "#8768A6", 
      "General-Prediction|d|f" = "gray50"),
    labels = c(
      "Female|a|e" = "Female",
      "Female-Prediction|b|f" = "Female-Prediction",
      "General Population|c|e" = "General Population",
      "General-Prediction|d|f" = "General-Prediction") ) +
  scale_shape_manual(name = "Type",
    values = c(
      "Female|a|e" = 16,
      "Female-Prediction|b|f" = 18,
      "General Population|c|e" = 16,
      "General-Prediction|d|f" = 18),
    labels = c(
      "Female|a|e" = "Female",
      "Female-Prediction|b|f" = "Female-Prediction",
      "General Population|c|e" = "General Population",
      "General-Prediction|d|f" = "General-Prediction" )) +
  scale_x_discrete(labels = custom_labels) +
  #scale_x_discrete(labels = function(x) {ifelse(grepl("Prediction", x), "Prediction", x)}) +
  guides(
    color = guide_legend(override.aes = list(size = 7))) +#,shape = guide_legend(override.aes = list(size = 7))
  theme_bw(base_family = "Times New Roman")+
  theme(
    plot.background = element_rect(fill = "White"),
    text = element_text(family = "Times New Roman"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.background = element_rect(fill = "White"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_markdown(size = 20, face = "bold", colour = "black", margin = margin(t = 8, unit = "pt"),
                               angle = 25, hjust = 0.9, vjust = 0.9),
    axis.text.y = element_text(size = 20, face = "bold", colour = "black", margin = margin(r = 8, unit = "pt")),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_line(linewidth = 1.2, color = "black"),
    axis.ticks.length.y = unit(0.2, "cm"),
    plot.tag = element_text(size = 30, face = "bold", family = "Times New Roman"),
    plot.tag.position = c(0.01, 1.05),
    plot.margin = margin(t = 10, r = 30, b = 10, l = 60, unit = "pt"),
    legend.position = c(0.02, 0.98),  
    legend.justification = c(0, 1),   
    legend.direction = "vertical",     
    legend.box = "vertical",  #horizontal
    legend.title = element_text(size = 28, face = "bold"),#, margin = margin(-10, 20, 0, 0)
    legend.text = element_text(size = 25, margin = margin(-2, 20, 0, 5)),
    legend.background =  element_blank())
print(p_scatter)

combined_plot <- grid.arrange(p_combined, p_scatter, nrow = 2, heights = c(1, 1), 
                              left = textGrob("GenX Concentration in Human Plasma (ng/mL)", x = unit(0.9, "npc"),
                                              gp = gpar(fontsize = 30, fontfamily = "Times New Roman", fontface = "bold"), rot = 90)  )
ggsave("Figure 4.tiff",scale = 1,
       plot = combined_plot,
       path = "D:/zs/PBPK/2025GenX/Artwork",
       width = 40, height = 40, units = "cm",dpi=320)

