## Load libraries
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

## Input mrgsolve-based PBPK Model
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

## NSC
Human.theta <- log(c(
  # Physiological parameters
  BW                             = 60,
  QCC                            = 16.42,
  QLC                            = 0.27,
  QLuC                           = 0.025,
  QKC                            = 0.17,
  Htc                            = 0.385,
  VPlasC                         = 0.04,
  VLC                            = 0.0233,
  VLuC                           = 0.007,
  VKC                            = 0.0046,
  VFilC                          = 0.00046,
  GEC                            = 3.51,
  GFRC                           = 27.28,
  FVBK                           = 0.160,
  
  # Chemical-specific parameters (final mean values)                       
  Free                           = 0.0704,                        
  PL                             = 1.0649,                        
  PK                             = 0.854,                         
  PLu                            = 0.331,
  PRest                          = 0.595,                         
  K0C                            = 0.135,                           
  Kabsc                          = 0.286,                      
  KunabsC                        = 0.0036,
  KurineC                        = 0.0165      
))

pred.Human <- function(pars, DOSE) {
  
  ## Get out of log domain
  pars <- lapply(pars, exp)## Return a list of exp (parametrs for gestational model) from log scale
  
  ## Exposure scenario for gestational exposure
  GBW          = 60                 ## Body weight during gestation 
  tinterval    = 24                 ## Time interval; 
  GTDOSE       = 365 * 30           ## Total dosing/Dose times; 
  GDOSE        = DOSE               ## Input oral dose  
  GDOSEoral    = GDOSE*GBW          ## Amount of oral dose
  
  # To create exposure scenario
  Gex.oral <- ev (ID   = rep(1, 24*365*30+1), 
                  time = 0,             ## Dosing start time
                  amt  = GDOSEoral,     ## Amount of dose 
                  ii   = tinterval,     ## Time interval
                  addl = GTDOSE - 1,    ## Additional doseing 
                  cmt  = "AST",         ## dosing: AST Stomach  
                  replicate = FALSE)    ## No replicate
  
  Gtsamp  = tgrid(0, tinterval*(GTDOSE - 1) + 24*1, 1) ## Simulation time 
  
  ## Simulation of exposure scenario
  Gout <- 
    mod %>%
    param (pars) %>%
    update(atol = 1E-3, maxsteps = 500000) %>%          
    mrgsim_d (data = Gex.oral, tgrid = Gtsamp)
  
  Goutdf = cbind.data.frame(Time      = Gout$time/(24 * 365), 
                            CPlas     = Gout$Plasma,
                            AUC_CPlas = Gout$AUC_CA,
                            AUC_CL    = Gout$AUC_CL,
                            AUC_CK    = Gout$AUC_CK,
                            AUC_CLu   = Gout$AUC_CLu)
  
  Goutdf <- Goutdf %>% filter (Time == 30) 
  return (list("G" = Goutdf))
  
}

NSC_func <- function (pars, Pred, DOSE) {
  nG <- length(pars)
  NSC_H    = matrix(NA, nrow = length(pars) ,ncol = 4)
  
  for (i in 1:nG) {
    pars.new      <- pars %>% replace(i, log(exp((pars[i]))*1.01))
    Rnew.G         <- Pred(pars.new, DOSE)
    R.G            <- Pred(pars, DOSE)
    delta.pars    <- exp(pars[i])/(exp(pars[i])*0.01)
    
    ## Estimated the AUC
    AUC.Plas.new       =  Rnew.G$G %>% select (AUC_CPlas)
    AUC.Plas.ori       =  R.G$G    %>% select (AUC_CPlas)
    AUC.L.new          =  Rnew.G$G %>% select (AUC_CL)
    AUC.L.ori          =  R.G$G    %>% select (AUC_CL)
    AUC.K.new          =  Rnew.G$G %>% select (AUC_CK)
    AUC.K.ori          =  R.G$G    %>% select (AUC_CK)
    AUC.Lu.new         =  Rnew.G$G %>% select (AUC_CLu)
    AUC.Lu.ori         =  R.G$G    %>% select (AUC_CLu)
    
    delta.AUC.Plas     =  AUC.Plas.new - AUC.Plas.ori
    delta.AUC.L        =  AUC.L.new - AUC.L.ori
    delta.AUC.K        =  AUC.K.new -  AUC.K.ori
    delta.AUC.Lu       =  AUC.Lu.new - AUC.Lu.ori
    
    NSC_H     [i, 1]   <- as.numeric((delta.AUC.Plas/AUC.Plas.ori) * delta.pars)
    NSC_H     [i, 2]   <- as.numeric((delta.AUC.L /AUC.L.ori) * delta.pars)
    NSC_H     [i, 3]   <- as.numeric((delta.AUC.K  /AUC.K.ori) * delta.pars)
    NSC_H     [i, 4]   <- as.numeric((delta.AUC.Lu  /AUC.Lu.ori) * delta.pars)
  }
  return (NSC_H)
}

## Gather all the data for plotting
B <- NSC_func (Human.theta, pred.Human, 3.074E-06)

rownames (B)  = names(Human.theta)
colnames (B)  = c("NSC_CPlas", "NSC_CL", "NSC_CK", "NSC_CLu")

NSC_H <- data.frame(B)

NSC_Clean <- NSC_H %>%
  mutate_all(~ifelse(. == 0, "<1e-5", 
                     ifelse(. > 0.1, round(., 2), 
                            formatC(., format = "e", digits = 2))))
write.csv(NSC_Clean, file = "NSC_H_F.csv", row.names = FALSE)

##################################### Circle barplot function ###############################################
## plot modifed from "R graph gallery: https://www.r-graph-gallery.com/297-circular-barplot-with-groups/ "  #
#############################################################################################################
melt.Plas        = melt(NSC_H[,1]) 
melt.Plas$group  = c("Plasma") 
melt.Liver       = melt(NSC_H[,2])
melt.Liver$group = c("Liver")
melt.Kidney      = melt(NSC_H[,3])
melt.Kidney$group= c("Kidney")
melt.Lung        = melt(NSC_H[,4])
melt.Lung$group  = c("Lung")

melt.data         = rbind (melt.Plas,melt.Liver,melt.Kidney,melt.Lung)
melt.data$par     = rep(rownames(NSC_H),4) 
data.NSC_H     = melt.data%>%filter(abs(value)>=0.2)
data.NSC_H$group <- factor(data.NSC_H$group)

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 4
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data.NSC_H$group), ncol(data.NSC_H)))
colnames(to_add) <- colnames(data.NSC_H)
to_add$group <- rep(levels(data.NSC_H$group), each=empty_bar)
data.NSC_H  <- rbind(data.NSC_H, to_add)
data.NSC_H  <- data.NSC_H %>% arrange(group)
data.NSC_H$id <- seq(1, nrow(data.NSC_H))

# Get the name and the y position of each label
label_data <- data.NSC_H 
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse(angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

#prepare a data frame for base lines
base_data <- data.NSC_H%>% 
  group_by(group)%>% 
  dplyr::summarize(start=min(id)-0.2, end=max(id) - empty_bar+0.2) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

windowsFonts(Times=windowsFont("Times New Roman"))

# Make the plot
p1 <- ggplot(data.NSC_H , aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  
  # Add a val=90/60/30 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 90, xend = start, yend = 90), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 30, xend = start, yend = 30), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data.NSC_H$id),3), y = c(30, 60, 90), label = c("30%", "60%", "90%")  , color="red", size=4, angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=abs(value*100), fill=group), stat="identity", alpha=0.5) +
  ylim(-100,180) +
  labs(tag = "(D)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    text= element_text (family = "Times"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm"), 
    plot.tag = element_text(size = 25, face = "bold", margin = margin(t = 150, l = 150)),  
    plot.tag.position = c(0, 1)  
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=abs(value*100)+10, label=par, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=5, angle= label_data$angle, inherit.aes = FALSE) +
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE) +
  geom_text(data=base_data, aes(x = title, y = -35, label=group), hjust=0.4, colour = "black", alpha=0.8, size=6, fontface="bold", inherit.aes = FALSE)
p1
p1 <- p1 + theme( plot.background = element_rect(fill = "white", color = NA))
ggsave("Figure S4-human-Female.tiff",scale = 1,
       plot = p1,
       path = "D:/zs/PBPK/2025GenX/Artwork",
       width = 24, height = 24, units = "cm",dpi=320)

##############################  MCMC  ##############################################
###normal
mean_value <- 27.28
sd_value   <-	8.184

lower_bound <- qnorm(0.025, mean = mean_value, sd = sd_value)
upper_bound <- qnorm(0.975, mean = mean_value, sd = sd_value)

print(lower_bound)
print(upper_bound)

###lognormal
u=0.03172224
CV=0.3 #or 0.5

mean_value <-log(u/sqrt((1+(CV)^2)))
sd_value <-	sqrt(log(1+CV^2))

lower_bound <- qnorm(0.025, mean = mean_value, sd = sd_value)
upper_bound <- qnorm(0.975, mean = mean_value, sd = sd_value)

print(exp(lower_bound))
print(exp(upper_bound))
print(mean_value)
print(sd_value)

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

N = 1000
PlotDat_A1   <- PBPK_H_pop (theta.int.H,DOSE_A,N = N)
write.csv(PlotDat_A1, file = "PlotDat_A2_F1000.csv", row.names = FALSE)
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
write.csv(H_CPlas, file = "H_CPlas_F.csv", row.names = FALSE)
