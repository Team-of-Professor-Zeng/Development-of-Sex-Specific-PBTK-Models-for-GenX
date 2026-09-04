##  Loading requried R package
library(mrgsolve)    ## R-package for Loading mrgsolve code into r via mcode from the 'mrgsolve' pckage
library(magrittr)    ## R-package for the pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)       ## R-package for transform and summarize tabular data with rows and columns; %>%
library(tidyverse)   ## R-package for transform and summarize tabular data with rows and columns; %>%
library(ggplot2)     ## R-package for GGplot
library(gridExtra)   ## R-package for grid
library(lattice)     ## for plotting the figure
library(FME)         ## R-package for MCMC simulation and model fitting
library(minpack.lm)  ## R-package for model fitting
library(reshape)     ## Package for melt function to reshape the table
library(truncnorm)   ## Package for the truncated normal distribution function   
library(EnvStats)    ## Package for Environmental Statistics, Including US EPA Guidance
library(invgamma)    ## Package for inverse gamma distribution function
library(foreach)     ## Package for parallel computing
library(doParallel)  ## Package for parallel computing
library(bayesplot)   ## Package for MCMC traceplot
library(patchwork)   ## Package for merge images
#fonts()
#loadfonts(device = "win")
windowsFonts("Times New Roman" = windowsFont("Times New Roman"))
windowsFonts(Times = windowsFont("Times New Roman"))
### Male
## Input mrgsolve-based PBPK Model
setwd("C:/Users/15960/Desktop/2025GenX/Modfit/Mice/Male")
source (file = "GenX MMod_M.R")

## Set working direction to the data files
## Build mrgsolve-based PBPK Model
mod <- mcode ("micepbpk", MicePBPK_M.code)

## input data set for model calibration/ oral
Data_A1    <- read.csv(file = "Male.csv")

# Model calibration for PBPK model based on the data of TK study #
# A1. : Mice, oral single dose of 10 mg/kg,   matrix: Plasma, Sampling time: 0.25,0.5,1,2,4,8,12,24,48,72,96,120,144,168h Gannon et al. (2015)
# A2. : Mice, oral single dose of 30 mg/kg,   matrix: Plasma, Sampling time: 0.25,0.5,1,2,4,8,12,24,48,72,96,120,144,168h Gannon et al. (2015)
# B1. : Mice, oral single dose of 12.5 mg/kg, matrix: serum , Sampling time: 0.25,0.5,1,2,4,8,12,24,36,48,72,96,168,336h  Hu et al. (2024)
# D1. : Mice, oral single dose of 0.01 mg/kg, matrix: serum , Sampling time: 6,12,24,36,48,72,96,168h  Zhang et al. (2023)
#===================================================================================================

## Read these datasets and later used in model calibration

OBS.A1  <- Data_A1 %>% filter(Study == 1 & Sample == "Plasma" & Dose == 10)   %>% select(Time = "Time", Plasma = "Conc") # TDosesA = 1
OBS.A2  <- Data_A1 %>% filter(Study == 1 & Sample == "Plasma" & Dose == 30)   %>% select(Time = "Time", Plasma = "Conc") # TDosesA = 1
OBS.B1  <- Data_A1 %>% filter(Study == 2 & Sample == "serum"  & Dose == 12.5) %>% select(Time = "Time", Plasma = "Conc") # TDosesA = 1
OBS.D1  <- Data_A1 %>% filter(Study == 4 & Sample == "serum"  & Dose == 0.01) %>% select(Time = "Time", Plasma = "Conc") # TDosesA = 1

###ADOSE:0.25mg
## Define the prediction function (for least squres fit using levenberg-marquart algorithm)
pred.Mice <- function(pars) {
  
  ## Get out of log domain
  pars %<>% lapply(exp)                 ## return a list of exp (parameters) from log domain
  
  ## Define the exposure scenario 
  BW          = 0.025                                     # Mice body weight
  tinterval   = 24                                        # Time interval
  TDosesA     = 1                                         # Total dosing/Dose times
  #TDosesB     = 28                                        # Total dosing/Dose times
  
  #A1
  PDOSEoral.A1 = 10                                       # Single oral dose from Gannon et al. (2015)
  DOSEoral.A1  = PDOSEoral.A1*BW                          # Amount of oral dose
  ex.oral.A1<- ev  (ID = 1,               ## One individual
                    amt  = DOSEoral.A1,   ## Amount of dose 
                    ii   = tinterval,     ## Time interval
                    addl = TDosesA-1,     ## Addtional doseing 
                    cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                    replicate = FALSE)    ## No replicate
  #A2
  PDOSEoral.A2 = 30                                       # Single oral dose from Gannon et al. (2015)
  DOSEoral.A2  = PDOSEoral.A2*BW                          # Amount of oral dose
  ex.oral.A2<- ev  (ID = 1,              ## One individual
                    amt  = DOSEoral.A2,  ## Amount of dose 
                    ii   = tinterval,    ## Time interval
                    addl = TDosesA-1,    ## Addtional doseing 
                    cmt  = "AST",        ## The dosing comaprtment: AST Stomach  
                    replicate = FALSE)   ## No replicate
  #B1
  PDOSEoral.B1 = 12.5                                     # Single oral dose from Hu et al. (2024)
  DOSEoral.B1  = PDOSEoral.B1*BW                          # Amount of oral dose
  ex.oral.B1<- ev  (ID = 1,              ## One individual
                    amt  = DOSEoral.B1,  ## Amount of dose 
                    ii   = tinterval,    ## Time interval
                    addl = TDosesA-1,    ## Addtional doseing 
                    cmt  = "AST",        ## The dosing comaprtment: AST Stomach  
                    replicate = FALSE)   ## No replicate
  #D1
  PDOSEoral.D1 = 0.01                                     # Single oral dose from Zhang et al. (2023)
  DOSEoral.D1  = PDOSEoral.D1*BW                          # Amount of oral dose
  ex.oral.D1<- ev  (ID = 1,              ## One individual
                    amt  = DOSEoral.D1,  ## Amount of dose 
                    ii   = tinterval,    ## Time interval
                    addl = TDosesA-1,    ## Addtional doseing 
                    cmt  = "AST",        ## The dosing comaprtment: AST Stomach  
                    replicate = FALSE)   ## No replicate
  
  ## set up the exposure time
  tsampA=tgrid(0,tinterval*(TDosesA-1)+24*15,0.5)
  
  ## Get a prediction
  out.A1 <- 
    mod %>%                                                 # model object
    param(pars) %>%                                         # to update the parameters in the model subject
    update(atol = 1E-20, maxsteps=500000) %>%                # solver setting, atol: Absolute tolerance parameter
    mrgsim_d(data = ex.oral.A1, tgrid=tsampA)               # Set up the simulation run
  out.A1<-cbind.data.frame(Time    =out.A1$time, 
                           Plasma  =out.A1$Plasma,
                           BAL     =out.A1$Balance)
  out.A1 <- out.A1 %>% filter (Time > 0)                    # filter the value at time = 0
  
  out.A2 <- 
    mod %>%                                                 # model object
    param(pars) %>%                                         # to update the parameters in the model subject
    update(atol = 1E-20, maxsteps=500000) %>%                # solver setting, atol: Absolute tolerance parameter
    mrgsim_d(data = ex.oral.A2, tgrid=tsampA)               # Set up the simulation run
  out.A2<-cbind.data.frame(Time    =out.A2$time, 
                           Plasma  =out.A2$Plasma,
                           BAL     =out.A2$Balance)
  out.A2 <- out.A2 %>% filter (Time > 0)                    # filter the value at time = 0
  
  out.B1 <- 
    mod %>%                                                 # model object
    param(pars) %>%                                         # to update the parameters in the model subject
    update(atol = 1E-20, maxsteps=500000) %>%                # solver setting, atol: Absolute tolerance parameter
    mrgsim_d(data = ex.oral.B1, tgrid=tsampA)               # Set up the simulation run
  out.B1<-cbind.data.frame(Time=out.B1$time, 
                           Plasma  =out.B1$Plasma,
                           BAL =out.B1$Balance)
  out.B1 <- out.B1 %>% filter (Time > 0)                    # filter the value at time = 0
  
  out.D1 <- 
    mod %>%                                                 # model object
    param(pars) %>%                                         # to update the parameters in the model subject
    update(atol = 1E-20, maxsteps=500000) %>%                # solver setting, atol: Absolute tolerance parameter
    mrgsim_d(data = ex.oral.D1, tgrid=tsampA)               # Set up the simulation run
  out.D1<-cbind.data.frame(Time=out.D1$time, 
                           Plasma  =out.D1$Plasma,
                           BAL =out.D1$Balance)
  out.D1 <- out.D1 %>% filter (Time > 0)                    # filter the value at time = 0
  
  return (list("out.A1"=out.A1, "out.A2"=out.A2,"out.B1"=out.B1, "out.D1"=out.D1))    # Return Dataframe  "out.A1"=out.A1, "out.A2"=out.A2,, "out.C1"=out.C1, "out.C2"=out.C2, "out.C3"=out.C3
}

## initial parmaeters
theta.int <- log(c(
  KurineC                        = 0.0628,                      ## Urinary elimination rate
  Free                           = 0.006,                       ## Free fraction in plasma
  PL                             = 1.339,                       ## Liver /plasma partition coefficient (PC)
  PK                             = 0.854,                       ## Kidney/plasma PC
  PLu                            = 0.431,                       ## Lung  /plasma PC
  PRest                          = 0.595,                       ## Rest of body/plasma PC
  K0C                            = 1,                           ## Rate of absorption of GenX in the stomach
  Kabsc                          = 2.12,                        ## Rate of absorption of GenX in the small intestines
  KunabsC                        = 0.00140                      ## Rate of unabobded dose to appear in feces
))
result <- pred.Mice(theta.int)

## Check mass balance and unfitted curve
plot(result$out.A1$Time,result$out.A1$Plasma,type="l",lwd=2,xlab="Time(hour)",ylab="GenX concentration in plasma")
plot(result$out.A1$Time,result$out.A1$BAL,type="l",lwd=2,xlab="Time(hour)",ylab="Mass Balance")

## Cost fuction (FME) 
## Estimate the model residual by modCost function
MCcost<-function (pars){
  outdf <- pred.Mice (pars)
  
  cost<- modCost  (model = outdf$out.A1, obs = OBS.A1, x ="Time" ,weight = "mean")
  cost<- modCost  (model = outdf$out.A2, obs = OBS.A2, x ="Time" ,weight = "mean",cost = cost)
  cost<- modCost  (model = outdf$out.B1, obs = OBS.B1, x ="Time" ,weight = "mean",cost = cost)
  cost<- modCost  (model = outdf$out.D1, obs = OBS.D1, x ="Time" ,weight = "mean",cost = cost)
  
  return(cost)
}

## Local sensitivity analysis

## Senstivity function (FME) 
## Check the senstive parameters in the model
SnsPlasma <- sensFun(func = MCcost, parms = theta.int, varscale = 1)
Sen       <- summary(SnsPlasma)
plot(Sen)

sen1 <- Sen %>% filter(Mean != 0)
SnsPlasma1 <- SnsPlasma[,row.names(sen1)]
col <-collin(SnsPlasma1)
col1 <-subset(col,collinearity<100)

## Selected senstive parameters
theta <- theta.int[abs(Sen$Mean) > 1.2*mean(abs(Sen$Mean))]
theta 

## Selected parameters
theta.int <- log(c(
  #KurineC                        = 0.0628,                      ## Urinary elimination rate
  Free                           = 0.006,                       ## Free fraction in plasma0.0006-0.06
  PL                             = 1.339,                       ## Liver /plasma partition coefficient (PC)
  #PK                             = 0.854,                       ## Kidney/plasma PC
  PLu                            = 0.431,                       ## Lung  /plasma PC
  PRest                          = 0.595,                       ## Rest of body/plasma PC 0.0595-0.595
  K0C                            = 1,                           ## Rate of absorption of GenX in the stomach
  Kabsc                          = 2.12                        ## Rate of absorption of GenX in the small intestines
  #KunabsC                        = 0.00140                      ## Rate of unabobded dose to appear in feces
))

## PBPK model calibration
## Least squres fit using method "Nelder-Mead" algorithm

Fit<- modFit(f=MCcost, p=theta.int, method ="Marq",
             control = nls.lm.control(nprint=1))
## Summary of fitting results
summary(Fit)               ## Summary of fit 
exp(Fit$par)               ## Get the arithmetic value out of the log domain
Cost <- MCcost(Fit$par)
####################################### Global fitting analysis ################################################################################# 
FDataA <- cbind.data.frame (name= Cost$residuals$name, OBS = Cost$residuals$obs, PRE = Cost$residuals$mod)

## Transformed the predicted and obseved values using log10-sacle to do the plot
FDataA %<>% mutate (Log.OBS = log(OBS,10), Log.PRE = log(PRE,10), Species = "Mice")

## Estimating the R-squared and goodness-of-fit using linear regression model
fit <- lm(Log.OBS ~Log.PRE, data = FDataA)
summary(fit)
adjusted_r_squared <- round(summary(fit)$adj.r.squared, digits = 2)
label_text <- paste("italic(R)^{2} == ", adjusted_r_squared)

######################################## log OPR = log PRE/ log OBS plot ###################################################################
FDataA %<>% mutate(res = residuals(fit), 
                   prediction = predict(fit), 
                   OPR = PRE/OBS,            ## OPR: the ratio of prediction value and observed data
                   log.OPR =  log(OPR,10)) 
write.csv(FDataA, file = 'FDataA_M.csv')
p1_M <- 
  ggplot(FDataA, aes(Log.OBS, Log.PRE)) + ## using log-sacle axis
  geom_abline (intercept = 0, 
               slope     = 1,
               color     ="black", linetype = "dashed",linewidth = 1, alpha = 0.8) +
  geom_point  (aes(shape   = as.factor(name)),size = 3, color = "#e3716e")  +
  annotation_logticks() +
  scale_y_continuous(limits = c(-8,8), labels = scales::math_format(10^.x))+
  scale_x_continuous(limits = c(-8,8),labels = scales::math_format(10^.x))+
  coord_cartesian(clip = "off")
p1_M <- p1_M + 
  theme (
    plot.background         = element_rect (fill="White"),
    text                    = element_text (family = "Times New Roman"),   # text front (Time new roman)
    panel.border            = element_rect (colour = "black", fill=NA, linewidth =2),
    panel.background        = element_rect (fill="White"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.text               = element_text (size   = 28, colour = "black", face = "bold"),    # tick labels along axes 
    axis.title              = element_text (size   = 30, colour = "black", face = "bold"),   # label of axes
    legend.position         =c(0.8, 0.2),
    legend.background = element_rect(fill = "white", color = "white"),
    legend.title = element_text(size = 30, face = "bold"),  
    legend.text = element_text(size = 28)) +
  labs (title = "", x = "Observed values (ug/ml)",  y = "Predicted values (ug/ml)", shape = "Sample")+
  annotate("text", x = -Inf, y = Inf, label = "(A)", hjust = 1.3, vjust = 0.8, size = 15, family = "Times New Roman", colour = "black") +
  annotate("text", x = -4, y = 5, label = label_text , parse = TRUE, size = 15, color = "black", family = "Times New Roman")
print(p1_M)
ggsave("Figure 2-Male.tiff", plot = p1_M,
       path = "D:/zs/PBPK/2025GenX/Artwork",
       width = 24,   height = 24,  units = "cm", dpi = 320)

pct_2e <- mean(FDataA$OPR >= 1/2 & FDataA$OPR <= 2) * 100
pct_3e <- mean(FDataA$OPR >= 1/3 & FDataA$OPR <= 3) * 100
#pct_5e <- mean(FDataA$OPR >= 1/5 & FDataA$OPR <= 5) * 100
#pct_10e <- mean(FDataA$OPR >= 1/10 & FDataA$OPR <= 10) * 100

p2_M <- ggplot(FDataA, aes(x = OPR)) +
  geom_histogram(aes(y = after_stat(count / sum(count) * 100)), 
                 bins = 30, fill = "steelblue", color = "black", alpha = 0.7) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  coord_cartesian(ylim = c(0,  35)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = c(0,0)) +
  coord_cartesian(clip = "off")+
  annotate("text", x = 0.005, y = 34, label = "(A)", 
           hjust = 2, vjust = 0.8, size = 15, family = "Times New Roman", colour = "black") +
  annotate("text", x = 0.005, y = 32, 
           label = paste0("% 2e: ", round(pct_2e, 1), "%"),
           hjust = 0, vjust = 1, size = 9, color = "black", fontface = "bold") +
  annotate("text", x = 0.005, y = 30, 
           label = paste0("% 3e: ", round(pct_3e, 1), "%"),
           hjust = 0, vjust = 1, size = 9, color = "grey", fontface = "bold") + 
  #annotate("text", x = 0.005, y = 28, label = paste0("% 5e: ", round(pct_5e, 1), "%"),
           #hjust = 0, vjust = 1, size = 9, color = "black", fontface = "bold") +
  #annotate("text", x = 0.005, y = 26, label = paste0("% 10e: ", round(pct_10e, 1), "%"),
           #hjust = 0, vjust = 1, size = 9, color = "black", fontface = "bold") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#e3716e", linewidth = 1) + 
  geom_vline(xintercept = c(1/2, 2), linetype = "dashed", color = "black", linewidth = 0.8) + 
  geom_vline(xintercept = c(1/3, 3), linetype = "dashed", color = "grey", linewidth = 0.8) + 
  labs(x = "Predicted / Observed", y = "Proportion (%)") +
  theme_bw(base_family = "Times New Roman") +
  theme(
    plot.background = element_rect(fill = "white"),
    text = element_text(family = "Times New Roman"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.text = element_text(size = 28, colour = "black", face = "bold"),
    axis.title = element_text(size = 30, colour = "black", face = "bold"),
    plot.margin = margin(0.8, 0.6, 0.3, 0.5, "cm"))
print(p2_M)
ggsave("Figure 2-Male1.tiff", plot = p2_M,
       path = "D:/zs/PBPK/2025GenX/Artwork",
       width = 24,   height = 24,  units = "cm", dpi = 320)

grid.arrange(p1,p2,ncol = 2)
##########################################################################################################
p2 <-
  ggplot(FDataA, aes(Log.PRE, log.OPR)) +
  geom_hline(yintercept = log10(2),  linetype = 3, colour = "black", linewidth = 1) +
  geom_hline(yintercept = log10(0.5), linetype = 3, colour = "black", linewidth = 1) +
  geom_hline(yintercept = log10(3),  linetype = 3, colour = "grey50", linewidth = 1) +
  geom_hline(yintercept = log10(0.33), linetype = 3, colour = "grey50", linewidth = 1) +
  annotation_logticks() +
  scale_y_continuous(limits = c(-4,4), labels = scales::math_format(10^.x), breaks = c(-4, -2, 0, 2, 4))+
  scale_x_continuous(limits = c(-7,7),labels = scales::math_format(10^.x))+ 
  annotate("text", x = -Inf, y = Inf, label = "Male (B)", hjust = -0.3, vjust = 1.5, size = 15, family = "Times", colour = "black") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = log10(0.5), ymax = log10(2), fill = "#C3E2EC", alpha = 0.4)+
  geom_point(color   = "steelblue4", aes(shape= as.factor(Species)),size = 3)+
  geom_smooth(se = FALSE, color   = "grey20") 

p2 <- p2 +
  theme (
    plot.background         = element_rect (fill ="White"),
    text                    = element_text (family = "Times"),   # text front (Time new roman)
    panel.border            = element_rect (colour = "black", fill=NA, linewidth =2),
    panel.background        = element_rect (fill ="White"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.text               = element_text (size   = 28, colour = "black", face = "bold"),    # tick labels along axes 
    axis.title              = element_text (size   = 30, colour = "black", face = "bold"),   # label of axes
    legend.position='none') +
  labs (title = "", x = "Predicted values (ug/ml)", y = "Predicted / Observed")
p3 <-ggMarginal(p2, type = "histogram", margins = "y",  yparams = list(binwidth = 0.1, fill = "steelblue4"))
print(p3)

############################################# Time-Conc Plot ################################################################
## Model calibration plot using ggplot2 
Sim.fit.A = pred.Mice(Fit$par)
df.sim.A1 <- cbind.data.frame(Time = Sim.fit.A$out.A1$Time, Plasma   = Sim.fit.A$out.A1$Plasma)
df.sim.A2 <- cbind.data.frame(Time = Sim.fit.A$out.A2$Time, Plasma   = Sim.fit.A$out.A2$Plasma)
df.sim.B1 <- cbind.data.frame(Time = Sim.fit.A$out.B1$Time, Plasma   = Sim.fit.A$out.B1$Plasma)
df.sim.D1 <- cbind.data.frame(Time = Sim.fit.A$out.D1$Time, Plasma   = Sim.fit.A$out.D1$Plasma)

## Setting an initial value
df.sim.A1 <- rbind(data.frame(Time = 0, Plasma = 1e-6),df.sim.A1)
df.sim.A2 <- rbind(data.frame(Time = 0, Plasma = 1e-6),df.sim.A2)
df.sim.B1 <- rbind(data.frame(Time = 0, Plasma = 1e-6),df.sim.B1)
df.sim.D1 <- rbind(data.frame(Time = 0, Plasma = 1e-6),df.sim.D1)

OBS.A1.1  <- Data_A1 %>% filter(Study == 1 & Sample == "Plasma" & Dose == 10)   %>% select(Time = "Time", Plasma = "Conc", SD = "SD") # TDosesA = 1
OBS.A2.1  <- Data_A1 %>% filter(Study == 1 & Sample == "Plasma" & Dose == 30)   %>% select(Time = "Time", Plasma = "Conc", SD = "SD") # TDosesA = 1
OBS.B1.1  <- Data_A1 %>% filter(Study == 2 & Sample == "serum"  & Dose == 12.5) %>% select(Time = "Time", Plasma = "Conc", SD = "SD") # TDosesA = 1
OBS.D1.1  <- Data_A1 %>% filter(Study == 4 & Sample == "serum"  & Dose == 0.01) %>% select(Time = "Time", Plasma = "Conc", SD = "SD") # TDosesA = 1

plot.A1 <- ggplot() +
  geom_line(data = df.sim.A1, aes(Time, Plasma), col = "#5A5A5A", lwd = 1.5) +
  geom_point(data = OBS.A1, aes(Time, Plasma), shape = 17, col = "#7FB2D5", size = 5) +
  geom_errorbar(data = OBS.A1.1, aes(x = Time, ymin = Plasma - SD, ymax = Plasma + SD), col = "#7FB2D5", width = 5, size = 1.5) +
  ylab("GenX Concentration in Plasma (μg/ml)")+
  xlab("Time (hours)") +
  xlim(c(0, 175)) +
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        plot.title = element_text(size = 30, face = "bold", hjust = 0.5, vjust = -8, margin = margin(b = 10) ),
        axis.text    = element_text(size   = 25, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 27),  
        axis.title.y = element_text(face = "bold", size = 27))+
  labs(title = "(A1) Plasma")
print(plot.A1)

plot.A2 <- ggplot() +
  geom_line(data = df.sim.A2, aes(Time, Plasma), col = "#5A5A5A", lwd = 1.5) +
  geom_point(data = OBS.A2, aes(Time, Plasma), shape = 17, col = "#7FB2D5", size = 5) +
  geom_errorbar(data = OBS.A2.1, aes(x = Time, ymin = Plasma - SD, ymax = Plasma + SD), col = "#7FB2D5", width = 5, size = 1.5) +
  ylab("GenX Concentration in Plasma (μg/ml)")+
  xlab("Time (hours)") +
  xlim(c(0, 175)) +
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        plot.title = element_text(size = 30, face = "bold", hjust = 0.5, vjust = -8, margin = margin(b = 10) ),
        axis.text    = element_text(size   = 25, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 27),  
        axis.title.y = element_text(face = "bold", size = 27))+
  labs(title = "(A2) Plasma")
print(plot.A2)

plot.B1 <- ggplot() +
  geom_line(data = df.sim.B1, aes(Time, Plasma), col = "#5A5A5A", lwd = 1.5) +
  geom_point(data = OBS.B1, aes(Time, Plasma), shape = 17, col = "#7FB2D5", size = 5) +
  geom_errorbar(data = OBS.B1.1, aes(x = Time, ymin = Plasma - SD, ymax = Plasma + SD), col = "#7FB2D5", width = 5, size = 1.5) +
  ylab("GenX Concentration in Plasma (μg/ml)")+
  xlab("Time (hours)") +
  xlim(c(0, 350)) +
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        plot.title = element_text(size = 30, face = "bold", hjust = 0.5, vjust = -8, margin = margin(b = 10) ),
        axis.text    = element_text(size   = 25, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 27),  
        axis.title.y = element_text(face = "bold", size = 27))+
  labs(title = "(A3) Plasma")
print(plot.B1)

plot.D1 <- ggplot() +
  geom_line(data = df.sim.D1, aes(Time, Plasma), col = "#5A5A5A", lwd = 1.5) +
  geom_point(data = OBS.D1, aes(Time, Plasma), shape = 17, col = "#7FB2D5", size = 5) +
  geom_errorbar(data = OBS.D1.1, aes(x = Time, ymin = Plasma - SD, ymax = Plasma + SD), col = "#7FB2D5", width = 5, size = 1.5) +
  ylab("GenX Concentration in Plasma (μg/ml)")+
  xlab("Time (hours)") +
  xlim(c(0, 175)) +
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        plot.title = element_text(size = 30, face = "bold", hjust = 0.5, vjust = -8, margin = margin(b = 10) ),
        axis.text    = element_text(size   = 25, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 27),  
        axis.title.y = element_text(face = "bold", size = 27))+
  labs(title = "(A4) Plasma")
print(plot.D1)

plot.A2 <- plot.A2 + theme(axis.title.y = element_blank())
plot.B1 <- plot.B1 + theme(axis.title.y = element_blank())
plot.D1 <- plot.D1 + theme(axis.title.y = element_blank())
ggsave("Figure S2-Male.tiff", grid.arrange(plot.A1, plot.A2, plot.B1, plot.D1, ncol = 4),
       path = "D:/zs/PBPK/2025GenX/Artwork",
       width = 80, height = 20, units = "cm", dpi = 320)

## Save the fitting results to RDS files
saveRDS(Fit, file = "Fit_R.rds") 
Fit_R <- readRDS("Fit_R.rds")

##ABSTRACT
plot.A2 <- ggplot() +
  geom_line(data = df.sim.A2, aes(Time, Plasma), col = "#5A5A5A", lwd = 2) +
  geom_point(data = OBS.A2, aes(Time, Plasma), col = "#e3716e", size = 7) +
  ylab("Concentration (μg/ml)")+
  xlab("Time (hours)") +
  xlim(c(0, 175)) +
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 3),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.text    = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(face = "bold", size = 30, margin = margin(t=15, b = 15)),  
        axis.title.y = element_text(face = "bold", size = 30, margin = margin(l=15, r = 15)))+
  labs(title = "")
print(plot.A2)
ggsave("Figure abstract.tiff", grid.arrange(plot.A2, ncol = 1),
       path = "D:/zs/PBPK/2025GenX/Artwork",
       width = 15, height = 15, units = "cm", dpi = 320)
#################################################################################################################################
###NSC                                                                                                                          #
#################################################################################################################################
mice.theta.G <- log(c(
  # Physiological parameters
  BW                             = 0.025,
  QCC                            = 16.5,
  QLC                            = 0.161,
  QLuC                           = 0.005,
  QKC                            = 0.091,
  Htc                            = 0.48,
  VPlasC                         = 0.049,
  VLC                            = 0.055,
  VLuC                           = 0.007,
  VKC                            = 0.017,
  VFilC                          = 0.0017,
  FVBK                           = 0.160,
  GFRC                           = 62.1,
  GEC                            = 0.54,
  
  # Chemical-specific parameters (final mean values)
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

pred.mice <- function(pars, DOSE) {#= c(10, 30, 12.5, 0.01)
  
  ## Get out of log domain
  pars <- exp(pars)                   ## Return a list of exp  from log scale
  
  ## Exposure scenario 
  BW               = 0.025                            ## kg, Rat body weight                                      
  tinterval        = 24                               ## hr, Time interval                                 
  TDoses           = 1                                ## Dose times                                    
  GDOSE             = DOSE                             ## Input oral dose
  DOSEoral         = GDOSE*BW                          ## Amount of oral dose
  
  # To create a exposure scenario
  ex.oral <- ev (ID   = 1,             ## One individual
                 amt  = DOSEoral,      ## Amount of dose 
                 ii   = tinterval,     ## Time interval
                 addl = TDoses - 1,    ## Addtional doseing 
                 cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                 replicate = TRUE)     ## No replicate
  
  # set up the exposure time
  tsamp = tgrid (0,tinterval*(TDoses-1)+24*15,0.5)       ## simulation 24*30 hr (30 days)
  
  out <- 
    mod %>%                                                 # model object
    param(pars) %>%                                         # to update the parameters in the model subject
    update(atol = 1E-20, maxsteps=500000) %>%               # solver setting, atol: Absolute tolerance parameter
    mrgsim_d(data = ex.oral, tgrid=tsamp) %>%
    filter(time!=0)
  
  outdfA = cbind.data.frame(Time      = out$time, 
                            Plasma    = out$Plasma,
                            AUC_CPlas = out$AUC_CA,
                            AUC_CL    = out$AUC_CL,
                            AUC_CK    = out$AUC_CK,
                            AUC_CLu   = out$AUC_CLu)   
  outdfA <- outdfA %>% filter (Time == 15*24) ## set up the output time
  return (list("G" = outdfA))
}
#result  <-  pred.mice(Fit$par)

par.G <- Fit_R$par
mice.theta.G[names(par.G)]  <- as.numeric(par.G) 

NSC_func <- function (pars, Pred, DOSE) {
  nG <- length(pars)
  NSC_GCA     = matrix(NA, nrow = length(pars) , ncol = 4)
  
  for (i in 1:nG) {
    pars.new      <- pars %>% replace(i, log(exp((pars[i]))*1.01))
    Mnew.G         <- Pred(pars.new, DOSE)
    M.G            <- Pred(pars, DOSE)
    delta.Gpars    <- exp(pars[i])/(exp(pars[i])*0.01)
    
    ## Estimated the AUC
    AUC.GPlas.new       =  Mnew.G$G %>% select (AUC_CPlas)
    AUC.GPlas.ori       =  M.G   $G %>% select (AUC_CPlas)
    AUC.GL.new          =  Mnew.G$G %>% select ( AUC_CL)
    AUC.GL.ori          =  M.G   $G %>% select ( AUC_CL)
    AUC.GK.new          =  Mnew.G$G %>% select ( AUC_CK)
    AUC.GK.ori          =  M.G   $G %>% select ( AUC_CK)
    AUC.GLu.new         =  Mnew.G$G %>% select ( AUC_CLu)
    AUC.GLu.ori         =  M.G   $G %>% select ( AUC_CLu)
    
    delta.AUC.GPlas     =  AUC.GPlas.new - AUC.GPlas.ori
    delta.AUC.GL        =  AUC.GL.new    -  AUC.GL.ori
    delta.AUC.GK        =  AUC.GK.new    -  AUC.GK.ori
    delta.AUC.GLu       =  AUC.GLu.new   -  AUC.GLu.ori
    
    NSC_GCA     [i, 1]   <- as.numeric((delta.AUC.GPlas/AUC.GPlas.ori) * delta.Gpars)
    NSC_GCA     [i, 2]   <- as.numeric((delta.AUC.GL   /AUC.GL.ori)    * delta.Gpars)
    NSC_GCA     [i, 3]   <- as.numeric((delta.AUC.GK   /AUC.GK.ori)    * delta.Gpars)
    NSC_GCA     [i, 4]   <- as.numeric((delta.AUC.GLu  /AUC.GLu.ori)   * delta.Gpars)
  }
  return (list(NSC_GCA = NSC_GCA))
}

# Model results
A <- NSC_func (mice.theta.G, pred.mice, 10)

rownames (A$NSC_GCA)  = names(mice.theta.G)
colnames (A$NSC_GCA)  = c("NSC_CPlas", "NSC_CL", "NSC_CK", "NSC_CLu")

NSC_GCA_M <- data.frame(A$NSC_GCA)

NSC_Mice <- NSC_GCA_M %>%
  rowwise() %>%
  filter(any(sapply(cols(), function(x) any(abs(x) >= 0.3)))) %>%
  ungroup()

NSC_GCA_M <- NSC_GCA_M %>%
  mutate_all(~replace(., . == 0, "<1e-5"))
write.csv(NSC_GCA_M, file = 'NSC.csv')
NSC_GCA_M <- read.csv(file = "NSC.csv", row.names = 1) 
##################################### Circle barplot function ###############################################
## plot modifed from "R graph gallery: https://www.r-graph-gallery.com/297-circular-barplot-with-groups/ "  #
#############################################################################################################
melt.Plas         = melt(NSC_GCA_M[,1]) 
melt.Plas$group   = c("Plasma") 
melt.Liver        = melt(NSC_GCA_M[,2])
melt.Liver$group  = c("Liver")
melt.Kidney       = melt(NSC_GCA_M[,3])
melt.Kidney$group = c("Kidney")
melt.Lung         = melt(NSC_GCA_M[,4])
melt.Lung$group   = c("Lung")

melt.data         = rbind (melt.Plas,melt.Liver,melt.Kidney,melt.Lung)
melt.data$par     = rep(rownames(NSC_GCA_M),4) 
str(melt.data$value)
melt.data$value <- as.numeric(as.character(melt.data$value))
data              = melt.data%>%filter(abs(value)>=0.3)
data$group        <- factor(data$group)

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 4
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)))
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data%>% arrange(group)
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

#prepare a data frame for base lines
base_data <- data %>% 
  group_by(group) %>% 
  summarize(start=min(id)-0.2, end=max(id) - empty_bar+0.2) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

windowsFonts(Times=windowsFont("Times New Roman"))

# Make the plot
p1 <- ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  
  # Add a val=90/60/30 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 90, xend = start, yend = 90), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 30, xend = start, yend = 30), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data$id),3), y = c(30, 60, 90), label = c("30%", "60%", "90%")  , color="red", size=4, angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=abs(value*100), fill=group), stat="identity", alpha=0.5) +
  ylim(-100,230) +
  labs(tag = "(A)") +
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
  geom_text(data=label_data, aes(x=id,  y=abs(value*100)+10, label=par, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=5, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -38, label=group), hjust=0.4, colour = "black", alpha=0.8, size=6, fontface="bold", inherit.aes = FALSE)
p1
p1 <- p1 + theme( plot.background = element_rect(fill = "white", color = NA))
ggsave("Figure S5-mice-Male.tiff",scale = 1,
       plot = p1,
       path = "D:/zs/PBPK/2025GenX/Artwork",
       width = 24, height = 24, units = "cm",dpi=320)

#############################################################################################################################
##MCMC                                                                                                                      #
## normal AND lognormal Calculation #########################################################################################
###normal
mean_value <- 62.1
sd_value   <-	18.63

lower_bound <- qnorm(0.025, mean = mean_value, sd = sd_value)
upper_bound <- qnorm(0.975, mean = mean_value, sd = sd_value)

print(lower_bound)
print(upper_bound)

###lognormal
u=0.07103239
CV=0.3  #or 0.5

mean_value <-log(u/sqrt((1+(CV)^2)))
sd_value <-	sqrt(log(1+CV^2))

lower_bound <- qnorm(0.025, mean = mean_value, sd = sd_value)
upper_bound <- qnorm(0.975, mean = mean_value, sd = sd_value)

print(exp(lower_bound))
print(exp(upper_bound))
print(mean_value)
print(sd_value)

## A1 MCMC ###################################################
Mice_A1 <- function (pars, N) {
  
  pars <- exp(pars)           ## Return a list of exp from log scale
  
  ## Exposure scenario 
  BW            = 0.025                            ## kg, Mice body weight                                      
  tinterval     = 24                               ## hr, Time interval                                 
  TDoses        = 1                                ## Dose times                                    
  DOSE          = 10                               ## Single oral dose from Gannon et al. (2015)
  DOSEoral      = DOSE*BW                          ## Amount of oral dose
  
  ## Amount of oral dose
  idata <- 
    tibble(ID = 1:N) %>% 
    mutate( BW        = rnormTrunc  (N, min = 0.0103, max = 0.0397, mean = 0.025,      sd = 0.0075),   #CV=0.3
            VKC       = rnormTrunc  (N, min = 0.007 , max = 0.0270, mean = 0.017,      sd = 0.0051),   #CV=0.3
            GFRC      = rnormTrunc  (N, min = 25.586, max = 98.614, mean = 62.1 ,      sd = 18.63) ,   #CV=0.3
            Free                           = 0.0106,                 ### fitting parameters
            PL                             = 0.753,                  ### fitting parameters
            PLu                            = 2.0399,                 ### fitting parameters
            PRest                          = 0.14 ,                  ### fitting parameters
            K0C                            = 1.838,                  ### fitting parameters               
            Kabsc                          = 1.0461,                 ### fitting parameters 
            DOSEoral  = DOSE*BW 
    )   
  
  ex.oral_1<- ev  (ID = 1:N,              ## One individual
                   time = 0,             ## Dosed strat time 
                   amt  = idata$DOSEoral,## Amount of dose 
                   ii   = tinterval,     ## Time interval
                   addl = TDoses-1,      ## Addtional doseing 
                   cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                   replicate = FALSE)    ## No replicate
  ex.oral_2 <- ev (ID   = 1,              ## One individual
                   time = 0,             ## Dosed strat time 
                   amt  = DOSEoral,      ## Amount of dose 
                   ii   = tinterval,     ## Time interval
                   addl = TDoses - 1,    ## Addtional doseing 
                   cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                   replicate = FALSE)    ## No replicate
  
  ex_1 <- ex.oral_1
  ex_2 <- ex.oral_2
  
  ## set up the exposure time
  tsamp = tgrid(0,tinterval*(TDoses-1)+24*8,0.5)
  
  # Combine data and run the simulation
  out_1 <- mod %>%data_set(ex_1) %>%
    idata_set(idata) %>% 
    update(atol = 1e-6, maxstep = 50000) %>%
    mrgsim(obsonly=TRUE, tgrid = tsamp)
  outdf_1 = cbind.data.frame ( ID     = out_1$ID,
                               Time   = out_1$time/24, 
                               CPlas  = out_1$Plasma, 
                               CL     = out_1$Liver,
                               CK     = out_1$Kidney,
                               CLu    = out_1$Lung)
  out_2 <- 
    mod %>%
    param (pars) %>%
    update(atol = 1E-6, maxsteps = 50000) %>%          
    mrgsim_d (data = ex_2, tgrid = tsamp)
  outdf_2 = cbind.data.frame(Time     = out_2$time/24, 
                             CPlas    = out_2$Plasma, 
                             CL       = out_2$Liver,
                             CK       = out_2$Kidney,
                             CLu      = out_2$Lung)
  
  return (list(outdf_1, outdf_2)) # Return outdf
}

R_Gpars <- Fit_R$par
N = 1000

PlotDat_A1     <- Mice_A1 (R_Gpars,N = N)[[1]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CPlas)  %>% mutate( Tissue = "Plasma")
PlotDat_A1_m   <- Mice_A1 (R_Gpars,N = N)[[2]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CPlas)  %>% mutate( Tissue = "Plasma")
OBS.A1  <- Data_A1 %>% filter(Study == 1 & Sample == "Plasma" & Dose == 10)   %>% select(Time = "Time", Plasma = "Conc", SD = "SD") # TDosesA = 1

PlotDat_A1_summary <- PlotDat_A1 %>%
  group_by(Time) %>%
  summarize(
    median_est = median(Conc, na.rm = TRUE),
    ci_q1 = quantile(Conc, probs = 0.25, names = FALSE, na.rm = TRUE),
    ci_q3 = quantile(Conc, probs = 0.75, names = FALSE, na.rm = TRUE),
    ci_10 = quantile(Conc, probs = 0.10, names = FALSE, na.rm = TRUE),
    ci_90 = quantile(Conc, probs = 0.90, names = FALSE, na.rm = TRUE),
    ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
    ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)
  ) %>%
  ungroup()

plot.A1_Plas <- 
  ggplot() + 
  geom_ribbon  (data = PlotDat_A1_summary, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), fill="#B0D9A5", alpha=0.3) +
  geom_ribbon  (data = PlotDat_A1_summary, aes(x = Time, ymin = ci_q1, ymax = ci_q3), fill="#1E803D", alpha = 0.3) +
  geom_line    (data = PlotDat_A1_m, aes(x = Time, y = Conc), colour = "#5A5A5A", lwd = 0.8) +
  geom_point   (data = OBS.A1, aes(x = Time/24, y = Plasma), colour = "#1E803D", size = 2.5) +
  geom_errorbar(data = OBS.A1, aes(x = Time/24, ymin= Plasma-SD, ymax = Plasma+SD), col = "#1E803D", width = 0.05, size = 0.8)+
  ylab("GenX Concentration in Plasma (μg/ml)")+
  xlab("Time (days)")+
  scale_x_continuous(trans = "log1p") +  # log Time
  scale_y_continuous(limits = c(-3, 70), expand = c(0, 0))+
  theme_bw() +
  theme(panel.border = element_rect(color = "grey10",fill = NA, linewidth = 1),
        strip.text   = element_text(size = rel(2), colour = "grey10"),
        axis.text    = element_text(size = rel(2), colour = "black"),
        axis.line    = element_line(color = "black", size = 0.3),
        axis.title   = element_text(size = 25, colour = "black", face = "bold"),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.5, vjust = -8, margin = margin(b = 10) ),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+
  labs(title = "(A1) Plasma")
plot.A1_Plas

## A2 MCMC ########################################################
Mice_A2 <- function (pars, N) {
  
  pars <- exp(pars)           ## Return a list of exp from log scale
  
  ## Exposure scenario 
  BW            = 0.025                            ## kg, Rat body weight                                      
  tinterval     = 24                               ## hr, Time interval                                 
  TDoses        = 1                                ## Dose times                                    
  DOSE          = 30                               ## Single oral dose from Gannon et al. (2015)
  DOSEoral      = DOSE*BW                          ## Amount of oral dose
  
  ## Amount of oral dose
  idata <- 
    tibble(ID = 1:N) %>% 
    mutate( BW        = rnormTrunc  (N, min = 0.0103, max = 0.0397, mean = 0.025,      sd = 0.0075),   #CV=0.3
            VKC       = rnormTrunc  (N, min = 0.007 , max = 0.0270, mean = 0.017,      sd = 0.0051),   #CV=0.3
            GFRC      = rnormTrunc  (N, min = 25.586, max = 98.614, mean = 62.1 ,      sd = 18.63) ,   #CV=0.3
            Free                           = 0.0106,                 ### fitting parameters
            PL                             = 0.753,                  ### fitting parameters
            PLu                            = 2.0399,                 ### fitting parameters
            PRest                          = 0.14 ,                  ### fitting parameters
            K0C                            = 1.838,                  ### fitting parameters               
            Kabsc                          = 1.0461,                 ### fitting parameters 
            DOSEoral  = DOSE*BW 
    )   
  
  ex.oral_1<- ev  (ID = 1:N,             ## One individual
                   time = 0,             ## Dossed strat time 
                   amt  = idata$DOSEoral,## Amount of dose 
                   ii   = tinterval,     ## Time interval
                   addl = TDoses-1,      ## Addtional doseing 
                   cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                   replicate = FALSE)    ## No replicate
  ex.oral_2 <- ev (ID   = 1,             ## One individual
                   time = 0,             ## Dossed strat time 
                   amt  = DOSEoral,      ## Amount of dose 
                   ii   = tinterval,     ## Time interval
                   addl = TDoses - 1,    ## Addtional doseing 
                   cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                   replicate = FALSE)    ## No replicate
  
  ex_1 <- ex.oral_1
  ex_2 <- ex.oral_2
  
  ## set up the exposure time
  tsamp = tgrid(0,tinterval*(TDoses-1)+24*8,0.5)
  
  # Combine data and run the simulation
  out_1 <- mod %>%data_set(ex_1) %>%
    idata_set(idata) %>% 
    update(atol = 1e-6, maxstep = 50000) %>%
    mrgsim(obsonly=TRUE, tgrid = tsamp)
  outdf_1 = cbind.data.frame ( ID     = out_1$ID,
                               Time   = out_1$time/24, 
                               CPlas  = out_1$Plasma, 
                               CL     = out_1$Liver,
                               CK     = out_1$Kidney,
                               CLu    = out_1$Lung)
  out_2 <- 
    mod %>%
    param (pars) %>%
    update(atol = 1E-6, maxsteps = 50000) %>%          
    mrgsim_d (data = ex_2, tgrid = tsamp)
  outdf_2 = cbind.data.frame(Time     = out_2$time/24, 
                             CPlas    = out_2$Plasma, 
                             CL       = out_2$Liver,
                             CK       = out_2$Kidney,
                             CLu      = out_2$Lung)
  
  return (list(outdf_1, outdf_2)) # Return outdf
}

R_Gpars <- Fit_R$par
N = 1000

PlotDat_A2     <- Mice_A2 (R_Gpars,N = N)[[1]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CPlas)  %>% mutate( Tissue = "Plasma")
PlotDat_A2_m   <- Mice_A2 (R_Gpars,N = N)[[2]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CPlas)  %>% mutate( Tissue = "Plasma")
OBS.A2  <- Data_A1 %>% filter(Study == 1 & Sample == "Plasma" & Dose == 30)   %>% select(Time = "Time", Plasma = "Conc", SD = "SD") # TDosesA = 1

PlotDat_A2_summary <- PlotDat_A2 %>%
  group_by(Time) %>%
  summarize(
    median_est = median(Conc, na.rm = TRUE),
    ci_q1 = quantile(Conc, probs = 0.25, names = FALSE, na.rm = TRUE),
    ci_q3 = quantile(Conc, probs = 0.75, names = FALSE, na.rm = TRUE),
    ci_10 = quantile(Conc, probs = 0.10, names = FALSE, na.rm = TRUE),
    ci_90 = quantile(Conc, probs = 0.90, names = FALSE, na.rm = TRUE),
    ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
    ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)
  ) %>%
  ungroup()

plot.A2_Plas <- 
  ggplot() + 
  geom_ribbon  (data = PlotDat_A2_summary, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), fill="#B0D9A5", alpha=0.3) +
  geom_ribbon  (data = PlotDat_A2_summary, aes(x = Time, ymin = ci_q1, ymax = ci_q3), fill="#1E803D", alpha = 0.3) +
  geom_line    (data = PlotDat_A2_m, aes(x = Time, y = Conc), colour = "#5A5A5A", size = 0.8) +
  geom_point   (data = OBS.A2, aes(x = Time/24, y = Plasma), colour = "#1E803D", size = 2.5) +
  geom_errorbar(data = OBS.A2, aes(x = Time/24, ymin= Plasma-SD, ymax = Plasma+SD), col = "#1E803D", width = 0.05, size = 0.8)+
  ylab("GenX Concentration in Plasma (μg/ml)")+
  xlab("Time (days)")+
  scale_x_continuous(trans = "log1p") +  # log Time
  scale_y_continuous(limits = c(-3, 170), expand = c(0, 0))+
  theme_bw() +
  theme(panel.border = element_rect(color = "grey10",fill = NA, linewidth = 1),
        strip.text   = element_text(size = rel(2), colour = "grey10"),
        axis.text    = element_text(size = rel(2), colour = "black"),
        axis.line    = element_line(color = "black", size = 0.3),
        axis.title   = element_text(size = 25, colour = "black", face = "bold"),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.5, vjust = -8, margin = margin(b = 10) ),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+
  labs(title = "(A2) Plasma")
plot.A2_Plas

## B1 MCMC #######################################################
Mice_B1 <- function (pars, N) {
  
  pars <- exp(pars)           ## Return a list of exp from log scale
  
  ## Exposure scenario 
  BW            = 0.025                            ## kg, Rat body weight                                      
  tinterval     = 24                               ## hr, Time interval                                 
  TDoses        = 1                                ## Dose times                                    
  DOSE          = 12.5                             ## Single oral dose from Hu et al. (2024)
  DOSEoral      = DOSE*BW                          ## Amount of oral dose
  
  ## Amount of oral dose
  idata <- 
    tibble(ID = 1:N) %>% 
    mutate( BW        = rnormTrunc  (N, min = 0.0103, max = 0.0397, mean = 0.025,      sd = 0.0075),   #CV=0.3
            VKC       = rnormTrunc  (N, min = 0.007 , max = 0.0270, mean = 0.017,      sd = 0.0051),   #CV=0.3
            GFRC      = rnormTrunc  (N, min = 25.586, max = 98.614, mean = 62.1 ,      sd = 18.63) ,   #CV=0.3
            Free                           = 0.0106,                 ### fitting parameters
            PL                             = 0.753,                  ### fitting parameters
            PLu                            = 2.0399,                 ### fitting parameters
            PRest                          = 0.14 ,                  ### fitting parameters
            K0C                            = 1.838,                  ### fitting parameters               
            Kabsc                          = 1.0461,                 ### fitting parameters 
            DOSEoral  = DOSE*BW 
    )   
  
  ex.oral_1<- ev  (ID = 1:N,             ## One individual
                   time = 0,             ## Dossed strat time 
                   amt  = idata$DOSEoral,## Amount of dose 
                   ii   = tinterval,     ## Time interval
                   addl = TDoses-1,      ## Addtional doseing 
                   cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                   replicate = FALSE)    ## No replicate
  ex.oral_2 <- ev (ID   = 1,             ## One individual
                   time = 0,             ## Dossed strat time 
                   amt  = DOSEoral,      ## Amount of dose 
                   ii   = tinterval,     ## Time interval
                   addl = TDoses - 1,    ## Addtional doseing 
                   cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                   replicate = FALSE)    ## No replicate
  
  ex_1 <- ex.oral_1
  ex_2 <- ex.oral_2
  
  ## set up the exposure time
  tsamp = tgrid(0,tinterval*(TDoses-1)+24*15,0.5)
  
  # Combine data and run the simulation
  out_1 <- mod %>%data_set(ex_1) %>%
    idata_set(idata) %>% 
    update(atol = 1e-6, maxstep = 50000) %>%
    mrgsim(obsonly=TRUE, tgrid = tsamp)
  outdf_1 = cbind.data.frame ( ID     = out_1$ID,
                               Time   = out_1$time/24, 
                               CPlas  = out_1$Plasma, 
                               CL     = out_1$Liver,
                               CK     = out_1$Kidney,
                               CLu    = out_1$Lung)
  out_2 <- 
    mod %>%
    param (pars) %>%
    update(atol = 1E-6, maxsteps = 50000) %>%          
    mrgsim_d (data = ex_2, tgrid = tsamp)
  outdf_2 = cbind.data.frame(Time     = out_2$time/24, 
                             CPlas    = out_2$Plasma, 
                             CL       = out_2$Liver,
                             CK       = out_2$Kidney,
                             CLu      = out_2$Lung)
  
  return (list(outdf_1, outdf_2)) # Return outdf
}

R_Gpars <- Fit_R$par
N = 1000

PlotDat_B1     <- Mice_B1 (R_Gpars,N = N)[[1]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CPlas)  %>% mutate( Tissue = "Plasma")
PlotDat_B1_m   <- Mice_B1 (R_Gpars,N = N)[[2]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CPlas)  %>% mutate( Tissue = "Plasma")
OBS.B1  <- Data_A1 %>% filter(Study == 2 & Sample == "serum"  & Dose == 12.5) %>% select(Time = "Time", Plasma = "Conc", SD = "SD") # TDosesA = 1

PlotDat_B1_summary <- PlotDat_B1 %>%
  group_by(Time) %>%
  summarize(
    median_est = median(Conc, na.rm = TRUE),
    ci_q1 = quantile(Conc, probs = 0.25, names = FALSE, na.rm = TRUE),
    ci_q3 = quantile(Conc, probs = 0.75, names = FALSE, na.rm = TRUE),
    ci_10 = quantile(Conc, probs = 0.10, names = FALSE, na.rm = TRUE),
    ci_90 = quantile(Conc, probs = 0.90, names = FALSE, na.rm = TRUE),
    ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
    ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)
  ) %>%
  ungroup()

plot.B1_Plas <- 
  ggplot() + 
  geom_ribbon  (data = PlotDat_B1_summary, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), fill="#B0D9A5", alpha=0.3) +
  geom_ribbon  (data = PlotDat_B1_summary, aes(x = Time, ymin = ci_q1, ymax = ci_q3), fill="#1E803D", alpha = 0.3) +
  geom_line    (data = PlotDat_B1_m, aes(x = Time, y = Conc), colour = "#5A5A5A", lwd = 0.8) +
  geom_point   (data = OBS.B1, aes(x = Time/24, y = Plasma), colour = "#1E803D", size = 2.5) +
  geom_errorbar(data = OBS.B1, aes(x = Time/24, ymin= Plasma-SD, ymax = Plasma+SD), col = "#1E803D", width = 0.05, size = 0.8)+
  ylab("GenX Concentration in Plasma (μg/ml)")+
  xlab("Time (days)")+
  scale_x_continuous(trans = "log1p") +  # log Time
  scale_y_continuous(limits = c(-3, 110), expand = c(0, 0))+
  theme_bw() +
  theme(panel.border = element_rect(color = "grey10",fill = NA, linewidth = 1),
        strip.text   = element_text(size = rel(2), colour = "grey10"),
        axis.text    = element_text(size = rel(2), colour = "black"),
        axis.line    = element_line(color = "black", size = 0.3),
        axis.title   = element_text(size = 25, colour = "black", face = "bold"),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.5, vjust = -8, margin = margin(b = 10) ),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+
  labs(title = "(A3) Plasma")
plot.B1_Plas

## D1 MCMC ########################################################
Mice_D1 <- function (pars, N) {
  
  pars <- exp(pars)           ## Return a list of exp from log scale
  
  ## Exposure scenario
  BW            = 0.025                            ## kg, Rat body weight                                      
  tinterval     = 24                               ## hr, Time interval                                 
  TDoses        = 1                                ## Dose times                                    
  DOSE          = 0.01                             ## Single oral dose from Zhang et al. (2023)
  DOSEoral      = DOSE*BW                          ## Amount of oral dose
  
  ## Amount of oral dose
  idata <- 
    tibble(ID = 1:N) %>% 
    mutate( BW        = rnormTrunc  (N, min = 0.0103, max = 0.0397, mean = 0.025,      sd = 0.0075),   #CV=0.3
            VKC       = rnormTrunc  (N, min = 0.007 , max = 0.0270, mean = 0.017,      sd = 0.0051),   #CV=0.3
            GFRC      = rnormTrunc  (N, min = 25.586, max = 98.614, mean = 62.1 ,      sd = 18.63) ,   #CV=0.3
            Free                           = 0.0106,                 ### fitting parameters
            PL                             = 0.753,                  ### fitting parameters
            PLu                            = 2.0399,                 ### fitting parameters
            PRest                          = 0.14 ,                  ### fitting parameters
            K0C                            = 1.838,                  ### fitting parameters               
            Kabsc                          = 1.0461,                 ### fitting parameters 
            DOSEoral  = DOSE*BW 
    )   
  
  ex.oral_1<- ev  (ID = 1:N,             ## One individual
                   time = 0,             ## Dossed strat time 
                   amt  = idata$DOSEoral,## Amount of dose 
                   ii   = tinterval,     ## Time interval
                   addl = TDoses-1,      ## Addtional doseing 
                   cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                   replicate = FALSE)    ## No replicate
  ex.oral_2 <- ev (ID   = 1,             ## One individual
                   time = 0,             ## Dossed strat time 
                   amt  = DOSEoral,      ## Amount of dose 
                   ii   = tinterval,     ## Time interval
                   addl = TDoses - 1,    ## Addtional doseing 
                   cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                   replicate = FALSE)    ## No replicate
  
  ex_1 <- ex.oral_1
  ex_2 <- ex.oral_2
  
  ## set up the exposure time
  tsamp = tgrid(0,tinterval*(TDoses-1)+24*8,0.5)
  
  # Combine data and run the simulation
  out_1 <- mod %>%data_set(ex_1) %>%
    idata_set(idata) %>% 
    update(atol = 1e-6, maxstep = 50000) %>%
    mrgsim(obsonly=TRUE, tgrid = tsamp)
  outdf_1 = cbind.data.frame ( ID     = out_1$ID,
                               Time   = out_1$time/24, 
                               CPlas  = out_1$Plasma, 
                               CL     = out_1$Liver,
                               CK     = out_1$Kidney,
                               CLu    = out_1$Lung)
  out_2 <- 
    mod %>%
    param (pars) %>%
    update(atol = 1E-6, maxsteps = 50000) %>%          
    mrgsim_d (data = ex_2, tgrid = tsamp)
  outdf_2 = cbind.data.frame(Time     = out_2$time/24, 
                             CPlas    = out_2$Plasma, 
                             CL       = out_2$Liver,
                             CK       = out_2$Kidney,
                             CLu      = out_2$Lung)
  
  return (list(outdf_1, outdf_2)) # Return outdf
}

R_Gpars <- Fit_R$par
N = 1000

PlotDat_D1     <- Mice_D1 (R_Gpars,N = N)[[1]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CPlas)  %>% mutate( Tissue = "Plasma")
PlotDat_D1_m   <- Mice_D1 (R_Gpars,N = N)[[2]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CPlas)  %>% mutate( Tissue = "Plasma")
OBS.D1  <- Data_A1 %>% filter(Study == 4 & Sample == "serum"  & Dose == 0.01) %>% select(Time = "Time", Plasma = "Conc", SD = "SD") # TDosesA = 1

PlotDat_D1_summary <- PlotDat_D1 %>%
  group_by(Time) %>%
  summarize(
    median_est = median(Conc, na.rm = TRUE),
    ci_q1 = quantile(Conc, probs = 0.25, names = FALSE, na.rm = TRUE),
    ci_q3 = quantile(Conc, probs = 0.75, names = FALSE, na.rm = TRUE),
    ci_10 = quantile(Conc, probs = 0.10, names = FALSE, na.rm = TRUE),
    ci_90 = quantile(Conc, probs = 0.90, names = FALSE, na.rm = TRUE),
    ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
    ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)
  ) %>%
  ungroup()

plot.D1_Plas <- 
  ggplot() + 
  geom_ribbon  (data = PlotDat_D1_summary, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), fill="#B0D9A5", alpha=0.3) +
  geom_ribbon  (data = PlotDat_D1_summary, aes(x = Time, ymin = ci_q1, ymax = ci_q3), fill="#1E803D", alpha = 0.3) +
  geom_line    (data = PlotDat_D1_m, aes(x = Time, y = Conc), colour = "#5A5A5A", lwd = 0.8) +
  geom_point   (data = OBS.D1, aes(x = Time/24, y = Plasma), colour = "#1E803D", size = 2.5) +
  geom_errorbar(data = OBS.D1, aes(x = Time/24, ymin= Plasma-SD, ymax = Plasma+SD), col = "#1E803D", width = 0.05, size = 0.8)+
  ylab("GenX Concentration in Plasma (μg/ml)")+
  xlab("Time (days)")+
  scale_x_continuous(trans = "log1p") +  # log Time
  scale_y_continuous(limits = c(-0.001, 0.07), expand = c(0, 0))+
  theme_bw() +
  theme(panel.border = element_rect(color = "grey10",fill = NA, linewidth = 1),
        strip.text   = element_text(size = rel(2), colour = "grey10"),
        axis.text    = element_text(size = rel(2), colour = "black"),
        axis.line    = element_line(color = "black", size = 0.3),
        axis.title   = element_text(size = 25, colour = "black", face = "bold"),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.5, vjust = -8, margin = margin(b = 10) ),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+
  labs(title = "(A4) Plasma")
plot.D1_Plas

plot.A2_Plas <- plot.A2_Plas + theme(axis.title.y = element_blank())
plot.B1_Plas <- plot.B1_Plas + theme(axis.title.y = element_blank())
plot.D1_Plas <- plot.D1_Plas + theme(axis.title.y = element_blank())
### Save figure #####
ggsave("Figure MCMC-Male.tiff",scale = 1,
       plot = grid.arrange(plot.A1_Plas, plot.A2_Plas, plot.B1_Plas, plot.D1_Plas, ncol = 4),
       path = "D:/zs/PBPK/2025GenX/Artwork",
       width = 80, height = 20, units = "cm",dpi=320)
warnings()
###############################################################################################################################
##Evaluation                                                                                                                  #
#Mice                                                                                                                         #
#oral daily dose to 0.0004, 0.002, 0.0032, 0.01, 0.32, 0.4, 1, 2, 10, 100 mg/kg-d for 1, 5, 14, 28,  29, 84, 126 and 140 days #
###############################################################################################################################
## Read the data and later used in model calibration and evaluation
evaluation <- read.csv(file = "Evalution-M.csv")

OBS.A_Pla <- evaluation %>% filter(Study == 1 & Sample == "Plasma") %>% select(Time = "Time", Plasma = "Conc", SD = "SD", Study = "Study", Dose = "Dose") # TDoses = 126
OBS.B_CL  <- evaluation %>% filter(Study == 2 & Sample == "Liver" ) %>% select(Time = "Time", CL     = "Conc", SD = "SD", Study = "Study", Dose = "Dose") # TDoses = 28
OBS.C_Pla <- evaluation %>% filter(Study == 3 & Sample == "Serum" ) %>% select(Time = "Time", Plasma = "Conc", SD = "SD", Study = "Study", Dose = "Dose") # TDoses = 140
OBS.C_CL  <- evaluation %>% filter(Study == 3 & Sample == "Liver" ) %>% select(Time = "Time", CL     = "Conc", SD = "SD", Study = "Study", Dose = "Dose") # TDoses = 140
OBS.D_Pla <- evaluation %>% filter(Study == 4 & Sample == "Serum" ) %>% select(Time = "Time", Plasma = "Conc", SD = "SD", Study = "Study", Dose = "Dose") # TDoses = 28
OBS.D_CL  <- evaluation %>% filter(Study == 4 & Sample == "Liver" ) %>% select(Time = "Time", CL     = "Conc", SD = "SD", Study = "Study", Dose = "Dose") # TDoses = 28
OBS.E_Pla <- evaluation %>% filter(Study == 5 & Sample == "Plasma") %>% select(Time = "Time", Plasma = "Conc", SD = "SD", Study = "Study", Dose = "Dose") # TDoses = 84
OBS.F_Pla <- evaluation %>% filter(Study == 6 & Sample == "Serum" ) %>% select(Time = "Time", Plasma = "Conc", SD = "SD", Study = "Study", Dose = "Dose") # TDoses = 28

pars.eva <- log(c(
  # Physiological parameters
  BW                             = 0.025,
  QCC                            = 16.5,
  QLC                            = 0.161,
  QLuC                           = 0.005,
  QKC                            = 0.091,
  Htc                            = 0.48,
  VPlasC                         = 0.049,
  VLC                            = 0.055,
  VLuC                           = 0.007,
  VKC                            = 0.017,
  VFilC                          = 0.0017,
  FVBK                           = 0.160,
  GFRC                           = 62.1,
  GEC                            = 0.54,
  
  # Chemical-specific parameters (final mean values)
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

pred.eva <- function(pars) { ## pars: input parameters, Dose: input dose, Dose regimen: 0.0004, 0.002, 0.0032, 0.01, 0.32, 0.4, 1, 2, 10, 100 mg/kg/day
  
  ## Get out of log domain
  pars <- exp(pars)                   ## Return a list of exp from log scale
  
  ## Exposure scenario 
  BW          = 0.025                  ## Body weight 
  tinterval   = 24                     ## Time interval; 
  
  ##1
  TDOSE1    = 126                    ## Total dosing/Dose times
  DOSE1     = 0.0032                 ## Input oral dose
  DOSEoral1 = DOSE1*BW               ## Amount of oral dose
  ex.oral1  <- ev(ID = 1, amt = DOSEoral1, ii = tinterval, addl = TDOSE1 - 1, cmt = "AST", replicate = FALSE)
  tsamp1    = tgrid(0, tinterval*(TDOSE1 - 1) + 24*7, 1) ## set up the output time 
  
  ## Simulation of exposure scenaior 
  out1 <- 
    mod %>%          # model object
    param    (pars) %>%
    Req      (Plasma)%>%
    update   (atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  ## Atol: Absolute tolerance parameter; maxsteps:maximum number of steps; mindt: simulation output time below which there model will assume to have not advanced          
    mrgsim_d (data = ex.oral1, tgrid = tsamp1)   
  
  ## Extract the "Time", "Plasma", "CL" from out; 
  outdf1 <- cbind.data.frame(Time   = out1$time/24, Plasma = out1$Plasma) 
  
  ##2
  TDOSE2    = 28                     
  DOSE2     = 100                    
  DOSEoral2 = DOSE2*BW               
  ex.oral2  <- ev(ID = 1, amt = DOSEoral2, ii = tinterval, addl = TDOSE2 - 1, cmt = "AST", replicate = FALSE)
  tsamp2    = tgrid(0, tinterval*(TDOSE2 - 1) + 24*7, 1) 
  
  out2 <- 
    mod %>%          
    param (pars) %>%
    Req   (Liver)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral2, tgrid = tsamp2)   
  
  outdf2 <- cbind.data.frame(Time = out2$time/24, CL = out2$Liver) 
  
  ##3
  TDOSE3      = 140                    
  tsamp3      = tgrid(0, tinterval*(TDOSE3 - 1) + 24*7, 1)
  
  #a
  DOSE3a      = 0.0004                 
  DOSEoral3a  = DOSE3a*BW              
  ex.oral3a   <- ev(ID = 1, amt = DOSEoral3a, ii = tinterval, addl = TDOSE3 - 1, cmt = "AST", replicate = FALSE)
  out3a <- 
    mod %>%          
    param (pars) %>%
    Req   (Plasma, Liver)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral3a, tgrid = tsamp3)   
  
  outdf3a <- cbind.data.frame(Time   = out3a$time/24, 
                              Plasma = out3a$Plasma,
                              CL     = out3a$Liver) 
  #b
  DOSE3b      = 0.002                 
  DOSEoral3b  = DOSE3b*BW              
  ex.oral3b   <- ev(ID = 1, amt = DOSEoral3b, ii = tinterval, addl = TDOSE3 - 1, cmt = "AST", replicate = FALSE)
  out3b <- 
    mod %>%          
    param (pars) %>%
    Req   (Plasma, Liver)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral3b, tgrid = tsamp3)   
  
  outdf3b <- cbind.data.frame(Time   = out3b$time/24, 
                              Plasma = out3b$Plasma,
                              CL     = out3b$Liver)
  #c
  DOSE3c      = 0.01                 
  DOSEoral3c  = DOSE3c*BW              
  ex.oral3c   <- ev(ID = 1, amt = DOSEoral3c, ii = tinterval, addl = TDOSE3 - 1, cmt = "AST", replicate = FALSE)
  out3c <- 
    mod %>%          
    param (pars) %>%
    Req   (Plasma, Liver)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral3c, tgrid = tsamp3)   
  
  outdf3c <- cbind.data.frame(Time   = out3c$time/24, 
                              Plasma = out3c$Plasma,
                              CL     = out3c$Liver)
  
  ##4
  TDOSE4      = 28                     
  tsamp4      = tgrid(0, tinterval*(TDOSE4 - 1) + 24*7, 1)
  
  #a
  DOSE4a      = 0.4                 
  DOSEoral4a  = DOSE4a*BW              
  ex.oral4a   <- ev(ID = 1, amt = DOSEoral4a, ii = tinterval, addl = TDOSE4 - 1, cmt = "AST", replicate = FALSE)
  out4a <- 
    mod %>%          
    param (pars) %>%
    Req   (Plasma, Liver)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral4a, tgrid = tsamp4)   
  
  outdf4a <- cbind.data.frame(Time   = out4a$time/24, 
                              Plasma = out4a$Plasma,
                              CL     = out4a$Liver) 
  #b
  DOSE4b      = 2                 
  DOSEoral4b  = DOSE4b*BW              
  ex.oral4b   <- ev(ID = 1, amt = DOSEoral4b, ii = tinterval, addl = TDOSE4 - 1, cmt = "AST", replicate = FALSE)
  out4b <- 
    mod %>%          
    param (pars) %>%
    Req   (Plasma, Liver)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral4b, tgrid = tsamp4)   
  
  outdf4b <- cbind.data.frame(Time   = out4b$time/24, 
                              Plasma = out4b$Plasma,
                              CL     = out4b$Liver)
  #c
  DOSE4c      = 10                 
  DOSEoral4c  = DOSE4c*BW              
  ex.oral4c   <- ev(ID = 1, amt = DOSEoral4c, ii = tinterval, addl = TDOSE4 - 1, cmt = "AST", replicate = FALSE)
  out4c <- 
    mod %>%          
    param (pars) %>%
    Req   (Plasma, Liver)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral4c, tgrid = tsamp4)   
  
  outdf4c <- cbind.data.frame(Time   = out4c$time/24, 
                              Plasma = out4c$Plasma,
                              CL     = out4c$Liver)
  
  ##5
  TDOSE5      = 84                      ## Total dosing/Dose times
  DOSE5       = 0.32                   ## Input oral dose
  DOSEoral5   = DOSE5*BW               ## Amount of oral dose
  ex.oral5    <- ev(ID = 1, amt = DOSEoral5, ii = tinterval, addl = TDOSE5 - 1, cmt = "AST", replicate = FALSE)
  tsamp5      = tgrid(0, tinterval*(TDOSE5 - 1) + 24*7, 1) 
  
  out5 <- 
    mod %>%          
    param (pars) %>%
    Req   (Plasma)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral5, tgrid = tsamp5)   
  
  outdf5 <- cbind.data.frame(Time = out5$time/24, Plasma = out5$Plasma) 
  
  ##6
  TDOSE6      = 28
  tsamp6      = tgrid(0, tinterval*(TDOSE6 - 1) + 24*7, 1)
  
  #a
  DOSE6a      = 1                 
  DOSEoral6a  = DOSE6a*BW              
  ex.oral6a   <- ev(ID = 1, amt = DOSEoral6a, ii = tinterval, addl = TDOSE6 - 1, cmt = "AST", replicate = FALSE)
  out6a <- 
    mod %>%          
    param (pars) %>%
    Req   (Plasma)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral6a, tgrid = tsamp6)   
  
  outdf6a <- cbind.data.frame(Time   = out6a$time/24, 
                              Plasma = out6a$Plasma) 
  #b
  DOSE6b      = 10                 
  DOSEoral6b  = DOSE6b*BW              
  ex.oral6b   <- ev(ID = 1, amt = DOSEoral6b, ii = tinterval, addl = TDOSE6 - 1, cmt = "AST", replicate = FALSE)
  out6b <- 
    mod %>%          
    param (pars) %>%
    Req   (Plasma)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral6b, tgrid = tsamp6)   
  
  outdf6b <- cbind.data.frame(Time   = out6b$time/24, 
                              Plasma = out6b$Plasma) 
  #c
  DOSE6c      = 100                 
  DOSEoral6c  = DOSE6c*BW              
  ex.oral6c   <- ev(ID = 1, amt = DOSEoral6c, ii = tinterval, addl = TDOSE6 - 1, cmt = "AST", replicate = FALSE)
  out6c <- 
    mod %>%          
    param (pars) %>%
    Req   (Plasma)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral6c, tgrid = tsamp6)   
  
  outdf6c <- cbind.data.frame(Time   = out6c$time/24, 
                              Plasma = out6c$Plasma) 
  
  return (list("outdf1"   = outdf1, 
               "outdf2"   = outdf2,
               "outdf3a"  = outdf3a, "outdf3b"  = outdf3b, "outdf3c"  = outdf3c,
               "outdf4a"  = outdf4a, "outdf4b"  = outdf4b, "outdf4c"  = outdf4c,
               "outdf5"   = outdf5,
               "outdf6a"  = outdf6a, "outdf6b"  = outdf6b, "outdf6c"  = outdf6c)) 
}

out_1  <- pred.eva(pars.eva)[[1]]  %>% mutate(SD = 0, Study = 1, Dose = 0.0032)%>% filter(Time > 0)
out_2  <- pred.eva(pars.eva)[[2]]  %>% mutate(SD = 0, Study = 2, Dose = 100   )%>% filter(Time > 0)
out_3a <- pred.eva(pars.eva)[[3]]  %>% mutate(SD = 0, Study = 3, Dose = 0.0004)%>% filter(Time > 0)
out_3b <- pred.eva(pars.eva)[[4]]  %>% mutate(SD = 0, Study = 3, Dose = 0.002 )%>% filter(Time > 0)
out_3c <- pred.eva(pars.eva)[[5]]  %>% mutate(SD = 0, Study = 3, Dose = 0.01  )%>% filter(Time > 0)
out_3  <- rbind.data.frame (out_3a, out_3b, out_3c)
out_4a <- pred.eva(pars.eva)[[6]]  %>% mutate(SD = 0, Study = 4, Dose = 0.4   )%>% filter(Time > 0)
out_4b <- pred.eva(pars.eva)[[7]]  %>% mutate(SD = 0, Study = 4, Dose = 2     )%>% filter(Time > 0)
out_4c <- pred.eva(pars.eva)[[8]]  %>% mutate(SD = 0, Study = 4, Dose = 10    )%>% filter(Time > 0)
out_4  <- rbind.data.frame (out_4a, out_4b, out_4c)
out_5  <- pred.eva(pars.eva)[[9]]  %>% mutate(SD = 0, Study = 5, Dose = 0.32  )%>% filter(Time > 0)
out_6a <- pred.eva(pars.eva)[[10]] %>% mutate(SD = 0, Study = 6, Dose = 1     )%>% filter(Time > 0)
out_6b <- pred.eva(pars.eva)[[11]] %>% mutate(SD = 0, Study = 6, Dose = 10    )%>% filter(Time > 0)
out_6c <- pred.eva(pars.eva)[[12]] %>% mutate(SD = 0, Study = 6, Dose = 100   )%>% filter(Time > 0)
out_6  <- rbind.data.frame (out_6a, out_6b, out_6c)

out_1_Pla <- out_1 %>% filter(Study == 1, Time == 126)
out_2_CL  <- out_2 %>% filter(Study == 2, Time == 29)
out_3_Pla <- out_3 %>% filter(Study == 3, Time == 140)%>% select(Time, Plasma, SD, Study, Dose)
out_3_CL  <- out_3 %>% filter(Study == 3, Time == 140)%>% select(Time, CL, SD, Study, Dose)
out_4_Pla <- out_4 %>% filter(Study == 4, Time == 28 )%>% select(Time, Plasma, SD, Study, Dose)
out_4_CL  <- out_4 %>% filter(Study == 4, Time == 28 )%>% select(Time, CL, SD, Study, Dose)
out_5_Pla <- out_5 %>% filter(Study == 5, Time == 84)
out_6_Pla <- out_6 %>% filter(Study == 6, Time %in% c(1, 5, 14, 28))

out_Pla <- rbind.data.frame (out_1_Pla, out_3_Pla, out_4_Pla, out_5_Pla, out_6_Pla)
out_Pla <- out_Pla %>% mutate(Matrix = c("Pre.Plasma"))
out_CL  <- rbind.data.frame (out_2_CL, out_3_CL, out_4_CL)
out_CL  <- out_CL  %>% mutate(Matrix = c("Pre.Liver"))

OBS_Pla <- rbind.data.frame (OBS.A_Pla, OBS.C_Pla, OBS.D_Pla, OBS.E_Pla, OBS.F_Pla)
OBS_Pla <- OBS_Pla %>% mutate(Matrix = c("Obs.Plasma"))
OBS_CL  <- rbind.data.frame (OBS.B_CL, OBS.C_CL, OBS.D_CL)
OBS_CL  <- OBS_CL  %>% mutate(Matrix = c("Obs.Liver"))

############################ Pre and Obs comparsion Plot Liver #####################################################################################################################################
combined_CL_M  <- rbind(out_CL , OBS_CL )
combined_CL1<- combined_CL %>% filter(Dose %in% c(0.4,2,10,100))
p_CL1<- 
  ggplot(combined_CL1, aes(x = as.factor(Dose), y = CL, fill = as.factor(Matrix))) + 
  geom_bar(stat = "identity", color = "black", size = 1.2, position = position_dodge(width = 0.8), width = 0.8) +
  geom_errorbar(aes(ymin = CL, ymax = CL+SD), width = 0.2, position = position_dodge(0.8),  size = 0.9) + 
  scale_fill_manual(values=c("#BFD0E1", "gray"))+
  labs(title = "(A2) Liver", x = "Dose (mg/kg/day)", y = "GenX Concentration in Liver (μg/g)", fill = "Legend") 

p_CL1 = p_CL1 + 
  scale_y_continuous(limits = c(-0.8, 40),breaks=c(10,20,30),expand = c(0,0))+
  theme_classic() +  
  theme(axis.text.x = element_text(size = 22),   
        axis.text.y = element_text(size = 22),   
        axis.title.x = element_text(face = "bold", size = 25), 
        axis.title.y = element_text(face = "bold", size = 25),
        legend.text  = element_text(size = 16),
        legend.title = element_text(size = 16),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.06, vjust = -8, margin = margin(b = 10) ), 
        legend.position = "none",  
        legend.justification = c(0, 1),   
        legend.background = element_rect(fill = "white", color = "black", size = 1))
p_CL1

combined_CL2<- combined_CL  %>% filter(Dose %in% c(0.002,0.01))
combined_CL2<- combined_CL2 %>% mutate(Dose = formatC(Dose, format = "f", digits = 4))
p_CL2<- 
  ggplot(combined_CL2, aes(x = as.factor(Dose), y = CL, fill = as.factor(Matrix))) + 
  geom_bar(stat = "identity", color = "black", size = 1.2, position = position_dodge(width = 0.6), width = 0.6) +
  geom_errorbar(aes(ymin = CL, ymax = CL+SD), width = 0.2, position = position_dodge(0.6),  size = 0.9) + 
  scale_fill_manual(values=c("#BFD0E1", "gray"))+##"#BFD0E1", "#7FB2D5"
  labs(title = "(A1) Liver", x = "Dose (mg/kg/day)", y = "GenX Concentration in Liver (μg/g)", fill = "Legend") 

p_CL2 = p_CL2 + 
  scale_y_continuous(limits = c(-0.001, 0.05),breaks=c(0.015,0.03,0.045),expand = c(0,0))+
  theme_classic() +  
  theme(axis.text.x = element_text(size = 22),   
        axis.text.y = element_text(size = 22),   
        axis.title.x = element_text(face = "bold", size = 25), 
        axis.title.y = element_text(face = "bold", size = 25),
        legend.text  = element_text(size = 16),
        legend.title = element_text(size = 16),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.06, vjust = -8, margin = margin(b = 10) ), 
        legend.position = "none",  
        legend.justification = c(0, 1),   
        legend.background = element_rect(fill = "white", color = "black", size = 1))
p_CL2
p_CL1 <- p_CL1 + theme(axis.title.y = element_blank())
ggsave("Figure 3-Male-liver.tiff",scale = 1,
       plot = grid.arrange(p_CL2, p_CL1 , ncol = 2),
       path = "D:/zs/PBPK/2025GenX/Artwork",
       width = 40, height = 20, units = "cm",dpi=320)

############################ Pre and Obs comparsion Plot Plasma #####################################################################################################################
combined_Pla_M <- rbind(out_Pla, OBS_Pla)
combined_Pla1<- combined_Pla %>% filter(Study==6)%>% filter(Dose==1)
p_Pla1<- 
  ggplot(combined_Pla1, aes(x = as.factor(Time), y = Plasma, fill = as.factor(Matrix))) + # do the histogram plot by the factor of matrix, so you have a bar for pre.plasma, a bar for obs.plasma, so on.
  geom_bar(stat = "identity", color = "black", size = 1.2, position = position_dodge(width = 0.8), width = 0.8) +
  geom_errorbar(aes(ymin = Plasma, ymax = Plasma+SD), width = 0.2, position = position_dodge(0.8),  size = 0.9) + 
  scale_fill_manual(values=c("#BFD0E1","gray"))+##"#F6B3AC","#F47F72"
  labs(title = "(A1) Plasma, Dose=1 mg/kg/day", x = "Time (day)", y = "GenX Concentration in Plasma (μg/ml)", fill = "Legend") 

p_Pla1 = p_Pla1 + 
  scale_y_continuous(limits = c(-0.06, 3),breaks=c(0.5,1,1.5,2,2.5),expand = c(0,0))+
  theme_classic() +  
  theme(axis.text.x = element_text(size = 22),   
        axis.text.y = element_text(size = 22),   
        axis.title.x = element_text(face = "bold", size = 25), 
        axis.title.y = element_text(face = "bold", size = 25),
        legend.text  = element_text(size = 16),
        legend.title = element_text(size = 16),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.06, vjust = -8, margin = margin(b = 10) ), 
        legend.position = "none",  
        legend.justification = c(0, 1),   
        legend.background = element_rect(fill = "white", color = "black", size = 1))
p_Pla1

combined_Pla2<- combined_Pla %>% filter(Study==6)%>% filter(Dose==10)
p_Pla2<- 
  ggplot(combined_Pla2, aes(x = as.factor(Time), y = Plasma, fill = as.factor(Matrix))) + # do the histogram plot by the factor of matrix, so you have a bar for pre.plasma, a bar for obs.plasma, so on.
  geom_bar(stat = "identity", color = "black", size = 1.2, position = position_dodge(width = 0.8), width = 0.8) +
  geom_errorbar(aes(ymin = Plasma, ymax = Plasma+SD), width = 0.2, position = position_dodge(0.8),  size = 0.9) + 
  scale_fill_manual(values=c("#BFD0E1","gray"))+
  labs(title = "(A2) Plasma, Dose=10 mg/kg/day", x = "Time (day)", y = "GenX Concentration in Plasma (μg/ml)", fill = "Legend") 

p_Pla2 = p_Pla2 + 
  scale_y_continuous(limits = c(-0.6, 30),breaks=c(5,10,15,20,25),expand = c(0,0))+
  theme_classic() +  
  theme(axis.text.x = element_text(size = 22),   
        axis.text.y = element_text(size = 22),   
        axis.title.x = element_text(face = "bold", size = 25), 
        axis.title.y = element_text(face = "bold", size = 25),
        legend.text  = element_text(size = 16),
        legend.title = element_text(size = 16),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.06, vjust = -8, margin = margin(b = 10) ), 
        legend.position = "none",  
        legend.justification = c(0, 1),   
        legend.background = element_rect(fill = "white", color = "black", size = 1))
p_Pla2

combined_Pla3<- combined_Pla %>% filter(Study==6)%>% filter(Dose==100)
p_Pla3<- 
  ggplot(combined_Pla3, aes(x = as.factor(Time), y = Plasma, fill = as.factor(Matrix))) + # do the histogram plot by the factor of matrix, so you have a bar for pre.plasma, a bar for obs.plasma, so on.
  geom_bar(stat = "identity", color = "black", size = 1.2, position = position_dodge(width = 0.8), width = 0.8) +
  geom_errorbar(aes(ymin = Plasma, ymax = Plasma+SD), width = 0.2, position = position_dodge(0.8),  size = 0.9) + 
  scale_fill_manual(values=c("#BFD0E1","gray"))+
  labs(title = "(A3) Plasma, Dose=100 mg/kg/day", x = "Time (day)", y = "GenX Concentration in Plasma (μg/ml)", fill = "Legend") 

p_Pla3 = p_Pla3 + 
  scale_y_continuous(limits = c(-4.8, 240),breaks=c(40,80,120,160,200),expand = c(0,0))+
  theme_classic() +  
  theme(axis.text.x = element_text(size = 22),   
        axis.text.y = element_text(size = 22),   
        axis.title.x = element_text(face = "bold", size = 25), 
        axis.title.y = element_text(face = "bold", size = 25),
        legend.text  = element_text(size = 16),
        legend.title = element_text(size = 16),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.06, vjust = -8, margin = margin(b = 10) ), 
        legend.position = "none",  
        legend.justification = c(0, 1),   
        legend.background = element_rect(fill = "white", color = "black", size = 1))
p_Pla3

p_Pla2 <- p_Pla2 + theme(axis.title.y = element_blank())
p_Pla3 <- p_Pla3 + theme(axis.title.y = element_blank())

combined_Pla4<- combined_Pla %>% filter(Dose %in% c(0.0004,0.002,0.0032))%>% mutate(Dose = formatC(Dose, format = "f", digits = 4))
p_Pla4<- 
  ggplot(combined_Pla4, aes(x = as.factor(Dose), y = Plasma, fill = as.factor(Matrix))) + # do the histogram plot by the factor of matrix, so you have a bar for pre.plasma, a bar for obs.plasma, so on.
  geom_bar(stat = "identity", color = "black", size = 1.2, position = position_dodge(width = 0.8), width = 0.8) +
  geom_errorbar(aes(ymin = Plasma, ymax = Plasma+SD), width = 0.2, position = position_dodge(0.8),  size = 0.9) + 
  scale_fill_manual(values=c("#BFD0E1","gray"))+
  labs(title = "(A4) Plasma", x = "Dose (mg/kg/day)", y = "GenX Concentration in Plasma (μg/ml)", fill = "Legend") 

p_Pla4 = p_Pla4 + 
  scale_y_continuous(limits = c(-0.0008, 0.04),breaks=c(0.01,0.02,0.03),expand = c(0,0))+
  theme_classic() +  
  theme(axis.text.x = element_text(size = 22),   
        axis.text.y = element_text(size = 22),   
        axis.title.x = element_text(face = "bold", size = 25), 
        axis.title.y = element_text(face = "bold", size = 25),
        legend.text  = element_text(size = 16),
        legend.title = element_text(size = 16),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.06, vjust = -8, margin = margin(b = 10) ), 
        legend.position = "none",  
        legend.justification = c(0, 1),   
        legend.background = element_rect(fill = "white", color = "black", size = 1))
p_Pla4

combined_Pla5<- combined_Pla %>% filter(Dose %in% c(0.01,0.32))
p_Pla5<- 
  ggplot(combined_Pla5, aes(x = as.factor(Dose), y = Plasma, fill = as.factor(Matrix))) + # do the histogram plot by the factor of matrix, so you have a bar for pre.plasma, a bar for obs.plasma, so on.
  geom_bar(stat = "identity", color = "black", size = 1.2, position = position_dodge(width = 0.5), width = 0.5) +
  geom_errorbar(aes(ymin = Plasma, ymax = Plasma+SD), width = 0.2, position = position_dodge(0.5),  size = 0.9) + 
  scale_fill_manual(values=c("#BFD0E1","gray"))+
  labs(title = "(A5) Plasma", x = "Dose (mg/kg/day)", y = "GenX Concentration in Plasma (μg/ml)", fill = "Legend") 

p_Pla5 = p_Pla5 + 
  scale_y_continuous(limits = c(-0.014, 0.7),breaks=c(0.15,0.3,0.45,0.6),expand = c(0,0))+
  theme_classic() +  
  theme(axis.text.x = element_text(size = 22),   
        axis.text.y = element_text(size = 22),   
        axis.title.x = element_text(face = "bold", size = 25), 
        axis.title.y = element_text(face = "bold", size = 25),
        legend.text  = element_text(size = 16),
        legend.title = element_text(size = 16),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.06, vjust = -8, margin = margin(b = 10) ), 
        legend.position = "none",  
        legend.justification = c(0, 1),   
        legend.background = element_rect(fill = "white", color = "black", size = 1))
p_Pla5
p_Pla5 <- p_Pla5 + theme(axis.title.y = element_blank())

combined_Pla6<- combined_Pla %>% filter(Study %in% c(4,5)) %>% filter(Dose %in% c(0.4,2,10))
p_Pla6<- 
  ggplot(combined_Pla6, aes(x = as.factor(Dose), y = Plasma, fill = as.factor(Matrix))) + # do the histogram plot by the factor of matrix, so you have a bar for pre.plasma, a bar for obs.plasma, so on.
  geom_bar(stat = "identity", color = "black", size = 1.2, position = position_dodge(width = 0.8), width = 0.8) +
  geom_errorbar(aes(ymin = Plasma, ymax = Plasma+SD), width = 0.2, position = position_dodge(0.8),  size = 0.9) + 
  scale_fill_manual(values=c("#BFD0E1","gray"))+
  labs(title = "(A6) Plasma", x = "Dose (mg/kg/day)", y = "GenX Concentration in Plasma (μg/ml)", fill = "Legend") 

p_Pla6 = p_Pla6 + 
  scale_y_continuous(limits = c(-1.3, 65),breaks=c(15,30,45,60),expand = c(0,0))+
  theme_classic() +  
  theme(axis.text.x = element_text(size = 22),   
        axis.text.y = element_text(size = 22),   
        axis.title.x = element_text(face = "bold", size = 25), 
        axis.title.y = element_text(face = "bold", size = 25),
        legend.text  = element_text(size = 16),
        legend.title = element_text(size = 16),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.06, vjust = -8, margin = margin(b = 10) ), 
        legend.position = "none",  
        legend.justification = c(0, 1),   
        legend.background = element_rect(fill = "white", color = "black", size = 1))
p_Pla6
p_Pla6 <- p_Pla6 + theme(axis.title.y = element_blank())

ggsave("Figure 3-Male.tiff",scale = 1,
       plot = grid.arrange(p_Pla1, p_Pla2, p_Pla3, p_Pla4, p_Pla5, p_Pla6, ncol = 3, nrow = 2),##, p_CL2, p_CL1
       path = "D:/zs/PBPK/2025GenX/Artwork",
       width = 60, height = 40, units = "cm",dpi=320)

write.csv(combined_CL , file = "combined_CL.csv" , row.names = FALSE)
write.csv(combined_Pla, file = "combined_Pla.csv", row.names = FALSE)

################################################################################################################
eva.outcome <- read.csv(file = "evadata.csv")
eva.outcome %<>% mutate (Log.OBS = log(OBS,10), Log.PRE = log(PRE,10))
fit <- lm(Log.OBS ~Log.PRE, data = eva.outcome)
summary(fit)
adjusted_r_squared <- round(summary(fit)$adj.r.squared, digits = 2)
label_text <- paste("italic(R)^{2} == ", adjusted_r_squared)
#######################################################################################
eva.outcome %<>% mutate(res = residuals(fit), 
                        prediction = predict(fit), 
                        OPR = PRE/OBS,            ## OPR: the ratio of prediction value and observed data
                        log.OPR =  log(OPR,10)) 
p <- 
  ggplot(eva.outcome, aes(Log.OBS, Log.PRE)) + ## using log-sacle axis
  geom_abline (intercept = 0, 
               slope     = 1,
               color     ="black", linetype = "dashed",linewidth = 1, alpha = 0.8) +
  geom_point  (aes(shape   = as.factor(Sample), color = as.factor(Sample)),size = 3)  +
  annotation_logticks() +
  scale_y_continuous(limits = c(-8,8), labels = scales::math_format(10^.x))+
  scale_x_continuous(limits = c(-8,8),labels = scales::math_format(10^.x)) +
  scale_color_manual(values = c("Liver" = "#7ac7e2", "Plasma" = "#e3716e")) +
  scale_shape_manual(values = c("Liver" = 16, "Plasma" = 17))+
  coord_cartesian(clip = "off")

p.eva <- p + theme_bw(base_family = "Times New Roman") +
  theme (
    plot.background         = element_rect (fill="White"),
    text                    = element_text (family = "Times New Roman"),   # text front (Time new roman)
    panel.border            = element_rect (colour = "black", fill=NA, linewidth =2),
    panel.background        = element_rect (fill="White"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.text               = element_text (size   = 28, colour = "black", face = "bold"),    # tick labels along axes 
    axis.title              = element_text (size   = 30, colour = "black", face = "bold"),   # label of axes
    legend.position         =c(0.8, 0.2),
    legend.background = element_rect(fill = "white", color = "white"),
    legend.title = element_text(size = 30, face = "bold"),  
    legend.text = element_text(size = 28) ) +
  labs (title = "", x = "Observed values (ug/ml or μg/g)",  y = "Predicted values (ug/ml or μg/g)", shape = "Sample", color = "Sample")+
  annotate("text", x = -Inf, y = Inf, label = "(A)", hjust = 1.3, vjust = 0.8, size = 15, family = "Times New Roman", colour = "black") +
  annotate("text", x = -3, y = 5, label = label_text , parse = TRUE, size = 15, color = "black", family = "Times New Roman")
print(p.eva)
ggsave("Figure S3-male.tiff",scale = 1,
       plot = p.eva,
       path = "D:/zs/PBPK/2025GenX/Artwork",
       width = 24, height = 24, units = "cm",dpi=320)
