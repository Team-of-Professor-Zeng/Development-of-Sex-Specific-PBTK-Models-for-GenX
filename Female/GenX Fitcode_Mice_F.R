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

setwd("D:/zs/PBPK/1-方案/GenX/小鼠/code/Mice/Female")
getwd()
### Female
## Input mrgsolve-based PBPK Model
source (file = "GenX MMod_F.R")

## Set working direction to the data files
## Build mrgsolve-based PBPK Model
mod <- mcode ("micepbpk", MicePBPK_F.code)

## input data set for model calibration/ oral
Data_A1    <- read.csv(file = "Female.csv")
# Model calibration for PBPK model based on the data of TK study #
# A1. : Mice, oral single dose of 10 mg/kg,   matrix: Plasma, Sampling time: 0.25,0.5,1,2,4,8,12,24,48,72,96,120,144,168h Gannon et al. (2015)
# A2. : Mice, oral single dose of 30 mg/kg,   matrix: Plasma, Sampling time: 0.25,0.5,1,2,4,8,12,24,48,72,96,120,144,168h Gannon et al. (2015)
# B1. : Mice, oral 0-14d dose of 1 mg/kg,     matrix: blood, Sampling time: 72,168,336,360,408,504h  Wen et al. (2022)
# B2. : Mice, oral 0-14d dose of 1 mg/kg,     matrix: Liver, Sampling time: 72,168,336,360,408,504h  Wen et al. (2022)
# B3. : Mice, oral 0-14d dose of 1 mg/kg,     matrix: Lung, Sampling time: 72,168,336,360,408,504h  Wen et al. (2022)
# B4. : Mice, oral 0-14d dose of 1 mg/kg,     matrix: Kidney, Sampling time: 72,168,336,360,408,504h  Wen et al. (2022)
#===================================================================================================

## Read these datasets and later used in model calibration

OBS.A1  <- Data_A1 %>% filter(Study == 1 & Sample == "Plasma" & Dose == 10)   %>% select(Time = "Time", Plasma = "Conc") # TDosesA = 1
OBS.A2  <- Data_A1 %>% filter(Study == 1 & Sample == "Plasma" & Dose == 30)   %>% select(Time = "Time", Plasma = "Conc") # TDosesA = 1
OBS.B1  <- Data_A1 %>% filter(Study == 2 & Sample == "blood"  & Dose == 1 )   %>% select(Time = "Time", Plasma = "Conc") # TDosesB = 14
OBS.B2  <- Data_A1 %>% filter(Study == 2 & Sample == "Liver"  & Dose == 1 )   %>% select(Time = "Time", Liver  = "Conc") # TDosesB = 14
OBS.B3  <- Data_A1 %>% filter(Study == 2 & Sample == "Lung"   & Dose == 1 )   %>% select(Time = "Time", Lung   = "Conc") # TDosesB = 14
OBS.B4  <- Data_A1 %>% filter(Study == 2 & Sample == "Kidney" & Dose == 1 )   %>% select(Time = "Time", Kidney = "Conc") # TDosesB = 14

###ADOSE:0.2mg
## Define the prediction function (for least squres fit using levenberg-marquart algorithm)
pred.Mice <- function(pars) {
  
  ## Get out of log domain
  pars %<>% lapply(exp)                 ## return a list of exp (parameters) from log domain
  
  ## Define the exposure scenario 
  BW          = 0.02                                      # Mice body weight
  tinterval   = 24                                        # Time interval
  TDosesA     = 1                                         # Total dosing/Dose times
  TDosesB     = 14                                        # Total dosing/Dose times
  
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
  #B
  PDOSEoral.B = 1                                        # Single oral dose from Wen et al. (2022)
  DOSEoral.B  = PDOSEoral.B*BW                          # Amount of oral dose
  ex.oral.B<- ev  (ID = 1,              ## One individual
                   amt  = DOSEoral.B,  ## Amount of dose 
                   ii   = tinterval,    ## Time interval
                   addl = TDosesB-1,    ## Addtional doseing 
                   cmt  = "AST",        ## The dosing comaprtment: AST Stomach  
                   replicate = FALSE)   ## No replicate
  
  ## set up the exposure time
  tsampA=tgrid(0,tinterval*(TDosesA-1)+24*14,0.5)
  tsampB=tgrid(0,tinterval*(TDosesB-1)+24*10,0.5)
  
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
  
  out.B <- 
    mod %>%                                                 # model object
    param(pars) %>%                                         # to update the parameters in the model subject
    update(atol = 1E-19, maxsteps=500000) %>%                # solver setting, atol: Absolute tolerance parameter
    mrgsim_d(data = ex.oral.B, tgrid=tsampB)               # Set up the simulation run
  out.B<-cbind.data.frame(Time  =out.B$time, 
                          Plasma =out.B$Plasma,
                          Liver  =out.B$Liver,
                          Lung   =out.B$Lung,
                          Kidney =out.B$Kidney,
                          BAL =out.B$Balance)
  out.B <- out.B %>% filter (Time > 0)                    # filter the value at time = 0
  
  return (list("out.A1"=out.A1, "out.A2"=out.A2, "out.B" =out.B))    # Return Dataframe
}

## initial parmaeters
theta.int <- log(c(
  KurineC                        = 0.122,                       ## Urinary elimination rate
  Free                           = 0.045,                       ## Free fraction in plasma
  PL                             = 1.339,                       ## Liver /plasma partition coefficient (PC)
  PK                             = 0.854,                       ## Kidney/plasma PC
  PLu                            = 0.431,                       ## Lung  /plasma PC
  PRest                          = 0.595,                       ## Rest of body/plasma PC
  K0C                            = 1,                           ## Rate of absorption of GenX in the stomach
  Kabsc                          = 2.12,                        ## Rate of absorption of GenX in the small intestines
  KunabsC                        = 0.0265                       ## Rate of unabobded dose to appear in feces
))
result <- pred.Mice(theta.int)

## Check mass balance and unfitted curve
plot(result$out.A1$Time,result$out.A1$Plasma ,type="l",lwd=2,xlab="Time(hour)",ylab="GenX concentration in plasma")
plot(result$out.A1$Time,result$out.A1$BAL,type="l",lwd=2,xlab="Time(hour)",ylab="Mass Balance")

## Cost fuction (FME) 
## Estimate the model residual by modCost function
MCcost<-function (pars){
  outdf <- pred.Mice (pars)
  cost<- modCost  (model = outdf$out.A1, obs = OBS.A1, x ="Time" ,weight = "mean")
  cost<- modCost  (model = outdf$out.A2, obs = OBS.A2, x ="Time" ,weight = "mean",cost = cost)
  cost<- modCost  (model = outdf$out.B , obs = OBS.B1, x ="Time" ,weight = "mean",cost = cost)
  cost<- modCost  (model = outdf$out.B , obs = OBS.B2, x ="Time" ,weight = "mean",cost = cost)
  cost<- modCost  (model = outdf$out.B , obs = OBS.B3, x ="Time" ,weight = "mean",cost = cost)
  cost<- modCost  (model = outdf$out.B , obs = OBS.B4, x ="Time" ,weight = "mean",cost = cost)
  return(cost)
}

## Local sensitivity analysis

## Senstivity function (FME) 
## Check the senstive parameters in the model
SnsPlasma <- sensFun(func = MCcost, parms = theta.int, varscale = 1)
Sen       <- summary(SnsPlasma)
plot(Sen)

## Selected senstive parameters
theta <- theta.int[abs(Sen$Mean) > 1.2*mean(abs(Sen$Mean))]
theta 

## Selected parameters
theta.int <- log(c(
  #KurineC                        = 0.122,                       ## Urinary elimination rate 0.0122-1.22
  Free                           = 0.045,                       ## Free fraction in plasma 0.0045-0.45
  PL                             = 1.339,                       ## Liver /plasma partition coefficient (PC) 0.1339-13.39
  #PK                             = 0.854,                       ## Kidney/plasma PC 0.0854-8.54
  PLu                            = 0.431                       ## Lung  /plasma PC 0.0431-4.31
  #PRest                          = 0.595                       ## Rest of body/plasma PC 0.0595-5.95
  #K0C                            = 1,                           ## Rate of absorption of GenX in the stomach 0.1-10
  #Kabsc                          = 2.12,                        ## Rate of absorption of GenX in the small intestines 0.212-21.2
  #KunabsC                        = 0.0265                       ## Rate of unabobded dose to appear in feces 0.00265-0.265
))

## PBPK model fitting 
## least squres fit using levenberg-marquart (method "Marq") algorithm
Fit<- modFit(f=MCcost, p=theta.int, method ="Marq",
             control = nls.lm.control(nprint=1))

summary(Fit)                                 ## Summary of fit 
exp(Fit$par)                                 ## Get the arithmetic value out of the log domain
Cost <- MCcost(Fit$par)
#############################Global fitting analysis##################################### 
FDataA <- cbind.data.frame (name= Cost$residuals$name,
                            OBS = Cost$residuals$obs, 
                            PRE = Cost$residuals$mod)

## Transformed the predicted and obseved values using log10-sacle to do the plot
FDataA %<>% mutate (Log.OBS = log(OBS,10), Log.PRE = log(PRE,10), Species = "Mice")

## Estimating the R-squared and goodness-of-fit using linear regression model
fit <- lm(Log.OBS ~Log.PRE, data = FDataA)
summary(fit)
adjusted_r_squared <- round(summary(fit)$adj.r.squared, digits = 2)
label_text <- paste("italic(R)^{2} == ", adjusted_r_squared)

#######################################################################################
FDataA %<>% mutate(res = residuals(fit), 
                   prediction = predict(fit), 
                   OPR = PRE/OBS,            ## OPR: the ratio of prediction value and observed data
                   log.OPR =  log(OPR,10)) 
write.csv(FDataA, file = 'FDataA_F.csv')
p <- 
  ggplot(FDataA, aes(Log.OBS, Log.PRE)) + 
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", linewidth = 1, alpha = 0.8) +
  geom_point(aes(shape = as.factor(name), color = as.factor(name)), size = 3) +
  annotation_logticks() +
  scale_y_continuous(limits = c(-6, 6), breaks = c(-5, 0, 5), labels = scales::math_format(10^.x)) +
  scale_x_continuous(limits = c(-6, 6), breaks = c(-5, 0, 5),  labels = scales::math_format(10^.x)) +
  scale_color_manual(values = c("Kidney" = "steelblue4", "Lung" = "#7A70B5", "Liver" = "#FCBB44", "Plasma" = "#e3716e"), labels = c("Kidney", "Lung", "Liver", "Plasma")) +
  scale_shape_manual(values = c("Kidney" = 7, "Lung" = 17, "Liver" = 16, "Plasma" = 19), labels = c("Kidney", "Lung", "Liver", "Plasma")) +
  labs(title = "", x = "Observed values (ug/ml)", y = "Predicted values (ug/ml)", color = "Sample", shape = "Sample") +
  guides(color = guide_legend(title = "Sample"), shape = guide_legend(title = "Sample"))+
  coord_cartesian(clip = "off")
p1 <- p + 
  theme(
    plot.background = element_rect(fill = "White"),
    text = element_text(family = "Times"),                                      # text front (Time new roman)
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.background = element_rect(fill = "White"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.text = element_text(size = 28, colour = "black", face = "bold"),       # tick labels along axes 
    axis.title = element_text(size = 30, colour = "black", face = "bold"),      # label of axes
    legend.position = c(0.8, 0.2),
    legend.background = element_rect(fill = "white", color = "white"),
    legend.title = element_text(size = 30, face = "bold"),
    legend.text = element_text(size = 28)
  ) +
  annotate("text", x = -Inf, y = Inf, label = "(B)", hjust = 1.3, vjust = 0.8, size = 15, family = "Times", colour = "black") +
  annotate("text", x = -3, y = 4, label = label_text , parse = TRUE, size = 15, color = "black", family = "Times")
print(p1)

ggsave("Figure 2-Female.tiff", plot = p1,
       path = "D:/zs/PBPK/2025GenX/Artwork",
       width = 24,   height = 24,  units = "cm", dpi = 320)

pct_2e <- mean(FDataA$OPR >= 1/2 & FDataA$OPR <= 2) * 100
pct_3e <- mean(FDataA$OPR >= 1/3 & FDataA$OPR <= 3) * 100
#pct_5e <- mean(FDataA$OPR >= 1/5 & FDataA$OPR <= 5) * 100
#pct_10e <- mean(FDataA$OPR >= 1/10 & FDataA$OPR <= 10) * 100
p2 <- ggplot(FDataA, aes(x = OPR)) +
  geom_histogram(aes(y = after_stat(count / sum(count) * 100)), 
                 bins = 30, fill = "steelblue", color = "black", alpha = 0.7) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  coord_cartesian(ylim = c(0,  35)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = c(0,0)) +
  coord_cartesian(clip = "off")+
  annotate("text", x = 0.005, y = 34, label = "(B)", 
           hjust = 2.1, vjust = 0.8, size = 15, family = "Times", colour = "black") +
  annotate("text", x = 0.005, y = 32, label = paste0("% 2e: ", round(pct_2e, 1), "%"),
           hjust = 0, vjust = 1, size = 9, color = "black", fontface = "bold") +
  annotate("text", x = 0.005, y = 30, label = paste0("% 3e: ", round(pct_3e, 1), "%"),
           hjust = 0, vjust = 1, size = 9, color = "grey", fontface = "bold") + 
  #annotate("text", x = 0.005, y = 28, label = paste0("% 5e: ", round(pct_5e, 1), "%"),
           #hjust = 0, vjust = 1, size = 9, color = "black", fontface = "bold") +
  #annotate("text", x = 0.005, y = 26, label = paste0("% 10e: ", round(pct_10e, 1), "%"),
           #hjust = 0, vjust = 1, size = 9, color = "black", fontface = "bold") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#e3716e", linewidth = 1) + 
  geom_vline(xintercept = c(1/2, 2), linetype = "dashed", color = "black", linewidth = 0.8) + 
  geom_vline(xintercept = c(1/3, 3), linetype = "dashed", color = "grey", linewidth = 0.8) + 
  labs(x = "Predicted / Observed", y = "Proportion (%)") +
  theme_bw(base_family = "Times") +
  theme(
    plot.background = element_rect(fill = "white"),
    text = element_text(family = "Times"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.text = element_text(size = 28, colour = "black", face = "bold"),
    axis.title = element_text(size = 30, colour = "black", face = "bold"),
    plot.margin = margin(0.8, 0.6, 0.3, 0.5, "cm"))
print(p2)

ggsave("Figure 2-Female1.tiff", p2,
       path = "D:/zs/PBPK/2025GenX/Artwork",
       width = 24,   height = 24,  units = "cm", dpi = 320)

grid.arrange(p1,p2,ncol = 2)
##############################################################################################################################################
p2 <-
  ggplot(FDataA, aes(Log.PRE, log.OPR)) +
  geom_hline(yintercept = log10(2),  linetype = 3, colour = "black", linewidth = 1) +
  geom_hline(yintercept = log10(0.5), linetype = 3, colour = "black", linewidth = 1) +
  geom_hline(yintercept = log10(3),  linetype = 3, colour = "grey50", linewidth = 1) +
  geom_hline(yintercept = log10(0.33), linetype = 3, colour = "grey50", linewidth = 1) +
  annotation_logticks() +
  scale_y_continuous(limits = c(-4,4), labels = scales::math_format(10^.x), breaks = c(-4, -2, 0, 2, 4))+
  scale_x_continuous(limits = c(-6,6),labels = scales::math_format(10^.x))+ 
  annotate("text", x = -Inf, y = Inf, label = "Female (B)", hjust = -0.25, vjust = 1.5, size = 15, family = "Times", colour = "black") +
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

## Model calibration plot using ggplot2 
Sim.fit.A = pred.Mice(Fit$par)
df.sim.A1 <- cbind.data.frame(Time  = Sim.fit.A$out.A1$Time, Plasma   = Sim.fit.A$out.A1$Plasma)
df.sim.A2 <- cbind.data.frame(Time  = Sim.fit.A$out.A2$Time, Plasma   = Sim.fit.A$out.A2$Plasma)
df.sim.B  <- cbind.data.frame(Time  = Sim.fit.A$out.B$Time,
                              Plasma= Sim.fit.A$out.B$Plasma,
                              Liver = Sim.fit.A$out.B$Liver,
                              Lung  = Sim.fit.A$out.B$Lung,
                              Kidney= Sim.fit.A$out.B$Kidney)

## Setting an initial value
df.sim.A1 <- rbind(data.frame(Time = 0, Plasma = 1e-6),df.sim.A1)
df.sim.A2 <- rbind(data.frame(Time = 0, Plasma = 1e-6),df.sim.A2)
df.sim.B  <- rbind(data.frame(Time = 0, Plasma = 1e-6, Liver = 1e-6, Lung = 1e-6, Kidney = 1e-6),df.sim.B) 

OBS.A1.1  <- Data_A1 %>% filter(Study == 1 & Sample == "Plasma" & Dose == 10)   %>% select(Time = "Time", Plasma = "Conc", SD = "SD") # TDosesA = 1
OBS.A2.1  <- Data_A1 %>% filter(Study == 1 & Sample == "Plasma" & Dose == 30)   %>% select(Time = "Time", Plasma = "Conc", SD = "SD") # TDosesA = 1
OBS.B1.1  <- Data_A1 %>% filter(Study == 2 & Sample == "blood"  & Dose == 1 )   %>% select(Time = "Time", Plasma = "Conc", SD = "SD") # TDosesB = 14
OBS.B2.1  <- Data_A1 %>% filter(Study == 2 & Sample == "Liver"  & Dose == 1 )   %>% select(Time = "Time", Liver  = "Conc", SD = "SD") # TDosesB = 14
OBS.B3.1  <- Data_A1 %>% filter(Study == 2 & Sample == "Lung"   & Dose == 1 )   %>% select(Time = "Time", Lung   = "Conc", SD = "SD") # TDosesB = 14
OBS.B4.1  <- Data_A1 %>% filter(Study == 2 & Sample == "Kidney" & Dose == 1 )   %>% select(Time = "Time", Kidney = "Conc", SD = "SD") # TDosesB = 14

plot.A1 <- ggplot() +
  geom_line(data = df.sim.A1, aes(Time, Plasma), col = "#5A5A5A", lwd = 1.5) +
  geom_point(data = OBS.A1, aes(Time, Plasma), col = "#e3716e", size = 5) +
  geom_errorbar(data = OBS.A1.1, aes(x = Time, ymin = Plasma - SD, ymax = Plasma + SD), col = "#e3716e", width = 5, size = 1.5) +
  ylab("GenX Concentration in Plasma (μg/ml)")+
  xlab("Time (hours)") +
  xlim(c(0, 130)) +
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        plot.title = element_text(size = 30, face = "bold", hjust = 0.5, vjust = -8, margin = margin(b = 10) ),
        axis.text    = element_text(size   = 25, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 27),  
        axis.title.y = element_text(face = "bold", size = 27))+
  labs(title = "(B1) Plasma")
print(plot.A1)

plot.A2 <- ggplot() +
  geom_line(data = df.sim.A2, aes(Time, Plasma), col = "#5A5A5A", lwd = 1.5) +
  geom_point(data = OBS.A2, aes(Time, Plasma), col = "#e3716e", size = 5) +
  geom_errorbar(data = OBS.A2.1, aes(x = Time, ymin = Plasma - SD, ymax = Plasma + SD), col = "#e3716e", width = 5, size = 1.5) +
  ylab("GenX Concentration in Plasma (μg/ml)")+
  xlab("Time (hours)") +
  xlim(c(0, 130)) +
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        plot.title = element_text(size = 30, face = "bold", hjust = 0.5, vjust = -8, margin = margin(b = 10) ),
        axis.text    = element_text(size   = 25, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 27),  
        axis.title.y = element_text(face = "bold", size = 27))+
  labs(title = "(B2) Plasma")
print(plot.A2)

plot.B1 <- ggplot() +
  geom_line(data = df.sim.B, aes(Time, Plasma), col = "#5A5A5A", lwd = 1.5) +
  geom_point(data = OBS.B1, aes(Time, Plasma), col = "#e3716e", size = 5) +
  geom_errorbar(data = OBS.B1.1, aes(x = Time, ymin = Plasma - SD, ymax = Plasma + SD), col = "#e3716e", width = 5, size = 1.5) +
  ylab("GenX Concentration in Plasma (μg/ml)")+
  xlab("Time (hours)") +
  xlim(c(0, 420)) +
  ylim(c(0, 2.1)) +
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        plot.title = element_text(size = 30, face = "bold", hjust = 0.5, vjust = -8, margin = margin(b = 10) ),
        axis.text    = element_text(size   = 25, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 27),  
        axis.title.y = element_text(face = "bold", size = 27))+
  labs(title = "(B3) Plasma")
print(plot.B1)

plot.B2 <- ggplot() +
  geom_line(data = df.sim.B, aes(Time, Liver), col = "#5A5A5A", lwd = 1.5) +
  geom_point(data = OBS.B2, aes(Time, Liver), col = "#e3716e", size = 5) +
  geom_errorbar(data = OBS.B2.1, aes(x = Time, ymin = Liver - SD, ymax = Liver + SD), col = "#e3716e", width = 5, size = 1.5) +
  ylab("GenX Concentration in liver (ug/g)")+
  xlab("Time (hours)") +
  xlim(c(0, 420)) +
  ylim(c(0, 6.7)) +
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        plot.title = element_text(size = 30, face = "bold", hjust = 0.5, vjust = -8, margin = margin(b = 10) ),
        axis.text    = element_text(size   = 25, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 27),  
        axis.title.y = element_text(face = "bold", size = 27))+
  labs(title = "(B4) Liver")
print(plot.B2)

plot.B3 <- ggplot() +
  geom_line(data = df.sim.B, aes(Time, Lung), col = "#5A5A5A", lwd = 1.5) +
  geom_point(data = OBS.B3, aes(Time, Lung), col = "#e3716e", size = 5) +
  geom_errorbar(data = OBS.B3.1, aes(x = Time, ymin = Lung - SD, ymax = Lung + SD), col = "#e3716e", width = 5, size = 1.5) +
  ylab("GenX Concentration in lung (ug/g)")+
  xlab("Time (hours)") +
  xlim(c(0, 420)) +
  ylim(c(0, 0.68)) +
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        plot.title = element_text(size = 30, face = "bold", hjust = 0.5, vjust = -8, margin = margin(b = 10) ),
        axis.text    = element_text(size   = 25, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 27),  
        axis.title.y = element_text(face = "bold", size = 27))+
  labs(title = "(B5) Lung")
print(plot.B3)

plot.B4 <- ggplot() +
  geom_line(data = df.sim.B, aes(Time, Kidney), col = "#5A5A5A", lwd = 1.5) +
  geom_point(data = OBS.B4, aes(Time, Kidney), col = "#e3716e", size = 5) +
  geom_errorbar(data = OBS.B4.1, aes(x = Time, ymin = Kidney - SD, ymax = Kidney + SD), col = "#e3716e", width = 5, size = 1.5) +
  ylab("GenX Concentration in kidney (ug/g)")+
  xlab("Time (hours)") +
  xlim(c(0, 420)) +
  ylim(c(0, 1.4)) +
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        plot.title = element_text(size = 30, face = "bold", hjust = 0.5, vjust = -8, margin = margin(b = 10) ),
        axis.text    = element_text(size   = 25, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 27),  
        axis.title.y = element_text(face = "bold", size = 27))+
  labs(title = "(B6) Kidney")
print(plot.B4)

plot.A2 <- plot.A2 + theme(axis.title.y = element_blank())
plot.B1 <- plot.B1 + theme(axis.title.y = element_blank())
plot.B2 <- plot.B2 + theme(axis.title.y = element_blank())
plot.B4 <- plot.B4 + theme(axis.title.y = element_blank())
tl3 <- rasterGrob(readTIFF("D:/zs/PBPK/2025GenX/Artwork/3tl.tif"))
ggsave("Figure S2-Female.tiff", grid.arrange(plot.A1, plot.A2, plot.B1, plot.B2 , plot.B3, plot.B4, tl3, ncol = 4,nrow = 2),
       path = "D:/zs/PBPK/2025GenX/Artwork",
       width = 80, height = 40, units = "cm", dpi = 320)

## Save the fitting results to RDS files
saveRDS(Fit, file = "Fit_R.rds") 
Fit_R <- readRDS("Fit_R.rds")

#################################################################################################################################
###NSC And Circle barplot function                                                                                             #
## plot modifed from "R graph gallery: https://www.r-graph-gallery.com/297-circular-barplot-with-groups/ "                     #
#################################################################################################################################
mice.theta.G <- log(c(
  # Physiological parameters
  BW                             = 0.02,
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
  GFRC                           = 41.04,
  GEC                            = 0.54,
  
  # Chemical-specific parameters (final mean values)
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
## combined NSC calculation and plot ##################################################################################
pred.mice <- function(pars, DOSE, TDoses = 1) {  
  
  ## Get out of log domain
  pars <- exp(pars)                   ## Return a list of exp from log scale
  
  ## Exposure scenario 
  BW               = 0.02                            ## kg, Mouse body weight                                      
  tinterval        = 24                               ## hr, Time interval                                 
  GDOSE            = DOSE                             ## Input oral dose
  DOSEoral         = GDOSE * BW                       ## Amount of oral dose
  
  # To create a exposure scenario
  ex.oral <- ev(ID   = 1,             ## One individual
                amt  = DOSEoral,      ## Amount of dose 
                ii   = tinterval,     ## Time interval
                addl = TDoses - 1,    ## Additional dosing 
                cmt  = "AST",         ## The dosing compartment: AST Stomach  
                replicate = FALSE)    ## No replicate
  
  # set up the exposure time
  # Single dose, 30 days；Repeat dose, exposure +7 days elimination
  if (TDoses == 1) { tsamp = tgrid(0, 24 * 14, 0.5)  
  } else { tsamp = tgrid(0, tinterval * (TDoses - 1) + 24 * 10, 0.5) }
  
  out <- 
    mod %>%                                                 # model object
    param(pars) %>%                                         # to update the parameters in the model subject
    update(atol = 1E-19, maxsteps = 500000) %>%             # solver setting
    mrgsim_d(data = ex.oral, tgrid = tsamp) %>%
    filter(time != 0)
  
  outdfA = cbind.data.frame(Time      = out$time, 
                            Plasma    = out$Plasma,
                            AUC_CPlas = out$AUC_CA,
                            AUC_CL    = out$AUC_CL,
                            AUC_CK    = out$AUC_CK,
                            AUC_CLu   = out$AUC_CLu,
                            Dose      = DOSE,
                            Regimen   = ifelse(TDoses == 1, "Single", "Repeated")) 
  return(list("G" = outdfA))
}

par.G <- Fit_R$par
mice.theta.G[names(par.G)]  <- as.numeric(par.G) 

NSC_func <- function(pars, Pred, doses, regimens) {
  nG <- length(pars)
  n_doses <- length(doses)
  
  NSC_GCA <- array(NA, dim = c(nG, 4, n_doses),
                   dimnames = list(names(pars), 
                                   c("NSC_CPlas", "NSC_CL", "NSC_CK", "NSC_CLu"),
                                   paste0("Dose_", doses, "mgkg_", regimens)))
  
  for (d in 1:n_doses) {
    dose <- doses[d]
    regimen <- regimens[d]
    TDoses <- ifelse(regimen == "Single", 1, 14)  
    
    for (i in 1:nG) {
      pars.new <- pars %>% replace(i, log(exp((pars[i])) * 1.01))
      Mnew.G <- Pred(pars.new, dose, TDoses)
      M.G <- Pred(pars, dose, TDoses)
      delta.Gpars <- exp(pars[i]) / (exp(pars[i]) * 0.01)
      
      AUC.GPlas.new <- tail(Mnew.G$G$AUC_CPlas, 1)
      AUC.GPlas.ori <- tail(M.G$G$AUC_CPlas, 1)
      AUC.GL.new <- tail(Mnew.G$G$AUC_CL, 1)
      AUC.GL.ori <- tail(M.G$G$AUC_CL, 1)
      AUC.GK.new <- tail(Mnew.G$G$AUC_CK, 1)
      AUC.GK.ori <- tail(M.G$G$AUC_CK, 1)
      AUC.GLu.new <- tail(Mnew.G$G$AUC_CLu, 1)
      AUC.GLu.ori <- tail(M.G$G$AUC_CLu, 1)
      
      delta.AUC.GPlas <- AUC.GPlas.new - AUC.GPlas.ori
      delta.AUC.GL <- AUC.GL.new - AUC.GL.ori
      delta.AUC.GK <- AUC.GK.new - AUC.GK.ori
      delta.AUC.GLu <- AUC.GLu.new - AUC.GLu.ori
      
      NSC_GCA[i, 1, d] <- as.numeric((delta.AUC.GPlas / AUC.GPlas.ori) * delta.Gpars)
      NSC_GCA[i, 2, d] <- as.numeric((delta.AUC.GL / AUC.GL.ori) * delta.Gpars)
      NSC_GCA[i, 3, d] <- as.numeric((delta.AUC.GK / AUC.GK.ori) * delta.Gpars)
      NSC_GCA[i, 4, d] <- as.numeric((delta.AUC.GLu / AUC.GLu.ori) * delta.Gpars)
    }
  }
  
  return(NSC_GCA)
}

doses <- c(10, 30, 1)  
regimens <- c("Single", "Single","Repeated")  

NSC_results <- NSC_func(mice.theta.G, pred.mice, doses, regimens)

process_and_save_NSC_results <- function(NSC_results, doses, regimens) {
  all_NSC_dfs <- list()
  
  for (i in 1:length(doses)) {
    dose <- doses[i]
    regimen <- regimens[i]
    NSC_df <- as.data.frame(NSC_results[, , i])
    rownames(NSC_df) <- names(mice.theta.G)
    NSC_df$Dose_Group <- paste0("Dose_", dose, "mgkg_", regimen)
    NSC_df$Parameter <- rownames(NSC_df)
    
    all_NSC_dfs[[i]] <- NSC_df
  }
  
  combined_NSC <- do.call(rbind, all_NSC_dfs)
  combined_NSC <- combined_NSC[, c("Parameter", "Dose_Group", "NSC_CPlas", "NSC_CL", "NSC_CK", "NSC_CLu")]
  numeric_cols <- c("NSC_CPlas", "NSC_CL", "NSC_CK", "NSC_CLu")
  for (col in numeric_cols) {
    combined_NSC[[col]] <- ifelse(abs(combined_NSC[[col]]) < 1e-5, "<1e-5", 
                                  format(combined_NSC[[col]], scientific = TRUE, digits = 4))
  }
  write.csv(combined_NSC, file = 'NSC_all_groups.csv', row.names = FALSE)
  
  return(combined_NSC)
}
NSC_combined <- process_and_save_NSC_results(NSC_results, doses, regimens)

NSC_significant <- NSC_combined %>%
  mutate_at(vars(NSC_CPlas, NSC_CL, NSC_CK, NSC_CLu), 
            ~as.numeric(ifelse(. == "<1e-5", 0, .))) %>%
  filter(if_any(c(NSC_CPlas, NSC_CL, NSC_CK, NSC_CLu), ~abs(.) >= 0.3))

write.csv(NSC_significant, file = 'NSC_significant.csv', row.names = FALSE)

##plot
NSC_GCA_M <- read.csv(file = "NSC.csv", row.names = 1)  
NSC_GCA_M[] <- lapply(NSC_GCA_M, function(x) {
  if(is.character(x)) { x <- gsub("<1e-5", "0", x)
  as.numeric(x)
  } else { x }})

print(head(NSC_GCA_M))
print(str(NSC_GCA_M))

melt.Plas <- data.frame(value = NSC_GCA_M$NSC_CPlas)
melt.Plas$group <- "Plasma"
melt.Plas$par <- rownames(NSC_GCA_M)

melt.Liver <- data.frame(value = NSC_GCA_M$NSC_CL)
melt.Liver$group <- "Liver"
melt.Liver$par <- rownames(NSC_GCA_M)

melt.Kidney <- data.frame(value = NSC_GCA_M$NSC_CK)
melt.Kidney$group <- "Kidney"
melt.Kidney$par <- rownames(NSC_GCA_M)

melt.Lung <- data.frame(value = NSC_GCA_M$NSC_CLu)
melt.Lung$group <- "Lung"
melt.Lung$par <- rownames(NSC_GCA_M)

melt.data <- rbind(melt.Plas, melt.Liver, melt.Kidney, melt.Lung)
melt.data$value <- as.numeric(melt.data$value)

melt.data <- melt.data[!is.na(melt.data$value), ]

data <- melt.data %>% filter(abs(value) >= 0.3)
if (nrow(data) == 0) {
  stop("not find NSC >=0.3")
}

data$group <- factor(data$group)

empty_bar <- 4
to_add <- data.frame(matrix(NA, empty_bar * nlevels(data$group), ncol(data)))
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each = empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))

label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id - 0.5) / number_of_bar
label_data$hjust <- ifelse(angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle + 180, angle)

base_data <- data %>% 
  group_by(group) %>% 
  summarize(start = min(id) - 0.2, end = max(id) - empty_bar + 0.2) %>% 
  rowwise() %>% 
  mutate(title = mean(c(start, end)))

grid_data <- base_data
grid_data$end <- grid_data$end[c(nrow(grid_data), 1:nrow(grid_data) - 1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1, ]

windowsFonts(Times = windowsFont("Times New Roman"))
p1 <- ggplot(data, aes(x = as.factor(id), y = value, fill = group)) +
  geom_bar(aes(x = as.factor(id), y = value, fill = group), stat = "identity", alpha = 0.5, na.rm = TRUE) +
  geom_segment(data = grid_data, aes(x = end, y = 90, xend = start, yend = 90), 
               colour = "grey", alpha = 1, size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = grid_data, aes(x = end, y = 60, xend = start, yend = 60), 
               colour = "grey", alpha = 1, size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = grid_data, aes(x = end, y = 30, xend = start, yend = 30), 
               colour = "grey", alpha = 1, size = 0.3, inherit.aes = FALSE) +
  
  annotate("text", x = rep(max(data$id, na.rm = TRUE), 3), y = c(30, 60, 90), 
           label = c("30%", "60%", "90%"), color = "red", size = 4, 
           angle = 0, fontface = "bold", hjust = 1) +
  
  geom_bar(aes(x = as.factor(id), y = abs(value * 100), fill = group), stat = "identity", alpha = 0.5, na.rm = TRUE) +
  ylim(-100, 230) +
  labs(tag = "(B)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(family = "Times"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1, 4), "cm"), 
    plot.tag = element_text(size = 25, face = "bold", margin = margin(t = 150, l = 150)),  
    plot.tag.position = c(0, 1) 
  ) +
  coord_polar() + 
  geom_text(data = label_data, aes(x = id, y = abs(value * 100) + 10, label = par, hjust = hjust), 
            color = "black", fontface = "bold", alpha = 0.6, size = 5, 
            angle = label_data$angle, inherit.aes = FALSE, na.rm = TRUE) +
  
  geom_segment(data = base_data, aes(x = start, y = -5, xend = end, yend = -5), 
               colour = "black", alpha = 0.8, size = 0.6, inherit.aes = FALSE) +
  geom_text(data = base_data, aes(x = title, y = -38, label = group), 
            hjust = 0.4, colour = "black", alpha = 0.8, size = 6, 
            fontface = "bold", inherit.aes = FALSE)
print(p1)
p1 <- p1 + theme( plot.background = element_rect(fill = "white", color = NA))
ggsave("Figure S6-mice-Female.tiff",scale = 1,
       plot = p1,
       path = "D:/zs/PBPK/2025GenX/Artwork",
       width = 24, height = 24, units = "cm",dpi=320)

#############################################################################################################################
##MCMC                                                                                                                      #
## normal AND lognormal Calculation #########################################################################################
###normal
mean_value <- 41.04
sd_value   <-	12.31

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
  BW            = 0.02                            ## kg, Mice body weight                                      
  tinterval     = 24                               ## hr, Time interval                                 
  TDoses        = 1                                ## Dose times                                    
  DOSE          = 10                               ## Single oral dose from Gannon et al. (2015)
  DOSEoral      = DOSE*BW                          ## Amount of oral dose
  
  ## Amount of oral dose
  idata <- 
    tibble(ID = 1:N) %>% 
    mutate( BW        = rnormTrunc  (N, min = 0.0082, max = 0.0318, mean = 0.02 ,      sd = 0.006 ),   #CV=0.3
            VKC       = rnormTrunc  (N, min = 0.007 , max = 0.0270, mean = 0.017,      sd = 0.0051),   #CV=0.3
            GFRC      = rnormTrunc  (N, min = 16.913, max = 65.167, mean = 41.04,      sd = 12.31 ) ,  #CV=0.3
            Free                           = 0.0704,                ### fitting parameters
            PL                             = 1.0649,                 ### fitting parameters
            PLu                            = 0.331,                 ### fitting parameters
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

PlotDat_A1     <- Mice_A1 (R_Gpars,N = N)[[1]]  %>% select (ID = ID, Time = Time, Conc = CPlas)  %>% mutate( Tissue = "Plasma")
PlotDat_A1_m   <- Mice_A1 (R_Gpars,N = N)[[2]]  %>% select (Time = Time, Conc = CPlas)  %>% mutate( Tissue = "Plasma")
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
  geom_ribbon  (data = PlotDat_A1_summary, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), fill="#CECCE5", alpha=0.3) +
  geom_ribbon  (data = PlotDat_A1_summary, aes(x = Time, ymin = ci_q1, ymax = ci_q3), fill="#715ea9", alpha = 0.3) +
  geom_line    (data = PlotDat_A1_m, aes(x = Time, y = Conc), colour = "#5A5A5A", lwd = 0.8) +
  geom_point   (data = OBS.A1, aes(x = Time/24, y = Plasma), colour = "#715ea9", size = 2.5) +
  geom_errorbar(data = OBS.A1, aes(x = Time/24, ymin= Plasma-SD, ymax = Plasma+SD), col = "#715ea9", width = 0.05, size = 0.8)+
  ylab("GenX Concentration in Plasma (μg/ml)")+
  xlab("Time (days)")+
  scale_x_continuous(trans = "log1p") +  # log Time
  scale_y_continuous(limits = c(-1, 65), expand = c(0, 0))+
  theme_bw() +
  theme(panel.border = element_rect(color = "grey10",fill = NA, linewidth = 1),
        strip.text   = element_text(size = rel(2), colour = "grey10"),
        axis.text    = element_text(size = rel(2), colour = "black"),
        axis.line    = element_line(color = "black", size = 0.3),
        axis.title   = element_text(size = 25, colour = "black", face = "bold"),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.5, vjust = -8, margin = margin(b = 10) ),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+
  labs(title = "(B1) Plasma")
plot.A1_Plas

## A2 MCMC ########################################################
Mice_A2 <- function (pars, N) {
  
  pars <- exp(pars)           ## Return a list of exp from log scale
  
  ## Exposure scenario 
  BW            = 0.02                            ## kg, Rat body weight                                      
  tinterval     = 24                               ## hr, Time interval                                 
  TDoses        = 1                                ## Dose times                                    
  DOSE          = 30                               ## Single oral dose from Gannon et al. (2015)
  DOSEoral      = DOSE*BW                          ## Amount of oral dose
  
  ## Amount of oral dose
  idata <- 
    tibble(ID = 1:N) %>% 
    mutate( BW        = rnormTrunc  (N, min = 0.0082, max = 0.0318, mean = 0.02 ,      sd = 0.006 ),   #CV=0.3
            VKC       = rnormTrunc  (N, min = 0.007 , max = 0.0270, mean = 0.017,      sd = 0.0051),   #CV=0.3
            GFRC      = rnormTrunc  (N, min = 16.913, max = 65.167, mean = 41.04,      sd = 12.31 ) ,  #CV=0.3
            Free                           = 0.0704,                ### fitting parameters
            PL                             = 1.0649,                 ### fitting parameters
            PLu                            = 0.331,                 ### fitting parameters
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

PlotDat_A2     <- Mice_A2 (R_Gpars,N = N)[[1]]  %>% select (ID = ID, Time = Time, Conc = CPlas)  %>% mutate( Tissue = "Plasma")
PlotDat_A2_m   <- Mice_A2 (R_Gpars,N = N)[[2]]  %>% select (Time = Time, Conc = CPlas)  %>% mutate( Tissue = "Plasma")
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
  geom_ribbon  (data = PlotDat_A2_summary, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), fill="#CECCE5", alpha=0.3) +
  geom_ribbon  (data = PlotDat_A2_summary, aes(x = Time, ymin = ci_q1, ymax = ci_q3), fill="#715ea9", alpha = 0.3) +
  geom_line    (data = PlotDat_A2_m, aes(x = Time, y = Conc), colour = "#5A5A5A", size = 0.8) +
  geom_point   (data = OBS.A2, aes(x = Time/24, y = Plasma), colour = "#715ea9", size = 2.5) +
  geom_errorbar(data = OBS.A2, aes(x = Time/24, ymin= Plasma-SD, ymax = Plasma+SD), col = "#715ea9", width = 0.05, size = 0.8)+
  ylab("GenX Concentration in Plasma (μg/ml)")+
  xlab("Time (days)")+
  scale_x_continuous(trans = "log1p") +  # log Time
  scale_y_continuous(limits = c(-3, 160), expand = c(0, 0))+
  theme_bw() +
  theme(panel.border = element_rect(color = "grey10",fill = NA, linewidth = 1),
        strip.text   = element_text(size = rel(2), colour = "grey10"),
        axis.text    = element_text(size = rel(2), colour = "black"),
        axis.line    = element_line(color = "black", size = 0.3),
        axis.title   = element_text(size = 25, colour = "black", face = "bold"),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.5, vjust = -8, margin = margin(b = 10) ),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+
  labs(title = "(B2) Plasma")
plot.A2_Plas

## B MCMC #######################################################
Mice_B1 <- function (pars, N) {
  
  pars <- exp(pars)           ## Return a list of exp from log scale
  
  ## Exposure scenario 
  BW            = 0.02                            ## kg, Rat body weight                                      
  tinterval     = 24                               ## hr, Time interval                                 
  TDoses        = 14                                ## Dose times                                    
  DOSE          = 1                              ## Repeat oral dose from Wen et al. (2022)
  DOSEoral      = DOSE*BW                          ## Amount of oral dose
  
  ## Amount of oral dose
  idata <- 
    tibble(ID = 1:N) %>% 
    mutate( BW        = rnormTrunc  (N, min = 0.0082, max = 0.0318, mean = 0.02 ,      sd = 0.006 ),   #CV=0.3
            VKC       = rnormTrunc  (N, min = 0.007 , max = 0.0270, mean = 0.017,      sd = 0.0051),   #CV=0.3
            GFRC      = rnormTrunc  (N, min = 16.913, max = 65.167, mean = 41.04,      sd = 12.31 ) ,  #CV=0.3
            Free                           = 0.0704,                ### fitting parameters
            PL                             = 1.0649,                 ### fitting parameters
            PLu                            = 0.331,                 ### fitting parameters
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
  tsamp = tgrid(0,tinterval*(TDoses-1)+24*7,0.05)
  
  # Combine data and run the simulation
  out_1 <- mod %>%data_set(ex_1) %>%
    idata_set(idata) %>% 
    update(atol = 1e-6, maxstep = 500000) %>%
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
    update(atol = 1E-6, maxsteps = 500000) %>%          
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

#B1
PlotDat_B1     <- Mice_B1 (R_Gpars,N = N)[[1]]  %>% select (ID = ID, Time = Time, Conc = CPlas)  %>% mutate( Tissue = "Plasma")
PlotDat_B1_m   <- Mice_B1 (R_Gpars,N = N)[[2]]  %>% select (Time = Time, Conc = CPlas)  %>% mutate( Tissue = "Plasma")
OBS.B1  <- Data_A1 %>% filter(Study == 2 & Sample == "blood"  & Dose == 1 )   %>% select(Time = "Time", Plasma = "Conc", SD = "SD") # TDosesB = 14

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
  geom_ribbon  (data = PlotDat_B1_summary, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), fill="#CECCE5", alpha=0.3) +
  geom_ribbon  (data = PlotDat_B1_summary, aes(x = Time, ymin = ci_q1, ymax = ci_q3), fill="#715ea9", alpha = 0.3) +
  geom_line    (data = PlotDat_B1_m, aes(x = Time, y = Conc), colour = "#5A5A5A", lwd = 0.8) +
  geom_point   (data = OBS.B1, aes(x = Time/24, y = Plasma), colour = "#715ea9", size = 2.5) +
  geom_errorbar(data = OBS.B1, aes(x = Time/24, ymin= Plasma-SD, ymax = Plasma+SD), col = "#715ea9", width = 0.05, size = 0.8)+
  ylab("GenX Concentration in Plasma (μg/ml)")+
  xlab("Time (days)")+
  scale_x_continuous(trans = "log1p") +  # log Time
  scale_y_continuous(limits = c(-0.5, 4), expand = c(0, 0))+
  theme_bw() +
  theme(panel.border = element_rect(color = "grey10",fill = NA, linewidth = 1),
        strip.text   = element_text(size = rel(2), colour = "grey10"),
        axis.text    = element_text(size = rel(2), colour = "black"),
        axis.line    = element_line(color = "black", size = 0.3),
        axis.title   = element_text(size = 25, colour = "black", face = "bold"),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.5, vjust = -8, margin = margin(b = 10) ),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+
  labs(title = "(B3) Plasma")
plot.B1_Plas

#B2 
PlotDat_B2     <- Mice_B1 (R_Gpars,N = N)[[1]]  %>% select (ID = ID, Time = Time, Conc = CL)  %>% mutate( Tissue = "Liver")
PlotDat_B2_m   <- Mice_B1 (R_Gpars,N = N)[[2]]  %>% select (Time = Time, Conc = CL)  %>% mutate( Tissue = "Liver")
OBS.B2  <- Data_A1 %>% filter(Study == 2 & Sample == "Liver"   & Dose == 1 )   %>% select(Time = "Time", Liver  = "Conc", SD = "SD") # TDosesB = 14

PlotDat_B2_summary <- PlotDat_B2 %>%
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

plot.B2_Liver <- 
  ggplot() + 
  geom_ribbon  (data = PlotDat_B2_summary, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), fill="#CECCE5", alpha=0.3) +
  geom_ribbon  (data = PlotDat_B2_summary, aes(x = Time, ymin = ci_q1, ymax = ci_q3), fill="#715ea9", alpha = 0.3) +
  geom_line    (data = PlotDat_B2_m, aes(x = Time, y = Conc), colour = "#5A5A5A", lwd = 0.8) +
  geom_point   (data = OBS.B2, aes(x = Time/24, y = Liver), colour = "#715ea9", size = 2.5) +
  geom_errorbar(data = OBS.B2, aes(x = Time/24, ymin= Liver-SD, ymax = Liver+SD), col = "#715ea9", width = 0.05, size = 0.8)+
  ylab("GenX Concentration in liver (ug/g)")+
  xlab("Time (days)")+
  scale_x_continuous(trans = "log1p") +  # log Time
  scale_y_continuous(limits = c(-0.5, 8.5), expand = c(0, 0))+
  theme_bw() +
  theme(panel.border = element_rect(color = "grey10",fill = NA, linewidth = 1),
        strip.text   = element_text(size = rel(2), colour = "grey10"),
        axis.text    = element_text(size = rel(2), colour = "black"),
        axis.line    = element_line(color = "black", size = 0.3),
        axis.title   = element_text(size = 25, colour = "black", face = "bold"),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.5, vjust = -8, margin = margin(b = 10) ),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+
  labs(title = "(B4) Liver")
plot.B2_Liver

#B3
PlotDat_B3     <- Mice_B1 (R_Gpars,N = N)[[1]]  %>% select (ID = ID, Time = Time, Conc = CLu)  %>% mutate( Tissue = "Lung")
PlotDat_B3_m   <- Mice_B1 (R_Gpars,N = N)[[2]]  %>% select (Time = Time, Conc = CLu)  %>% mutate( Tissue = "Lung")
OBS.B3  <- Data_A1 %>% filter(Study == 2 & Sample == "Lung"   & Dose == 1 )   %>% select(Time = "Time", Lung   = "Conc", SD = "SD") # TDosesB = 14

PlotDat_B3_summary <- PlotDat_B3 %>%
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

plot.B3_Lung <- 
  ggplot() + 
  geom_ribbon  (data = PlotDat_B3_summary, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), fill="#CECCE5", alpha=0.3) +
  geom_ribbon  (data = PlotDat_B3_summary, aes(x = Time, ymin = ci_q1, ymax = ci_q3), fill="#715ea9", alpha = 0.3) +
  geom_line    (data = PlotDat_B3_m, aes(x = Time, y = Conc), colour = "#5A5A5A", lwd = 0.8) +
  geom_point   (data = OBS.B3, aes(x = Time/24, y = Lung), colour = "#715ea9", size = 2.5) +
  geom_errorbar(data = OBS.B3, aes(x = Time/24, ymin= Lung-SD, ymax = Lung+SD), col = "#715ea9", width = 0.05, size = 0.8)+
  ylab("GenX Concentration in Lung (ug/g)")+
  xlab("Time (days)")+
  scale_x_continuous(trans = "log1p") +  # log Time
  scale_y_continuous(limits = c(-0.2, 1.2), expand = c(0, 0))+
  theme_bw() +
  theme(panel.border = element_rect(color = "grey10",fill = NA, linewidth = 1),
        strip.text   = element_text(size = rel(2), colour = "grey10"),
        axis.text    = element_text(size = rel(2), colour = "black"),
        axis.line    = element_line(color = "black", size = 0.3),
        axis.title   = element_text(size = 25, colour = "black", face = "bold"),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.5, vjust = -8, margin = margin(b = 10) ),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+
  labs(title = "(B5) Lung")
plot.B3_Lung

#B4
PlotDat_B4     <- Mice_B1 (R_Gpars,N = N)[[1]]  %>% select (ID = ID, Time = Time, Conc = CK)  %>% mutate( Tissue = "Kidney")
PlotDat_B4_m   <- Mice_B1 (R_Gpars,N = N)[[2]]  %>% select (Time = Time, Conc = CK)  %>% mutate( Tissue = "Kidney")
OBS.B4  <- Data_A1 %>% filter(Study == 2 & Sample == "Kidney" & Dose == 1 )   %>% select(Time = "Time", Kidney = "Conc", SD = "SD") # TDosesB = 14

PlotDat_B4_summary <- PlotDat_B4 %>%
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

plot.B4_Kidney <- 
  ggplot() + 
  geom_ribbon  (data = PlotDat_B4_summary, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est), fill="#CECCE5", alpha=0.3) +
  geom_ribbon  (data = PlotDat_B4_summary, aes(x = Time, ymin = ci_q1, ymax = ci_q3), fill="#715ea9", alpha = 0.3) +
  geom_line    (data = PlotDat_B4_m, aes(x = Time, y = Conc), colour = "#5A5A5A", lwd = 0.8) +
  geom_point   (data = OBS.B4, aes(x = Time/24, y = Kidney), colour = "#715ea9", size = 2.5) +
  geom_errorbar(data = OBS.B4, aes(x = Time/24, ymin= Kidney-SD, ymax = Kidney+SD), col = "#715ea9", width = 0.05, size = 0.8)+
  ylab("GenX Concentration in Kidney (ug/g)")+
  xlab("Time (days)")+
  scale_x_continuous(trans = "log1p") +  # log Time
  scale_y_continuous(limits = c(-0.5, 3), expand = c(0, 0))+
  theme_bw() +
  theme(panel.border = element_rect(color = "grey10",fill = NA, linewidth = 1),
        strip.text   = element_text(size = rel(2), colour = "grey10"),
        axis.text    = element_text(size = rel(2), colour = "black"),
        axis.line    = element_line(color = "black", size = 0.3),
        axis.title   = element_text(size = 25, colour = "black", face = "bold"),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.5, vjust = -8, margin = margin(b = 10) ),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+
  labs(title = "(B6) Kidney")
plot.B4_Kidney

grid.arrange(plot.A1_Plas, plot.A2_Plas, plot.B1_Plas, plot.B2_Liver, plot.B3_Lung, plot.B4_Kidney, ncol = 4, nrow = 2)

plot.A2_Plas   <- plot.A2_Plas + theme(axis.title.y = element_blank())
plot.B1_Plas   <- plot.B1_Plas + theme(axis.title.y = element_blank())
plot.B2_Liver  <- plot.B2_Liver + theme(axis.title.y = element_blank())
plot.B4_Kidney <- plot.B4_Kidney + theme(axis.title.y = element_blank())
ggsave("Figure MCMC-Female.tiff", grid.arrange(plot.A1_Plas, plot.A2_Plas, plot.B1_Plas, plot.B2_Liver, plot.B3_Lung, plot.B4_Kidney, ncol = 4, nrow = 2),
       path = "D:/zs/PBPK/2025GenX/Artwork",
       width = 80, height = 40, units = "cm", dpi = 320)

###############################################################################################################################
##Evaluation                                                                                                                  #
#Mice                                                                                                                         #
#oral daily dose to 0.32, 1, 10, 100 mg/kg-d for 1, 5, 12, 14, 28, 29 and 84 days #
###############################################################################################################################
## Read the data and later used in model calibration and evaluation
evaluation <- read.csv(file = "Evalution-F.csv")

OBS.A_CL  <- evaluation %>% filter(Study == 1 & Sample == "Liver" ) %>% select(Time = "Time", CL     = "Conc", SD = "SD", Study = "Study", Dose = "Dose") # TDoses = 28
OBS.B_Pla <- evaluation %>% filter(Study == 2 & Sample == "Plasma") %>% select(Time = "Time", Plasma = "Conc", SD = "SD", Study = "Study", Dose = "Dose") # TDoses = 84
OBS.C_Pla <- evaluation %>% filter(Study == 3 & Sample == "serum" ) %>% select(Time = "Time", Plasma = "Conc", SD = "SD", Study = "Study", Dose = "Dose") # TDoses = 28

pars.eva <- log(c(
  # Physiological parameters
  BW                             = 0.02,
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
  GFRC                           = 41.04,
  GEC                            = 0.54,
  
  # Chemical-specific parameters (final mean values)
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

pred.eva <- function(pars) { ## pars: input parameters, Dose: input dose, Dose regimen: 0.0004, 0.002, 0.0032, 0.01, 0.32, 0.4, 1, 2, 10, 100 mg/kg/day
  
  ## Get out of log domain
  pars <- exp(pars)                   ## Return a list of exp from log scale
  
  ## Exposure scenario 
  BW          = 0.02                  ## Body weight 
  tinterval   = 24                     ## Time interval; 
  
  ##1
  TDOSE1    = 28                     
  DOSE1     = 100                    
  DOSEoral1 = DOSE1*BW               
  ex.oral1  <- ev(ID = 1, amt = DOSEoral1, ii = tinterval, addl = TDOSE1 - 1, cmt = "AST", replicate = FALSE)
  tsamp1    = tgrid(0, tinterval*(TDOSE1 - 1) + 24*7, 1) 
  
  out1 <- 
    mod %>%          
    param (pars) %>%
    Req   (Liver)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral1, tgrid = tsamp1)   
  
  outdf1 <- cbind.data.frame(Time = out1$time/24, CL = out1$Liver) 
  
  ##2
  TDOSE2      = 84                      ## Total dosing/Dose times
  DOSE2       = 0.32                   ## Input oral dose
  DOSEoral2   = DOSE2*BW               ## Amount of oral dose
  ex.oral2    <- ev(ID = 1, amt = DOSEoral2, ii = tinterval, addl = TDOSE2 - 1, cmt = "AST", replicate = FALSE)
  tsamp2      = tgrid(0, tinterval*(TDOSE2 - 1) + 24*7, 1) 
  
  out2 <- 
    mod %>%          
    param (pars) %>%
    Req   (Plasma)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral2, tgrid = tsamp2)   
  
  outdf2 <- cbind.data.frame(Time = out2$time/24, Plasma = out2$Plasma) 
  
  ##3
  TDOSE3      = 28
  tsamp3      = tgrid(0, tinterval*(TDOSE3 - 1) + 24*7, 1)
  
  #a
  DOSE3a      = 1                 
  DOSEoral3a  = DOSE3a*BW              
  ex.oral3a   <- ev(ID = 1, amt = DOSEoral3a, ii = tinterval, addl = TDOSE3 - 1, cmt = "AST", replicate = FALSE)
  out3a <- 
    mod %>%          
    param (pars) %>%
    Req   (Plasma)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral3a, tgrid = tsamp3)   
  
  outdf3a <- cbind.data.frame(Time   = out3a$time/24, 
                              Plasma = out3a$Plasma) 
  #b
  DOSE3b      = 10                 
  DOSEoral3b  = DOSE3b*BW              
  ex.oral3b   <- ev(ID = 1, amt = DOSEoral3b, ii = tinterval, addl = TDOSE3 - 1, cmt = "AST", replicate = FALSE)
  out3b <- 
    mod %>%          
    param (pars) %>%
    Req   (Plasma)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral3b, tgrid = tsamp3)   
  
  outdf3b <- cbind.data.frame(Time   = out3b$time/24, 
                              Plasma = out3b$Plasma) 
  #c
  DOSE3c      = 100                 
  DOSEoral3c  = DOSE3c*BW              
  ex.oral3c   <- ev(ID = 1, amt = DOSEoral3c, ii = tinterval, addl = TDOSE3 - 1, cmt = "AST", replicate = FALSE)
  out3c <- 
    mod %>%          
    param (pars) %>%
    Req   (Plasma)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral3c, tgrid = tsamp3)   
  
  outdf3c <- cbind.data.frame(Time   = out3c$time/24, 
                              Plasma = out3c$Plasma) 
  
  return (list("outdf1"   = outdf1, 
               "outdf2"   = outdf2,
               "outdf3a"  = outdf3a, "outdf3b"  = outdf3b, "outdf3c"  = outdf3c)) 
}

out_1  <- pred.eva(pars.eva)[[1]]  %>% mutate(SD = 0, Study = 1, Dose = 100 )%>% filter(Time > 0)
out_2  <- pred.eva(pars.eva)[[2]]  %>% mutate(SD = 0, Study = 2, Dose = 0.32)%>% filter(Time > 0)
out_3a <- pred.eva(pars.eva)[[3]]  %>% mutate(SD = 0, Study = 3, Dose = 1   )%>% filter(Time > 0)
out_3b <- pred.eva(pars.eva)[[4]]  %>% mutate(SD = 0, Study = 3, Dose = 10  )%>% filter(Time > 0)
out_3c <- pred.eva(pars.eva)[[5]]  %>% mutate(SD = 0, Study = 3, Dose = 100 )%>% filter(Time > 0)
out_3  <- rbind.data.frame (out_3a, out_3b, out_3c)

out_1_CL  <- out_1 %>% filter(Study == 1, Time == 29)
out_2_Pla <- out_2 %>% filter(Study == 2, Time == 84)
out_3_Pla <- out_3 %>% filter(Study == 3, Time %in% c(1, 5, 14, 28))

out_Pla <- rbind.data.frame (out_2_Pla, out_3_Pla)
out_Pla <- out_Pla %>% mutate(Matrix = c("Pre.Plasma"))
out_CL  <- rbind.data.frame (out_1_CL)
out_CL  <- out_CL  %>% mutate(Matrix = c("Pre.Liver"))

OBS_Pla <- rbind.data.frame (OBS.B_Pla, OBS.C_Pla)
OBS_Pla <- OBS_Pla %>% mutate(Matrix = c("Obs.Plasma"))
OBS_CL  <- rbind.data.frame (OBS.A_CL)
OBS_CL  <- OBS_CL  %>% mutate(Matrix = c("Obs.Liver"))

############################ Pre and Obs comparsion Plot Liver #####################################################################################################################################
combined_CL  <- rbind(out_CL , OBS_CL )
p_CL1<- 
  ggplot(combined_CL, aes(x = as.factor(Dose), y = CL, fill = as.factor(Matrix))) + 
  geom_bar(stat = "identity", color = "black", size = 1.2, position = position_dodge(width = 0.4), width = 0.4) +
  geom_errorbar(aes(ymin = CL, ymax = CL+SD), width = 0.2, position = position_dodge(0.4),  size = 0.9) + 
  scale_fill_manual(values=c("#F6B3AC", "gray"))+##"#BFD0E1", "#7FB2D5"
  labs(title = "(B) Liver, Time=29d", x = "Dose (mg/kg/day)", y = "GenX Concentration in liver (ug/g)", fill = "Data Type") 

p_CL1 = p_CL1 + 
  scale_y_continuous(limits = c(-0.09, 4.5),breaks=c(1,2,3,4),expand = c(0,0))+
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
tl <- rasterGrob(readTIFF("D:/zs/PBPK/2025GenX/Artwork/4tl.tif"))
ggsave("Figure 3-Female-liver.tiff",scale = 1,
       plot = grid.arrange(p_CL1, tl , ncol = 2),
       path = "D:/zs/PBPK/2025GenX/Artwork",
       width = 40, height = 20, units = "cm",dpi=320)
############################ Pre and Obs comparsion Plot Plasma #####################################################################################################################
combined_Pla <- rbind(out_Pla, OBS_Pla)

combined_Pla1<- combined_Pla %>% filter(Study==3)%>% filter(Dose==1)
p_Pla1<- 
  ggplot(combined_Pla1, aes(x = as.factor(Time), y = Plasma, fill = as.factor(Matrix))) + # do the histogram plot by the factor of matrix, so you have a bar for pre.plasma, a bar for obs.plasma, so on.
  geom_bar(stat = "identity", color = "black", size = 1.2, position = position_dodge(width = 0.8), width = 0.8) +
  geom_errorbar(aes(ymin = Plasma, ymax = Plasma+SD), width = 0.2, position = position_dodge(0.8),  size = 0.9) + 
  scale_fill_manual(values=c("#F6B3AC","gray"))+##"#F6B3AC","#F47F72"
  labs(title = "(B1) Plasma, Dose=1 mg/kg/day", x = "Time (day)", y = "GenX Concentration in Plasma (μg/ml)", fill = "Data Type") 

p_Pla1 = p_Pla1 + 
  scale_y_continuous(limits = c(-0.05, 2.5),breaks=c(0.5,1,1.5,2),expand = c(0,0))+
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

combined_Pla2<- combined_Pla %>% filter(Study==3)%>% filter(Dose==10)
p_Pla2<- 
  ggplot(combined_Pla2, aes(x = as.factor(Time), y = Plasma, fill = as.factor(Matrix))) + # do the histogram plot by the factor of matrix, so you have a bar for pre.plasma, a bar for obs.plasma, so on.
  geom_bar(stat = "identity", color = "black", size = 1.2, position = position_dodge(width = 0.8), width = 0.8) +
  geom_errorbar(aes(ymin = Plasma, ymax = Plasma+SD), width = 0.2, position = position_dodge(0.8),  size = 0.9) + 
  scale_fill_manual(values=c("#F6B3AC","gray"))+
  labs(title = "(B2) Plasma, Dose=10 mg/kg/day", x = "Time (day)", y = "GenX Concentration in Plasma (μg/ml)", fill = "Data Type") 

p_Pla2 = p_Pla2 + 
  scale_y_continuous(limits = c(-0.13, 6.5),breaks=c(1.5,3,4.5,6),expand = c(0,0))+
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

combined_Pla3<- combined_Pla %>% filter(Study==3)%>% filter(Dose==100)
p_Pla3<- 
  ggplot(combined_Pla3, aes(x = as.factor(Time), y = Plasma, fill = as.factor(Matrix))) + # do the histogram plot by the factor of matrix, so you have a bar for pre.plasma, a bar for obs.plasma, so on.
  geom_bar(stat = "identity", color = "black", size = 1.2, position = position_dodge(width = 0.8), width = 0.8) +
  geom_errorbar(aes(ymin = Plasma, ymax = Plasma+SD), width = 0.2, position = position_dodge(0.8),  size = 0.9) + 
  scale_fill_manual(values=c("#F6B3AC","gray"))+
  labs(title = "(B3) Plasma, Dose=100 mg/kg/day", x = "Time (day)", y = "GenX Concentration in Plasma (μg/ml)", fill = "Data Type") 

p_Pla3 = p_Pla3 + 
  scale_y_continuous(limits = c(-0.7, 35),breaks=c(10,20,30),expand = c(0,0))+
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

combined_Pla4<- combined_Pla %>% filter(Time==84)
p_Pla4<- 
  ggplot(combined_Pla4, aes(x = as.factor(Time), y = Plasma, fill = as.factor(Matrix))) + # do the histogram plot by the factor of matrix, so you have a bar for pre.plasma, a bar for obs.plasma, so on.
  geom_bar(stat = "identity", color = "black", size = 1.2, position = position_dodge(width = 0.5), width = 0.5) +
  geom_errorbar(aes(ymin = Plasma, ymax = Plasma+SD), width = 0.2, position = position_dodge(0.5),  size = 0.9) + 
  scale_fill_manual(values=c("#F6B3AC","gray"))+
  labs(title = "(B4) Plasma, Dose=0.32 mg/kg/day", x = "Time (day)", y = "GenX Concentration in Plasma (μg/ml)", fill = "Data Type") 

p_Pla4 = p_Pla4 + 
  scale_y_continuous(limits = c(-0.008, 0.4),breaks=c(0.1,0.2,0.3),expand = c(0,0))+
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
library(gridExtra)
library(tiff)
library(grid)
tl <- rasterGrob(readTIFF("D:/zs/PBPK/2025GenX/Artwork/4tl.tif"))
ggsave("Figure 3-Female.tiff",scale = 1,
       plot = grid.arrange(p_Pla1, p_Pla2, p_Pla3, p_Pla4, tl, ncol = 3, nrow = 2),
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
formatted_r2 <- formatC(adjusted_r_squared, digits = 2, format = "f")
label_text <- paste0("italic(R)^{2} == '", formatted_r2, "'")

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
  scale_y_continuous(limits = c(-5,5), breaks = c(-4, 0, 4),labels = scales::math_format(10^.x))+
  scale_x_continuous(limits = c(-5,5),labels = scales::math_format(10^.x)) +
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
  annotate("text", x = -Inf, y = Inf, label = "(B)", hjust = 1.3, vjust = 0.8, size = 15, family = "Times New Roman", colour = "black") +
  annotate("text", x = -2, y = 3, label = label_text , parse = TRUE, size = 15, color = "black", family = "Times New Roman")
print(p.eva)
ggsave("Figure S4-female.tiff",scale = 1,
       plot = p.eva,
       path = "D:/zs/PBPK/2025GenX/Artwork",
       width = 24, height = 24, units = "cm",dpi=320)
