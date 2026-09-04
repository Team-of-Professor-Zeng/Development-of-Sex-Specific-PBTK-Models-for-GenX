
## Male ##########################################################################################
setwd("D:/zs/PBPK/1-方案/GenX/小鼠/code/Mice/Male")
evaluation <- read.csv(file = "Evalution-M.csv")

OBS.A_Pla <- evaluation %>% filter(Study == 1 & Sample == "Plasma") %>% select(Time = "Time", Plasma = "Conc", SD = "SD", Study = "Study", Dose = "Dose") # TDoses = 126
OBS.B_CL  <- evaluation %>% filter(Study == 2 & Sample == "Liver" ) %>% select(Time = "Time", CL     = "Conc", SD = "SD", Study = "Study", Dose = "Dose") # TDoses = 28
OBS.C_Pla <- evaluation %>% filter(Study == 3 & Sample == "Serum" ) %>% select(Time = "Time", Plasma = "Conc", SD = "SD", Study = "Study", Dose = "Dose") # TDoses = 140
OBS.C_CL  <- evaluation %>% filter(Study == 3 & Sample == "Liver" ) %>% select(Time = "Time", CL     = "Conc", SD = "SD", Study = "Study", Dose = "Dose") # TDoses = 140
OBS.D_Pla <- evaluation %>% filter(Study == 4 & Sample == "Serum" ) %>% select(Time = "Time", Plasma = "Conc", SD = "SD", Study = "Study", Dose = "Dose") # TDoses = 28
OBS.D_CL  <- evaluation %>% filter(Study == 4 & Sample == "Liver" ) %>% select(Time = "Time", CL     = "Conc", SD = "SD", Study = "Study", Dose = "Dose") # TDoses = 28
OBS.E_Pla <- evaluation %>% filter(Study == 5 & Sample == "Plasma") %>% select(Time = "Time", Plasma = "Conc", SD = "SD", Study = "Study", Dose = "Dose") # TDoses = 84
OBS.F_Pla <- evaluation %>% filter(Study == 6 & Sample == "Serum" ) %>% select(Time = "Time", Plasma = "Conc", SD = "SD", Study = "Study", Dose = "Dose") # TDoses = 28

OBS_Pla <- rbind.data.frame (OBS.A_Pla, OBS.C_Pla, OBS.D_Pla, OBS.E_Pla, OBS.F_Pla)
OBS_Pla <- OBS_Pla %>% mutate(Matrix = c("Obs.Plasma"))
OBS_CL  <- rbind.data.frame (OBS.B_CL, OBS.C_CL, OBS.D_CL)
OBS_CL  <- OBS_CL  %>% mutate(Matrix = c("Obs.Liver"))

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

##M1
pred.eva <- function(pars, N) { 
  pars <- exp(pars)                   
  BW          = 0.025                
  tinterval   = 24                   
  TDOSE1    = 126                    
  DOSE1     = 0.0032                 
  idata1 <- tibble(ID = 1:N) %>%  mutate(BW        = rnormTrunc  (N, min = 0.0103, max = 0.0397, mean = 0.025,      sd = 0.0075),   #CV=0.3
                                         VKC       = rnormTrunc  (N, min = 0.007 , max = 0.0270, mean = 0.017,      sd = 0.0051),   #CV=0.3
                                         GFRC      = rnormTrunc  (N, min = 25.586, max = 98.614, mean = 62.1 ,      sd = 18.63) ,   #CV=0.3
                                         Free                           = 0.0106,                 ### fitting parameters
                                         PL                             = 0.753,                  ### fitting parameters
                                         PLu                            = 2.0399,                 ### fitting parameters
                                         PRest                          = 0.14 ,                  ### fitting parameters
                                         K0C                            = 1.838,                  ### fitting parameters               
                                         Kabsc                          = 1.0461,                 ### fitting parameters 
                                         DOSEoral1 = DOSE1*BW)   
  DOSEoral1 = DOSE1*BW               
  ex.oral1  <- ev(ID = 1, amt = DOSEoral1, ii = tinterval, addl = TDOSE1 - 1, cmt = "AST", replicate = FALSE)
  ex.oral1_1  <- ev(ID = 1:N, amt = idata1$DOSEoral1, ii = tinterval, addl = TDOSE1 - 1, cmt = "AST", replicate = FALSE)
  tsamp1    = tgrid(0, tinterval*(TDOSE1 - 1) + 24*7, 1) ## set up the output time 
  
  out1 <- 
    mod %>%          
    param    (pars) %>%
    Req      (Plasma)%>%
    update   (atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral1, tgrid = tsamp1)  
  outdf1 <- cbind.data.frame(Time   = out1$time/24, Plasma = out1$Plasma) 
  
  out1_1 <- mod %>%data_set(ex.oral1_1) %>%
    idata_set(idata1) %>% 
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%
    mrgsim(obsonly=TRUE, tgrid = tsamp1)
  outdf1_1 = cbind.data.frame (ID     = out1_1$ID, Time   = out1_1$time/24, Plasma = out1_1$Plasma)
  
  return (list(outdf1, outdf1_1)) 
}
N = 1000

M1_m   <- pred.eva (pars.eva,N = N)[[1]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = Plasma) %>% mutate( Sex = "Male", Study= "M1", Dose = 0.0032,Tissue = "Plasma", Matrix = c("Pre.Plasma"))
M1_m1  <- M1_m %>% filter(Time == 126)    
M1     <- pred.eva (pars.eva,N = N)[[2]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = Plasma) 
M1_1   <- M1%>% filter(Time == 126)
M1_summary <- M1_1 %>% summarize(Mean = mean(Conc, na.rm = TRUE),
                                 SD = sd(Conc, na.rm = TRUE),
                                 Median = median(Conc, na.rm = TRUE),
                                 P25 = quantile(Conc, probs = 0.25, na.rm = TRUE),
                                 P75 = quantile(Conc, probs = 0.75, na.rm = TRUE),
                                 ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
                                 ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)) 
M1_summary <-M1_summary %>% mutate( Time = 126, Sex = "Male", Study= "M1", Dose = 0.0032, Tissue = "Plasma", Matrix = c("Pre.Plasma"))

##M2
pred.eva <- function(pars, N) { 
  pars <- exp(pars)                   
  BW          = 0.025                  
  tinterval   = 24                     
  TDOSE2    = 28                     
  DOSE2     = 100 
  idata2 <- tibble(ID = 1:N) %>%  mutate(BW        = rnormTrunc  (N, min = 0.0103, max = 0.0397, mean = 0.025,      sd = 0.0075),   #CV=0.3
                                         VKC       = rnormTrunc  (N, min = 0.007 , max = 0.0270, mean = 0.017,      sd = 0.0051),   #CV=0.3
                                         GFRC      = rnormTrunc  (N, min = 25.586, max = 98.614, mean = 62.1 ,      sd = 18.63) ,   #CV=0.3
                                         Free                           = 0.0106,                 ### fitting parameters
                                         PL                             = 0.753,                  ### fitting parameters
                                         PLu                            = 2.0399,                 ### fitting parameters
                                         PRest                          = 0.14 ,                  ### fitting parameters
                                         K0C                            = 1.838,                  ### fitting parameters               
                                         Kabsc                          = 1.0461,                 ### fitting parameters 
                                         DOSEoral2 = DOSE2*BW)   
  DOSEoral2 = DOSE2*BW               
  ex.oral2  <- ev(ID = 1, amt = DOSEoral2, ii = tinterval, addl = TDOSE2 - 1, cmt = "AST", replicate = FALSE)
  ex.oral2_1  <- ev(ID = 1:N, amt = idata2$DOSEoral2, ii = tinterval, addl = TDOSE2 - 1, cmt = "AST", replicate = FALSE)
  tsamp2    = tgrid(0, tinterval*(TDOSE2 - 1) + 24*7, 1) 
  
  out2 <- 
    mod %>%          
    param (pars) %>%
    Req   (Liver)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral2, tgrid = tsamp2)   
  outdf2 <- cbind.data.frame(Time = out2$time/24, CL = out2$Liver) 
  
  out2_1 <- mod %>%data_set(ex.oral2_1) %>%
    idata_set(idata2) %>% 
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%
    mrgsim(obsonly=TRUE, tgrid = tsamp2)
  outdf2_1 = cbind.data.frame (ID     = out2_1$ID, Time   = out2_1$time/24, CL = out2_1$Liver)
  return (list(outdf2, outdf2_1))
}
N = 1000

M2_m   <- pred.eva (pars.eva,N = N)[[1]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CL) %>% mutate( Sex = "Male", Study= "M2", Dose = 100,Tissue = "Liver", Matrix = c("Pre.Liver"))
M2_m1  <- M2_m  %>% filter(Time == 29)    
M2     <- pred.eva (pars.eva,N = N)[[2]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CL) 
M2_1   <- M2    %>% filter(Time == 29)
M2_summary <- M2_1 %>% summarize(Mean = mean(Conc, na.rm = TRUE),
                                 SD = sd(Conc, na.rm = TRUE),
                                 Median = median(Conc, na.rm = TRUE),
                                 P25 = quantile(Conc, probs = 0.25, na.rm = TRUE),
                                 P75 = quantile(Conc, probs = 0.75, na.rm = TRUE),
                                 ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
                                 ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)) 
M2_summary <-M2_summary %>% mutate( Time = 29, Sex = "Male", Study= "M2", Dose = 100, Tissue = "Liver", Matrix = c("Pre.Liver"))

##M3
pred.eva <- function(pars, N) { 
  pars <- exp(pars)                   
  BW          = 0.025                  
  tinterval   = 24                     
  TDOSE3      = 140                    
  tsamp3      = tgrid(0, tinterval*(TDOSE3 - 1) + 24*7, 1)
  
  #a
  DOSE3a      = 0.0004                 
  DOSEoral3a  = DOSE3a*BW 
  idata3a <- tibble(ID = 1:N) %>%  mutate(BW        = rnormTrunc  (N, min = 0.0103, max = 0.0397, mean = 0.025,      sd = 0.0075),   #CV=0.3
                                         VKC       = rnormTrunc  (N, min = 0.007 , max = 0.0270, mean = 0.017,      sd = 0.0051),   #CV=0.3
                                         GFRC      = rnormTrunc  (N, min = 25.586, max = 98.614, mean = 62.1 ,      sd = 18.63) ,   #CV=0.3
                                         Free                           = 0.0106,                 ### fitting parameters
                                         PL                             = 0.753,                  ### fitting parameters
                                         PLu                            = 2.0399,                 ### fitting parameters
                                         PRest                          = 0.14 ,                  ### fitting parameters
                                         K0C                            = 1.838,                  ### fitting parameters               
                                         Kabsc                          = 1.0461,                 ### fitting parameters 
                                         DOSEoral3a  = DOSE3a*BW)  
  ex.oral3a   <- ev(ID = 1, amt = DOSEoral3a, ii = tinterval, addl = TDOSE3 - 1, cmt = "AST", replicate = FALSE)
  ex.oral3a_1  <- ev(ID = 1:N, amt = idata3a$DOSEoral3a, ii = tinterval, addl = TDOSE3 - 1, cmt = "AST", replicate = FALSE)
  out3a <- 
    mod %>%          
    param (pars) %>%
    Req   (Plasma, Liver)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral3a, tgrid = tsamp3) 
  outdf3a <- cbind.data.frame(Time   = out3a$time/24,  Plasma = out3a$Plasma, CL = out3a$Liver) 
  
  out3a_1 <- mod %>%data_set(ex.oral3a_1) %>%
    idata_set(idata3a) %>% 
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%
    mrgsim(obsonly=TRUE, tgrid = tsamp3)
  outdf3a_1 = cbind.data.frame (ID     = out3a_1$ID, Time   = out3a_1$time/24, Plasma = out3a_1$Plasma, CL = out3a_1$Liver)
  
  #b
  DOSE3b      = 0.002                 
  DOSEoral3b  = DOSE3b*BW  
  idata3b <- tibble(ID = 1:N) %>%  mutate(BW        = rnormTrunc  (N, min = 0.0103, max = 0.0397, mean = 0.025,      sd = 0.0075),   #CV=0.3
                                          VKC       = rnormTrunc  (N, min = 0.007 , max = 0.0270, mean = 0.017,      sd = 0.0051),   #CV=0.3
                                          GFRC      = rnormTrunc  (N, min = 25.586, max = 98.614, mean = 62.1 ,      sd = 18.63) ,   #CV=0.3
                                          Free                           = 0.0106,                 ### fitting parameters
                                          PL                             = 0.753,                  ### fitting parameters
                                          PLu                            = 2.0399,                 ### fitting parameters
                                          PRest                          = 0.14 ,                  ### fitting parameters
                                          K0C                            = 1.838,                  ### fitting parameters               
                                          Kabsc                          = 1.0461,                 ### fitting parameters 
                                          DOSEoral3b  = DOSE3b*BW)  
  ex.oral3b   <- ev(ID = 1, amt = DOSEoral3b, ii = tinterval, addl = TDOSE3 - 1, cmt = "AST", replicate = FALSE)
  ex.oral3b_1  <- ev(ID = 1:N, amt = idata3b$DOSEoral3b, ii = tinterval, addl = TDOSE3 - 1, cmt = "AST", replicate = FALSE)
  out3b <- 
    mod %>%          
    param (pars) %>%
    Req   (Plasma, Liver)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral3b, tgrid = tsamp3)   
  outdf3b <- cbind.data.frame(Time   = out3b$time/24,  Plasma = out3b$Plasma, CL     = out3b$Liver)
  out3b_1 <- mod %>%data_set(ex.oral3b_1) %>%
    idata_set(idata3b) %>% 
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%
    mrgsim(obsonly=TRUE, tgrid = tsamp3)
  outdf3b_1 = cbind.data.frame (ID     = out3b_1$ID, Time   = out3b_1$time/24, Plasma = out3b_1$Plasma, CL = out3b_1$Liver)
  
  #c
  DOSE3c      = 0.01                 
  DOSEoral3c  = DOSE3c*BW 
  idata3c <- tibble(ID = 1:N) %>%  mutate(BW        = rnormTrunc  (N, min = 0.0103, max = 0.0397, mean = 0.025,      sd = 0.0075),   #CV=0.3
                                          VKC       = rnormTrunc  (N, min = 0.007 , max = 0.0270, mean = 0.017,      sd = 0.0051),   #CV=0.3
                                          GFRC      = rnormTrunc  (N, min = 25.586, max = 98.614, mean = 62.1 ,      sd = 18.63) ,   #CV=0.3
                                          Free                           = 0.0106,                 ### fitting parameters
                                          PL                             = 0.753,                  ### fitting parameters
                                          PLu                            = 2.0399,                 ### fitting parameters
                                          PRest                          = 0.14 ,                  ### fitting parameters
                                          K0C                            = 1.838,                  ### fitting parameters               
                                          Kabsc                          = 1.0461,                 ### fitting parameters 
                                          DOSEoral3c  = DOSE3c*BW)  
  ex.oral3c   <- ev(ID = 1, amt = DOSEoral3c, ii = tinterval, addl = TDOSE3 - 1, cmt = "AST", replicate = FALSE)
  ex.oral3c_1  <- ev(ID = 1:N, amt = idata3c$DOSEoral3c, ii = tinterval, addl = TDOSE3 - 1, cmt = "AST", replicate = FALSE)
  out3c <- 
    mod %>%          
    param (pars) %>%
    Req   (Plasma, Liver)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral3c, tgrid = tsamp3)  
  outdf3c <- cbind.data.frame(Time   = out3c$time/24, Plasma = out3c$Plasma, CL     = out3c$Liver)
  out3c_1 <- mod %>%data_set(ex.oral3c_1) %>%
    idata_set(idata3c) %>% 
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%
    mrgsim(obsonly=TRUE, tgrid = tsamp3)
  outdf3c_1 = cbind.data.frame (ID     = out3c_1$ID, Time   = out3c_1$time/24, Plasma = out3c_1$Plasma, CL = out3c_1$Liver)
  
  return (list(outdf3a, outdf3a_1, outdf3b, outdf3b_1, outdf3c, outdf3c_1))
}
N = 1000

##a
M3a_m_Pla   <- pred.eva (pars.eva,N = N)[[1]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = Plasma)  %>% mutate(  Sex = "Male", Study= "M3", Dose = 0.0004, Tissue = "Plasma", Matrix = c("Pre.Plasma"))
M3a_m_CL    <- pred.eva (pars.eva,N = N)[[1]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CL)  %>% mutate(  Sex = "Male", Study= "M3", Dose = 0.0004, Tissue = "Liver", Matrix = c("Pre.Liver"))
M3a_m       <- rbind(M3a_m_Pla, M3a_m_CL)
M3a_m1      <- M3a_m  %>% filter(Time == 140)

M3a_Pla     <- pred.eva (pars.eva,N = N)[[2]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = Plasma)
M3a_1_Pla   <- M3a_Pla %>% filter(Time == 140)
M3a_summary_Pla <- M3a_1_Pla %>% summarize(Mean = mean(Conc, na.rm = TRUE),
                                 SD = sd(Conc, na.rm = TRUE),
                                 Median = median(Conc, na.rm = TRUE),
                                 P25 = quantile(Conc, probs = 0.25, na.rm = TRUE),
                                 P75 = quantile(Conc, probs = 0.75, na.rm = TRUE),
                                 ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
                                 ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)) 
M3a_summary_Pla <-M3a_summary_Pla %>% mutate( Time = 140, Sex = "Male", Study= "M3", Dose = 0.0004, Tissue = "Plasma", Matrix = c("Pre.Plasma"))
M3a_CL      <- pred.eva (pars.eva,N = N)[[2]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CL)  
M3a_1_CL    <- M3a_CL %>% filter(Time == 140)
M3a_summary_CL <- M3a_1_CL %>% summarize(Mean = mean(Conc, na.rm = TRUE),
                                           SD = sd(Conc, na.rm = TRUE),
                                           Median = median(Conc, na.rm = TRUE),
                                           P25 = quantile(Conc, probs = 0.25, na.rm = TRUE),
                                           P75 = quantile(Conc, probs = 0.75, na.rm = TRUE),
                                           ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
                                           ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)) 
M3a_summary_CL <-M3a_summary_CL %>% mutate( Time = 140, Sex = "Male", Study= "M3", Dose = 0.0004, Tissue = "Liver", Matrix = c("Pre.Liver"))
M3a        <- rbind(M3a_summary_Pla, M3a_summary_CL)

##b
M3b_m_Pla   <- pred.eva (pars.eva,N = N)[[3]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = Plasma)  %>% mutate( Sex = "Male", Study= "M3", Dose = 0.002, Tissue = "Plasma", Matrix = c("Pre.Plasma"))
M3b_m_CL    <- pred.eva (pars.eva,N = N)[[3]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CL)  %>% mutate( Sex = "Male", Study= "M3", Dose = 0.002, Tissue = "Liver", Matrix = c("Pre.Liver"))
M3b_m       <- rbind(M3b_m_Pla, M3b_m_CL)
M3b_m1      <- M3b_m  %>% filter(Time == 140) 

M3b_Pla     <- pred.eva (pars.eva,N = N)[[4]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = Plasma)
M3b_1_Pla   <- M3b_Pla %>% filter(Time == 140)
M3b_summary_Pla <- M3b_1_Pla %>% summarize(Mean = mean(Conc, na.rm = TRUE),
                                           SD = sd(Conc, na.rm = TRUE),
                                           Median = median(Conc, na.rm = TRUE),
                                           P25 = quantile(Conc, probs = 0.25, na.rm = TRUE),
                                           P75 = quantile(Conc, probs = 0.75, na.rm = TRUE),
                                           ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
                                           ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)) 
M3b_summary_Pla <-M3b_summary_Pla %>% mutate( Time = 140, Sex = "Male", Study= "M3", Dose = 0.002, Tissue = "Plasma", Matrix = c("Pre.Plasma"))
M3b_CL      <- pred.eva (pars.eva,N = N)[[4]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CL)  
M3b_1_CL    <- M3b_CL %>% filter(Time == 140)
M3b_summary_CL <- M3b_1_CL %>% summarize(Mean = mean(Conc, na.rm = TRUE),
                                         SD = sd(Conc, na.rm = TRUE),
                                         Median = median(Conc, na.rm = TRUE),
                                         P25 = quantile(Conc, probs = 0.25, na.rm = TRUE),
                                         P75 = quantile(Conc, probs = 0.75, na.rm = TRUE),
                                         ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
                                         ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)) 
M3b_summary_CL <-M3b_summary_CL %>% mutate( Time = 140, Sex = "Male", Study= "M3", Dose = 0.002, Tissue = "Liver", Matrix = c("Pre.Liver"))
M3b        <- rbind(M3b_summary_Pla, M3b_summary_CL)

##c
M3c_m_Pla   <- pred.eva (pars.eva,N = N)[[5]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = Plasma)  %>% mutate(Sex = "Male", Study= "M3", Dose = 0.01, Tissue = "Plasma", Matrix = c("Pre.Plasma"))
M3c_m_CL    <- pred.eva (pars.eva,N = N)[[5]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CL)  %>% mutate(Sex = "Male", Study= "M3", Dose = 0.01, Tissue = "Liver", Matrix = c("Pre.Liver"))
M3c_m       <- rbind(M3c_m_Pla, M3c_m_CL)
M3c_m1      <- M3c_m  %>% filter(Time == 140) 

M3c_Pla     <- pred.eva (pars.eva,N = N)[[6]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = Plasma)
M3c_1_Pla   <- M3c_Pla %>% filter(Time == 140)
M3c_summary_Pla <- M3c_1_Pla %>% summarize(Mean = mean(Conc, na.rm = TRUE),
                                           SD = sd(Conc, na.rm = TRUE),
                                           Median = median(Conc, na.rm = TRUE),
                                           P25 = quantile(Conc, probs = 0.25, na.rm = TRUE),
                                           P75 = quantile(Conc, probs = 0.75, na.rm = TRUE),
                                           ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
                                           ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)) 
M3c_summary_Pla <-M3c_summary_Pla %>% mutate( Time = 140, Sex = "Male", Study= "M3", Dose = 0.01, Tissue = "Plasma", Matrix = c("Pre.Plasma"))
M3c_CL      <- pred.eva (pars.eva,N = N)[[6]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CL)  
M3c_1_CL    <- M3c_CL %>% filter(Time == 140)
M3c_summary_CL <- M3c_1_CL %>% summarize(Mean = mean(Conc, na.rm = TRUE),
                                         SD = sd(Conc, na.rm = TRUE),
                                         Median = median(Conc, na.rm = TRUE),
                                         P25 = quantile(Conc, probs = 0.25, na.rm = TRUE),
                                         P75 = quantile(Conc, probs = 0.75, na.rm = TRUE),
                                         ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
                                         ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)) 
M3c_summary_CL <-M3c_summary_CL %>% mutate( Time = 140, Sex = "Male", Study= "M3", Dose = 0.01, Tissue = "Liver", Matrix = c("Pre.Liver"))
M3c        <- rbind(M3c_summary_Pla, M3c_summary_CL)

M3_m1 <- rbind(M3a_m1, M3b_m1, M3c_m1)
M3_summary <- rbind(M3a, M3b, M3c)

##M4
pred.eva <- function(pars, N) { 
  pars <- exp(pars)                   
  BW          = 0.025                  
  tinterval   = 24                     
  TDOSE4      = 28                     
  tsamp4      = tgrid(0, tinterval*(TDOSE4 - 1) + 24*7, 1)
  
  #a
  DOSE4a      = 0.4                 
  DOSEoral4a  = DOSE4a*BW  
  idata4a <- tibble(ID = 1:N) %>%  mutate(BW        = rnormTrunc  (N, min = 0.0103, max = 0.0397, mean = 0.025,      sd = 0.0075),   #CV=0.3
                                          VKC       = rnormTrunc  (N, min = 0.007 , max = 0.0270, mean = 0.017,      sd = 0.0051),   #CV=0.3
                                          GFRC      = rnormTrunc  (N, min = 25.586, max = 98.614, mean = 62.1 ,      sd = 18.63) ,   #CV=0.3
                                          Free                           = 0.0106,                 ### fitting parameters
                                          PL                             = 0.753,                  ### fitting parameters
                                          PLu                            = 2.0399,                 ### fitting parameters
                                          PRest                          = 0.14 ,                  ### fitting parameters
                                          K0C                            = 1.838,                  ### fitting parameters               
                                          Kabsc                          = 1.0461,                 ### fitting parameters 
                                          DOSEoral4a  = DOSE4a*BW)  
  ex.oral4a   <- ev(ID = 1, amt = DOSEoral4a, ii = tinterval, addl = TDOSE4 - 1, cmt = "AST", replicate = FALSE)
  ex.oral4a_1  <- ev(ID = 1:N, amt = idata4a$DOSEoral4a, ii = tinterval, addl = TDOSE4 - 1, cmt = "AST", replicate = FALSE)
  out4a <- 
    mod %>%          
    param (pars) %>%
    Req   (Plasma, Liver)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral4a, tgrid = tsamp4)
  outdf4a <- cbind.data.frame(Time   = out4a$time/24, Plasma = out4a$Plasma, CL     = out4a$Liver) 
  
  out4a_1 <- mod %>%data_set(ex.oral4a_1) %>%
    idata_set(idata4a) %>% 
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%
    mrgsim(obsonly=TRUE, tgrid = tsamp4)
  outdf4a_1 = cbind.data.frame (ID     = out4a_1$ID, Time   = out4a_1$time/24, Plasma = out4a_1$Plasma, CL = out4a_1$Liver)
  
  #b
  DOSE4b      = 2                 
  DOSEoral4b  = DOSE4b*BW            
  idata4b <- tibble(ID = 1:N) %>%  mutate(BW        = rnormTrunc  (N, min = 0.0103, max = 0.0397, mean = 0.025,      sd = 0.0075),   #CV=0.3
                                          VKC       = rnormTrunc  (N, min = 0.007 , max = 0.0270, mean = 0.017,      sd = 0.0051),   #CV=0.3
                                          GFRC      = rnormTrunc  (N, min = 25.586, max = 98.614, mean = 62.1 ,      sd = 18.63) ,   #CV=0.3
                                          Free                           = 0.0106,                 ### fitting parameters
                                          PL                             = 0.753,                  ### fitting parameters
                                          PLu                            = 2.0399,                 ### fitting parameters
                                          PRest                          = 0.14 ,                  ### fitting parameters
                                          K0C                            = 1.838,                  ### fitting parameters               
                                          Kabsc                          = 1.0461,                 ### fitting parameters 
                                          DOSEoral4b  = DOSE4b*BW)  
  ex.oral4b   <- ev(ID = 1, amt = DOSEoral4b, ii = tinterval, addl = TDOSE4 - 1, cmt = "AST", replicate = FALSE)
  ex.oral4b_1  <- ev(ID = 1:N, amt = idata4b$DOSEoral4b, ii = tinterval, addl = TDOSE4 - 1, cmt = "AST", replicate = FALSE)
  out4b <- 
    mod %>%          
    param (pars) %>%
    Req   (Plasma, Liver)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral4b, tgrid = tsamp4)  
  outdf4b <- cbind.data.frame(Time   = out4b$time/24,  Plasma = out4b$Plasma, CL     = out4b$Liver)
  
  out4b_1 <- mod %>%data_set(ex.oral4b_1) %>%
    idata_set(idata4b) %>% 
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%
    mrgsim(obsonly=TRUE, tgrid = tsamp4)
  outdf4b_1 = cbind.data.frame (ID     = out4b_1$ID, Time   = out4b_1$time/24, Plasma = out4b_1$Plasma, CL = out4b_1$Liver)
  
  #c
  DOSE4c      = 10                 
  DOSEoral4c  = DOSE4c*BW         
  idata4c <- tibble(ID = 1:N) %>%  mutate(BW        = rnormTrunc  (N, min = 0.0103, max = 0.0397, mean = 0.025,      sd = 0.0075),   #CV=0.3
                                          VKC       = rnormTrunc  (N, min = 0.007 , max = 0.0270, mean = 0.017,      sd = 0.0051),   #CV=0.3
                                          GFRC      = rnormTrunc  (N, min = 25.586, max = 98.614, mean = 62.1 ,      sd = 18.63) ,   #CV=0.3
                                          Free                           = 0.0106,                 ### fitting parameters
                                          PL                             = 0.753,                  ### fitting parameters
                                          PLu                            = 2.0399,                 ### fitting parameters
                                          PRest                          = 0.14 ,                  ### fitting parameters
                                          K0C                            = 1.838,                  ### fitting parameters               
                                          Kabsc                          = 1.0461,                 ### fitting parameters 
                                          DOSEoral4c  = DOSE4c*BW)  
  ex.oral4c   <- ev(ID = 1, amt = DOSEoral4c, ii = tinterval, addl = TDOSE4 - 1, cmt = "AST", replicate = FALSE)
  ex.oral4c_1  <- ev(ID = 1:N, amt = idata4c$DOSEoral4c, ii = tinterval, addl = TDOSE4 - 1, cmt = "AST", replicate = FALSE)
  out4c <- 
    mod %>%          
    param (pars) %>%
    Req   (Plasma, Liver)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral4c, tgrid = tsamp4)  
  outdf4c <- cbind.data.frame(Time   = out4c$time/24, Plasma = out4c$Plasma, CL     = out4c$Liver)
  
  out4c_1 <- mod %>%data_set(ex.oral4c_1) %>%
    idata_set(idata4c) %>% 
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%
    mrgsim(obsonly=TRUE, tgrid = tsamp4)
  outdf4c_1 = cbind.data.frame (ID     = out4c_1$ID, Time   = out4c_1$time/24, Plasma = out4c_1$Plasma, CL = out4c_1$Liver)
  
  return (list(outdf4a, outdf4a_1, outdf4b, outdf4b_1, outdf4c, outdf4c_1))
}
N = 1000

##a
M4a_m_Pla   <- pred.eva (pars.eva,N = N)[[1]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = Plasma)  %>% mutate(  Sex = "Male", Study= "M4", Dose = 0.4, Tissue = "Plasma", Matrix = c("Pre.Plasma"))
M4a_m_CL    <- pred.eva (pars.eva,N = N)[[1]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CL)  %>% mutate(  Sex = "Male", Study= "M4", Dose = 0.4, Tissue = "Liver", Matrix = c("Pre.Liver"))
M4a_m       <- rbind(M4a_m_Pla, M4a_m_CL)
M4a_m1      <- M4a_m  %>% filter(Time == 28)

M4a_Pla     <- pred.eva (pars.eva,N = N)[[2]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = Plasma)
M4a_1_Pla   <- M4a_Pla %>% filter(Time == 28)
M4a_summary_Pla <- M4a_1_Pla %>% summarize(Mean = mean(Conc, na.rm = TRUE),
                                           SD = sd(Conc, na.rm = TRUE),
                                           Median = median(Conc, na.rm = TRUE),
                                           P25 = quantile(Conc, probs = 0.25, na.rm = TRUE),
                                           P75 = quantile(Conc, probs = 0.75, na.rm = TRUE),
                                           ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
                                           ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)) 
M4a_summary_Pla <-M4a_summary_Pla %>% mutate( Time = 28, Sex = "Male", Study= "M4", Dose = 0.4, Tissue = "Plasma", Matrix = c("Pre.Plasma"))
M4a_CL      <- pred.eva (pars.eva,N = N)[[2]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CL)  
M4a_1_CL    <- M4a_CL %>% filter(Time == 28)
M4a_summary_CL <- M4a_1_CL %>% summarize(Mean = mean(Conc, na.rm = TRUE),
                                         SD = sd(Conc, na.rm = TRUE),
                                         Median = median(Conc, na.rm = TRUE),
                                         P25 = quantile(Conc, probs = 0.25, na.rm = TRUE),
                                         P75 = quantile(Conc, probs = 0.75, na.rm = TRUE),
                                         ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
                                         ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)) 
M4a_summary_CL <-M4a_summary_CL %>% mutate( Time = 28, Sex = "Male", Study= "M4", Dose = 0.4, Tissue = "Liver", Matrix = c("Pre.Liver"))
M4a        <- rbind(M4a_summary_Pla, M4a_summary_CL)

##b
M4b_m_Pla   <- pred.eva (pars.eva,N = N)[[3]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = Plasma)  %>% mutate( Sex = "Male", Study= "M4", Dose = 2, Tissue = "Plasma", Matrix = c("Pre.Plasma"))
M4b_m_CL    <- pred.eva (pars.eva,N = N)[[3]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CL)  %>% mutate( Sex = "Male", Study= "M4", Dose = 2, Tissue = "Liver", Matrix = c("Pre.Liver"))
M4b_m       <- rbind(M4b_m_Pla, M4b_m_CL)
M4b_m1      <- M4b_m  %>% filter(Time == 28) 

M4b_Pla     <- pred.eva (pars.eva,N = N)[[4]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = Plasma)
M4b_1_Pla   <- M4b_Pla %>% filter(Time == 28)
M4b_summary_Pla <- M4b_1_Pla %>% summarize(Mean = mean(Conc, na.rm = TRUE),
                                           SD = sd(Conc, na.rm = TRUE),
                                           Median = median(Conc, na.rm = TRUE),
                                           P25 = quantile(Conc, probs = 0.25, na.rm = TRUE),
                                           P75 = quantile(Conc, probs = 0.75, na.rm = TRUE),
                                           ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
                                           ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)) 
M4b_summary_Pla <- M4b_summary_Pla %>% mutate( Time = 28, Sex = "Male", Study= "M4", Dose = 2, Tissue = "Plasma", Matrix = c("Pre.Plasma"))
M4b_CL      <- pred.eva (pars.eva,N = N)[[4]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CL)  
M4b_1_CL    <- M4b_CL %>% filter(Time == 28)
M4b_summary_CL <- M4b_1_CL %>% summarize(Mean = mean(Conc, na.rm = TRUE),
                                         SD = sd(Conc, na.rm = TRUE),
                                         Median = median(Conc, na.rm = TRUE),
                                         P25 = quantile(Conc, probs = 0.25, na.rm = TRUE),
                                         P75 = quantile(Conc, probs = 0.75, na.rm = TRUE),
                                         ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
                                         ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)) 
M4b_summary_CL <- M4b_summary_CL %>% mutate( Time = 28, Sex = "Male", Study= "M4", Dose = 2, Tissue = "Liver", Matrix = c("Pre.Liver"))
M4b        <- rbind(M4b_summary_Pla, M4b_summary_CL)

##c
M4c_m_Pla   <- pred.eva (pars.eva,N = N)[[5]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = Plasma)  %>% mutate(Sex = "Male", Study= "M4", Dose = 10, Tissue = "Plasma", Matrix = c("Pre.Plasma"))
M4c_m_CL    <- pred.eva (pars.eva,N = N)[[5]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CL)  %>% mutate(Sex = "Male", Study= "M4", Dose = 10, Tissue = "Liver", Matrix = c("Pre.Liver"))
M4c_m       <- rbind(M4c_m_Pla, M4c_m_CL)
M4c_m1      <- M4c_m  %>% filter(Time == 28) 

M4c_Pla     <- pred.eva (pars.eva,N = N)[[6]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = Plasma)
M4c_1_Pla   <- M4c_Pla %>% filter(Time == 28)
M4c_summary_Pla <- M4c_1_Pla %>% summarize(Mean = mean(Conc, na.rm = TRUE),
                                           SD = sd(Conc, na.rm = TRUE),
                                           Median = median(Conc, na.rm = TRUE),
                                           P25 = quantile(Conc, probs = 0.25, na.rm = TRUE),
                                           P75 = quantile(Conc, probs = 0.75, na.rm = TRUE),
                                           ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
                                           ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)) 
M4c_summary_Pla <- M4c_summary_Pla %>% mutate( Time = 28, Sex = "Male", Study= "M4", Dose = 10, Tissue = "Plasma", Matrix = c("Pre.Plasma"))
M4c_CL      <- pred.eva (pars.eva,N = N)[[6]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CL)  
M4c_1_CL    <- M4c_CL %>% filter(Time == 28)
M4c_summary_CL <- M4c_1_CL %>% summarize(Mean = mean(Conc, na.rm = TRUE),
                                         SD = sd(Conc, na.rm = TRUE),
                                         Median = median(Conc, na.rm = TRUE),
                                         P25 = quantile(Conc, probs = 0.25, na.rm = TRUE),
                                         P75 = quantile(Conc, probs = 0.75, na.rm = TRUE),
                                         ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
                                         ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)) 
M4c_summary_CL <- M4c_summary_CL %>% mutate( Time = 28, Sex = "Male", Study= "M4", Dose = 10, Tissue = "Liver", Matrix = c("Pre.Liver"))
M4c        <- rbind(M4c_summary_Pla, M4c_summary_CL)

M4_m1 <- rbind(M4a_m1, M4b_m1, M4c_m1)
M4_summary <- rbind(M4a, M4b, M4c)

##M5
pred.eva <- function(pars, N) { 
  pars <- exp(pars)                   
  BW          = 0.025                
  tinterval   = 24                   
  TDOSE5      = 84                     
  DOSE5       = 0.32   
  DOSEoral5   = DOSE5*BW
  idata5 <- tibble(ID = 1:N) %>%  mutate(BW        = rnormTrunc  (N, min = 0.0103, max = 0.0397, mean = 0.025,      sd = 0.0075),   #CV=0.3
                                         VKC       = rnormTrunc  (N, min = 0.007 , max = 0.0270, mean = 0.017,      sd = 0.0051),   #CV=0.3
                                         GFRC      = rnormTrunc  (N, min = 25.586, max = 98.614, mean = 62.1 ,      sd = 18.63) ,   #CV=0.3
                                         Free                           = 0.0106,                 ### fitting parameters
                                         PL                             = 0.753,                  ### fitting parameters
                                         PLu                            = 2.0399,                 ### fitting parameters
                                         PRest                          = 0.14 ,                  ### fitting parameters
                                         K0C                            = 1.838,                  ### fitting parameters               
                                         Kabsc                          = 1.0461,                 ### fitting parameters 
                                         DOSEoral5   = DOSE5*BW)   
  ex.oral5    <- ev(ID = 1, amt = DOSEoral5, ii = tinterval, addl = TDOSE5 - 1, cmt = "AST", replicate = FALSE)
  ex.oral5_1  <- ev(ID = 1:N, amt = idata5$DOSEoral5, ii = tinterval, addl = TDOSE5 - 1, cmt = "AST", replicate = FALSE)
  tsamp5      = tgrid(0, tinterval*(TDOSE5 - 1) + 24*7, 1) 
  
  out5 <- 
    mod %>%          
    param (pars) %>%
    Req   (Plasma)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral5, tgrid = tsamp5)   
  outdf5 <- cbind.data.frame(Time = out5$time/24, Plasma = out5$Plasma) 
  
  out5_1 <- mod %>%data_set(ex.oral5_1) %>%
    idata_set(idata5) %>% 
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%
    mrgsim(obsonly=TRUE, tgrid = tsamp5)
  outdf5_1 = cbind.data.frame (ID     = out5_1$ID, Time   = out5_1$time/24, Plasma = out5_1$Plasma)
  
  return (list(outdf5, outdf5_1)) 
}
N = 1000

M5_m   <- pred.eva (pars.eva,N = N)[[1]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = Plasma) %>% mutate( Sex = "Male", Study= "M5", Dose = 0.32,Tissue = "Plasma", Matrix = c("Pre.Plasma"))
M5_m1  <- M5_m %>% filter(Time == 84)    
M5     <- pred.eva (pars.eva,N = N)[[2]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = Plasma) 
M5_1   <- M5 %>% filter(Time == 84)
M5_summary <- M5_1 %>% summarize(Mean = mean(Conc, na.rm = TRUE),
                                 SD = sd(Conc, na.rm = TRUE),
                                 Median = median(Conc, na.rm = TRUE),
                                 P25 = quantile(Conc, probs = 0.25, na.rm = TRUE),
                                 P75 = quantile(Conc, probs = 0.75, na.rm = TRUE),
                                 ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
                                 ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)) 
M5_summary <-M5_summary %>% mutate( Time = 84, Sex = "Male", Study= "M5", Dose = 0.32, Tissue = "Plasma", Matrix = c("Pre.Plasma"))

##M6
pred.eva <- function(pars, N) { 
  pars <- exp(pars)                   
  BW          = 0.025                  
  tinterval   = 24                     
  TDOSE6      = 28
  tsamp6      = tgrid(0, tinterval*(TDOSE6 - 1) + 24*7, 1)
  
  #a
  DOSE6a      = 1                 
  DOSEoral6a  = DOSE6a*BW  
  idata6a <- tibble(ID = 1:N) %>%  mutate(BW        = rnormTrunc  (N, min = 0.0103, max = 0.0397, mean = 0.025,      sd = 0.0075),   #CV=0.3
                                          VKC       = rnormTrunc  (N, min = 0.007 , max = 0.0270, mean = 0.017,      sd = 0.0051),   #CV=0.3
                                          GFRC      = rnormTrunc  (N, min = 25.586, max = 98.614, mean = 62.1 ,      sd = 18.63) ,   #CV=0.3
                                          Free                           = 0.0106,                 ### fitting parameters
                                          PL                             = 0.753,                  ### fitting parameters
                                          PLu                            = 2.0399,                 ### fitting parameters
                                          PRest                          = 0.14 ,                  ### fitting parameters
                                          K0C                            = 1.838,                  ### fitting parameters               
                                          Kabsc                          = 1.0461,                 ### fitting parameters 
                                          DOSEoral6a  = DOSE6a*BW)  
  ex.oral6a   <- ev(ID = 1, amt = DOSEoral6a, ii = tinterval, addl = TDOSE6 - 1, cmt = "AST", replicate = FALSE)
  ex.oral6a_1  <- ev(ID = 1:N, amt = idata6a$DOSEoral6a, ii = tinterval, addl = TDOSE6 - 1, cmt = "AST", replicate = FALSE)
  out6a <- 
    mod %>%          
    param (pars) %>%
    Req   (Plasma)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral6a, tgrid = tsamp6)   
  outdf6a <- cbind.data.frame(Time   = out6a$time/24, Plasma = out6a$Plasma) 
  
  out6a_1 <- mod %>%data_set(ex.oral6a_1) %>%
    idata_set(idata6a) %>% 
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%
    mrgsim(obsonly=TRUE, tgrid = tsamp6)
  outdf6a_1 = cbind.data.frame (ID     = out6a_1$ID, Time   = out6a_1$time/24, Plasma = out6a_1$Plasma)
  
  #b
  DOSE6b      = 10                 
  DOSEoral6b  = DOSE6b*BW           
  idata6b <- tibble(ID = 1:N) %>%  mutate(BW        = rnormTrunc  (N, min = 0.0103, max = 0.0397, mean = 0.025,      sd = 0.0075),   #CV=0.3
                                          VKC       = rnormTrunc  (N, min = 0.007 , max = 0.0270, mean = 0.017,      sd = 0.0051),   #CV=0.3
                                          GFRC      = rnormTrunc  (N, min = 25.586, max = 98.614, mean = 62.1 ,      sd = 18.63) ,   #CV=0.3
                                          Free                           = 0.0106,                 ### fitting parameters
                                          PL                             = 0.753,                  ### fitting parameters
                                          PLu                            = 2.0399,                 ### fitting parameters
                                          PRest                          = 0.14 ,                  ### fitting parameters
                                          K0C                            = 1.838,                  ### fitting parameters               
                                          Kabsc                          = 1.0461,                 ### fitting parameters 
                                          DOSEoral6b  = DOSE6b*BW)  
  ex.oral6b   <- ev(ID = 1, amt = DOSEoral6b, ii = tinterval, addl = TDOSE6 - 1, cmt = "AST", replicate = FALSE)
  ex.oral6b_1  <- ev(ID = 1:N, amt = idata6b$DOSEoral6b, ii = tinterval, addl = TDOSE6 - 1, cmt = "AST", replicate = FALSE)
  out6b <- 
    mod %>%          
    param (pars) %>%
    Req   (Plasma)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral6b, tgrid = tsamp6)   
  outdf6b <- cbind.data.frame(Time   = out6b$time/24, Plasma = out6b$Plasma) 
  
  out6b_1 <- mod %>%data_set(ex.oral6b_1) %>%
    idata_set(idata6b) %>% 
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%
    mrgsim(obsonly=TRUE, tgrid = tsamp6)
  outdf6b_1 = cbind.data.frame (ID     = out6b_1$ID, Time   = out6b_1$time/24, Plasma = out6b_1$Plasma)
  
  #c
  DOSE6c      = 100                 
  DOSEoral6c  = DOSE6c*BW         
  idata6c <- tibble(ID = 1:N) %>%  mutate(BW        = rnormTrunc  (N, min = 0.0103, max = 0.0397, mean = 0.025,      sd = 0.0075),   #CV=0.3
                                          VKC       = rnormTrunc  (N, min = 0.007 , max = 0.0270, mean = 0.017,      sd = 0.0051),   #CV=0.3
                                          GFRC      = rnormTrunc  (N, min = 25.586, max = 98.614, mean = 62.1 ,      sd = 18.63) ,   #CV=0.3
                                          Free                           = 0.0106,                 ### fitting parameters
                                          PL                             = 0.753,                  ### fitting parameters
                                          PLu                            = 2.0399,                 ### fitting parameters
                                          PRest                          = 0.14 ,                  ### fitting parameters
                                          K0C                            = 1.838,                  ### fitting parameters               
                                          Kabsc                          = 1.0461,                 ### fitting parameters 
                                          DOSEoral6c  = DOSE6c*BW)  
  ex.oral6c   <- ev(ID = 1, amt = DOSEoral6c, ii = tinterval, addl = TDOSE6 - 1, cmt = "AST", replicate = FALSE)
  ex.oral6c_1  <- ev(ID = 1:N, amt = idata6c$DOSEoral6c, ii = tinterval, addl = TDOSE6 - 1, cmt = "AST", replicate = FALSE)
  out6c <- 
    mod %>%          
    param (pars) %>%
    Req   (Plasma)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral6c, tgrid = tsamp6)   
  outdf6c <- cbind.data.frame(Time   = out6c$time/24,  Plasma = out6c$Plasma)
  
  out6c_1 <- mod %>%data_set(ex.oral6c_1) %>%
    idata_set(idata6c) %>% 
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%
    mrgsim(obsonly=TRUE, tgrid = tsamp6)
  outdf6c_1 = cbind.data.frame (ID     = out6c_1$ID, Time   = out6c_1$time/24, Plasma = out6c_1$Plasma)
  
  return (list(outdf6a, outdf6a_1, outdf6b, outdf6b_1, outdf6c, outdf6c_1))
}
N = 1000

##a
M6a_m       <- pred.eva (pars.eva,N = N)[[1]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = Plasma)  %>% mutate(Sex = "Male", Study= "M6", Dose = 1, Tissue = "Plasma", Matrix = c("Pre.Plasma"))
M6a_m1      <- M6a_m  %>% filter(Time %in% c(1, 5, 14, 28))

M6a         <- pred.eva (pars.eva,N = N)[[2]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = Plasma)
M6a_summary <- M6a %>% 
  filter(Time %in% c(1, 5, 14, 28)) %>%
  group_by(Time) %>%
  summarize(
    Mean = mean(Conc, na.rm = TRUE),
    SD = sd(Conc, na.rm = TRUE),
    Median = median(Conc, na.rm = TRUE),
    P25 = quantile(Conc, probs = 0.25, na.rm = TRUE),
    P75 = quantile(Conc, probs = 0.75, na.rm = TRUE),
    ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
    ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)
  ) %>% mutate( Sex = "Male", Study= "M6", Dose = 1, Tissue = "Plasma", Matrix = c("Pre.Plasma"))

##b
M6b_m       <- pred.eva (pars.eva,N = N)[[3]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = Plasma)  %>% mutate(Sex = "Male", Study= "M6", Dose = 10, Tissue = "Plasma", Matrix = c("Pre.Plasma"))
M6b_m1      <- M6b_m  %>% filter(Time %in% c(1, 5, 14, 28))

M6b         <- pred.eva (pars.eva,N = N)[[4]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = Plasma)
M6b_summary <- M6b %>% 
  filter(Time %in% c(1, 5, 14, 28)) %>%
  group_by(Time) %>%
  summarize(
    Mean = mean(Conc, na.rm = TRUE),
    SD = sd(Conc, na.rm = TRUE),
    Median = median(Conc, na.rm = TRUE),
    P25 = quantile(Conc, probs = 0.25, na.rm = TRUE),
    P75 = quantile(Conc, probs = 0.75, na.rm = TRUE),
    ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
    ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)
  ) %>% mutate( Sex = "Male", Study= "M6", Dose = 10, Tissue = "Plasma", Matrix = c("Pre.Plasma"))

##c
M6c_m       <- pred.eva (pars.eva,N = N)[[5]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = Plasma)  %>% mutate(Sex = "Male", Study= "M6", Dose = 100, Tissue = "Plasma", Matrix = c("Pre.Plasma"))
M6c_m1      <- M6c_m  %>% filter(Time %in% c(1, 5, 14, 28))

M6c         <- pred.eva (pars.eva,N = N)[[6]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = Plasma)
M6c_summary <- M6c %>% 
  filter(Time %in% c(1, 5, 14, 28)) %>%
  group_by(Time) %>%
  summarize(
    Mean = mean(Conc, na.rm = TRUE),
    SD = sd(Conc, na.rm = TRUE),
    Median = median(Conc, na.rm = TRUE),
    P25 = quantile(Conc, probs = 0.25, na.rm = TRUE),
    P75 = quantile(Conc, probs = 0.75, na.rm = TRUE),
    ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
    ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)
  ) %>% mutate( Sex = "Male", Study= "M6", Dose = 100, Tissue = "Plasma", Matrix = c("Pre.Plasma"))

M6a_curve <- M6a %>% 
  group_by(Time) %>%
  summarize(
    Mean = mean(Conc, na.rm = TRUE),
    SD = sd(Conc, na.rm = TRUE),
    Median = median(Conc, na.rm = TRUE),
    P25 = quantile(Conc, probs = 0.25, na.rm = TRUE),
    P75 = quantile(Conc, probs = 0.75, na.rm = TRUE),
    ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
    ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)
  ) %>% ungroup()%>% mutate( Sex = "Male", Dose = 1)

M6b_curve <- M6b %>% 
  group_by(Time) %>%
  summarize(
    Mean = mean(Conc, na.rm = TRUE),
    SD = sd(Conc, na.rm = TRUE),
    Median = median(Conc, na.rm = TRUE),
    P25 = quantile(Conc, probs = 0.25, na.rm = TRUE),
    P75 = quantile(Conc, probs = 0.75, na.rm = TRUE),
    ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
    ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)
  ) %>% ungroup()%>% mutate( Sex = "Male", Dose = 10)

M6c_curve <- M6c %>% 
  group_by(Time) %>%
  summarize(
    Mean = mean(Conc, na.rm = TRUE),
    SD = sd(Conc, na.rm = TRUE),
    Median = median(Conc, na.rm = TRUE),
    P25 = quantile(Conc, probs = 0.25, na.rm = TRUE),
    P75 = quantile(Conc, probs = 0.75, na.rm = TRUE),
    ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
    ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)
  ) %>% ungroup()%>% mutate( Sex = "Male", Dose = 100)

M6_curve   <-rbind(M6a_curve, M6b_curve, M6c_curve)
write.csv(M6_curve, file = 'M6_curve.csv')
M6_m1      <- rbind(M6a_m1, M6b_m1, M6c_m1)
M6_summary <- rbind(M6a_summary, M6b_summary, M6c_summary)%>% 
  relocate(Time, .after = 8)

Male_eva <- rbind(M1_summary, M2_summary, M3_summary,M4_summary,M5_summary,M6_summary)
write.csv(Male_eva, file = 'Male_eva.csv')
Pla_M  <- Male_eva %>% filter(Tissue == "Plasma")
CL_M   <- Male_eva %>% filter(Tissue == "Liver")
write.csv(Pla_M, file = 'Pla_M.csv')
write.csv(CL_M , file = 'CL_M.csv')
Male_eva_fix <- rbind(M1_m1, M2_m1,M3_m1,M4_m1,M5_m1,M6_m1)
write.csv(Male_eva_fix, file = 'Male_eva_fix.csv')

##Female ############################################################################################################################
setwd("D:/zs/PBPK/1-方案/GenX/小鼠/code/Mice/Female")
evaluation <- read.csv(file = "Evalution-F.csv")

OBS.A_CL  <- evaluation %>% filter(Study == 1 & Sample == "Liver" ) %>% select(Time = "Time", CL     = "Conc", SD = "SD", Study = "Study", Dose = "Dose") # TDoses = 28
OBS.B_Pla <- evaluation %>% filter(Study == 2 & Sample == "Plasma") %>% select(Time = "Time", Plasma = "Conc", SD = "SD", Study = "Study", Dose = "Dose") # TDoses = 84
OBS.C_Pla <- evaluation %>% filter(Study == 3 & Sample == "serum" ) %>% select(Time = "Time", Plasma = "Conc", SD = "SD", Study = "Study", Dose = "Dose") # TDoses = 28
OBS_Pla <- rbind.data.frame (OBS.B_Pla, OBS.C_Pla)
OBS_Pla <- OBS_Pla %>% mutate(Matrix = c("Obs.Plasma"))
OBS_CL  <- rbind.data.frame (OBS.A_CL)
OBS_CL  <- OBS_CL  %>% mutate(Matrix = c("Obs.Liver"))

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

##F1
pred.eva <- function(pars, N) { 
  pars <- exp(pars)                   
  BW          = 0.02                  
  tinterval   = 24                     
  TDOSE1    = 28                     
  DOSE1     = 100                    
  DOSEoral1 = DOSE1*BW  
  idata1 <- tibble(ID = 1:N) %>%  mutate( BW        = rnormTrunc  (N, min = 0.0082, max = 0.0318, mean = 0.02 ,      sd = 0.006 ),   #CV=0.3
                                          VKC       = rnormTrunc  (N, min = 0.007 , max = 0.0270, mean = 0.017,      sd = 0.0051),   #CV=0.3
                                          GFRC      = rnormTrunc  (N, min = 16.913, max = 65.167, mean = 41.04,      sd = 12.31 ) ,  #CV=0.3
                                          Free                           = 0.0704,                ### fitting parameters
                                          PL                             = 1.0649,                 ### fitting parameters
                                          PLu                            = 0.331,                 ### fitting parameters
                                          DOSEoral1 = DOSE1*BW)   
  ex.oral1  <- ev(ID = 1, amt = DOSEoral1, ii = tinterval, addl = TDOSE1 - 1, cmt = "AST", replicate = FALSE)
  ex.oral1_1  <- ev(ID = 1:N, amt = idata1$DOSEoral1, ii = tinterval, addl = TDOSE1 - 1, cmt = "AST", replicate = FALSE)
  tsamp1    = tgrid(0, tinterval*(TDOSE1 - 1) + 24*7, 1)  
  
  out1 <- 
    mod %>%          
    param (pars) %>%
    Req   (Liver)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral1, tgrid = tsamp1) 
  outdf1 <- cbind.data.frame(Time = out1$time/24, CL = out1$Liver) 
  
  out1_1 <- mod %>%data_set(ex.oral1_1) %>%
    idata_set(idata1) %>% 
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%
    mrgsim(obsonly=TRUE, tgrid = tsamp1)
  outdf1_1 = cbind.data.frame (ID     = out1_1$ID, Time   = out1_1$time/24, CL = out1_1$Liver)
  return (list(outdf1, outdf1_1))
}
N = 1000

F1_m   <- pred.eva (pars.eva,N = N)[[1]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CL) %>% mutate( Sex = "Female", Study= "F1", Dose = 100,Tissue = "Liver", Matrix = c("Pre.Liver"))
F1_m1  <- F1_m  %>% filter(Time == 29)    
F1     <- pred.eva (pars.eva,N = N)[[2]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CL) 
F1_1   <- F1    %>% filter(Time == 29)
F1_summary <- F1_1 %>% summarize(Mean = mean(Conc, na.rm = TRUE),
                                 SD = sd(Conc, na.rm = TRUE),
                                 Median = median(Conc, na.rm = TRUE),
                                 P25 = quantile(Conc, probs = 0.25, na.rm = TRUE),
                                 P75 = quantile(Conc, probs = 0.75, na.rm = TRUE),
                                 ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
                                 ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)) 
F1_summary <- F1_summary %>% mutate( Time = 29, Sex = "Female", Study= "F1", Dose = 100, Tissue = "Liver", Matrix = c("Pre.Liver"))

##F2
pred.eva <- function(pars, N) { 
  pars <- exp(pars)                   
  BW          = 0.02                  
  tinterval   = 24                     
  TDOSE2      = 84                      ## Total dosing/Dose times
  DOSE2       = 0.32                   ## Input oral dose
  DOSEoral2   = DOSE2*BW  
  idata2 <- tibble(ID = 1:N) %>%  mutate( BW        = rnormTrunc  (N, min = 0.0082, max = 0.0318, mean = 0.02 ,      sd = 0.006 ),   #CV=0.3
                                          VKC       = rnormTrunc  (N, min = 0.007 , max = 0.0270, mean = 0.017,      sd = 0.0051),   #CV=0.3
                                          GFRC      = rnormTrunc  (N, min = 16.913, max = 65.167, mean = 41.04,      sd = 12.31 ) ,  #CV=0.3
                                          Free                           = 0.0704,                ### fitting parameters
                                          PL                             = 1.0649,                 ### fitting parameters
                                          PLu                            = 0.331,                 ### fitting parameters
                                          DOSEoral2   = DOSE2*BW)   
  ex.oral2    <- ev(ID = 1, amt = DOSEoral2, ii = tinterval, addl = TDOSE2 - 1, cmt = "AST", replicate = FALSE)
  ex.oral2_1  <- ev(ID = 1:N, amt = idata2$DOSEoral2, ii = tinterval, addl = TDOSE2 - 1, cmt = "AST", replicate = FALSE)
  tsamp2      = tgrid(0, tinterval*(TDOSE2 - 1) + 24*7, 1) 
  
  out2 <- 
    mod %>%          
    param (pars) %>%
    Req   (Plasma)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral2, tgrid = tsamp2)  
  outdf2 <- cbind.data.frame(Time = out2$time/24, Plasma = out2$Plasma) 
  
  out2_1 <- mod %>%data_set(ex.oral2_1) %>%
    idata_set(idata2) %>% 
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%
    mrgsim(obsonly=TRUE, tgrid = tsamp2)
  outdf2_1 = cbind.data.frame (ID     = out2_1$ID, Time   = out2_1$time/24, Plasma = out2_1$Plasma)
  return (list(outdf2, outdf2_1))
}
N = 1000

F2_m   <- pred.eva (pars.eva,N = N)[[1]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = Plasma) %>% mutate( Sex = "Female", Study= "F2", Dose = 0.32,Tissue = "Plasma", Matrix = c("Pre.Plasma"))
F2_m1  <- F2_m %>% filter(Time == 84)    
F2     <- pred.eva (pars.eva,N = N)[[2]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = Plasma) 
F2_1   <- F2 %>% filter(Time == 84)
F2_summary <- F2_1 %>% summarize(Mean = mean(Conc, na.rm = TRUE),
                                 SD = sd(Conc, na.rm = TRUE),
                                 Median = median(Conc, na.rm = TRUE),
                                 P25 = quantile(Conc, probs = 0.25, na.rm = TRUE),
                                 P75 = quantile(Conc, probs = 0.75, na.rm = TRUE),
                                 ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
                                 ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)) 
F2_summary <- F2_summary %>% mutate( Time = 84, Sex = "Female", Study= "F2", Dose = 0.32, Tissue = "Plasma", Matrix = c("Pre.Plasma"))

##F3
pred.eva <- function(pars, N) { 
  pars <- exp(pars)                   
  BW          = 0.02                  
  tinterval   = 24                     
  TDOSE3      = 28
  tsamp3      = tgrid(0, tinterval*(TDOSE3 - 1) + 24*7, 1)
  
  #a
  DOSE3a      = 1                 
  DOSEoral3a  = DOSE3a*BW    
  idata3a <- tibble(ID = 1:N) %>%  mutate( BW        = rnormTrunc  (N, min = 0.0082, max = 0.0318, mean = 0.02 ,      sd = 0.006 ),   #CV=0.3
                                           VKC       = rnormTrunc  (N, min = 0.007 , max = 0.0270, mean = 0.017,      sd = 0.0051),   #CV=0.3
                                           GFRC      = rnormTrunc  (N, min = 16.913, max = 65.167, mean = 41.04,      sd = 12.31 ) ,  #CV=0.3
                                           Free                           = 0.0704,                ### fitting parameters
                                           PL                             = 1.0649,                 ### fitting parameters
                                           PLu                            = 0.331,                 ### fitting parameters
                                           DOSEoral3a  = DOSE3a*BW)  
  ex.oral3a   <- ev(ID = 1, amt = DOSEoral3a, ii = tinterval, addl = TDOSE3 - 1, cmt = "AST", replicate = FALSE)
  ex.oral3a_1  <- ev(ID = 1:N, amt = idata3a$DOSEoral3a, ii = tinterval, addl = TDOSE3 - 1, cmt = "AST", replicate = FALSE)
  out3a <- 
    mod %>%          
    param (pars) %>%
    Req   (Plasma)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral3a, tgrid = tsamp3) 
  outdf3a <- cbind.data.frame(Time   = out3a$time/24, Plasma = out3a$Plasma) 
  
  out3a_1 <- mod %>%data_set(ex.oral3a_1) %>%
    idata_set(idata3a) %>% 
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%
    mrgsim(obsonly=TRUE, tgrid = tsamp3)
  outdf3a_1 = cbind.data.frame (ID     = out3a_1$ID, Time   = out3a_1$time/24, Plasma = out3a_1$Plasma)
  
  #b
  DOSE3b      = 10                 
  DOSEoral3b  = DOSE3b*BW           
  idata3b <- tibble(ID = 1:N) %>%  mutate( BW        = rnormTrunc  (N, min = 0.0082, max = 0.0318, mean = 0.02 ,      sd = 0.006 ),   #CV=0.3
                                           VKC       = rnormTrunc  (N, min = 0.007 , max = 0.0270, mean = 0.017,      sd = 0.0051),   #CV=0.3
                                           GFRC      = rnormTrunc  (N, min = 16.913, max = 65.167, mean = 41.04,      sd = 12.31 ) ,  #CV=0.3
                                           Free                           = 0.0704,                ### fitting parameters
                                           PL                             = 1.0649,                 ### fitting parameters
                                           PLu                            = 0.331,                 ### fitting parameters
                                           DOSEoral3b  = DOSE3b*BW)  
  ex.oral3b   <- ev(ID = 1, amt = DOSEoral3b, ii = tinterval, addl = TDOSE3 - 1, cmt = "AST", replicate = FALSE)
  ex.oral3b_1  <- ev(ID = 1:N, amt = idata3b$DOSEoral3b, ii = tinterval, addl = TDOSE3 - 1, cmt = "AST", replicate = FALSE)
  out3b <- 
    mod %>%          
    param (pars) %>%
    Req   (Plasma)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral3b, tgrid = tsamp3) 
  outdf3b <- cbind.data.frame(Time   = out3b$time/24, Plasma = out3b$Plasma) 
  
  out3b_1 <- mod %>%data_set(ex.oral3b_1) %>%
    idata_set(idata3b) %>% 
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%
    mrgsim(obsonly=TRUE, tgrid = tsamp3)
  outdf3b_1 = cbind.data.frame (ID     = out3b_1$ID, Time   = out3b_1$time/24, Plasma = out3b_1$Plasma)
  
  #c
  DOSE3c      = 100                 
  DOSEoral3c  = DOSE3c*BW         
  idata3c <- tibble(ID = 1:N) %>%  mutate( BW        = rnormTrunc  (N, min = 0.0082, max = 0.0318, mean = 0.02 ,      sd = 0.006 ),   #CV=0.3
                                           VKC       = rnormTrunc  (N, min = 0.007 , max = 0.0270, mean = 0.017,      sd = 0.0051),   #CV=0.3
                                           GFRC      = rnormTrunc  (N, min = 16.913, max = 65.167, mean = 41.04,      sd = 12.31 ) ,  #CV=0.3
                                           Free                           = 0.0704,                ### fitting parameters
                                           PL                             = 1.0649,                 ### fitting parameters
                                           PLu                            = 0.331,                 ### fitting parameters
                                           DOSEoral3c  = DOSE3c*BW )  
  ex.oral3c   <- ev(ID = 1, amt = DOSEoral3c, ii = tinterval, addl = TDOSE3 - 1, cmt = "AST", replicate = FALSE)
  ex.oral3c_1  <- ev(ID = 1:N, amt = idata3c$DOSEoral3c, ii = tinterval, addl = TDOSE3 - 1, cmt = "AST", replicate = FALSE)
  out3c <- 
    mod %>%          
    param (pars) %>%
    Req   (Plasma)%>%
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%  
    mrgsim_d (data = ex.oral3c, tgrid = tsamp3)   
  outdf3c <- cbind.data.frame(Time   = out3c$time/24,  Plasma = out3c$Plasma) 
  
  out3c_1 <- mod %>%data_set(ex.oral3c_1) %>%
    idata_set(idata3c) %>% 
    update(atol = 1E-5, rtol= 1e-5, maxsteps = 2000) %>%
    mrgsim(obsonly=TRUE, tgrid = tsamp3)
  outdf3c_1 = cbind.data.frame (ID     = out3c_1$ID, Time   = out3c_1$time/24, Plasma = out3c_1$Plasma)
  
  return (list(outdf3a, outdf3a_1, outdf3b, outdf3b_1, outdf3c, outdf3c_1))
}
N = 1000

##a
F3a_m       <- pred.eva (pars.eva,N = N)[[1]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = Plasma)  %>% mutate(Sex = "Female", Study= "F3", Dose = 1, Tissue = "Plasma", Matrix = c("Pre.Plasma"))
F3a_m1      <- F3a_m  %>% filter(Time %in% c(1, 5, 14, 28))

F3a         <- pred.eva (pars.eva,N = N)[[2]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = Plasma)
F3a_summary <- F3a %>% 
  filter(Time %in% c(1, 5, 14, 28)) %>%
  group_by(Time) %>%
  summarize(
    Mean = mean(Conc, na.rm = TRUE),
    SD = sd(Conc, na.rm = TRUE),
    Median = median(Conc, na.rm = TRUE),
    P25 = quantile(Conc, probs = 0.25, na.rm = TRUE),
    P75 = quantile(Conc, probs = 0.75, na.rm = TRUE),
    ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
    ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)
  ) %>% mutate( Sex = "Female", Study= "F3", Dose = 1, Tissue = "Plasma", Matrix = c("Pre.Plasma"))

##b
F3b_m       <- pred.eva (pars.eva,N = N)[[3]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = Plasma)  %>% mutate(Sex = "Female", Study= "F3", Dose = 10, Tissue = "Plasma", Matrix = c("Pre.Plasma"))
F3b_m1      <- F3b_m  %>% filter(Time %in% c(1, 5, 14, 28))

F3b         <- pred.eva (pars.eva,N = N)[[4]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = Plasma)
F3b_summary <- F3b %>% 
  filter(Time %in% c(1, 5, 14, 28)) %>%
  group_by(Time) %>%
  summarize(
    Mean = mean(Conc, na.rm = TRUE),
    SD = sd(Conc, na.rm = TRUE),
    Median = median(Conc, na.rm = TRUE),
    P25 = quantile(Conc, probs = 0.25, na.rm = TRUE),
    P75 = quantile(Conc, probs = 0.75, na.rm = TRUE),
    ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
    ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)
  ) %>% mutate( Sex = "Female", Study= "F3", Dose = 10, Tissue = "Plasma", Matrix = c("Pre.Plasma"))

##c
F3c_m       <- pred.eva (pars.eva,N = N)[[5]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = Plasma)  %>% mutate(Sex = "Female", Study= "F3", Dose = 100, Tissue = "Plasma", Matrix = c("Pre.Plasma"))
F3c_m1      <- F3c_m  %>% filter(Time %in% c(1, 5, 14, 28))

F3c         <- pred.eva (pars.eva,N = N)[[6]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = Plasma)
F3c_summary <- F3c %>% 
  filter(Time %in% c(1, 5, 14, 28)) %>%
  group_by(Time) %>%
  summarize(
    Mean = mean(Conc, na.rm = TRUE),
    SD = sd(Conc, na.rm = TRUE),
    Median = median(Conc, na.rm = TRUE),
    P25 = quantile(Conc, probs = 0.25, na.rm = TRUE),
    P75 = quantile(Conc, probs = 0.75, na.rm = TRUE),
    ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
    ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)
  ) %>% mutate( Sex = "Female", Study= "F3", Dose = 100, Tissue = "Plasma", Matrix = c("Pre.Plasma"))

F3a_curve <- F3a %>% 
  group_by(Time) %>%
  summarize(
    Mean = mean(Conc, na.rm = TRUE),
    SD = sd(Conc, na.rm = TRUE),
    Median = median(Conc, na.rm = TRUE),
    P25 = quantile(Conc, probs = 0.25, na.rm = TRUE),
    P75 = quantile(Conc, probs = 0.75, na.rm = TRUE),
    ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
    ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)
  ) %>% ungroup()%>% mutate( Sex = "Female", Dose = 1)

F3b_curve <- F3b %>% 
  group_by(Time) %>%
  summarize(
    Mean = mean(Conc, na.rm = TRUE),
    SD = sd(Conc, na.rm = TRUE),
    Median = median(Conc, na.rm = TRUE),
    P25 = quantile(Conc, probs = 0.25, na.rm = TRUE),
    P75 = quantile(Conc, probs = 0.75, na.rm = TRUE),
    ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
    ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)
  ) %>% ungroup()%>% mutate( Sex = "Female", Dose = 10)

F3c_curve <- F3c %>% 
  group_by(Time) %>%
  summarize(
    Mean = mean(Conc, na.rm = TRUE),
    SD = sd(Conc, na.rm = TRUE),
    Median = median(Conc, na.rm = TRUE),
    P25 = quantile(Conc, probs = 0.25, na.rm = TRUE),
    P75 = quantile(Conc, probs = 0.75, na.rm = TRUE),
    ci_lower_est = quantile(Conc, probs = 0.025, names = FALSE, na.rm = TRUE),
    ci_upper_est = quantile(Conc, probs = 0.975, names = FALSE, na.rm = TRUE)
  ) %>% ungroup()%>% mutate( Sex = "Female", Dose = 100)

F3_curve   <- rbind(F3a_curve, F3b_curve, F3c_curve)
write.csv(F3_curve, file = 'F3_curve.csv')
F3_m1      <- rbind(F3a_m1, F3b_m1, F3c_m1)
F3_summary <- rbind(F3a_summary, F3b_summary, F3c_summary) %>% 
  relocate(Time, .after = 8)

Female_eva <- rbind(F1_summary, F2_summary, F3_summary)
write.csv(Female_eva, file = 'Female_eva.csv')
Pla_F  <- Female_eva %>% filter(Tissue == "Plasma")
CL_F   <- Female_eva %>% filter(Tissue == "Liver")
write.csv(Pla_F, file = 'Pla_F.csv')
write.csv(CL_F , file = 'CL_F.csv')
Female_eva_fix <- rbind(F1_m1, F2_m1,F3_m1)
write.csv(Female_eva_fix, file = 'Female_eva_fix.csv')

## Plot ###########################################################################################################################
setwd("D:/zs/PBPK/2025GenX/Modfit/Mice/Figure 3")
Figure_3_data<- read.csv(file = "Figure 3_data.csv")
CL_M    <- read.csv(file = "CL_M.csv")
Pla_M    <- read.csv(file = "Pla_M.csv")
Male_eva <- read.csv(file = "Male_eva.csv")
M6_curve   <- read.csv(file = "M6_curve.csv")%>% select(-1)
Male_eva_fix <- read.csv(file = "Male_eva_fix.csv")

CL_F       <- read.csv(file = "CL_F.csv")
Pla_F      <- read.csv(file = "Pla_F.csv")
Female_eva <- read.csv(file = "Female_eva.csv")
F3_curve   <- read.csv(file = "F3_curve.csv") %>% select(-1)
Female_eva_fix <- read.csv(file = "Female_eva_fix.csv")

combined_Pla <- rbind(Pla_M , Pla_F )%>% select(-1)
write.csv(combined_Pla, file = 'combined_Pla.csv')
combined_CL <- rbind(CL_M , CL_F )%>% select(-1)
write.csv(combined_CL , file = 'combined_CL.csv')

Figure_3_A<- read.csv(file = "Figure 3_A.csv")
combined_curve <- rbind(M6_curve , F3_curve ) %>% select(Time, Mean, SD, Sex, Dose)

## Sex ##################################################################################################
A1_OBS <- Figure_3_A %>% filter(Sex=="Male")
A1<- combined_curve %>% filter(Sex=="Male")
plot.A1 <- 
  ggplot(A1, aes(x = Time, y = Mean, color = factor(Dose), fill = factor(Dose))) + 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.3, colour = NA) +
  geom_line(size = 0.8) +
  geom_point(data = A1_OBS, aes(x = Time, y = Mean, color = factor(Dose)), size = 2.5) +
  geom_errorbar(data = A1_OBS, aes(x = Time, ymin = Mean - SD, ymax = Mean + SD, color = factor(Dose)), 
                width = 0.3, size = 1) +
  ylab("GenX Concentration in Plasma (μg/ml)") +
  xlab("Time (days)") +
  scale_color_manual(name = "Dose (mg/kg/day)", values = c("1" = "steelblue", "10" = "#1E803D", "100" = "#715ea9")) +
  scale_fill_manual(name = "Dose (mg/kg/day)", values = c("1" = "lightblue", "10" = "#B0D9A5", "100" = "#CECCE5")) +
  scale_y_continuous(limits = c(-100, 900), expand = c(0, 0)) +
  theme_bw() +
  theme(
    panel.border = element_rect(color = "grey10", fill = NA, linewidth = 1),
    strip.text = element_text(size = rel(2), colour = "grey10"),
    axis.text = element_text(size = rel(2), colour = "black"),
    axis.line = element_line(color = "black", size = 0.3),
    axis.title = element_text(size = 25, colour = "black", face = "bold"),
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5, vjust = -5, margin = margin(b = 10)),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.position = c(0.92, 0.92)  
  ) +
  labs(title = "(A1) Male")
print(plot.A1)

## Dose ######################################################################################################################
A1_OBS <- Figure_3_A %>% filter(Dose==1)
A1<- combined_curve %>% filter(Dose==1)
plot.A1 <- 
  ggplot(A1, aes(x = Time, y = Mean, color = Sex, fill = Sex)) + 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.3, colour = NA) +
  geom_line(size = 0.8) +
  geom_point(data = A1_OBS, aes(x = Time, y = Mean, color = Sex), size = 3.5) +
  geom_errorbar(data = A1_OBS, aes(x = Time, ymin= Mean-SD, ymax = Mean+SD, color = Sex), width = 0.3, size = 1) +
  ylab("GenX Concentration in Plasma (μg/ml)") +
  xlab("Time (days)") +
  scale_color_manual(values = c("Male" = "steelblue", "Female" = "#E9687A")) +
  scale_fill_manual(values = c("Male" = "lightblue", "Female" = "#F6B3AC")) +
  scale_y_continuous(limits = c(-0.4, 9), expand = c(0, 0)) +
  theme_bw(base_family = "Times New Roman") +
  theme(
    panel.border = element_rect(color = "grey10", fill = NA, linewidth = 1),
    strip.text = element_text(size = rel(2), colour = "grey10"),
    axis.text = element_text(size = rel(2), colour = "black"),
    axis.line = element_line(color = "black", size = 0.3),
    axis.title = element_text(size = 25, colour = "black", face = "bold"),
    plot.title = element_text(size = 25, face = "bold", hjust = 0, vjust = 0, margin = margin(b = 10)),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    legend.text  = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.background = element_rect(fill = NA, colour = NA),
    legend.position = c(0.89, 0.9)  
  ) +
  labs(title = "(A1)")
print(plot.A1)

A2<- combined_curve %>% filter(Dose==10)
A2_OBS <- Figure_3_A %>% filter(Dose==10)
plot.A2 <- 
  ggplot(A2, aes(x = Time, y = Mean, color = Sex, fill = Sex)) + 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.3, colour = NA) +
  geom_line(size = 0.8) +
  geom_point(data = A2_OBS, aes(x = Time, y = Mean, color = Sex), size = 3.5) +
  geom_errorbar(data = A2_OBS, aes(x = Time, ymin= Mean-SD, ymax = Mean+SD, color = Sex), width = 0.25, size = 1) +
  ylab("GenX Concentration in Plasma (μg/ml)") +
  xlab("Time (days)") +
  scale_color_manual(values = c("Male" = "steelblue", "Female" = "#E9687A")) +
  scale_fill_manual(values = c("Male" = "lightblue", "Female" = "#F6B3AC")) +
  scale_y_continuous(limits = c(-4, 90), expand = c(0, 0)) +
  theme_bw(base_family = "Times New Roman") +
  theme(
    panel.border = element_rect(color = "grey10", fill = NA, linewidth = 1),
    strip.text = element_text(size = rel(2), colour = "grey10"),
    axis.text = element_text(size = rel(2), colour = "black"),
    axis.line = element_line(color = "black", size = 0.3),
    axis.title = element_text(size = 25, colour = "black", face = "bold"),
    plot.title = element_text(size = 25, face = "bold", hjust = 0, vjust = 0, margin = margin(b = 10)),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    legend.text  = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.background = element_rect(fill = NA, colour = NA),
    legend.position = c(0.89, 0.9) 
  ) +
  labs(title = "(A2)")
print(plot.A2)

A3<- combined_curve %>% filter(Dose==100)
A3_OBS <- Figure_3_A %>% filter(Dose==100)
plot.A3 <- 
  ggplot(A3, aes(x = Time, y = Mean, color = Sex, fill = Sex)) + 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.3, colour = NA) +
  geom_line(size = 0.8) +
  geom_point(data = A3_OBS, aes(x = Time, y = Mean, color = Sex), size = 3.5) +
  geom_errorbar(data = A3_OBS, aes(x = Time, ymin= Mean-SD, ymax = Mean+SD, color = Sex), width = 0.25, size = 1) +
  ylab("GenX Concentration in Plasma (μg/ml)") +
  xlab("Time (days)") +
  scale_color_manual(values = c("Male" = "steelblue", "Female" = "#E9687A")) +
  scale_fill_manual(values = c("Male" = "lightblue", "Female" = "#F6B3AC")) +
  scale_y_continuous(limits = c(-35, 900), expand = c(0, 0)) +
  theme_bw(base_family = "Times New Roman") +
  theme(
    panel.border = element_rect(color = "grey10", fill = NA, linewidth = 1),
    strip.text = element_text(size = rel(2), colour = "grey10"),
    axis.text = element_text(size = rel(2), colour = "black"),
    axis.line = element_line(color = "black", size = 0.3),
    axis.title = element_text(size = 25, colour = "black", face = "bold"),
    plot.title = element_text(size = 25, face = "bold", hjust = 0, vjust = 0, margin = margin(b = 10)),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    legend.text  = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.background = element_rect(fill = NA, colour = NA),
    legend.position = c(0.89, 0.9) 
  ) +
  labs(title = "(A3)")
print(plot.A3)

plot.A2 <- plot.A2 + theme(axis.title.y = element_blank())
plot.A3 <- plot.A3 + theme(axis.title.y = element_blank())
ggsave("Figure 3A.tiff",scale = 1,
       plot = grid.arrange(plot.A1, plot.A2, plot.A3, ncol = 3),
       path = "D:/zs/PBPK/2025GenX/Artwork",
       width = 60, height = 20, units = "cm",dpi=320)
