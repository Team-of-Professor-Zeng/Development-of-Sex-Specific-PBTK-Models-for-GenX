MicePBPK_F.code <- '
$PARAM @annotated

//Physiological parameters		

BW                  :   0.02  : kg,                  Bodyweight (Brown 1997)
QCC                 :   16.5  : L/h/kg^0.75,         Cardiac output (Brown 1997)
QLC                 :   0.161 : Unitless,            Fraction blood flow to liver (Brown 1997)
QLuC                :   0.005 : Unitless,            Fraction blood flow to lung (Brown 1997)
QKC                 :   0.091	: Unitless,            Fraction blood flow to kidney (Brown 1997)
Htc                 :   0.48  : Unitless,            Hematocrit for Mice (Hejtmancik 2002)
VPlasC              :   0.049 : L/kg BW,             Fractional plasma (Brown 1997)
VLC                 :   0.055 : Unitless,            Fractional liver tissue (Brown 1997)
VLuC                :   0.007 : Unitless,            Fractional lung tissue (Brown 1997)
VKC                 :   0.017 : Unitless,            Fractional kidney tissue (Brown 1997)
VfilC               :  0.0017 : L/kg BW,             Fraction vol. of filtrate; 10% of Kidney volume; ((Brown 1997)
FVBK                :   0.160 : Unitless,            Blood volume fraction of kidney (Brown, 1997) 

//Chemical-specific parameters 		

Free                :   0.045 : Unitless,            Free fraction; (Loccisano 2012)  
PL                  :   1.339 : Unitless,            Liver/ plasma PC; (Wen 2022) 
PLu                 :   0.431 : Unitless,            Lung/ plasma PC; (Wen 2022) 
PK                  :   0.854 : Unitless,            Kidney/ plasma PC; (Wen 2022) 
PRest               :   0.595 : Unitless,            Restofbody/ plasma PC; (Wen 2022)

GFRC                :  41.04  : L/hr/kg kiney,       Glomerular filtration rate (female) (Corley, 2005)
GEC                 :   0.54  : 1/(h*BW^0.25),       Gastric emptying time  (Yang et al., 2015)

Kabsc               :   2.12   : 1/(h*BW^-0.25),      Rate of absorption of GenX  from small intestine to liver; (initial value assumed the same as PFOA (2.12) from Worley and Fisher, 2015 and then re-fitting)    
K0C                 :   1.0    : 1/(h*BW^-0.25),      Rate of uptake from the stomach into the liver (initial value assumed the same as PFOA (1) from Worley and Fisher, 2015 and then re-fitting)
KunabsC             :  0.0265  : 1/(h*BW^-0.25),      Rate of unabsorbed dose to appear in feces; ) 
KurineC             :  0.122   : 1/(h*BW^-0.25),      Rate of urine elimination from urine storage; (Gannon 2015)

$MAIN
double QC  = QCC*pow(BW, 0.75)*(1-Htc);              // L/h, Cardiac output (adjusted for plasma)
double QK  = QKC*QC;                                 // L/h, Plasma flow to kidney
double QL  = QLC*QC;                                 // L/h, Plasma flow to liver
double QLu = QLuC*QC;                                // L/h, Plasma flow to lung
double QRest = QC-QK-QL-QLu;                         // L/h, Plasma flow to the rest of body

double VL = VLC*BW;                                  // L,   Volume of liver 
double VLu= VLuC*BW;                                 // L,   Volume of lung
double VPlas = VPlasC*BW;                            // L,   Volume of plasma
double VK = VKC*BW;                                  // L,   Volume of kidney 
double Vfil = VfilC*BW;                              // L,   Volume of filtrate  
double VKb = VK*FVBK;                                // L,   Volume of blood in the kidney; fraction blood volume of kidney (0.16) from Brown, 1997
double VRest = (0.93*BW)-VL-VLu-VK-VPlas-Vfil;            // L,   Rest of body; volume of remaining tissue (L); Revised original equation (VR = (0.93*BW) - VPlas - VPTC - Vfil - VL) from Worley and Fisher, 2015  
double MK = VK*1000;                                 // g,   kidney weight in gram

double Kurine = KurineC*pow(BW,(-0.25));             // 1/h, Urinary elimination; 
double GFR = GFRC*(MK/1000);                         // L/h, Glomerular filtration rate, scaled to mass of kidney 

//GI tract parameters
double Kabs = Kabsc*pow(BW,(-0.25));                 // 1/h, Rate of absorption of GenX  from small intestine to liver
double Kunabs = KunabsC*pow(BW,(-0.25));             // 1/h, Rate of unabsorbed dose to appear in feces
double GE = GEC*pow(BW,(-0.25));                     // 1/h, Gastric emptying time 
double K0 = K0C*pow(BW,(-0.25));                     // 1/h, Rate of uptake from the stomach into the liver

$CMT APlas_free AUCCA_free ARest AUCCRest AST AabsST ASI AabsSI Afeces Aurine AL AUCCL ALu AUCCLu AKb AUCCK AFil AUCCfil 

$INIT @annotated
ADOSE:0.2: mg, Amount of input dose; assumed a virtual compartment for validating the model mass balance 

$ODE

// Concentrations in plasma
double CA_free  = APlas_free/VPlas;                  // mg/L, Free GenX  concentration in the plasma
double CA       = CA_free/Free;                      // mg/L, Concentration of total GenX in the plasma

// Concentrations in liver
double CL = AL/VL;                                   // mg/L, Concentration of GenX  in the liver compartment
double CVL = CL/PL;                                  // mg/L, Concentration of GenX  in venous plasma leaving liver

// Concentrations in kidney
double CKb = AKb/VKb;                                // mg/L, Concetraitons of GenX  in venous plasma leaving kidney 
double CVK = CKb;                                    // mg/L, Concentration of GenX  in plasma leaving kidney
double CK  = CVK*PK;                                 // mg/L, Concetraitons of GenX  in Kidney compartment

// Concentrations in lung
double CLu = ALu/VLu;                                // mg/L, Concentration of GenX  in the lung compartment
double CVLu= CLu/PLu;                                // mg/L, Concentration of GenX  in venous plasma leaving lung

// Concentrations in Rest of body
double CRest = ARest/VRest;                          // mg/L, Concentration of GenX  in the rest of the body
double CVRest = ARest/(VRest*PRest);                 // mg/L, Concentration of GenX  in the venous palsma leaving the rest of body

// Kidney compartment plus 1 subcompartment (Filtrate: Fil)
// Concentration in kidney and fil

double Cfil = AFil/Vfil;                             // mg/L, Concentraitons of GenX  in Fil

// {GenX  distribution in each compartment}
// Free GenX  in plasma
double RPlas_free = (QRest*CVRest*Free)+(QK*CVK*Free)+(QL*CVL*Free)+(QLu*CVLu*Free)-(QC*CA*Free);  // mg/h,   Rate of change in the plasma 
dxdt_APlas_free = RPlas_free;                                                                      // mg,     Amount of free GenX  in the plasma
dxdt_AUCCA_free = CA_free;                                                                         // mg*h/L, Area under curve of free GenX  in plasma compartment

// Proximal Tubule Lumen/ Filtrate (Fil)
double Rfil = CA*GFR*Free - AFil*Kurine;                                                    // mg/h,   Rate of change in the Fil
dxdt_AFil = Rfil;                                                                           // mg,     Amount in the Fil
dxdt_AUCCfil = Cfil;                                                                        // mg*h/L, Area under curve of GenX  in the compartment of Fil

// Urine elimination
double Rurine = Kurine*AFil;                                                                // mg/h,   Rate of change in urine
dxdt_Aurine = Rurine;                                                                       // mg,     Amount in urine

// Kidney compartment
double RKb = QK*(CA-CVK)*Free - CA*GFR*Free ;                                               // mg/h,   Rate of change in Kidney compartment
dxdt_AKb = RKb;                                                                             // mg,     Amount in kidney compartment
dxdt_AUCCK= CK;                                                                             // mg*h/L, Area under curve of GenX  in the Kidney compartment

// GenX  in the compartment of rest of body, flow-limited model
double RRest = QRest*(CA-CVRest)*Free;                                                      // mg/h,   Rate of change in rest of body
dxdt_ARest = RRest;                                                                         // mg,     Amount in rest of body 
dxdt_AUCCRest = CRest;                                                                      // mg*h/L, Area under curve of GenX  in the compartment of rest of body

// GenX  in the compartment of Lung, flow-limited model
double RLu = QLu*(CA-CVLu)*Free;                                                            // mg/h,   Rate of change in Lung
dxdt_ALu   = RLu;                                                                           // mg,     Amount in Lung 
dxdt_AUCCLu= CLu;                                                                           // mg*h/L, Area under curve of GenX  in the compartment of Lung
    
// Gastrointestinal (GI) tract
// Stomach compartment
double RST = - K0*AST - GE*AST;                                                             // mg/h,   Rate of chagne in stomach compartment
dxdt_AST = RST;                                                                             // mg,     Amount in Stomach
double RabsST = K0*AST;                                                                     // mg/h,   Rate of absorption in the stomach
dxdt_AabsST = RabsST;                                                                       // mg,     Amount absorbed in the stomach

// Small intestine compartment
double RSI = GE*AST - Kabs*ASI - Kunabs*ASI;                                                // mg/h,   Rate of chagne in small intestine compartment
dxdt_ASI = RSI;                                                                             // mg,     Amount in small intestine
double RabsSI = Kabs*ASI;                                                                   // mg/h,   Rate of absorption in the small intestine
dxdt_AabsSI = RabsSI;                                                                       // mg,     Amount absorbed in the small intestine
double Total_oral_uptake = AabsSI + AabsST;                                                 // mg,     Total oral uptake in the GI

// Feces compartment
double Rfeces = Kunabs*ASI ;                                                                // mg/h,   Rate of change in feces compartment 
dxdt_Afeces = Rfeces;                                                                       // mg,     Amount of the feces compartment

// GenX  in liver compartment, flow-limited model
double RL = QL*(CA-CVL)*Free + Kabs*ASI + K0*AST;                                           // mg/h,   Rate of change in liver compartment
dxdt_AL = RL;                                                                               // mg,     Amount in liver compartment
dxdt_AUCCL = CL;                                                                            // mg*h/L, Area under curve of GenX  in liver compartment

// #+ Virtural compartment; input dose
dxdt_ADOSE         = 0; 

// {Mass balance equations}
double Qbal  = QC-QL-QK-QRest-QLu;
double Tmass = APlas_free + ARest + AKb + AFil + AL + AST + ASI + ALu;
double Loss  = Aurine + Afeces;
double Input = ADOSE;
double Bal   = Input- Tmass - Loss;

$TABLE  
capture Plasma = CA_free/Free;
capture Liver  = AL/VL;
capture Kidney = CK;
capture Lung   = ALu/VLu;
capture Rest   = ARest/VRest;
capture AUC_CA = AUCCA_free;
capture AUC_CL = AUCCL;
capture AUC_CK = AUCCK;
capture AUC_CLu= AUCCLu;
capture Balance= Bal;
capture QB     = Qbal;
'
