# load libraries
library(psych)
library(tidyverse)
library(dplyr)
library(reshape2)
library(viridis)
library(hrbrthemes)
library(ggpmisc)
library(cowplot)
library(ggpubr)
library(lavaan)
library(lcsm)
library(OpenMx)
library(rrr)
library(doBy)
library(Hmisc)
library(corrplot)
library(mice)
packages <- c("corrplot", "foreign", "gridExtra", "knitr", 
              "lavaan", "mice", "plyr", "tidyverse", "magrittr", "reshape2", 
              "naniar", "Hmisc", "devtools", 
              "psych")
invisible(lapply(packages, library, character.only = TRUE))


# load data
setwd("~/")
LCSM_table <- read_csv("W:/LCSM/LCSM_table.csv")

# subset adversity score and brain data
LCSM_DATA <- LCSM_table[,c("subcode", "all_site", "sex",  "BL_age", "FU3_age", "BL_TLFB_CigCB", "FU2_TLFB_CigCB", "FU3_TLFB_CigCB",
                           "BL_TLFB_alc", "FU2_TLFB_alc", "FU3_TLFB_alc","BL_TLFB_drug",  
                           "FU2_TLFB_drug","FU3_TLFB_drug","BL_SDQ_hyper", "FU2_SDQ_hyper", 
                           "FU3_SDQ_hyper","BL_SDQ_inattends","FU2_SDQ_inattends","FU3_SDQ_inattends",
                           "BL_ICV", "FU2_ICV", "BL_SSRT", "FU2_SSRT", "FU3_SSRT", "FU3_ICV",
                           "BL_SA_GT_PosNetStr", "FU2_SA_GT_PosNetStr", "FU3_SA_GT_PosNetStr",
                           "BL_SA_GT_NegNetStr", "FU2_SA_GT_NegNetStr", "FU3_SA_GT_NegNetStr",
                           "BL_SA_SS_PosNetStr", "FU2_SA_SS_PosNetStr", "FU3_SA_SS_PosNetStr",
                           "BL_SA_SS_NegNetStr", "FU2_SA_SS_NegNetStr", "FU3_SA_SS_NegNetStr",
                           "BL_IC_GT_PosNetStr", "FU2_IC_GT_PosNetStr", "FU3_IC_GT_PosNetStr",
                           "BL_IC_GT_NegNetStr", "FU2_IC_GT_NegNetStr", "FU3_IC_GT_NegNetStr",
                           "BL_IC_SS_PosNetStr", "FU2_IC_SS_PosNetStr", "FU3_IC_SS_PosNetStr",
                           "BL_IC_SS_NegNetStr", "FU2_IC_SS_NegNetStr", "FU3_IC_SS_NegNetStr",
                           "mode_centered_pds",  "BL_meanFD", "FU2_meanFD", "FU3_meanFD")]


LCSM_DATA %<>% replace_with_na_at(.vars = c("all_site","sex"),
                                  condition = ~.x == -100)

# scale all relevant numeric variables
LCSM_DATA[, -c(1,2,3)] <- scale(LCSM_DATA[, -c(1,2,3)])



################## Cigarette and cannabis use #################

# Latent change score model: ICV & CigCB

icv <- '

#####     The following lines specify the core assumptions of the LCS 
#####     and should not generally be modified

# Fixes regressions between CigCB measurements at each timepoint to 1
FU2_TLFB_CigCB ~ 1*BL_TLFB_CigCB
FU3_TLFB_CigCB ~ 1*FU2_TLFB_CigCB

FU2_ICV ~ 1*BL_ICV     
FU3_ICV ~ 1*FU2_ICV

# Fixes fator loadings on latent change scores to 1
dCIGCB1 =~ 1*FU2_TLFB_CigCB      
dCIGCB2 =~ 1*FU3_TLFB_CigCB

dNEU1 =~ 1*FU2_ICV      
dNEU2 =~ 1*FU3_ICV

# Constrains intercepts at t2 and t3 to 0 
FU2_TLFB_CigCB ~ 0*1            
FU3_TLFB_CigCB ~ 0*1

FU2_ICV ~ 0*1           
FU3_ICV ~ 0*1

 # This fixes the variance of T2 and T3 measures to 0 
FU2_TLFB_CigCB ~~ 0*FU2_TLFB_CigCB      
FU3_TLFB_CigCB ~~ 0*FU3_TLFB_CigCB

FU2_ICV ~~ 0*FU2_ICV       
FU3_ICV ~~ 0*FU3_ICV

###### The following parameters will be estimated in the model.

# Estimates the intercept of CigCB/COGNITION at T3 (What are mean scores at T1?)
BL_TLFB_CigCB ~ 1      
BL_ICV ~ 1 

# This estimates the variance of the change scores (How much variation is there in improvement?)
dCIGCB1 ~~  vardCigCB1*dCIGCB1        #  variances constrained to equality across waves
dCIGCB2 ~~  vardCigCB1*dCIGCB2
dNEU1 ~~  vardNEU1*dNEU1       
dNEU2 ~~  vardNEU1*dNEU2

# Estimates intercepts for latent change scores (equality constraints across timepoint with labels)
dCIGCB1 ~ CigCBint1 * 1   
dCIGCB2 ~ CigCBint2 * 1 # exception to equality constraint

dNEU1 ~ NEUint1 * 1   
dNEU2 ~ NEUint2 * 1  # exception to equality constraint

# This estimates the variance of T1 scores (How much variation is there at T1?)
BL_TLFB_CigCB ~~ BL_TLFB_CigCB      
BL_ICV ~~ BL_ICV 

# The following estimates the self-feedback parameters (Does change in cognition depend on cognition starting point?) and coupling parameters (Does change in cognition depend on block design starting point?)

# Feedback parameters allowed to vary across timepoints:

dCIGCB1 ~ selfCigCB1*BL_TLFB_CigCB + couplinga*BL_ICV  
dCIGCB2 ~ selfCigCB2*FU2_TLFB_CigCB + couplinga*FU2_ICV
dNEU1 ~ selfNEU1*BL_ICV + couplingb*BL_TLFB_CigCB
dNEU2 ~ selfNEU2*FU2_ICV + couplingb*FU2_TLFB_CigCB

BL_TLFB_CigCB ~~ BL_ICV         # Covariance between cognition and CigCB at T1 

#Covariance between latent change scores at each timepoint (set to equivalence)
dCIGCB1 ~~ a * dNEU1
dCIGCB2 ~~ a * dNEU2


# add sites of recruitment

BL_ICV~all_site
FU2_ICV~all_site
FU3_ICV~all_site
BL_TLFB_CigCB~all_site
FU2_TLFB_CigCB~all_site
FU3_TLFB_CigCB~all_site
all_site~~all_site
all_site~1

# add sex 
BL_ICV~sex
FU2_ICV~sex
FU3_ICV~sex
BL_TLFB_CigCB~sex
FU2_TLFB_CigCB~sex
FU3_TLFB_CigCB~sex
sex~~sex
sex~1

# add age
BL_ICV~BL_age
FU2_ICV~BL_age
FU3_ICV~BL_age
BL_TLFB_CigCB~BL_age
FU2_TLFB_CigCB~BL_age
FU3_TLFB_CigCB~BL_age
BL_age~~BL_age
BL_age~1
'
icv_out <- growth(icv, data=LCSM_DATA, estimator='mlr',fixed.x=FALSE,missing='fiml')
fitMeasures(icv_out, c("CFI","rmsea","srmr"))
summary(icv_out, fit.measures=TRUE) 




# Latent change score model: sustained attention network strength of positive network derived from Go trials & CigCB

SA_GT_PosNetStr <- '

#####     The following lines specify the core assumptions of the LCS 
#####     and should not generally be modified

# Fixes regressions between CigCB measurements at each timepoint to 1
FU2_TLFB_CigCB ~ 1*BL_TLFB_CigCB
FU3_TLFB_CigCB ~ 1*FU2_TLFB_CigCB

FU2_SA_GT_PosNetStr ~ 1*BL_SA_GT_PosNetStr     
FU3_SA_GT_PosNetStr ~ 1*FU2_SA_GT_PosNetStr

# Fixes fator loadings on latent change scores to 1
dCIGCB1 =~ 1*FU2_TLFB_CigCB      
dCIGCB2 =~ 1*FU3_TLFB_CigCB

dNEU1 =~ 1*FU2_SA_GT_PosNetStr      
dNEU2 =~ 1*FU3_SA_GT_PosNetStr

# Constrains intercepts at t2 and t3 to 0 
FU2_TLFB_CigCB ~ 0*1            
FU3_TLFB_CigCB ~ 0*1

FU2_SA_GT_PosNetStr ~ 0*1           
FU3_SA_GT_PosNetStr ~ 0*1

 # This fixes the variance of T2 and T3 measures to 0 
FU2_TLFB_CigCB ~~ 0*FU2_TLFB_CigCB      
FU3_TLFB_CigCB ~~ 0*FU3_TLFB_CigCB

FU2_SA_GT_PosNetStr ~~ 0*FU2_SA_GT_PosNetStr       
FU3_SA_GT_PosNetStr ~~ 0*FU3_SA_GT_PosNetStr

###### The following parameters will be estimated in the model.

# Estimates the intercept of CigCB/NEU at T1 (What are mean scores at T1?)
BL_TLFB_CigCB ~ 1      
BL_SA_GT_PosNetStr ~ 1 

# This estimates the variance of the change scores (How much variation is there in improvement?)
dCIGCB1 ~~  vardCigCB1*dCIGCB1        #  variances constrained to equality across waves
dCIGCB2 ~~  vardCigCB1*dCIGCB2
dNEU1 ~~  vardNEU1*dNEU1       
dNEU2 ~~  vardNEU1*dNEU2

# Estimates intercepts for latent change scores (equality constraints across timepoint with labels)
dCIGCB1 ~ CigCBint1 * 1   
dCIGCB2 ~ CigCBint2 * 1 # exception to equality constraint

dNEU1 ~ NEUint1 * 1   
dNEU2 ~ NEUint2 * 1  # exception to equality constraint

# This estimates the variance of T1 scores (How much variation is there at T1?)
BL_TLFB_CigCB ~~ BL_TLFB_CigCB      
BL_SA_GT_PosNetStr ~~ BL_SA_GT_PosNetStr 

# The following estimates the self-feedback parameters (Does change in NEU depend on NEU starting point?) and coupling parameters (Does change in NEU depend on CigCB starting point?)

# Feedback parameters allowed to vary across timepoints:

dCIGCB1 ~ selfCigCB1*BL_TLFB_CigCB + couplinga*BL_SA_GT_PosNetStr  
dCIGCB2 ~ selfCigCB2*FU2_TLFB_CigCB + couplinga*FU2_SA_GT_PosNetStr
dNEU1 ~ selfNEU1*BL_SA_GT_PosNetStr + couplingb*BL_TLFB_CigCB
dNEU2 ~ selfNEU2*FU2_SA_GT_PosNetStr + couplingb*FU2_TLFB_CigCB

BL_TLFB_CigCB ~~ BL_SA_GT_PosNetStr         # Covariance between NEU and CigCB at T1 

#Covariance between latent change scores at each timepoint (set to equivalence)
dCIGCB1 ~~ a * dNEU1
dCIGCB2 ~~ a * dNEU2


# add sites of recruitment

BL_SA_GT_PosNetStr~all_site
FU2_SA_GT_PosNetStr~all_site
FU3_SA_GT_PosNetStr~all_site
BL_TLFB_CigCB~all_site
FU2_TLFB_CigCB~all_site
FU3_TLFB_CigCB~all_site
all_site~~all_site
all_site~1

# add sex 
BL_SA_GT_PosNetStr~sex
FU2_SA_GT_PosNetStr~sex
FU3_SA_GT_PosNetStr~sex
BL_TLFB_CigCB~sex
FU2_TLFB_CigCB~sex
FU3_TLFB_CigCB~sex
sex~~sex
sex~1

# add age
BL_SA_GT_PosNetStr~BL_age
FU2_SA_GT_PosNetStr~BL_age
FU3_SA_GT_PosNetStr~BL_age
BL_TLFB_CigCB~BL_age
FU2_TLFB_CigCB~BL_age
FU3_TLFB_CigCB~BL_age
BL_age~~BL_age
BL_age~1

'
SA_GT_PosNetStr_out <- growth(SA_GT_PosNetStr, data=LCSM_DATA, estimator='mlr',fixed.x=FALSE,missing='fiml')
fitMeasures(SA_GT_PosNetStr_out, c("CFI","rmsea","srmr"))
summary(SA_GT_PosNetStr_out, fit.measures=TRUE) 


# Latent change score model: sustained attention network strength of negative network derived from Go trials & CigCB

SA_GT_NegNetStr <- '

#####     The following lines specify the core assumptions of the LCS 
#####     and should not generally be modified

# Fixes regressions between CigCB measurements at each timepoint to 1
FU2_TLFB_CigCB ~ 1*BL_TLFB_CigCB
FU3_TLFB_CigCB ~ 1*FU2_TLFB_CigCB

FU2_SA_GT_NegNetStr ~ 1*BL_SA_GT_NegNetStr     
FU3_SA_GT_NegNetStr ~ 1*FU2_SA_GT_NegNetStr

# Fixes fator loadings on latent change scores to 1
dCIGCB1 =~ 1*FU2_TLFB_CigCB      
dCIGCB2 =~ 1*FU3_TLFB_CigCB

dNEU1 =~ 1*FU2_SA_GT_NegNetStr      
dNEU2 =~ 1*FU3_SA_GT_NegNetStr

# Constrains intercepts at t2 and t3 to 0 
FU2_TLFB_CigCB ~ 0*1            
FU3_TLFB_CigCB ~ 0*1

FU2_SA_GT_NegNetStr ~ 0*1           
FU3_SA_GT_NegNetStr ~ 0*1

 # This fixes the variance of T2 and T3 measures to 0 
FU2_TLFB_CigCB ~~ 0*FU2_TLFB_CigCB      
FU3_TLFB_CigCB ~~ 0*FU3_TLFB_CigCB

FU2_SA_GT_NegNetStr ~~ 0*FU2_SA_GT_NegNetStr       
FU3_SA_GT_NegNetStr ~~ 0*FU3_SA_GT_NegNetStr

###### The following parameters will be estimated in the model.

# Estimates the intercept of CigCB/NEU at T1 (What are mean scores at T1?)
BL_TLFB_CigCB ~ 1      
BL_SA_GT_NegNetStr ~ 1 

# This estimates the variance of the change scores (How much variation is there in improvement?)
dCIGCB1 ~~  vardCigCB1*dCIGCB1        #  variances constrained to equality across waves
dCIGCB2 ~~  vardCigCB1*dCIGCB2
dNEU1 ~~  vardNEU1*dNEU1       
dNEU2 ~~  vardNEU1*dNEU2

# Estimates intercepts for latent change scores (equality constraints across timepoint with labels)
dCIGCB1 ~ CigCBint1 * 1   
dCIGCB2 ~ CigCBint2 * 1 # exception to equality constraint

dNEU1 ~ NEUint1 * 1   
dNEU2 ~ NEUint2 * 1  # exception to equality constraint

# This estimates the variance of T1 scores (How much variation is there at T1?)
BL_TLFB_CigCB ~~ BL_TLFB_CigCB      
BL_SA_GT_NegNetStr ~~ BL_SA_GT_NegNetStr 

# The following estimates the self-feedback parameters (Does change in NEU depend on NEU starting point?) and coupling parameters (Does change in NEU depend on CigCB starting point?)

# Feedback parameters allowed to vary across timepoints:

dCIGCB1 ~ selfCigCB1*BL_TLFB_CigCB + couplinga*BL_SA_GT_NegNetStr  
dCIGCB2 ~ selfCigCB2*FU2_TLFB_CigCB + couplinga*FU2_SA_GT_NegNetStr
dNEU1 ~ selfNEU1*BL_SA_GT_NegNetStr + couplingb*BL_TLFB_CigCB
dNEU2 ~ selfNEU2*FU2_SA_GT_NegNetStr + couplingb*FU2_TLFB_CigCB

BL_TLFB_CigCB ~~ BL_SA_GT_NegNetStr         # Covariance between NEU and CigCB at T1 

#Covariance between latent change scores at each timepoint (set to equivalence)
dCIGCB1 ~~ a * dNEU1
dCIGCB2 ~~ a * dNEU2


# add sites of recruitment

BL_SA_GT_NegNetStr~all_site
FU2_SA_GT_NegNetStr~all_site
FU3_SA_GT_NegNetStr~all_site
BL_TLFB_CigCB~all_site
FU2_TLFB_CigCB~all_site
FU3_TLFB_CigCB~all_site
all_site~~all_site
all_site~1

# add sex 
BL_SA_GT_NegNetStr~sex
FU2_SA_GT_NegNetStr~sex
FU3_SA_GT_NegNetStr~sex
BL_TLFB_CigCB~sex
FU2_TLFB_CigCB~sex
FU3_TLFB_CigCB~sex
sex~~sex
sex~1

# add age
BL_SA_GT_NegNetStr~BL_age
FU2_SA_GT_NegNetStr~BL_age
FU3_SA_GT_NegNetStr~BL_age
BL_TLFB_CigCB~BL_age
FU2_TLFB_CigCB~BL_age
FU3_TLFB_CigCB~BL_age
BL_age~~BL_age
BL_age~1

'
SA_GT_NegNetStr_out <- growth(SA_GT_NegNetStr, data=LCSM_DATA, estimator='mlr',fixed.x=FALSE,missing='fiml')
fitMeasures(SA_GT_NegNetStr_out, c("CFI","rmsea","srmr"))
summary(SA_GT_NegNetStr_out, fit.measures=TRUE) 


# Latent change score model: sustained attention network strength of positive network derived from Successful Stop trials & CigCB

SA_SS_PosNetStr <- '

#####     The following lines specify the core assumptions of the LCS 
#####     and should not generally be modified

# Fixes regressions between CigCB measurements at each timepoint to 1
FU2_TLFB_CigCB ~ 1*BL_TLFB_CigCB
FU3_TLFB_CigCB ~ 1*FU2_TLFB_CigCB

FU2_SA_SS_PosNetStr ~ 1*BL_SA_SS_PosNetStr     
FU3_SA_SS_PosNetStr ~ 1*FU2_SA_SS_PosNetStr

# Fixes fator loadings on latent change scores to 1
dCIGCB1 =~ 1*FU2_TLFB_CigCB      
dCIGCB2 =~ 1*FU3_TLFB_CigCB

dNEU1 =~ 1*FU2_SA_SS_PosNetStr      
dNEU2 =~ 1*FU3_SA_SS_PosNetStr

# Constrains intercepts at t2 and t3 to 0 
FU2_TLFB_CigCB ~ 0*1            
FU3_TLFB_CigCB ~ 0*1

FU2_SA_SS_PosNetStr ~ 0*1           
FU3_SA_SS_PosNetStr ~ 0*1

 # This fixes the variance of T2 and T3 measures to 0 
FU2_TLFB_CigCB ~~ 0*FU2_TLFB_CigCB      
FU3_TLFB_CigCB ~~ 0*FU3_TLFB_CigCB

FU2_SA_SS_PosNetStr ~~ 0*FU2_SA_SS_PosNetStr       
FU3_SA_SS_PosNetStr ~~ 0*FU3_SA_SS_PosNetStr

###### The following parameters will be estimated in the model.

# Estimates the intercept of CigCB/NEU at T1 (What are mean scores at T1?)
BL_TLFB_CigCB ~ 1      
BL_SA_SS_PosNetStr ~ 1 

# This estimates the variance of the change scores (How much variation is there in improvement?)
dCIGCB1 ~~  vardCigCB1*dCIGCB1        #  variances constrained to equality across waves
dCIGCB2 ~~  vardCigCB1*dCIGCB2
dNEU1 ~~  vardNEU1*dNEU1       
dNEU2 ~~  vardNEU1*dNEU2

# Estimates intercepts for latent change scores (equality constraints across timepoint with labels)
dCIGCB1 ~ CigCBint1 * 1   
dCIGCB2 ~ CigCBint2 * 1 # exception to equality constraint

dNEU1 ~ NEUint1 * 1   
dNEU2 ~ NEUint2 * 1  # exception to equality constraint

# This estimates the variance of T1 scores (How much variation is there at T1?)
BL_TLFB_CigCB ~~ BL_TLFB_CigCB      
BL_SA_SS_PosNetStr ~~ BL_SA_SS_PosNetStr 

# The following estimates the self-feedback parameters (Does change in NEU depend on NEU starting point?) and coupling parameters (Does change in NEU depend on CigCB starting point?)

# Feedback parameters allowed to vary across timepoints:

dCIGCB1 ~ selfCigCB1*BL_TLFB_CigCB + couplinga*BL_SA_SS_PosNetStr  
dCIGCB2 ~ selfCigCB2*FU2_TLFB_CigCB + couplinga*FU2_SA_SS_PosNetStr
dNEU1 ~ selfNEU1*BL_SA_SS_PosNetStr + couplingb*BL_TLFB_CigCB
dNEU2 ~ selfNEU2*FU2_SA_SS_PosNetStr + couplingb*FU2_TLFB_CigCB

BL_TLFB_CigCB ~~ BL_SA_SS_PosNetStr         # Covariance between NEU and CigCB at T1 

#Covariance between latent change scores at each timepoint (set to equivalence)
dCIGCB1 ~~ a * dNEU1
dCIGCB2 ~~ a * dNEU2



# add sites of recruitment

BL_SA_SS_PosNetStr~all_site
FU2_SA_SS_PosNetStr~all_site
FU3_SA_SS_PosNetStr~all_site
BL_TLFB_CigCB~all_site
FU2_TLFB_CigCB~all_site
FU3_TLFB_CigCB~all_site
all_site~~all_site
all_site~1

# add sex 
BL_SA_SS_PosNetStr~sex
FU2_SA_SS_PosNetStr~sex
FU3_SA_SS_PosNetStr~sex
BL_TLFB_CigCB~sex
FU2_TLFB_CigCB~sex
FU3_TLFB_CigCB~sex
sex~~sex
sex~1

# add age
BL_SA_SS_PosNetStr~BL_age
FU2_SA_SS_PosNetStr~BL_age
FU3_SA_SS_PosNetStr~BL_age
BL_TLFB_CigCB~BL_age
FU2_TLFB_CigCB~BL_age
FU3_TLFB_CigCB~BL_age
BL_age~~BL_age
BL_age~1

'
SA_SS_PosNetStr_out <- growth(SA_SS_PosNetStr, data=LCSM_DATA, estimator='mlr',fixed.x=FALSE,missing='fiml')
fitMeasures(SA_SS_PosNetStr_out, c("CFI","rmsea","srmr"))
summary(SA_SS_PosNetStr_out, fit.measures=TRUE) 


# Latent change score model: sustained attention network strength of negative network derived from Successful Stop trials & CigCB


SA_SS_NegNetStr <- '

#####     The following lines specify the core assumptions of the LCS 
#####     and should not generally be modified

# Fixes regressions between CigCB measurements at each timepoint to 1
FU2_TLFB_CigCB ~ 1*BL_TLFB_CigCB
FU3_TLFB_CigCB ~ 1*FU2_TLFB_CigCB

FU2_SA_SS_NegNetStr ~ 1*BL_SA_SS_NegNetStr     
FU3_SA_SS_NegNetStr ~ 1*FU2_SA_SS_NegNetStr

# Fixes fator loadings on latent change scores to 1
dCIGCB1 =~ 1*FU2_TLFB_CigCB      
dCIGCB2 =~ 1*FU3_TLFB_CigCB

dNEU1 =~ 1*FU2_SA_SS_NegNetStr      
dNEU2 =~ 1*FU3_SA_SS_NegNetStr

# Constrains intercepts at t2 and t3 to 0 
FU2_TLFB_CigCB ~ 0*1            
FU3_TLFB_CigCB ~ 0*1

FU2_SA_SS_NegNetStr ~ 0*1           
FU3_SA_SS_NegNetStr ~ 0*1

 # This fixes the variance of T2 and T3 measures to 0 
FU2_TLFB_CigCB ~~ 0*FU2_TLFB_CigCB      
FU3_TLFB_CigCB ~~ 0*FU3_TLFB_CigCB

FU2_SA_SS_NegNetStr ~~ 0*FU2_SA_SS_NegNetStr       
FU3_SA_SS_NegNetStr ~~ 0*FU3_SA_SS_NegNetStr

###### The following parameters will be estimated in the model.

# Estimates the intercept of CigCB/NEU at T1 (What are mean scores at T1?)
BL_TLFB_CigCB ~ 1      
BL_SA_SS_NegNetStr ~ 1 

# This estimates the variance of the change scores (How much variation is there in improvement?)
dCIGCB1 ~~  vardCigCB1*dCIGCB1        #  variances constrained to equality across waves
dCIGCB2 ~~  vardCigCB1*dCIGCB2
dNEU1 ~~  vardNEU1*dNEU1       
dNEU2 ~~  vardNEU1*dNEU2

# Estimates intercepts for latent change scores (equality constraints across timepoint with labels)
dCIGCB1 ~ CigCBint1 * 1   
dCIGCB2 ~ CigCBint2 * 1 # exception to equality constraint

dNEU1 ~ NEUint1 * 1   
dNEU2 ~ NEUint2 * 1  # exception to equality constraint

# This estimates the variance of T1 scores (How much variation is there at T1?)
BL_TLFB_CigCB ~~ BL_TLFB_CigCB      
BL_SA_SS_NegNetStr ~~ BL_SA_SS_NegNetStr 

# The following estimates the self-feedback parameters (Does change in NEU depend on NEU starting point?) and coupling parameters (Does change in NEU depend on CigCB starting point?)

# Feedback parameters allowed to vary across timepoints:

dCIGCB1 ~ selfCigCB1*BL_TLFB_CigCB + couplinga*BL_SA_SS_NegNetStr  
dCIGCB2 ~ selfCigCB2*FU2_TLFB_CigCB + couplinga*FU2_SA_SS_NegNetStr
dNEU1 ~ selfNEU1*BL_SA_SS_NegNetStr + couplingb*BL_TLFB_CigCB
dNEU2 ~ selfNEU2*FU2_SA_SS_NegNetStr + couplingb*FU2_TLFB_CigCB

BL_TLFB_CigCB ~~ BL_SA_SS_NegNetStr         # Covariance between NEU and CigCB at T1 

#Covariance between latent change scores at each timepoint (set to equivalence)
dCIGCB1 ~~ a * dNEU1
dCIGCB2 ~~ a * dNEU2



# add sites of recruitment

BL_SA_SS_NegNetStr~all_site
FU2_SA_SS_NegNetStr~all_site
FU3_SA_SS_NegNetStr~all_site
BL_TLFB_CigCB~all_site
FU2_TLFB_CigCB~all_site
FU3_TLFB_CigCB~all_site
all_site~~all_site
all_site~1

# add sex 
BL_SA_SS_NegNetStr~sex
FU2_SA_SS_NegNetStr~sex
FU3_SA_SS_NegNetStr~sex
BL_TLFB_CigCB~sex
FU2_TLFB_CigCB~sex
FU3_TLFB_CigCB~sex
sex~~sex
sex~1

# add age
BL_SA_SS_NegNetStr~BL_age
FU2_SA_SS_NegNetStr~BL_age
FU3_SA_SS_NegNetStr~BL_age
BL_TLFB_CigCB~BL_age
FU2_TLFB_CigCB~BL_age
FU3_TLFB_CigCB~BL_age
BL_age~~BL_age
BL_age~1

'
SA_SS_NegNetStr_out <- growth(SA_SS_NegNetStr, data=LCSM_DATA, estimator='mlr',fixed.x=FALSE,missing='fiml')
fitMeasures(SA_SS_NegNetStr_out, c("CFI","rmsea","srmr"))
summary(SA_SS_NegNetStr_out, fit.measures=TRUE) 




################## alcohol use #################


# Latent change score model: ICV & alcohol use
icv <- '

#####     The following lines specify the core assumptions of the LCS 
#####     and should not generally be modified

# Fixes regressions between CigCB measurements at each timepoint to 1
FU2_TLFB_alc ~ 1*BL_TLFB_alc
FU3_TLFB_alc ~ 1*FU2_TLFB_alc

FU2_ICV ~ 1*BL_ICV     
FU3_ICV ~ 1*FU2_ICV

# Fixes fator loadings on latent change scores to 1
dCIGCB1 =~ 1*FU2_TLFB_alc      
dCIGCB2 =~ 1*FU3_TLFB_alc

dNEU1 =~ 1*FU2_ICV      
dNEU2 =~ 1*FU3_ICV

# Constrains intercepts at t2 and t3 to 0 
FU2_TLFB_alc ~ 0*1            
FU3_TLFB_alc ~ 0*1

FU2_ICV ~ 0*1           
FU3_ICV ~ 0*1

 # This fixes the variance of T2 and T3 measures to 0 
FU2_TLFB_alc ~~ 0*FU2_TLFB_alc      
FU3_TLFB_alc ~~ 0*FU3_TLFB_alc

FU2_ICV ~~ 0*FU2_ICV       
FU3_ICV ~~ 0*FU3_ICV

###### The following parameters will be estimated in the model.

# Estimates the intercept of CigCB/COGNITION at T3 (What are mean scores at T1?)
BL_TLFB_alc ~ 1      
BL_ICV ~ 1 

# This estimates the variance of the change scores (How much variation is there in improvement?)
dCIGCB1 ~~  vardCigCB1*dCIGCB1        #  variances constrained to equality across waves
dCIGCB2 ~~  vardCigCB1*dCIGCB2
dNEU1 ~~  vardNEU1*dNEU1       
dNEU2 ~~  vardNEU1*dNEU2

# Estimates intercepts for latent change scores (equality constraints across timepoint with labels)
dCIGCB1 ~ CigCBint1 * 1   
dCIGCB2 ~ CigCBint2 * 1 # exception to equality constraint

dNEU1 ~ NEUint1 * 1   
dNEU2 ~ NEUint2 * 1  # exception to equality constraint

# This estimates the variance of T1 scores (How much variation is there at T1?)
BL_TLFB_alc ~~ BL_TLFB_alc      
BL_ICV ~~ BL_ICV 

# The following estimates the self-feedback parameters (Does change in cognition depend on cognition starting point?) and coupling parameters (Does change in cognition depend on block design starting point?)

# Feedback parameters allowed to vary across timepoints:

dCIGCB1 ~ selfCigCB1*BL_TLFB_alc + couplinga*BL_ICV  
dCIGCB2 ~ selfCigCB2*FU2_TLFB_alc + couplinga*FU2_ICV
dNEU1 ~ selfNEU1*BL_ICV + couplingb*BL_TLFB_alc
dNEU2 ~ selfNEU2*FU2_ICV + couplingb*FU2_TLFB_alc

BL_TLFB_alc ~~ BL_ICV         # Covariance between cognition and CigCB at T1 

#Covariance between latent change scores at each timepoint (set to equivalence)
dCIGCB1 ~~ a * dNEU1
dCIGCB2 ~~ a * dNEU2



# add sites of recruitment

BL_ICV~all_site
FU2_ICV~all_site
FU3_ICV~all_site
BL_TLFB_alc~all_site
FU2_TLFB_alc~all_site
FU3_TLFB_alc~all_site
all_site~~all_site
all_site~1

# add sex 
BL_ICV~sex
FU2_ICV~sex
FU3_ICV~sex
BL_TLFB_alc~sex
FU2_TLFB_alc~sex
FU3_TLFB_alc~sex
sex~~sex
sex~1

# add age
BL_ICV~BL_age
FU2_ICV~BL_age
FU3_ICV~BL_age
BL_TLFB_alc~BL_age
FU2_TLFB_alc~BL_age
FU3_TLFB_alc~BL_age
BL_age~~BL_age
BL_age~1

'
icv_out <- growth(icv, data=LCSM_DATA, estimator='mlr',fixed.x=FALSE,missing='fiml')
fitMeasures(icv_out, c("CFI","rmsea","srmr"))
summary(icv_out, fit.measures=TRUE) 



# Latent change score model: sustained attention network strength of positive network derived from Go trials & alcohol use
SA_GT_PosNetStr <- '

#####     The following lines specify the core assumptions of the LCS 
#####     and should not generally be modified

# Fixes regressions between CigCB measurements at each timepoint to 1
FU2_TLFB_alc ~ 1*BL_TLFB_alc
FU3_TLFB_alc ~ 1*FU2_TLFB_alc

FU2_SA_GT_PosNetStr ~ 1*BL_SA_GT_PosNetStr     
FU3_SA_GT_PosNetStr ~ 1*FU2_SA_GT_PosNetStr

# Fixes fator loadings on latent change scores to 1
dCIGCB1 =~ 1*FU2_TLFB_alc      
dCIGCB2 =~ 1*FU3_TLFB_alc

dNEU1 =~ 1*FU2_SA_GT_PosNetStr      
dNEU2 =~ 1*FU3_SA_GT_PosNetStr

# Constrains intercepts at t2 and t3 to 0 
FU2_TLFB_alc ~ 0*1            
FU3_TLFB_alc ~ 0*1

FU2_SA_GT_PosNetStr ~ 0*1           
FU3_SA_GT_PosNetStr ~ 0*1

 # This fixes the variance of T2 and T3 measures to 0 
FU2_TLFB_alc ~~ 0*FU2_TLFB_alc      
FU3_TLFB_alc ~~ 0*FU3_TLFB_alc

FU2_SA_GT_PosNetStr ~~ 0*FU2_SA_GT_PosNetStr       
FU3_SA_GT_PosNetStr ~~ 0*FU3_SA_GT_PosNetStr

###### The following parameters will be estimated in the model.

# Estimates the intercept of CigCB/NEU at T1 (What are mean scores at T1?)
BL_TLFB_alc ~ 1      
BL_SA_GT_PosNetStr ~ 1 

# This estimates the variance of the change scores (How much variation is there in improvement?)
dCIGCB1 ~~  vardCigCB1*dCIGCB1        #  variances constrained to equality across waves
dCIGCB2 ~~  vardCigCB1*dCIGCB2
dNEU1 ~~  vardNEU1*dNEU1       
dNEU2 ~~  vardNEU1*dNEU2

# Estimates intercepts for latent change scores (equality constraints across timepoint with labels)
dCIGCB1 ~ CigCBint1 * 1   
dCIGCB2 ~ CigCBint2 * 1 # exception to equality constraint

dNEU1 ~ NEUint1 * 1   
dNEU2 ~ NEUint2 * 1  # exception to equality constraint

# This estimates the variance of T1 scores (How much variation is there at T1?)
BL_TLFB_alc ~~ BL_TLFB_alc      
BL_SA_GT_PosNetStr ~~ BL_SA_GT_PosNetStr 

# The following estimates the self-feedback parameters (Does change in NEU depend on NEU starting point?) and coupling parameters (Does change in NEU depend on CigCB starting point?)

# Feedback parameters allowed to vary across timepoints:

dCIGCB1 ~ selfCigCB1*BL_TLFB_alc + couplinga*BL_SA_GT_PosNetStr  
dCIGCB2 ~ selfCigCB2*FU2_TLFB_alc + couplinga*FU2_SA_GT_PosNetStr
dNEU1 ~ selfNEU1*BL_SA_GT_PosNetStr + couplingb*BL_TLFB_alc
dNEU2 ~ selfNEU2*FU2_SA_GT_PosNetStr + couplingb*FU2_TLFB_alc

BL_TLFB_alc ~~ BL_SA_GT_PosNetStr         # Covariance between NEU and CigCB at T1 

#Covariance between latent change scores at each timepoint (set to equivalence)
dCIGCB1 ~~ a * dNEU1
dCIGCB2 ~~ a * dNEU2



# add sites of recruitment

BL_SA_GT_PosNetStr~all_site
FU2_SA_GT_PosNetStr~all_site
FU3_SA_GT_PosNetStr~all_site
BL_TLFB_alc~all_site
FU2_TLFB_alc~all_site
FU3_TLFB_alc~all_site
all_site~~all_site
all_site~1

# add sex 
BL_SA_GT_PosNetStr~sex
FU2_SA_GT_PosNetStr~sex
FU3_SA_GT_PosNetStr~sex
BL_TLFB_alc~sex
FU2_TLFB_alc~sex
FU3_TLFB_alc~sex
sex~~sex
sex~1

# add age
BL_SA_GT_PosNetStr~BL_age
FU2_SA_GT_PosNetStr~BL_age
FU3_SA_GT_PosNetStr~BL_age
BL_TLFB_alc~BL_age
FU2_TLFB_alc~BL_age
FU3_TLFB_alc~BL_age
BL_age~~BL_age
BL_age~1

'
SA_GT_PosNetStr_out <- growth(SA_GT_PosNetStr, data=LCSM_DATA, estimator='mlr',fixed.x=FALSE,missing='fiml')
fitMeasures(SA_GT_PosNetStr_out, c("CFI","rmsea","srmr"))
summary(SA_GT_PosNetStr_out, fit.measures=TRUE) 

# Latent change score model: sustained attention network strength of negative network derived from Go trials & alcohol use

SA_GT_NegNetStr <- '

#####     The following lines specify the core assumptions of the LCS 
#####     and should not generally be modified

# Fixes regressions between CigCB measurements at each timepoint to 1
FU2_TLFB_alc ~ 1*BL_TLFB_alc
FU3_TLFB_alc ~ 1*FU2_TLFB_alc

FU2_SA_GT_NegNetStr ~ 1*BL_SA_GT_NegNetStr     
FU3_SA_GT_NegNetStr ~ 1*FU2_SA_GT_NegNetStr

# Fixes fator loadings on latent change scores to 1
dCIGCB1 =~ 1*FU2_TLFB_alc      
dCIGCB2 =~ 1*FU3_TLFB_alc

dNEU1 =~ 1*FU2_SA_GT_NegNetStr      
dNEU2 =~ 1*FU3_SA_GT_NegNetStr

# Constrains intercepts at t2 and t3 to 0 
FU2_TLFB_alc ~ 0*1            
FU3_TLFB_alc ~ 0*1

FU2_SA_GT_NegNetStr ~ 0*1           
FU3_SA_GT_NegNetStr ~ 0*1

 # This fixes the variance of T2 and T3 measures to 0 
FU2_TLFB_alc ~~ 0*FU2_TLFB_alc      
FU3_TLFB_alc ~~ 0*FU3_TLFB_alc

FU2_SA_GT_NegNetStr ~~ 0*FU2_SA_GT_NegNetStr       
FU3_SA_GT_NegNetStr ~~ 0*FU3_SA_GT_NegNetStr

###### The following parameters will be estimated in the model.

# Estimates the intercept of CigCB/NEU at T1 (What are mean scores at T1?)
BL_TLFB_alc ~ 1      
BL_SA_GT_NegNetStr ~ 1 

# This estimates the variance of the change scores (How much variation is there in improvement?)
dCIGCB1 ~~  vardCigCB1*dCIGCB1        #  variances constrained to equality across waves
dCIGCB2 ~~  vardCigCB1*dCIGCB2
dNEU1 ~~  vardNEU1*dNEU1       
dNEU2 ~~  vardNEU1*dNEU2

# Estimates intercepts for latent change scores (equality constraints across timepoint with labels)
dCIGCB1 ~ CigCBint1 * 1   
dCIGCB2 ~ CigCBint2 * 1 # exception to equality constraint

dNEU1 ~ NEUint1 * 1   
dNEU2 ~ NEUint2 * 1  # exception to equality constraint

# This estimates the variance of T1 scores (How much variation is there at T1?)
BL_TLFB_alc ~~ BL_TLFB_alc      
BL_SA_GT_NegNetStr ~~ BL_SA_GT_NegNetStr 

# The following estimates the self-feedback parameters (Does change in NEU depend on NEU starting point?) and coupling parameters (Does change in NEU depend on CigCB starting point?)

# Feedback parameters allowed to vary across timepoints:

dCIGCB1 ~ selfCigCB1*BL_TLFB_alc + couplinga*BL_SA_GT_NegNetStr  
dCIGCB2 ~ selfCigCB2*FU2_TLFB_alc + couplinga*FU2_SA_GT_NegNetStr
dNEU1 ~ selfNEU1*BL_SA_GT_NegNetStr + couplingb*BL_TLFB_alc
dNEU2 ~ selfNEU2*FU2_SA_GT_NegNetStr + couplingb*FU2_TLFB_alc

BL_TLFB_alc ~~ BL_SA_GT_NegNetStr         # Covariance between NEU and CigCB at T1 

#Covariance between latent change scores at each timepoint (set to equivalence)
dCIGCB1 ~~ a * dNEU1
dCIGCB2 ~~ a * dNEU2


# add sites of recruitment

BL_SA_GT_NegNetStr~all_site
FU2_SA_GT_NegNetStr~all_site
FU3_SA_GT_NegNetStr~all_site
BL_TLFB_alc~all_site
FU2_TLFB_alc~all_site
FU3_TLFB_alc~all_site
all_site~~all_site
all_site~1

# add sex 
BL_SA_GT_NegNetStr~sex
FU2_SA_GT_NegNetStr~sex
FU3_SA_GT_NegNetStr~sex
BL_TLFB_alc~sex
FU2_TLFB_alc~sex
FU3_TLFB_alc~sex
sex~~sex
sex~1

# add age
BL_SA_GT_NegNetStr~BL_age
FU2_SA_GT_NegNetStr~BL_age
FU3_SA_GT_NegNetStr~BL_age
BL_TLFB_alc~BL_age
FU2_TLFB_alc~BL_age
FU3_TLFB_alc~BL_age
BL_age~~BL_age
BL_age~1

'
SA_GT_NegNetStr_out <- growth(SA_GT_NegNetStr, data=LCSM_DATA, estimator='mlr',fixed.x=FALSE,missing='fiml')
fitMeasures(SA_GT_NegNetStr_out, c("CFI","rmsea","srmr"))
summary(SA_GT_NegNetStr_out, fit.measures=TRUE) 


# Latent change score model: sustained attention network strength of positive network derived from Successful stop trials & alcohol use

SA_SS_PosNetStr <- '

#####     The following lines specify the core assumptions of the LCS 
#####     and should not generally be modified

# Fixes regressions between CigCB measurements at each timepoint to 1
FU2_TLFB_alc ~ 1*BL_TLFB_alc
FU3_TLFB_alc ~ 1*FU2_TLFB_alc

FU2_SA_SS_PosNetStr ~ 1*BL_SA_SS_PosNetStr     
FU3_SA_SS_PosNetStr ~ 1*FU2_SA_SS_PosNetStr

# Fixes fator loadings on latent change scores to 1
dCIGCB1 =~ 1*FU2_TLFB_alc      
dCIGCB2 =~ 1*FU3_TLFB_alc

dNEU1 =~ 1*FU2_SA_SS_PosNetStr      
dNEU2 =~ 1*FU3_SA_SS_PosNetStr

# Constrains intercepts at t2 and t3 to 0 
FU2_TLFB_alc ~ 0*1            
FU3_TLFB_alc ~ 0*1

FU2_SA_SS_PosNetStr ~ 0*1           
FU3_SA_SS_PosNetStr ~ 0*1

 # This fixes the variance of T2 and T3 measures to 0 
FU2_TLFB_alc ~~ 0*FU2_TLFB_alc      
FU3_TLFB_alc ~~ 0*FU3_TLFB_alc

FU2_SA_SS_PosNetStr ~~ 0*FU2_SA_SS_PosNetStr       
FU3_SA_SS_PosNetStr ~~ 0*FU3_SA_SS_PosNetStr

###### The following parameters will be estimated in the model.

# Estimates the intercept of CigCB/NEU at T1 (What are mean scores at T1?)
BL_TLFB_alc ~ 1      
BL_SA_SS_PosNetStr ~ 1 

# This estimates the variance of the change scores (How much variation is there in improvement?)
dCIGCB1 ~~  vardCigCB1*dCIGCB1        #  variances constrained to equality across waves
dCIGCB2 ~~  vardCigCB1*dCIGCB2
dNEU1 ~~  vardNEU1*dNEU1       
dNEU2 ~~  vardNEU1*dNEU2

# Estimates intercepts for latent change scores (equality constraints across timepoint with labels)
dCIGCB1 ~ CigCBint1 * 1   
dCIGCB2 ~ CigCBint2 * 1 # exception to equality constraint

dNEU1 ~ NEUint1 * 1   
dNEU2 ~ NEUint2 * 1  # exception to equality constraint

# This estimates the variance of T1 scores (How much variation is there at T1?)
BL_TLFB_alc ~~ BL_TLFB_alc      
BL_SA_SS_PosNetStr ~~ BL_SA_SS_PosNetStr 

# The following estimates the self-feedback parameters (Does change in NEU depend on NEU starting point?) and coupling parameters (Does change in NEU depend on CigCB starting point?)

# Feedback parameters allowed to vary across timepoints:

dCIGCB1 ~ selfCigCB1*BL_TLFB_alc + couplinga*BL_SA_SS_PosNetStr  
dCIGCB2 ~ selfCigCB2*FU2_TLFB_alc + couplinga*FU2_SA_SS_PosNetStr
dNEU1 ~ selfNEU1*BL_SA_SS_PosNetStr + couplingb*BL_TLFB_alc
dNEU2 ~ selfNEU2*FU2_SA_SS_PosNetStr + couplingb*FU2_TLFB_alc

BL_TLFB_alc ~~ BL_SA_SS_PosNetStr         # Covariance between NEU and CigCB at T1 

#Covariance between latent change scores at each timepoint (set to equivalence)
dCIGCB1 ~~ a * dNEU1
dCIGCB2 ~~ a * dNEU2



# add sites of recruitment

BL_SA_SS_PosNetStr~all_site
FU2_SA_SS_PosNetStr~all_site
FU3_SA_SS_PosNetStr~all_site
BL_TLFB_alc~all_site
FU2_TLFB_alc~all_site
FU3_TLFB_alc~all_site
all_site~~all_site
all_site~1

# add sex 
BL_SA_SS_PosNetStr~sex
FU2_SA_SS_PosNetStr~sex
FU3_SA_SS_PosNetStr~sex
BL_TLFB_alc~sex
FU2_TLFB_alc~sex
FU3_TLFB_alc~sex
sex~~sex
sex~1

# add age
BL_SA_SS_PosNetStr~BL_age
FU2_SA_SS_PosNetStr~BL_age
FU3_SA_SS_PosNetStr~BL_age
BL_TLFB_alc~BL_age
FU2_TLFB_alc~BL_age
FU3_TLFB_alc~BL_age
BL_age~~BL_age
BL_age~1

'

SA_SS_PosNetStr_out <- growth(SA_SS_PosNetStr, data=LCSM_DATA, estimator='mlr',fixed.x=FALSE,missing='fiml')
fitMeasures(SA_SS_PosNetStr_out, c("CFI","rmsea","srmr"))
summary(SA_SS_PosNetStr_out, fit.measures=TRUE) 


# Latent change score model: sustained attention network strength of negative network derived from Successful stop trials & alcohol use

SA_SS_NegNetStr <- '

#####     The following lines specify the core assumptions of the LCS 
#####     and should not generally be modified

# Fixes regressions between CigCB measurements at each timepoint to 1
FU2_TLFB_alc ~ 1*BL_TLFB_alc
FU3_TLFB_alc ~ 1*FU2_TLFB_alc

FU2_SA_SS_NegNetStr ~ 1*BL_SA_SS_NegNetStr     
FU3_SA_SS_NegNetStr ~ 1*FU2_SA_SS_NegNetStr

# Fixes fator loadings on latent change scores to 1
dCIGCB1 =~ 1*FU2_TLFB_alc      
dCIGCB2 =~ 1*FU3_TLFB_alc

dNEU1 =~ 1*FU2_SA_SS_NegNetStr      
dNEU2 =~ 1*FU3_SA_SS_NegNetStr

# Constrains intercepts at t2 and t3 to 0 
FU2_TLFB_alc ~ 0*1            
FU3_TLFB_alc ~ 0*1

FU2_SA_SS_NegNetStr ~ 0*1           
FU3_SA_SS_NegNetStr ~ 0*1

 # This fixes the variance of T2 and T3 measures to 0 
FU2_TLFB_alc ~~ 0*FU2_TLFB_alc      
FU3_TLFB_alc ~~ 0*FU3_TLFB_alc

FU2_SA_SS_NegNetStr ~~ 0*FU2_SA_SS_NegNetStr       
FU3_SA_SS_NegNetStr ~~ 0*FU3_SA_SS_NegNetStr

###### The following parameters will be estimated in the model.

# Estimates the intercept of CigCB/NEU at T1 (What are mean scores at T1?)
BL_TLFB_alc ~ 1      
BL_SA_SS_NegNetStr ~ 1 

# This estimates the variance of the change scores (How much variation is there in improvement?)
dCIGCB1 ~~  vardCigCB1*dCIGCB1        #  variances constrained to equality across waves
dCIGCB2 ~~  vardCigCB1*dCIGCB2
dNEU1 ~~  vardNEU1*dNEU1       
dNEU2 ~~  vardNEU1*dNEU2

# Estimates intercepts for latent change scores (equality constraints across timepoint with labels)
dCIGCB1 ~ CigCBint1 * 1   
dCIGCB2 ~ CigCBint2 * 1 # exception to equality constraint

dNEU1 ~ NEUint1 * 1   
dNEU2 ~ NEUint2 * 1  # exception to equality constraint

# This estimates the variance of T1 scores (How much variation is there at T1?)
BL_TLFB_alc ~~ BL_TLFB_alc      
BL_SA_SS_NegNetStr ~~ BL_SA_SS_NegNetStr 

# The following estimates the self-feedback parameters (Does change in NEU depend on NEU starting point?) and coupling parameters (Does change in NEU depend on CigCB starting point?)

# Feedback parameters allowed to vary across timepoints:

dCIGCB1 ~ selfCigCB1*BL_TLFB_alc + couplinga*BL_SA_SS_NegNetStr  
dCIGCB2 ~ selfCigCB2*FU2_TLFB_alc + couplinga*FU2_SA_SS_NegNetStr
dNEU1 ~ selfNEU1*BL_SA_SS_NegNetStr + couplingb*BL_TLFB_alc
dNEU2 ~ selfNEU2*FU2_SA_SS_NegNetStr + couplingb*FU2_TLFB_alc

BL_TLFB_alc ~~ BL_SA_SS_NegNetStr         # Covariance between NEU and CigCB at T1 

#Covariance between latent change scores at each timepoint (set to equivalence)
dCIGCB1 ~~ a * dNEU1
dCIGCB2 ~~ a * dNEU2



# add sites of recruitment

BL_SA_SS_NegNetStr~all_site
FU2_SA_SS_NegNetStr~all_site
FU3_SA_SS_NegNetStr~all_site
BL_TLFB_alc~all_site
FU2_TLFB_alc~all_site
FU3_TLFB_alc~all_site
all_site~~all_site
all_site~1

# add sex 
BL_SA_SS_NegNetStr~sex
FU2_SA_SS_NegNetStr~sex
FU3_SA_SS_NegNetStr~sex
BL_TLFB_alc~sex
FU2_TLFB_alc~sex
FU3_TLFB_alc~sex
sex~~sex
sex~1

# add age
BL_SA_SS_NegNetStr~BL_age
FU2_SA_SS_NegNetStr~BL_age
FU3_SA_SS_NegNetStr~BL_age
BL_TLFB_alc~BL_age
FU2_TLFB_alc~BL_age
FU3_TLFB_alc~BL_age
BL_age~~BL_age
BL_age~1

'
SA_SS_NegNetStr_out <- growth(SA_SS_NegNetStr, data=LCSM_DATA, estimator='mlr',fixed.x=FALSE,missing='fiml')
fitMeasures(SA_SS_NegNetStr_out, c("CFI","rmsea","srmr"))
summary(SA_SS_NegNetStr_out, fit.measures=TRUE) 




