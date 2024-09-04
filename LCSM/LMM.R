rm(list=ls(all=TRUE))
# Load packages
library(lme4)
library(data.table)
library(ggplot2)
library(emmeans)
library(lmerTest)

####################################
# 1) load data
df = read.csv(file="W:/Behavior_LMM/ICV_Behavior_table.csv");

df$age_c <- (df$BL_age-mean(df$BL_age,na.rm=TRUE))/365  # (age-mean)/365
df$sex.f <- df$sex
df$sex.f[df$sex.f==0]<- "0Male"
df$sex.f[df$sex.f==1]<- "Female"  
df$sex.f <- factor(df$sex.f)
df$site.f <- df$sites 
df$site.f <- factor(df$site.f)
df$timepoint.f <- df$timepoint
df$timepoint.f <- factor(df$timepoint.f)
df$subject.f <- df$ID
df$subject.f <- factor(df$subject.f)
df$ICV <- (df$ICV - mean(df$ICV,na.rm=TRUE))/sd(df$ICV,na.rm=TRUE)

####################################
## 2) Fit model
## LMM
md <- lmer(ICV~timepoint.f+sex.f+site.f+pds+age_c+(1|subject.f), data = df, REML = TRUE)

md_summary<- summary(md,type = 3)
(md_anova<- anova(md,type=3))
(emmeans(md, specs = pairwise ~ timepoint.f, mode = "satterthwaite"))
md_emameans <- emmeans(md, specs = pairwise ~ timepoint.f, mode = "satterthwaite")

####################################
## 3) Save data
outputfile = "W:/Behavior_LMM/"
output_filename<- paste0(outputfile, "LMM_results_ICV.txt")

output_filehandle <- file(output_filename, "a+")
cat("\n #################################################", file = output_filehandle)
cat("\n Fitted variable ","ICV~timepoint.f+sex.f+site.f+pds+age_c+(1|subject.f)", "\n", file = output_filehandle)
capture.output(md_summary, file = output_filehandle)
cat("\n", file = output_filehandle)
capture.output(md_anova, file = output_filehandle)
capture.output(md_emameans, file = output_filehandle)
close(output_filehandle)  

