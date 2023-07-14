# Rowan K. O'Hara
# 07.09.2023
# Thesis Project
#
# PRSxE
# NOTE (if run on Fenn):
# Need to run interactively on Fenn using R in the terminal
# Push to terminal using Option + Command + Return
################################################################################

# Run on Fenn

cd /vcu_gpfs2/home/oharark/S4S/PRS/dep/

R

# LOAD PACKAGES & READ IN DATA #################################################

library(tidyverse)
library(data.table)

pheno<- fread("thesis_s4s_clean.csv")
colnames(pheno)

afr <- fread("output/S4S_AFR_dep_scores.sscore")
colnames(afr) <- c("FID", "unique_id", "ALLELE_CT",	"NAMED_ALLELE_DOSAGE_SUM",	"depPRSAFRmeta")
afr <- subset(afr, select=c(unique_id,depPRSAFRmeta))

amr <- fread("output/S4S_AMR_dep_scores.sscore")
colnames(amr) <- c("FID", "unique_id", "ALLELE_CT",	"NAMED_ALLELE_DOSAGE_SUM",	"depPRSAMRmeta")
amr <- subset(amr, select=c(unique_id,depPRSAMRmeta))

eas <- fread("output/S4S_EAS_dep_scores.sscore")
colnames(eas) <- c("FID", "unique_id", "ALLELE_CT",	"NAMED_ALLELE_DOSAGE_SUM",	"depPRSEASmeta")
eas <- subset(eas, select=c(unique_id,depPRSEASmeta))

eur <- fread("output/S4S_EUR_dep_scores.sscore")
colnames(eur) <- c("FID", "unique_id", "ALLELE_CT",	"NAMED_ALLELE_DOSAGE_SUM",	"depPRSEURmeta")
eur <- subset(eur, select=c(unique_id,depPRSEURmeta))

sas <- fread("output/S4S_SAS_dep_scores.sscore")
colnames(sas) <- c("FID", "unique_id", "ALLELE_CT",	"NAMED_ALLELE_DOSAGE_SUM",	"depPRSSASmeta")
sas <- subset(sas, select=c(unique_id,depPRSSASmeta))


# put all data frames into list
df_list <- list(pheno, 
                afr,
                #amr,
                eas,
                eur,
                sas)

# merge all data frames in list
phenoPRS<- df_list %>% reduce(full_join, by='unique_id')

AFRphenoPRS<-subset(phenoPRS, ancestry=="AFR") 
AMRphenoPRS<-subset(phenoPRS, ancestry=="AMR") 
EASphenoPRS<-subset(phenoPRS, ancestry=="EAS") 
EURphenoPRS<-subset(phenoPRS, ancestry=="EUR") 
SASphenoPRS<-subset(phenoPRS, ancestry=="SAS") 


# LINEAR REGRESSION ############################################################

## EUR ----
fit.alt <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRSEURmeta + ever_ipv + depPRSEURmeta*ever_ipv, data = EURphenoPRS)
fit.alt1 <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRSEURmeta + ever_ipv + biosex*ever_ipv + depPRSEURmeta*ever_ipv, data = EURphenoPRS)

summary(fit.alt)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             6.541e+00  2.202e-01  29.702  < 2e-16 ***
#   PC1                    -2.660e+01  4.164e+00  -6.389 1.88e-10 ***
#   PC2                    -1.033e+01  3.952e+00  -2.614 0.008987 ** 
#   PC3                     5.387e+00  4.196e+00   1.284 0.199302    
#   PC4                    -1.198e+01  3.848e+00  -3.113 0.001866 ** 
#   PC5                    -1.768e+00  3.896e+00  -0.454 0.649898    
#   biosex                  1.084e+00  1.243e-01   8.725  < 2e-16 ***
#   depPRSEURmeta           1.639e+06  4.429e+05   3.701 0.000218 ***
#   ever_ipv                1.538e+00  1.630e-01   9.433  < 2e-16 ***
#   depPRSEURmeta:ever_ipv -7.997e+05  6.904e+05  -1.158 0.246868    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3.63 on 3690 degrees of freedom
# Multiple R-squared:  0.08492,   Adjusted R-squared:  0.08269 
# F-statistic: 38.05 on 9 and 3690 DF,  p-value: < 2.2e-16

summary(fit.alt1)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             6.593e+00  2.663e-01  24.758  < 2e-16 ***
#   PC1                    -2.657e+01  4.165e+00  -6.380 1.99e-10 ***
#   PC2                    -1.030e+01  3.953e+00  -2.607 0.009175 ** 
#   PC3                     5.370e+00  4.197e+00   1.279 0.200870    
#   PC4                    -1.196e+01  3.849e+00  -3.108 0.001898 ** 
#   PC5                    -1.756e+00  3.896e+00  -0.451 0.652279    
#   biosex                  1.051e+00  1.564e-01   6.723 2.06e-11 ***
#   depPRSEURmeta           1.636e+06  4.430e+05   3.692 0.000226 ***
#   ever_ipv                1.392e+00  4.508e-01   3.089 0.002024 ** 
#   biosex:ever_ipv         8.873e-02  2.564e-01   0.346 0.729291    
#   depPRSEURmeta:ever_ipv -8.017e+05  6.906e+05  -1.161 0.245758    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3.631 on 3689 degrees of freedom
# Multiple R-squared:  0.08495,   Adjusted R-squared:  0.08247 
# F-statistic: 34.25 on 10 and 3689 DF,  p-value: < 2.2e-16

fit.alt1 <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRSEURmeta + ever_ipv + biosex*ever_ipv + depPRSEURmeta*ever_ipv, data = EURphenoPRS)
fit.noSex <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + depPRSEURmeta + ever_ipv + depPRSEURmeta*ever_ipv, data = EURphenoPRS)
fit.noSxE <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRSEURmeta + ever_ipv + depPRSEURmeta*ever_ipv, data = EURphenoPRS)
fit.noG <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + ever_ipv + biosex*ever_ipv, data = EURphenoPRS)
fit.noGxE <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRSEURmeta + ever_ipv + biosex*ever_ipv, data = EURphenoPRS)
fit.noE <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRSEURmeta, data = EURphenoPRS)

rSex     <-summary(fit.alt1)$r.squared - summary(fit.noSex)$r.squared # due to Sex  0.01890889 (1.89%)
rE       <-summary(fit.alt1)$r.squared - summary(fit.noE)$r.squared # due to IPV  0.0331644 (3.32%)
rSexIPV  <-summary(fit.alt1)$r.squared - summary(fit.noSxE)$r.squared # due to SexByIPV  2.971171e-05
rPRS     <-summary(fit.alt1)$r.squared - summary(fit.noG)$r.squared # due to PRS  0.0039983 (0.4%)
rPRSxIPV <-summary(fit.alt1)$r.squared - summary(fit.noGxE)$r.squared # due to PRSxIPV  0.0003342919


## AFR ----

fit.alt <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRSAFRmeta + ever_ipv + depPRSAFRmeta*ever_ipv, data = AFRphenoPRS)
fit.alt1 <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRSAFRmeta + ever_ipv + biosex*ever_ipv + depPRSAFRmeta*ever_ipv, data = AFRphenoPRS)


summary(fit.alt)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             6.576e+00  3.793e-01  17.337  < 2e-16 ***
#   PC1                     2.918e+01  4.040e+00   7.222 8.10e-13 ***
#   PC2                    -3.123e+00  3.935e+00  -0.794   0.4275    
#   PC3                    -8.214e+00  4.158e+00  -1.976   0.0484 *  
#   PC4                     8.700e+00  4.136e+00   2.104   0.0356 *  
#   PC5                     5.109e+00  3.847e+00   1.328   0.1844    
#   biosex                  8.610e-01  2.092e-01   4.116 4.07e-05 ***
#   depPRSAFRmeta          -1.130e+05  3.952e+05  -0.286   0.7750    
#   ever_ipv                1.489e+00  1.933e-01   7.705 2.35e-14 ***
#   depPRSAFRmeta:ever_ipv -1.205e+05  6.419e+05  -0.188   0.8511    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3.591 on 1515 degrees of freedom
# Multiple R-squared:  0.08732,   Adjusted R-squared:  0.0819 
# F-statistic: 16.11 on 9 and 1515 DF,  p-value: < 2.2e-16

summary(fit.alt1)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             6.859e+00  4.656e-01  14.732  < 2e-16 ***
#   PC1                     2.904e+01  4.042e+00   7.183 1.06e-12 ***
#   PC2                    -2.970e+00  3.938e+00  -0.754  0.45075    
#   PC3                    -8.215e+00  4.158e+00  -1.976  0.04835 *  
#   PC4                     8.625e+00  4.136e+00   2.085  0.03722 *  
#   PC5                     5.259e+00  3.849e+00   1.366  0.17206    
#   biosex                  6.963e-01  2.616e-01   2.662  0.00785 ** 
#   depPRSAFRmeta          -1.192e+05  3.952e+05  -0.302  0.76295    
#   ever_ipv                6.919e-01  7.848e-01   0.882  0.37809    
#   biosex:ever_ipv         4.569e-01  4.359e-01   1.048  0.29470    
#   depPRSAFRmeta:ever_ipv -1.113e+05  6.420e+05  -0.173  0.86237    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3.591 on 1514 degrees of freedom
# Multiple R-squared:  0.08798,   Adjusted R-squared:  0.08196 
# F-statistic: 14.61 on 10 and 1514 DF,  p-value: < 2.2e-16

fit.alt1 <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRSAFRmeta + ever_ipv + biosex*ever_ipv + depPRSAFRmeta*ever_ipv, data = AFRphenoPRS)
fit.noSex <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + depPRSAFRmeta + ever_ipv + depPRSAFRmeta*ever_ipv, data = AFRphenoPRS)
fit.noSxE <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRSAFRmeta + ever_ipv + depPRSAFRmeta*ever_ipv, data = AFRphenoPRS)
fit.noG <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + ever_ipv + biosex*ever_ipv, data = AFRphenoPRS)
fit.noGxE <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRSAFRmeta + ever_ipv + biosex*ever_ipv, data = AFRphenoPRS)
fit.noE <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRSAFRmeta, data = AFRphenoPRS)

rSex     <-summary(fit.alt1)$r.squared - summary(fit.noSex)$r.squared # due to Sex  0.01086638 (1.09%)
rE       <-summary(fit.alt1)$r.squared - summary(fit.noE)$r.squared # due to IPV  0.0381159 (3.81%)
rSexIPV  <-summary(fit.alt1)$r.squared - summary(fit.noSxE)$r.squared # due to SexByIPV  0.0006619024
rPRS     <-summary(fit.alt1)$r.squared - summary(fit.noG)$r.squared # due to PRS  0.0001802638
rPRSxIPV <-summary(fit.alt1)$r.squared - summary(fit.noGxE)$r.squared # due to PRSxIPV  1.810919e-05


## AMR ----

fit.alt <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRSAMRmeta + ever_ipv + depPRSAMRmeta*ever_ipv, data = AMRphenoPRS)
fit.alt1 <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRSAMRmeta + ever_ipv + biosex*ever_ipv + depPRSAMRmeta*ever_ipv, data = AMRphenoPRS)

summary(fit.alt)
summary(fit.alt1)


fit.alt1 <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRSAMRmeta + ever_ipv + biosex*ever_ipv + depPRSAMRmeta*ever_ipv, data = AMRphenoPRS)
fit.noSex <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + depPRSAMRmeta + ever_ipv + depPRSAMRmeta*ever_ipv, data = AMRphenoPRS)
fit.noSxE <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRSAMRmeta + ever_ipv + depPRSAMRmeta*ever_ipv, data = AMRphenoPRS)
fit.noG <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + ever_ipv + biosex*ever_ipv, data = AMRphenoPRS)
fit.noGxE <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRSAMRmeta + ever_ipv + biosex*ever_ipv, data = AMRphenoPRS)
fit.noE <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRSAMRmeta, data = AMRphenoPRS)

rSex     <-summary(fit.alt1)$r.squared - summary(fit.noSex)$r.squared #due to Sex 
rE       <-summary(fit.alt1)$r.squared - summary(fit.noE)$r.squared #due to IPV 
rSexIPV  <-summary(fit.alt1)$r.squared - summary(fit.noSxE)$r.squared #due to SexByIPV 
rPRS     <-summary(fit.alt1)$r.squared - summary(fit.noG)$r.squared #due to PRS 
rPRSxIPV <-summary(fit.alt1)$r.squared - summary(fit.noGxE)$r.squared #due to PRSxIPV 


## EAS ----

fit.alt <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRSEASmeta + ever_ipv + depPRSEASmeta*ever_ipv, data = EASphenoPRS)
fit.alt1 <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRSEASmeta + ever_ipv + biosex*ever_ipv + depPRSEASmeta*ever_ipv, data = EASphenoPRS)


summary(fit.alt)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             6.102e+00  5.801e-01  10.519  < 2e-16 ***
#   PC1                     8.208e+00  4.039e+00   2.032   0.0425 *  
#   PC2                    -1.749e+00  4.051e+00  -0.432   0.6661    
#   PC3                    -1.013e+00  4.200e+00  -0.241   0.8095    
#   PC4                    -3.117e+00  4.123e+00  -0.756   0.4500    
#   PC5                     1.942e+00  3.964e+00   0.490   0.6244    
#   biosex                  1.597e+00  2.850e-01   5.603 2.99e-08 ***
#   depPRSEASmeta           1.558e+05  8.639e+05   0.180   0.8569    
#   ever_ipv                1.377e+00  6.279e-01   2.193   0.0286 *  
#   depPRSEASmeta:ever_ipv  1.688e+06  1.519e+06   1.112   0.2667    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3.684 on 724 degrees of freedom
# Multiple R-squared:  0.1074,    Adjusted R-squared:  0.09632 
# F-statistic: 9.681 on 9 and 724 DF,  p-value: 4.762e-14

summary(fit.alt1)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             6.145e+00  6.506e-01   9.444  < 2e-16 ***
#   PC1                     8.197e+00  4.043e+00   2.028    0.043 *  
#   PC2                    -1.760e+00  4.054e+00  -0.434    0.664    
#   PC3                    -9.569e-01  4.220e+00  -0.227    0.821    
#   PC4                    -3.128e+00  4.127e+00  -0.758    0.449    
#   PC5                     1.942e+00  3.966e+00   0.490    0.625    
#   biosex                  1.571e+00  3.350e-01   4.691 3.25e-06 ***
#   depPRSEASmeta           1.512e+05  8.650e+05   0.175    0.861    
#   ever_ipv                1.215e+00  1.283e+00   0.946    0.344    
#   biosex:ever_ipv         9.221e-02  6.357e-01   0.145    0.885    
#   depPRSEASmeta:ever_ipv  1.714e+06  1.530e+06   1.120    0.263    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3.687 on 723 degrees of freedom
# Multiple R-squared:  0.1074,    Adjusted R-squared:  0.0951 
# F-statistic: 8.703 on 10 and 723 DF,  p-value: 1.528e-13

fit.alt1 <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRSEASmeta + ever_ipv + biosex*ever_ipv + depPRSEASmeta*ever_ipv, data = EASphenoPRS)
fit.noSex <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + depPRSEASmeta + ever_ipv + depPRSEASmeta*ever_ipv, data = EASphenoPRS)
fit.noSxE <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRSEASmeta + ever_ipv + depPRSEASmeta*ever_ipv, data = EASphenoPRS)
fit.noG <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + ever_ipv + biosex*ever_ipv, data = EASphenoPRS)
fit.noGxE <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRSEASmeta + ever_ipv + biosex*ever_ipv, data = EASphenoPRS)
fit.noE <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRSEASmeta, data = EASphenoPRS)

rSex     <-summary(fit.alt1)$r.squared - summary(fit.noSex)$r.squared # due to Sex 0.03873258 (3.87%)
rE       <-summary(fit.alt1)$r.squared - summary(fit.noE)$r.squared # due to IPV 0.05583847 (5.58%)
rSexIPV  <-summary(fit.alt1)$r.squared - summary(fit.noSxE)$r.squared # due to SexByIPV  2.597601e-05
rPRS     <-summary(fit.alt1)$r.squared - summary(fit.noG)$r.squared # due to PRS  0.002750196 (0.28%)
rPRSxIPV <-summary(fit.alt1)$r.squared - summary(fit.noGxE)$r.squared # due to PRSxIPV  0.001548858 (0.15%)


## SAS ----

fit.alt <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRSSASmeta + ever_ipv + depPRSSASmeta*ever_ipv, data = SASphenoPRS)
fit.alt1 <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRSSASmeta + ever_ipv + biosex*ever_ipv + depPRSSASmeta*ever_ipv, data = SASphenoPRS)

summary(fit.alt)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             5.901e+00  5.077e-01  11.623  < 2e-16 ***
#   PC1                    -8.957e+00  4.241e+00  -2.112   0.0351 *  
#   PC2                     7.075e+00  4.027e+00   1.757   0.0794 .  
#   PC3                     1.790e+01  4.271e+00   4.190 3.21e-05 ***
#   PC4                    -5.704e+00  4.137e+00  -1.379   0.1684    
#   PC5                    -1.480e+00  4.243e+00  -0.349   0.7274    
#   biosex                  1.491e+00  2.996e-01   4.975 8.50e-07 ***
#   depPRSSASmeta          -3.389e+05  3.077e+05  -1.102   0.2711    
#   ever_ipv                1.498e+00  3.288e-01   4.557 6.27e-06 ***
#   depPRSSASmeta:ever_ipv  3.111e+05  5.146e+05   0.605   0.5457    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3.673 on 612 degrees of freedom
# Multiple R-squared:  0.127,     Adjusted R-squared:  0.1141 
# F-statistic: 9.888 on 9 and 612 DF,  p-value: 3.164e-14

summary(fit.alt1)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             6.105e+00  6.011e-01  10.157  < 2e-16 ***
#   PC1                    -8.821e+00  4.248e+00  -2.077 0.038262 *  
#   PC2                     7.051e+00  4.029e+00   1.750 0.080608 .  
#   PC3                     1.782e+01  4.275e+00   4.168 3.51e-05 ***
#   PC4                    -5.820e+00  4.143e+00  -1.405 0.160577    
#   PC5                    -1.338e+00  4.251e+00  -0.315 0.753029    
#   biosex                  1.362e+00  3.621e-01   3.760 0.000186 ***
#   depPRSSASmeta          -3.431e+05  3.079e+05  -1.114 0.265591    
#   ever_ipv                8.514e-01  1.070e+00   0.796 0.426404    
#   biosex:ever_ipv         4.099e-01  6.451e-01   0.635 0.525466    
#   depPRSSASmeta:ever_ipv  3.083e+05  5.149e+05   0.599 0.549595    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3.674 on 611 degrees of freedom
# Multiple R-squared:  0.1275,    Adjusted R-squared:  0.1133 
# F-statistic: 8.931 on 10 and 611 DF,  p-value: 8.614e-14


fit.alt1 <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRSSASmeta + ever_ipv + biosex*ever_ipv + depPRSSASmeta*ever_ipv, data = SASphenoPRS)
fit.noSex <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + depPRSSASmeta + ever_ipv + depPRSSASmeta*ever_ipv, data = SASphenoPRS)
fit.noSxE <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRSSASmeta + ever_ipv + depPRSSASmeta*ever_ipv, data = SASphenoPRS)
fit.noG <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + ever_ipv + biosex*ever_ipv, data = SASphenoPRS)
fit.noGxE <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRSSASmeta + ever_ipv + biosex*ever_ipv, data = SASphenoPRS)
fit.noE <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRSSASmeta, data = SASphenoPRS)

rSex     <-summary(fit.alt1)$r.squared - summary(fit.noSex)$r.squared # due to Sex 0.035882 (3.59%)
rE       <-summary(fit.alt1)$r.squared - summary(fit.noE)$r.squared # due to IPV  0.03389128 (3.39%)
rSexIPV  <-summary(fit.alt1)$r.squared - summary(fit.noSxE)$r.squared # due to SexByIPV 0.0005763366
rPRS     <-summary(fit.alt1)$r.squared - summary(fit.noG)$r.squared # due to PRS  0.00177977 (0.18%)
rPRSxIPV <-summary(fit.alt1)$r.squared - summary(fit.noGxE)$r.squared # due to PRSxIPV  0.0005118311


## complete model ----

phenoPRS2 <- phenoPRS
phenoPRS2$depPRS <- apply(phenoPRS2[, c(18:21)], 1, function(x) as.numeric(paste(x[!is.na(x)])))
phenoPRS2 <- phenoPRS2 %>%
  select(-c(ends_with('meta')))

fit.alt <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRS + ever_ipv + depPRS*ever_ipv, data = phenoPRS2)
fit.alt1 <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRS + ever_ipv + biosex*ever_ipv + depPRS*ever_ipv, data = phenoPRS2)

summary(fit.alt)
summary(fit.alt1)


fit.alt1 <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRS + ever_ipv + biosex*ever_ipv + depPRS*ever_ipv, data = phenoPRS2)
fit.noSex <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + depPRS + ever_ipv + depPRS*ever_ipv, data = phenoPRS2)
fit.noSxE <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRS + ever_ipv + depPRS*ever_ipv, data = phenoPRS2)
fit.noG <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + ever_ipv + biosex*ever_ipv, data = phenoPRS2)
fit.noGxE <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRS + ever_ipv + biosex*ever_ipv, data = phenoPRS2)
fit.noE <- lm(depScore~PC1 + PC2 + PC3 + PC4 + PC5 + biosex + depPRS, data = phenoPRS)

rSex     <-summary(fit.alt1)$r.squared - summary(fit.noSex)$r.squared # due to Sex 
rE       <-summary(fit.alt1)$r.squared - summary(fit.noE)$r.squared # due to IPV 
rSexIPV  <-summary(fit.alt1)$r.squared - summary(fit.noSxE)$r.squared # due to SexByIPV 
rPRS     <-summary(fit.alt1)$r.squared - summary(fit.noG)$r.squared # due to PRS 
rPRSxIPV <-summary(fit.alt1)$r.squared - summary(fit.noGxE)$r.squared # due to PRSxIPV 

q()