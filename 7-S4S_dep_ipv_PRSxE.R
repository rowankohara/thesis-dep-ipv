# Rowan K. O'Hara
# 07.09.2023
# Thesis Project - PRSxE script
# Run depression-PRS association models, complete PRS x IPV interaction models,
# & meta-analyze to estimate cross-population effects
#
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
library(car)
library(meta)

# read in phenotypic data
pheno <- fread("thesis_s4s_clean.csv")
colnames(pheno)
# [1] "unique_id" "PC1"       "PC2"       "PC3"       "PC4"       "PC5"      
# [7] "PC6"       "PC7"       "PC8"       "PC9"       "PC10"      "ancestry" 
# [13] "biosex"    "depScore"  "ever_ipv"  "phys_ipv"  "sex_ipv"

# make factor variables as factors
pheno <- pheno %>%
  mutate(biosex = as.factor(biosex),
         ever_ipv = as.factor(ever_ipv),
         phys_ipv = as.factor(phys_ipv),
         sex_ipv = as.factor(sex_ipv))

# read in PRS scores
PRS <- fread("output/All_dep_S4S_PRS.summary")
PRS <- rename(PRS, 
              "unique_id" = "IID",
              "depPRS" = "PRS_SCORE")
colnames(PRS)
# [1] "unique_id" "depPRS"

# combine data frames
phenoPRS <- inner_join(pheno, PRS, by = "unique_id")
colnames(phenoPRS)
# [1] "unique_id" "PC1"       "PC2"       "PC3"       "PC4"       "PC5"      
# [7] "PC6"       "PC7"       "PC8"       "PC9"       "PC10"      "ancestry" 
# [13] "biosex"    "depScore"  "ever_ipv"  "phys_ipv"  "sex_ipv"   "depPRS"

# subset by ancestry & center and scale PRS within ancestry as well
AFRphenoPRS <- phenoPRS %>%
  subset(ancestry=="AFR") %>%
  mutate(across(c(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, depPRS), scale, .names = "{.col}.centered"))
dim(AFRphenoPRS)
# [1] 1525   29

AMRphenoPRS <- phenoPRS %>%
  subset(ancestry=="AMR") %>%
  mutate(across(c(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, depPRS), scale, .names = "{.col}.centered"))
dim(AMRphenoPRS)
# [1] 921  29

EASphenoPRS <- phenoPRS %>%
  subset(ancestry=="EAS") %>%
  mutate(across(c(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, depPRS), scale, .names = "{.col}.centered"))
dim(EASphenoPRS)
# [1] 734  29

EURphenoPRS <- phenoPRS %>%
  subset(ancestry=="EUR") %>%
  mutate(across(c(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, depPRS), scale, .names = "{.col}.centered"))
dim(EURphenoPRS)
# [1] 3700   29

SASphenoPRS <- phenoPRS %>%
  subset(ancestry=="SAS") %>%
  mutate(across(c(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, depPRS), scale, .names = "{.col}.centered"))
dim(SASphenoPRS)
# [1] 622  29


# LINEAR REGRESSION ############################################################

## EUR ----
fit.alt1 <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipv + biosex*ever_ipv + depPRS.centered*ever_ipv, data = EURphenoPRS)

car::vif(fit.alt1, type = 'predictor')
#                       GVIF Df GVIF^(1/(2*Df))          Interacts With
#   PC1.centered    1.067479  1        1.033189                    --  
#   PC2.centered    1.018356  1        1.009136                    --  
#   PC3.centered    1.080563  1        1.039501                    --  
#   PC4.centered    1.003733  1        1.001865                    --  
#   PC5.centered    1.054042  1        1.026666                    --  
#   PC6.centered    1.138225  1        1.066876                    --  
#   PC7.centered    1.108622  1        1.052911                    --  
#   PC8.centered    1.284960  1        1.133561                    --  
#   PC9.centered    1.004110  1        1.002053                    --  
#   PC10.centered   1.006508  1        1.003249                    --  
#   biosex          1.023133  3        1.003819                ever_ipv
#   depPRS.centered 2.888449  3        1.193376                ever_ipv
#   ever_ipv        1.022775  5        1.002254 biosex, depPRS.centered

summary(fit.alt1) # 8% of variance explained, sex + ipv + PRS are sig
# Coefficients:
#                         Estimate    Std. Error   t value   Pr(>|t|)
# (Intercept)                7.89926    0.11745  67.255 0.00e+00 ***
# PC1.centered              -0.36129    0.06170  -5.856 5.16e-09 ***
# PC2.centered              -0.16700    0.06026  -2.771 0.005611 ** 
# PC3.centered               0.10061    0.06207   1.621 0.105140    
# PC4.centered              -0.18337    0.05983  -3.065 0.002192 ** 
# PC5.centered              -0.04665    0.06131  -0.761 0.446783    
# PC6.centered              -0.05214    0.06371  -0.818 0.413140    
# PC7.centered               0.09344    0.06287   1.486 0.137329    
# PC8.centered               0.08945    0.06769   1.322 0.186412    
# PC9.centered               0.04079    0.05984   0.682 0.495503    
# PC10.centered             -0.01970    0.05991  -0.329 0.742286    
# biosex2                    1.05003    0.15650   6.710 2.25e-11 ***
# depPRS.centered            0.27733    0.07813   3.550 0.000391 ***
# ever_ipv1                  1.35244    0.20580   6.572 5.67e-11 ***
# biosex2:ever_ipv1          0.09426    0.25660   0.367 0.713380    
# depPRS.centered:ever_ipv1 -0.14823    0.12183  -1.217 0.223814    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3.632 on 3684 degrees of freedom
# Multiple R-squared:  0.08557,   Adjusted R-squared:  0.08184 
# F-statistic: 22.98 on 15 and 3684 DF,  p-value: < 2.2e-16


fit.noSex <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + depPRS.centered + ever_ipv + depPRS.centered*ever_ipv, data = EURphenoPRS)
fit.noSxE <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipv + depPRS.centered*ever_ipv, data = EURphenoPRS)
fit.noG <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + ever_ipv + biosex*ever_ipv, data = EURphenoPRS)
fit.noGxE <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipv + biosex*ever_ipv, data = EURphenoPRS)
fit.noE <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered, data = EURphenoPRS)

rSex     <-summary(fit.alt1)$r.squared - summary(fit.noSex)$r.squared # due to Sex - 0.01893371 (1.89%)
rE       <-summary(fit.alt1)$r.squared - summary(fit.noE)$r.squared # due to IPV - 0.03292689 (3.29%)
rSexIPV  <-summary(fit.alt1)$r.squared - summary(fit.noSxE)$r.squared # due to SexByIPV - 3.349582e-05
rPRS     <-summary(fit.alt1)$r.squared - summary(fit.noG)$r.squared # due to PRS - 0.003600713 (0.4%)
rPRSxIPV <-summary(fit.alt1)$r.squared - summary(fit.noGxE)$r.squared # due to PRSxIPV - 0.0003674234




## AFR ----

fit.alt1 <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipv + biosex*ever_ipv + depPRS.centered*ever_ipv, data = AFRphenoPRS)

car::vif(fit.alt1, type = 'predictor')
#                       GVIF Df GVIF^(1/(2*Df))          Interacts With
#   PC1.centered    1.005875  1        1.002933                    --  
#   PC2.centered    1.004619  1        1.002307                    --  
#   PC3.centered    1.006885  1        1.003437                    --  
#   PC4.centered    1.002911  1        1.001454                    --  
#   PC5.centered    1.002834  1        1.001416                    --  
#   PC6.centered    1.004558  1        1.002276                    --  
#   PC7.centered    1.002022  1        1.001011                    --  
#   PC8.centered    1.005181  1        1.002587                    --  
#   PC9.centered    1.003952  1        1.001974                    --  
#   PC10.centered   1.003549  1        1.001773                    --  
#   biosex          1.023666  3        1.003906                ever_ipv
#   depPRS.centered 4.041200  3        1.262075                ever_ipv
#   ever_ipv        1.033032  5        1.003255 biosex, depPRS.centered

summary(fit.alt1) # 8% of variance explained, biosex + ipv are sig
# Coefficients:
#                           Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                7.57966    0.22232  34.094  2.18e-189 ***
# PC1.centered               0.66217    0.09225   7.178  1.1e-12 ***
# PC2.centered              -0.06872    0.09219  -0.745  0.45612    
# PC3.centered              -0.18086    0.09229  -1.960  0.05022 .  
# PC4.centered               0.18946    0.09211   2.057  0.03987 *  
# PC5.centered               0.12657    0.09211   1.374  0.16959    
# PC6.centered              -0.13973    0.09219  -1.516  0.12980    
# PC7.centered              -0.03040    0.09207  -0.330  0.74127    
# PC8.centered              -0.03094    0.09222  -0.336  0.73727    
# PC9.centered               0.12951    0.09216   1.405  0.16015    
# PC10.centered              0.02158    0.09214   0.234  0.81486    
# biosex2                    0.71493    0.26197   2.729  0.00643 ** 
# depPRS.centered           -0.02485    0.11738  -0.212  0.83235    
# ever_ipv1                  1.16730    0.37790   3.089  0.00205 ** 
# biosex2:ever_ipv1          0.43697    0.43720   0.999  0.31773    
# depPRS.centered:ever_ipv1 -0.07412    0.19040  -0.389  0.69711    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3.591 on 1509 degrees of freedom
# Multiple R-squared:  0.09098,   Adjusted R-squared:  0.08195 
# F-statistic: 10.07 on 15 and 1509 DF,  p-value: < 2.2e-16

fit.noSex <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + depPRS.centered + ever_ipv + depPRS.centered*ever_ipv, data = AFRphenoPRS)
fit.noSxE <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipv + depPRS.centered*ever_ipv, data = AFRphenoPRS)
fit.noG <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + ever_ipv + biosex*ever_ipv, data = AFRphenoPRS)
fit.noGxE <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipv + biosex*ever_ipv, data = AFRphenoPRS)
fit.noE <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered, data = AFRphenoPRS)

rSex     <-summary(fit.alt1)$r.squared - summary(fit.noSex)$r.squared # due to Sex - 0.01106739 (1.11%)
rE       <-summary(fit.alt1)$r.squared - summary(fit.noE)$r.squared # due to IPV - 0.0378609 (3.79%)
rSexIPV  <-summary(fit.alt1)$r.squared - summary(fit.noSxE)$r.squared # due to SexByIPV - 0.0006017464
rPRS     <-summary(fit.alt1)$r.squared - summary(fit.noG)$r.squared # due to PRS - 0.000289854
rPRSxIPV <-summary(fit.alt1)$r.squared - summary(fit.noGxE)$r.squared # due to PRSxIPV - 9.129761e-05


## AMR ----

fit.alt1 <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipv + biosex*ever_ipv + depPRS.centered*ever_ipv, data = AMRphenoPRS)

car::vif(fit.alt1, type = 'predictor')
#                       GVIF Df GVIF^(1/(2*Df))          Interacts With
#   PC1.centered    1.650267  1        1.284627                    --  
#   PC2.centered    1.034974  1        1.017337                    --  
#   PC3.centered    1.009994  1        1.004985                    --  
#   PC4.centered    2.737991  1        1.654688                    --  
#   PC5.centered    1.298140  1        1.139359                    --  
#   PC6.centered    2.092731  1        1.446627                    --  
#   PC7.centered    1.914006  1        1.383476                    --  
#   PC8.centered    1.100298  1        1.048951                    --  
#   PC9.centered    1.018413  1        1.009164                    --  
#   PC10.centered   1.023570  1        1.011716                    --  
#   biosex          1.056471  3        1.009198                ever_ipv
#   depPRS.centered 3.312123  3        1.220912                ever_ipv
#   ever_ipv        1.051776  5        1.005061 biosex, depPRS.centered

summary(fit.alt1) # 9.8% of variance explained, sex + PRS + ipv are sig
# Coefficients:
#                             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                7.70313    0.26315  29.273 4.54e-133 ***
# PC1.centered               0.01440    0.16164   0.089 0.929036    
# PC2.centered               0.02010    0.12800   0.157 0.875280    
# PC3.centered               0.17599    0.12645   1.392 0.164339    
# PC4.centered               0.27057    0.20820   1.300 0.194068    
# PC5.centered               0.03020    0.14336   0.211 0.833191    
# PC6.centered              -0.04906    0.18202  -0.270 0.787572    
# PC7.centered              -0.03275    0.17407  -0.188 0.850787    
# PC8.centered              -0.13028    0.13198  -0.987 0.323841    
# PC9.centered              -0.19843    0.12698  -1.563 0.118467    
# PC10.centered             -0.08452    0.12730  -0.664 0.506873    
# biosex2                    1.14545    0.33719   3.397 0.000711 ***
# depPRS.centered            0.45430    0.16528   2.749 0.006102 ** 
# ever_ipv1                  2.18716    0.46155   4.739  2.5e-06 ***
# biosex2:ever_ipv1         -0.17503    0.55916  -0.313 0.754341    
# depPRS.centered:ever_ipv1 -0.27520    0.25838  -1.065 0.287118    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3.816 on 905 degrees of freedom
# Multiple R-squared:  0.1129,    Adjusted R-squared:  0.09823 
# F-statistic: 7.681 on 15 and 905 DF,  p-value: 2.272e-16

fit.noSex <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + depPRS.centered + ever_ipv + depPRS.centered*ever_ipv, data = AMRphenoPRS)
fit.noSxE <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipv + depPRS.centered*ever_ipv, data = AMRphenoPRS)
fit.noG <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + ever_ipv + biosex*ever_ipv, data = AMRphenoPRS)
fit.noGxE <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipv + biosex*ever_ipv, data = AMRphenoPRS)
fit.noE <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered, data = AMRphenoPRS)

rSex     <-summary(fit.alt1)$r.squared - summary(fit.noSex)$r.squared # due to Sex - 0.01583512 (1.58%)
rE       <-summary(fit.alt1)$r.squared - summary(fit.noE)$r.squared # due to IPV - 0.06264464 (6.26%)
rSexIPV  <-summary(fit.alt1)$r.squared - summary(fit.noSxE)$r.squared # due to SexByIPV - 9.603702e-05
rPRS     <-summary(fit.alt1)$r.squared - summary(fit.noG)$r.squared # due to PRS - 0.00819108
rPRSxIPV <-summary(fit.alt1)$r.squared - summary(fit.noGxE)$r.squared # due to PRSxIPV - 0.00111194


## EAS ----

fit.alt1 <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipv + biosex*ever_ipv + depPRS.centered*ever_ipv, data = EASphenoPRS)

car::vif(fit.alt1, type = 'predictor')
#                       GVIF Df GVIF^(1/(2*Df))          Interacts With
#   PC1.centered    1.028382  1        1.014091                    --  
#   PC2.centered    1.036144  1        1.017912                    --  
#   PC3.centered    1.012710  1        1.006335                    --  
#   PC4.centered    1.002751  1        1.001375                    --  
#   PC5.centered    1.008198  1        1.004091                    --  
#   PC6.centered    1.014123  1        1.007037                    --  
#   PC7.centered    1.110003  1        1.053567                    --  
#   PC8.centered    1.012492  1        1.006227                    --  
#   PC9.centered    1.068021  1        1.033451                    --  
#   PC10.centered   1.013633  1        1.006793                    --  
#   biosex          1.087284  3        1.014045                ever_ipv
#   depPRS.centered 3.187488  3        1.213132                ever_ipv
#   ever_ipv        1.097252  5        1.009324 biosex, depPRS.centered

summary(fit.alt1) # 9.2% of varaince explained, sex + ipv are sig
# Coefficients:
#                           Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                7.777825   0.262577  29.621  1.4e-126 ***
# PC1.centered               0.281302   0.138287   2.034   0.0423 *  
# PC2.centered              -0.044880   0.138808  -0.323   0.7465    
# PC3.centered              -0.035352   0.137229  -0.258   0.7968    
# PC4.centered              -0.099749   0.136553  -0.730   0.4653    
# PC5.centered               0.072941   0.136923   0.533   0.5944    
# PC6.centered              -0.046006   0.137325  -0.335   0.7377    
# PC7.centered              -0.018495   0.143670  -0.129   0.8976    
# PC8.centered              -0.081284   0.137214  -0.592   0.5538    
# PC9.centered              -0.207588   0.140927  -1.473   0.1412    
# PC10.centered             -0.009306   0.137291  -0.068   0.9460    
# biosex2                    1.569346   0.336627   4.662 3.73e-06 ***
# depPRS.centered            0.057434   0.167569   0.343   0.7319    
# ever_ipv1                  1.950299   0.521714   3.738   0.0002 ***
# biosex2:ever_ipv1          0.059701   0.639661   0.093   0.9257    
# depPRS.centered:ever_ipv1  0.281829   0.301616   0.934   0.3504    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3.692 on 718 degrees of freedom
# Multiple R-squared:  0.111,     Adjusted R-squared:  0.09243 
# F-statistic: 5.977 on 15 and 718 DF,  p-value: 7.603e-12

fit.noSex <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + depPRS.centered + ever_ipv + depPRS.centered*ever_ipv, data = EASphenoPRS)
fit.noSxE <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipv + depPRS.centered*ever_ipv, data = EASphenoPRS)
fit.noG <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + ever_ipv + biosex*ever_ipv, data = EASphenoPRS)
fit.noGxE <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipv + biosex*ever_ipv, data = EASphenoPRS)
fit.noE <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered, data = EASphenoPRS)

rSex     <-summary(fit.alt1)$r.squared - summary(fit.noSex)$r.squared # due to Sex - 0.03800151 (3.80%)
rE       <-summary(fit.alt1)$r.squared - summary(fit.noE)$r.squared # due to IPV - 0.05603426 (5.60%)
rSexIPV  <-summary(fit.alt1)$r.squared - summary(fit.noSxE)$r.squared # due to SexByIPV - 1.078562e-05
rPRS     <-summary(fit.alt1)$r.squared - summary(fit.noG)$r.squared # due to PRS - 0.002453416
rPRSxIPV <-summary(fit.alt1)$r.squared - summary(fit.noGxE)$r.squared # due to PRSxIPV - 0.00108103


## SAS ----

fit.alt1 <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipv + biosex*ever_ipv + depPRS.centered*ever_ipv, data = SASphenoPRS)

car::vif(fit.alt1, type = 'predictor')
#                       GVIF Df GVIF^(1/(2*Df))          Interacts With
#   PC1.centered    1.030075  1        1.014926                    --  
#   PC2.centered    1.012834  1        1.006397                    --  
#   PC3.centered    1.036517  1        1.018095                    --  
#   PC4.centered    1.023994  1        1.011926                    --  
#   PC5.centered    1.027658  1        1.013735                    --  
#   PC6.centered    1.007937  1        1.003961                    --  
#   PC7.centered    1.091602  1        1.044798                    --  
#   PC8.centered    1.031497  1        1.015627                    --  
#   PC9.centered    1.074890  1        1.036769                    --  
#   PC10.centered   1.015932  1        1.007934                    --  
#   biosex          1.074611  3        1.012065                ever_ipv
#   depPRS.centered 2.588256  3        1.171749                ever_ipv
#   ever_ipv        1.122868  5        1.011656 biosex, depPRS.centered

summary(fit.alt1) # 12% of variance explained, only sex is sig
# Coefficients:
#                           Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                7.45332    0.27552  27.052 2.7e-106 ***
# PC1.centered              -0.37074    0.14903  -2.488 0.013125 *  
# PC2.centered               0.27774    0.14778   1.879 0.060658 .  
# PC3.centered               0.65776    0.14949   4.400 1.28e-05 ***
# PC4.centered              -0.19694    0.14859  -1.325 0.185528    
# PC5.centered              -0.01151    0.14885  -0.077 0.938405    
# PC6.centered              -0.05546    0.14742  -0.376 0.706895    
# PC7.centered               0.32007    0.15341   2.086 0.037370 *  
# PC8.centered              -0.04882    0.14913  -0.327 0.743482    
# PC9.centered              -0.19472    0.15224  -1.279 0.201367    
# PC10.centered              0.23194    0.14800   1.567 0.117606    
# biosex2                    1.40880    0.36223   3.889 0.000112 ***
# depPRS.centered            0.27915    0.17730   1.574 0.115899    
# ever_ipv1                  1.34615    0.49155   2.739 0.006352 ** 
# biosex2:ever_ipv1          0.35033    0.64617   0.542 0.587910    
# depPRS.centered:ever_ipv1 -0.10764    0.32988  -0.326 0.744320    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3.659 on 606 degrees of freedom
# Multiple R-squared:  0.1418,    Adjusted R-squared:  0.1206 
# F-statistic: 6.676 on 15 and 606 DF,  p-value: 2.055e-13

fit.noSex <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + depPRS.centered + ever_ipv + depPRS.centered*ever_ipv, data = SASphenoPRS)
fit.noSxE <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipv + depPRS.centered*ever_ipv, data = SASphenoPRS)
fit.noG <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + ever_ipv + biosex*ever_ipv, data = SASphenoPRS)
fit.noGxE <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipv + biosex*ever_ipv, data = SASphenoPRS)
fit.noE <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered, data = SASphenoPRS)

rSex     <-summary(fit.alt1)$r.squared - summary(fit.noSex)$r.squared # due to Sex - 0.03669247 (3.67%)
rE       <-summary(fit.alt1)$r.squared - summary(fit.noE)$r.squared # due to IPV - 0.03348438 (3.35%)
rSexIPV  <-summary(fit.alt1)$r.squared - summary(fit.noSxE)$r.squared # due to SexByIPV - 0.000416257
rPRS     <-summary(fit.alt1)$r.squared - summary(fit.noG)$r.squared # due to PRS - 0.004065141
rPRSxIPV <-summary(fit.alt1)$r.squared - summary(fit.noGxE)$r.squared # due to PRSxIPV - 0.0001507681


## complete model ----

# wrangle data to get the n, mean, and standard deviation of the depression symptom scores
# separated by those who had experienced IPV and no IPV exposure
afr <- AFRphenoPRS %>%
  group_by(ever_ipv) %>%
  summarise(n = n(), m = mean(depScore), sd = sd(depScore)) %>%
  pivot_wider(names_from = ever_ipv, values_from = c(n, m, sd)) %>%
  add_column(ancestry = "AFR", .before = 1) %>%
  select(c(ancestry, n_1, m_1, sd_1, n_0, m_0, sd_0))

amr <- AMRphenoPRS %>%
  group_by(ever_ipv) %>%
  summarise(n = n(), m = mean(depScore), sd = sd(depScore)) %>%
  pivot_wider(names_from = ever_ipv, values_from = c(n, m, sd)) %>%
  add_column(ancestry = "AMR", .before = 1) %>%
  select(c(ancestry, n_1, m_1, sd_1, n_0, m_0, sd_0))

eas <- EASphenoPRS %>%
  group_by(ever_ipv) %>%
  summarise(n = n(), m = mean(depScore), sd = sd(depScore)) %>%
  pivot_wider(names_from = ever_ipv, values_from = c(n, m, sd)) %>%
  add_column(ancestry = "EAS", .before = 1) %>%
  select(c(ancestry, n_1, m_1, sd_1, n_0, m_0, sd_0))

eur <- EURphenoPRS %>%
  group_by(ever_ipv) %>%
  summarise(n = n(), m = mean(depScore), sd = sd(depScore)) %>%
  pivot_wider(names_from = ever_ipv, values_from = c(n, m, sd)) %>%
  add_column(ancestry = "EUR", .before = 1) %>%
  select(c(ancestry, n_1, m_1, sd_1, n_0, m_0, sd_0))

sas <- SASphenoPRS %>%
  group_by(ever_ipv) %>%
  summarise(n = n(), m = mean(depScore), sd = sd(depScore)) %>%
  pivot_wider(names_from = ever_ipv, values_from = c(n, m, sd)) %>%
  add_column(ancestry = "SAS", .before = 1) %>%
  select(c(ancestry, n_1, m_1, sd_1, n_0, m_0, sd_0))

df.meta <- rbind(afr, amr, eas, eur, sas)

# use the metacont function for the continuous outcome for meta-analysis
# and the forest function to create a visual table/plot
m.cont <- metacont(n.e = n_1,
         mean.e = m_1,
         sd.e = sd_1,
         n.c = n_0,
         mean.c = m_0,
         sd.c = sd_0,
         studlab = ancestry,
         data = df.meta)
forest(m.cont)
summary(m.cont)
# MD           95%-CI %W(common) %W(random)
# AFR 1.5288 [1.1424; 1.9152]       21.1       23.8
# AMR 2.2701 [1.7436; 2.7966]       11.4       17.0
# EAS 2.0792 [1.4672; 2.6912]        8.4       13.9
# EUR 1.6228 [1.3763; 1.8694]       51.9       32.9
# SAS 1.7745 [1.1124; 2.4365]        7.2       12.5
# 
# Number of studies: k = 5
# Number of observations: o = 7502
# 
# MD           95%-CI     z  p-value
# Common effect model  1.7259 [1.5484; 1.9035] 19.05 < 0.0001
# Random effects model 1.7926 [1.5175; 2.0678] 12.77 < 0.0001
# 
# Quantifying heterogeneity:
#   tau^2 = 0.0441 [0.0000; 0.7451]; tau = 0.2099 [0.0000; 0.8632]
# I^2 = 43.5% [0.0%; 79.2%]; H = 1.33 [1.00; 2.19]
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 7.08    4  0.1319
# 
# Details on meta-analytical method:
#   - Inverse variance method
# - Restricted maximum-likelihood estimator for tau^2
# - Q-Profile method for confidence interval of tau^2 and tau

q()
