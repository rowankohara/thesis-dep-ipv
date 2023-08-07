# Rowan K. O'Hara
# 07.09.2023
# Thesis Project - PRSxE script
# Run depression-PRS association models, complete PRS x ipt interaction models,
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
library(sensemakr)
library(patchwork)

# read in phenotypic data
pheno <- fread("thesis_s4s_clean.csv")
colnames(pheno)
# [1] "unique_id" "PC1"       "PC2"       "PC3"       "PC4"       "PC5"      
# [7] "PC6"       "PC7"       "PC8"       "PC9"       "PC10"      "ancestry" 
# [13] "biosex"    "depScore"  "ever_ipv"  "phys_ipt"  "sex_ipt"

# make factor variables as factors
pheno <- pheno %>%
  mutate(biosex = as.factor(biosex),
         ever_ipv = as.factor(ever_ipv),
         phys_ipv = as.factor(phys_ipv),
         sex_ipv = as.factor(sex_ipv))
pheno <- rename(pheno,
                "ever_ipt" = "ever_ipv",
                "phys_ipt" = "phys_ipv",
                "sex_ipt" = "sex_ipv")

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
# [13] "biosex"    "depScore"  "ever_ipt"  "phys_ipt"  "sex_ipt"   "depPRS"

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
fit.alt1 <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipt + biosex*ever_ipt + depPRS.centered*ever_ipt, data = EURphenoPRS)

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
#   biosex          1.023133  3        1.003819                ever_ipt
#   depPRS.centered 2.888449  3        1.193376                ever_ipt
#   ever_ipt        1.022775  5        1.002254 biosex, depPRS.centered

summary(fit.alt1) # 8% of variance explained, sex + ipt + PRS are sig
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
# ever_ipt1                  1.35244    0.20580   6.572 5.67e-11 ***
# biosex2:ever_ipt1          0.09426    0.25660   0.367 0.713380    
# depPRS.centered:ever_ipt1 -0.14823    0.12183  -1.217 0.223814    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3.632 on 3684 degrees of freedom
# Multiple R-squared:  0.08557,   Adjusted R-squared:  0.08184 
# F-statistic: 22.98 on 15 and 3684 DF,  p-value: < 2.2e-16


fit.noSex <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + depPRS.centered + ever_ipt + depPRS.centered*ever_ipt, data = EURphenoPRS)
fit.noSxE <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipt + depPRS.centered*ever_ipt, data = EURphenoPRS)
fit.noG <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + ever_ipt + biosex*ever_ipt, data = EURphenoPRS)
fit.noGxE <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipt + biosex*ever_ipt, data = EURphenoPRS)
fit.noE <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered, data = EURphenoPRS)

rSex     <-summary(fit.alt1)$r.squared - summary(fit.noSex)$r.squared # due to Sex - 0.01893371 (1.89%)
rE       <-summary(fit.alt1)$r.squared - summary(fit.noE)$r.squared # due to ipt - 0.03292689 (3.29%)
rSexipt  <-summary(fit.alt1)$r.squared - summary(fit.noSxE)$r.squared # due to SexByipt - 3.349582e-05
rPRS     <-summary(fit.alt1)$r.squared - summary(fit.noG)$r.squared # due to PRS - 0.003600713 (0.4%)
rPRSxipt <-summary(fit.alt1)$r.squared - summary(fit.noGxE)$r.squared # due to PRSxipt - 0.0003674234

fit.nul <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipt + biosex * ever_ipt, data = EURphenoPRS)
summary(fit.nul)
# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        7.89861    0.11746  67.246  < 2e-16 ***
# PC1.centered      -0.36282    0.06169  -5.882 4.42e-09 ***
# PC2.centered      -0.16760    0.06026  -2.781 0.005443 ** 
# PC3.centered       0.10031    0.06208   1.616 0.106190    
# PC4.centered      -0.18611    0.05979  -3.113 0.001867 ** 
# PC5.centered      -0.04517    0.06130  -0.737 0.461206    
# PC6.centered      -0.05231    0.06371  -0.821 0.411659    
# PC7.centered       0.09472    0.06287   1.507 0.131973    
# PC8.centered       0.08920    0.06769   1.318 0.187659    
# PC9.centered       0.04167    0.05984   0.696 0.486203    
# PC10.centered     -0.01832    0.05990  -0.306 0.759744    
# biosex2            1.04586    0.15647   6.684 2.67e-11 ***
# depPRS.centered    0.21637    0.05995   3.609 0.000312 ***
# ever_ipt1          1.35099    0.20581   6.564 5.96e-11 ***
# biosex2:ever_ipt1  0.09198    0.25661   0.358 0.720036    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3.632 on 3685 degrees of freedom
# Multiple R-squared:  0.0852,	Adjusted R-squared:  0.08172 
# F-statistic: 24.51 on 14 and 3685 DF,  p-value: < 2.2e-16

partial_r2(fit.nul)
# (Intercept)      PC1.centered      PC2.centered      PC3.centered      PC4.centered 
# 5.509934e-01      9.300553e-03      2.094669e-03      7.081293e-04      2.622756e-03 
# PC5.centered      PC6.centered      PC7.centered      PC8.centered      PC9.centered 
# 1.473539e-04      1.829126e-04      6.156720e-04      4.710244e-04      1.316014e-04 
# PC10.centered           biosex2   depPRS.centered         ever_ipt1 biosex2:ever_ipt1 
# 2.538279e-05      1.197873e-02      3.521975e-03      1.155826e-02      3.486379e-05 

eurBeta <- summary(fit.alt1)$coefficient['depPRS.centered',1]
eurCI_lower <- confint(fit.alt1, 'depPRS.centered', level=0.95)[1,1]
eurCI_upper <- confint(fit.alt1, 'depPRS.centered', level=0.95)[1,2]

eurBeta_ipt <- summary(fit.alt1)$coefficient['ever_ipt1',1]
eurCI_lower_ipt <- confint(fit.alt1, 'ever_ipt1', level=0.95)[1,1]
eurCI_upper_ipt <- confint(fit.alt1, 'ever_ipt1', level=0.95)[1,2]

eurBeta_sex <- summary(fit.alt1)$coefficient['biosex2',1]
eurCI_lower_sex <- confint(fit.alt1, 'biosex2', level=0.95)[1,1]
eurCI_upper_sex <- confint(fit.alt1, 'biosex2', level=0.95)[1,2]


## AFR ----

fit.alt1 <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipt + biosex*ever_ipt + depPRS.centered*ever_ipt, data = AFRphenoPRS)

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
#   biosex          1.023666  3        1.003906                ever_ipt
#   depPRS.centered 4.041200  3        1.262075                ever_ipt
#   ever_ipt        1.033032  5        1.003255 biosex, depPRS.centered

summary(fit.alt1) # 8% of variance explained, biosex + ipt are sig
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
# ever_ipt1                  1.16730    0.37790   3.089  0.00205 ** 
# biosex2:ever_ipt1          0.43697    0.43720   0.999  0.31773    
# depPRS.centered:ever_ipt1 -0.07412    0.19040  -0.389  0.69711    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3.591 on 1509 degrees of freedom
# Multiple R-squared:  0.09098,   Adjusted R-squared:  0.08195 
# F-statistic: 10.07 on 15 and 1509 DF,  p-value: < 2.2e-16

fit.noSex <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + depPRS.centered + ever_ipt + depPRS.centered*ever_ipt, data = AFRphenoPRS)
fit.noSxE <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipt + depPRS.centered*ever_ipt, data = AFRphenoPRS)
fit.noG <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + ever_ipt + biosex*ever_ipt, data = AFRphenoPRS)
fit.noGxE <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipt + biosex*ever_ipt, data = AFRphenoPRS)
fit.noE <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered, data = AFRphenoPRS)

rSex     <-summary(fit.alt1)$r.squared - summary(fit.noSex)$r.squared # due to Sex - 0.01106739 (1.11%)
rE       <-summary(fit.alt1)$r.squared - summary(fit.noE)$r.squared # due to ipt - 0.0378609 (3.79%)
rSexipt  <-summary(fit.alt1)$r.squared - summary(fit.noSxE)$r.squared # due to SexByipt - 0.0006017464
rPRS     <-summary(fit.alt1)$r.squared - summary(fit.noG)$r.squared # due to PRS - 0.000289854
rPRSxipt <-summary(fit.alt1)$r.squared - summary(fit.noGxE)$r.squared # due to PRSxipt - 9.129761e-05

fit.nul <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipt + biosex * ever_ipt, data = AFRphenoPRS)
summary(fit.nul)
# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        7.57953    0.22225  34.103  < 2e-16 ***
# PC1.centered       0.66249    0.09222   7.184 1.06e-12 ***
# PC2.centered      -0.07035    0.09207  -0.764   0.4449    
# PC3.centered      -0.17943    0.09219  -1.946   0.0518 .  
# PC4.centered       0.19030    0.09206   2.067   0.0389 *  
# PC5.centered       0.12621    0.09208   1.371   0.1707    
# PC6.centered      -0.13984    0.09216  -1.517   0.1294    
# PC7.centered      -0.03106    0.09203  -0.337   0.7358    
# PC8.centered      -0.03105    0.09219  -0.337   0.7363    
# PC9.centered       0.12914    0.09213   1.402   0.1612    
# PC10.centered      0.02286    0.09206   0.248   0.8039    
# biosex2            0.71362    0.26188   2.725   0.0065 ** 
# depPRS.centered   -0.05304    0.09236  -0.574   0.5659    
# ever_ipt1          1.16391    0.37770   3.082   0.0021 ** 
# biosex2:ever_ipt1  0.43951    0.43703   1.006   0.3147    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3.59 on 1510 degrees of freedom
# Multiple R-squared:  0.09089,	Adjusted R-squared:  0.08246 
# F-statistic: 10.78 on 14 and 1510 DF,  p-value: < 2.2e-16

partial_r2(fit.nul)
# (Intercept)      PC1.centered      PC2.centered      PC3.centered      PC4.centered 
# 4.350954e-01      3.304982e-02      3.865695e-04      2.502204e-03      2.821847e-03 
# PC5.centered      PC6.centered      PC7.centered      PC8.centered      PC9.centered 
# 1.242680e-03      1.522392e-03      7.540719e-05      7.512129e-05      1.299611e-03 
# PC10.centered           biosex2   depPRS.centered         ever_ipt1 biosex2:ever_ipt1 
# 4.082903e-05      4.893582e-03      2.183600e-04      6.249602e-03      6.693225e-04

afrBeta <- summary(fit.alt1)$coefficient['depPRS.centered',1]
afrCI_lower <- confint(fit.alt1, 'depPRS.centered', level=0.95)[1,1]
afrCI_upper <- confint(fit.alt1, 'depPRS.centered', level=0.95)[1,2]

afrBeta_ipt <- summary(fit.alt1)$coefficient['ever_ipt1',1]
afrCI_lower_ipt <- confint(fit.alt1, 'ever_ipt1', level=0.95)[1,1]
afrCI_upper_ipt <- confint(fit.alt1, 'ever_ipt1', level=0.95)[1,2]

afrBeta_sex <- summary(fit.alt1)$coefficient['biosex2',1]
afrCI_lower_sex <- confint(fit.alt1, 'biosex2', level=0.95)[1,1]
afrCI_upper_sex <- confint(fit.alt1, 'biosex2', level=0.95)[1,2]


## AMR ----

fit.alt1 <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipt + biosex*ever_ipt + depPRS.centered*ever_ipt, data = AMRphenoPRS)

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
#   biosex          1.056471  3        1.009198                ever_ipt
#   depPRS.centered 3.312123  3        1.220912                ever_ipt
#   ever_ipt        1.051776  5        1.005061 biosex, depPRS.centered

summary(fit.alt1) # 9.8% of variance explained, sex + PRS + ipt are sig
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
# ever_ipt1                  2.18716    0.46155   4.739  2.5e-06 ***
# biosex2:ever_ipt1         -0.17503    0.55916  -0.313 0.754341    
# depPRS.centered:ever_ipt1 -0.27520    0.25838  -1.065 0.287118    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3.816 on 905 degrees of freedom
# Multiple R-squared:  0.1129,    Adjusted R-squared:  0.09823 
# F-statistic: 7.681 on 15 and 905 DF,  p-value: 2.272e-16

fit.noSex <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + depPRS.centered + ever_ipt + depPRS.centered*ever_ipt, data = AMRphenoPRS)
fit.noSxE <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipt + depPRS.centered*ever_ipt, data = AMRphenoPRS)
fit.noG <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + ever_ipt + biosex*ever_ipt, data = AMRphenoPRS)
fit.noGxE <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipt + biosex*ever_ipt, data = AMRphenoPRS)
fit.noE <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered, data = AMRphenoPRS)

rSex     <-summary(fit.alt1)$r.squared - summary(fit.noSex)$r.squared # due to Sex - 0.01583512 (1.58%)
rE       <-summary(fit.alt1)$r.squared - summary(fit.noE)$r.squared # due to ipt - 0.06264464 (6.26%)
rSexipt  <-summary(fit.alt1)$r.squared - summary(fit.noSxE)$r.squared # due to SexByipt - 9.603702e-05
rPRS     <-summary(fit.alt1)$r.squared - summary(fit.noG)$r.squared # due to PRS - 0.00819108
rPRSxipt <-summary(fit.alt1)$r.squared - summary(fit.noGxE)$r.squared # due to PRSxipt - 0.00111194

fit.nul <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipt + biosex * ever_ipt, data = AMRphenoPRS)
summary(fit.nul)
# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        7.69439    0.26304  29.252  < 2e-16 ***
# PC1.centered       0.01499    0.16165   0.093 0.926117    
# PC2.centered       0.02244    0.12799   0.175 0.860848    
# PC3.centered       0.17684    0.12646   1.398 0.162338    
# PC4.centered       0.26322    0.20810   1.265 0.206237    
# PC5.centered       0.03060    0.14337   0.213 0.831053    
# PC6.centered      -0.05414    0.18197  -0.298 0.766129    
# PC7.centered      -0.02810    0.17403  -0.161 0.871766    
# PC8.centered      -0.13404    0.13194  -1.016 0.309958    
# PC9.centered      -0.19897    0.12698  -1.567 0.117497    
# PC10.centered     -0.08033    0.12725  -0.631 0.528029    
# biosex2            1.14365    0.33721   3.391 0.000725 ***
# depPRS.centered    0.34195    0.12725   2.687 0.007337 ** 
# ever_ipt1          2.19305    0.46155   4.751 2.35e-06 ***
# biosex2:ever_ipt1 -0.19896    0.55875  -0.356 0.721859    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3.817 on 906 degrees of freedom
# Multiple R-squared:  0.1118,	Adjusted R-squared:  0.0981 
# F-statistic: 8.148 on 14 and 906 DF,  p-value: < 2.2e-16

partial_r2(fit.nul)
# (Intercept)      PC1.centered      PC2.centered      PC3.centered      PC4.centered 
# 4.857120e-01      9.496506e-06      3.393467e-05      2.153739e-03      1.762826e-03 
# PC5.centered      PC6.centered      PC7.centered      PC8.centered      PC9.centered 
# 5.026807e-05      9.770016e-05      2.877359e-05      1.137790e-03      2.702435e-03 
# PC10.centered           biosex2   depPRS.centered         ever_ipt1 biosex2:ever_ipt1 
# 4.396418e-04      1.253621e-02      7.907406e-03      2.431290e-02      1.399334e-04

amrBeta <- summary(fit.alt1)$coefficient['depPRS.centered',1]
amrCI_lower <- confint(fit.alt1, 'depPRS.centered', level=0.95)[1,1]
amrCI_upper <- confint(fit.alt1, 'depPRS.centered', level=0.95)[1,2]

amrBeta_ipt <- summary(fit.alt1)$coefficient['ever_ipt1',1]
amrCI_lower_ipt <- confint(fit.alt1, 'ever_ipt1', level=0.95)[1,1]
amrCI_upper_ipt <- confint(fit.alt1, 'ever_ipt1', level=0.95)[1,2]

amrBeta_sex <- summary(fit.alt1)$coefficient['biosex2',1]
amrCI_lower_sex <- confint(fit.alt1, 'biosex2', level=0.95)[1,1]
amrCI_upper_sex <- confint(fit.alt1, 'biosex2', level=0.95)[1,2]


## EAS ----

fit.alt1 <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipt + biosex*ever_ipt + depPRS.centered*ever_ipt, data = EASphenoPRS)

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
#   biosex          1.087284  3        1.014045                ever_ipt
#   depPRS.centered 3.187488  3        1.213132                ever_ipt
#   ever_ipt        1.097252  5        1.009324 biosex, depPRS.centered

summary(fit.alt1) # 9.2% of varaince explained, sex + ipt are sig
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
# ever_ipt1                  1.950299   0.521714   3.738   0.0002 ***
# biosex2:ever_ipt1          0.059701   0.639661   0.093   0.9257    
# depPRS.centered:ever_ipt1  0.281829   0.301616   0.934   0.3504    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3.692 on 718 degrees of freedom
# Multiple R-squared:  0.111,     Adjusted R-squared:  0.09243 
# F-statistic: 5.977 on 15 and 718 DF,  p-value: 7.603e-12

fit.noSex <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + depPRS.centered + ever_ipt + depPRS.centered*ever_ipt, data = EASphenoPRS)
fit.noSxE <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipt + depPRS.centered*ever_ipt, data = EASphenoPRS)
fit.noG <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + ever_ipt + biosex*ever_ipt, data = EASphenoPRS)
fit.noGxE <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipt + biosex*ever_ipt, data = EASphenoPRS)
fit.noE <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered, data = EASphenoPRS)

rSex     <-summary(fit.alt1)$r.squared - summary(fit.noSex)$r.squared # due to Sex - 0.03800151 (3.80%)
rE       <-summary(fit.alt1)$r.squared - summary(fit.noE)$r.squared # due to ipt - 0.05603426 (5.60%)
rSexipt  <-summary(fit.alt1)$r.squared - summary(fit.noSxE)$r.squared # due to SexByipt - 1.078562e-05
rPRS     <-summary(fit.alt1)$r.squared - summary(fit.noG)$r.squared # due to PRS - 0.002453416
rPRSxipt <-summary(fit.alt1)$r.squared - summary(fit.noGxE)$r.squared # due to PRSxipt - 0.00108103

fit.nul <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipt + biosex * ever_ipt, data = EASphenoPRS)
summary(fit.nul)
# Coefficients:
#                   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        7.77574    0.26254  29.617  < 2e-16 ***
# PC1.centered       0.28202    0.13827   2.040   0.0418 *  
# PC2.centered      -0.03447    0.13835  -0.249   0.8033    
# PC3.centered      -0.03934    0.13715  -0.287   0.7743    
# PC4.centered      -0.09753    0.13652  -0.714   0.4752    
# PC5.centered       0.06898    0.13685   0.504   0.6143    
# PC6.centered      -0.04703    0.13731  -0.342   0.7321    
# PC7.centered      -0.03521    0.14254  -0.247   0.8050    
# PC8.centered      -0.08226    0.13720  -0.600   0.5490    
# PC9.centered      -0.21526    0.14067  -1.530   0.1264    
# PC10.centered     -0.01241    0.13724  -0.090   0.9280    
# biosex2            1.58037    0.33639   4.698 3.15e-06 ***
# depPRS.centered    0.14571    0.13839   1.053   0.2927    
# ever_ipt1          2.02248    0.51592   3.920 9.69e-05 ***
# biosex2:ever_ipt1 -0.01593    0.63446  -0.025   0.9800    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3.692 on 719 degrees of freedom
# Multiple R-squared:  0.1099,	Adjusted R-squared:  0.09259 
# F-statistic: 6.343 on 14 and 719 DF,  p-value: 4.256e-12

partial_r2(fit.nul)
# (Intercept)      PC1.centered      PC2.centered      PC3.centered      PC4.centered 
# 5.495437e-01      5.752567e-03      8.633780e-05      1.144187e-04      7.093025e-04 
# PC5.centered      PC6.centered      PC7.centered      PC8.centered      PC9.centered 
# 3.533060e-04      1.631078e-04      8.484672e-05      4.997179e-04      3.246015e-03 
# PC10.centered           biosex2   depPRS.centered         ever_ipt1 biosex2:ever_ipt1 
# 1.137985e-05      2.978311e-02      1.539500e-03      2.092645e-02      8.767752e-07

easBeta <- summary(fit.alt1)$coefficient['depPRS.centered',1]
easCI_lower <- confint(fit.alt1, 'depPRS.centered', level=0.95)[1,1]
easCI_upper <- confint(fit.alt1, 'depPRS.centered', level=0.95)[1,2]

easBeta_ipt <- summary(fit.alt1)$coefficient['ever_ipt1',1]
easCI_lower_ipt <- confint(fit.alt1, 'ever_ipt1', level=0.95)[1,1]
easCI_upper_ipt <- confint(fit.alt1, 'ever_ipt1', level=0.95)[1,2]

easBeta_sex <- summary(fit.alt1)$coefficient['biosex2',1]
easCI_lower_sex <- confint(fit.alt1, 'biosex2', level=0.95)[1,1]
easCI_upper_sex <- confint(fit.alt1, 'biosex2', level=0.95)[1,2]


## SAS ----

fit.alt1 <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipt + biosex*ever_ipt + depPRS.centered*ever_ipt, data = SASphenoPRS)

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
#   biosex          1.074611  3        1.012065                ever_ipt
#   depPRS.centered 2.588256  3        1.171749                ever_ipt
#   ever_ipt        1.122868  5        1.011656 biosex, depPRS.centered

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
# ever_ipt1                  1.34615    0.49155   2.739 0.006352 ** 
# biosex2:ever_ipt1          0.35033    0.64617   0.542 0.587910    
# depPRS.centered:ever_ipt1 -0.10764    0.32988  -0.326 0.744320    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3.659 on 606 degrees of freedom
# Multiple R-squared:  0.1418,    Adjusted R-squared:  0.1206 
# F-statistic: 6.676 on 15 and 606 DF,  p-value: 2.055e-13

fit.noSex <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + depPRS.centered + ever_ipt + depPRS.centered*ever_ipt, data = SASphenoPRS)
fit.noSxE <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipt + depPRS.centered*ever_ipt, data = SASphenoPRS)
fit.noG <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + ever_ipt + biosex*ever_ipt, data = SASphenoPRS)
fit.noGxE <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipt + biosex*ever_ipt, data = SASphenoPRS)
fit.noE <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered, data = SASphenoPRS)

rSex     <-summary(fit.alt1)$r.squared - summary(fit.noSex)$r.squared # due to Sex - 0.03669247 (3.67%)
rE       <-summary(fit.alt1)$r.squared - summary(fit.noE)$r.squared # due to ipt - 0.03348438 (3.35%)
rSexipt  <-summary(fit.alt1)$r.squared - summary(fit.noSxE)$r.squared # due to SexByipt - 0.000416257
rPRS     <-summary(fit.alt1)$r.squared - summary(fit.noG)$r.squared # due to PRS - 0.004065141
rPRSxipt <-summary(fit.alt1)$r.squared - summary(fit.noGxE)$r.squared # due to PRSxipt - 0.0001507681

fit.nul <- lm(depScore~PC1.centered + PC2.centered + PC3.centered + PC4.centered + PC5.centered + PC6.centered + PC7.centered + PC8.centered + PC9.centered + PC10.centered + biosex + depPRS.centered + ever_ipt + biosex * ever_ipt, data = SASphenoPRS)
summary(fit.nul)
# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        7.45463    0.27529  27.080  < 2e-16 ***
# PC1.centered      -0.37098    0.14892  -2.491 0.012998 *  
# PC2.centered       0.27964    0.14755   1.895 0.058540 .  
# PC3.centered       0.65739    0.14938   4.401 1.27e-05 ***
# PC4.centered      -0.19524    0.14839  -1.316 0.188748    
# PC5.centered      -0.01533    0.14828  -0.103 0.917672    
# PC6.centered      -0.05597    0.14730  -0.380 0.704091    
# PC7.centered       0.32148    0.15324   2.098 0.036328 *  
# PC8.centered      -0.05330    0.14839  -0.359 0.719550    
# PC9.centered      -0.18996    0.15142  -1.254 0.210153    
# PC10.centered      0.23257    0.14788   1.573 0.116307    
# biosex2            1.40538    0.36182   3.884 0.000114 ***
# depPRS.centered    0.24783    0.14896   1.664 0.096673 .  
# ever_ipt1          1.34810    0.49115   2.745 0.006235 ** 
# biosex2:ever_ipt1  0.34172    0.64516   0.530 0.596532    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3.656 on 607 degrees of freedom
# Multiple R-squared:  0.1417,	Adjusted R-squared:  0.1219 
# F-statistic: 7.155 on 14 and 607 DF,  p-value: 7.814e-14

partial_r2(fit.nul)
# (Intercept)      PC1.centered      PC2.centered      PC3.centered      PC4.centered 
# 5.471164e-01      1.012060e-02      5.882485e-03      3.091986e-02      2.844030e-03 
# PC5.centered      PC6.centered      PC7.centered      PC8.centered      PC9.centered 
# 1.761695e-05      2.378128e-04      7.198455e-03      2.125486e-04      2.585883e-03 
# PC10.centered           biosex2   depPRS.centered         ever_ipt1 biosex2:ever_ipt1 
# 4.058280e-03      2.425277e-02      4.539680e-03      1.225925e-02      4.619850e-04

sasBeta <- summary(fit.alt1)$coefficient['depPRS.centered',1]
sasCI_lower <- confint(fit.alt1, 'depPRS.centered', level=0.95)[1,1]
sasCI_upper <- confint(fit.alt1, 'depPRS.centered', level=0.95)[1,2]

sasBeta_ipt <- summary(fit.alt1)$coefficient['ever_ipt1',1]
sasCI_lower_ipt <- confint(fit.alt1, 'ever_ipt1', level=0.95)[1,1]
sasCI_upper_ipt <- confint(fit.alt1, 'ever_ipt1', level=0.95)[1,2]

sasBeta_sex <- summary(fit.alt1)$coefficient['biosex2',1]
sasCI_lower_sex <- confint(fit.alt1, 'biosex2', level=0.95)[1,1]
sasCI_upper_sex <- confint(fit.alt1, 'biosex2', level=0.95)[1,2]


## meta-analysis ----

# meta-analyzed the estimates (beta) and standard error of each ancestry using METASOFT
meta.dep <- fread("../../data/S4S_dep_ipt_meta.csv")

metaBeta <- 0.20910100
metaCI_lower <- 0.20910100 - 1.96 * 0.0541985
metaCI_upper <- 0.20910100 + 1.96 * 0.0541985

metaBeta_ipt <- 1.32042000
metaCI_lower_ipt <- 1.32042000 - 1.96 * 0.3236910
metaCI_upper_ipt <- 1.32042000 + 1.96 * 0.3236910

metaBeta_sex <- 1.09063000
metaCI_lower_sex <- 1.09063000 - 1.96 * 0.1113570
metaCI_upper_sex <- 1.09063000 + 1.96 * 0.1113570


# Forest plot ----

forest.df <- data.frame(ancestry = c("META", "SAS", "EAS", "AMR", "AFR", "EUR"),
                        index=1:6,
                        effect = c(metaBeta, sasBeta, easBeta, amrBeta, afrBeta, eurBeta),
                        lower = c(metaCI_lower, sasCI_lower, easCI_lower, amrCI_lower, afrCI_lower, eurCI_lower),
                        upper = c(metaCI_upper, sasCI_upper, easCI_upper, amrCI_upper, afrCI_upper, eurCI_upper))

# create forest plot
p <- forest.df %>%
  ggplot(aes(y = index, x = effect, xmin = lower, xmax = upper)) +
  geom_point() + 
  geom_vline(xintercept = 0, linetype = "dashed", 
             color = "gray", linewidth = 1) +
  geom_errorbarh(height = 0.1) +
  scale_y_continuous(name = "Genetic Ancestry", breaks=1:nrow(forest.df), labels=forest.df$ancestry) +
  theme_minimal() +
  labs(title = "Within-Ancestry Major Depression-PRS Effect Sizes",
       x = "Major Depression-PRS Effect Size") +
  xlim(-0.3, 3.25) +
  theme(text = element_text(family = "Times",
                            size = 12),
        plot.title = element_text(hjust = 0.5))

## IPT ----

forest.df1 <- data.frame(ancestry = c("META", "SAS", "EAS", "AMR", "AFR", "EUR"),
                        index=1:6,
                        effect = c(metaBeta_ipt, sasBeta_ipt, easBeta_ipt, amrBeta_ipt, afrBeta_ipt, eurBeta_ipt),
                        lower = c(metaCI_lower_ipt, sasCI_lower_ipt, easCI_lower_ipt, amrCI_lower_ipt, afrCI_lower_ipt, eurCI_lower_ipt),
                        upper = c(metaCI_upper_ipt, sasCI_upper_ipt, easCI_upper_ipt, amrCI_upper_ipt, afrCI_upper_ipt, eurCI_upper_ipt))

# create forest plot
p1 <- forest.df1 %>%
  ggplot(aes(y = index, x = effect, xmin = lower, xmax = upper)) +
  geom_point() + 
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray", linewidth = 1) +
  geom_errorbarh(height = 0.1) +
  scale_y_continuous(name = "Genetic Ancestry", breaks=1:nrow(forest.df1), labels=forest.df1$ancestry) +
  theme_minimal() +
  labs(title = "Within-Ancestry IPT Exposure Effect Sizes",
       x = "IPT Exposure Effect Sizes") +
  xlim(-0.3, 3.25) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(family = "Times New Roman"))

## Biological Sex ----

forest.df2 <- data.frame(ancestry = c("META", "SAS", "EAS", "AMR", "AFR", "EUR"),
                         index=1:6,
                         effect = c(metaBeta_sex, sasBeta_sex, easBeta_sex, amrBeta_sex, afrBeta_sex, eurBeta_sex),
                         lower = c(metaCI_lower_sex, sasCI_lower_sex, easCI_lower_sex, amrCI_lower_sex, afrCI_lower_sex, eurCI_lower_sex),
                         upper = c(metaCI_upper_sex, sasCI_upper_sex, easCI_upper_sex, amrCI_upper_sex, afrCI_upper_sex, eurCI_upper_sex))

# create forest plot
p2 <- forest.df2 %>%
  ggplot(aes(y = index, x = effect, xmin = lower, xmax = upper)) +
  geom_point() + 
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray", size = 1) +
  geom_errorbarh(height = 0.1) +
  scale_y_continuous(name = "Genetic Ancestry", breaks=1:nrow(forest.df2), labels=forest.df2$ancestry) +
  theme_minimal() +
  labs(title = "Within-Ancestry Biological Sex Effect Sizes",
       x = "Biological Sex Effect Sizes") +
  xlim(-0.3, 3.25) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(family = "Times New Roman"))

patchwork <- p / p1 / p2
patchwork +
  plot_annotation(title = "Within-Ancestry Depressive Symptoms-PRS, IPT Exposure, and Sex Effect Sizes",
                  tag_levels = 'A',
                  theme = theme(plot.title = element_text(hjust = 0.5),
                                text = element_text(family = "Times New Roman")))

q()

