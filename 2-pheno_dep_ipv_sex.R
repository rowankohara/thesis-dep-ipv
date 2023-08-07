# Rowan K. O'Hara
# 04.05.2023
# Thesis Project
#
# This script investigates the phenotypic associations between depression symptoms,
# interpersonal violence, and sex.
################################################################################


# LOAD PACKAGES & READ IN DATA #################################################

library(tidyverse)
library(report)
library(modelsummary)
library(sensemakr)

s4s <- read_csv("../../data/thesis_s4s_clean.csv")
# N = 7502


# STATISTICS ###################################################################

### t-test ####

##### depression ----
female_dep <- s4s %>%
  filter(biosex == 2) %>%
  select(depScore) %>%
  as.vector()

mean(female_dep$depScore)
# mean = 9.50

male_dep <- s4s %>%
  filter(biosex == 1) %>%
  select(depScore) %>%
  as.vector()

mean(male_dep$depScore)
# mean = 8.21

dep_ttest <- t.test(x = female_dep$depScore,
                    y = male_dep$depScore, 
                    alternative = "greater",
                    paired = FALSE)
# p value = 3.22e-46
# t = 14.35

sdDepF <- s4s %>%
  filter(biosex == 2) %>%
  summarise(depScore = sd(depScore))

sdDepF <- as.numeric(sdDepF)
# sd = 3.86

sdDepM <- s4s %>%
  filter(biosex == 1) %>%
  summarise(depScore = sd(depScore))

sdDepM <- as.numeric(sdDepM)
# sd = 3.64

report(dep_ttest)
glance(dep_ttest)


##### ipv ----

no_ipv_dep <- s4s %>%
  filter(ever_ipv == 0) %>%
  select(depScore) %>%
  as.vector()

mean(no_ipv_dep$depScore)
# mean = 8.38

ipv_dep <- s4s %>%
  filter(ever_ipv == 1) %>%
  select(depScore) %>%
  as.vector()

mean(ipv_dep$depScore)
# mean = 10.11

dep_ttest <- t.test(x = ipv_dep$depScore,
                    y = no_ipv_dep$depScore, 
                    alternative = "greater",
                    paired = FALSE)
# p value = 1.15e-78
# t = 19.035

sdDepN <- s4s %>%
  filter(ever_ipv == 0) %>%
  summarise(depScore = sd(depScore))

sdDepN <- as.numeric(sdDepN)
# sd = 3.61

sdDepIPV <- s4s %>%
  filter(ever_ipv == 1) %>%
  summarise(depScore = sd(depScore))

sdDepIPV <- as.numeric(sdDepIPV)
# sd = 3.95

report(dep_ttest)
glance(dep_ttest)


### prop-test ####

##### ipv ----
s4s_table <- s4s %>%
  mutate(ever_ipv = as.factor(ever_ipv),
         sex_ipv = as.factor(sex_ipv),
         phys_ipv = as.factor(phys_ipv)) %>%
  mutate(ever_ipv = fct_recode(ever_ipv, Yes = "1", No = "0"),
         sex_ipv = fct_recode(sex_ipv, Yes = "1", No = "0"),
         phys_ipv = fct_recode(phys_ipv, Yes = "1", No = "0"),
         biosex = as.factor(biosex)) %>%
  mutate(biosex = fct_recode(biosex, Male = "1", Female = "2"))

# Size of female and male samples
female_n <- s4s_table %>%
  filter(biosex == "Female") %>%
  nrow() %>%
  as.numeric()
# N = 4829

male_n <- s4s_table %>%
  filter(biosex == "Male") %>%
  nrow() %>%
  as.numeric()
# N = 2673

# Size of those reporting IPV by sex
female_ipv <- s4s_table %>%
  filter(biosex == "Female" & ever_ipv == "Yes") %>%
  nrow() %>%
  as.numeric()
# n = 2018

male_ipv <- s4s_table %>%
  filter(biosex == "Male" & ever_ipv == "Yes") %>%
  nrow() %>%
  as.numeric()
# n = 857

s4s_table %>%
  count(ever_ipv)
# A tibble: 2 × 2
# ever_ipv     n
# <fct>    <int>
# 1 No        4627
# 2 Yes       2875


# Two proportion z-test
ipv_proptest <- prop.test(x = c(female_ipv, male_ipv),
                          n = c(female_n, male_n),
                          alternative = "greater")
# 42% vs 32%
# p value = 6.43e-17
# chi squared = 68.5

glance(ipv_proptest)


##### ipv pa ----

# Size of those reporting IPV by sex
female_pa <- s4s_table %>%
  filter(biosex == "Female" & phys_ipv == "Yes") %>%
  nrow() %>%
  as.numeric()
# n = 1039

male_pa <- s4s_table %>%
  filter(biosex == "Male" & phys_ipv == "Yes") %>%
  nrow() %>%
  as.numeric()
# n = 726


# Two proportion z-test
pa_proptest <- prop.test(x = c(female_pa, male_pa),
                         n = c(female_n, male_n),
                         alternative = "less")
# 22% vs 27%
# p value = 0.0000000199
# chi squared = 30.2

glance(pa_proptest)


##### ipv sa ----

# Size of those reporting IPV by sex
female_sa <- s4s_table %>%
  filter(biosex == "Female" & sex_ipv == "Yes") %>%
  nrow() %>%
  as.numeric()
# n = 1564

male_sa <- s4s_table %>%
  filter(biosex == "Male" & sex_ipv == "Yes") %>%
  nrow() %>%
  as.numeric()
# n = 291


# Two proportion z-test
sa_proptest <- prop.test(x = c(female_sa, male_sa),
                         n = c(female_n, male_n),
                         alternative = "greater")
# 32% vs 10%
# p value = 5.45e-95
# chi squared = 426

glance(sa_proptest)


## plots ----

ggplot(s4s_table, aes(x = ever_ipv, y = depScore, fill = biosex)) +
  geom_boxplot(alpha = 0.7) +
  stat_summary(fun = mean, geom = "text", col = "black", position = position_dodge(0.75), vjust = -0.1, aes(label = round(after_stat(y), digits = 2))) +
  theme_minimal() +
  labs(title = "Sex Differences in IPT Exposure and Depressive Symptoms Score",
       x = "IPT Exposure",
       y = "Depressive Symptoms Score",
       fill = "Biological Sex") +
  theme(text = element_text(family = "Times",
                            size = 16)) +
  scale_fill_manual(values = c("#40B0A6", "#E1BE6A"))
# Color-blind friendly colors

ggplot(s4s_table, aes(x = depScore, y = phys_ipv, fill = biosex)) +
  geom_boxplot(alpha = 0.7) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Sex Differences in IPV Exposure and Depression Symptoms",
       x = "Depression Symptom Score",
       y = "Physical Assault IPV Exposure",
       fill = "Biological Sex") +
  theme(text = element_text(family = "Times",
                            size = 16)) +
  scale_fill_manual(values = c("#40B0A6", "#E1BE6A"))

ggplot(s4s_table, aes(x = depScore, y = sex_ipv, fill = biosex)) +
  geom_boxplot(alpha = 0.7) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Sex Differences in IPV Exposure and Depression Symptoms",
       x = "Depression Symptom Score",
       y = "Sexual Assault IPV Exposure",
       fill = "Biological Sex") +
  theme(text = element_text(family = "Times",
                            size = 16)) +
  scale_fill_manual(values = c("#40B0A6", "#E1BE6A"))


# LINEAR REGRESSION ############################################################
s4s_model1 <- lm(depScore ~ biosex + ever_ipv + biosex*ever_ipv, data = s4s_table)

summary(s4s_model1)

# Run these to get info for table
glance(s4s_model1)

summary(s4s_model1)$coefficient

report(s4s_model1)

# Odds Ratio
exp(coef(s4s_model1))
# biosexFemale              2.867680
# ever_ipvYes               4.380554
# biosexFemale:ever_ipvYes  1.239828

# for confidence intervals
exp(confint(s4s_model1))

partial_r2(s4s_model1)
# biosexFemale              0.011773133
# ever_ipvYes               0.012209921
# biosexFemale:ever_ipvYes  0.000175003


### exploratory analyses ----
# Look at the effect of sex on physical and sexual ipv exposure
s4s_model2 <- glm(depScore ~ biosex + phys_ipv + biosex*phys_ipv, data = s4s_table)
s4s_model3 <- glm(depScore ~ biosex + sex_ipv + biosex*sex_ipv, data = s4s_table)

summary(s4s_model2)
summary(s4s_model3)

glance(s4s_model2)
glance(s4s_model3)

summary(s4s_model2)$coefficient
summary(s4s_model3)$coefficient

report(s4s_model2)
report(s4s_model3)


### plots ----
s4s_table %>%
  ggplot(aes(x=biosex, y=depScore)) + 
    geom_jitter(aes(color = biosex)) +
    stat_smooth(method="lm", color="black", se=FALSE) +
  theme_minimal() +
  labs(title = "Sex Differences in IPT Exposure and Depression Symptoms",
       x = "Depression Symptom Score",
       y = "IPT Exposure",
       color = "Biological Sex") +
  theme(text = element_text(family = "Times",
                            size = 16)) +
  scale_color_manual(values = c("#40B0A6", "#E1BE6A")) +
  facet_grid(ever_ipv ~ .)






