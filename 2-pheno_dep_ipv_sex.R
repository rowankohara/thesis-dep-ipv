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

s4s <- read_csv("../../data/thesis_s4s_clean.csv")
# N = 7502


# STATISTICS ###################################################################

### t-test ####

##### depression ----
female_dep <- s4s %>%
  filter(biosex == 2) %>%
  select(depScore) %>%
  as.vector()
# mean = 9.50

male_dep <- s4s %>%
  filter(biosex == 1) %>%
  select(depScore) %>%
  as.vector()
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


### prop-test ####

##### ipv ----
s4s_table <- s4s %>%
  mutate(ever_ipv = as.factor(ever_ipv),
         SAscore = as.factor(if_else(SAscore > 0 | USEscore > 0, 1, 0)),
         PAscore = as.factor(PAscore)) %>%
  mutate(ever_ipv = fct_recode(ever_ipv, Yes = "1", No = "0"),
         SAscore = fct_recode(SAscore, Yes = "1", No = "0"),
         PAscore = fct_recode(PAscore, Yes = "1", No = "0"),
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
# N = 2875


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
  filter(biosex == "Female" & PAscore == "Yes") %>%
  nrow() %>%
  as.numeric()
# n = 1039

male_pa <- s4s_table %>%
  filter(biosex == "Male" & PAscore == "Yes") %>%
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
  filter(biosex == "Female" & SAscore == "Yes") %>%
  nrow() %>%
  as.numeric()
# n = 1564

male_sa <- s4s_table %>%
  filter(biosex == "Male" & SAscore == "Yes") %>%
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

ggplot(s4s_table, aes(x = depScore, y = ever_ipv, fill = biosex)) +
  geom_boxplot(alpha = 0.7) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Sex Differences in IPV Exposure and Depression Symptoms",
       x = "Depression Symptom Score",
       y = "IPV Exposure",
       fill = "Biological Sex") +
  theme(text = element_text(family = "Times",
                            size = 16)) +
  scale_fill_manual(values = c("#40B0A6", "#E1BE6A"))
# Color-blind friendly colors

ggplot(s4s_table, aes(x = depScore, y = PAscore, fill = biosex)) +
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

ggplot(s4s_table, aes(x = depScore, y = SAscore, fill = biosex)) +
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


# LOGISTIC REGRESSION ##########################################################
s4s_model1 <- glm(depScore ~ biosex, data = s4s_table)
s4s_model2 <- glm(depScore ~ ever_ipv, data = s4s_table)
s4s_model3 <- glm(depScore ~ biosex*ever_ipv, data = s4s_table)


models <- list(
  "Model 1"     = s4s_model1,
  "Model 2"     = s4s_model2,
  "Model 3"     = s4s_model3
)

modelsummary(models, 
             "markdown", 
             stars = TRUE,
             title = "Model Comparison")

summary(s4s_model1)
summary(s4s_model2)
summary(s4s_model3)


# Run these to get info for table
glance(s4s_model1)
glance(s4s_model2)
glance(s4s_model3)



summary(s4s_model1)$coefficient
summary(s4s_model2)$coefficient
summary(s4s_model3)$coefficient


report(s4s_model1)
report(s4s_model2)
report(s4s_model3)


### exploratory analyses ----
# Look at the effect of sex on physical and sexual ipv exposure
s4s_model4 <- glm(depScore ~ biosex + PAscore, data = s4s_table)
s4s_model5 <- glm(depScore ~ biosex + SAscore, data = s4s_table)

summary(s4s_model4)
summary(s4s_model5)

glance(s4s_model4)
glance(s4s_model5)

summary(s4s_model4)$coefficient
summary(s4s_model5)$coefficient

report(s4s_model4)
report(s4s_model5)


### plots ----
s4s_table %>%
  mutate(ever_ipv = as.numeric(ever_ipv)) %>%
  mutate(ever_ipv = if_else(ever_ipv == 1, 0, 1)) %>%
  ggplot(aes(x=depScore, y=ever_ipv)) + 
    geom_point(aes(color = biosex)) +
    stat_smooth(method="glm", color="black", se=FALSE,
                method.args = list(family=binomial)) +
  theme_minimal() +
  labs(title = "Sex Differences in IPV Exposure and Depression Symptoms",
       x = "Depression Symptom Score",
       y = "IPV Exposure",
       color = "Biological Sex") +
  theme(text = element_text(family = "Times",
                            size = 16)) +
  scale_color_manual(values = c("#40B0A6", "#E1BE6A"))






