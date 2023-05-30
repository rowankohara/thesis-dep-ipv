# Rowan K. O'Hara
# 03.25.2023
# Thesis Project
#
# This script cleans the raw Spit for Science data using the exclusions outlined  
# in the preregistration.
################################################################################


# LOAD PACKAGES & READ IN DATA #################################################

library(tidyverse)
library(bazar)
library(skimr)

## raw Spit for Science data ---------------------------------------------------
s4s <- read_csv("../../data/S4S_rawExport_withCohort5&Upperclassman2018&Soph2019_081420.csv")
# N = 12358


## FAM & EVEC files by ancestry ------------------------------------------------
# AMR = Admixed American ancestry
# AFR = African Ancestry
# EAS = East Asian ancestry
# EUR = European ancestry
# SAS = South Asian ancestry

famAMR <- read.table("../../data/S4S_cohort_15_AMR_info8_maf05_hwe6_geno025_mind025_final.fam", 
                     col.names= c("famid", "iid", "dad", "mom", "biosex", "phenoval"))

famAFR <- read.table("../../data/S4S_cohort_15_AFR_info8_maf05_hwe6_geno025_mind025_final.fam",
                     col.names= c("famid", "iid", "dad", "mom", "biosex", "phenoval"))

famEAS <- read.table("../../data/S4S_cohort_15_EAS_info8_maf05_hwe6_geno025_mind025_final.fam",
                     col.names= c("famid", "iid", "dad", "mom", "biosex", "phenoval"))

famEUR <- read.table("../../data/S4S_cohort_15_EUR_info8_maf05_hwe6_geno025_mind025_final.fam",
                     col.names= c("famid", "iid", "dad", "mom", "biosex", "phenoval"))

famSAS <- read.table("../../data/S4S_cohort_15_SAS_info8_maf05_hwe6_geno025_mind025_final.fam",
                     col.names= c("famid", "iid", "dad", "mom", "biosex", "phenoval"))

evecAMR <- read.table("../../data/Cohort_15_AMR_PCs.txt",
                      header = TRUE)

evecAFR <- read.table("../../data/Cohort_15_AFR_PCs.txt",
                      header = TRUE)

evecEAS <- read.table("../../data/Cohort_15_EAS_PCs.txt",
                      header = TRUE)

evecEUR <- read.table("../../data/Cohort_15_EUR_PCs.txt",
                      header = TRUE)

evecSAS <- read.table("../../data/Cohort_15_SAS_PCs.txt",
                      header = TRUE)


# GENETIC EXCLUSIONS ###########################################################
# Only include participants with genetic data available who passed QC

## combine EVEC and FAM files --------------------------------------------------
afr <- left_join(famAFR, evecAFR, by = c("famid" = "FID"))
amr <- left_join(famAMR, evecAMR, by = c("famid" = "FID"))
eas <- left_join(famEAS, evecEAS, by = c("famid" = "FID"))
eur <- left_join(famEUR, evecEUR, by = c("famid" = "FID"))
sas <- left_join(famSAS, evecSAS, by = c("famid" = "FID"))

# Add prefix to unique_id of S4S df (to make it easier to join with fam and evec dfs)
s4s$unique_id <- paste0("S4S_", s4s$unique_id)

# Combine QC first
afr$fam <- "AFR"
amr$fam <- "AMR"
eas$fam <- "EAS"
eur$fam <- "EUR"
sas$fam <- "SAS"

all <- bind_rows(afr, amr, eas, eur, sas)


## join genetic data to Spit for Science data ----------------------------------
s4s_qc <- right_join(s4s, all, by = c("unique_id" = "famid"))
# N = 9588
# Removed 2770

# FAM #s are correct
s4s_qc %>%
  count(fam)
# A tibble: 5 × 2
# fam       n
# <chr> <int>
# 1 AFR    2040
# 2 AMR    1198
# 3 EAS     929
# 4 EUR    4594
# 5 SAS     827


# CLEAN DEPRESSION VARAIBLES ###################################################

## exclusion -------------------------------------------------------------------
dep_s4s <- s4s_qc %>%
  rowwise() %>%               # https://dplyr.tidyverse.org/articles/rowwise.html
  mutate(dep_noAnswer = sum(c_across(c(y1f_hea_1c,y1f_hea_1d, y1f_hea_1e, y1f_hea_1g)) == -99), # https://stackoverflow.com/questions/51783095/compute-row-wise-counts-in-subsets-of-columns-in-dplyr
         dep_skip = sum(c_across(c(y1f_hea_1c,y1f_hea_1d, y1f_hea_1e, y1f_hea_1g)) == -9),
         dep_filter = sum(c_across(dep_noAnswer:dep_skip))) %>%
  filter(dep_filter < 3)
# N = 7732
# Removed 1856


## calculate depressive symptom scores -----------------------------------------
# These scores are pro-rated sum scores; so after filtering out
# participants who did not answer the threshold number of questions (2) for each variable, the skip
# choose not to answer values are converted to zeros. Then, the filter variable, dep_filter,
# is converted into a counter of the number of questions participant answered by subtracting the amount of
# questions they did not answer from the total questions (4). The symptom scores are calculated by summing all
# of the questions for depression. If the participant did not answer all questions (filter variable
# does not equal 4), their score is pro-rated by finding the average of answered questions and multiplying by
# the total number of questions (4).

dep_s4s <- dep_s4s %>%
  mutate(y1f_hea_1c = if_else(y1f_hea_1c == -9 | y1f_hea_1c == -99, 0, y1f_hea_1c),
         y1f_hea_1d = if_else(y1f_hea_1d == -9 | y1f_hea_1d == -99, 0, y1f_hea_1d),
         y1f_hea_1e = if_else(y1f_hea_1e == -9 | y1f_hea_1e == -99, 0, y1f_hea_1e),
         y1f_hea_1g = if_else(y1f_hea_1g == -9 | y1f_hea_1g == -99, 0, y1f_hea_1g)) %>%  # https://dplyr.tidyverse.org/reference/if_else.html
  mutate(dep_filter = 4 - dep_filter) %>%
  rowwise() %>%
  mutate(depScore = sum(c_across(c(y1f_hea_1c,y1f_hea_1d, y1f_hea_1e, y1f_hea_1g)))) %>%  # https://dplyr.tidyverse.org/reference/c_across.html
  mutate(depScore = if_else(dep_filter != 4, (depScore / dep_filter) * 4, depScore))

# Check scores
summary(dep_s4s$depScore)
# min = 4
# 1st Q = 6
# median = 8
# mean = 9.059
# 3rd Q = 11
# max = 20


# CLEAN IPV VARIABLES ##########################################################

## create interpersonal violence scores ----------------------------------------
# Replace missing values with NAs
dat <- dep_s4s %>%
  mutate_at(names(dep_s4s[,grep("_str_1[b-d]", names(dep_s4s), ignore.case = T)]),
            funs(na_if(., -99)))

# Do rowwise summation of questions for each category (physical assault, sexual 
# assault, & other uncomfortable sexual events), ignoring missing data.
# 0 = No, 1 = Yes, NA = Missing data 
ipv_dat <- dat %>%
  mutate(phys_ipv = pmax(y1f_str_1b_before12, y1f_str_1b_past12, Y1F_str_1b_beforeVCU, na.rm= T),
         sex_ipv = pmax(y1f_str_1c_before12, y1f_str_1c_past12, Y1F_str_1c_beforeVCU, na.rm= T),
         other_ipv = pmax(y1f_str_1d_before12, y1f_str_1d_past12, Y1F_str_1d_beforeVCU, na.rm= T))

# Do a rowwise summation of all sumscores to get the total score of at least one occurrence of assault in any category.
ipv_dat$ever_ipv <- apply(ipv_dat[,c("phys_ipv", "sex_ipv", "other_ipv")], 1, sumNA, na.rm=T)

ipv_dat %>%
  count(ever_ipv)
# A tibble: 5 × 2
# Rowwise: 
# ever_ipv     n
# <dbl> <int>
# 1        0  4749
# 2        1  1883
# 3        2   709
# 4        3   342
# 5       NA    49


## calculate negating scores ---------------------------------------------------
# Since the questions were "check one or more of the boxes," we need to check if 
# any participants selected both that they did have an experience and that they 
# never have (in the same question). Something else to consider is that the spring 
# enrollees have a different `never` variable. Here, I collapsed them into a single 
# variable for both fall and spring enrollees.
# 1 = Never happened, 0 = Yes, it happened, NA = Missing data
ipv_calc <- ipv_dat %>%
  rowwise() %>%
  mutate(y1f_str_1b_total = if_else(is.na(y1f_str_1b_never), Y1S_str_1b_never, y1f_str_1b_never)) %>%
  mutate(y1f_str_1c_total = if_else(is.na(y1f_str_1c_never), Y1S_str_1c_never, y1f_str_1c_never)) %>%
  mutate(y1f_str_1d_total = if_else(is.na(y1f_str_1d_never), Y1S_str_1d_never, y1f_str_1d_never))

ipv_calc %>% 
  count(y1f_str_1b_total, y1f_str_1c_total, y1f_str_1d_total)


## check congruency ------------------------------------------------------------
# If their answers are impossible (incongruent), mark as missing
# Roll sexual assault and other uncomfortable experiends into one variable
# Count how many missing variables there are
# Make ever_ipv into a binary variable
# Mark participants who will not be included in our final sample if they are missing all 3 ipv variables
# and/or they report no ipv and have any missingness.
ipv_qc <- ipv_calc %>%
  rowwise() %>%
  mutate(phys_ipv = if_else(phys_ipv == 1 & y1f_str_1b_total == 1, NA_real_, phys_ipv),
         sex_ipv = if_else(sex_ipv == 1 & y1f_str_1c_total == 1, NA_real_, sex_ipv),
         other_ipv = if_else(other_ipv == 1 & y1f_str_1d_total == 1, NA_real_, other_ipv)) %>%
  mutate(sex_ipv = if_else(sex_ipv == 1 | other_ipv == 1, 1, sex_ipv)) %>%
  mutate(ipv_NA = sum(is.na(c(phys_ipv, sex_ipv, other_ipv)))) %>%
  mutate(ipv_total = sumNA(c(phys_ipv, sex_ipv, other_ipv), na.rm = T)) %>%
  mutate(ever_ipv = if_else(ipv_total > 0, 1, 0)) %>%
  mutate(ipv_include = if_else(ipv_NA == 3, 0, 1)) %>%
  mutate(ipv_include = if_else(ever_ipv == 0 & ipv_NA != 0, 0, ipv_include)) %>%
  filter(ipv_include == 1)

ipv_qc %>%
  count(phys_ipv, sex_ipv, other_ipv, ipv_NA, ipv_total, ever_ipv, ipv_include)
# A tibble: 13 × 8
# Rowwise: 
# phys_ipv sex_ipv other_ipv ipv_NA ipv_total ever_ipv ipv_include     n
# <dbl>   <dbl>     <dbl>  <int>     <dbl>    <dbl>       <dbl> <int>
# 1        0       0         0      0         0        0           1  4627
# 2        0       1         0      0         1        1           1    56
# 3        0       1         1      0         2        1           1  1008
# 4        0       1        NA      1         1        1           1     7
# 5        1       0         0      0         1        1           1   957
# 6        1       1         0      0         2        1           1    48
# 7        1       1         1      0         3        1           1   686
# 8        1       1        NA      1         2        1           1    11
# 9        1      NA         0      1         1        1           1    12
# 10       1      NA        NA      2         1        1           1    51
# 11      NA       1         0      1         1        1           1     4
# 12      NA       1         1      1         2        1           1    29
# 13      NA       1        NA      2         1        1           1     6

ipv_qc %>%
  count(ipv_include, ever_ipv)
# A tibble: 2 × 3
# Rowwise: 
# ipv_include ever_ipv     n
# <dbl>    <dbl> <int>
# 1           1        0  4627
# 2           1        1  2875
#
# N = 7502
# Removed 230

ipv_qc %>%
  count(biosex)
# A tibble: 2 × 2
# Rowwise: 
# biosex     n
# <int> <int>
# 1      1  2673
# 2      2  4829

ipv_qc %>%
  count(biosex, ever_ipv)
# A tibble: 4 × 3
# Rowwise: 
# biosex ever_ipv     n
# <int>    <dbl> <int>
# 1      1        0  1816
# 2      1        1   857
# 3      2        0  2811
# 4      2        1  2018


# SUBSET #######################################################################

final_dat <- ipv_qc[, c("unique_id", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6",
                          "PC7", "PC8", "PC9", "PC10", "fam", "biosex", "depScore",
                        "ever_ipv", "phys_ipv", "sex_ipv")]

final_dat <- final_dat %>%
  rename("ancestry" = "fam")


# OUTPUT #######################################################################

write_csv(final_dat, "../../data/thesis_s4s_clean.csv",
          append = FALSE, 
          col_names = TRUE)

