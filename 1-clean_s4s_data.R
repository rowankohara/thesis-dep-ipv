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
# EAS = Easy Asian ancestry
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

# FAM #s are correct
s4s_qc %>%
  count(fam)


# CLEAN DEPRESSION VARAIBLES ###################################################

## exclusion -------------------------------------------------------------------
dep_s4s <- s4s_qc %>%
  rowwise() %>%               # https://dplyr.tidyverse.org/articles/rowwise.html
  mutate(dep_noAnswer = sum(c_across(c(y1f_hea_1c,y1f_hea_1d, y1f_hea_1e, y1f_hea_1g)) == -99), # https://stackoverflow.com/questions/51783095/compute-row-wise-counts-in-subsets-of-columns-in-dplyr
         dep_skip = sum(c_across(c(y1f_hea_1c,y1f_hea_1d, y1f_hea_1e, y1f_hea_1g)) == -9),
         dep_filter = sum(c_across(dep_noAnswer:dep_skip))) %>%
  filter(dep_filter < 3)
# N = 7732


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
dep_s4s %>%
  summarize(min_score = min(dep_s4s$depScore), 
          mean_score = mean(dep_s4s$depScore), 
          med_score = median(dep_s4s$depScore), 
          max_score = max(dep_s4s$depScore))


# CLEAN IPV VARIABLES ##########################################################

## create interpersonal violence scores ----------------------------------------
# Replace missing values with NAs
ipv_s4s <- dep_s4s %>%
  mutate(y1f_str_1b_before12 = if_else(y1f_str_1b_before12 == -99, NA_real_, y1f_str_1b_before12),
         y1f_str_1b_past12 = if_else(y1f_str_1b_past12 == -99, NA_real_, y1f_str_1b_past12),
         Y1F_str_1b_beforeVCU = if_else(Y1F_str_1b_beforeVCU == -99, NA_real_, Y1F_str_1b_beforeVCU),
         y1s_str_1b_sincevcu = if_else(y1s_str_1b_sincevcu == -99, NA_real_, y1s_str_1b_sincevcu),
         y1f_str_1c_before12 = if_else(y1f_str_1c_before12 == -99, NA_real_, y1f_str_1c_before12),
         y1f_str_1c_past12 = if_else(y1f_str_1c_past12 == -99, NA_real_, y1f_str_1c_past12),
         Y1F_str_1c_beforeVCU = if_else(Y1F_str_1c_beforeVCU == -99, NA_real_, Y1F_str_1c_beforeVCU),
         y1s_str_1c_sincevcu = if_else(y1s_str_1c_sincevcu == -99, NA_real_, y1s_str_1c_sincevcu),
         y1f_str_1d_before12 = if_else(y1f_str_1d_before12 == -99, NA_real_, y1f_str_1d_before12),
         y1f_str_1d_past12 = if_else(y1f_str_1d_past12 == -99, NA_real_, y1f_str_1d_past12),
         Y1F_str_1d_beforeVCU = if_else(Y1F_str_1d_beforeVCU == -99, NA_real_, Y1F_str_1d_beforeVCU),
         y1s_str_1d_sincevcu = if_else(y1s_str_1d_sincevcu == -99, NA_real_, y1s_str_1d_sincevcu),
         y1f_str_1b_never = if_else(y1f_str_1b_never == -99, NA_real_, y1f_str_1b_never),
         Y1S_str_1b_never = if_else(Y1S_str_1b_never == -99, NA_real_, Y1S_str_1b_never),
         y1f_str_1c_never = if_else(y1f_str_1c_never == -99, NA_real_, y1f_str_1c_never),
         Y1S_str_1c_never = if_else(Y1S_str_1c_never == -99, NA_real_, Y1S_str_1c_never),
         y1f_str_1d_never = if_else(y1f_str_1d_never == -99, NA_real_, y1f_str_1d_never),
         Y1S_str_1d_never = if_else(Y1S_str_1d_never == -99, NA_real_, Y1S_str_1d_never))


# Sum each category (physical assault, sexual assault, & other uncomfortable sexual events), ignoring missing data.
ipv_calc_s4s <- ipv_s4s %>%
  rowwise() %>%
  mutate(sumscore_PA = sum(y1f_str_1b_before12, y1f_str_1b_past12),
         sumscore_SA = sum(y1f_str_1c_before12, y1f_str_1c_past12),
         sumscore_other = sum(y1f_str_1d_before12, y1f_str_1d_past12)) %>%
  mutate(sumscore_PA = if_else(is.na(sumscore_PA), Y1F_str_1b_beforeVCU, sumscore_PA),
         sumscore_SA = if_else(is.na(sumscore_SA), Y1F_str_1c_beforeVCU, sumscore_SA),
         sumscore_other = if_else(is.na(sumscore_other), Y1F_str_1d_beforeVCU, sumscore_other))

# Looking at any occurrence of assault, replace any score above one with one
# 0 = No  
# 1 = Yes  
# NA = Missing data  
ipv_calc <- ipv_calc_s4s %>%
  mutate(sumscore_PA = if_else(sumscore_PA > 1, 1, sumscore_PA),
         sumscore_SA = if_else(sumscore_SA > 1, 1, sumscore_SA),
         sumscore_other = if_else(sumscore_other > 1, 1, sumscore_other))

# Do a rowwise summation of all sumscores to get the total score of at least one 
# occurrence of assault in any category.
ipv_calc$sumscore_total <- apply(ipv_calc[,c("sumscore_PA", "sumscore_SA", "sumscore_other")], 1, sumNA, na.rm=T)

# 1 = Never happened  
# 0 = Yes, it happened
ipv_calc %>%
  count(sumscore_total)


## calculate negating scores ---------------------------------------------------

never_dat <- ipv_calc %>%
  rowwise() %>%
  mutate(y1f_str_1b_never = if_else(is.na(y1f_str_1b_never), Y1S_str_1b_never, y1f_str_1b_never)) %>%
  mutate(y1f_str_1c_never = if_else(is.na(y1f_str_1c_never), Y1S_str_1c_never, y1f_str_1c_never)) %>%
  mutate(y1f_str_1d_never = if_else(is.na(y1f_str_1d_never), Y1S_str_1d_never, y1f_str_1d_never))


## check congruency ------------------------------------------------------------
# Participants were able to check more than one option for the IPV exposure questions, 
# and it was possible to check both “yes” and “never happened to me” to an IPV exposure.
# Participants who marked both “yes” and “never happened to me” are called “incongruent” for that question. 
# Participants will be excluded if all 3 answered IPV questions were incongruent. For IPV, 
# participants that had only one incongruent response were not excluded, but their incongruent 
# response was marked as missing. Participants that have congruent “never happened to me” answers and 
# any missingness will be excluded due to the uncertainty if they would have responded “yes” to the questions skipped.

# If `align` variable == 1, then `sumscore` and `never` variable do not equal (are congruent -- good!).
# If `align` variable == 0, then `sumscore` and `never` variable are equal (are incongruent -- bad!). 
align_dat <- never_dat %>%
  rowwise() %>%
  mutate(align_PA = if_else(sumscore_PA != y1f_str_1b_never, 1, 0),
         align_SA = if_else(sumscore_SA != y1f_str_1c_never, 1, 0),
         align_other = if_else(sumscore_other != y1f_str_1d_never, 1, 0))

align_dat$align_total <- apply(align_dat[,c("align_PA", "align_SA", "align_other")], 1, sumNA, na.rm=T)

align_dat %>%
  count(align_total) %>%
  knitr::kable()

# Count NAs for each participant
align_dat <- align_dat %>%
  rowwise() %>%
  mutate(ipv_NA = sum(is.na(c(sumscore_PA, sumscore_SA, sumscore_other))))

align_dat %>%
  count(ipv_NA) %>%
  knitr::kable()

# Remove participants that are incongruent for all questions answered and then 
# anyone whose `sumscore_total` equals zero and have any missingness.
ipv_final <- align_dat %>%
  rowwise() %>%
  mutate(sumscore_PA = if_else(sumscore_PA != y1f_str_1b_never, sumscore_PA, NA_real_),
         sumscore_SA = if_else(sumscore_SA != y1f_str_1c_never, sumscore_SA, NA_real_),
         sumscore_other = if_else(sumscore_other != y1f_str_1d_never, sumscore_other, NA_real_),
         ever_ipv = if_else(sumscore_total > 0, 1, 0)) %>%
  mutate(ipv_include = if_else(ipv_NA == 3, 0, 1)) %>%
  mutate(ipv_include = if_else(sumscore_total == 0 & ipv_NA != 0, 0, ipv_include)) %>%
  filter(ipv_include == 1)
# N = 7561

ipv_final %>%
  group_by(biosex) %>%
  count()
# Males: N = 2689
# Females: N = 4872


# OUTPUT #######################################################################

write_csv(ipv_final, "../../data/thesis_s4s_clean.csv",
          append = FALSE, 
          col_names = TRUE)




