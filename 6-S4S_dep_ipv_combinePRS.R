# Rowan K. O'Hara
# 07.09.2023
# Thesis Project
#
# Combine PRS scores
#
# NOTE (if run on Fenn):
# Need to run interactively on Fenn using R in the terminal
# Push to terminal using Option + Command + Return
################################################################################


# Run on Fenn

cd /vcu_gpfs2/home/oharark/S4S/PRS/dep/output/

R

# LOAD PACKAGES & READ IN DATA #################################################

library(tidyverse)
library(readr)
library(data.table)

datadir <- "/home/oharark/S4S/PRS/dep/output/" # set datadir
j <- "dep"
z <- "_S4S_" ## dataset


# CHANGE FILE STRUCTURES #######################################################

# merging all the files into one to get the sum of the scores
PRS_list  <- list.files(".", full.names=F, recursive = FALSE) ## read the files and list them

PRS_list <- PRS_list[grepl(".sscore", PRS_list)] ## grep the scores

PRS_list <- PRS_list[!grepl(".vars", PRS_list)] ## get rid of the variants file


df <-  PRS_list %>% ## name them df
  set_names(.) %>% ##set the names for them 
  map_df(read_table2, .id = "FileName") ## read them

colnames(df) <- c("FileName", "FID" , "IID","NMISS_ALLELE_CT","NAMED_ALLELE_DOSAGE","SCORE1_AVG") # rename the columns
##getting the sum of the scores

summary_ALL <- df %>% ## combine the scores
  group_by(IID) %>%  ## group them by ID for each individual in each chr file
  summarise(num = n(), ## summarize them 
            PRS_SCORE = sum(SCORE1_AVG)) ## based on their score

summary_ALL <- summary_ALL[,c(1,3)] #get only the ID and the score

outFile <- paste0("All_",j,z,"PRS.summary",sep="") # name the output file

write.table(summary_ALL,outFile,quote = F,row.names = F) # write it. 

q()
