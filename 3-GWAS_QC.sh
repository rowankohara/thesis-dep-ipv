# Rowan K. O'Hara
# 07.07.2023
# Thesis Project
#
# GWAS Summary Statistics QC for PRS-CSx
# NOTE (if run on Fenn):
# Need to run interactively on Fenn using R in the terminal
# Push to terminal using Option + Command + Return
################################################################################

# Run on Fenn

cd /vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats
  
R

# LOAD PACKAGES & READ IN DATA #################################################

library(tidyverse)
library(data.table)

afr <- fread("/vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/MDD.AFR.MVP.NatNeuro2021")
eas <- fread("/vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/jamapsy_Giannakopoulou_2021_exclude_whi_23andMe_ukb.txt")
eur <- fread("/vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/MDD_meta_AA_FIN_QC.txt")


# CHANGE FILE STRUCTURES #######################################################

# eur has already been QC'd
# SNPS - SNP rsid, A1 - allele 1, A2 - allele 2, BETA - effect, P - p value

afr <- afr %>%
  select(c(rsid, A1, A2, EFFECT, P)) %>%
  rename("SNP" = "rsid",
         "BETA" = "EFFECT")
dim(afr) # 7,410,896 SNPs

eas <- eas %>%
  select(c(MarkerName, Allele1, Allele2, Effect, P.SE)) %>%
  rename("SNP" = "MarkerName",
         "A1" = "Allele1",
         "A2" = "Allele2",
         "BETA" = "Effect",
         "P" = "P.SE")
dim(eas) # 7,440,922 SNPs

eur <- eur %>%
  select(c(SNPS, A1, A2, BETA, P)) %>%
  rename("SNP" = "SNPS")
dim(eur) # 5,326,424 SNPs


# QC ###########################################################################

# filter  sumafr based on MAF (MAF <0.01), p-vals <0 or >1, SNPs with missing values,
# Remove variants that are not SNPs or are strand-ambiguous, info < 0.9, SNPs with duplicated rs

## remove p-vals <0 or >1 ----

afr.r <- afr[!afr$P < 0 & !afr$P > 1,]
eas.r <- eas[!eas$P < 0 & !eas$P > 1,]
eur.r <- eur[!eur$P < 0 & !eur$P > 1,]

dim(afr)[1]-dim(afr.r)[1] # 0 dropped
dim(eas)[1]-dim(eas.r)[1] # 0 dropped
dim(eur)[1]-dim(eur.r)[1] # 0 dropped


## remove missing SNPs ----

afr.r2 <- afr.r[is.na(afr.r$P) == F, ]
eas.r2 <- eas.r[is.na(eas.r$P) == F, ]
eur.r2 <- eur.r[is.na(eur.r$P) == F, ]

dim(afr.r)[1]-dim(afr.r2)[1] # 0 dropped
dim(eas.r)[1]-dim(eas.r2)[1] # 0 dropped
dim(eur.r)[1]-dim(eur.r2)[1] # 0 dropped


## remove variants that are not SNPs ----

afr.r3 <- afr.r2[grep("rs", afr.r2$SNP),]
eas.r3 <- eas.r2[grep("rs", eas.r2$SNP),]
eur.r3 <- eur.r2[grep("rs", eur.r2$SNP),]

dim(afr.r2)[1]-dim(afr.r3)[1] # 14,514 dropped
dim(eas.r2)[1]-dim(eas.r3)[1] # 11,864 dropped
dim(eur.r2)[1]-dim(eur.r3)[1] # 0 dropped

## remove duplicated SNPs ----
afr.final <- afr.r3[duplicated(afr.r3$SNP) == F,]
eas.final <- eas.r3[duplicated(eas.r3$SNP) == F,]
eur.final <- eur.r3[duplicated(eur.r3$SNP) == F,]

dim(afr.r3)[1]-dim(afr.final)[1] # 1,825 dropped
dim(eas.r3)[1]-dim(eas.final)[1] # 3,383 dropped
dim(eur.r3)[1]-dim(eur.final)[1] # 0 dropped

# ----
eas.final$A1 <- toupper(eas.final$A1)
eas.final$A2 <- toupper(eas.final$A2)

dim(afr.final) # 7,394,557
dim(eas.final) # 7,425,675
dim(eur.final) # 5,326,424

## write output ----
write.table(afr.final,"/vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/MDD.AFR.MVP.NatNeuro2021_half_QCd",
            row.names=F,col.names=T,quote=F,sep="\t")
write.table(eur.final,"/vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/MDD_meta_AA_FIN_QC_half_QCd",
            row.names=F,col.names=T,quote=F,sep="\t")
write.table(eas.final,"/vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/jamapsy_Giannakopoulou_2021_exclude_whi_23andMe_ukb_half_QCd",
            row.names=F,col.names=T,quote=F,sep="\t")

q()


################################################################################
# in bash


# remove ambiguous SNPs 
awk '!( ($2=="A" && $3=="T") || \
        ($2=="T" && $3=="A") || \
        ($2=="G" && $3=="C") || \
        ($2=="C" && $3=="G")) {print}' \
        /vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/MDD.AFR.MVP.NatNeuro2021_half_QCd > \
        /vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/MDD.AFR.MVP.NatNeuro2021_QCd_final.txt

# remove ambiguous SNPs 
awk '!( ($2=="A" && $3=="T") || \
        ($2=="T" && $3=="A") || \
        ($2=="G" && $3=="C") || \
        ($2=="C" && $3=="G")) {print}' \
        /vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/jamapsy_Giannakopoulou_2021_exclude_whi_23andMe_ukb_half_QCd > \
        /vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/jamapsy_Giannakopoulou_2021_exclude_whi_23andMe_ukb_QCd_final.txt
        
# remove ambiguous SNPs 
awk '!( ($2=="A" && $3=="T") || \
        ($2=="T" && $3=="A") || \
        ($2=="G" && $3=="C") || \
        ($2=="C" && $3=="G")) {print}' \
        /vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/MDD_meta_AA_FIN_QC_half_QCd > \
        /vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/MDD_meta_AA_FIN_QCd_final.txt
# check #s removed



wc -l /vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/MDD.AFR.MVP.NatNeuro2021_half_QCd
# 7394558
wc -l /vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/MDD.AFR.MVP.NatNeuro2021_QCd_final.txt
# 6356949
wc -l /vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/MDD.AFR.MVP.NatNeuro2021
# 7410898 originally

wc -l /vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/jamapsy_Giannakopoulou_2021_exclude_whi_23andMe_ukb_half_QCd
# 7425676
wc -l /vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/jamapsy_Giannakopoulou_2021_exclude_whi_23andMe_ukb_QCd_final.txt
# 6386187
wc -l /vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/jamapsy_Giannakopoulou_2021_exclude_whi_23andMe_ukb.txt
# 7440923 originally

wc -l /vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/MDD_meta_AA_FIN_QC_half_QCd
# 5326425
wc -l /vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/MDD_meta_AA_FIN_QCd_final.txt
# 5326425
wc -l /vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/MDD_meta_AA_FIN_QC.txt
# 5326425 originally



# EOF ----