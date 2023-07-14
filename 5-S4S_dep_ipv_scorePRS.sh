#!/usr/bin/env bash

######################################################
# execute script in current directory
#$ -cwd
# want any .e/.o stuff to show up here too
#$ -e ./logs/
#$ -o ./logs/
# shell for qsub to use:
#$ -S /bin/bash
# label for the job; displayed by qstat
#$ -N ipv_scorePRS
######################################################


date

# combine all chromosome weights

cat /vcu_gpfs2/home/oharark/S4S/PRS/dep/output/S4S-EUR_dep*_META* > /vcu_gpfs2/home/oharark/S4S/PRS/dep/output/S4S-EUR_dep_META_pst_eff_a1_b0.5_phi1e-02_ALLchr.txt
cat /vcu_gpfs2/home/oharark/S4S/PRS/dep/output/S4S-SAS_dep*_META* > /vcu_gpfs2/home/oharark/S4S/PRS/dep/output/S4S-SAS_dep_META_pst_eff_a1_b0.5_phi1e-02_ALLchr.txt
cat /vcu_gpfs2/home/oharark/S4S/PRS/dep/output/S4S-EAS_dep*_META* > /vcu_gpfs2/home/oharark/S4S/PRS/dep/output/S4S-EAS_dep_META_pst_eff_a1_b0.5_phi1e-02_ALLchr.txt
cat /vcu_gpfs2/home/oharark/S4S/PRS/dep/output/S4S-AFR_dep*_META* > /vcu_gpfs2/home/oharark/S4S/PRS/dep/output/S4S-AFR_dep_META_pst_eff_a1_b0.5_phi1e-02_ALLchr.txt
cat /vcu_gpfs2/home/oharark/S4S/PRS/dep/output/S4S-AMR_dep*_META* > /vcu_gpfs2/home/oharark/S4S/PRS/dep/output/S4S-AMR_dep_META_pst_eff_a1_b0.5_phi1e-02_ALLchr.txt


# prep score file
# keep variant ID $2, ref allele $4, and score $6 columns only

cut -f 2,4,6 -d$'\t' /vcu_gpfs2/home/oharark/S4S/PRS/dep/output/S4S-AMR_dep_META_pst_eff_a1_b0.5_phi1e-02_ALLchr.txt > /vcu_gpfs2/home/oharark/S4S/PRS/dep/output/S4S-AMR_dep_META_pst_eff_a1_b0.5_phi1e-02_ALLchr.score
cut -f 2,4,6 -d$'\t' /vcu_gpfs2/home/oharark/S4S/PRS/dep/output/S4S-AFR_dep_META_pst_eff_a1_b0.5_phi1e-02_ALLchr.txt > /vcu_gpfs2/home/oharark/S4S/PRS/dep/output/S4S-AFR_dep_META_pst_eff_a1_b0.5_phi1e-02_ALLchr.score
cut -f 2,4,6 -d$'\t' /vcu_gpfs2/home/oharark/S4S/PRS/dep/output/S4S-EUR_dep_META_pst_eff_a1_b0.5_phi1e-02_ALLchr.txt > /vcu_gpfs2/home/oharark/S4S/PRS/dep/output/S4S-EUR_dep_META_pst_eff_a1_b0.5_phi1e-02_ALLchr.score
cut -f 2,4,6 -d$'\t' /vcu_gpfs2/home/oharark/S4S/PRS/dep/output/S4S-EAS_dep_META_pst_eff_a1_b0.5_phi1e-02_ALLchr.txt > /vcu_gpfs2/home/oharark/S4S/PRS/dep/output/S4S-EAS_dep_META_pst_eff_a1_b0.5_phi1e-02_ALLchr.score
cut -f 2,4,6 -d$'\t' /vcu_gpfs2/home/oharark/S4S/PRS/dep/output/S4S-SAS_dep_META_pst_eff_a1_b0.5_phi1e-02_ALLchr.txt > /vcu_gpfs2/home/oharark/S4S/PRS/dep/output/S4S-SAS_dep_META_pst_eff_a1_b0.5_phi1e-02_ALLchr.score


# score each ancestry ---

/vcu_gpfs2/home/oharark/Utils/plink2  --bfile /vcu_gpfs2/home/S4S/S4S_cohorts_15_work/S4S_cohort_15_AFR_info8_maf05_hwe6_geno025_mind025_final \
--bim /vcu_gpfs2/home/S4S/S4S_cohorts_15_work/S4S_cohort_15_AFR_info8_maf05_hwe6_geno025_mind025_final.bim --rm-dup exclude-all \
--score /vcu_gpfs2/home/oharark/S4S/PRS/dep/output/S4S-AFR_dep_META_pst_eff_a1_b0.5_phi1e-02_ALLchr.score list-variants \
--out /vcu_gpfs2/home/oharark/S4S/PRS/dep/output/S4S_AFR_dep_scores

/vcu_gpfs2/home/oharark/Utils/plink2  --bfile /vcu_gpfs2/home/S4S/S4S_cohorts_15_work/S4S_cohort_15_AMR_info8_maf05_hwe6_geno025_mind025_final \
--bim /vcu_gpfs2/home/S4S/S4S_cohorts_15_work/S4S_cohort_15_AMR_info8_maf05_hwe6_geno025_mind025_final.bim --rm-dup exclude-all \
--score /vcu_gpfs2/home/oharark/S4S/PRS/dep/output/S4S-AMR_dep_META_pst_eff_a1_b0.5_phi1e-02_ALLchr.score list-variants \
--out /vcu_gpfs2/home/oharark/S4S/PRS/dep/output/S4S_AMR_dep_scores

/vcu_gpfs2/home/oharark/Utils/plink2  --bfile /vcu_gpfs2/home/S4S/S4S_cohorts_15_work/S4S_cohort_15_EUR_info8_maf05_hwe6_geno025_mind025_final \
--bim /vcu_gpfs2/home/S4S/S4S_cohorts_15_work/S4S_cohort_15_EUR_info8_maf05_hwe6_geno025_mind025_final.bim --rm-dup exclude-all \
--score /vcu_gpfs2/home/oharark/S4S/PRS/dep/output/S4S-EUR_dep_META_pst_eff_a1_b0.5_phi1e-02_ALLchr.score list-variants \
--out /vcu_gpfs2/home/oharark/S4S/PRS/dep/output/S4S_EUR_dep_scores

/vcu_gpfs2/home/oharark/Utils/plink2  --bfile /vcu_gpfs2/home/S4S/S4S_cohorts_15_work/S4S_cohort_15_EAS_info8_maf05_hwe6_geno025_mind025_final \
--bim /vcu_gpfs2/home/S4S/S4S_cohorts_15_work/S4S_cohort_15_EAS_info8_maf05_hwe6_geno025_mind025_final.bim --rm-dup exclude-all \
--score /vcu_gpfs2/home/oharark/S4S/PRS/dep/output/S4S-EAS_dep_META_pst_eff_a1_b0.5_phi1e-02_ALLchr.score list-variants \
--out /vcu_gpfs2/home/oharark/S4S/PRS/dep/output/S4S_EAS_dep_scores

/vcu_gpfs2/home/oharark/Utils/plink2  --bfile /vcu_gpfs2/home/S4S/S4S_cohorts_15_work/S4S_cohort_15_SAS_info8_maf05_hwe6_geno025_mind025_final \
--bim /vcu_gpfs2/home/S4S/S4S_cohorts_15_work/S4S_cohort_15_SAS_info8_maf05_hwe6_geno025_mind025_final.bim --rm-dup exclude-all \
--score /vcu_gpfs2/home/oharark/S4S/PRS/dep/output/S4S-SAS_dep_META_pst_eff_a1_b0.5_phi1e-02_ALLchr.score list-variants \
--out /vcu_gpfs2/home/oharark/S4S/PRS/dep/output/S4S_SAS_dep_scores