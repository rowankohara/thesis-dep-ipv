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
#$ -N ipv_mkbim
######################################################

# combine the bims into one file for all ancestries
# then break up by chromosome to be able to parallelize prscsx


echo -e "/vcu_gpfs2/home/S4S/S4S_cohorts_15_work/S4S_cohort_15_EAS_info8_maf05_hwe6_geno025_mind025_final\n\
/vcu_gpfs2/home/S4S/S4S_cohorts_15_work/S4S_cohort_15_AMR_info8_maf05_hwe6_geno025_mind025_final\n\
/vcu_gpfs2/home/S4S/S4S_cohorts_15_work/S4S_cohort_15_EUR_info8_maf05_hwe6_geno025_mind025_final\n\
/vcu_gpfs2/home/S4S/S4S_cohorts_15_work/S4S_cohort_15_SAS_info8_maf05_hwe6_geno025_mind025_final\n\
/vcu_gpfs2/home/S4S/S4S_cohorts_15_work/S4S_cohort_15_AFR_info8_maf05_hwe6_geno025_mind025_final" > \
/vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/mergelist.txt



for i in {1..22}
do
/vcu_gpfs2/home/oharark/Utils/plink --merge-list /vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/mergelist.txt \
--chr $i \
--make-just-bim \
--out /vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/S4S_cohort_15_All_info8_maf05_hwe6_geno025_mind025_final_chr$i
done


#EOF