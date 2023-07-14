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
#$ -N ipv_mkPRS
######################################################

# Parallel computation for the 22 autosomes is recommended.
# Default is iterating through 22 autosomes (can be time-consuming).

# EAS ---

for i in {1..22}
do

python /vcu_gpfs2/home/oharark/S4S/PRS/PRScsx/PRScsx.py \
 --ref_dir=/vcu_gpfs2/home/oharark/S4S/PRS/refs \
 --bim_prefix=/vcu_gpfs2/home/S4S/S4S_cohorts_15_work/S4S_cohort_15_EAS_info8_maf05_hwe6_geno025_mind025_final \
 --sst_file=/vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/MDD_meta_AA_FIN_QCd_final.txt,/vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/MDD.AFR.MVP.NatNeuro2021_QCd_final.txt,/vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/jamapsy_Giannakopoulou_2021_exclude_whi_23andMe_ukb_QCd_final.txt \
 --n_gwas=1307228,59600,98502 \
 --pop=EUR,AFR,EAS \
 --chrom=$i \
 --phi=1e-2 \
 --meta=True \
 --out_dir=/vcu_gpfs2/home/oharark/S4S/PRS/dep/output \
 --out_name=S4S-EAS_dep$i

 
done

date



# EUR ---


for i in {1..22}
do

python /vcu_gpfs2/home/oharark/S4S/PRS/PRScsx/PRScsx.py \
 --ref_dir=/vcu_gpfs2/home/oharark/S4S/PRS/refs \
 --bim_prefix=/vcu_gpfs2/home/S4S/S4S_cohorts_15_work/S4S_cohort_15_EUR_info8_maf05_hwe6_geno025_mind025_final \
 --sst_file=/vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/MDD_meta_AA_FIN_QCd_final.txt,/vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/MDD.AFR.MVP.NatNeuro2021_QCd_final.txt,/vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/jamapsy_Giannakopoulou_2021_exclude_whi_23andMe_ukb_QCd_final.txt \
 --n_gwas=1307228,59600,98502 \
 --pop=EUR,AFR,EAS \
 --chrom=$i \
 --phi=1e-2 \
 --meta=True \
 --out_dir=/vcu_gpfs2/home/oharark/S4S/PRS/dep/output \
 --out_name=S4S-EUR_dep$i

 
done

date



# AFR ---


for i in {1..22}
do

python /vcu_gpfs2/home/oharark/S4S/PRS/PRScsx/PRScsx.py \
 --ref_dir=/vcu_gpfs2/home/oharark/S4S/PRS/refs \
 --bim_prefix=/vcu_gpfs2/home/S4S/S4S_cohorts_15_work/S4S_cohort_15_AFR_info8_maf05_hwe6_geno025_mind025_final \
 --sst_file=/vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/MDD_meta_AA_FIN_QCd_final.txt,/vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/MDD.AFR.MVP.NatNeuro2021_QCd_final.txt,/vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/jamapsy_Giannakopoulou_2021_exclude_whi_23andMe_ukb_QCd_final.txt \
 --n_gwas=1307228,59600,98502 \
 --pop=EUR,AFR,EAS \
 --chrom=$i \
 --phi=1e-2 \
 --meta=True \
 --out_dir=/vcu_gpfs2/home/oharark/S4S/PRS/dep/output \
 --out_name=S4S-AFR_dep$i

 
done

date



# SAS ---


for i in {1..22}
do

python /vcu_gpfs2/home/oharark/S4S/PRS/PRScsx/PRScsx.py \
 --ref_dir=/vcu_gpfs2/home/oharark/S4S/PRS/refs \
 --bim_prefix=/vcu_gpfs2/home/S4S/S4S_cohorts_15_work/S4S_cohort_15_SAS_info8_maf05_hwe6_geno025_mind025_final \
 --sst_file=/vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/MDD_meta_AA_FIN_QCd_final.txt,/vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/MDD.AFR.MVP.NatNeuro2021_QCd_final.txt,/vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/jamapsy_Giannakopoulou_2021_exclude_whi_23andMe_ukb_QCd_final.txt \
 --n_gwas=1307228,59600,98502 \
 --pop=EUR,AFR,EAS \
 --chrom=$i \
 --phi=1e-2 \
 --meta=True \
 --out_dir=/vcu_gpfs2/home/oharark/S4S/PRS/dep/output \
 --out_name=S4S-SAS_dep$i

 
done

date



# AMR ---


for i in {1..22}
do

python /vcu_gpfs2/home/oharark/S4S/PRS/PRScsx/PRScsx.py \
 --ref_dir=/vcu_gpfs2/home/oharark/S4S/PRS/refs \
 --bim_prefix=/vcu_gpfs2/home/S4S/S4S_cohorts_15_work/S4S_cohort_15_AMR_info8_maf05_hwe6_geno025_mind025_final \
 --sst_file=/vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/MDD_meta_AA_FIN_QCd_final.txt,/vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/MDD.AFR.MVP.NatNeuro2021_QCd_final.txt,/vcu_gpfs2/home/oharark/S4S/PRS/GWAS_SumStats/jamapsy_Giannakopoulou_2021_exclude_whi_23andMe_ukb_QCd_final.txt \
 --n_gwas=1307228,59600,98502 \
 --pop=EUR,AFR,EAS \
 --chrom=$i \
 --phi=1e-2 \
 --meta=True \
 --out_dir=/vcu_gpfs2/home/oharark/S4S/PRS/dep/output \
 --out_name=S4S-AMR_dep$i

 
done

date