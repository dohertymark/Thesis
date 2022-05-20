#!/bin/bash
###########################################################################################################
###########################################################################################################

# email: dohertm7@tcd.ie

# This is a script to perform two EPACTS SKAT tests on hard filtered SNPs and INDELs
# TO RUN:
#		Copy the SNPs and INDELs files to the directory you are working in and define the variables below 
#		Run this as Cuban_EPACTS.sh
# OUTPUT:
#		Two files:
#					1. ${project}.chr_rem.sites.anno.LOF.grp
#					2. ${project}.chr_rem.sites.anno.snps.grp
# AFTER:
# 		These files are analysed in Chapter4_Results_Plots.R

###########################################################################################################
###########################################################################################################
### DEFINE VARIABLES AND PATHS 
###########################################################################################################
###########################################################################################################

# Copy VCF to Current Dir - basefile should be the name WITHOUT vcf at the end 
snps_file="Cuban_results_hard_filtered_snps.tidy.vcf"
indels_file="Cuban_results_hard_filtered_indels.tidy.vcf"
ped="/path/to/cuban_pedigree.ped"
project="Cubans_ALS"

###########################################################################################################
###########################################################################################################
# Merge VCF Files 
###########################################################################################################
###########################################################################################################

more ${snps_file}| bgzip -c >  ${snps_file}.gz
tabix -p vcf ${snps_file}.gz
more ${indels_file}| bgzip -c >  ${indels_file}.gz
tabix -p vcf ${indels_file}.gz

bcftools concat -a ${snps_file}.gz ${indels_file}.gz  -o ${project}.vcf

# Remove chr to line up with EPACTS references 
more ${project}.vcf | sed 's/chr//g' | bgzip -c > ${project}.chr_rem.vcf.gz

tabix -p vcf ${project}.chr_rem.vcf.gz

###########################################################################################################
###########################################################################################################
# Annotate File
###########################################################################################################
###########################################################################################################

# We are going to annotate all of our variants with gene, functinoal impact etc. - first lets pull out each variant 
gunzip -c ${project}.chr_rem.vcf.gz | cut -f 1-8 | bgzip -c > ${project}.chr_rem.sites.vcf.gz

# Annotate all sites 
$epacts anno --in ${project}.chr_rem.sites.vcf.gz --out ${project}.chr_rem.sites.anno.vcf.gz 

###########################################################################################################
###########################################################################################################
# Analysis 1 : Missense Variants 
###########################################################################################################
###########################################################################################################

# Collapse to list of variant present grouped in each gene
$epacts make-group --vcf ${project}.chr_rem.sites.anno.vcf.gz --out ${project}.chr_rem.sites.anno.snps.grp -format epacts --type Missense --type Nonsynonymous  -pass
$epacts group --vcf ${project}.chr_rem.vcf.gz --ped sample.ped --max-maf 0.05 --groupf ${project}.chr_rem.sites.anno.snps.grp --pheno DISEASE --test skat --out results.geno.missense.skat --skat-o --run 5

###########################################################################################################
###########################################################################################################
# Analysis 2 : LOF Variants 
###########################################################################################################
###########################################################################################################

# Collapse to list of variant present grouped in each gene
$epacts make-group --vcf ${project}.chr_rem.sites.anno.vcf.gz --out ${project}.chr_rem.sites.anno.LOF.grp -format epacts --type StructuralVariation  --type Stop_Gain  --type Stop_Loss  --type Start_Gain  --type Start_Loss  --type Frameshift  --type CodonGain  --type CodonLoss  --type CodonRegion  --type Insertion  --type Deletion  --type Essential_Splice_Site  --type Nonsense  -pass
$epacts group --vcf ${project}.chr_rem.vcf.gz --ped sample.ped --max-maf 0.05 --groupf ${project}.chr_rem.sites.anno.LOF.grp --pheno DISEASE --test skat --out results.geno.LOF.skat --skat-o --run 5
###########################################################################################################
###########################################################################################################
exit 0
