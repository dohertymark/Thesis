# Thesis Chapter 3 
# Commands for Running Repeat Expansion Tools


##########################################
# Set Working Directory
##########################################
# Set The base folder to work in here
base_directory= 
# Set Directory To House Repeat Calls 
Repeat_Calls_Dir=${base_directory}/Chapter3/Repeat_Calls/

##########################################
# Locations of repeat files for each tool
##########################################

EH2_json_dir=~/programs/ExpansionHunter-v2.5.5-linux_x86_64/data/repeat-specs/hg19/
EH3_json_file=~/programs/ExpansionHunter-v3.2.2-linux_x86_64/variant_catalog/hg19/variant_catalog.json
GangSTR_beds=~/programs/GangSTR-2.4.4/target_files/
exSTRa_repeats=~/programs/exSTRa/repeat_database/hg19_extra_repeats
HipSTR_bed_location=~/programs/HipSTR/str_directory/
HipSTR_beds_autosome=~/programs/HipSTR/str_directory/chr1.22.hg19.hipstr_reference.bed
HipSTR_beds_X=~/programs/HipSTR/str_directory/chrX.hg19.hipstr_reference.bed
HipSTR_beds_Y=~/programs/HipSTR/str_directory/chrY.hg19.hipstr_reference.bed
STRetch_repeats=~/programs/STRetch/reference-data/hg19.simpleRepeat_period1-6_dedup.sorted.bed
exSTRa_repeats_targeted=~/programs/exSTRa/repeat_database/repeat_expansion_disorders_hg19.txt
Repeat_Seq_Repeats=~/programs/repeatseq/regions/RepeatSeq_Target_Repeats

##########################################
# Other Required Files & Locations 
##########################################

# Location of reference genome etc. 
resources=${base_directory}/Resources/Chapter3/Chapter3_Resources/
# File with the sex of each sample														
sex_file=${resources}/WGS_WES_ALL_STUDIES_ALL_SAMPLES_SEX
# File denoting the case / control status of our ALS WGS data 
ALS_WGS_cc=${resources}/pmid_vs_cc_status
# Bed file for SureSelect Target Regions (This was used for Irish ALS WES)
SureSelect_bed=${resources}/SureSelect_Targets.bed
# Bed file for MedExome Target Regions (This was used for RCSI WES)
MedExome_bed=${resources}/MedExome_hg19_capture_targets.bed
# Combined exome bed file (bedtools merge the ones above)
Combined_Exome_Targets_bed=${resources}/Combined_Exome_Targets.bed
# a bed file of the target disease beds and have added 20bp buffer to each region
Target_Disease_Beds=${resources}/All_Disease_Target_Beds.bed


##########################################
# Create Directories We Need 
##########################################

# We Have Five Different Sets of Data That Require Repeat Calls 
mkdir ${Repeat_Calls_Dir}/RCSI_WGS_PCR_FREE
mkdir ${Repeat_Calls_Dir}/RCSI_WGS_PCR
mkdir ${Repeat_Calls_Dir}/RCSI_WES
mkdir ${Repeat_Calls_Dir}/ALS_WGS_PCR_FREE
mkdir ${Repeat_Calls_Dir}/ALS_WES

for folder in "RCSI_WGS_PCR_FREE" "RCSI_WGS_PCR" "RCSI_WES" "ALS_WGS_PCR_FREE" "ALS_WES" ; do
	mkdir ${Repeat_Calls_Dir}/${folder}/EH_v2
	mkdir ${Repeat_Calls_Dir}/${folder}/EH_v2/Calls/
	mkdir ${Repeat_Calls_Dir}/${folder}/EH_v3
	mkdir ${Repeat_Calls_Dir}/${folder}/EH_v3/Calls/
	mkdir ${Repeat_Calls_Dir}/${folder}/ExpansionHunterDenovo
	mkdir ${Repeat_Calls_Dir}/${folder}/ExpansionHunterDenovo/Calls/
	mkdir ${Repeat_Calls_Dir}/${folder}/ExpansionHunterDenovo/Calls/str-profiles/
	mkdir ${Repeat_Calls_Dir}/${folder}/exSTRa/
	mkdir ${Repeat_Calls_Dir}/${folder}/GangSTR
 	mkdir ${Repeat_Calls_Dir}/${folder}/GangSTR/GangSTR_Target_Mode_With_Disease_Targets/
	mkdir ${Repeat_Calls_Dir}/${folder}/GangSTR/GangSTR_Target_Mode_With_Disease_Targets/Calls
	mkdir ${Repeat_Calls_Dir}/${folder}/GangSTR/GangSTR_NonTarget_Mode_With_Disease_Targets
	mkdir ${Repeat_Calls_Dir}/${folder}/GangSTR/GangSTR_NonTarget_Mode_With_Disease_Targets/Calls/
	mkdir ${Repeat_Calls_Dir}/${folder}/tredparse
	mkdir ${Repeat_Calls_Dir}/${folder}/tredparse/Calls/
	mkdir ${Repeat_Calls_Dir}/${folder}/RepeatSeq/
	mkdir ${Repeat_Calls_Dir}/${folder}/RepeatSeq/Calls/
	mkdir ${Repeat_Calls_Dir}/${folder}/STRetch/
	mkdir ${Repeat_Calls_Dir}/${folder}/STRetch/Calls/
    	mkdir ${Repeat_Calls_Dir}/${folder}/STRetch/Cases
    	mkdir ${Repeat_Calls_Dir}/${folder}/STRetch/Controls 
	mkdir ${Repeat_Calls_Dir}/${folder}/HipSTR/
	mkdir ${Repeat_Calls_Dir}/${folder}/HipSTR/Calls/
done

##########################################
# Call Each Tool 
##########################################

################################################################################
################################################################################
#########################    ExpansionHunter v2.5.5    ######################### 
################################################################################
################################################################################

# WGS Samples 
for folder in "ALS_WGS_PCR_FREE" "RCSI_WGS_PCR" "RCSI_WGS_PCR_FREE" ; do 
	data_folder=${base_directory}/Working_Data/${folder}/
	bams=${data_folder}/*bam
	for bam in ${bams} ; do 
		filename=`echo $( basename $bam) | sed 's/.bam//g'`
		sex=`more ${sex_file} | grep -w ${filename} | cut -f 2`
		if [ "$sex" == "F" ] ; then
			ExpansionHunter2 --bam ${sample} --ref-fasta ${resources}/canonical.hg19.fasta --repeat-spec ~/programs/ExpansionHunter-v2.5.5-linux_x86_64/data/repeat-specs/hg19/ --vcf ${Repeat_Calls_Dir}/${folder}/EH_v2/Calls/EH_v2_${filename}.vcf --json ${Repeat_Calls_Dir}/${folder}/EH_v2/Calls/EH_v2_${filename}.json --log ${Repeat_Calls_Dir}/${folder}/EH_v2/Calls/EH_v2_${filename}.log --sex female 
		else 
			ExpansionHunter2 --bam ${sample} --ref-fasta ${resources}/canonical.hg19.fasta --repeat-spec ~/programs/ExpansionHunter-v2.5.5-linux_x86_64/data/repeat-specs/hg19/ --vcf ${Repeat_Calls_Dir}/${folder}/EH_v2/Calls/EH_v2_${filename}.vcf --json ${Repeat_Calls_Dir}/${folder}/EH_v2/Calls/EH_v2_${filename}.json --log ${Repeat_Calls_Dir}/${folder}/EH_v2/Calls/EH_v2_${filename}.log --sex male
		fi
	done
done

# Exome Samples 
for folder in "ALS_WES" "RCSI_WES" ; do 
	data_folder=${base_directory}/Working_Data/${folder}/
	bams=${data_folder}/*bam
	for bam in ${bams} ; do 
		filename=`echo $( basename $bam) | sed 's/.bam//g'`
		sex=`more ${sex_file} | grep -w ${filename} | cut -f 2`
		coverage=`more ${coverage_directory}/${folder}/DOC_${filename}*
		if [ "$sex" == "M" ] ; then
			ExpansionHunter2 --bam ${bam} --ref-fasta ${resources}/canonical.hg19.fasta --repeat-spec ~/programs/ExpansionHunter-v2.5.5-linux_x86_64/data/repeat-specs/hg19/ --vcf ${Repeat_Calls_Dir}/${folder}/EH_v2/Calls/EH_v2_${filename}.vcf --json ${Repeat_Calls_Dir}/${folder}/EH_v2/Calls/EH_v2_${filename}.json --log ${Repeat_Calls_Dir}/${folder}/EH_v2/Calls/EH_v2_${filename}.log --read-depth $coverage --sex male
		else 
			ExpansionHunter2 --bam ${bam} --ref-fasta ${resources}/canonical.hg19.fasta --repeat-spec ~/programs/ExpansionHunter-v2.5.5-linux_x86_64/data/repeat-specs/hg19/ --vcf ${Repeat_Calls_Dir}/${folder}/EH_v2/Calls/EH_v2_${filename}.vcf --json ${Repeat_Calls_Dir}/${folder}/EH_v2/Calls/EH_v2_${filename}.json --log ${Repeat_Calls_Dir}/${folder}/EH_v2/Calls/EH_v2_${filename}.log --read-depth $coverage --sex female 
		fi
	done
done

################################################################################
################################################################################
#########################    ExpansionHunter v3.2.2    ######################### 
################################################################################
################################################################################

for folder in "ALS_WGS_PCR_FREE" "ALS_WES" "RCSI_WGS_PCR_FREE" "RCSI_WGS_PCR" "RCSI_WES" ; do 
	data_folder=${base_directory}/Working_Data/${folder}/
	bams=$data_folder/*bam
	for bam in $bams ; do
		filename=`echo $( basename $bam) | sed 's/.bam//g'`
		sex=`more ${sex_file} | grep -w ${filename} | cut -f 2`
		if [ "$sex" == "M" ] ; then
			ExpansionHunter --reads ${bam} --reference ${resources}/canonical.hg19.fasta --output-prefix ${Repeat_Calls_Dir}/${folder}/EH_v3/Calls/EH_v3_`echo $( basename $bam) | sed 's/.bam//g'` --variant-catalog ~/programs/ExpansionHunter-v3.2.2-linux_x86_64/variant_catalog/hg19/variant_catalog.json --sex male 
		else 
			ExpansionHunter --reads ${bam} --reference ${resources}/canonical.hg19.fasta --output-prefix ${Repeat_Calls_Dir}/${folder}/EH_v3/Calls/EH_v3_`echo $( basename $bam) | sed 's/.bam//g'` --variant-catalog ~/programs/ExpansionHunter-v3.2.2-linux_x86_64/variant_catalog/hg19/variant_catalog.json --sex female
		fi
	done
	####
	# Create REViewer files 
	####
	# ATXN8OS is a combination of different repeats (CTA(N)CTG(N))
	# EH gives the output for each repeat separately but doesn't phase the data so we can't reconstruct the entire repeat (don't know which CTA allele goes with which CTG allele)
	# We need to run another program (called REViewer) on the EH results to get a phased version of this information
	# First we need to sort and index the EH generated bam files 
	for bam in ${Repeat_Calls_Dir}/${folder}/EH_v3/Calls/EH_*bam ; do
		filename=`echo $( basename $bam) | sed 's/.bam//g'`
		# Sort file 
		samtools sort -l 9 -T ${Repeat_Calls_Dir}/${folder}/EH_v3/Calls/temp.${filename}.bam -@ 2 -o ${Repeat_Calls_Dir}/${folder}/EH_v3/Calls/sorted.${filename}.bam ${bam}
		# Index file 
		samtools index ${Repeat_Calls_Dir}/${folder}/EH_v3/Calls/sorted.${filename}.bam
		# Run REViewer to get the phased information 
		REViewer --reads ${Repeat_Calls_Dir}/${folder}/EH_v3/Calls/sorted.${filename}.bam --vcf `echo ${bam} | sed 's/.bam//g' | sed 's/_realigned//g'`.vcf --reference ${resources}/canonical.hg19.fasta --catalog ~/programs/ExpansionHunter-v3.2.2-linux_x86_64/variant_catalog/hg19/variant_catalog.json --locus ATXN8OS --output-prefix ${Repeat_Calls_Dir}/${folder}/EH_v3/Calls/REViewer_ATXN8OS_${filename} --output-phasing-info TRUE
	done
done

################################################################################
################################################################################
#########################    ExpansionHunter de novo   ######################### 
################################################################################
################################################################################


folder="ALS_WGS_PCR_FREE"
data_folder=${base_directory}/Working_Data/${folder}/
bams=${data_folder}/*bam


################################################################################
# Create STR profile for each sample 
################################################################################

for bam in ${bams} ; do 
	filename=`echo $( basename $bam) | sed 's/.bam//g'`
	ExpansionHunterDenovo profile --reads ${bam} --reference ${resources}/canonical.hg19.fasta --output-prefix ${Repeat_Calls_Dir}/${folder}/ExpansionHunterDenovo/Calls/str-profiles/${filename} --min-anchor-mapq 50 --max-irr-mapq 40
done

################################################################################
# Create Manifest for whole dataset 
################################################################################

for file in ${Repeat_Calls_Dir}/${folder}/ExpansionHunterDenovo/Calls/str-profiles/*json ; do 
	filename=`echo $( basename $file) | sed 's/.str_profile.json//g'`
	echo ${filename} >> ${Repeat_Calls_Dir}/${folder}/ExpansionHunterDenovo/Calls/temp1
	# case/control by case/control status 
	more ${ALS_WGS_cc} | grep ${filename} | cut -f 2 | sed 's/1/case/g' | sed 's/0/control/g' >> ${Repeat_Calls_Dir}/${folder}/ExpansionHunterDenovo/Calls/temp2
	echo "str-profiles/${filename}.str_profile.json" >> ${Repeat_Calls_Dir}/${folder}/ExpansionHunterDenovo/Calls/temp4
done

cd ${Repeat_Calls_Dir}/${folder}/ExpansionHunterDenovo/Calls
# case/control by case/control status
paste temp1 temp2 | paste - temp4 > manifest.case_control.tsv
# Doing a test removing all non positive ALS samples 
more $ALS_WGS_cc | awk '$5 ==0' | cut -f 1 > negative_cases
paste temp1 temp3 | paste - temp4 | grep -vf negative_cases > manifest.C9orf72_test.tsv

################################################################################
# Merge Files 
################################################################################

# case/control by case/control status
ExpansionHunterDenovo merge --reference ${resources}/canonical.hg19.fasta --manifest manifest.case_control.tsv --output-prefix ExpansionHunterDeNovo.case_control
# case/control by C9orf72 status (not actual case/control just to test)
ExpansionHunterDenovo merge --reference ${resources}/canonical.hg19.fasta --manifest manifest.C9orf72_test.tsv --output-prefix ExpansionHunterDeNovo.C9orf72_test

################################################################################
# Case / Control Analysis
################################################################################

# Test Each Locus
# case/control by case/control status
~/programs/ExpansionHunterDenovo/scripts/casecontrol.py locus --manifest manifest.case_control.tsv --multisample-profile ExpansionHunterDeNovo.case_control.multisample_profile.json --output ExpansionHunterDeNovo.case_control.locus.tsv 
# case/control by C9orf72 status (only using controls and C9 cases)
~/programs/ExpansionHunterDenovo/scripts/casecontrol.py locus --manifest manifest.C9orf72_test.tsv --multisample-profile ExpansionHunterDeNovo.C9orf72_test.multisample_profile.json --output ExpansionHunterDeNovo.C9orf72_test.locus.tsv 

################################################################################
################################################################################
#########################            exSTRa            ######################### 
################################################################################
################################################################################ 

# Not Running exSTRa for ALS_WES or RCSI_WGS_PCR_FREE as we don't have appropriate controls
for folder in "RCSI_WGS_PCR" "RCSI_WES" "ALS_WGS_PCR_FREE" ; do 
	data_folder=$RCSI_WGS_PCR_dir
	exSTRa_score.pl ${resources}/canonical.hg19.fasta ${exSTRa_repeats_targeted} ${data_folder}/*bam > ${Repeat_Calls_Dir}/${folder}/exSTRa/${folder}_exSTRa_scores_targeted
	Rscript exSTRa_p_values ${Repeat_Calls_Dir}/${folder}/exSTRa/ RCSI_WGS_PCR_exSTRa_scores_Targeted ${folder}_exSTRa_scores_targeted "regex to recognise controls"
done

################################################################################
################################################################################
#########################            GangSTR           ######################### 
################################################################################
################################################################################ 

# merge all target beds into one 
cat ${GangSTR_beds}/*bed | bedtools sort -faidx ${resources}/canonical.hg19.fasta.fai > ${GangSTR_beds}/all_repeats.bed

# PCR FREE WGS Data 
for folder in "ALS_WGS_PCR_FREE" "RCSI_WGS_PCR_FREE" ; do 
	data_folder=${base_directory}/Working_Data/${folder}/
	bams_csv=`echo ${data_folder}/*bam | sed 's/ /,/g'`
	files_csv=`echo $bams_csv | tr ',' '\n' | sed 's/^.*sorted_//g' | sed 's/.bam//g' | tr '\n' ',' | sed 's/,$//g'`
	# GangSTR_Target_Mode_With_Disease_Targets
	GangSTR --bam ${bams_csv} --bam-samps ${files_csv} --targeted --ref ${resources}/canonical.hg19.fasta --regions ${GangSTR_beds}/all_repeats.bed 	--output-readinfo --out ${Repeat_Calls_Dir}/${folder}/GangSTR/GangSTR_Target_Mode_With_Disease_Targets/Calls/GangSTR_Target_Mode_With_Disease_Targets  --include-ggl
	# GangSTR_NonTarget_Mode_With_Disease_Targets
	GangSTR --bam ${bams_csv} --bam-samps ${files_csv}            --ref ${resources}/canonical.hg19.fasta --regions ${GangSTR_beds}/all_repeats.bed 	--output-readinfo --out ${Repeat_Calls_Dir}/${folder}/GangSTR/GangSTR_NonTarget_Mode_With_Disease_Targets/Calls/GangSTR_NonTarget_Mode_With_Disease_Targets --include-ggl
done

# Add nonuniform flag for PCR WGS Data 
folder="RCSI_WGS_PCR" 
data_folder=${base_directory}/Working_Data/${folder}/
bams_csv=`echo ${data_folder}/*bam | sed 's/ /,/g'`
files_csv=`echo $bams_csv | tr ',' '\n' | sed 's/^.*sorted_//g' | sed 's/.bam//g' | tr '\n' ',' | sed 's/,$//g'`
# GangSTR_Target_Mode_With_Disease_Targets
GangSTR --bam ${bams_csv} --bam-samps ${files_csv} --nonuniform --targeted --ref ${resources}/canonical.hg19.fasta --regions ${GangSTR_beds}/all_repeats.bed 	--output-readinfo --out ${Repeat_Calls_Dir}/${folder}/GangSTR/GangSTR_Target_Mode_With_Disease_Targets/Calls/GangSTR_Target_Mode_With_Disease_Targets --include-ggl
# GangSTR_NonTarget_Mode_With_Disease_Targets
GangSTR --bam ${bams_csv} --bam-samps ${files_csv}  --nonuniform          --ref ${resources}/canonical.hg19.fasta --regions ${GangSTR_beds}/all_repeats.bed 	--output-readinfo --out ${Repeat_Calls_Dir}/${folder}/GangSTR/GangSTR_NonTarget_Mode_With_Disease_Targets/Calls/GangSTR_NonTarget_Mode_With_Disease_Targets --include-ggl


# Exome Data 
# For exome wide screen going to use average overall sample coverage - the neccessary files have previously been computed 

for folder in "ALS_WGS_PCR_FREE" "RCSI_WGS_PCR_FREE" ; do 
	data_folder=$RCSI_WES_dir
	folder="RCSI_WES"
	bams_csv=`echo ${data_folder}/*bam | sed 's/ /,/g'`
	data_folder=${base_directory}/Working_Data/${folder}/
	bams_csv=`echo ${data_folder}/*bam | sed 's/ /,/g'`
	files_csv=`echo $bams_csv | tr ',' '\n' | sed 's/^.*sorted_//g' | sed 's/.bam//g' | tr '\n' ',' | sed 's/,$//g'`
	coverage_csv=`cat ${coverage_directory}/${folder}/DOC_*  | tr '\n' ',' |  sed 's/,$//g'` 
	# There is a sample with coverage too low to compute fragment length and standrad deviation - will take the first 10000 lines from each bam file and get the fragment length and stdev from this 
	# Choosing fragments lengths above 0 (to avoid mates) and below 1000 to avoid misalignments 
	# Note : this can be avoided if all samples are above threshold
	for bam in ${data_folder}/*bam ; do 
		samtools view $bam | cut -f 9 | head -n 1000 | awk '$1 > 0' | awk '$1 <10000'>> ${Repeat_Calls_Dir}/${folder}/GangSTR/GangSTR_Target_Mode_With_Disease_Targets/fragment_lengths
	done 
	fragment_length_mean=`more ${Repeat_Calls_Dir}/${folder}/GangSTR/GangSTR_Target_Mode_With_Disease_Targets/fragment_lengths | awk '{ sum += $1 } END { print sum / NR }'`
	fragment_length_stdev=`more ${Repeat_Calls_Dir}/${folder}/GangSTR/GangSTR_Target_Mode_With_Disease_Targets/fragment_lengths | awk '{sum+=$1; sumsq+=$1*$1} END {print sqrt(sumsq/NR - (sum/NR)**2)}'`
	if [ $folder == "RCSI_WES" ] ; then 
		# This requires the specified fragment length
		# GangSTR_ non Target_Mode_With_Disease_Targets 
		GangSTR --bam ${bams_csv} --bam-samps ${files_csv} --nonuniform            --ref ${resources}/canonical.hg19.fasta --regions ${GangSTR_beds}/all_repeats.bed --coverage ${coverage_csv}	--output-readinfo --out ${Repeat_Calls_Dir}/${folder}/GangSTR/GangSTR_NonTarget_Mode_With_Disease_Targets/Calls/GangSTR_NonTarget_Mode_With_Disease_Targets --insertmean ${fragment_length_mean} --insertsdev ${fragment_length_stdev} --include-ggl
		# GangSTR_Target_Mode_With_Disease_Targets 
		GangSTR --bam ${bams_csv} --bam-samps ${files_csv} --nonuniform --targeted --ref ${resources}/canonical.hg19.fasta --regions ${GangSTR_beds}/all_repeats.bed --coverage ${coverage_csv}	--output-readinfo --out ${Repeat_Calls_Dir}/${folder}/GangSTR/GangSTR_Target_Mode_With_Disease_Targets/Calls/GangSTR_Target_Mode_With_Disease_Targets --insertmean ${fragment_length_mean} --insertsdev ${fragment_length_stdev} --include-ggl
	else 
		# GangSTR_Target_Mode_With_Disease_Targets 
		GangSTR --bam ${bams_csv} --bam-samps ${files_csv} --nonuniform --targeted --ref ${resources}/canonical.hg19.fasta --regions ${GangSTR_beds}/all_repeats.bed --coverage ${coverage_csv}	--output-readinfo --out ${Repeat_Calls_Dir}/${folder}/GangSTR/GangSTR_Target_Mode_With_Disease_Targets/Calls/GangSTR_Target_Mode_With_Disease_Targets  --include-ggl
		# GangSTR non Target mode with disease targets  
		GangSTR --bam ${bams_csv} --bam-samps ${files_csv} --nonuniform             --ref ${resources}/canonical.hg19.fasta --regions ${GangSTR_beds}/all_repeats.bed --coverage ${coverage_csv}	--output-readinfo --out ${Repeat_Calls_Dir}/${folder}/GangSTR/GangSTR_NonTarget_Mode_With_Disease_Targets/Calls/GangSTR_NonTarget_Mode_With_Disease_Targets --include-ggl
	fi
done

################################################################################
################################################################################
#########################            HipSTR            ######################### 
################################################################################
################################################################################

for folder in "RCSI_WGS_PCR_FREE" "ALS_WES" "RCSI_WGS_PCR" "ALS_WGS_PCR_FREE" "RCSI_WES" ; do 
	data_folder=${base_directory}/Working_Data/${folder}/
	bams=${data_folder}/*bam
	# Need to get a list of which samples are male and female 
	echo $bams | sed "s|$data_folder||g" | sed 's/\///g' | sed 's/.bam//g' | tr ' ' '\n'  > ${Repeat_Calls_Dir}/${folder}/HipSTR/all_samples_parsed
	ls $data_folder/*bam > ${Repeat_Calls_Dir}/${folder}/HipSTR/all_samples
	more $sex_file | grep -wf ${Repeat_Calls_Dir}/${folder}/HipSTR/all_samples | grep -w "M"  | cut -f 1 > ${Repeat_Calls_Dir}/${folder}/HipSTR/temp_males_1
	more ${Repeat_Calls_Dir}/${folder}/HipSTR/all_samples | grep -f ${Repeat_Calls_Dir}/${folder}/HipSTR/temp_males_1 > ${Repeat_Calls_Dir}/${folder}/HipSTR/male_samples
	more $sex_file | grep -wf ${Repeat_Calls_Dir}/${folder}/HipSTR/all_samples | grep -w "F" | cut -f 1 > ${Repeat_Calls_Dir}/${folder}/HipSTR/temp_females_1
	more ${Repeat_Calls_Dir}/${folder}/HipSTR/all_samples | grep -f ${Repeat_Calls_Dir}/${folder}/HipSTR/temp_females_1 > ${Repeat_Calls_Dir}/${folder}/HipSTR/females_samples
	for chrom in `seq 1 22` ; do 
		$HipSTR --bam-files ${Repeat_Calls_Dir}/${folder}/HipSTR/all_samples --fasta ${resources}/canonical.hg19.fasta --regions ${HipSTR_bed_location}/disease.autosome.hipstr_reference.bed --str-vcf ${Repeat_Calls_Dir}/${folder}/HipSTR/chr${chrom}_HipSTR_calls.vcf.gz --log ${Repeat_Calls_Dir}/${folder}/HipSTR/Disease_Genes/chr${chrom}.log.txt --min-reads 5 --output-filters --chrom chr${chrom}  & 
	done 
	# Genotype chr X in males 
	$HipSTR --bam-files ${Repeat_Calls_Dir}/${folder}/HipSTR/male_samples --fasta ${resources}/canonical.hg19.fasta --regions ${HipSTR_bed_location}/disease.chrX.hipstr_reference.bed --haploid-chrs chrX --str-vcf ${Repeat_Calls_Dir}/${folder}/HipSTR/chr.males_X_HipSTR_calls.vcf.gz --log ${Repeat_Calls_Dir}/${folder}/HipSTR/Disease_Genes/chr.males.xlog.txt --min-reads 5 --chrom chrX --output-filters &
	# Genotype chrX in females 
	$HipSTR --bam-files ${Repeat_Calls_Dir}/${folder}/HipSTR/females_samples --fasta ${resources}/canonical.hg19.fasta --regions ${HipSTR_bed_location}/disease.chrX.hipstr_reference.bed --str-vcf ${Repeat_Calls_Dir}/${folder}/HipSTR/chr.females_X_HipSTR_calls.vcf.gz --log ${Repeat_Calls_Dir}/${folder}/HipSTR/Disease_Genes/chr.females.x.log.txt --min-reads 5 --chrom chrX --output-filters &
done

################################################################################
################################################################################
#########################           RepeatSeq          ######################### 
################################################################################
################################################################################

for folder in "ALS_WGS_PCR_FREE" "RCSI_WGS_PCR_FREE" "RCSI_WGS_PCR" "RCSI_WES" "ALS_WES" ; do 
	cd ${Repeat_Calls_Dir}/${folder}/RepeatSeq/Calls/
	data_folder=${base_directory}/Working_Data/${folder}/
	for bam in ${data_folder}/*bam ; do 
		repeatseq -calls -repeatseq  ${bam} ${resources}/canonical.hg19.fasta  ${Repeat_Seq_Repeats} & 		
	done
done

################################################################################
################################################################################
#########################            STRetch           ######################### 
################################################################################
################################################################################

# Step 1 : Calculate STRetch Files for All Samples
# WGS 
# Not running this on RCSI WGS PCR FREE as there are no appropriate controls
for folder in "ALS_WGS_PCR_FREE" "RCSI_WGS_PCR" ; do 
	data_folder=${base_directory}/Working_Data/${folder}/
	cd ${Repeat_Calls_Dir}/${folder}/STRetch/Calls/
	bams=${data_folder}/*bam
	~/programs/bpipe-0.9.9.9/bin/bpipe run -p input_regions=~programs/STRetch/reference-data/hg19.simpleRepeat_period1-6_dedup.sorted.bed ~/programs/STRetch/pipelines/STRetch_wgs_bam_pipeline.groovy  ${bams}
done

# Exome 
for folder in "RCSI_WES"  ; do 
	data_folder=${base_directory}/Working_Data/${folder}/
	bams=${data_folder}/*bam
	cd ${Repeat_Calls_Dir}/${folder}/STRetch/Calls/
	~/programs/bpipe-0.9.9.9/bin/bpipe run -p input_regions=~/programs/STRetch/reference-data/hg19.simpleRepeat_period1-6_dedup.sorted.bed ~/programs/STRetch/pipelines/STRetch_exome_bam_pipeline.groovy ${bams}
done

# Step 2 : Give each dataset appropriate controls
for folder in "ALS_WGS_PCR_FREE" "RCSI_WGS_PCR" "RCSI_WES" ; do 
    # For ALS Treating Controls as Controls (because controls are controls)
	# For RCSI treating Parents as Proxy for Controls 
	# Move all control files from ${Repeat_Calls_Dir}/${folder}/STRetch/Calls to ${Repeat_Calls_Dir}/${folder}/STRetch/Controls
	# Move all case files from ${Repeat_Calls_Dir}/${folder}/STRetch/Calls to ${Repeat_Calls_Dir}/${folder}/STRetch/Cases
	cd ${Repeat_Calls_Dir}/${folder}/STRetch/Controls
	# Creat control database from this 
	~/programs/STRetch/tools/bin/python ~/programs/STRetch/scripts/estimateSTR.py --model ~/programs/STRetch/scripts/STRcov.model.csv --locus_counts *.locus_counts --STR_counts *.STR_counts --median_cov *.median_cov --emit ${Repeat_Calls_Dir}/${folder}/STRetch/Controls/${folder}_Controls.tsv
	# In pipeline_config.groovy replace control file with ${Repeat_Calls_Dir}/${folder}/STRetch/Controls/${folder}_Controls.tsv
	# Now check cases against these controls
	cd ${Repeat_Calls_Dir}/${folder}/STRetch/Cases
	~/programs/STRetch/tools/bin/python ~/programs/STRetch/scripts/estimateSTR.py --model ~/programs/STRetch/scripts/STRcov.model.csv --locus_counts *.locus_counts --STR_counts *.STR_counts --median_cov *.median_cov
done


################################################################################
################################################################################
#########################           TREDPARSE          ######################### 
################################################################################
################################################################################

for folder in "ALS_WGS_PCR_FREE" "RCSI_WGS_PCR_FREE" "RCSI_WGS_PCR" "RCSI_WES" "ALS_WES" ; do 
	data_folder=${base_directory}/Working_Data/${folder}/
	for bam in ${data_folder}/*bam ; do 
		filename=`echo $( basename $bam) | sed 's/.bam//g'`
		sex=`more ${sex_file} | grep -w ${filename} | cut -f 2`
		if [ "$sex" == "M" ] ; then 
			$tredparse $sample --sample_id ${filename} --workdir ${Repeat_Calls_Dir}/${folder}/tredparse/Calls/ --ref hg19 --haploid chrX --useclippedreads --fullsearch 
		else
			$tredparse $sample --sample_id ${filename} --workdir ${Repeat_Calls_Dir}/${folder}/tredparse/Calls/ --ref hg19 --useclippedreads --fullsearch 
		fi
	done
done
