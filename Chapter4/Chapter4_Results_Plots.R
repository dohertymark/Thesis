# Chapter 4: Generate Series of Figures and Some Analysis for Results 
# To Run: 
#    	1. Ensure Directories are set correctly 
#		2. Make sure you are in the main thesis directory i.e. /path/to/Thesis
#		3. Rscript /path/to/Chapter_2_Results_Plots.R 
# Figures Generated:
#		1. Reanalysis of journALS AOO including new variants 
# 		2. Proportion of cases explained in Cuba, North Am. and South Am. 
# 		3. A page of figures showing AAFs in cases and controls, EPACTS results, and testing for oligogenic results
# 		4. ATXN2 Plots and comparisons of different tools 
# email: dohertm7@tcd.ie

################################################################################
################################################################################

# Set Directory 

################################################################################
################################################################################

setwd("/path/to/Thesis")
output_directory=""


################################################################################
################################################################################

# Clear Memory

################################################################################
################################################################################

rm(list = ls())

################################################################################
################################################################################

#Load Packages

################################################################################
################################################################################

library(data.table)
library(rcompanion)
library(tidyr)
library(binom)
library(berryFunctions)
library(stringr)
library(kinship2)



################################################################################
################################################################################

# Read in Data 

################################################################################
################################################################################

#####
# Cuban Individuals Data 
#####

cuban_population_studies      <<- fread("/path/to/cuban_population_studies.txt",header=TRUE,sep='\t',na.strings=c("","NA"))
cuban_population_people_dataset <<- fread("/path/to/cuban_population_people_dataset.txt",header=TRUE,sep='\t',na.strings=c("","NA"))

#####
# journALS Data 
#####

# Read in journALS Population Studies 
population_studies      <<- fread("/path/to/journALS_population_studies_long_format.tsv",header=TRUE,sep='\t',na.strings=c("","NA"))
# Read in journALS Data 
lit_review_pen<- fread("/path/to/journALS_post_analysis.tsv.gz",header=TRUE,sep='\t',na.strings=c("","NA"))
# Age P-Value Cutoff
age_p_cutoff <<-  data.frame(fread("/path/to/journALS_p_cutoff_df.tsv",header=TRUE,sep='\t',na.strings=c("","NA")))
age_p_cutoff <<- age_p_cutoff$p_cutoff
# All ages df
all_ages_df <<- data.frame(fread("/path/to/journALS_all_ages_df.tsv",header=TRUE,sep='\t',na.strings=c("","NA")))
# These are needed for regions analysis plot
# Find list of pathogenic variants
pathogenic_variants<<-lit_review_pen$HGVS[lit_review_pen$acmg_literal=="Pathogenic" & !is.na(lit_review_pen$acmg_literal)]
# Add Cuban Variant 
pathogenic_variants<<-c(pathogenic_variants,"FUS:c.1512_1513delAG(p.[G505fs])")
#
# Find list of pathogenic or likely variants
p_or_lp_variants<<-lit_review_pen$HGVS[(lit_review_pen$acmg_literal=="Pathogenic" | lit_review_pen$acmg_literal=="Likely Pathogenic") & !is.na(lit_review_pen$acmg_literal)]
# Add Cuban Variant
p_or_lp_variants<<-c(p_or_lp_variants,"FUS:c.1512_1513delAG(p.[G505fs])")
# Reported variants in P or LP genes 
p_or_lp_genes<<-unique(lit_review_pen$gene[(lit_review_pen$acmg_literal=="Pathogenic" | lit_review_pen$acmg_literal=="Likely Pathogenic") & !is.na(lit_review_pen$acmg_literal)])
variants_in_p_or_lp_genes<<-unique(c(lit_review_pen$HGVS[lit_review_pen$population_carriers_count>0 & !is.na(lit_review_pen$population_carriers_count) & lit_review_pen$gene %in% p_or_lp_genes & !is.na(lit_review_pen$gene)],cuban_population_people_dataset$HGVS))
list_of_pathogenic_variants<<-list(pathogenic_variants,p_or_lp_variants,variants_in_p_or_lp_genes)

#####
# ATXN2 Data 
#####

ATXN2_Coverage <- data.frame(fread("/path/to/ATXN2_Coverage_Report",header=TRUE,sep='\t',na.strings=c("","NA")))
HipSTR <- data.frame(fread("/path/to/ATXN2_Genotyping/HipSTR/HipSTR_output",header=TRUE,sep='\t',na.strings=c("","NA")))
colnames(HipSTR) <- c("SAMPLE_ID","ATXN2_a1","ATXN2_a2")
tredparse <- data.frame(fread("/path/to/ATXN2_Genotyping/tredparse/tredparse_output",header=TRUE,sep='\t',na.strings=c("","NA")))
colnames(tredparse) <- c("SAMPLE_ID","ATXN2_a1","ATXN2_a2")
cc <- data.frame(fread("/path/to/cubans.ped",header=FALSE,sep='\t',na.strings=c("","NA")))
cc <- subset(cc,select=c("V2","V6"))
colnames(cc) <- c("SAMPLE_ID","cc_status")
cc$cc_status <- gsub("2","Case",cc$cc_status)
cc$cc_status <- gsub("1","Control",cc$cc_status)


#####
# Data for AAF Plot 
#####

# Need to Run These Commands First 

#gemini query -q "SELECT chrom,start,ref,alt,aaf_gnomad_all,gene,impact from variants where filter is NULL" Cubans_hard_filtered_snps.db > snps_aafs
#gemini query -q "SELECT chrom,start,ref,alt,aaf_gnomad_all,gene,impact from variants where filter is NULL" Cubans_hard_filtered_indels.db >> indels_aafs
#gemini query -q "SELECT chrom,start,ref,alt,aaf_gnomad_all,gene,impact from variants where filter is NULL" --show-samples Cubans_hard_filtered_snps.db > snps_aafs_with_cases
#gemini query -q "SELECT chrom,start,ref,alt,aaf_gnomad_all,gene,impact from variants where filter is NULL" --show-samples Cubans_hard_filtered_indels.db >> indels_aafs_with_cases

snps_aafs <- fread("/path/to/Gemini_Analysis/snps_aafs",header=FALSE,sep='\t',na.strings=c("","NA"))
indels_aafs <- fread("/path/to/Gemini_Analysis/indels_aafs",header=FALSE,sep='\t',na.strings=c("","NA"))
snps_aafs_with_cases <- fread("/path/to/Gemini_Analysis/snps_aafs_with_cases",header=FALSE,sep='\t',na.strings=c("","NA"))
indels_aafs_with_cases <- fread("/path/to/Gemini_Analysis/indels_aafs_with_cases",header=FALSE,sep='\t',na.strings=c("","NA"))
variant_sample_count.indels.cases <- fread("/path/to/Gemini_Analysis/variant_sample_count.indels.cases",header=TRUE,sep='\t',na.strings=c("","NA"))
variant_sample_count.indels.controls <- fread("/path/to/Gemini_Analysis/variant_sample_count.indels.controls",header=TRUE,sep='\t',na.strings=c("","NA"))
variant_sample_count.snps.controls <- fread("/path/to/Gemini_Analysis/variant_sample_count.snps.controls",header=TRUE,sep='\t',na.strings=c("","NA"))
variant_sample_count.snps.cases <- fread("/path/to/Gemini_Analysis/variant_sample_count.snps.cases",header=TRUE,sep='\t',na.strings=c("","NA"))

#####
# Coverage Data For Oligogenic Filtering 
#####

coverage_plate_1_1 <- read.csv("/path/to/PLATE_1/bams/PLATE_1_coverage",header=FALSE,sep='\t')
coverage_plate_2_1 <- read.csv("/path/to/PLATE_2/bams/PLATE_2_coverage",header=FALSE,sep='\t')
coverage_plate_2_2 <- read.csv("/path/to/PLATE_2_2/bams/PLATE_2_2_coverage",header=FALSE,sep='\t')
coverage_plate_3_1 <- read.csv("/path/to/PLATE_3_1/bams/PLATE_3_1_coverage",header=FALSE,sep='\t')
coverage_plate_4_1 <- read.csv("/path/to/PLATE_4/bams/PLATE_4_coverage",header=FALSE,sep='\t')

#####
# EPACTS burden testing results 
#####

epacts_missense <- fread("/path/to/EPACTS/results.geno.missense.skat.epacts",header=TRUE,sep='\t',na.strings=c("","NA"))
epacts_lof <- fread("/path/to/EPACTS/results.geno.LOF.skat.epacts",header=TRUE,sep='\t',na.strings=c("","NA"))

#####
# Pedigree
#####

ped <- read.table("/path/to/Pedigree_2302.ped.txt", header=T, fileEncoding="UTF-8")

#####
# VQSR
#####

#
# Create files 
#
# Before Opening these files need to cd into Cuban_PLS_Discordant_Exomes/Variant_Calling/ and run the following:
#more Exomes_recalibrate_SNP_plots.R | grep data..- | cut -d' ' -f 3- | sed 's/c(//g'  | sed 's/)//g' | head -n 1 | tail -n 1 > VQSR.MQ.QD.vector.txt 
#more Exomes_recalibrate_SNP_plots.R | grep data..- | cut -d' ' -f 3- | sed 's/c(//g'  | sed 's/)//g' | head -n 2 | tail -n 1 > VQSR.MQ.MQRankSum.vector.txt
#more Exomes_recalibrate_SNP_plots.R | grep data..- | cut -d' ' -f 3- | sed 's/c(//g'  | sed 's/)//g' | head -n 3 | tail -n 1 > VQSR.MQ.SOR.vector.txt
#more Exomes_recalibrate_SNP_plots.R | grep data..- | cut -d' ' -f 3- | sed 's/c(//g'  | sed 's/)//g' | head -n 4 | tail -n 1 > VQSR.MQ.FS.vector.txt
#more Exomes_recalibrate_SNP_plots.R | grep data..- | cut -d' ' -f 3- | sed 's/c(//g'  | sed 's/)//g' | head -n 5 | tail -n 1 > VQSR.MQ.ReadPosRankSum.vector.txt

MQ.FS <- scan("/path/to/Cuban_PLS_Discordant_Exomes/Variant_Calling/VQSR.MQ.FS.vector.txt",sep=",")
MQ.MQRankSum <- scan("/path/to/Cuban_PLS_Discordant_Exomes/Variant_Calling/VQSR.MQ.MQRankSum.vector.txt",sep=",")
MQ.QD <- scan("/path/to/Cuban_PLS_Discordant_Exomes/Variant_Calling/VQSR.MQ.QD.vector.txt",sep=",")
MQ.ReadPosRankSum <- scan("/path/to/Cuban_PLS_Discordant_Exomes/Variant_Calling/VQSR.MQ.ReadPosRankSum.vector.txt",sep=",")
MQ.SOR <- scan("/path/to/Cuban_PLS_Discordant_Exomes/Variant_Calling/VQSR.MQ.SOR.vector.txt",sep=",")

#####
# Exomes Duplication and Coverage Rate 
#####

duplication_rate <- read.csv("/path/to/Cuban_PLS_Discordant_Exomes/Exomes_duplication_rate",header=F,sep='\t')
coverage_rate <- read.csv("/path/to/Cuban_PLS_Discordant_Exomes/Variant_Calling/all_samples_mean_coverage",header=F,sep='\t')


################################################################################
################################################################################

# Define Useful Colours

################################################################################
################################################################################

# Define Colours - Not neccessaily all used here 
nicegrey        <<- "grey30"
niceblack       <<- "#3b444b"
nicedarkblue    <<- "#146D97"
nicered         <<- "#C02906"
niceorange      <<- "#FF7400"
nicelightblue   <<- "#1E96B9"
reallydarkblue  <<- "#01184E"
reallylightblue <<- "#CFF6F5"
darkgreen       <<- "#CFAA28"
lightgreen      <<- "#EDDC58"
pink<<- "#f4ddf1"
green<<-"#4897a5"
lime<<-"#c1ef37"
medium_green<<-"#5bcfa5"
brown<<-"#7f4827"
#bright_yellow<<-"#f4ff00"
bright_yellow<<-"#f5db37"
lavendar<<-"#8da3d3"
algae<<-"#526708"
yellow1<<-"#cad874"
sunset<<-"#c24802"
purple<<-"#523386"
peach<<-"#dca995"
light_purple<<-"#624f94"
pink_purple<<-"#9783dc"
hull<<-"#CAA575"
peach2<<-"#d45150"

# Colour fades
darkgreen_rgb_fade    <<- rgb(0.8117647,0.6666667,0.1568627,alpha=0.2)
nicegrey_rgb_fade     <<- rgb(0.3019608,0.3019608,0.3019608,alpha=0.2)
nicered_rgb_fade      <<- rgb(0.7529412, 0.1607843, 0.02352941,alpha=0.2)
nicedarkblue_rgb_fade <<- rgb(0.07843137,0.427451,0.5921569,alpha=0.2)
peach_rgb_fade <<- rgb(0.863,0.663,0.584,alpha=0.2)
reallydarkblue_rgb_fade <<- rgb(0.004,0.094,0.306,alpha=0.2)
niceorange_rgb_fade <<- rgb(1,0.455,0,alpha=0.2)


vangogh_palette<<-c(nicered,hull,peach,sunset,niceorange,lightgreen,darkgreen,reallylightblue,nicelightblue,nicedarkblue,reallydarkblue,pink,green,lime,brown,bright_yellow,lavendar,
algae,yellow1,purple,light_purple,pink_purple)

################################################################################
################################################################################

# Plot 1 - New AOO Plots for FUS Variants 

################################################################################
################################################################################


#####
# Define Functions 
#####

# Select What to Plot in Legend 
selections_legend=function(phenotype,sex,family_history,height){
	if (sex=="F"){
		sex <- "Female"
	}
	if (sex=="M"){
		sex <- "Male"
	}
	if (family_history=="F"){
		family_history <- "Familial"
	}
	if (family_history=="S"){
		family_history <- "Sporadic"
	}
	text(0,height,paste("Phenotype: ",phenotype,sep=""),col=nicegrey,adj=0,cex=0.45)
	text(30,height,paste("Sex: ",sex,sep=""),col=nicegrey,adj=0,cex=0.45)
	text(60,height,paste("Family History: ",family_history,sep=""),col=nicegrey,adj=0,cex=0.45)
}

# Calculate median CI 
median_ci=function(data,range){
	tempdf<-data.frame(data=data,group=1)
	tempdf_median<- groupwiseMedian(data ~ group,data=tempdf,conf = 0.95,R= 500,percentile = TRUE,bca = FALSE,digits = 30)
	if(range=="lower"){
		return(tempdf_median$Percentile.lower)
	} 
	else if (range == "upper"){
		return(tempdf_median$Percentile.upper)
	}
}

# Create Plots 
age_plot_cumulative_publication = function (all_ages_x,gene_ages_x,variant_ages_x,variant,gene,phenotype,sex,family_history,p_gene,plot_letter){
	### select medians ###
	### calculate kruskal wallis p value ###
	if (length(variant_ages_x >0) & length(all_ages_x >0)){
		kruskal_p <- kruskal_test(variant_ages_x,all_ages_x)
	}
	if (length(variant_ages_x)>0){
		variant_ecdf<-ecdf(variant_ages_x)
		variant_ages_y<- variant_ecdf(variant_ages_x)
	}
	if (length(gene_ages_x)>0){
		gene_ecdf<-ecdf(gene_ages_x)
		gene_ages_y<- gene_ecdf(gene_ages_x)
	}
	if (length(all_ages_x)>0){
		dist_funct<-ecdf(all_ages_x)
	}
	all_ages_x <- sort(all_ages_x)
	all_ages_y <- dist_funct(all_ages_x)
	smooth_all_ages<-smooth.spline(all_ages_x, all_ages_y, spar=0.35)
  	#plot_margins_variable
  	plot(
  		smooth_all_ages,     
  		lwd=1,
  		col=darkgreen,
  		type="l",
  		xaxt="n",
  		yaxt="n",
  		las=1,
    		xlab="",
		ylab="",
		main=variant,
		cex.lab=0.77,
		cex.main=0.66,
		col.main=nicegrey, 
		col.lab=nicegrey,
		xlim=c(0,100),
		ylim=c(0,1.3)
  	)
  	title(ylab="Cumulative Proportion", xlab="Age",line=1.5, cex.lab=0.77)
  	#calculate and plot median 95% CI using bootstraps
  	# Plot median and CI for all carriers except P and LP of interest
  	if(length(sort(unique(all_ages_x)))>1){
  		median_all_ages <- median(all_ages_x)
  		all_ages_lower<- median_ci(all_ages_x,"lower")
  		all_ages_upper<- median_ci(all_ages_x,"upper")
  		lines(x=c(median_all_ages,median_all_ages),y=c(0,1),lty=2,lwd=1,col=darkgreen)
  		rect(all_ages_lower, 0, all_ages_upper, 1,border=NA,col=darkgreen_rgb_fade)
  	}
  	# Plot median and CI for VUS in gene of interest 
  	if(length(sort(unique(gene_ages_x)))>1){
  		this_median     <- median(gene_ages_x)
  		lines(x=c(this_median,this_median),y=c(0,1),col=nicedarkblue,lty=2,lwd=1)
  		ages_lower<- median_ci(gene_ages_x,"lower")
  		ages_upper<- median_ci(gene_ages_x,"upper")
  		rect(ages_lower, 0, ages_upper, 1,border=NA,col=nicedarkblue_rgb_fade)
  	}
  	# Plot median and CI for P and LP variants of interest 
  	if(length(sort(unique(variant_ages_x)))>1){
  		this_median     <- median(variant_ages_x)
  		lines(x=c(this_median,this_median),y=c(0,1),col=nicered,lty=2,lwd=1)
  		ages_lower<- median_ci(variant_ages_x,"lower")
  		ages_upper<- median_ci(variant_ages_x,"upper")
  		rect(ages_lower, 0, ages_upper, 1,border=NA,col=nicered_rgb_fade)
  	}
  	# Replot all ages line over CIs
  	lines(smooth_all_ages,lwd=2,col=darkgreen)
  	# Plot all age line in legend 
  	lines(x=c(-1,0.5),y=c(1.3*0.85,1.3*0.85),lwd=2,col=darkgreen)
  	# Plot all age text in legend 
  	text(5,(1.3*0.85),"All Cases Distribution",col=nicegrey,adj=0,cex=0.45)
  	# Plot kruskal test 
  	if (length(variant_ages_x >0)){
  		points(46.5,(1.3*0.85),pch=21,col="white",bg=darkgreen,cex=0.75,lwd=0.5)
  		points(45,(1.3*0.85),pch=21,col="white",bg=nicered,cex=0.75,lwd=0.5)
  		text(50,(1.3*0.85),paste("Test For Difference: p=",kruskal_p,sep=""),col=nicegrey,adj=0,cex=0.45)
  	}
  	# plot gene and variant points if there are any
  	if (length(gene_ages_x) >0){
  		points(
  			x=gene_ages_x,
  			y=gene_ages_y,
  			cex=0.75,
  			pch=21,
  			col="white",
  			bg=nicedarkblue,
  			lwd=0.5
  		)
  		points(x=0,y=(1.3*0.9),cex=0.75,pch=21,col="white",bg=nicedarkblue,lwd=0.5)
  		text(5,(1.3*0.9),gene,col=nicegrey,adj=0,cex=0.45)
  	}
  	else{
  		text(5,(1.3*0.9),paste("No observed ages for",gene," excluding ",variant," carriers",sep=""),col=nicegrey,adj=0,cex=0.45)
  	}
  	if (length(variant_ages_x)>0){
  		points(
  			x=variant_ages_x,
  			y=variant_ages_y,
  			cex=0.75,
  			pch=21,
  			col="white",
  			bg=nicered,
  			lwd=0.5
  		)
  		points(x=0,y=(1.3*0.95),cex=0.75,pch=21,col="white",bg=nicered,lwd=0.5)
  		text(5,(1.3*0.95),variant,col=nicegrey,adj=0,cex=0.45)
  	}
  	else {
  		text(5,(1.3*0.95),paste("No observed ages for ",variant," carriers",sep=""),col=nicegrey,adj=0,cex=0.45)
  	}
  	selections_legend(phenotype,sex,family_history,1.3)
  	axis(side=1,lwd=1,cex.axis=0.6,col.axis=nicegrey,col=nicegrey,at = seq(0, 100, by = 20),mgp = c(3,0.1,0))
  	axis(side=2,lwd=1,las=1,cex.axis=0.6,col.axis=nicegrey,col=nicegrey,at = seq(0, 1, by = 0.2),mgp = c(3,0.8,0))
  	box(lwd=2,col=nicegrey)
}

kruskal_test=function(data1,data2){
	tempdf<-data.frame(data=data1,group=1)
	tempdf<-rbind(tempdf,data.frame(data=data2,group=2))
	kruskal<-kruskal.test(data ~ group, data = tempdf)$p.value
	if (kruskal >= 0.001){
		return(round(kruskal,4))
	} 
	else {
		return(format(kruskal,digits=3,scientific=T))
	}
}

plot_letter_label=function(x_min,x_max,y_min,y_max,letter){
	y_val <- y_min+1.15*(y_max-y_min)
	x_val <- x_min-0.25*(x_max-x_min)
	text(x_val,y_val,letter,col=nicegrey,adj=0,cex=1.2,xpd=NA,font=2)
}


# PDF = 210x297 - margins are 35,20,20,20 mm 
# PDF = 8.3x11.7 - margins are 1.38,0.79,0.79,0.79 inches

pdf(paste(output_directory,"/Chapter2_FUS_AOO_Reanalysis_Results.pdf"),width=(8.3-1.38-0.8),height=(11.7-0.8-0.8)/4)
#par(mfrow=c(1,2), mar = c(3,3,1.5,1.5),mgp=c(0.2,0.1,0),las=0)
par(mfrow=c(1,2),mar=c(3,3,1.5,1))
# FUS:c.1512_1513delAG(p.[G505fs])
variant_ages_x <- all_ages_df$all_carriers_aoo[all_ages_df$HGVS=="FUS:c.1512_1513delAG(p.[G505fs])" & !is.na(all_ages_df$HGVS)]
variant_ages_x <- c(variant_ages_x,37)
gene_ages_x <- all_ages_df$all_carriers_aoo[all_ages_df$gene=="FUS" & !is.na(all_ages_df$gene) & all_ages_df$HGVS != "FUS:c.1512_1513delAG(p.[G505fs])" & all_ages_df$all_carriers_primary_phenotype!="FTD" & !is.na(all_ages_df$all_carriers_primary_phenotype)]
all_ages_x <- all_ages_df$all_carriers_aoo[all_ages_df$HGVS != "FUS:c.1512_1513delAG(p.[G505fs])" & all_ages_df$all_carriers_primary_phenotype!="FTD" & !is.na(all_ages_df$all_carriers_primary_phenotype)]
age_plot_cumulative_publication(all_ages_x,gene_ages_x,variant_ages_x,"FUS:c.1512_1513delAG(p.[G505fs])","FUS","ALS","All","All","FUS","A")
plot_letter_label(0,100,0,1.3,"A")
# FUS:c.684_686dupCGG(p.[G229dup])
variant_ages_x <- all_ages_df$all_carriers_aoo[all_ages_df$HGVS=="FUS:c.684_686dupCGG(p.[G229dup])" & !is.na(all_ages_df$HGVS)]
variant_ages_x <- c(variant_ages_x,55,66,49)
gene_ages_x <- all_ages_df$all_carriers_aoo[all_ages_df$gene=="FUS" & !is.na(all_ages_df$gene) & all_ages_df$HGVS != "FUS:c.684_686dupCGG(p.[G229dup])" & all_ages_df$all_carriers_primary_phenotype!="FTD" & !is.na(all_ages_df$all_carriers_primary_phenotype)]
all_ages_x <- all_ages_df$all_carriers_aoo[all_ages_df$HGVS != "FUS:c.684_686dupCGG(p.[G229dup])" & all_ages_df$all_carriers_primary_phenotype!="FTD" & !is.na(all_ages_df$all_carriers_primary_phenotype)]
age_plot_cumulative_publication(all_ages_x,gene_ages_x,variant_ages_x,"FUS:c.684_686dupCGG(p.[G229dup])","FUS","ALS","All","All","FUS","A")
plot_letter_label(0,100,0,1.3,"B")
dev.off()




################################################################################
################################################################################

# Plot 2 - Proportion Explained - Proportions explained in Cuba, North Am. South Am.

################################################################################
################################################################################


#####
# Define Functions 
#####


population_people_dataset_function=function(lit_review_pen){
	# Function to pull out population information individuals from lit_review_pen and split this based on '|'
  	population_people_dataset=subset(lit_review_pen[!is.na(lit_review_pen$population_carriers_count) &  lit_review_pen$population_carriers_count>0,],select=c("gene","HGVS","population_carriers_primary_phenotype","population_carriers_continent","population_carriers_nationality","population_carriers_family_history","population_carriers_pmid")) 
  	population_people_dataset=separate_rows(population_people_dataset,c("population_carriers_primary_phenotype","population_carriers_continent","population_carriers_nationality","population_carriers_family_history","population_carriers_pmid"),sep="\\|")
  	population_people_dataset$population_carriers_family_history<-gsub("F","Familial",population_people_dataset$population_carriers_family_history)
  	population_people_dataset$population_carriers_family_history<-gsub("S","Sporadic",population_people_dataset$population_carriers_family_history)
  	return(population_people_dataset)
}

blank_plot=function(xlim,ylim){
  	# Function to create blank 
  	plot(xlim/2,ylim/2,
  		col="white",       
  		pch=21,
  		bg="white",      
  		type="l",
  		xaxt="n",
  		yaxt="n",
  		xlab="",
  		ylab="",
    		#main="",
    		cex.lab=1.4,
    		cex.main=1.4,
    		main="",
    		col.main=nicegrey,
    		col.lab=nicegrey,
    		xlim=c(0,xlim),
   		ylim=c(0,ylim),
    		bty='n'
	)
}

scale_bar_plot_linear=function(left,right){
  	# Plot side by side linear scale bars 
  	# Horizontal Lines
  	lines(x=c(left,right),y=c(10,10),col=nicegrey,lty=1,lwd=0.6)
  	# Vertical Lines 
  	for (number in seq(0,1,0.2)){
  		lines(x=c(left+number*(right-left),left+number*(right-left)),y=c(10,8),col=nicegrey,lty=1,lwd=0.6)
  	}
  	# text for scale bar 
  	text(c(left,left+0.2*(right-left),left+0.4*(right-left),left+0.6*(right-left),left+0.8*(right-left),right),5,c("0","20","40","60","80","100%"),col=nicegrey,adj=0.5,cex=0.8)
}

log_position=function(numberofbreaks,value){
  	# Function to convert a value to it's position on a log scale 
  	100*(log10(value)+2)/numberofbreaks
}

scale_bar_plot_log=function(left,right){
  	# Plot side by side log scale bars 
  	# Horizontal Lines
  	# Need to create axis break 
  	lines(x=c(left,left+4),y=c(10,10),col=nicegrey,lty=1,lwd=0.6)
  	lines(x=c(left+7,right),y=c(10,10),col=nicegrey,lty=1,lwd=0.6)
  	# Vertical Lines 
  	for (number in seq(0,1,0.25)){
  		lines(x=c(left+number*(right-left),left+number*(right-left)),y=c(10,8),col=nicegrey,lty=1,lwd=0.6)
  	}
  	for (number in c(seq(0.01,0.09,0.01),seq(0.1,0.9,0.1),seq(1,9,1),seq(10,90,10))){
  		lines(x=c(left+log_position(4,number),left+log_position(4,number)),y=c(10,8),col=nicegrey,lty=1,lwd=0.6)
  	}
  	# text for scale bar 
  	text(c(left,left+0.25*(right-left),left+0.5*(right-left),left+0.75*(right-left),right),5,c("0","0.1","1","10","100%") ,col=nicegrey,adj=0.5,cex=0.8)
}

#####
# Parse Data 
#####

population_people_dataset <<- data.frame(population_people_dataset_function(lit_review_pen))
# The current existing Cuba data in journALS is a subset of this data 
population_people_dataset <<- population_people_dataset[population_people_dataset$population_carriers_nationality!="Cuba",]
population_studies <<- population_studies[population_studies$Country != "Cuba",]
# Get gene for Cuban Data 
cuban_population_people_dataset$gene <- gsub(":.*","",cuban_population_people_dataset$HGVS)
# Merge Back in new Cuban Data 
population_studies <<- rbind(population_studies,cuban_population_studies)
cuban_population_people_dataset <- subset(cuban_population_people_dataset,select=colnames(population_people_dataset))
population_people_dataset <<- rbind(population_people_dataset,cuban_population_people_dataset)
population_people_dataset<-data.frame(population_people_dataset)

#####
# Create Starting Variables 
#####

# Set Starting Plotting Variables
gap_between   <- 5
end_pos       <- 100
top           <- 90
bottom        <- 85
als_ftd_gap   <- 45
left.1        <- 75
right.1       <- left.1 + 100
left.2        <- right.1 + als_ftd_gap 
right.2       <- left.2 + 100
box_dimension <- 4

# Define North and South American Countries 
NAM_countries <- unique(population_studies$Country[population_studies$Continent=="NAM"])[unique(population_studies$Country[population_studies$Continent=="NAM"]) !="Cuba"]
SAM_countries <- unique(population_studies$Country[population_studies$Continent=="SAM"])[unique(population_studies$Country[population_studies$Continent=="SAM"]) !="Cuba"]

#####
# Create Plot
#####

pdf(paste(output_directory,"/Chapter2_Proportion_of_Cases_Explained_Results.pdf"),height=11.7,width=(8.3),onefile=T)

for (filter in list("Cuba",NAM_countries,SAM_countries)){


  	####
  	# OVERVIEW PLOT 
  	###
  	par(mfrow=c(4,1), mai = c(0.1,0.1,0.1,0.1))

  	####
  	# Will Create 3 2x2 plots - Pathogenic, Pathogenic and Likely Pathogenic, Reported Variants
  	###
  	for (regions_pathogenicity_selection in c(1,3)){ # 1 is P variants, 2 is P or LP variants, 3 is all variants 
  		regions_pathogenicity_selection1=list_of_pathogenic_variants[[as.numeric(regions_pathogenicity_selection)]]
    		# Create Data Frame of Mean Proportions 
    		row_count<-length(unique(population_studies$Gene[population_studies$Gene %in% p_or_lp_genes]))
    		mean_gene_dataframe<-data.frame(ALS_Familial=rep(NA,row_count),FTD_Familial=rep(NA,row_count),ALS_Sporadic=rep(NA,row_count),FTD_Sporadic=rep(NA,row_count))
    		rownames(mean_gene_dataframe) <- unique(population_studies$Gene[population_studies$Gene %in% p_or_lp_genes])
    		mean_gene_dataframe[is.na(mean_gene_dataframe)]<-NA
    		for(gene in rownames(mean_gene_dataframe)){
    			for (regions_phenotype_selection in c("ALS","ALS")){
    				for (regions_history_selection in c("Familial","Sporadic")){
    					total=sum(population_studies$Count[!is.na(population_studies$Gene) & population_studies$Gene==gene & !is.na(population_studies$Phenotype) & population_studies$Phenotype==regions_phenotype_selection & !is.na(population_studies$History) & population_studies$History==regions_history_selection & (population_studies$Country %in% filter)])
    					cases=nrow(population_people_dataset[!is.na(population_people_dataset$gene) & population_people_dataset$gene==gene & !is.na(population_people_dataset$HGVS) & population_people_dataset$HGVS %in% regions_pathogenicity_selection1 & !is.na(population_people_dataset$population_carriers_primary_phenotype) & population_people_dataset$population_carriers_primary_phenotype ==regions_phenotype_selection & !is.na(population_people_dataset$population_carriers_family_history) & population_people_dataset$population_carriers_family_history==regions_history_selection & (population_people_dataset$population_carriers_nationality %in% filter),])
    					mean_gene_dataframe[gene,][[paste(regions_phenotype_selection,regions_history_selection,sep="_")]] <- 100*binom.confint(x=cases,n=total,method='wilson')$mean
    				}
    			}
    		}
   	for (regions_history_selection in c("Familial","Sporadic")){
       		# Create blank plot 
       		par(mai = c(0,0.1,0.2,0.1))
       		if (regions_history_selection=="Familial"){
       			blank_plot(375,125)}
       			if (regions_pathogenicity_selection==1){
       				if (filter == "Cuba"){
       					text((left.1+right.2)/2,120,"Cuba: Pathogenic Variants",col=nicegrey,adj=0.5,cex=1.1,xpd=NA)
       					text(left.1-(0.12*(right.2-left.1)),120,col=nicegrey,"A",adj=0,cex=1.2,xpd=NA,font=2)
       				}
       				if (filter == NAM_countries){
       					text((left.1+right.2)/2,120,"North America: Pathogenic Variants",col=nicegrey,adj=0.5,cex=1.1,xpd=NA)
       					text(left.1-(0.12*(right.2-left.1)),120,"B",col=nicegrey,adj=0,cex=1.2,xpd=NA,font=2)
       				}
       				if (filter == SAM_countries){
       					text((left.1+right.2)/2,120,"South America: Pathogenic Variants",col=nicegrey,adj=0.5,cex=1.1,xpd=NA)
       					text(left.1-(0.12*(right.2-left.1)),120,col=nicegrey,"C",adj=0,cex=1.2,xpd=NA,font=2)
       				}
       			}
       			if (regions_pathogenicity_selection==3){
       				if (filter == "Cuba"){
       					text((left.1+right.2)/2,120,"Cuba: Reported Variants in Genes with Pathogenic or Likely Pathogenic Variants",col=nicegrey,adj=0.5,cex=1.1,xpd=NA)
       				}
       				if (filter == NAM_countries){
       					text((left.1+right.2)/2,120,"North America: Reported Variants in Genes with Pathogenic or Likely Pathogenic Variants",col=nicegrey,adj=0.5,cex=1.1,xpd=NA)
       				}
       				if (filter == SAM_countries){
       					text((left.1+right.2)/2,120,"South America: Reported Variants in Genes with Pathogenic or Likely Pathogenic Variants",col=nicegrey,adj=0.5,cex=1.1,xpd=NA)
       				}
       			}
       			if (regions_history_selection=="Familial"){
       				scale_bar_plot_log(left.1,right.1)
       			}
       			else {
       				scale_bar_plot_log(left.2,right.2)
       			}
       			regions_phenotype_selection="ALS"
        		# Find order to plot genes y reordering dataframe 
        		mean_gene_dataframe<-mean_gene_dataframe[order(-mean_gene_dataframe[[paste(regions_phenotype_selection,regions_history_selection,sep="_")]],as.character(rownames(mean_gene_dataframe))),]
        		# height of bottom box / text 
        		height=13
        		# Need to cycle through each gene 
        		for (gene in rev(rownames(mean_gene_dataframe))){
        			if (regions_history_selection == "Familial"){
        				left_pos=left.1
        			} 
        			else {
        				left_pos=left.2
        			}
        			colour<-vangogh_palette[match(gene,p_or_lp_genes)]
        			total=sum(population_studies$Count[!is.na(population_studies$Gene) & population_studies$Gene==gene & !is.na(population_studies$Phenotype) & population_studies$Phenotype==regions_phenotype_selection & !is.na(population_studies$History) & population_studies$History==regions_history_selection & (population_studies$Country %in% filter)])
        			cases=nrow(population_people_dataset[!is.na(population_people_dataset$gene) & population_people_dataset$gene==gene & !is.na(population_people_dataset$HGVS) & population_people_dataset$HGVS %in% regions_pathogenicity_selection1 & !is.na(population_people_dataset$population_carriers_primary_phenotype) & population_people_dataset$population_carriers_primary_phenotype ==regions_phenotype_selection & !is.na(population_people_dataset$population_carriers_family_history) & population_people_dataset$population_carriers_family_history==regions_history_selection  & (population_people_dataset$population_carriers_nationality %in% filter),])
        			mean_position=100*binom.confint(x=cases,n=total,method='wilson')$mean
          			# if mean_position is zero set to 0.01 which will plot at zero 
          			if (mean_position == 0 & !is.na(mean_position)){
          				mean_position=0.01
          			}
          			lower_position=100*binom.confint(x=cases,n=total,method='wilson')$lower
          			# if lower_position is zero set to 0.01 which will plot at zero 
          			if ((lower_position <= 0 | log_position(4,lower_position) < 0) & !is.na(lower_position)) {
          				lower_position=0.01
          			}
          			# if upper_position is zero set to 0.01 which will plot at zero 
          			upper_position=100*binom.confint(x=cases,n=total,method='wilson')$upper
          			# if mean_position is zero set to 0.01 which will plot at zero 
          			if ((upper_position <= 0 | log_position(4,upper_position)<0) & !is.na(upper_position)) {
          				upper_position=0.01
          			}
          			# Plot vertical lines
          			for (number in c(0.1,1,10,100)){
          				lines(x=c(left_pos+log_position(4,number),left_pos+log_position(4,number)),y=c(height,height+1),col="grey70",lty=2,lwd=0.6)
          				lines(x=c(left.2+log_position(4,number),left.2+log_position(4,number)),y=c(height,height+1),col="grey70",lty=2,lwd=0.6)
          			}
          			# Plot confidence interval
          			lines(x=c(left_pos+log_position(4,lower_position),left_pos+log_position(4,upper_position)),y=c(height,height),col=nicegrey,lty=1,lwd=0.6)
          			# Plot rectangle 
          			roundedRect(left_pos+log_position(4,mean_position)-box_dimension/2,height-box_dimension/2,left_pos+log_position(4,mean_position)+box_dimension/2,height+box_dimension/2,border=NA,col=colour, bothsame=TRUE, aspcorrect=TRUE,corfactor=6)
          			# Plot gene name on first pass
          			text(left_pos-5,height,gene,col=nicegrey,adj=1,cex=0.8,font=3)
          			# Move up to next box position 
          			height=height+5
          		}
          		text((left.1+right.1)/2,110,"ALS",col=nicegrey,adj=0.5,cex=0.8)
          		text((left.2+right.2)/2,110,"ALS",col=nicegrey,adj=0.5,cex=0.8)
          		if (regions_history_selection=="Familial"){
          			text(left.1-30,70,regions_history_selection,col=nicegrey,adj=0.5,cex=0.8,srt=90)
          		}
          		else {
          			text(left.2-30,70,regions_history_selection,col=nicegrey,adj=0.5,cex=0.8,srt=90)
          		}
          	}
      	}
  }
  dev.off()

################################################################################
################################################################################

# Plot 3 : AAFs in cases and controls, EPACTS results, and testing for oligogenic results

################################################################################
################################################################################

################################################################################
# EPACTS - This is plotted below with AAF 
################################################################################

#####
# Parse Data 
#####

# Define genes in study to ignore alignments to overlapping etc. 
list_of_study_genes <-c("ALS2","ANG","ATXN2","C21orf2","CHCHD10","CHMP2B","DAO","DCTN1","ELP3","ERBB4","FIG4","FUS","GRN","HNRNPA1","LMNB1","MAPT","MATR3","NEFH","NEK1","OPTN","PFN1","PRPH","PSEN1","PSEN2","SARM1","SETX","SIGMAR1","SOD1","SPAST","SPG11","SQSTM1","TAF15","TARDBP","TBK1","UBQLN2","UNC13A","VAPB","VCP")

# Subset data 
epacts_missense<- subset(epacts_missense,select=c(MARKER_ID,PVALUE))
epacts_lof<- subset(epacts_lof,select=c(MARKER_ID,PVALUE))

# Trim To Necessary Columns 
epacts_missense$MARKER_ID <- gsub('^.*_',"",epacts_missense$MARKER_ID)
epacts_lof$MARKER_ID <- gsub('^.*_',"",epacts_lof$MARKER_ID)

# Remove off targets -missense
epacts_missense<-epacts_missense[epacts_missense$MARKER_ID %in% list_of_study_genes,]
epacts_missense<-epacts_missense[epacts_missense$MARKER_ID != "ATXN2",]

# Remove off targets - LOF
epacts_lof<-epacts_lof[epacts_lof$MARKER_ID %in% list_of_study_genes,]
epacts_lof<-epacts_lof[epacts_lof$MARKER_ID != "ATXN2",]

# While all genes are in missense not all are in lof data - need to correct this for plot 
all_genes <- unique(c(epacts_lof$MARKER_ID,epacts_missense$MARKER_ID))
tempdf<-data.frame(MARKER_ID=all_genes,PVALUE=rep(NA,length(all_genes)))
tempdf <- tempdf[!(tempdf$MARKER_ID %in% epacts_lof$MARKER_ID),]
epacts_lof <- merge(epacts_lof,tempdf,all.x=T,all.y=T)

#####
# Create EPACTS Plot Function 
#####

epacts_plot_function=function(dataset,test_type,plot_letter){
	#Add Colour 
	dataset$colour <- "grey"
	dataset$colour <- ifelse(dataset$MARKER_ID %in% p_or_lp_genes,vangogh_palette[match(dataset$MARKER_ID,p_or_lp_genes)],dataset$colour)
	dataset$MARKER_ID <- factor(dataset$MARKER_ID,levels=rev(dataset$MARKER_ID))
	plot(
		dataset$PVALUE,
		dataset$MARKER_ID,
		xaxt="n",
		yaxt="n",
		xlim=c(1,0.0001),
		log="x",
		bty="n",
		xlab="",
		ylab="",
		bg=dataset$colour,
		pch=21,
		col=dataset$colour,
		lwd=0.01,
		cex=1.2
		)
	title(main=paste("Burden Testing (",test_type,")",sep=""),line=0.3, cex.main=1,col.main=nicegrey)
	title(xlab=expression('-log'[10]*' (p-value)'),line=2.2, cex.lab=1,col.lab=nicegrey,cex.main=1,col.main=nicegrey)
	# Plot main x axis points and labels 
	axis(side=1,lwd=0.6,las=1,cex.axis=1,col.axis=nicegrey,col=nicegrey,tick=c(rev(seq(0.2,1,by=0.1)),rev(seq(0.02,0.1,by=0.01)),rev(seq(0.002,0.01,by=0.001)),rev(seq(0.0001,0.001,by=0.0001))),at=c(1,0.1,0.01,0.001,0.0001),labels=c("0","1","2","3","4"))
	# Plot intermediate x axis ticks 
	axis(side=1,lwd=0.6,las=1,cex.axis=1,col.axis=nicegrey,col=nicegrey,at=c(rev(seq(0.2,1,by=0.1)),rev(seq(0.02,0.1,by=0.01)),rev(seq(0.002,0.01,by=0.001)),rev(seq(0.0001,0.001,by=0.0001))),labels=c("",rep("",8),"",rep("",8),"",rep("",8),"",rep("",8),""))
	# Plot grey swimming lane rectangles
	rect(1,c(seq(0.5,32.5,by=2)),0.0001,c(seq(1.5,33.5,by=2)),col=peach_rgb_fade ,border=NA)
	# Write gene names
	text(1.1,c(seq(1,33,by=1)),rev(dataset$MARKER_ID),col=nicegrey,adj=0,cex=0.8,pos=2,xpd=NA,font=3)
	# Draw p-value threshold
	lines(x=c(0.05/nrow(dataset),0.05/nrow(dataset)),y=c(0.5,33.5),lty=1,lwd=1,col=nicered)
	text(0.05/(nrow(dataset)+23),30.5,"p-value threshold",col=nicered,adj=0,cex=1,pos=1,xpd=NA,srt=90)
	# Replot points over grey boxes
	points(
		dataset$PVALUE,
		dataset$MARKER_ID,
		bg=dataset$colour,
		pch=21,
		col=dataset$colour,
		lwd=0.01,
		cex=1.2
		)
	# Plot Letter 
	plot_letter_label(1,-30,0,30.5,plot_letter)
}

################################################################################
# AAF - Compare AAFs in cases and controls 
################################################################################

# Define Potentially Pathogenic Impacts
functional_impacts <-c("structural_interaction_variant","missense_variant","exon_loss_variant","disruptive_inframe_insertion","conservative_inframe_insertion","conservative_inframe_deletion","splice_acceptor_variant","stop_gained","frameshift_variant","initiator_codon_variant","splice_donor_variant","start_lost","stop_lost","disruptive_inframe_deletion")
# Transparency of points for plot
point_transparency=0.5

#####
# Parse Data 
#####

# Get in standard chr,pos,ref,alt format (gemini outputs in 0 base so converting to 1 base)
colnames(snps_aafs) <- c("chr","pos","ref","alt","aaf_gnomad_all","gene","impact")
snps_aafs$pos <- snps_aafs$pos+1
snps_aafs$identifier <- paste(gsub("chr","",snps_aafs$chr),snps_aafs$pos,snps_aafs$ref,snps_aafs$alt,sep=":")
snps_aafs <- subset(snps_aafs,select=c("aaf_gnomad_all","gene","identifier","impact"))

# Get in standard chr,pos,ref,alt format (gemini outputs in 0 base so converting to 1 base)
colnames(indels_aafs) <- c("chr","pos","ref","alt","aaf_gnomad_all","gene","impact")
indels_aafs$pos <- indels_aafs$pos+1
indels_aafs$identifier <- paste(gsub("chr","",indels_aafs$chr),indels_aafs$pos,indels_aafs$ref,indels_aafs$alt,sep=":")
indels_aafs <- subset(indels_aafs,select=c("aaf_gnomad_all","gene","identifier","impact"))
gnomAD_aafs <- unique(rbind(snps_aafs,indels_aafs))

# Get in standard chr,pos,ref,alt format (gemini outputs in 0 base so converting to 1 base)
variant_sample_count.indels.cases$start <- variant_sample_count.indels.cases$start +1 
variant_sample_count.indels.cases$identifier <- paste(gsub("chr","",variant_sample_count.indels.cases$chr),variant_sample_count.indels.cases$start,variant_sample_count.indels.cases$ref,variant_sample_count.indels.cases$alt,sep=":")
variant_sample_count.indels.cases$case_aaf <- (variant_sample_count.indels.cases$num_het+2*variant_sample_count.indels.cases$num_hom_alt) / (variant_sample_count.indels.cases$num_het+2*variant_sample_count.indels.cases$num_hom_alt+2*variant_sample_count.indels.cases$num_hom_ref)
variant_sample_count.indels.cases <- subset(variant_sample_count.indels.cases,select=c("identifier","case_aaf"))

# Get in standard chr,pos,ref,alt format (gemini outputs in 0 base so converting to 1 base)
variant_sample_count.snps.cases$start <- variant_sample_count.snps.cases$start +1 
variant_sample_count.snps.cases$identifier <- paste(gsub("chr","",variant_sample_count.snps.cases$chr),variant_sample_count.snps.cases$start,variant_sample_count.snps.cases$ref,variant_sample_count.snps.cases$alt,sep=":")
variant_sample_count.snps.cases$case_aaf <- (variant_sample_count.snps.cases$num_het+2*variant_sample_count.snps.cases$num_hom_alt) / (variant_sample_count.snps.cases$num_het+2*variant_sample_count.snps.cases$num_hom_alt+2*variant_sample_count.snps.cases$num_hom_ref)
variant_sample_count.snps.cases <- subset(variant_sample_count.snps.cases,select=c("identifier","case_aaf"))
variant_sample_count.cases <- rbind(variant_sample_count.snps.cases,variant_sample_count.indels.cases)

# Get in standard chr,pos,ref,alt format (gemini outputs in 0 base so converting to 1 base)
variant_sample_count.indels.controls$start <- variant_sample_count.indels.controls$start +1 
variant_sample_count.indels.controls$identifier <- paste(gsub("chr","",variant_sample_count.indels.controls$chr),variant_sample_count.indels.controls$start,variant_sample_count.indels.controls$ref,variant_sample_count.indels.controls$alt,sep=":")
variant_sample_count.indels.controls$control_aaf <- (variant_sample_count.indels.controls$num_het+2*variant_sample_count.indels.controls$num_hom_alt) / (variant_sample_count.indels.controls$num_het+2*variant_sample_count.indels.controls$num_hom_alt+2*variant_sample_count.indels.controls$num_hom_ref)
variant_sample_count.indels.controls <- subset(variant_sample_count.indels.controls,select=c("identifier","control_aaf"))

# Get in standard chr,pos,ref,alt format (gemini outputs in 0 base so converting to 1 base)
variant_sample_count.snps.controls$start <- variant_sample_count.snps.controls$start +1 
variant_sample_count.snps.controls$identifier <- paste(gsub("chr","",variant_sample_count.snps.controls$chr),variant_sample_count.snps.controls$start,variant_sample_count.snps.controls$ref,variant_sample_count.snps.controls$alt,sep=":")
variant_sample_count.snps.controls$control_aaf <- (variant_sample_count.snps.controls$num_het+2*variant_sample_count.snps.controls$num_hom_alt) / (variant_sample_count.snps.controls$num_het+2*variant_sample_count.snps.controls$num_hom_alt+2*variant_sample_count.snps.controls$num_hom_ref)
variant_sample_count.snps.controls <- subset(variant_sample_count.snps.controls,select=c("identifier","control_aaf"))
variant_sample_count.controls <- rbind(variant_sample_count.snps.controls,variant_sample_count.indels.controls)
variant_sample_count <- merge(variant_sample_count.controls,variant_sample_count.cases,by="identifier",all.x=T,all.y=T)
variant_sample_count <- merge(gnomAD_aafs,variant_sample_count,by="identifier",all.x=T,all.y=T)

# Want to use transparent colours rather than full
variant_sample_count$red <- col2rgb(vangogh_palette[match(variant_sample_count$gene,p_or_lp_genes)])["red",]/(col2rgb(vangogh_palette[match(variant_sample_count$gene,p_or_lp_genes)])["red",]+col2rgb(vangogh_palette[match(variant_sample_count$gene,p_or_lp_genes)])["blue",]+col2rgb(vangogh_palette[match(variant_sample_count$gene,p_or_lp_genes)])["green",])
variant_sample_count$blue <- col2rgb(vangogh_palette[match(variant_sample_count$gene,p_or_lp_genes)])["blue",]/(col2rgb(vangogh_palette[match(variant_sample_count$gene,p_or_lp_genes)])["red",]+col2rgb(vangogh_palette[match(variant_sample_count$gene,p_or_lp_genes)])["blue",]+col2rgb(vangogh_palette[match(variant_sample_count$gene,p_or_lp_genes)])["green",])
variant_sample_count$green <- col2rgb(vangogh_palette[match(variant_sample_count$gene,p_or_lp_genes)])["green",]/(col2rgb(vangogh_palette[match(variant_sample_count$gene,p_or_lp_genes)])["red",]+col2rgb(vangogh_palette[match(variant_sample_count$gene,p_or_lp_genes)])["blue",]+col2rgb(vangogh_palette[match(variant_sample_count$gene,p_or_lp_genes)])["green",])
variant_sample_count$colour <- rgb(variant_sample_count$red,variant_sample_count$blue,variant_sample_count$green,alpha=point_transparency)
grey_red <- col2rgb("grey")["red",]/(col2rgb("grey")["red",]+col2rgb("grey")["blue",]+col2rgb("grey")["green",])
grey_blue <- col2rgb("grey")["blue",]/(col2rgb("grey")["red",]+col2rgb("grey")["blue",]+col2rgb("grey")["green",])
grey_green <- col2rgb("grey")["green",]/(col2rgb("grey")["red",]+col2rgb("grey")["blue",]+col2rgb("grey")["green",])
variant_sample_count$colour <- ifelse(is.na(variant_sample_count$colour)==T,rgb(grey_red,grey_blue,grey_green,alpha=point_transparency),variant_sample_count$colour)

# Shuffle DF to avoid clumps 
set.seed(42)
rows <- sample(nrow(variant_sample_count))
variant_sample_count <- variant_sample_count[rows, ]



################################################################################
# OLIGOGENIC BINOMIAL TEST 
################################################################################

#####
# Define Functions to create oligogenic p-value DFs and create plots 
#####
Create_Binomial_DF=function(dataset,gnomAD_AF,Case_total,Control_total){
	# Filter to gnomAD AF 
	dataset <- dataset[variant_sample_count_with_cases_filter$aaf_gnomad_all <= gnomAD_AF,]
	# One line per sample 
	dataset <- data.frame(separate_rows(dataset,c("samples")))

	# Get count of each sample
	dataset <-data.frame(table(dataset$samples))
	colnames(dataset) <- c("SAMPLE_ID","Variant_Count")

	# Filter to Just Samples and merge with cc status 
	dataset <- merge(dataset,cc,all.x=T,all.y=F)

	# Get Counts 
	Control_above_1 <- nrow(dataset[dataset$cc_status=="Control" & dataset$Variant_Count>1,])
	Case_above_1 <- nrow(dataset[dataset$cc_status=="Case" & dataset$Variant_Count>1,])
	Control_above_or_equal_to_1 <- nrow(dataset[dataset$cc_status=="Control" & dataset$Variant_Count>=1,])
	Case_above_or_equal_to_1 <- nrow(dataset[dataset$cc_status=="Case" & dataset$Variant_Count>=1,])


	# Perform Binomial Test 
	binom_test <- binom.test(Case_above_1,Case_total,((Case_above_or_equal_to_1/Case_total)*(Control_above_or_equal_to_1/Control_total)),alternative="greater")
	return(binom_test$p.value)
}

Create_Binomial_Plot=function(dataset,title,plot_letter){
	dataset$gnomAD_AF <- factor(dataset$gnomAD_AF,levels=rev(dataset$gnomAD_AF))
	plot(
		NULL,
		xaxt="n",
		yaxt="n",
		ylim=c(1,0.001),
		log="y",
		bty="n",
		xlab="",
		ylab="",
		bg=reallydarkblue,
		pch=21,
		col=reallydarkblue,
		type="p",
		lwd=0.1,
		cex=1.2,
		xlim=c(1,length(levels(dataset$gnomAD_AF)))
	)
  	points(
  		dataset$gnomAD_AF,
  		dataset$binoms,
  		cex=1.6,
  		pch=21,
  		bg=nicedarkblue,
  		col=nicedarkblue
  		)
  	# X - axis 
  	axis(side=1,lwd=1,cex.axis=0.8,col.axis=nicegrey,col=nicegrey,at=c(1,2,3,4,5),labels=c("1e-2","1e03","1e-4","1e-5","0"))#,at = seq(0, 100, by = 20),mgp = c(3,0.1,0))
  	# y axis 
  	axis(side=2,lwd=1,las=1,cex.axis=0.8,col.axis=nicegrey,col=nicegrey,at=c(1,0.1,0.01,0.001),labels=c("0","1","2","3"))#,at = seq(0, 1, by = 0.2),mgp = c(3,0.8,0))
  	box(lwd=2,col=nicegrey)
  	# Create Axis Titles (need some extra room with these)
  	title(ylab=expression('-log'[10]*' (p-value)'), xlab="gnomAD AF Cutoff",line=2.2, cex.lab=1,col.lab=nicegrey,cex.main=1,col.main=nicegrey)
  	# Create Plot Title (Needs less room)
  	title(main=title,line=1, cex.lab=1,col.lab=nicegrey,cex.main=1,col.main=nicegrey)
  	# Draw p-value threshold
	lines(x=c(10,0.000000000000001),y=c(0.05/5,0.05/5),lty=1,lwd=1,col=nicered)
	if(plot_letter=="E"){
		text(3.5,0.004,"p-value threshold",col=nicered,adj=0,cex=1,pos=1,xpd=NA)
  		plot_letter_label(1,7,1,0.1307,plot_letter)
	}
}


#####
# Parse SNPs Data 
#####

colnames(snps_aafs_with_cases) <- c("chr","pos","ref","alt","aaf_gnomad_all","gene","impact","samples","het_samples","hom_alt_samples")
snps_aafs_with_cases$pos <- snps_aafs_with_cases$pos+1
snps_aafs_with_cases$identifier <- paste(gsub("chr","",snps_aafs_with_cases$chr),snps_aafs_with_cases$pos,snps_aafs_with_cases$ref,snps_aafs_with_cases$alt,sep=":")
snps_aafs_with_cases <- subset(snps_aafs_with_cases,select=c("aaf_gnomad_all","gene","identifier","impact","samples","het_samples","hom_alt_samples"))

#####
# Parse INDELs Data 
#####

colnames(indels_aafs_with_cases) <- c("chr","pos","ref","alt","aaf_gnomad_all","gene","impact","samples","het_samples","hom_alt_samples")
indels_aafs_with_cases$pos <- indels_aafs_with_cases$pos+1
indels_aafs_with_cases$identifier <- paste(gsub("chr","",indels_aafs_with_cases$chr),indels_aafs_with_cases$pos,indels_aafs_with_cases$ref,indels_aafs_with_cases$alt,sep=":")
indels_aafs_with_cases <- subset(indels_aafs_with_cases,select=c("aaf_gnomad_all","gene","identifier","impact","samples","het_samples","hom_alt_samples"))

#####
# Merge SNPs and INDELs 
#####

gnomAD_aafs_with_cases <- unique(rbind(snps_aafs_with_cases,indels_aafs_with_cases))

#####
# Merge SNPs and INDELs with Variant Data Generated Above 
#####

variant_sample_count_with_cases <- merge(variant_sample_count.controls,variant_sample_count.cases,by="identifier",all.x=T,all.y=T)
variant_sample_count_with_cases <- merge(gnomAD_aafs_with_cases,variant_sample_count_with_cases,by="identifier",all.x=T,all.y=T)

#####
# PARSE COVERAGE DATA 
#####

# Subset and aggregate each plate - this combines all regions for each sample
coverage_plate_1_1 <- subset(coverage_plate_1_1,select=c("V1","V2"))
coverage_plate_1_1 <- aggregate(.~ V2,data=coverage_plate_1_1,mean)

coverage_plate_2_1 <- subset(coverage_plate_2_1,select=c("V1","V2"))
coverage_plate_2_1 <- aggregate(.~ V2,data=coverage_plate_2_1,mean)

coverage_plate_2_2 <- subset(coverage_plate_2_2,select=c("V1","V2"))
coverage_plate_2_2 <- aggregate(.~ V2,data=coverage_plate_2_2,mean)

coverage_plate_3_1 <- subset(coverage_plate_3_1,select=c("V1","V2"))
coverage_plate_3_1 <- aggregate(.~ V2,data=coverage_plate_3_1,mean)

coverage_plate_4_1 <- subset(coverage_plate_4_1,select=c("V1","V2"))
coverage_plate_4_1 <- aggregate(.~ V2,data=coverage_plate_4_1,mean)

# Merge all coverage plates 
coverage_df <- rbind(coverage_plate_1_1,coverage_plate_2_1)
coverage_df <- rbind(coverage_df,coverage_plate_2_2)
coverage_df <- rbind(coverage_df,coverage_plate_3_1)
coverage_df <- rbind(coverage_df,coverage_plate_4_1)

# Aggregate plates - this combines samples sequenced across different plates 
coverage_df <- aggregate(.~ V2,data=coverage_df,sum)

# Tidy names
coverage_df$V2 <- gsub(".bam","",coverage_df$V2)
colnames(coverage_df) <- c("SAMPLE_ID","Coverage")
coverage_df <- merge(coverage_df,cc,all.x=T,all.y=F)
samples_above_20X_coverage <- coverage_df$SAMPLE_ID[coverage_df$Coverage>20]


#####
# APPLY FILTERS TO DATA
#####

Case_total <- nrow(cc[cc$cc_status=="Case",])
Control_total <- nrow(cc[cc$cc_status=="Control",]) - 9 # 9 controls were excluded for coverage below 5X - this is used for the first two plots (E and F)

# Filter to functional variants 
variant_sample_count_with_cases_filter <- variant_sample_count_with_cases[variant_sample_count_with_cases$impact %in% functional_impacts,]

#####
# Create Three Datasets 
#####

##
# DF1 : this is filtered to functional variants across a range of gnomAD AFs 
##
binom_df1<-data.frame(gnomAD_AF=c(0,0.00001,0.0001,0.001,0.01),binoms=rep(NA,length(c(0,0.00001,0.0001,0.001,0.01))))
for (i in 1:nrow(binom_df1)){
	gnomAD_AF <- binom_df1$gnomAD_AF[i]
	binom_df1$binoms[i] <- Create_Binomial_DF(variant_sample_count_with_cases_filter,gnomAD_AF,Case_total,Control_total)
}

##
# DF2 : This is the same as DF1 but restricitng analysis to genes with P or LP variants 
##
variant_sample_count_with_cases_filter <- variant_sample_count_with_cases_filter[variant_sample_count_with_cases_filter$gene %in% p_or_lp_genes,]

binom_df2<-data.frame(gnomAD_AF=c(0,0.00001,0.0001,0.001,0.01),binoms=rep(NA,length(c(0,0.00001,0.0001,0.001,0.01))))
for (i in 1:nrow(binom_df2)){
	gnomAD_AF <- binom_df2$gnomAD_AF[i]
	binom_df2$binoms[i] <- Create_Binomial_DF(variant_sample_count_with_cases_filter,gnomAD_AF,Case_total,Control_total)
}

##
# DF3 : The same as DFs 1 and 20 but just testing samples with coverage above 20X (therefore total cases and controls is different) 
##
variant_sample_count_with_cases_filter <- variant_sample_count_with_cases_filter[variant_sample_count_with_cases_filter$samples %in% samples_above_20X_coverage,]
Case_total <- nrow(coverage_df[coverage_df$cc_status=="Case" & coverage_df$Coverage>20,])
Control_total <- nrow(coverage_df[coverage_df$cc_status=="Control" & coverage_df$Coverage>20,])
binom_df3<-data.frame(gnomAD_AF=c(0,0.00001,0.0001,0.001,0.01),binoms=rep(NA,length(c(0,0.00001,0.0001,0.001,0.01))))
for (i in 1:nrow(binom_df3)){
	gnomAD_AF <- binom_df3$gnomAD_AF[i]
	binom_df3$binoms[i] <- Create_Binomial_DF(variant_sample_count_with_cases_filter,gnomAD_AF,Case_total,Control_total)
}


################################################################################
# Create Plot 3 : AAFs in cases and controls, EPACTS results, and testing for oligogenic results
################################################################################


pdf("path/to/Chapter2_AAF_EPACTs_Oligogenic_Results.pdf",width=(8.3-1.38-0.8),height=(11.7-0.8-0.8))
layout(matrix(c(1,4,4,6,1,4,4,6,2,4,4,7,2,5,5,7,3,5,5,8,3,5,5,8), nrow=4,byrow = FALSE))
par(mar=c(3.2,3.2,2,1))

#######
# Plot 1
#######

blank_plot(1,1)
abline(0,1)
# Plot Ugly Grey Points
points(
	x=variant_sample_count$case_aaf[variant_sample_count$impact %in% functional_impacts & !(variant_sample_count$gene %in% p_or_lp_genes)],
	y=variant_sample_count$control_aaf[variant_sample_count$impact %in% functional_impacts & !(variant_sample_count$gene %in% p_or_lp_genes)],
	col="white",
    	pch=21, # good 
    	bg=variant_sample_count$colour,
    	cex=1.8,
    	lwd=0.1
    )
# Plot Pretty Colour Points 
points(
	x=variant_sample_count$case_aaf[variant_sample_count$impact %in% functional_impacts & variant_sample_count$gene %in% p_or_lp_genes],
	y=variant_sample_count$control_aaf[variant_sample_count$impact %in% functional_impacts & variant_sample_count$gene %in% p_or_lp_genes],
	col="white",
    	pch=21, 
    	bg=variant_sample_count$colour,
    	cex=1.8,
    	lwd=0.1
	)
axis(side=1,lwd=1,cex.axis=0.8,col.axis=nicegrey,col=nicegrey)#,at = seq(0, 100, by = 20),mgp = c(3,0.1,0))
axis(side=2,lwd=1,las=1,cex.axis=0.8,col.axis=nicegrey,col=nicegrey)#,at = seq(0, 1, by = 0.2),mgp = c(3,0.8,0))
box(lwd=2,col=nicegrey)
# Create Axis Titles (need some extra room with these)
title(ylab="Control AF", xlab="Case AF",line=2.2, cex.lab=1,col.lab=nicegrey,cex.main=1,col.main=nicegrey)
# Create Plot Title (Needs less room)
title(main="Observed AFs",line=1, cex.lab=1,col.lab=nicegrey,cex.main=1,col.main=nicegrey)
plot_letter_label(0,1,0,1,"A")

#######
# Plot 2
#######
blank_plot(0.05,0.05)
abline(0,1)
# Plot Ugly Grey Points
points(
	x=variant_sample_count$case_aaf[variant_sample_count$impact %in% functional_impacts & !(variant_sample_count$gene %in% p_or_lp_genes)],
	y=variant_sample_count$control_aaf[variant_sample_count$impact %in% functional_impacts & !(variant_sample_count$gene %in% p_or_lp_genes)],
	col="white",
    	pch=21, # good 
    	bg=variant_sample_count$colour,
    	cex=1.8,
    	lwd=0.1
	)
# Plot Pretty Colour Points
points(
	x=variant_sample_count$case_aaf[variant_sample_count$impact %in% functional_impacts & variant_sample_count$gene %in% p_or_lp_genes],
	y=variant_sample_count$control_aaf[variant_sample_count$impact %in% functional_impacts & variant_sample_count$gene %in% p_or_lp_genes],
	col="white",
    	pch=21, # good 
    	bg=variant_sample_count$colour,
    	cex=1.8,
    	lwd=0.1
)
axis(side=1,lwd=1,cex.axis=0.8,col.axis=nicegrey,col=nicegrey)#,at = seq(0, 100, by = 20),mgp = c(3,0.1,0))
axis(side=2,lwd=1,las=1,cex.axis=0.8,col.axis=nicegrey,col=nicegrey)#,at = seq(0, 1, by = 0.2),mgp = c(3,0.8,0))
box(lwd=2,col=nicegrey)
# Create Axis Titles (need some extra room with these)
title(ylab="Control AF", xlab="Case AF",line=2.2, cex.lab=1,col.lab=nicegrey,cex.main=1,col.main=nicegrey)
# Create Plot Title (Needs less room)
title(main="Observed AFs (<0.05)",line=1, cex.lab=1,col.lab=nicegrey,cex.main=1,col.main=nicegrey)
plot_letter_label(0,0.05,0,0.05,"B")

#########
# Plot Legend 
#########
blank_plot(100,100)
points(x=rep(0,11),y=seq(0,100,10),cex=1.8,col="white",bg=vangogh_palette[match(sort(p_or_lp_genes)[1:11],p_or_lp_genes)],pch=21)
text(x=rep(5,11),y=seq(0,100,10),sort(p_or_lp_genes)[1:11],col=nicegrey,adj=0,cex=0.8,font=3)
points(x=rep(40,10),y=seq(10,100,10),cex=1.8,col="white",bg=vangogh_palette[match(sort(p_or_lp_genes)[12:21],p_or_lp_genes)],pch=21)
text(x=rep(45,10),y=seq(10,100,10),sort(p_or_lp_genes)[12:21],col=nicegrey,adj=0,cex=0.8,font=3)
points(x=40,y=0,cex=1.8,col="white",bg="grey",pch=21)
text(x=45,y=0,"Other",col=nicegrey,adj=0,cex=0.8)

#########
# EPACTS PLOTS  
#########

par(mar=c(4,4,3,1))
#layout(matrix(c(1,2), nrow=1,byrow = TRUE))
epacts_plot_function(epacts_missense,"missense","C")
epacts_plot_function(epacts_lof,"LOF","D")

#########
# EPACTS PLOTS  
#########
Create_Binomial_Plot(binom_df1,"Oligogenic Testing","E")
Create_Binomial_Plot(binom_df2,"P or LP Genes","F")
Create_Binomial_Plot(binom_df3,"20X Coverage Filter","G")

dev.off()

################################################################################
################################################################################

# Plot 4 : ATXN2  - PLOTS AND ANALYSIS 

################################################################################
################################################################################

##########################################
# Neccessary Functions
##########################################
dataset_shuffle <- function(dataset){
	# Randomly shuffle dataframe to avoid confusing clumps of blue and red 
	set.seed(42)
	rows<-sample(nrow(dataset))
	dataset<-dataset[rows,]
	return(dataset)
}

##########################################
# Plotting Functions - For Later 
##########################################

rmsd_create_plot<-function(rmsd_df,letter){
	# A function to plot the RMSD cutoff 
	# create blank plot
	plot(
		x=1,  
		y=1,   
		col="white",
		bg="white",
		xaxt="n",
		yaxt="n",
		las=1,
		pch=21,
		cex=2,
		las=1,
		xlab="",
		ylab="",
		main="",
		cex.lab=0.8,
		cex.main=1.2,
		col.main=nicegrey,
		col.lab=nicegrey,
		xlim=c(0,15),
		ylim=c(0.6,1.6)
		)
	# axis titles 
	title(ylab="RMSD", xlab="Coverage",line=2.2, cex.lab=1,col.lab=nicegrey)
	title(main="ATXN2 Predictions",line=0.3, cex.main=1,col.main=nicegrey)
	 # Plot x axis 
	 axis(side=1,lwd=1,cex.axis=1,col.axis=nicegrey,col=nicegrey,at = seq(0,15,by=5))
	# Plot y axis
	axis(side=2,lwd=1,las=1,cex.axis=1,col.axis=nicegrey,col=nicegrey,at = seq(0.6,1.6,by=0.2))
	# Plot cutoff line 
	abline(v=2,col=niceorange,lwd=0.8,lty=2)
	# Plot rmsd line
	lines(rmsd_df$coverage_cutoff,rmsd_df$rmsd,col=reallydarkblue,lwd=2.2)
	# Plot outline 
	box(lwd=2,col=nicegrey)
	# Plot Letter 
	plot_letter_label(0,15,0.6,1.6,letter)
}

program_comparison_plot <- function(dataset_plotdf,x_axis_title,y_axis_title,plot_title,legend_print,letter){
	# A function to create plots B and C 
	##########################################
	# Format DF layout 
	##########################################
	colnames(dataset_plotdf)[2] <- "a1"
	colnames(dataset_plotdf)[3] <- "a2"
	##########################################
	# Create Variables and Columns for Plotting
	##########################################
	# create transparent fill colours 
	red_rgb_fade      <<- rgb(0.7529412, 0.1607843, 0.02352941,alpha=0.6)
	blue_rgb_fade <<- rgb(0.07843137,0.427451,0.5921569,alpha=0.6)
	dataset_plotdf$fill<-"white"
	dataset_plotdf$fill<-ifelse(dataset_plotdf$cc=="Control",blue_rgb_fade ,dataset_plotdf$fill)
	dataset_plotdf$fill<-ifelse(dataset_plotdf$cc=="Case",red_rgb_fade,dataset_plotdf$fill)
	# Create outline colours 
	dataset_plotdf$outline<-"white"
	dataset_plotdf$outline<-ifelse(dataset_plotdf$cc=="Control",nicedarkblue,dataset_plotdf$outline)
	dataset_plotdf$outline<-ifelse(dataset_plotdf$cc=="Case",nicered,dataset_plotdf$outline)
	##########################################
	# Shuffle 
	##########################################
	dataset_plotdf<-dataset_shuffle(dataset_plotdf)
	##########################################
	# Create Plot Area
	##########################################
	# create blank plot
	plot(
		x=1,  
		y=1,   
		col="white",
		bg="white",
		xaxt="n",
		yaxt="n",
		las=1,
		pch=21,
		cex=2,
		las=1,
		xlab="",
		ylab="",
		main="",
		cex.lab=0.8,
		cex.main=1.2,
		col.main=nicegrey,
		col.lab=nicegrey,
		xlim=c(9,35),
		ylim=c(9,41)
		)
	# axis titles 
	title(ylab=y_axis_title, xlab=x_axis_title,line=2.2, cex.lab=1,col.lab=nicegrey)
	title(main=plot_title,line=0.3, cex.main=1,col.main=nicegrey)
	# Plot ab line 
	abline(0,1,lwd=1.4,col="grey",lty=1)
	# Plot points 
	points(
		x=dataset_plotdf$a1,  
		y=dataset_plotdf$a2,   
		col=dataset_plotdf$outline,
		bg=dataset_plotdf$fill,
		xaxt="n",
		yaxt="n",
		las=1,
		pch=21,
		cex=1.4,
		las=1,
		lwd=0.6
		)
	# Plot x axis 
	axis(side=1,lwd=1,cex.axis=1,col.axis=nicegrey,col=nicegrey,at = seq(10,34,by=8))
	# Plot y axis
	axis(side=2,lwd=1,las=1,cex.axis=1,col.axis=nicegrey,col=nicegrey,at = seq(10,34,by=8))
	# Plot outline 
	box(lwd=2,col=nicegrey)
	if (legend_print == "legend_yes"){
		# Plot legend
		points(11,40,col=nicered,bg=red_rgb_fade,pch=21,cex=1.4,lwd=0.6)
		points(11,37,col=nicedarkblue,bg=blue_rgb_fade,pch=21,cex=1.4,lwd=0.6)
		text(13,40,"Case Allele",col=nicegrey,adj=0,cex=1)
		text(13,37,"Control Allele",col=nicegrey,adj=0,cex=1)
	}
	# Plot Letter 
	plot_letter_label(9,35,9,41,letter)
}


create_atxn2_length_plot=function(dataset,column_name,program,legend_print,letter){
	##########################################
	# This is a roundabout way of getting the data in the right format for barplot()
	##########################################

	# Create column counting longer alleles
	dataset$A1_Count <- str_count(dataset$ID,"A1")
	# Filter to longer alleles and coverage above 2
	dataset <- dataset[dataset$A1_Count==0 & dataset$ATXN2_Coverage>=2,]
	# Convert To Factors 
	dataset[[column_name]] <- as.factor(dataset[[column_name]])
	dataset$cc <- as.factor(dataset$cc)
	# Collapse Data Frame 
	dataset_table <- subset(dataset,select=c(column_name,"cc"))
	dataset_table <- table(dataset_table)
	dataset_table.1<- data.frame(allele=names(dataset_table[,"Case"]),Count=as.vector(dataset_table[,"Case"]),cc=rep("Case",length(as.vector(dataset_table[,"Case"]))))
	dataset_table.1$Count <- dataset_table.1$Count*100 / sum(dataset_table.1$Count)
	dataset_table.2<- data.frame(allele=names(dataset_table[,"Control"]),Count=as.vector(dataset_table[,"Control"]),cc=rep("Control",length(as.vector(dataset_table[,"Control"]))))
	dataset_table.2$Count <- dataset_table.2$Count*100 / sum(dataset_table.2$Count)
	dataset_table <- rbind(dataset_table.1,dataset_table.2)
	dataset_table <- reshape(dataset_table,                        
		idvar = "cc",
		timevar = "allele",
		direction = "wide"
		)
	# Rename rows and columns
	row.names(dataset_table) <- dataset_table$cc
	dataset_table <- dataset_table[ , 2:ncol(dataset_table)]
	colnames(dataset_table) <- levels(dataset[[column_name]])
	dataset_table <-as.matrix(dataset_table)
	##########################################
	# Plot 
	##########################################
	# Plot barplot
	bp <- barplot(dataset_table,
		beside=T,
		ylim=c(0,100),					
		col=c(nicered,nicedarkblue),
		lwd=2,
		yaxt="n",
		space=c(0,0.35),
		cex.names=1,
		col.axis=nicegrey,
		xlab="",
		cex.lab=1,
		col.lab=nicegrey,
		border=c(nicered,nicedarkblue)
		)
	# y axis 1 
	title(ylab=expression(paste(italic("ATXN2")," allele carrier frequency,",sep="")),cex.lab=1,col.lab=nicegrey,line=4.5)
	# y axis 2
	title(ylab="longer allele (%)",cex.lab=1,col.lab=nicegrey,line=3.5)
	# main title 
	title(main=paste(program," Allele Length",sep=""),cex.main=1,col.main=nicegrey,line=0.6)
	# grid lines 
	abline(h=c(0,25,50,75,100),lty=2,lwd=0.8,col="grey")
	# axis
	axis(2,at=c(0,25,50,75,100),lab=c(0,25,50,75,100),las=1,lwd=2,col=nicegrey,col.axis=nicegrey,cex.axis=1)
	axis(1,at=bp[1],lab="CAG repeats: ",las=1,col=nicegrey,col.axis=nicegrey,cex.axis=1,hadj=1,line=NA,tick=F)
	# replot bars over grid lines
	barplot(dataset_table,
		beside=T,
		ylim=c(0,100),					
		col=c(nicered,nicedarkblue),
		lwd=2,
		yaxt="n",
		space=c(0,0.35),
		cex.names=1,
		col.axis=nicegrey,
		xlab="",
		cex.lab=1,
		col.lab=nicegrey,
		border=c(nicered,nicedarkblue),
		add=T
		)
	# Plot Legend 
	if (legend_print == "legend_yes"){
		points(x=c(bp[14],bp[14]),y=c(95,85),col=c(nicered,nicedarkblue),bg=c(nicered,nicedarkblue),pch=21,cex=1.4,lwd=0.6)
		text(x=c(bp[15],bp[15]),y=c(95,85),c("Case","Control"),col=nicegrey,adj=0,cex=1)
	}
	# Plot Letter 
	plot_letter_label(min(bp),max(bp),0,100,letter)
}
##########################################
# Data Parsing etc.
##########################################

##########################################
# Convert case / control to long format  
##########################################
cc.1<-data.frame(ID=paste(cc$SAMPLE_ID,"A1",sep="."),cc=cc$cc_status)
cc.2<-data.frame(ID=paste(cc$SAMPLE_ID,"A2",sep="."),cc=cc$cc_status)
cc <- rbind(cc.1,cc.2)

##########################################
# Convert coverage file to long format  
##########################################

ATXN2_Coverage.1 <- data.frame(ID=paste(ATXN2_Coverage$SAMPLE_ID,"A1",sep="."),ATXN2_Coverage=ATXN2_Coverage$ATXN2_Coverage)
ATXN2_Coverage.2 <- data.frame(ID=paste(ATXN2_Coverage$SAMPLE_ID,"A2",sep="."),ATXN2_Coverage=ATXN2_Coverage$ATXN2_Coverage)
ATXN2_Coverage <- rbind(ATXN2_Coverage.1,ATXN2_Coverage.2)

##########################################
# Filter to Just Called Sites  
##########################################

HipSTR <- HipSTR[!is.na(HipSTR$ATXN2_a1),]
tredparse <- tredparse[tredparse$ATXN2_a1 != -1,]

##########################################
# Convert to a1 being the smaller allele 
##########################################
for (number in 1:nrow(tredparse)){
	tredparse$a1[number]<-min(tredparse$ATXN2_a1[number],tredparse$ATXN2_a2[number])
	tredparse$a2[number]<-max(tredparse$ATXN2_a1[number],tredparse$ATXN2_a2[number])
}
for (number in 1:nrow(HipSTR)){
	HipSTR$a1[number]<-min(HipSTR$ATXN2_a1[number],HipSTR$ATXN2_a2[number])
	HipSTR$a2[number]<-max(HipSTR$ATXN2_a1[number],HipSTR$ATXN2_a2[number])
}
##########################################
# Convert to long format DFs 
##########################################
# Tredparse
tredparse.parse.1<-data.frame(ID=paste(tredparse$SAMPLE_ID,"A1",sep="."),tred.allele=tredparse$a1)
tredparse.parse.2<-data.frame(ID=paste(tredparse$SAMPLE_ID,"A2",sep="."),tred.allele=tredparse$a2)
tredparse.merge<-rbind(tredparse.parse.1,tredparse.parse.2)
# HipSTR
hipstr.parse.1<-data.frame(ID=paste(HipSTR$SAMPLE_ID,"A1",sep="."),hip.allele=HipSTR$a1)
hipstr.parse.2<-data.frame(ID=paste(HipSTR$SAMPLE_ID,"A2",sep="."),hip.allele=HipSTR$a2)
hipstr.merge<-rbind(hipstr.parse.1,hipstr.parse.2)

##########################################
# Merge 
##########################################

tredparse.hipstr <- merge(tredparse.merge,hipstr.merge,by="ID",all.x=F,all.y=F)
tredparse.hipstr <- merge(tredparse.hipstr,cc,by="ID",all.x=F,all.y=F)
tredparse.hipstr <- merge(tredparse.hipstr,ATXN2_Coverage,by="ID",all.x=T,all.y=F)
print(paste(table(tredparse.hipstr$cc)['Case']/2," cases and ",table(tredparse.hipstr$cc)['Control']/2," controls had ATXN2 genotypes called by both HipSTR and TREDPARSE",sep=""))


##########################################
# Calculate RMSD
##########################################

# Get Difference 
tredparse.hipstr$diff <- tredparse.hipstr$tred.allele -tredparse.hipstr$hip.allele
rmsd<- (sum(tredparse.hipstr$diff^2)/nrow(tredparse.hipstr))^(1/2)
rmsd_df<-data.frame(coverage_cutoff=0:15)
rmsd_df$rmsd <- NA
for (i in 0:15){
	rmsd_df$rmsd[i+1] <- (sum(tredparse.hipstr[tredparse.hipstr$ATXN2_Coverage>=i,]$diff^2)/nrow(tredparse.hipstr[tredparse.hipstr$ATXN2_Coverage>=i,]))^(1/2)
}

print(paste(table(tredparse.hipstr$cc[tredparse.hipstr$ATXN2_Coverage>=2])['Case']/2," cases and ",table(tredparse.hipstr$cc[tredparse.hipstr$ATXN2_Coverage>=2])['Control']/2," controls had ATXN2 genotypes called by both HipSTR and TREDPARSE with coverage above 2X as per the chosen RMSD filter",sep=""))
print(paste("Applying a coverage filter of 2X yields a rmsd of ",rmsd_df$rmsd[rmsd_df$coverage_cutoff==2],sep=""))
##########################################
# Calculate odds ratio
##########################################
e<<-2.718281828459045


calculate_atxn2_or=function(dataset,column,cutoff){
	# Apply Coverage Filter 
	or_df <- dataset[dataset$ATXN2_Coverage>2,]
	# Only Retain Longest Allele 
	or_df$A1_count <- str_count(or_df$ID,"A1")
	or_df <- or_df[or_df$A1_count != 1,]
	# Filter to just column of interest 
	ort<-table(or_df[[column]],or_df$cc)
	# Count cases and controls falling above and below cutoff 
	cases_below <- sum(ort[rownames(ort)[as.numeric(rownames(ort))<cutoff],][,1])
	cases_above <- sum(ort[rownames(ort)[as.numeric(rownames(ort))>=cutoff],][,1])
	controls_below <- sum(ort[rownames(ort)[as.numeric(rownames(ort))<cutoff],][,2])
	controls_above <- sum(ort[rownames(ort)[as.numeric(rownames(ort))>=cutoff],][,2])
	# Create dataframe that we don't actually use 
	ort1 <- data.frame(cases=c(cases_above,cases_below),controls=c(controls_above,controls_below))
	rownames(ort1) <- c("27+","0-26")
	# Calculate odds ratio
	or<- (cases_above/controls_above)/(cases_below/controls_below)
	# Calculate upper CI 
	or_uci <- e^(log(or)+1.96*sqrt((1/cases_above)+(1/controls_above)+(1/cases_below)+(1/controls_below)))
    # Calculate lower CI 
    or_lci <- e^(log(or)-1.96*sqrt((1/cases_above)+(1/controls_above)+(1/cases_below)+(1/controls_below)))
    print(paste("For ",column," the observed OR is: ",round(or,2)," (95% CI: ",round(or_uci,2),"-",round(or_lci,2),")",sep=""))
}

calculate_atxn2_or(tredparse.hipstr,"tred.allele",27)
calculate_atxn2_or(tredparse.hipstr,"hip.allele",27)


# PDF = 210x297 - margins are 35,20,20,20 mm 
# PDF = 8.3x11.7 - margins are 1.38,0.79,0.79,0.79 inches

pdf(paste(output_directory,"/Chapter2_ATXN2_Results.pdf"),width=(8.3-1.38-0.8),height=(11.7-0.8-0.8)*2/5)
par(mar=c(4,4,3,1))
layout(matrix(c(1,1,2,2,3,3,4,4,4,5,5,5), nrow=2,byrow = TRUE))
rmsd_create_plot(rmsd_df,"A")
program_comparison_plot(tredparse.hipstr[tredparse.hipstr$ATXN2_Coverage>=0,],"TREDPARSE Prediction","HipSTR Prediction","No Coverage Filter","legend_no","B")
program_comparison_plot(tredparse.hipstr[tredparse.hipstr$ATXN2_Coverage>=2,],"TREDPARSE Prediction","HipSTR Prediction","2X Coverage Filter","legend_yes","C")
par(mar=c(4,6,3,1))
create_atxn2_length_plot(tredparse.hipstr,"tred.allele","TREDPARSE","legend_no","D")
create_atxn2_length_plot(tredparse.hipstr,"hip.allele","HipSTR","legend_yes","E")
dev.off()

################################################################################
################################################################################

# Plot 5 : Plot Pedigree and Labels 

################################################################################
################################################################################

# Have to attach the pedigree 
attach(ped)
# Convert it to usable pedigree format 
pedi <- pedigree(ID, Father, Mother, Sex.1.Male.2.Female., AffectionStatus.NA.Unknown.Uncertain.1.Unaffected.2.Affected.)#, startage=AgeOrAgeOfOnset, endage=Duration(months))
# Tell it how to print (want NAs to print blank)
strid=gsub("NA","",paste(AgeOrAgeOfOnset_.NA.Unknown.,Duration.months.,Genotype.1.0.heterozygous.1.1.homozygous.....ungenotyped.,Genotype2,Genotype3, sep="\n"))
# Call PDF 
pdf(paste(output_directory,"/Chapter2_2302_Pedigree.pdf"),width=(11.7-0.8-0.8),height=(8.3-1.38-0.8)*4/5)
par(mar=c(0,0,0,0))
plot(pedi, id=strid, cex=0.5, symbolsize=1, branch=0.8, density = c(-1, 35, 55, 25), packed=F)
dev.off()

################################################################################
################################################################################

# Plot 6 : VQSR  - PLOTS AND ANALYSIS 

################################################################################
################################################################################

call_VQSR_plot=function(dataset,x_min,x_max,y_min,y_max,y_variable,plot_letter){
	dataset <- matrix(dataset,ncol=5,byrow=T)
	dataset <- data.frame(x=dataset[,1],y=dataset[,2],retained=dataset[,3])
	dataset$colour <- reallydarkblue_rgb_fade
	dataset$colour[dataset$retained== -1 ] <- niceorange_rgb_fade
	blank_plot2(x_min,x_max,y_min,y_max)
	points(
		x=dataset$x,
		y=dataset$y,
		pch=21,
		cex=1,
		col=NA,
		bg=dataset$colour
	)
	axis(side=1,lwd=1,cex.axis=1,col.axis=nicegrey,col=nicegrey,at=c(40,60,80),labels=c(40,60,80))
	axis(side=2,lwd=1,las=1,cex.axis=1,col.axis=nicegrey,col=nicegrey,at=c(y_min,((y_max^2)^0.5-(y_min^2)^0.5)/2,y_max),labels=c(y_min,((y_max^2)^0.5-(y_min^2)^0.5)/2,y_max))
	box(lwd=2,col=nicegrey)
	title(main="VQSR",line=0.3, cex.main=1,col.main=nicegrey)
	title(xlab="MQ",ylab=y_variable,line=2.2, cex.lab=1,col.lab=nicegrey,cex.main=1,col.main=nicegrey)
	plot_letter_label(40,90,y_min,y_max,plot_letter)
	if(y_variable=="QD"){
		text(c(69,69),c(38,35),c("Filtered","Retained"),col=nicegrey,adj=0,cex=0.8,xpd=NA)
		points(c(67,67),c(38,35),pch=21,col=NA,bg=c(niceorange,reallydarkblue),cex=1)
	}
}

plot_letter_label=function(x_min,x_max,y_min,y_max,letter){
	y_val <- y_min+1.15*(y_max-y_min)
	x_val <- x_min-0.25*(x_max-x_min)
	text(x_val,y_val,letter,col=nicegrey,adj=0,cex=1.2,xpd=NA,font=2)
}

blank_plot2=function(x_min,x_max,y_min,y_max){
	# Function to create blank 
	plot((x_max-x_min)/2,(y_max-y_min)/2,
		col="white",       
		pch=21,
		bg="white",      
		type="l",
		xaxt="n",
		yaxt="n",
		xlab="",
		ylab="",
		#main="",
		cex.lab=1.4,
		cex.main=1.4,
		main="",
		col.main=nicegrey,
		col.lab=nicegrey,
		xlim=c(x_min,x_max),
		ylim=c(y_min,y_max),
		bty='n'
	)
}
pdf(paste(output_directory,"/Chapter2_VQSR_Filtering.pdf"),width=(8.3-1.38-0.8),height=(11.7-0.8-0.8)*2/5)
#par(mfrow=c(1,2), mar = c(3,3,1.5,1.5),mgp=c(0.2,0.1,0),las=0)
par(mfrow=c(2,3),mar=c(3.5,3.5,2,1))
call_VQSR_plot(MQ.QD,40,80,0,40,"QD","A")
call_VQSR_plot(MQ.MQRankSum,40,80,-7.5,7.5,"MQRankSum","B")
call_VQSR_plot(MQ.SOR,40,80,0,3,"SOR","C")
call_VQSR_plot(MQ.FS,40,80,0,30,"FS","D")
call_VQSR_plot(MQ.ReadPosRankSum,40,80,-7.5,7.5,"ReadPosRankSum","E")
dev.off()


################################################################################
################################################################################

# Exomes Duplication Rate, Exomes Coverage Rate & Filtering Numbers 

################################################################################
################################################################################

# No plots but useful statistics

mean_duplication_rate <- mean(duplication_rate$V2)
diff <- qnorm(.975)*(sd(duplication_rate$V2)/sqrt(nrow(duplication_rate)))
print(paste("In the exome study samples had a mean duplication rate of ",round(100*mean_duplication_rate,2),"% (95% CI: ",round(100*(mean_duplication_rate+diff),2),"-",round(100*(mean_duplication_rate-diff),2),"%)",sep=""))

mean_coverage_rate <- mean(coverage_rate$V3)
diff <- qnorm(.975)*(sd(coverage_rate$V3)/sqrt(nrow(coverage_rate)))
print(paste("In the exome study samples had a mean coverage of ",round(mean_coverage_rate,2),"X (95% CI: ",round(mean_coverage_rate-diff,2),"-",round(mean_coverage_rate+diff,2),"X)",sep=""))

# Get Filtering Numbers for Variants remaining in Cuban Pedigree 

print("Present in any family member: ")
system("gemini query -q \"SELECT * FROM variants\" --gt-filter \"(gt_types).(project='CUBA').(!=HOM_REF).(any)\" Cuban_PLS_Discordant_Exomes/Gemini_Database/Cub_PLS_Disc_Exomes.db | wc -l")
print("Present in any family member and passing sequencing filters: ")
system("gemini query -q \"SELECT * FROM variants where FILTER is NULL\" --gt-filter \"(gt_types).(project='CUBA').(!=HOM_REF).(any)\" Cuban_PLS_Discordant_Exomes/Gemini_Database/Cub_PLS_Disc_Exomes.db | wc -l")
print("Present in all family members and passing sequencing filters: ")
system("gemini query -q \"SELECT * FROM variants where FILTER is NULL\" --gt-filter \"(gt_types).(project='CUBA').(!=HOM_REF).(all)\" Cuban_PLS_Discordant_Exomes/Gemini_Database/Cub_PLS_Disc_Exomes.db | wc -l")
print("Present in all family members and passing sequencing filters and at a frequency below 0.1% AF in gnomAD: ")
system("gemini query -q \"SELECT * FROM variants where FILTER is NULL and aaf_gnomad_all < 0.001\" --gt-filter \"(gt_types).(project='CUBA').(!=HOM_REF).(all)\" Cuban_PLS_Discordant_Exomes/Gemini_Database/Cub_PLS_Disc_Exomes.db | wc -l")
print("Present in all family members and passing sequencing filters and at a frequency below 0.1% AF in gnomAD: ")
system("gemini query -q \"SELECT * FROM variants where FILTER is NULL and aaf_gnomad_all < 0.001 and (impact == 'conservative_inframe_deletion' or impact == 'conservative_inframe_insertion' or impact == 'disruptive_inframe_deletion disruptive_inframe_insertion' or impact == 'frameshift_variant' or impact == 'missense_variant' or impact == 'splice_donor_variant' or impact == 'splice_region_variant' or impact == 'stop_gained' or impact == 'stop_lost' or impact == 'structural_interaction_variant')\" --gt-filter \"(gt_types).(project='CUBA').(!=HOM_REF).(all)\" Cuban_PLS_Discordant_Exomes/Gemini_Database/Cub_PLS_Disc_Exomes.db | wc -l")
print("Present in all family members and passing sequencing filters and at a frequency below 0.1% AF in gnomAD and functional and not present in more than 10% of PLS samples: ")
system("gemini query -q \"SELECT * FROM variants where FILTER is NULL and aaf_gnomad_all < 0.001 and (impact == 'conservative_inframe_deletion' or impact == 'conservative_inframe_insertion' or impact == 'disruptive_inframe_deletion disruptive_inframe_insertion' or impact == 'frameshift_variant' or impact == 'missense_variant' or impact == 'splice_donor_variant' or impact == 'splice_region_variant' or impact == 'stop_gained' or impact == 'stop_lost' or impact == 'structural_interaction_variant')\" --gt-filter \"(gt_types).(project='CUBA').(!=HOM_REF).(all) and (gt_types).(project='PLS').(!=HOM_REF).(count<2)\" Cuban_PLS_Discordant_Exomes/Gemini_Database/Cub_PLS_Disc_Exomes.db | wc -l")
