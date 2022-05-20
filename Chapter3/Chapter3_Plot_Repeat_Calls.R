# Chapter 3 Create analysis plots and output useful statistics
# This a script to create all plots and analyses for chapter 3 of my thesis
# The script requires reading in parsed results files from in silico genotyping tools, see 'Read in Data'
################################################################################
################################################################################

rm(list=ls())

################################################################################
################################################################################

# Define Directories 

################################################################################
################################################################################

# Set output directory here 
output_directory=

################################################################################
################################################################################

# Load packages

################################################################################
################################################################################

library(stringr)
library(qqman) # for ehdenovo manhanhattan plot 
library(beeswarm)
library(plotrix) # For axis break

################################################################################
################################################################################

# Colour Variables

################################################################################
################################################################################


# Some nice colours 
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
# Colour fades
darkgreen_rgb_fade_light   <<- rgb(0.8117647,0.6666667,0.1568627,alpha=0.2)
nicegrey_rgb_fade_light     <<- rgb(0.3019608,0.3019608,0.3019608,alpha=0.2)
nicered_rgb_fade_light      <<- rgb(0.7529412, 0.1607843, 0.02352941,alpha=0.2) 
nicedarkblue_rgb_fade_light <<- rgb(0.07843137,0.427451,0.5921569,alpha=0.2) 
darkgreen_rgb_fade_heavy    <<- rgb(0.8117647,0.6666667,0.1568627,alpha=0.8)
nicegrey_rgb_fade_heavy     <<- rgb(0.3019608,0.3019608,0.3019608,alpha=0.8)
nicered_rgb_fade_heavy      <<- rgb(0.7529412, 0.1607843, 0.02352941,alpha=0.8)
nicedarkblue_rgb_fade_heavy <<- rgb(0.07843137,0.427451,0.5921569,alpha=0.8)
# Extra colours for use in regions analysis plot 
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
# This pallete is used for regions analysis plot 
vangogh_palette<<-  c(  nicered, hull,  peach,    sunset, niceorange,lightgreen,darkgreen,reallylightblue,nicelightblue,nicedarkblue,reallydarkblue,pink,green,lime,brown,bright_yellow,lavendar,algae,yellow1,purple,light_purple,pink_purple)

################################################################################
################################################################################

# Create Useful Variables

################################################################################
################################################################################

# These are the genes tested by ExpansionHunter3
gene_list_quotes<<-  c(	"AFF2",	"AR",	   "ATN1", "ATXN1",	"ATXN10",	"ATXN2",	"ATXN3",	"ATXN7",	"ATXN8OS",		"C9orf72",  "CACNA1A",  "CBL",	"CNBP",		"CSTB",			    "DIP2B",	"DMPK",	"FMR1",	"FXN",	"GIPC1",	"GLS",	"HTT",	"JPH3",	"NIPA1",	"NOP56",	"PABPN1",	"PHOX2B",	"PPP2R2B", "RFC1",	"TBP",	"TCF4", "HOXD13",        "ARX", "PAPBN1", "SOX3",         "ZIC2",        "HOXA13", "FOXL2",  "RUNX2")
repeat_motifs <<-    c(	"GCC",	"CAG",	 "CAG",	 "CAG",		"ATTCT",	"CAG",		"CAG",		"CAG",		"CTG/CTA",		"GGGGCC",	  "CAG",		  "CCG",	"CCTG",		"CCCCGCCCCGCG",	"CGG",		"CTG",	"CGG",	"GAA",	"GGC",		"GCA",	"CAG",	"CAG",	"GCG",		"GGCCTG",	"GCG",		"GCN",		"CAG",		 "AAGGG",	"CAG",	"CTG",  "GCN",   "CCG", "GCG",    "GCN",  "GCN",  "GCN",    "GCN",  "GCN")

cutoff_genes <-   c("AFF2", "AR", "ARX",  "ATN1", "ATXN1",  "ATXN10",   "ATXN2",  "ATXN3",  "ATXN7",  "ATXN8OS","C9orf72","CACNA1A",  "CBL",  "CNBP", "CSTB", "DIP2B",  "DMPK", "FMR1", "FOXL2","FXN",  "GIPC1","GLS",  "HOXA13", "HOXD13", "HTT",  "JPH3", "LRP12","MARCHF6","NIPA1","NOP56","NOTCH2NLA","NUTM2B", "PABPN1", "PHOX2B", "PPP2R2B","RAPGEF2","RFC1", "RUNX2","SAMD12", "SOX3","STARD7","TBP","TCF4", "TNRC6A","YEATS2","ZIC2")
cutoff_cutoffs <- c(200,    38,   23,     48,     39,       280,        32,       55,       36,       74,        30,      20,         100,    50,     30,     150,      50,     200,    19,     66,     97,     680,    24,       22,       35,     41,     90,     660,      NA,      650,    61,         40,       12,       24,       51,       22,       400,    27,     440,      26,    40,      47,   70,     NA,       800,     25)

################################################################################
################################################################################

# Useful Functions

################################################################################
################################################################################

# A function to place a letter in the upper left hand corner of a plot 
plot_letter_label=function(x_min,x_max,y_min,y_max,letter){
  y_val <- y_min+1.15*(y_max-y_min)
  x_val <- x_min-0.17*(x_max-x_min)
  text(x_val,y_val,letter,col=nicegrey,adj=0,cex=1.2,xpd=NA,font=2)
}

# A function to create colour fades 
see_through_colour=function(hexcode,fade_proportion){
  rgb(red=as.vector(col2rgb(hexcode))[1]/255, green=as.vector(col2rgb(hexcode))[2]/255, blue=as.vector(col2rgb(hexcode))[3]/255,alpha=fade_proportion)
}

# A function to create a blank plot. This is useful for have greater control of exact elements you want to plot.
blank_plot=function(xmin,xmax,ymin,ymax){
  # Function to create blank 
  plot((xmax-xmin)/2,(ymax-ymin)/2,
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
    xlim=c(xmin,xmax),
    ylim=c(ymin,ymax),
    bty='n')
}

plot_letter_label_b=function(x_min,x_max,y_min,y_max,letter){
  # This differs from original plot_letter_label function by giving a slighter wider x value 
  y_val <- y_min+1.15*(y_max-y_min)
  x_val <- x_min-0.26*(x_max-x_min)
  text(x_val,y_val,letter,col=nicegrey,adj=0,cex=1.2,xpd=NA,font=2)
}

################################################################################
################################################################################

# Read In Data 

################################################################################
################################################################################

##############################
# Read in case/control information
##############################

# Read in case / control information
pmid_vs_cc_status <- read.table("pmid_vs_cc_status",sep='\t',stringsAsFactors=F,header=T) 
# Rename columns
colnames(pmid_vs_cc_status) <- gsub("dna_bank","ID",colnames(pmid_vs_cc_status)) 
# Retain useful columns
ALS_CC_STATUS <- subset(pmid_vs_cc_status,select=c(Sample_well,cc,c9orf72_status_clean))
# Final column rename
colnames(ALS_CC_STATUS) <- c("Sample","CC","C9orf72_Status")

##############################
# Read in Results of in silico genotyping tools
##############################

tredparse <- read.table("tredparse_all_sample_results",sep='\t',stringsAsFactors=F,header=F)
hipstr <- read.table("HipSTR_all_sample_results",sep='\t',stringsAsFactors=F,header=F)
repeatseq <- read.table("RepeatSeq_all_sample_results",sep='\t',stringsAsFactors=F,header=F)
EHV3 <- read.table("EHv3_all_sample_results",sep='\t',stringsAsFactors=F,header=F)
EHV2 <- read.table("EHv2_all_sample_results",sep='\t',stringsAsFactors=F,header=F,na.strings=c(""," ","NA"))
STRetch <- read.table("STRetch_all_sample_results",sep='\t',stringsAsFactors=F,header=F)
GangSTR_NonTargeted <- read.table("GangSTR_NonTargeted_all_sample_results",sep='\t',stringsAsFactors=F,header=F)
GangSTR_Targeted <- read.table("GangSTR_Targeted_all_sample_results",sep='\t',stringsAsFactors=F,header=F)
EHdN_CC <- read.table("ExpansionHunterDeNovo.case_control.locus.tsv",sep='\t',stringsAsFactors=F,header=T)
EHdN_C9 <- read.table("ExpansionHunterDeNovo.C9orf72_test.locus.tsv",sep='\t',stringsAsFactors=F,header=T)
exSTRa <- read.table("exSTRa_all_sample_results",sep='\t',stringsAsFactors=F,header=F)


##############################
# Read in the Gold Standard PCR data 
##############################

PCR_Genotypes <- read.table("Gold_Standard_PCR_Genotypes.txt",sep='\t',stringsAsFactors=F,header=T)

##############################
# Read in Data For Plot Comparing Coverage Between Exome Samples  
##############################

# Repeat coverage of RCSI exome data 
RCSI_Repeat_Coverage <- read.csv("RCSI_Repeat_Coverage",sep='\t',header=F,stringsAsFactors=F)
# Read in ALS Coverage Stats
ALS_Repeat_Coverage <- read.csv("ALS_Repeat_Coverage",sep='\t',header=F,stringsAsFactors=F)
# Read in regions where bed files overlap repeat regions 
SeqCap_Regions <- read.csv("SeqCap_Exonic_Repeats_Intersection.bed",sep='\t',header=F,stringsAsFactors=F)  
SureSelect_Regions <- read.csv("SureSelect_Exonic_Repeats_Intersection.bed",sep='\t',header=F,stringsAsFactors=F) 
# Read in coordinates of exonic repeats 
Exonic_Repeats_Coordinates <- read.csv("Exonic_Repeats_Coordinates.bed",sep='\t',header=F,stringsAsFactors=F)
# File with structures of all genes
gene_coordinates      <<- read.csv("genomic_regions_for_coverage_genes",header=TRUE,sep='\t',na.strings=c("","NA"))
# Overall Exome Coverage Files 
ALS_DOC    <<- read.csv("ALS_DOC",header=F,sep='\t',na.strings=c("","NA"))
RCSI_DOC   <<- read.csv("RCSI_DOC",header=F,sep='\t',na.strings=c("","NA"))

##############################
# Data needed for follow up verification of putative positive expansions 
##############################

# TCF4 Coverage 
TCF4_coverage <- read.table("TCF4_WES_Coverage_Comparison",sep='\t',stringsAsFactors=F,header=F)
TCF4_exSTRa_WES <- read.table("RCSI_WES_exSTRa_Parsed",sep='\t',stringsAsFactors=F,header=T)
# NOTCH2 Coverage 
WES_NOTCH2_Coverage <- read.table("RCSI_WES_NOTCH2_Coverage",sep='\t',stringsAsFactors=F,header=F)

# This is used for C9orf72 plot otherwise we just use the p values 
ALS_WGS_PCR_FREE_exSTRa_scores_Targeted <- read.table("ALS_WGS_PCR_FREE_exSTRa_scores_Targeted",sep='\t',stringsAsFactors=F,header=T)

################################################################################
################################################################################

# Data Parsing 

################################################################################
################################################################################

##############################
# Functions for data parsing
##############################
filter_genes_failing_WES_WGS_comparison=function(primary_dataset,primary_dataset_doubles){
  # When comparing the results of WES and WGS sequencing some genes may be fine while others may give poor results due to the context around the gene in the WES data
  # This is a function to find those genes which give poor comparisons (rmsd>1) and remove these genes from the exome and PCR WGS data 
  # Create empty vector 
  WES_fail_genes <-c()
  # Cycle through each gene calculating the rmsd for that gene
  for (gene in sort(unique(primary_dataset_doubles$Gene[primary_dataset_doubles$Dataset=="RCSI_WES"]))){
    # Create working dataset 
    Working_Dataset=primary_dataset_doubles
    # Filter to gene 
    Working_Dataset <- Working_Dataset[Working_Dataset$Gene==gene,]
    # Find samples that appear twice (not all samples will have a call in both exome and genome for many reasons)
    retain_samples <- names(table(Working_Dataset$Sample)[table(Working_Dataset$Sample)>1])
    if (length(retain_samples)>0){
      # Filter to samples with both WES and WGS result 
      Working_Dataset <- Working_Dataset[Working_Dataset$Sample %in% retain_samples,]
      # Convert dataframe format 
      Working_Dataset_1 <- data.frame(Sample=c(paste(retain_samples,"_A1",sep=""),paste(retain_samples,"_A2",sep="")),WES=NA,WGS=NA)
      for (sample in retain_samples){
        Working_Dataset_1$WES[Working_Dataset_1$Sample==paste(sample,"_A1",sep="")] <- Working_Dataset$A1[Working_Dataset$Sample==sample & Working_Dataset$Dataset=="RCSI_WES"]
        Working_Dataset_1$WES[Working_Dataset_1$Sample==paste(sample,"_A2",sep="")] <- Working_Dataset$A2[Working_Dataset$Sample==sample & Working_Dataset$Dataset=="RCSI_WES"]
        Working_Dataset_1$WGS[Working_Dataset_1$Sample==paste(sample,"_A1",sep="")] <- Working_Dataset$A1[Working_Dataset$Sample==sample & Working_Dataset$Dataset=="RCSI_WGS_PCR_FREE"]
        Working_Dataset_1$WGS[Working_Dataset_1$Sample==paste(sample,"_A2",sep="")] <- Working_Dataset$A2[Working_Dataset$Sample==sample & Working_Dataset$Dataset=="RCSI_WGS_PCR_FREE"]
      }
      ###
      # Calculate rmsd
      ###
      Working_Dataset_1$WGS <- as.numeric(Working_Dataset_1$WGS)
      Working_Dataset_1$WES <- as.numeric(Working_Dataset_1$WES)
      Working_Dataset_1$diff <- ((Working_Dataset_1$WES-Working_Dataset_1$WGS)^2)^0.5
      rmsd <- (sum(Working_Dataset_1$diff^2)/nrow(Working_Dataset_1))^(1/2)
      if (rmsd >= 1) {
        WES_fail_genes <-c(WES_fail_genes,gene)
      }
    } else {
      WES_fail_genes <-c(WES_fail_genes,gene)
    }
  }
  # Exclude these genes in WES samples 
  primary_dataset <- primary_dataset[!((primary_dataset$Dataset=="RCSI_WES" | primary_dataset$Dataset=="ALS_WES") & primary_dataset$Gene %in% WES_fail_genes),]
  return(list(dataset=primary_dataset,fail_genes=WES_fail_genes))
}


##############################
# PARSE EACH DATASET 
##############################

##############################
# Data for depth of coverage (DOC) Plots
##############################

# Get max difference from 5' to 3' (use this later to ensure all plots are the same relative size )
max_offset=max(c(RCSI_Repeat_Coverage$V5,ALS_Repeat_Coverage$V5))
# Convert 'position' to correct base (easier when plotting gene coordinates etc)
RCSI_Repeat_Coverage$V5<-RCSI_Repeat_Coverage$V2+RCSI_Repeat_Coverage$V5-1
ALS_Repeat_Coverage$V5<-ALS_Repeat_Coverage$V2+ALS_Repeat_Coverage$V5-1
# Filter to columns of interest
RCSI_Repeat_Coverage_Trim <- subset(RCSI_Repeat_Coverage,select=c(V4,V5,V6))
ALS_Repeat_Coverage_Trim <- subset(ALS_Repeat_Coverage,select=c(V4,V5,V6))
# Calculate SD for coverage at each base
RCSI_Repeat_Coverage_SD <- aggregate(RCSI_Repeat_Coverage_Trim$V6,list(RCSI_Repeat_Coverage_Trim$V5,RCSI_Repeat_Coverage_Trim$V4),sd)
ALS_Repeat_Coverage_SD <- aggregate(ALS_Repeat_Coverage_Trim$V6,list(ALS_Repeat_Coverage_Trim$V5,ALS_Repeat_Coverage_Trim$V4),sd)
# Calculate Mean for coverage at each base 
RCSI_Repeat_Coverage_Mean <- aggregate(RCSI_Repeat_Coverage_Trim$V6,list(RCSI_Repeat_Coverage_Trim$V5,RCSI_Repeat_Coverage_Trim$V4),mean)
ALS_Repeat_Coverage_Mean <- aggregate(ALS_Repeat_Coverage_Trim$V6,list(ALS_Repeat_Coverage_Trim$V5,ALS_Repeat_Coverage_Trim$V4),mean)
# Combine mean and SD 
RCSI_Repeat_Coverage_Mean$SD <- RCSI_Repeat_Coverage_SD$x
ALS_Repeat_Coverage_Mean$SD <- ALS_Repeat_Coverage_SD$x

# Rename columns 
colnames(RCSI_Repeat_Coverage_Mean) <- c("Position","Gene","Mean_Coverage","Coverage_SD")
colnames(ALS_Repeat_Coverage_Mean) <- c("Position","Gene","Mean_Coverage","Coverage_SD")
# Calculate min and max values for x axis 
RCSI_Repeat_Coverage_Mean$Mean_plus_SD <- RCSI_Repeat_Coverage_Mean$Mean_Coverage + RCSI_Repeat_Coverage_Mean$Coverage_SD
RCSI_Repeat_Coverage_Mean$Mean_minus_SD <- pmax(0,(RCSI_Repeat_Coverage_Mean$Mean_Coverage - RCSI_Repeat_Coverage_Mean$Coverage_SD))
ALS_Repeat_Coverage_Mean$Mean_plus_SD <-  ALS_Repeat_Coverage_Mean$Mean_Coverage +  ALS_Repeat_Coverage_Mean$Coverage_SD
ALS_Repeat_Coverage_Mean$Mean_minus_SD <- pmax(0,( ALS_Repeat_Coverage_Mean$Mean_Coverage -  ALS_Repeat_Coverage_Mean$Coverage_SD))

##############################
# PCR Genotypes
##############################

# parse PCR Gene Names
PCR_Genotypes$Gene <- gsub(" ","",PCR_Genotypes$Gene)
PCR_Genotypes$Gene <- gsub("C9ORF72","C9orf72",PCR_Genotypes$Gene)
PCR_Genotypes$Gene <- gsub("-AS1","",PCR_Genotypes$Gene)
PCR_Genotypes$Gene <- gsub("NLC","",PCR_Genotypes$Gene)
# Retain Columns Of Interest
PCR_Genotypes <- unique(subset(PCR_Genotypes,select=c(ID,Gene,Allele1_repeat,Allele2_repeat)))
PCR_Genotypes$Sample_Gene <- paste(PCR_Genotypes$ID,PCR_Genotypes$Gene,sep="_")
# Remove Conflicting Calls 
PCR_Genotypes <- PCR_Genotypes[PCR_Genotypes$Sample_Gene %in% names(table(PCR_Genotypes$Sample_Gene)[table(PCR_Genotypes$Sample_Gene)==1]),]
# Merge With LP IDs
lp_subset <- subset(pmid_vs_cc_status,select=c(Sample_well,ID))
colnames(lp_subset) <- c("Sample","ID")
PCR_Genotypes <- merge(lp_subset,PCR_Genotypes,all.x=F,all.Y=F)
PCR_Genotypes <- subset(PCR_Genotypes,select=c(Sample,Gene,Allele1_repeat,Allele2_repeat))
colnames(PCR_Genotypes) <- c("Sample","Gene","PCR_A1","PCR_A2")
# Filter to Just Called Genotypes 
PCR_Genotypes <- PCR_Genotypes[PCR_Genotypes$PCR_A1>0 & PCR_Genotypes$PCR_A2>0 & !is.na(PCR_Genotypes$PCR_A1) & !is.na(PCR_Genotypes$PCR_A2),]

##############################
# GangSTR NonTargeted 
##############################


colnames(GangSTR_NonTargeted) <- c("Sample","Dataset","Gene","Alleles","Filter")
gangstr_loci <- c(  "chr2_191745599", "chr3_63898362",  "chr3_128891420", "chr4_3076604", "chr4_39350045",  "chr5_146258292", "chr6_16327867",  "chr6_170870996", "chr8_105601281", "chr9_27573527",  "chr9_71652203",  "chr10_81586140", "chr11_119077000",  "chr12_7045892",  "chr12_50898785", "chr12_112036755",  "chr13_70713516", "chr14_23790682", "chr14_92537355", "chr15_23086365", "chr16_87637894", "chr19_13318673", "chr19_14606854", "chr19_46273463", "chr20_2633382",  "chr21_45196324", "chr22_46191235", "chrX_25031777",  "chrX_66765160",  "chrX_146993569", "chrX_147582158")
gangstr_genes <- c( "GLS",            "ATXN7",          "CNBP",           "HTT",          "RFC1",           "PPP2R2B",        "ATXN1",          "TBP",            "LRP12",          "C9orf72",        "FXN",            "NUTM2B",         "CBL",              "ATN1",           "DIP2B",          "ATXN2",            "ATXN8OS",        "PABPN1",         "ATXN3",          "NIPA1",          "JPH3",           "CACNA1A",        "GIPC1",          "DM1",            "NOP56",          "CSTB",           "ATXN10",         "ARX",            "AR",             "FMR1",           "AFF2")
gangstr_fails <- c("GangSTRCallMinDepth","GangSTRCallMinDepth,GangSTRCallBadCI","GangSTRCallMinDepth,GangSTRCallSpanBoundOnly","GangSTRCallSpanBoundOnly","NOCALL")

# Step 1 :
# Filter
GangSTR_NonTargeted <- GangSTR_NonTargeted[!(GangSTR_NonTargeted$Filter %in% gangstr_fails) & !is.na(GangSTR_NonTargeted$Filter) & GangSTR_NonTargeted$Alleles !=".",]

# Convert Genes to standardised gene name 
GangSTR_NonTargeted$Gene <- gangstr_genes[match(GangSTR_NonTargeted$Gene,gangstr_loci)]

# Step 2 : 
# Split alleles and make A1 smallest A2 largest 
GangSTR_NonTargeted$Allele1 <- as.numeric(gsub(",.*","",GangSTR_NonTargeted$Alleles))
GangSTR_NonTargeted$Allele2 <- as.numeric(gsub("^.*,","",GangSTR_NonTargeted$Alleles))
GangSTR_NonTargeted$A1 <- ifelse(GangSTR_NonTargeted$Allele1 <= GangSTR_NonTargeted$Allele2, GangSTR_NonTargeted$Allele1,GangSTR_NonTargeted$Allele2)
GangSTR_NonTargeted$A2 <- ifelse(GangSTR_NonTargeted$Allele1 >= GangSTR_NonTargeted$Allele2, GangSTR_NonTargeted$Allele1,GangSTR_NonTargeted$Allele2)
GangSTR_NonTargeted <- unique(GangSTR_NonTargeted)
GangSTR_NonTargeted$A1 <- as.numeric(GangSTR_NonTargeted$A1)
GangSTR_NonTargeted$A2 <- as.numeric(GangSTR_NonTargeted$A2)

# Step 3 :
# The same individual may be included across multiple RCSI datasets - for RepeatSeq analysis want to retain in priority 1) WGS PCR FREE 2) WGS PCR 3) Exome 
# Find individuals in more than once 
data_subset <- unique(subset(GangSTR_NonTargeted,select=c(Sample,Dataset)))
# Find individuals present > once 
doubles<-names(table(data_subset$Sample)[table(data_subset$Sample)!=1])
# Create Sub Dataset Just of the double people 
GangSTR_NonTargeted_doubles <- GangSTR_NonTargeted[GangSTR_NonTargeted$Sample %in% doubles,]
# First remove those in the exome sequencing as this is lowest priority
GangSTR_NonTargeted<- GangSTR_NonTargeted[!(GangSTR_NonTargeted$Sample %in% doubles & GangSTR_NonTargeted$Dataset=="RCSI_WES"),]
# Find doubles still remaining 
# Find individuals in more than once 
data_subset <- unique(subset(GangSTR_NonTargeted,select=c(Sample,Dataset)))
# Find individuals present > once 
doubles<-names(table(data_subset$Sample)[table(data_subset$Sample)!=1])
# If there are still people to be removed
if (length(doubles)>0){
  # Then they should be removed from the lower priority WGS data 
  GangSTR_NonTargeted <- GangSTR_NonTargeted[!(GangSTR_NonTargeted$Sample %in% doubles & GangSTR_NonTargeted$Dataset=="RCSI_WGS_PCR"),]
}

# Step 4 :
# Add ALS cc status 
GangSTR_NonTargeted <- merge(GangSTR_NonTargeted,ALS_CC_STATUS,all.x=T,all.y=F)
# Need to count B and C to get parent samples 
GangSTR_NonTargeted$B_count <- str_count(GangSTR_NonTargeted$Sample,"B")
GangSTR_NonTargeted$C_count <- str_count(GangSTR_NonTargeted$Sample,"C")
GangSTR_NonTargeted$BC_count <- GangSTR_NonTargeted$B_count + GangSTR_NonTargeted$C_count
GangSTR_NonTargeted$Parent <- ifelse(GangSTR_NonTargeted$BC==1 & !is.na(GangSTR_NonTargeted$BC) & (GangSTR_NonTargeted$Dataset=="RCSI_WES" | GangSTR_NonTargeted$Dataset=="RCSI_WGS_PCR" | GangSTR_NonTargeted$Dataset=="RCSI_WGS_PCR_FREE"),"Yes","No")
GangSTR_NonTargeted <- subset(GangSTR_NonTargeted,select=c("Sample","Dataset","A1","A2","Gene","Parent","CC"))
# Get master Group
GangSTR_NonTargeted$Dataset_Master <- str_count(GangSTR_NonTargeted$Dataset,"RCSI")
GangSTR_NonTargeted$Dataset_Master <- ifelse(GangSTR_NonTargeted$Dataset_Master>=1,"RCSI","ALS")

# Step 5 : 
# Some genes may be acceptable to use the WES + WGS PCR data, others are not, taking a rmsd >= 1 as the cutoff and want to remove these 
GangSTR_NonTargeted_pre_exclusion <- GangSTR_NonTargeted
GangSTR_NonTargeted <- filter_genes_failing_WES_WGS_comparison(GangSTR_NonTargeted,GangSTR_NonTargeted_doubles)$dataset
GangSTR_NonTargeted_fail_genes <- filter_genes_failing_WES_WGS_comparison(GangSTR_NonTargeted,GangSTR_NonTargeted_doubles)$fail_genes

##############################
# GangSTR Targeted
##############################


colnames(GangSTR_Targeted) <- c("Sample","Dataset","Gene","Alleles","Filter")                                      

# Step 1 :
# Filter
GangSTR_Targeted <- GangSTR_Targeted[!(GangSTR_Targeted$Filter %in% gangstr_fails) & !is.na(GangSTR_Targeted$Filter) & GangSTR_Targeted$Alleles !=".",]
#GangSTR_Targeted <- GangSTR_Targeted[!(GangSTR_Targeted$Filter != "PASS"),]

# Convert Genes
GangSTR_Targeted$Gene <- gangstr_genes[match(GangSTR_Targeted$Gene,gangstr_loci)]

# Step 2 : 
# Split alleles and make A1 smallest A2 largest 
GangSTR_Targeted$Allele1 <- as.numeric(gsub(",.*","",GangSTR_Targeted$Alleles))
GangSTR_Targeted$Allele2 <- as.numeric(gsub("^.*,","",GangSTR_Targeted$Alleles))
GangSTR_Targeted$A1 <- ifelse(GangSTR_Targeted$Allele1 <= GangSTR_Targeted$Allele2, GangSTR_Targeted$Allele1,GangSTR_Targeted$Allele2)
GangSTR_Targeted$A2 <- ifelse(GangSTR_Targeted$Allele1 >= GangSTR_Targeted$Allele2, GangSTR_Targeted$Allele1,GangSTR_Targeted$Allele2)
GangSTR_Targeted <- unique(GangSTR_Targeted)
GangSTR_Targeted$A1 <- as.numeric(GangSTR_Targeted$A1)
GangSTR_Targeted$A2 <- as.numeric(GangSTR_Targeted$A2)

# Step 3 :
# The same individual may be included across multiple RCSI datasets - for RepeatSeq analysis want to retain in priority 1) WGS PCR FREE 2) WGS PCR 3) Exome 
# Find individuals in more than once 
data_subset <- unique(subset(GangSTR_Targeted,select=c(Sample,Dataset)))
# Find individuals present > once 
doubles<-names(table(data_subset$Sample)[table(data_subset$Sample)!=1])
# Create Sub Dataset Just of the double people 
GangSTR_Targeted_doubles <- GangSTR_Targeted[GangSTR_Targeted$Sample %in% doubles,]
# First remove those in the exome sequencing as this is lowest priority
GangSTR_Targeted<- GangSTR_Targeted[!(GangSTR_Targeted$Sample %in% doubles & GangSTR_Targeted$Dataset=="RCSI_WES"),]
# Find doubles still remaining 
# Find individuals in more than once 
data_subset <- unique(subset(GangSTR_Targeted,select=c(Sample,Dataset)))
# Find individuals present > once 
doubles<-names(table(data_subset$Sample)[table(data_subset$Sample)!=1])
# If there are still people to be removed
if (length(doubles)>0){
  # Then they should be removed from the lower priority WGS data 
  GangSTR_Targeted <- GangSTR_Targeted[!(GangSTR_Targeted$Sample %in% doubles & GangSTR_Targeted$Dataset=="RCSI_WGS_PCR"),]
}

# Step 4 :
# Add ALS cc status 
GangSTR_Targeted <- merge(GangSTR_Targeted,ALS_CC_STATUS,all.x=T,all.y=F)
# Need to count B and C to get parent samples 
GangSTR_Targeted$B_count <- str_count(GangSTR_Targeted$Sample,"B")
GangSTR_Targeted$C_count <- str_count(GangSTR_Targeted$Sample,"C")
GangSTR_Targeted$BC_count <- GangSTR_Targeted$B_count + GangSTR_Targeted$C_count
GangSTR_Targeted$Parent <- ifelse(GangSTR_Targeted$BC==1 & !is.na(GangSTR_Targeted$BC) & (GangSTR_Targeted$Dataset=="RCSI_WES" | GangSTR_Targeted$Dataset=="RCSI_WGS_PCR" | GangSTR_Targeted$Dataset=="RCSI_WGS_PCR_FREE"),"Yes","No")
GangSTR_Targeted <- subset(GangSTR_Targeted,select=c("Sample","Dataset","A1","A2","Gene","Parent","CC"))
# Get master Group
GangSTR_Targeted$Dataset_Master <- str_count(GangSTR_Targeted$Dataset,"RCSI")
GangSTR_Targeted$Dataset_Master <- ifelse(GangSTR_Targeted$Dataset_Master>=1,"RCSI","ALS")

# Step 5 : 
# Some genes may be acceptable to use the WES + WGS PCR data, others are not, taking a rmsd >= 1 as the cutoff and want to remove these 
GangSTR_Targeted_pre_exclusion <- GangSTR_Targeted
GangSTR_Targeted <- filter_genes_failing_WES_WGS_comparison(GangSTR_Targeted,GangSTR_Targeted_doubles)$dataset
GangSTR_Targeted_fail_genes <- filter_genes_failing_WES_WGS_comparison(GangSTR_Targeted,GangSTR_Targeted_doubles)$fail_genes


##############################
# STRetch
##############################

# Step 1:
# Convert to gene names 

STRetch_loci <- c("chr3_128891419_128891502", "chr4_3076665_3076695", "chr4_3076603_3076667", "chr4_39350044_39350103", "chr6_16327864_16327955", "chr9_27573482_27573544", "chr9_71652202_71652220", "chr2_191745598_191745646", "chr20_2633378_2633421",  "chr19_14606853_14606897",  "chr19_46273462_46273524",  "chr18_53253384_53253460",  "chr19_13318672_13318712",  "chr13_70713515_70713561",  "chr14_92537354_92537396",  "chr12_112036753_112036823",  "chr13_70713483_70713517",  "chrX_146993554_146993629", "chrX_147582124_147582273", "chrX_66765158_66765261", "chr22_46191234_46191304")
STRetch_genes <- c("CNBP",                     "HTT_exclude",          "HTT",                  "RFC1",                   "ATXN1",                  "C9orf72",                "FXN",                    "GLS",                      "NOP56",                  "GIPC1",                     "DM1",                      "TCF4",                     "CACNA1A",                  "ATXN8OS_exclude",                  "ATXN3",                    "ATXN2",              "ATXN8OS",                  "FMR1",                     "AFF2",                     "AR",                     "ATXN10")

colnames(STRetch) <- c("Dataset","Gene","Sample","Adjusted_PValue","A2")
# Convert to Gene Name 
STRetch$Gene <- STRetch_genes[match(STRetch$Gene,STRetch_loci)]
STRetch <- unique(STRetch)
STRetch$A2 <- round(STRetch$A2,0)

# Step 4 :
# Add ALS cc status 
STRetch <- merge(STRetch,ALS_CC_STATUS,all.x=T,all.y=F)
# Get master Group
STRetch$Dataset_Master <- str_count(STRetch$Dataset,"RCSI")
STRetch$Dataset_Master <- ifelse(STRetch$Dataset_Master>=1,"RCSI","ALS")

# The full list of genes searched by STRetch
STRetch_Full_Gene_List <- c("AFF2","AR","ARX","ATN1","ATXN1","ATXN10","ATXN2","ATXN3","ATXN7","ATXN8OS","C9orf72","CACNA1A","CNBP","DIP2B","DMPK","FMR1","FXN","GIPC1","GLS","HTT","JPH3","LRP12","NOP56","PAPBN1","PHOX2B","PPP2R2B","RFC1","TBP","TCF4")



##############################
# exSTRa
##############################
exSTRa_loci  <-  c("CANVAS","DBQD2",  "DM1",  "DM2",  "DRPLA",  "EPM1A","FAME1",  "FAME2",  "FAME3",    "FAME4",  "FAME6",  "FAME7",    "FECD3",  "FRAXA","FRAXE","FRDA", "FTDALS1",  "GDPAG","HD",   "HDL2", "NIID",   "OPDM1",  "OPML1",  "SBMA", "SCA1",   "SCA10",  "SCA12",   "SCA17", "SCA2",   "SCA3",   "SCA36",  "SCA6",   "SCA7", "SCA8")
exSTRa_genes <-  c("RFC1",  "XYLT1",  "DMPK", "CNBP", "ATN1",   "CSTB", "SAMD12", "STARD7", "MARCHF6",  "YEATS2", "TNRC6A", "RAPGEF2",  "TCF4",   "FMR1", "AFF2", "FXN",  "C9orf72",  "GLS",   "HTT",  "JPH3", "NOTCH2", "LRP12",  "NUTM2B", "AR",   "ATXN1",  "ATXN10", "PPP2R2B", "TBP",   "ATXN2",  "ATXN3",  "NOP56",  "CACNA1A","ATXN7","ATXN8OS")

# Rename Columns
colnames(exSTRa) <- c("Sample","Gene","P_Value","Signif","Dataset")
# Convert gene names
exSTRa$Gene <- exSTRa_genes[match(exSTRa$Gene,exSTRa_loci)]
# Filter missing calls 
exSTRa <- exSTRa[!is.na(exSTRa$P_Value),]
# Add CC Status 
exSTRa <- merge(exSTRa,ALS_CC_STATUS,all.x=T,all.y=F)

# Add ALS cc status 
exSTRa_ALS_WGS <- exSTRa[exSTRa$Dataset=="ALS_WGS_PCR_FREE",]
# Assume controls are C9 negative 
exSTRa_ALS_WGS$C9orf72_Status <- ifelse(is.na(exSTRa_ALS_WGS$C9orf72_Status)==T,0,exSTRa_ALS_WGS$C9orf72_Status)
# Filter C9 Data 
ALS_WGS_PCR_FREE_exSTRa_scores_Targeted$locus <- exSTRa_genes[match(ALS_WGS_PCR_FREE_exSTRa_scores_Targeted$locus,exSTRa_loci)]
ALS_WGS_PCR_FREE_exSTRa_scores_Targeted <- ALS_WGS_PCR_FREE_exSTRa_scores_Targeted[ALS_WGS_PCR_FREE_exSTRa_scores_Targeted$locus=="C9orf72",]
ALS_WGS_PCR_FREE_exSTRa_scores_Targeted <- subset(ALS_WGS_PCR_FREE_exSTRa_scores_Targeted,select=c(locus,rep,sample))
colnames(ALS_WGS_PCR_FREE_exSTRa_scores_Targeted) <- c("Gene","rep","Sample")
ALS_WGS_PCR_FREE_exSTRa_scores_Targeted <- merge(ALS_WGS_PCR_FREE_exSTRa_scores_Targeted,ALS_CC_STATUS,all.x=T,all.y=F)
ALS_WGS_PCR_FREE_exSTRa_scores_Targeted$C9orf72_Status <- ifelse(is.na(ALS_WGS_PCR_FREE_exSTRa_scores_Targeted$C9orf72_Status)==T,0,ALS_WGS_PCR_FREE_exSTRa_scores_Targeted$C9orf72_Status)

# This is TCF4 data because tredparse finds some significants that are to be explored 
# Convert locus to gene name 
TCF4_exSTRa_WES$locus <- exSTRa_genes[match(TCF4_exSTRa_WES$locus,exSTRa_loci)]
# Retain just columns of interest
TCF4_exSTRa_WES<- subset(TCF4_exSTRa_WES,select=c(locus,rep,sample))
colnames(TCF4_exSTRa_WES) <- c("Gene","rep","Sample")

##############################
# ExpansionHunter v2
##############################

colnames(EHV2) <- c("Sample","Dataset","Filter","Gene","Allele1","Allele2")
# Step 1 :
# Filter
EHV2 <- EHV2[EHV2$Filter=="PASS" & !is.na(EHV2$Filter),]
# Convert C9
EHV2$Gene <- gsub("C9ORF72","C9orf72",EHV2$Gene)

# Step 2 : 
# Split alleles and make A1 smallest A2 largest 

EHV2$A1 <- ifelse(EHV2$Allele1 <= EHV2$Allele2, EHV2$Allele1,EHV2$Allele2)
EHV2$A2 <- ifelse(EHV2$Allele1 >= EHV2$Allele2, EHV2$Allele1,EHV2$Allele2)
# Step 3 :
# The same individual may be included across multiple RCSI datasets - for RepeatSeq analysis want to retain in priority 1) WGS PCR FREE 2) WGS PCR 3) Exome 
# Find individuals in more than once 
data_subset <- unique(subset(EHV2,select=c(Sample,Dataset)))
# Find individuals present > once 
doubles<-names(table(data_subset$Sample)[table(data_subset$Sample)!=1])
# Create Sub Dataset Just of the double people 
EHV2_doubles <- EHV2[EHV2$Sample %in% doubles,]
# First remove those in the exome sequencing as this is lowest priority
EHV2<- EHV2[!(EHV2$Sample %in% doubles & EHV2$Dataset=="RCSI_WES"),]
# Find doubles still remaining 
# Find individuals in more than once 
data_subset <- unique(subset(EHV2,select=c(Sample,Dataset)))
# Find individuals present > once 
doubles<-names(table(data_subset$Sample)[table(data_subset$Sample)!=1])
# If there are still people to be removed
if (length(doubles)>0){
  # Then they should be removed from the lower priority WGS data 
  EHV2 <- EHV2[!(EHV2$Sample %in% doubles & EHV2$Dataset=="RCSI_WGS_PCR"),]
}

# Step 4 :
# Add ALS cc status 
EHV2 <- merge(EHV2,ALS_CC_STATUS,all.x=T,all.y=F)
# Need to count B and C to get parent samples 
EHV2$B_count <- str_count(EHV2$Sample,"B")
EHV2$C_count <- str_count(EHV2$Sample,"C")
EHV2$BC_count <- EHV2$B_count + EHV2$C_count
EHV2$Parent <- ifelse(EHV2$BC==1 & !is.na(EHV2$BC) & (EHV2$Dataset=="RCSI_WES" | EHV2$Dataset=="RCSI_WGS_PCR" | EHV2$Dataset=="RCSI_WGS_PCR_FREE"),"Yes","No")
EHV2 <- subset(EHV2,select=c("Sample","Dataset","A1","A2","Gene","Parent","CC"))
# Get master Group
EHV2$Dataset_Master <- str_count(EHV2$Dataset,"RCSI")
EHV2$Dataset_Master <- ifelse(EHV2$Dataset_Master>=1,"RCSI","ALS")

# Step 5 : 
# Some genes may be acceptable to use the WES + WGS PCR data, others are not, taking a rmsd >= 1 as the cutoff and want to remove these 
EHV2_pre_exclusion <- EHV2
EHV2 <- filter_genes_failing_WES_WGS_comparison(EHV2,EHV2_doubles)$dataset
EHV2_fail_genes <- filter_genes_failing_WES_WGS_comparison(EHV2,EHV2_doubles)$fail_genes

##############################
# ExpansionHunter v3 
##############################

colnames(EHV3) <- c("Sample","Dataset","Filter","Gene","Allele")

# Step 1 :
# Remove genotype calls we don't need (these are repeats near the pathogenic repeats but from the literature aren't involved in pathogenesis)  
EHV3_remove <-c("ATXN8OS","ATXN8OS_CTA","ATXN7_GCC","CNBP_CA","CNBP_CAGA","HTT_CCG","FXN_A","NOP56_CGCCTG")
EHV3 <- EHV3[!(EHV3$Gene %in% EHV3_remove),]
EHV3$Gene <- gsub("C9ORF72","C9orf72",EHV3$Gene)
# Convert ATXN8OS
EHV3$Gene <- gsub("_Combined","",EHV3$Gene)

# Step 2 :
# Filter failing genotypes 
EHV3 <- EHV3[EHV3$Filter=="PASS",]

# Step 3 
# Split alleles and make A1 smallest A2 largest 
EHV3$A1_temp <- as.numeric(gsub("/.*","",EHV3$Allele))
EHV3$A2_temp <- as.numeric(gsub("^.*/","",EHV3$Allele))
EHV3$A1 <- ifelse(EHV3$A1_temp <= EHV3$A2_temp, EHV3$A1_temp,EHV3$A2_temp)
EHV3$A2 <- ifelse(EHV3$A1_temp >= EHV3$A2_temp, EHV3$A1_temp,EHV3$A2_temp)

# Step 4 :
# The same individual may be included across multiple RCSI datasets - for RepeatSeq analysis want to retain in priority 1) WGS PCR FREE 2) WGS PCR 3) Exome 
# Find individuals in more than once 
data_subset <- unique(subset(EHV3,select=c(Sample,Dataset)))
# Find individuals present > once 
doubles<-names(table(data_subset$Sample)[table(data_subset$Sample)!=1])
# Create Sub Dataset Just of the double people 
EHV3_doubles <- EHV3[EHV3$Sample %in% doubles,]
# First remove those in the exome sequencing as this is lowest priority
EHV3<- EHV3[!(EHV3$Sample %in% doubles & EHV3$Dataset=="RCSI_WES"),]
# Find doubles still remaining 
# Find individuals in more than once 
data_subset <- unique(subset(EHV3,select=c(Sample,Dataset)))
# Find individuals present > once 
doubles<-names(table(data_subset$Sample)[table(data_subset$Sample)!=1])
# If there are still people to be removed
if (length(doubles)>0){
  # Then they should be removed from the lower priority WGS data 
  EHV3 <- EHV3[!(EHV3$Sample %in% doubles & EHV3$Dataset=="RCSI_WGS_PCR"),]
}


# Step 5 :
# Add ALS cc status 
EHV3 <- merge(EHV3,ALS_CC_STATUS,all.x=T,all.y=F)
# Need to count B and C to get parent samples 
EHV3$B_count <- str_count(EHV3$Sample,"B")
EHV3$C_count <- str_count(EHV3$Sample,"C")
EHV3$BC_count <- EHV3$B_count + EHV3$C_count
EHV3$Parent <- ifelse(EHV3$BC==1 & !is.na(EHV3$BC) & (EHV3$Dataset=="RCSI_WES" | EHV3$Dataset=="RCSI_WGS_PCR" | EHV3$Dataset=="RCSI_WGS_PCR_FREE"),"Yes","No")
EHV3 <- subset(EHV3,select=c("Sample","Dataset","A1","A2","Gene","Parent","CC"))
# Get master Group
EHV3$Dataset_Master <- str_count(EHV3$Dataset,"RCSI")
EHV3$Dataset_Master <- ifelse(EHV3$Dataset_Master>=1,"RCSI","ALS")

# Step 6 : 
# Some genes may be acceptable to use the WES + WGS PCR data, others are not, taking a rmsd >= 1 as the cutoff and want to remove these 
EHV3_pre_exclusion <- EHV3
EHV3 <- filter_genes_failing_WES_WGS_comparison(EHV3,EHV3_doubles)$dataset
EHV3_fail_genes <- filter_genes_failing_WES_WGS_comparison(EHV3,EHV3_doubles)$fail_genes



##############################
# RepeatSeq
##############################

colnames(repeatseq) <- c("Sample","Dataset","Gene","Allele")

# Step 1 :
# Convert position to gene name   
repeatseq_regions <-    c("chr3:63898361-63898392", "chr4:3076604-3076667", "chr5:146258291-146258322", "chr6:16327865-16327955", "chr6:170870995-170871105", "chr9:27573483-27573544", "chr12:112036754-112036823",  "chr13:70713516-70713561",  "chr14:92537355-92537396",  "chr19:13318673-13318712",  "chr19:46273463-46273524",  "chrX:146993555-146993629", "chr2:191745599-191745646", "chr3:128891420-128891499", "chr4:39350045-39350103", "chr4:41747993-41748039", "chr9:71652201-71652220", "chr11:119077000-119077033",  "chr12:7045880-7045938",  "chr12:50898785-50898807",  "chr14:23790681-23790701",  "chr15:23086365-23086396",  "chr16:87637889-87637935",  "chr18:53253385-53253460",  "chr19:14606854-14606897",  "chrX:25031771-25031814", "chrX:66765159-66765261", "chrX:147582125-147582273")
repeatseq_gene <-       c("ATXN7",                  "HTT",                  "PPP2R2B",                  "ATXN1",                  "TBP",                       "C9orf72",               "ATXN2",                      "ATXN8OS",                  "ATXN3",                    "CACNA1A",                  "DMPK",                     "FMR1",                     "GLS",                      "CNBP",                     "RFC1",                   "PHOX2B",                 "FXN",                    "CBL",                        "ATN1",                   "DIP2B",                    "PABPN1",                   "NIPA1",                    "JPH3",                     "TCF4",                     "GIPC1",                    "ARX",                    "AR",                     "AFF2")
 
# Convert from disease/gene to just gene 
repeatseq$Gene <- gsub("^.*/","",repeatseq$Gene)
# Convert position to gene 
repeatseq$Gene <- repeatseq_gene[match(repeatseq$Gene,repeatseq_regions)]


# Step 2 
# Remove samples and genes that weren't called 
repeatseq <- repeatseq[!is.na(repeatseq$Allele),]

# Step 3 
# Split alleles and calculate length 
repeatseq$Allele1 <- gsub("L.*","",repeatseq$Allele) # RepeatSeq output ist komisch, some of them have homozygous allele length followed by quality e.g. 56L:50 - in this case removing the L:50 which is the quality 
repeatseq$Allele1 <- gsub("h.*","",repeatseq$Allele1) # Heterozygotes are separated by a small h
repeatseq$Allele2 <- gsub("L.*","",repeatseq$Allele)
repeatseq$Allele2 <- gsub("^.*h","",repeatseq$Allele2)
# Get repeat units 
repeatseq$repeat_unit <- repeat_motifs[match(repeatseq$Gene,gene_list_quotes)]
repeat_length<-nchar(gsub("/.*","",repeatseq$repeat_unit))
repeatseq$Allele1 <- as.numeric(round(as.numeric(repeatseq$Allele1)/repeat_length,0))
repeatseq$Allele2 <- as.numeric(round(as.numeric(repeatseq$Allele2)/repeat_length,0))
repeatseq$A1 <- ifelse(repeatseq$Allele1 <= repeatseq$Allele2, repeatseq$Allele1,repeatseq$Allele2)
repeatseq$A2 <- ifelse(repeatseq$Allele1 >= repeatseq$Allele2, repeatseq$Allele1,repeatseq$Allele2)




# Step 4 :
# The same individual may be included across multiple RCSI datasets - for RepeatSeq analysis want to retain in priority 1) WGS PCR FREE 2) WGS PCR 3) Exome 
# Find individuals in more than once 
data_subset <- unique(subset(repeatseq,select=c(Sample,Dataset)))
# Find individuals present > once 
doubles<-names(table(data_subset$Sample)[table(data_subset$Sample)!=1])
# Create Sub Dataset Just of the double people 
repeatseq_doubles <- repeatseq[repeatseq$Sample %in% doubles,]
# First remove those in the exome sequencing as this is lowest priority
repeatseq<- repeatseq[!(repeatseq$Sample %in% doubles & repeatseq$Dataset=="RCSI_WES"),]
# Find doubles still remaining 
# Find individuals in more than once 
data_subset <- unique(subset(repeatseq,select=c(Sample,Dataset)))
# Find individuals present > once 
doubles<-names(table(data_subset$Sample)[table(data_subset$Sample)!=1])
# If there are still people to be removed
if (length(doubles)>0){
  # Then they should be removed from the lower priority WGS data 
  repeatseq <- repeatseq[!(repeatseq$Sample %in% doubles & repeatseq$Dataset=="RCSI_WGS_PCR"),]
}


# Step 5 :
# Add ALS cc status 
repeatseq <- merge(repeatseq,ALS_CC_STATUS,all.x=T,all.y=F)
# Need to count B and C to get parent samples 
repeatseq$B_count <- str_count(repeatseq$Sample,"B")
repeatseq$C_count <- str_count(repeatseq$Sample,"C")
repeatseq$BC_count <- repeatseq$B_count + repeatseq$C_count
repeatseq$Parent <- ifelse(repeatseq$BC==1 & !is.na(repeatseq$BC) & (repeatseq$Dataset=="RCSI_WES" | repeatseq$Dataset=="RCSI_WGS_PCR" | repeatseq$Dataset=="RCSI_WGS_PCR_FREE"),"Yes","No")
repeatseq <- subset(repeatseq,select=c("Sample","Dataset","A1","A2","Gene","Parent","CC"))
# Get master Group
repeatseq$Dataset_Master <- str_count(repeatseq$Dataset,"RCSI")
repeatseq$Dataset_Master <- ifelse(repeatseq$Dataset_Master>=1,"RCSI","ALS")

# Step 6 :

# Some genes may be acceptable to use the WES + WGS PCR data, others are not, taking a rmsd >= 1 as the cutoff and want to remove these 
repeatseq_pre_exclusion <- repeatseq
repeatseq <- filter_genes_failing_WES_WGS_comparison(repeatseq,repeatseq_doubles)$dataset
repeatseq_fail_genes <- filter_genes_failing_WES_WGS_comparison(repeatseq,repeatseq_doubles)$fail_genes

##############################
# HIPSTR
##############################
colnames(hipstr) <- c("Sample","Dataset","Gene","Ref_Allele","Filter","Alleles")
# Step 1 :
# HipSTR output doesn't give proper gene codes so need to convert these  
hipstr_disease <-         c("Human_STR_173941",  "Human_STR_266121", "Expansion_ATN1/DRPLA", "Expansion_FRA12A_MR/DIP2B",  "Expansion_SCA2/ATXN2", "Expansion_SCA8/ATXN8", "Expansion_OPMD/PAPBN1",  "Expansion_HDL2/JPH3","Expansion_SCA6/CACNA1A", "Human_STR_668860", "Expansion_DM1/DMPK", "Human_STR_886260", "Human_STR_803302", "Expansion_SCA7/ATXN7", "Expansion_HD/HTT", "Human_STR_1036600",  "Expansion_SCA12/PPP2R2B", "Expansion_SCA1/ATXN1","Expansion_FXS/FMR1",  "Human_STR_395660", "Human_STR_833719", "Human_STR_909208", "Human_STR_796394", "Human_STR_1475502",  "Human_STR_1542603",  "Human_STR_1597777")
hipstr_gene <-            c("NUTM2B",            "CBL",              "ATN1",                 "DIP2B",                      "ATXN2",                "ATXN8OS",              "PAPBN1",                 "JPH3",               "CACNA1A",                "GIPC1",            "DMPK",               "CSTB",             "GLS",              "ATXN7",                "HTT",              "RFC1",               "PPP2R2B",                 "ATXN1",               "FMR1",                "ZIC2",             "NOP56",            "ATXN10",           "HOXD13",           "C9orf72",            "ARX",                "SOX3")
hipstr_repeat_unit <-     c("CGG",               "CGG",              "CAG",                  "CGG",                        "CAG",                  "CTG/CTA",              "GCG",                    "CAG",                "CAG",                    "CCG",              "CTG",                "CGGGG",            "CAG",              "CAG",                  "CAG",              "AAAAG",              "CAG",                     "CAG",                 "CGG",                 "CAG/CTG/CGG",      "GGCCTG",           "ATTCT",            "CAG/CTG/CGG",      "GGGGCC",             "CCG",                "GCN") 


# Convert from disease/gene to just gene 
# Convert disease to gene 
hipstr$Gene <- hipstr_gene[match(hipstr$Gene,hipstr_disease)]

# Step 2 :
# Remove samples that fail filtering 
hipstr <- hipstr[hipstr$Filter == "PASS",]

# Step 3 :
# Calculate allele lengths 
# Calculate Reference Length
hipstr$ref_length <- nchar(hipstr$Ref_Allele)
hipstr$diff1 <- gsub("\\|.*","",hipstr$Alleles)
hipstr$diff2 <- gsub("^.*\\|","",hipstr$Alleles)
hipstr$A1_length <- hipstr$ref_length+as.numeric(hipstr$diff1)
hipstr$A2_length <- hipstr$ref_length+as.numeric(hipstr$diff2)
hipstr$A1 <- ifelse(hipstr$A1_length<=hipstr$A2_length,hipstr$A1_length,hipstr$A2_length)
hipstr$A2 <- ifelse(hipstr$A1_length>=hipstr$A2_length,hipstr$A1_length,hipstr$A2_length)


# Get Repeat Unit 
hipstr$repeat_unit <- hipstr$Gene
# Convert gene to repeat unit 
hipstr$repeat_unit <- hipstr_repeat_unit[match(hipstr$Gene,hipstr_gene)]

repeat_length<-nchar(gsub("/.*","",hipstr$repeat_unit))
hipstr$A1 <- round(hipstr$A1/repeat_length,0)
hipstr$A2 <- round(hipstr$A2/repeat_length,0)
hipstr <- subset(hipstr,select=c(-Alleles))
hipstr <- unique(hipstr)
# hipstr HTT genotype includes seconday repeat that isn't counted in other software 
hipstr$A2[hipstr$Gene=="HTT"] <- hipstr$A2[hipstr$Gene=="HTT"]-16
hipstr$A1[hipstr$Gene=="HTT"] <- hipstr$A1[hipstr$Gene=="HTT"]-16

# Step 4 :
# The same individual may be included across multiple RCSI datasets - for HipSTR analysis want to retain in priority 1) WGS PCR FREE 2) WGS PCR 3) Exome 
# Find individuals in more than once 
data_subset <- unique(subset(hipstr,select=c(Sample,Dataset)))
# Find individuals present > once 
doubles<-names(table(data_subset$Sample)[table(data_subset$Sample)!=1])
# Create Sub Dataset Just of the double people 
hipstr_doubles <- hipstr[hipstr$Sample %in% doubles,]
# First remove those in the exome sequencing as this is lowest priority
hipstr <- hipstr[!(hipstr$Sample %in% doubles & hipstr$Dataset=="RCSI_WES"),]
# Find doubles still remaining 
# Find individuals in more than once 
data_subset <- unique(subset(hipstr,select=c(Sample,Dataset)))
# Find individuals present > once 
doubles<-names(table(data_subset$Sample)[table(data_subset$Sample)!=1])
# If there are still people to be removed
if (length(doubles)>0){
  # Then they should be removed from the lower priority WGS data 
  hipstr <- hipstr[!(hipstr$Sample %in% doubles & hipstr$Dataset=="RCSI_WGS_PCR"),]
}

# Step 5 :
# Add ALS cc status 
hipstr <- merge(hipstr,ALS_CC_STATUS,all.x=T,all.y=F)
# Need to count B and C to get parent samples 
hipstr$B_count <- str_count(hipstr$Sample,"B")
hipstr$C_count <- str_count(hipstr$Sample,"C")
hipstr$BC_count <- hipstr$B_count + hipstr$C_count
hipstr$Parent <- ifelse(hipstr$BC==1 & !is.na(hipstr$BC) & (hipstr$Dataset=="RCSI_WES" | hipstr$Dataset=="RCSI_WGS_PCR" | hipstr$Dataset=="RCSI_WGS_PCR_FREE"),"Yes","No")
hipstr <- subset(hipstr,select=c("Sample","Dataset","A1","A2","Gene","Parent","CC"))
# Get master Group
hipstr$Dataset_Master <- str_count(hipstr$Dataset,"RCSI")
hipstr$Dataset_Master <- ifelse(hipstr$Dataset_Master>=1,"RCSI","ALS")

# Step 6 :
 
# Some genes may be acceptable to use the WES + WGS PCR data, others are not, taking a rmsd >= 1 as the cutoff and want to remove these 
hipstr_pre_exclusion <- hipstr
hipstr <- filter_genes_failing_WES_WGS_comparison(hipstr,hipstr_doubles)$dataset
hipstr_fail_genes <- filter_genes_failing_WES_WGS_comparison(hipstr,hipstr_doubles)$fail_genes

##############################
# TREDPARSE
##############################

# Step 1 :
# Tredpase output gives disease codes rather than gene names so lets convert these 
tredparse_disease <-  c("ALS",    "AR", "BPES", "CCD",  "CCHS",   "DM1",  "DM2",  "DRPLA","EIEE1","FRAXE",  "FRDA", "FXS",  "FXTAS",  "HD", "HDL",  "HFG",    "HPE5","OPMD",  "SBMA", "SCA1", "SCA10",  "SCA12",  "SCA17","SCA2", "SCA3", "SCA36","SCA6",   "SCA7", "SCA8",   "SD5",    "ULD",  "XLMR", "NOTCH2NLA", "NUTM2B", "CBL",  "DIP2B",  "PABPN1","NIPA1","TCF4","GIPC1","GLS","RFC1","LRP12","FXN")
tredparse_gene <-     c("C9orf72","AR", "FOXL2","RUNX2","PHOX2B", "DMPK", "CNBP", "ATN1", "ARX",  "AFF2",   "FXN",  "FMR1", "FMR1",   "HTT","JPH3", "HOXA13", "ZIC2","PABPN1","AR",   "ATXN1","ATXN10", "PPP2R2B","TBP",  "ATXN2","ATXN3","NOP56","CACNA1A","ATXN7","ATXN8OS","HOXD13", "CSTB", "SOX3", "NOTCH2",    "NUTM2B", "CBL",  "DIP2B",  "PABPN1","NIPA1","TCF4","GIPC1","GLS","RFC1","LRP12","FXN")

tredparse$V3 <- tredparse_gene[match(tredparse$V3,tredparse_disease)]

tredparse <- unique(tredparse)

# Step 2 :
# Make the longer allele A2 and the shorter allele A1
# Remove Sites that have not been called (e.g. non exonic repeats in exome samples)
tredparse <- tredparse[tredparse$V4 != "-1" & tredparse$V5 != "-1",]
tredparse$A1 <- ifelse(tredparse$V4 <= tredparse$V5, tredparse$V4,tredparse$V5) 
tredparse$A2 <- ifelse(tredparse$V4 >= tredparse$V5, tredparse$V4,tredparse$V5) 
# Rename Columns 
colnames(tredparse) <- c("Sample","Dataset","Gene","V4","V5","A1","A2")

# Step 3 :
# The same individual may be included across multiple RCSI datasets - for tredparse analysis want to retain in priority 1) WGS PCR FREE 2) WGS PCR 3) Exome 
# Find individuals in more than once 
data_subset <- unique(subset(tredparse,select=c(Sample,Dataset)))
# Find individuals present > once 
doubles<-names(table(data_subset$Sample)[table(data_subset$Sample)!=1])
# Create Sub Dataset Just of the double people 
tredparse_doubles <- tredparse[tredparse$Sample %in% doubles,]
# First remove those in the exome sequencing as this is lowest priority
tredparse<- tredparse[!(tredparse$Sample %in% doubles & tredparse$Sample=="RCSI_WES"),]
# Find doubles still remaining 
# Find individuals in more than once 
data_subset <- unique(subset(tredparse,select=c(Sample,Dataset)))
# Find individuals present > once 
doubles<-names(table(data_subset$Sample)[table(data_subset$Sample)!=1])
# If there are still people to be removed
if (length(doubles)>0){
  # Then they should be removed from the lower priority WGS data 
  tredparse<- tredparse[!(tredparse$Sample%in% doubles & tredparse$Dataset=="RCSI_WGS_PCR"),]
}

# Step 5 :
# Add ALS cc status 
tredparse <- merge(tredparse,ALS_CC_STATUS,all.x=T,all.y=F)
# Need to count B and C to get parent samples 
tredparse$B_count <- str_count(tredparse$Sample,"B")
tredparse$C_count <- str_count(tredparse$Sample,"C")
tredparse$BC_count <- tredparse$B_count + tredparse$C_count
tredparse$Parent <- ifelse(tredparse$BC==1 & !is.na(tredparse$BC) & (tredparse$Dataset=="RCSI_WES" | tredparse$Dataset=="RCSI_WGS_PCR" | tredparse$Dataset=="RCSI_WGS_PCR_FREE"),"Yes","No")
tredparse <- subset(tredparse,select=c("Sample","Dataset","A1","A2","Gene","Parent","CC"))

# Get master Group
tredparse$Dataset_Master <- str_count(tredparse$Dataset,"RCSI")
tredparse$Dataset_Master <- ifelse(tredparse$Dataset_Master>=1,"RCSI","ALS")

# Step 6 :
 
# Some genes may be acceptable to use the WES + WGS PCR data, others are not, taking a rmsd >= 1 as the cutoff and want to remove these 
tredparse_pre_exclusion <- tredparse
tredparse <- filter_genes_failing_WES_WGS_comparison(tredparse,tredparse_doubles)$dataset
tredparse_fail_genes <- filter_genes_failing_WES_WGS_comparison(tredparse,tredparse_doubles)$fail_genes


################################################################################
################################################################################

# Plots  

################################################################################
################################################################################

################################################################################
################################################################################

# Plotting Functions 

################################################################################
################################################################################

create_gene_table=function(temp_df,gene){
  # Create a table from a dataframe  
  temp_df <- temp_df[temp_df$Gene==gene,]
  # Retain neccessary columns
  temp_df <- subset(temp_df,select=c(Status,A2))
  temp_df$A2 <- as.numeric(temp_df$A2)
  # Convert to table 
   temp_df <- table(temp_df)
  # Convert to percentages
  temp_df["Patient",] <- temp_df["Patient",]*100/sum(temp_df["Patient",])
  temp_df["Control",] <- temp_df["Control",]*100/sum(temp_df["Control",])
  return(temp_df)
}

create_bp=function(gene){
# Function to create barplot to get coordinates later 
bp1 <- barplot(get(paste(gene,"_working_dataframe",sep="")),
    beside=T,
    ylim=c(0,50),         
    col=c(nicedarkblue,nicered),
    #lwd=2,
    yaxt="n",
    space=c(0,0.35),
    cex.names=1,
    col.axis=nicegrey,
    xlab="",
    cex.lab=1,
    col.lab=nicegrey,
    #border=c("white","white"),
    #border=c(nicedarkblue,nicered)
    )
}

bring_dfs_to_max_columns=function(gene,max_columns){
  con_pat <- c("Control","Patient")
  # This is a function to bring all dataframes to the same number of columns, creating blank columns with 0,0 entries
  # Finc the number of columns in the original dataframe
  N_original_col <- length(get(paste(gene,"_working_dataframe",sep="")))/2
  # Find the number of new columns to be added 
  N_new_col <- (max_columns-length(get(paste(gene,"_working_dataframe",sep="")))/2)
  if (N_new_col > 0){ # For any other than the already longest df 
    # Create the vector of blank names
    A2 <- c()
    this<-c(" "," ")
    for (i in 1:N_new_col){
      A2 <- c(A2,this)
      this <- gsub("^ ","  ",this)
    }
    # Create DF of blank names and patient status 
    temp_df <- data.frame("Status"=c(rep(con_pat,N_new_col)),"A2"=A2)
    # Merge with original table 
    temp_df<-cbind(get(paste(gene,"_working_dataframe",sep="")),table(temp_df))
    # Replace new entries with 0,0
    replacement<-c(0,0)
    for (i in seq(N_original_col+1,N_original_col+N_new_col)){
      temp_df[,i] <- replacement
    }
  return(temp_df)
  } else {
    return(get(paste(gene,"_working_dataframe",sep="")))
  }
}

plot_letter_label=function(x_min,x_max,y_min,y_max,letter){
  y_val <- y_min+1.15*(y_max-y_min)
  x_val <- x_min-0.17*(x_max-x_min)
  text(x_val,y_val,letter,col=nicegrey,adj=0,cex=1.2,xpd=NA,font=2)
}

analyse_data_function_manhattan_plot=function(gwasResults){
  gwasResults$CHR <- gwasResults$contig
  # remove chr 
  gwasResults$CHR <- gsub("chr","",gwasResults$CHR)
  # remove X, Y and MT (neccessary for qqman)
  gwasResults <- gwasResults[gwasResults$CHR != "X" & gwasResults$CHR != "Y",]
  gwasResults$CHR <- as.numeric(gsub("chr","",gwasResults$CHR))
  gwasResults$BP <- gwasResults$start
  gwasResults$P <- gwasResults$pvalue
  gwasResults$colour<- NA
  gwasResults$colour<- ifelse(gwasResults$CHR %% 2 == 0,nicedarkblue,nicered)
  gwasResults$SNP <- paste(gwasResults$contig,gwasResults$start,sep="_")
  return(gwasResults)
}

# Make function to create Manhattan and QQ plot 
manhattan_and_qq = function(dataset,title){
  dataset<-analyse_data_function_manhattan_plot(dataset)
  manhattan(dataset,suggestiveline=FALSE,genomewideline=FALSE,main =title ,col=c(nicered,nicedarkblue),ylim=c(0,20),col.main="grey20", col.lab="grey20", cex.lab=1, cex.names=1, border="grey20", col.axis="grey20", cex.axis=1 , chrlabs = c(1,"",3,"",5,"",7,"",9,"",11,"",13,"",15,"",17,"",19,"",21,""))
  abline(h=-log10(0.05/nrow(dataset)),b=0,col="grey20",lty=2,lwd=0.8)
  # Not happy with how qq plot looks so removing it
  # making blank plot for now to avoid resetting up plot space 
  blank_plot(0,10,0,10)
}

vector_of_datasets <-         c("EHV2",               "EHV3",               "GangSTR_Targeted",     "GangSTR_NonTargeted",    "tredparse",  "repeatseq",  "hipstr", "STRetch")
vector_of_software_names <-   c("ExpansionHunter_v2", "ExpansionHunter_v3", "GangSTR_Target_Mode",  "GangSTR_NonTarget_Mode", "Tredparse",  "RepeatSeq",  "HipSTR", "STRetch","exSTRa")
vector_of_software_colours <- c(nicered,              nicedarkblue,         niceorange,             nicelightblue,            darkgreen,    algae,         lavendar, hull,light_purple)
vector_of_datasets2 <-         c("EHV2",               "EHV3",               "GangSTR_Targeted",     "GangSTR_NonTargeted",    "tredparse",  "repeatseq",  "hipstr", "STRetch","exSTRa")

################################################################################
################################################################################

# Figure 3.5: Differential Repeat Coverage with Alternative Exome Targets 

################################################################################
################################################################################

########
# Create Colour Variables for beeswarm plot 
########

ALS_DOC$V2 <- "ALS"
ALS_DOC$Fill <- nicedarkblue_rgb_fade_light 
ALS_DOC$Outline <- nicedarkblue
RCSI_DOC$V2 <- "RCSI"
RCSI_DOC$Fill <- nicered_rgb_fade_light
RCSI_DOC$Outline <- nicered
Combined_DOC<- rbind(ALS_DOC,RCSI_DOC)

# Create beeswarm plot 
# This is need to get x positions later
beeswarm_variable <- beeswarm(V1 ~ V2, data=Combined_DOC,corral="wrap",ylim=c(0,150),pwcol = Outline,pwbg = Fill,pch = 21,lwd=2,axes=F,xlab="",ylab="",cex=2)

########
# Data wrangling for coverage comparisons plots 
########

Ratio_DF <- data.frame(Gene=sort(unique(Exonic_Repeats_Coordinates$V4)),Ratio=NA)
for (gene in sort(unique(Exonic_Repeats_Coordinates$V4))){
  start_pos <- Exonic_Repeats_Coordinates$V2[Exonic_Repeats_Coordinates$V4==gene]
  end_pos <- Exonic_Repeats_Coordinates$V3[Exonic_Repeats_Coordinates$V4==gene]
  RCSI_repeat_mean <-mean(RCSI_Repeat_Coverage_Mean$Mean_Coverage[RCSI_Repeat_Coverage_Mean$Position>start_pos & RCSI_Repeat_Coverage_Mean$Position<end_pos])
  ALS_repeat_mean <-mean(ALS_Repeat_Coverage_Mean$Mean_Coverage[ALS_Repeat_Coverage_Mean$Position>start_pos & ALS_Repeat_Coverage_Mean$Position<end_pos])
  Ratio_DF$Ratio[Ratio_DF$Gene==gene] <- ALS_repeat_mean/RCSI_repeat_mean
}
Overall_Ratio <- mean(ALS_DOC$V1) / mean(RCSI_DOC$V1)
# We have a point that is significantly higher than all others, going to manually adjust this and create an axis break
Ratio_DF$Ratio <- ifelse(Ratio_DF$Ratio>12,Ratio_DF$Ratio-19,Ratio_DF$Ratio)


####################
# Create Plot  
####################

# PDF = 210x297 - margins are 35,20,20,20 mm 
# PDF = 8.3x11.7 - margins are 1.38,0.79,0.79,0.79 inches
pdf(paste(output_directory,"Chapter3_Exonic_Expansions_Coverage_in_Exome_Data.pdf",sep=""),width=(8.3-1.38-0.8),height=(11.7-0.8-0.8))
par(mfrow=c(3,1),mar=c(5.5,8,3,1))# PDF = 210x297 - margins are 35,20,20,20 mm 


####################
# beeswarm plot (A)
####################

beeswarm(V1 ~ V2, data=Combined_DOC,corral="wrap",ylim=c(0,150),pwcol = Outline,pwbg = Fill,pch = 21,lwd=2,axes=F,xlab="",ylab="",cex=2)
#Plot central mean lines
x_pos=1
lines(x=c(x_pos-0.3,x_pos+0.3),y= c(mean(Combined_DOC$V1[Combined_DOC$V2=="ALS"]),mean(Combined_DOC$V1[Combined_DOC$V2=="ALS"])),col=nicedarkblue,lwd=1.5)
x_pos=2
lines(x=c(x_pos-0.3,x_pos+0.3),y= c(mean(Combined_DOC$V1[Combined_DOC$V2=="RCSI"]),mean(Combined_DOC$V1[Combined_DOC$V2=="RCSI"])),col=nicered,lwd=1.5)

# Plot Legend
points(x=c(2.2,2.2),y=c(130,140),bg=c(nicered,nicedarkblue),col="white",pch=21,cex=1.4,lwd=0)
text(2.25,130,"SeqCap",col=nicegrey,adj=0,cex=1,xpd=NA)
text(2.25,140,"SureSelect",col=nicegrey,adj=0,cex=1,xpd=NA)
# Plot X Axis 
axis(side=1,lwd=1,cex.axis=1,col.axis=nicegrey,col=nicegrey,at=c(1,2),labels=c("SureSelect (ALS Exome Samples)","SeqCap (Epilepsy Samples)"))
#Plot Y Axis
axis(side=2,lwd=1,las=1,cex.axis=1.2,col.axis=nicegrey,col=nicegrey)

# Plot Title
title(
  main="Overall Sample Coverage in Target Region", 
  xlab="Samples", 
  ylab="Mean Coverage (X)",
  col.main=nicegrey,
  cex.main=1.2,
  col.lab=nicegrey,
  cex.lab=1.2
  )

# Plot Letter Label
plot_letter_label_b(x_min=0,x_max=max(beeswarm_variable$x),y_min=0,y_max=150,"A")

####################
# Ratio Plot (B)
####################

blank_plot(xmin=0,xmax=8,ymin=0,ymax=11)
# Plot median line 
lines(x=c(0,8),y=c(Overall_Ratio,Overall_Ratio),col=niceorange,lwd=2)
# Plot points 
points(0:8,Ratio_DF$Ratio,cex=1.5,col=nicegrey,pch=21,bg=nicegrey)
# Plot Y axis 
axis(side=2,lwd=1,las=1,cex.axis=1.2,col.axis=nicegrey,col=nicegrey,at=0:11,labels=c(0:9,29,30))
axis.break(2,breakpos=9.5,breakcol=nicegrey)
# Plot X axis 
axis(side=1,lwd=1,las=1,cex.axis=1,col.axis=nicegrey,col=nicegrey,at=0:8,labels=sort(unique(Exonic_Repeats_Coordinates$V4)),font=3)
# Plot explanation text 
text(0,2.7,"Total Mean SureSelect Coverage / Total Mean SeqCap Coverage",col=niceorange,adj=0,cex=0.7,xpd=NA)
# Add Arrows 
arrows(x0=1,y0=Overall_Ratio+1,x1=1,y1=Overall_Ratio+2,length=0.1,angle=30,lwd=1.2,col=lavendar)
arrows(x0=1,y0=Overall_Ratio-1,x1=1,y1=Overall_Ratio-2,length=0.1,angle=30,lwd=1.2,col=lavendar)
# Add Arrow Text
text(1.2,Overall_Ratio+2,"SeqCap Below Expected Coverage",col=lavendar,adj=0,cex=0.7,xpd=NA)
text(1.2,Overall_Ratio-2,"SeqCap Above Expected Coverage",col=lavendar,adj=0,cex=0.7,xpd=NA)
title(
  main="Exome Target Kit Coverage Ratios", 
  xlab="", 
  ylab="SureSelect Coverage / SeqCap Coverage",
  col.main=nicegrey,
  cex.main=1.2,
  col.lab=nicegrey,
  cex.lab=1.2
  )
plot_letter_label_b(x_min=0,x_max=8,y_min=0,y_max=11,"B")
i=2
for (gene in unique(ALS_Repeat_Coverage_Mean$Gene)) {
  # Get Letter To Plot 
  i=i+1
  plot_letter=LETTERS[i]
  RCSI_Working<-RCSI_Repeat_Coverage_Mean[RCSI_Repeat_Coverage_Mean$Gene==gene,]
  ALS_Working<-ALS_Repeat_Coverage_Mean[ALS_Repeat_Coverage_Mean$Gene==gene,]
  # Define repeat motif
  repeat_motif=repeat_motifs[match(gene,gene_list_quotes)]
  #
  # Create Plot 
  #
  # Create blank plot 
  # plot height (1.1 x the combined heights of two plots )
  plot_height <- 1.1*(max(RCSI_Working$Mean_plus_SD[RCSI_Working$Gene==gene])+max(ALS_Working$Mean_plus_SD[ALS_Working$Gene==gene]))
  blank_plot(min(RCSI_Working$Position),min(RCSI_Working$Position)+max_offset,0,plot_height)
  # PLOT GENE AND VERTICAL LINES 
  # This is the positions of the target regions for RCSI exomes

  # Plot Regions that are targeted by exome kits 
  SC_x_starts<-SeqCap_Regions$V2[SeqCap_Regions$V4==gene]
  SC_x_ends<-SeqCap_Regions$V3[SeqCap_Regions$V4==gene]
  rect(xleft=SC_x_starts,ybottom=rep(-0.02*plot_height,length(SC_x_starts)),xright=SC_x_ends,ytop=rep(-0.01*plot_height,length(SC_x_starts)),border=NA,col=nicered,xpd=NA)
  SS_x_starts<-SureSelect_Regions$V2[SureSelect_Regions$V4==gene]
  SS_x_ends<-SureSelect_Regions$V3[SureSelect_Regions$V4==gene]
  rect(xleft=SS_x_starts,ybottom=rep(max(RCSI_Working$Mean_plus_SD) + 0.1*plot_height-0.02*plot_height,length(SS_x_starts)),xright=SS_x_ends,ytop=rep(max(RCSI_Working$Mean_plus_SD) + 0.1*plot_height-0.01*plot_height,length(SS_x_starts)),border=NA,col=nicedarkblue,xpd=NA)
  # Give labels to regions targeted by exome kits 
  text(x=max(c(RCSI_Working$Position,ALS_Working$Position)),y=-0.05*plot_height,"SeqCap Targets",col=nicered,adj=1,cex=0.7,xpd=NA)
  text(x=max(c(RCSI_Working$Position,ALS_Working$Position)),y=max(RCSI_Working$Mean_plus_SD) + 0.1*plot_height-0.05*plot_height,"SureSelect Targets",col=nicedarkblue,adj=1,cex=0.8,xpd=NA)
  # Get exon coordinates 
  this_gene_coordinates <- tail(head(gene_coordinates[gene_coordinates$name2==gene,],3),1)
  # 5' coordinate
  x_min<-this_gene_coordinates$txStart
  # 3' coordinate
  x_max<-this_gene_coordinates$txEnd
  # Find exon coordinates
  left<-as.numeric(unlist(strsplit(as.character(this_gene_coordinates$exonStarts), split=",")))
  right<-as.numeric(unlist(strsplit(as.character(this_gene_coordinates$exonEnds), split=",")))
  # Find coding exon coordinates
  coding_left<-left
  coding_right<-right
  remove_exons<-c()
  for (number in 1:length(coding_left)){
    # Find completely non coding exons
    if ((coding_left[number] < this_gene_coordinates$cdsStart & coding_right[number] < this_gene_coordinates$cdsStart) | (coding_left[number] > this_gene_coordinates$cdsEnd & coding_right[number]> this_gene_coordinates$cdsEnd)) {
      remove_exons<-append(remove_exons,number)
    }
    # Correct first exon
    if (coding_left[number]<this_gene_coordinates$cdsStart & coding_right[number]>this_gene_coordinates$cdsStart & coding_right[number]<this_gene_coordinates$cdsEnd){
      coding_left[number]<-this_gene_coordinates$cdsStart
    }
    # Correct last exon
    if (coding_left[number]>this_gene_coordinates$cdsStart & coding_right[number]>this_gene_coordinates$cdsEnd){
      coding_right[number] <- this_gene_coordinates$cdsEnd
    }
    # If only one exon 
    if (coding_left[number]<this_gene_coordinates$cdsStart & coding_right[number]>this_gene_coordinates$cdsEnd){
      coding_right[number] <- this_gene_coordinates$cdsEnd
      coding_left[number]<-this_gene_coordinates$cdsStart
    }
  }
  # Remove non coding exons 
  if (length(remove_exons)>0){
    coding_left<-coding_left[-remove_exons]
    coding_right<-coding_right[-remove_exons]
  }
  # Plot central line 
  lines(x=c(min(RCSI_Working$Position),max(RCSI_Working$Position)),y=c(-0.1*plot_height,-0.1*plot_height),col="grey",lty=1,lwd=1.2,xpd=NA)
  # Plot All Exons
  rect(xleft=left,ybottom=rep(-0.11*plot_height,length(left)),xright=right,ytop=rep(-0.09*plot_height,length(left)),border=NA,col="grey",xpd=NA)
  # Plot Coding Exons
  rect(xleft=coding_left,ybottom=rep(-0.125*plot_height,length(coding_left)),xright=coding_right,ytop=rep(-0.075*plot_height,length(coding_left)),border=NA,col="grey",xpd=NA)
  # Plot Repeat Region 
  repeat_start <- Exonic_Repeats_Coordinates$V2[Exonic_Repeats_Coordinates$V4==gene]
  repeat_end <- Exonic_Repeats_Coordinates$V3[Exonic_Repeats_Coordinates$V4==gene]
  rect(xleft=repeat_start,ybottom=-0.125*plot_height,xright=repeat_end,ytop=-0.075*plot_height,border=NA,col=niceorange,xpd=NA)
  # Plot blank squares to cover gene that went outside plot area 
  # Left 
  rect(xleft=min(ALS_Working$Position)-(max(ALS_Working$Position)-min(ALS_Working$Position)),ybottom=-0.125*plot_height,xright=min(ALS_Working$Position),ytop=-0.075*plot_height,border=NA,col="white",xpd=NA)
  # Right
  rect(xleft=max(ALS_Working$Position),ybottom=-0.125*plot_height,xright=max(ALS_Working$Position)+(max(ALS_Working$Position)-min(ALS_Working$Position)),ytop=-0.075*plot_height,border=NA,col="white",xpd=NA)
  
  # RCSI PLOT
  # Plot Mean Polygon
  polygon(c(min(RCSI_Working$Position),RCSI_Working$Position,max(RCSI_Working$Position)),c(0,RCSI_Working$Mean_Coverage,0),col=nicered_rgb_fade_heavy,border=NA)
  #Plot SD Polygon
  polygon(c(RCSI_Working$Position,rev(RCSI_Working$Position),min(RCSI_Working$Position)),c(RCSI_Working$Mean_plus_SD,rev(RCSI_Working$Mean_minus_SD),RCSI_Working$Mean_plus_SD[RCSI_Working$Position==min(RCSI_Working$Position)]),col=nicered_rgb_fade_light,border=NA)

  # ALS PLOT
  # Need to add the height of the RCSI to place the plot correctly 
  ALS_Working$Mean_plus_SD <- ALS_Working$Mean_plus_SD + max(RCSI_Working$Mean_plus_SD) + 0.1*plot_height
  ALS_Working$Mean_minus_SD <- ALS_Working$Mean_minus_SD + max(RCSI_Working$Mean_plus_SD) + 0.1*plot_height
  ALS_Working$Mean_Coverage <- ALS_Working$Mean_Coverage + max(RCSI_Working$Mean_plus_SD) + 0.1*plot_height

  # Plot Mean Polygon
  polygon(c(min(ALS_Working$Position),ALS_Working$Position,max(ALS_Working$Position)),c(max(RCSI_Working$Mean_plus_SD) + 0.1*plot_height,ALS_Working$Mean_Coverage,max(RCSI_Working$Mean_plus_SD) + 0.1*plot_height),col=nicedarkblue_rgb_fade_heavy,border=NA)
  #Plot SD Polygon
  polygon(c(ALS_Working$Position,rev(ALS_Working$Position),min(ALS_Working$Position)),c(ALS_Working$Mean_plus_SD,rev(ALS_Working$Mean_minus_SD),ALS_Working$Mean_plus_SD[ALS_Working$Position==min(ALS_Working$Position)]),col=nicedarkblue_rgb_fade_light,border=NA)
  # Plot left legend lines 
  # RCSI
  lines(x=c(min(c(RCSI_Working$Position,ALS_Working$Position))-0.01*(max(c(RCSI_Working$Position,ALS_Working$Position))-min(c(RCSI_Working$Position,ALS_Working$Position))),min(c(RCSI_Working$Position,ALS_Working$Position))-0.01*(max(c(RCSI_Working$Position,ALS_Working$Position))-min(c(RCSI_Working$Position,ALS_Working$Position)))),y=c(0,max(RCSI_Working$Mean_plus_SD)),col=nicegrey,lty=1,lwd=1,xpd=NA)
  lines(x=c(min(c(RCSI_Working$Position,ALS_Working$Position))-0.01*(max(c(RCSI_Working$Position,ALS_Working$Position))-min(c(RCSI_Working$Position,ALS_Working$Position))),min(c(RCSI_Working$Position,ALS_Working$Position))-0.02*(max(c(RCSI_Working$Position,ALS_Working$Position))-min(c(RCSI_Working$Position,ALS_Working$Position)))),y=c(0,0),col=nicegrey,lty=1,lwd=1,xpd=NA)
  lines(x=c(min(c(RCSI_Working$Position,ALS_Working$Position))-0.01*(max(c(RCSI_Working$Position,ALS_Working$Position))-min(c(RCSI_Working$Position,ALS_Working$Position))),min(c(RCSI_Working$Position,ALS_Working$Position))-0.02*(max(c(RCSI_Working$Position,ALS_Working$Position))-min(c(RCSI_Working$Position,ALS_Working$Position)))),y=c(max(RCSI_Working$Mean_plus_SD),max(RCSI_Working$Mean_plus_SD)),col=nicegrey,lty=1,lwd=1,xpd=NA)

  # ALS
  lines(x=c(min(c(RCSI_Working$Position,ALS_Working$Position))-0.01*(max(c(RCSI_Working$Position,ALS_Working$Position))-min(c(RCSI_Working$Position,ALS_Working$Position))),min(c(RCSI_Working$Position,ALS_Working$Position))-0.01*(max(c(RCSI_Working$Position,ALS_Working$Position))-min(c(RCSI_Working$Position,ALS_Working$Position)))),y=c(max(RCSI_Working$Mean_plus_SD) + 0.1*plot_height,max(ALS_Working$Mean_plus_SD)),col=nicegrey,lty=1,lwd=1,xpd=NA)
  lines(x=c(min(c(RCSI_Working$Position,ALS_Working$Position))-0.01*(max(c(RCSI_Working$Position,ALS_Working$Position))-min(c(RCSI_Working$Position,ALS_Working$Position))),min(c(RCSI_Working$Position,ALS_Working$Position))-0.02*(max(c(RCSI_Working$Position,ALS_Working$Position))-min(c(RCSI_Working$Position,ALS_Working$Position)))),y=c(max(RCSI_Working$Mean_plus_SD) + 0.1*plot_height,max(RCSI_Working$Mean_plus_SD) + 0.1*plot_height),col=nicegrey,lty=1,lwd=1,xpd=NA)
  lines(x=c(min(c(RCSI_Working$Position,ALS_Working$Position))-0.01*(max(c(RCSI_Working$Position,ALS_Working$Position))-min(c(RCSI_Working$Position,ALS_Working$Position))),min(c(RCSI_Working$Position,ALS_Working$Position))-0.02*(max(c(RCSI_Working$Position,ALS_Working$Position))-min(c(RCSI_Working$Position,ALS_Working$Position)))),y=c(max(ALS_Working$Mean_plus_SD),max(ALS_Working$Mean_plus_SD)),col=nicegrey,lty=1,lwd=1,xpd=NA)

  # Add coverage labels 
  # RCSI
  text(min(c(RCSI_Working$Position,ALS_Working$Position))-0.03*(max(c(RCSI_Working$Position,ALS_Working$Position))-min(c(RCSI_Working$Position,ALS_Working$Position))),0,"0",col=nicegrey,adj=1,cex=1.2,xpd=NA)
  text(min(c(RCSI_Working$Position,ALS_Working$Position))-0.03*(max(c(RCSI_Working$Position,ALS_Working$Position))-min(c(RCSI_Working$Position,ALS_Working$Position))),max(RCSI_Working$Mean_plus_SD),round(max(RCSI_Working$Mean_plus_SD),0),col=nicegrey,adj=1,cex=1.2,xpd=NA)
  # ALS
  text(min(c(RCSI_Working$Position,ALS_Working$Position))-0.03*(max(c(RCSI_Working$Position,ALS_Working$Position))-min(c(RCSI_Working$Position,ALS_Working$Position))),max(RCSI_Working$Mean_plus_SD) + 0.1*plot_height,"0",col=nicegrey,adj=1,cex=1.2,xpd=NA)
  text(min(c(RCSI_Working$Position,ALS_Working$Position))-0.03*(max(c(RCSI_Working$Position,ALS_Working$Position))-min(c(RCSI_Working$Position,ALS_Working$Position))),max(ALS_Working$Mean_plus_SD),round(max(ALS_Working$Mean_plus_SD),0),col=nicegrey,adj=1,cex=1.2,xpd=NA)

  # Add y axis labels 
  text(min(c(RCSI_Working$Position,ALS_Working$Position))-0.15*(max(c(RCSI_Working$Position,ALS_Working$Position))-min(c(RCSI_Working$Position,ALS_Working$Position))),min(RCSI_Working$Mean_plus_SD)+0.5*(max(RCSI_Working$Mean_plus_SD)-min(RCSI_Working$Mean_plus_SD)),"SeqCap \n Coverage (X)",col=nicegrey,adj=0.5,cex=1.2,xpd=NA,font=3)
  text(x=min(c(RCSI_Working$Position,ALS_Working$Position))-0.15*(max(c(RCSI_Working$Position,ALS_Working$Position))-min(c(RCSI_Working$Position,ALS_Working$Position))),y=max(RCSI_Working$Mean_plus_SD)+0.1*(plot_height)+0.5*(max(ALS_Working$Mean_plus_SD)-min(ALS_Working$Mean_plus_SD)),"SureSelect \n Coverage (X)",col=nicegrey,adj=0.5,cex=1.2,xpd=NA,font=3)

  # Plot vertical dashed lines 

  lines(x=c(repeat_start,repeat_start),y=c(-0.075*plot_height,plot_height),col=niceorange,lty=2,lwd=1,xpd=NA)
  lines(x=c(repeat_end,repeat_end),y=c(-0.075*plot_height,plot_height),col=niceorange,lty=2,lwd=1,xpd=NA)

  # Write gene name
  text(min(c(RCSI_Working$Position,ALS_Working$Position))-0.03*(max(c(RCSI_Working$Position,ALS_Working$Position))-min(c(RCSI_Working$Position,ALS_Working$Position))),-0.1*plot_height,gene,col=nicegrey,adj=1,cex=1.2,xpd=NA,font=3)
  # Write repeat motif 
  text(repeat_start+(repeat_end-repeat_start)/2,-0.16*plot_height,bquote(paste(.(repeat_motif)['N'])),col=niceorange,adj=0.5,cex=1.2,xpd=NA)
  plot_letter_label_b(x_min=min(RCSI_Working$Position),x_max=max(RCSI_Working$Position),y_min=0,y_max=plot_height,plot_letter)
}

dev.off()


################################################################################
################################################################################

# Table 3.4 : RMSD, Sensitivity, Specificity of in silico Genotyping Tools Relative to PCR Data and Between WES and WGS : Print Statistics

################################################################################
################################################################################

# Different tools start counting the alleles at certain loci from a different point to PCR, using a correction to bring them in concurrence with the PCR 
#
EHV2_PCR_Gene_Correction_Gene <-c(  "AR", "ATN1", "ATXN3","ATXN7",  "GIPC1","HTT","NOP56","PPP2R2B","TBP")
EHV2_PCR_Gene_Correction_Value <-c(  -1,    -4,     3,     -5,       -2,      -3,     2,    -2,        -1)
EHV3_PCR_Gene_Correction_Gene <-c(  "AR", "ATN1", "ATXN3","ATXN7",  "GIPC1","HTT","NOP56","PPP2R2B","TBP")
EHV3_PCR_Gene_Correction_Value <-c(  -1,    -4,     3,     -5,             -2,      -3,     2,    -2,        -1)
GangSTR_NonTargeted_PCR_Gene_Correction_Gene <-c(  "AR", "ATN1", "ATXN3","ATXN7",  "GIPC1","HTT","NOP56","PPP2R2B","TBP")
GangSTR_NonTargeted_PCR_Gene_Correction_Value <-c(  0,    0,     6,     -5,         -2,      -3,     4,    -2,        -1)
GangSTR_Targeted_PCR_Gene_Correction_Gene <-c(  "AR", "ATN1", "ATXN3","ATXN7",  "GIPC1","HTT","NOP56","PPP2R2B","TBP")
GangSTR_Targeted_PCR_Gene_Correction_Value <-c(  0,    0,     6,     -5,            -2,      -3,     4,    -2,        -1)
hipstr_PCR_Gene_Correction_Gene <-c(  "AR", "ATN1", "ATXN3","ATXN7",  "GIPC1","HTT","NOP56","PPP2R2B","TBP","JPH3")
hipstr_PCR_Gene_Correction_Value <-c(  3,    -4,     3,     -5,            -6,  0,     1,    -3,        -1,-2)
tredparse_PCR_Gene_Correction_Gene <- c("AR", "ATN1", "ATXN3","ATXN7",  "GIPC1","HTT","NOP56","PPP2R2B","TBP")
tredparse_PCR_Gene_Correction_Value <- c(  1,    0,     6,     -5,       -2,      -3,     1,    -2, -2     )
STRetch_PCR_Gene_Correction_Gene <-c(  "AR", "ATN1", "ATXN3","ATXN7",  "GIPC1","HTT","NOP56","PPP2R2B","TBP")
STRetch_PCR_Gene_Correction_Value <-c(  -1,    -4,     3,     -5,             -2,      -3,     2,    -2,        -1)
repeatseq_PCR_Gene_Correction_Gene <-c(  "AR", "ATN1", "ATXN3","ATXN7", "GIPC1","HTT","NOP56","PPP2R2B","TBP")
repeatseq_PCR_Gene_Correction_Value <-c(  0,    -4,     3,     0,      -2,      -7,     2,    -2,        -1)

# All Genes Plot 
for (primary_dataset in vector_of_datasets2){
  # STEP 1: 
  # Find Called and True Positives for C9orf72 (note: this could be done a lot more efficiently but reusing code from elsewhere)
  #
  Working_Dataset <-get(primary_dataset)[get(primary_dataset)$Gene=="C9orf72" & get(primary_dataset)$Dataset=="ALS_WGS_PCR_FREE",]
  # Merge 
  Working_Dataset <- merge(Working_Dataset,ALS_CC_STATUS,all.x=T,all.y=F)
  Working_Dataset$C9orf72_Status <- ifelse(is.na(Working_Dataset$C9orf72_Status)==T,0,Working_Dataset$C9orf72_Status)
  Working_Dataset$PCR_Expanded <- NA
  Working_Dataset$PCR_Expanded<- ifelse(Working_Dataset$C9orf72_Status==1,"TRUE_POSITIVE","TRUE_NEGATIVE")
  Working_Dataset$Software_Expanded <- NA
  if (primary_dataset!="exSTRa" & primary_dataset!="STRetch"){
    Working_Dataset$Software_Expanded<- ifelse(Working_Dataset$A2 >=30, "GENOTYPING_POSITIVE","GENOTYPING_NEGATIVE")
  } else if (primary_dataset=="STRetch"){
    Working_Dataset$Software_Expanded<- ifelse(Working_Dataset$Adjusted_PValue <0.05, "GENOTYPING_POSITIVE","GENOTYPING_NEGATIVE")
  }else {
    Working_Dataset$Software_Expanded<- ifelse(Working_Dataset$P_Value <=(0.05/34), "GENOTYPING_POSITIVE","GENOTYPING_NEGATIVE")
  }
  Working_Dataset$Final_Call <- "TRUE_NEGATIVE"
  # Get True Positives (i.e. PCR and genotyping positive)
  Working_Dataset$Final_Call <- ifelse(Working_Dataset$PCR_Expanded=="TRUE_POSITIVE" & Working_Dataset$Software_Expanded=="GENOTYPING_POSITIVE","TRUE_POSITIVE",Working_Dataset$Final_Call)
  # Get False Positives 
  Working_Dataset$Final_Call <- ifelse(Working_Dataset$PCR_Expanded=="TRUE_NEGATIVE" & Working_Dataset$Software_Expanded=="GENOTYPING_POSITIVE","FALSE_POSITIVE",Working_Dataset$Final_Call)
  # Get False Negative
  Working_Dataset$Final_Call <- ifelse(Working_Dataset$PCR_Expanded=="TRUE_POSITIVE" & Working_Dataset$Software_Expanded=="GENOTYPING_NEGATIVE","FALSE_NEGATIVE",Working_Dataset$Final_Call)
  final_calls <- Working_Dataset$Final_Call
  sensitivity <- length(final_calls[final_calls=="TRUE_POSITIVE"])*100/(length(final_calls[final_calls=="TRUE_POSITIVE"])+length(final_calls[final_calls=="FALSE_NEGATIVE"]))
  specificity <- length(final_calls[final_calls=="TRUE_NEGATIVE"])*100/(length(final_calls[final_calls=="TRUE_NEGATIVE"])+length(final_calls[final_calls=="FALSE_POSITIVE"]))
  print(paste(primary_dataset,": C9orf72 Sensitivity: ",round(sensitivity,2),sep=""))
  print(paste(primary_dataset,": C9orf72 Specificity: ",round(specificity,2),sep=""))

  # STEP 2: 
  # All genes except C9orf72 
  #
  # Filter to PCR FREE genotypes
  Working_Dataset <- get(primary_dataset)[get(primary_dataset)$Dataset=="ALS_WGS_PCR_FREE",]
  # Merge 
  Working_Dataset <- merge(Working_Dataset,PCR_Genotypes,all.x=T,all.y=F)
  Working_Dataset <- Working_Dataset[!is.na(Working_Dataset$PCR_A1) & !is.na(Working_Dataset$PCR_A2),]
  # Currently a problem with the SCA1 & DIP2B PCR genotypes and C9orf72 has been looked at elsewhere
  Working_Dataset <- Working_Dataset[Working_Dataset$Gene!="ATXN1" & Working_Dataset$Gene!="C9orf72" & Working_Dataset$Gene!="DIP2B",]
  # Remove TNRC6A and NIPA1 as they don't have defined cutoffs (see methods)
  Working_Dataset <- Working_Dataset[Working_Dataset$Gene!="TNRC6A" & Working_Dataset$Gene!="NIPA1",]
  # Find Called and True Positives for All Genes Except C9
  Working_Dataset$PCR_Expanded <- NA 
  Working_Dataset$Software_Expanded <- NA 
  if (primary_dataset!="exSTRa" & primary_dataset != "STRetch"){
    # Apply Correction
    for (entry in get(paste(primary_dataset,"_PCR_Gene_Correction_Gene",sep=""))){
      Working_Dataset$PCR_A2[Working_Dataset$Gene==entry] <- Working_Dataset$PCR_A2[Working_Dataset$Gene==entry]-as.numeric(get(paste(primary_dataset,"_PCR_Gene_Correction_Value",sep=""))[match(entry,get(paste(primary_dataset,"_PCR_Gene_Correction_Gene",sep="")))])
    }
    for (entry in unique(Working_Dataset$Gene)){
      Working_Dataset$PCR_Expanded[Working_Dataset$Gene==entry & Working_Dataset$PCR_A2 > cutoff_cutoffs[match(entry,cutoff_genes)]] <- "TRUE_POSITIVE"
      Working_Dataset$PCR_Expanded[Working_Dataset$Gene==entry & Working_Dataset$PCR_A2 <= cutoff_cutoffs[match(entry,cutoff_genes)]] <- "TRUE_NEGATIVE"   
      Working_Dataset$Software_Expanded[Working_Dataset$Gene==entry & Working_Dataset$A2 > cutoff_cutoffs[match(entry,cutoff_genes)]] <- "GENOTYPING_POSITIVE"
      Working_Dataset$Software_Expanded[Working_Dataset$Gene==entry & Working_Dataset$A2 <= cutoff_cutoffs[match(entry,cutoff_genes)]] <- "GENOTYPING_NEGATIVE"    
    }
  } else if (primary_dataset != "STRetch"){
    for (entry in unique(Working_Dataset$Gene)){
      Working_Dataset$PCR_Expanded[Working_Dataset$Gene==entry & Working_Dataset$PCR_A2 > cutoff_cutoffs[match(entry,cutoff_genes)]] <- "TRUE_POSITIVE"
      Working_Dataset$PCR_Expanded[Working_Dataset$Gene==entry & Working_Dataset$PCR_A2 <= cutoff_cutoffs[match(entry,cutoff_genes)]] <- "TRUE_NEGATIVE"
    }
    Working_Dataset$Software_Expanded[Working_Dataset$P_Value<(0.05/34)]  <- "GENOTYPING_POSITIVE"
    Working_Dataset$Software_Expanded[Working_Dataset$P_Value>=(0.05/34)]  <- "GENOTYPING_NEGATIVE"
  } else if (primary_dataset == "STRetch"){
    Working_Dataset <- get(primary_dataset)[get(primary_dataset)$Dataset=="ALS_WGS_PCR_FREE",]
    # For STRetch want to find the full list of genes which overlap with the PCR genotyping 
    # Find all genes in PCR
    genes_in_pcr <- unique(PCR_Genotypes$Gene)
    # Find genes that overlap between PCR and STRetch calls 
    STRetch_pcr_overlap_genes <- STRetch_Full_Gene_List[STRetch_Full_Gene_List %in% genes_in_pcr]
    # Filter to genes that STRetch has studied 
    PCR_Genotypes_Working <- PCR_Genotypes[(PCR_Genotypes$Gene %in% STRetch_pcr_overlap_genes),]
    Working_Dataset <- merge(PCR_Genotypes_Working,Working_Dataset,all.x=T,all.y=F)
    # Find true+false positives and negatives 
    Working_Dataset$PCR_Expanded <- NA 
    Working_Dataset$Software_Expanded <- "GENOTYPING_NEGATIVE" 
    for (entry in unique(Working_Dataset$Gene)){
      Working_Dataset$PCR_Expanded[Working_Dataset$Gene==entry & Working_Dataset$PCR_A2 > cutoff_cutoffs[match(entry,cutoff_genes)]] <- "TRUE_POSITIVE"
      Working_Dataset$PCR_Expanded[Working_Dataset$Gene==entry & Working_Dataset$PCR_A2 <= cutoff_cutoffs[match(entry,cutoff_genes)]] <- "TRUE_NEGATIVE"   
      Working_Dataset$Software_Expanded[Working_Dataset$Gene==entry & Working_Dataset$Adjusted_PValue < 0.05  & !is.na(Working_Dataset$Adjusted_PValue)] <- "GENOTYPING_POSITIVE"
    }
  }
  Working_Dataset$Final_Call <- "TRUE_NEGATIVE"
  # Get True Positives (i.e. PCR and genotyping positive)
  Working_Dataset$Final_Call <- ifelse(Working_Dataset$PCR_Expanded=="TRUE_POSITIVE" & Working_Dataset$Software_Expanded=="GENOTYPING_POSITIVE","TRUE_POSITIVE",Working_Dataset$Final_Call)
  # Get False Positives 
  Working_Dataset$Final_Call <- ifelse(Working_Dataset$PCR_Expanded=="TRUE_NEGATIVE" & Working_Dataset$Software_Expanded=="GENOTYPING_POSITIVE","FALSE_POSITIVE",Working_Dataset$Final_Call)
  # Get False Negative
  Working_Dataset$Final_Call <- ifelse(Working_Dataset$PCR_Expanded=="TRUE_POSITIVE" & Working_Dataset$Software_Expanded=="GENOTYPING_NEGATIVE","FALSE_NEGATIVE",Working_Dataset$Final_Call)
  # Add Final Calls to C9 calls 
  final_calls <- c(final_calls,Working_Dataset$Final_Call)
  sensitivity <- length(final_calls[final_calls=="TRUE_POSITIVE"])*100/(length(final_calls[final_calls=="TRUE_POSITIVE"])+length(final_calls[final_calls=="FALSE_NEGATIVE"]))
  specificity <- length(final_calls[final_calls=="TRUE_NEGATIVE"])*100/(length(final_calls[final_calls=="TRUE_NEGATIVE"])+length(final_calls[final_calls=="FALSE_POSITIVE"]))
  print(paste(primary_dataset,": All Gene Sensitivity: ",round(sensitivity,2),sep=""))
  print(paste(primary_dataset,": All Gene Specificity: ",round(specificity,2),sep=""))
  # Calculate RMSD - per gene per software and report rmsd + CI 
  if (primary_dataset=="STRetch"){
    Working_Dataset<-Working_Dataset[!is.na(Working_Dataset$A2) & !is.na(Working_Dataset$PCR_A2),]
    EH_alleles <-c(Working_Dataset$A2)
    PCR_alleles <- c(Working_Dataset$PCR_A2)
  } else if (primary_dataset != "exSTRa"){
    EH_alleles <-c(Working_Dataset$A1,Working_Dataset$A2)
    PCR_alleles <- c(Working_Dataset$PCR_A1,Working_Dataset$PCR_A2)
  } else {
    EH_alleles <- NA
    PCR_alleles <- NA
  }
  difference <- c(PCR_alleles-EH_alleles)
  rmsd <- (sum(difference^2)/length(difference))^(1/2)
  print(paste(primary_dataset,": RMSD: ",rmsd))
}
print("Sensitivity/Specificity/RMSD Calculations Complete")
 
 


################################################################################
################################################################################

# Figure 3.6 + Supplementary figures S3.XXX: in silico genotyping of the C9orf72 Repeat Expansion 

################################################################################
################################################################################

#########################################
# C9orf72 Plot - Compare software ability to correctly genotype the C9orf72 repeat expansion 
#########################################
i=1
# PDF = 210x297 - margins are 35,20,20,20 mm 
# PDF = 8.3x11.7 - margins are 1.38,0.79,0.79,0.79 inches
pdf(paste(output_directory,"Chapter3_C9orf72_Benchmarking.pdf",sep=""),width=(11.7-0.8-0.8),height=(8.3-1.38-0.8))
par(mfrow=c(3,1),mar=c(3.5,7.5,2.5,1.5))
for (primary_dataset in vector_of_datasets){
  software=vector_of_software_names[match(primary_dataset,vector_of_datasets)]
  # Call Dataset 
  Working_Dataset <-get(primary_dataset)[get(primary_dataset)$Gene=="C9orf72" & get(primary_dataset)$Dataset=="ALS_WGS_PCR_FREE",]
  # Merge 
  Working_Dataset <- merge(Working_Dataset,ALS_CC_STATUS,all.x=T,all.y=F)
  Working_Dataset$C9orf72_Status <- ifelse(is.na(Working_Dataset$C9orf72_Status)==T,0,Working_Dataset$C9orf72_Status)
  Working_Dataset$CC_C9 <- ifelse(Working_Dataset$C9orf72_Status==1,"C9_Positive_Patient",ifelse(Working_Dataset$CC==1,"C9_Negative_Patient","Control"))
  # Calculate statistics (don't move this and leave as EHV3 regardless of current software )
  if(primary_dataset != "STRetch"){
    # total 408 cases and controls 
    total_sample_count=408
  } else {
    # STRetch only gives output for cases 
    total_sample_count=272
  }
  ###
  # Calculate false pos / neg
  ###
  Percentage_of_all_samples_genotyped=nrow(Working_Dataset)*100/total_sample_count
  Positive_Samples<- pmid_vs_cc_status$Sample_well[pmid_vs_cc_status$c9orf72_status_clean==1 & !is.na(pmid_vs_cc_status$c9orf72_status_clean) & (pmid_vs_cc_status$Sample_well %in% unique(EHV3$Sample[EHV3$Dataset=="ALS_WGS_PCR_FREE"]))]
  Negative_Samples <- pmid_vs_cc_status$Sample_well[!(pmid_vs_cc_status$Sample_well %in% Positive_Samples)]
  positive_count<-length(Positive_Samples) # leave this as EHV3 it will be right for all
  negative_count<-total_sample_count - positive_count
  positives_genotyped_count <- nrow(Working_Dataset[Working_Dataset$Sample %in% Positive_Samples & !is.na(Working_Dataset$A2) & Working_Dataset$CC_C9=="C9_Positive_Patient",])
  negatives_genotyped_count <- nrow(Working_Dataset[Working_Dataset$Sample %in% Negative_Samples & !is.na(Working_Dataset$A2),])
  Percentage_of_positive_samples_genotyped=positives_genotyped_count*100/positive_count
  if(primary_dataset != "STRetch"){
    Sensitivity=nrow(Working_Dataset[Working_Dataset$A2>29 & !is.na(Working_Dataset$A2) & Working_Dataset$CC_C9=="C9_Positive_Patient",])*100/positives_genotyped_count
    Specificity=nrow(Working_Dataset[Working_Dataset$A2<30 & !is.na(Working_Dataset$A2) & (Working_Dataset$CC_C9=="C9_Negative_Patient" | Working_Dataset$CC_C9=="Control"),])*100/negatives_genotyped_count
  } else if (primary_dataset=="STRetch"){
    Sensitivity=nrow(Working_Dataset[Working_Dataset$Adjusted_PValue<0.05 & !is.na(Working_Dataset$A2) & Working_Dataset$CC_C9=="C9_Positive_Patient",])*100/positives_genotyped_count
    Specificity=nrow(Working_Dataset[Working_Dataset$Adjusted_PValue>=0.05 & !is.na(Working_Dataset$A2) & (Working_Dataset$CC_C9=="C9_Negative_Patient" | Working_Dataset$CC_C9=="Control"),])*100/negatives_genotyped_count
  }
  # All predicted repeats above 30 are being grouped as >30 
  Working_Dataset <- subset(Working_Dataset,select=c(A2,CC_C9))
  Working_Dataset$A2 <- as.numeric(Working_Dataset$A2)
  Working_Dataset <- Working_Dataset[,c(2,1)]
  C9orf72_predicitons_above_30 <- sort(Working_Dataset$A2[Working_Dataset$A2>=30])
  C9orf72_predicitons_above_30_not_sort <- Working_Dataset$A2[Working_Dataset$A2>=30]
  C9orf72_predicitons_above_30_colours_not_sort <- Working_Dataset$CC_C9[Working_Dataset$A2>=30]
  C9orf72_predicitons_above_30_colours <-C9orf72_predicitons_above_30_colours_not_sort[match(C9orf72_predicitons_above_30_not_sort,C9orf72_predicitons_above_30)]
  C9orf72_predicitons_above_30_colours <- gsub("Control",nicedarkblue,C9orf72_predicitons_above_30_colours)
  C9orf72_predicitons_above_30_colours <- gsub("C9_Negative_Patient",niceorange,C9orf72_predicitons_above_30_colours)
  C9orf72_predicitons_above_30_colours <- gsub("C9_Positive_Patient",nicered,C9orf72_predicitons_above_30_colours)
  Working_Dataset$A2 <- ifelse(Working_Dataset$A2>30,30,Working_Dataset$A2)
  # Adding columns so all numbers between 0 and 30 are represented then adding extra columns to write things beside plot (31-40)
  Working_Dataset$CC_C9 <- factor(Working_Dataset$CC_C9,levels=c("C9_Positive_Patient","C9_Negative_Patient","Control"))
  Working_Dataset$A2 <- factor(Working_Dataset$A2,levels=1:40)
  Working_Dataset<- table(Working_Dataset)
  Working_Dataset["C9_Negative_Patient",] <- Working_Dataset["C9_Negative_Patient",]*100/sum(Working_Dataset["C9_Negative_Patient",])
  Working_Dataset["C9_Positive_Patient",] <- Working_Dataset["C9_Positive_Patient",]*100/sum(Working_Dataset["C9_Positive_Patient",])
  Working_Dataset["Control",] <- Working_Dataset["Control",]*100/sum(Working_Dataset["Control",])
  # This controls the thickness of bar borders 
  opar <- par(lwd = 0.01)
  bp <- barplot(Working_Dataset,
    beside=T,
    ylim=c(0,100),         
    col=c(niceorange,nicered,nicedarkblue),
    yaxt="n",
    xaxt="n",
    space=c(0,0.6),
    cex.names=0.6,
    col.axis=nicegrey,
    xlab="",
    cex.lab=0.1,
    col.lab=nicegrey,
    border=NA
  )
  text(x=1,y=110,bquote(bold(.(software) ~ ": " ~ italic(.("C9orf72")) ~ " Allele Predicition")),adj=0,cex=1,col=nicegrey,xpd=TRUE,font=2)
  # Grid Lines:  only want to plot these over the original distance of the plot, not the fake distance that keeps all the bars uniform 
  lines(x=c(-5,max(bp[,30])),y=c(25,25),lty=2,lwd=0.8,col="grey")
  lines(x=c(-5,max(bp[,30])),y=c(50,50),lty=2,lwd=0.8,col="grey")
  lines(x=c(-5,max(bp[,30])),y=c(75,75),lty=2,lwd=0.8,col="grey")
  lines(x=c(-5,max(bp[,30])),y=c(100,100),lty=2,lwd=0.8,col="grey")
  # Replot bars over grid lines
  barplot(Working_Dataset,
    beside=T,
    ylim=c(0,100),         
    col=c(nicered,niceorange,nicedarkblue),
    yaxt="n",
    xaxt="n",
    space=c(0,0.6),
    cex.names=0.6,
    col.axis=nicegrey,
    xlab="",
    cex.lab=0.1,
    col.lab=nicegrey,
    lwd=0.01,
    add=T,
    border=NA
  )
  # white line at 0 to hide 0 values (that otherwise plot as tiny thin lines)
  abline(h=-0.3,lty=1,lwd=1.2,col="white",xpd=NA)
  # Title for y-axis 1 
  title(ylab="Software Predicted Alleles (%)",cex.lab=1,col.lab=nicegrey,line=3.5,font=3)
  # Plot axis
  axis(2,at=c(0,25,50,75,100),lab=c(0,25,50,75,100),las=1,lwd=1,col=nicegrey,col.axis=nicegrey,cex.axis=1)
  # Plot "repeats:"
  axis(1,at=-1,lab="repeats: ",las=1,col=nicegrey,col.axis=nicegrey,cex.axis=1,hadj=1,line=NA,tick=F)
  # Plot Coloured Points
  points(x=c(1,1,1),y=c(93,82,71),col=c(niceorange,nicered,nicedarkblue),bg=c(niceorange,nicered,nicedarkblue),pch=21,cex=1,lwd=0.6,xpd=NA)
  # Plot Case/Control Labels
  text(x=c(2,2,2),y=c(93,82,71),c(expression(paste(italic("C9orf72")," Negative Patients",sep="")),expression(paste(italic("C9orf72")," Positive Patients",sep="")),"Controls"),col=nicegrey,adj=0,cex=0.8,xpd=NA)
  # Repeat count label
  # plot lines above labels 
  for (number in 1:30){
    if (number==30){
      label=">29"
    } else {
      label=number
    }
    text(x=bp[,number][2],y=-16,label,col=nicegrey,adj=0.5,cex=0.8,xpd=NA)
    lines(x=c(bp[,number][1],bp[,number][3]),y=c(-12,-12),lty=1,col=nicegrey,lwd=0.5,xpd=NA)
  }
  # Plot Letter 
  if (i <= 26){
    plot_letter=LETTERS[i]
  } else {
    j=i
    j=j-26
    plot_letter=paste(LETTERS[j],".2",sep="")
  }
  plot_letter_label(min(bp),max(bp)*0.9,0,100,plot_letter)
  # Information to right of plot 
  # Allele lengths > 29
  if (length(C9orf72_predicitons_above_30)>0){
    #text(x=bp[,31][1],y=10+as.numeric(Working_Dataset[,30][1]),c("Predicted Repeat Lengths >29"),col=nicegrey,adj=0,cex=0.8,xpd=NA)
    text(x=bp[,31][1],y=90,c("Predicted Repeat Lengths >29"),col=nicegrey,adj=0,cex=0.8,xpd=NA)
    height=80
    allele_pos<-as.numeric(bp[,31][1])
    k=1
    for (entry in rev(sort(C9orf72_predicitons_above_30))){
      if (height <0){
        height=80
        allele_pos<-as.numeric(bp[,33][1])-0.02*bp[,40][1]
        lines(x=c(as.numeric(bp[,32][1]),as.numeric(bp[,32][1])),y=c(height,0),col=nicegrey,lwd=0.6,xpd=NA)
      }
      text(x=allele_pos,y=height,entry,col=C9orf72_predicitons_above_30_colours[k],adj=0,cex=0.6,xpd=NA)
      height=height-7
      k=k+1
    }
  }
  # Sensitivity specificity etc
  text(x=bp[,34][1],y=79,software,col=nicegrey,adj=0,cex=0.8,xpd=NA,font=2)
  text(x=bp[,34][1],y=65,paste("Percentage of all alleles genotyped: ",round(Percentage_of_all_samples_genotyped,2),"%",sep=""),col=nicegrey,adj=0,cex=0.8,xpd=NA)
  text(x=bp[,34][1],y=58,paste("Percentage of positive alleles genotyped: ",round(Percentage_of_positive_samples_genotyped,2),"%",sep=""),col=nicegrey,adj=0,cex=0.8,xpd=NA)
  text(x=bp[,34][1],y=51,paste("Sensitivity: ",round(Sensitivity,2),"%",sep=""),col=nicegrey,adj=0,cex=0.8,xpd=NA)
  text(x=bp[,34][1],y=44,paste("Specificity: ",round(Specificity,2),"%",sep=""),col=nicegrey,adj=0,cex=0.8,xpd=NA)
  rect(xleft=bp[,34][1]-0.01*bp[,40][1],xright=bp[,40][1]+0.075*bp[,40][1],ytop=70,ybottom=39,col=NA,border="grey",lwd=0.8,xpd=NA)
  i=i+1
}
######
# exSTRa
######
#exSTRa plot 1
# Technically doing two plots on one to avoid redoing margins and page setup 
blank_plot(
  xmin=0
  ,xmax=450
  ,ymin=0
  ,ymax=1
)
# Randomise plotting order to avoid clumps of blue and red 
for (sample in sample(unique(ALS_WGS_PCR_FREE_exSTRa_scores_Targeted$Sample))){
  if (ALS_WGS_PCR_FREE_exSTRa_scores_Targeted$C9orf72_Status[ALS_WGS_PCR_FREE_exSTRa_scores_Targeted$Sample==sample]==1){
    colour_to_plot <- nicered
  } else if (ALS_WGS_PCR_FREE_exSTRa_scores_Targeted$CC[ALS_WGS_PCR_FREE_exSTRa_scores_Targeted$Sample==sample]==1){
    colour_to_plot <- niceorange
  } else {
    colour_to_plot <- nicedarkblue
  }
  sample_ecdf <- ecdf(sort(ALS_WGS_PCR_FREE_exSTRa_scores_Targeted$rep[ALS_WGS_PCR_FREE_exSTRa_scores_Targeted$Sample==sample]))
  # Plot gridlines
  lines(x=c(-17,150),y=c(25,25),lty=2,lwd=0.8,col="grey",xpd=NA)
  lines(x=c(-17,150),y=c(50,50),lty=2,lwd=0.8,col="grey",xpd=NA)
  lines(x=c(-17,150),y=c(75,75),lty=2,lwd=0.8,col="grey",xpd=NA)
  lines(x=c(-17,150),y=c(100,100),lty=2,lwd=0.8,col="grey",xpd=NA)
  plot(sample_ecdf, verticals=TRUE, do.points=FALSE,add=TRUE,col=colour_to_plot,lwd=0.2,col.01line = NULL)
}
# plot white rectangle to cover edges
rect(xleft=155,xright=200,ytop=1.1,ybottom=0.9,col="white",border=NA,lwd=0.8,xpd=NA)
rect(xleft=-25,xright=-5,ytop=0.1,ybottom=-0.1,col="white",border=NA,lwd=0.8,xpd=NA)
lines(x=c(-17,-17),y=c(0,1),lwd=1,col=nicegrey,xpd=NA)
lines(x=c(0,150),y=c(-0.1,-0.1),lwd=1,col=nicegrey,xpd=NA)
# plot axis ticks 
lines(x=c(-20,-17),y=c(0,0),lwd=1,col=nicegrey,xpd=NA)
lines(x=c(-20,-17),y=c(0.25,0.25),lwd=1,col=nicegrey,xpd=NA)
lines(x=c(-20,-17),y=c(0.5,0.5),lwd=1,col=nicegrey,xpd=NA)
lines(x=c(-20,-17),y=c(0.75,0.75),lwd=1,col=nicegrey,xpd=NA)
lines(x=c(-20,-17),y=c(1,1),lwd=1,col=nicegrey,xpd=NA)
lines(x=c(0,0),y=c(-0.1,-0.1-3/150),lwd=1,col=nicegrey,xpd=NA)
lines(x=c(150,150),y=c(-0.1,-0.1-3/150),lwd=1,col=nicegrey,xpd=NA)
# Plot axis labels
text(x=c(-21,-21,-21,-21,-21),y=c(0,0.25,0.5,0.75,1),c(0,0.25,0.5,0.75,1),cex=1,col=nicegrey,adj=1,xpd=NA)
text(x=c(0,150),y=c(-0.1-15/150,-0.1-15/150),c(0,150),cex=1,col=nicegrey,adj=0.5,xpd=NA)
text(x=-45,y=0.5,"ecdf(Repeated Bases)",srt=90,adj=0.5,cex=1,col=nicegrey,xpd=NA)
text(x=75,y=-0.3,"Repeated Bases",,adj=0.5,cex=1,col=nicegrey,xpd=NA)
# plot title
# This is left aligned 
text(x=10,y=1.10,bquote(bold(.("exSTRa") ~ ": " ~ italic(.("C9orf72")) ~ " Allele Predicition")),adj=0,cex=1,col=nicegrey,xpd=TRUE,font=2)

####
# exSTRa plot 2 
###
# Going to compare exSTRa to EHV2 (the most accurate to see how well it classifies the output)
EHV2_C9_ALS <- subset(EHV2[EHV2$Dataset=="ALS_WGS_PCR_FREE" & EHV2$Gene=="C9orf72",],select=c(Sample,A2))
ALS_WGS_PCR_FREE_exSTRA_pvalues <- exSTRa[exSTRa$Dataset=="ALS_WGS_PCR_FREE",]
exSTRA_ALS_C9 <- subset(ALS_WGS_PCR_FREE_exSTRA_pvalues[ALS_WGS_PCR_FREE_exSTRA_pvalues$Gene=="C9orf72" & ALS_WGS_PCR_FREE_exSTRA_pvalues$Sample %in% unique(EHV2_C9_ALS$Sample),],select=c(Sample,P_Value))
exSTRA_ALS_C9 <- merge(exSTRA_ALS_C9,EHV2_C9_ALS,all.x=F,all.y=F)
# Want to plot this beside other exSTRa plots so normalising to be same units as that (but give it correct labels later)
exSTRA_ALS_C9_normalise <- exSTRA_ALS_C9
exSTRA_ALS_C9_normalise$A2 <- 200+exSTRA_ALS_C9_normalise$A2/4
exSTRA_ALS_C9_normalise$P_Value <- (-log10(exSTRA_ALS_C9_normalise$P_Value))/3
exSTRA_ALS_C9_normalise <- merge(exSTRA_ALS_C9_normalise,ALS_CC_STATUS,all.x=T,all.y=F)
exSTRA_ALS_C9_normalise$C9orf72_Status <- ifelse(exSTRA_ALS_C9_normalise$C9orf72_Status==1,see_through_colour(nicered,0.80),see_through_colour(niceorange,0.80))
exSTRA_ALS_C9_normalise$C9orf72_Status <- ifelse(is.na(exSTRA_ALS_C9_normalise$C9orf72_Status)==T,see_through_colour(niceorange,0.80),exSTRA_ALS_C9_normalise$C9orf72_Status)
# plot p-value cutoff
lines(x=c(200,350),y=c((-log10(0.05/34))/3,(-log10(0.05/34))/3),col="grey",lwd=0.7,xpd=NA)
text(x=200,y=1,"p-value threshold",cex=0.8,col=nicegrey,adj=0,xpd=NA)
points(exSTRA_ALS_C9_normalise$A2,exSTRA_ALS_C9_normalise$P_Value,col=NA,bg=exSTRA_ALS_C9_normalise$C9orf72_Status,pch=21,cex=1.2,xpd=NA,lwd=0.3)
# Plot axes
# x axis
lines(x=c(200,350),y=c(-0.1,-0.1),lwd=1,col=nicegrey,xpd=NA)
# y axis
lines(x=c(187,187),y=c(0,1),lwd=1,col=nicegrey,xpd=NA)
# plot axis ticks 
lines(x=c(184,187),y=c(0,0),lwd=1,col=nicegrey,xpd=NA)
lines(x=c(184,187),y=c(0.3333333,0.3333333),lwd=1,col=nicegrey,xpd=NA)
lines(x=c(184,187),y=c(0.6666666,0.6666666),lwd=1,col=nicegrey,xpd=NA)
lines(x=c(184,187),y=c(1,1),lwd=1,col=nicegrey,xpd=NA)
lines(x=c(200,200),y=c(-0.1,-0.1-3/150),lwd=1,col=nicegrey,xpd=NA)
lines(x=c(350,350),y=c(-0.1,-0.1-3/150),lwd=1,col=nicegrey,xpd=NA)
lines(x=c(500,500),y=c(-0.1,-0.1-3/150),lwd=1,col=nicegrey,xpd=NA)
lines(x=c(650,650),y=c(-0.1,-0.1-3/150),lwd=1,col=nicegrey,xpd=NA)
# Plot axis labels
text(x=c(183,183,183,183),y=c(0,0.3333333,0.6666666,1),c(0,1,2,3),cex=1,col=nicegrey,adj=1,xpd=NA)
text(x=c(200,237.5,275,312.5,350),y=c(-0.1-15/150,-0.1-15/150,-0.1-15/150,-0.1-15/150,-0.1-15/150),c(0,150,300,450,600),cex=1,col=nicegrey,adj=0.5,xpd=NA)
# 
text(x=170,y=0.5,expression('-log'[10]*' (exSTRa p-value)'),srt=90,adj=0.5,cex=1,col=nicegrey,xpd=NA)
text(x=264,y=-0.3,"ExpansionHunterv2 Allele Length Prediction (bp)",,adj=0.5,cex=1,col=nicegrey,xpd=NA)
# exSTRa Stats 
exSTRa_sensitivity <- nrow(exSTRA_ALS_C9[exSTRA_ALS_C9$P_Value < 0.05/34,])*100/positive_count
# Have manually checked these stats 
exSTRa_specificity <- 100
exSTRa_Percentage_of_all_samples_genotyped <- 100
exSTRa_Percentage_of_positive_samples_genotyped <- 100
# Sensitivity specificity etc
text(x=375,y=0.79,"exSTRa",col=nicegrey,adj=0,cex=0.8,xpd=NA,font=2)
text(x=375,y=0.65,paste("Percentage of all alleles genotyped: ",round(exSTRa_Percentage_of_all_samples_genotyped,2),"%",sep=""),col=nicegrey,adj=0,cex=0.8,xpd=NA)
text(x=375,y=0.58,paste("Percentage of positive alleles genotyped: ",round(exSTRa_Percentage_of_positive_samples_genotyped,2),"%",sep=""),col=nicegrey,adj=0,cex=0.8,xpd=NA)
text(x=375,y=0.51,paste("Sensitivity: ",round(exSTRa_sensitivity,2),"%",sep=""),col=nicegrey,adj=0,cex=0.8,xpd=NA)
text(x=375,y=0.44,paste("Specificity: ",round(exSTRa_specificity,2),"%",sep=""),col=nicegrey,adj=0,cex=0.8,xpd=NA)
rect(xleft=370,xright=473,ytop=0.70,ybottom=0.39,col=NA,border="grey",lwd=0.8,xpd=NA)
# Plot Letter 
if (i <= 26){
  plot_letter=LETTERS[i]
} else {
  j=i
  j=j-26
  plot_letter=paste(LETTERS[j],".2",sep="")
}
plot_letter_label(0,450*0.9,0,1,plot_letter)
i=i+1

####
# Expansion Hunter de Novo Plot 
###
# create plot 
# Manually set up plot zone
layout(matrix(c(1,2,3,4,5,6),nrow=3,byrow=T),heights=c(1,1,1,1,1,1),widths=c(3,1,3,1,3,1))
par(lwd=2,fg="grey30",mar=c(6,9,2,1))
i=i+1
manhattan_and_qq(EHdN_CC,expression("ExpansionHunter de novo: Manhattan plot (272 cases, 136 Controls)"))
manhattan_and_qq(EHdN_C9,expression(paste("ExpansionHunter de novo:Manhattan plot (26 ",italic("C9orf72")," Positive Cases, 136 Controls)",sep="")))
# Rather than be awkward and try to put plot letters in our manhattan plots going to make a blank plot below and place plot letters 
blank_plot(0,100,0,100)
# Plot Letter 
if (i <= 26){
  plot_letter=LETTERS[i]
} else {
  j=i
  j=j-26
  plot_letter=paste(LETTERS[j],".2",sep="")
}
plot_letter_label(0,100*0.9,0,480,plot_letter)
i=i+1
# Plot Letter 
if (i <= 26){
  plot_letter=LETTERS[i]
} else {
  j=i
  j=j-26
  plot_letter=paste(LETTERS[j],".2",sep="")
}
plot_letter_label(0,100*0.9,0,300,plot_letter)
dev.off()
print(" Plot Chapter3_C9orf72_Benchmarking.pdf Complete")

################################################################################
################################################################################

# Figure 3.7 : C9orf72 in silico Genotyping, Specificity as a Function of Sensitivity

################################################################################
################################################################################

#########################################
# PLOT OF SENSITIVITY AND SPECIFICITY GRADIENT ACROSS C9orf72 
#########################################
sens_spec_calc_c9 <- function(Working_Dataset,primary_dataset){
  # A function to calculate the specificity of a tool as sensitivity ranges from 0-100
  # Create empty variables
  Working_Sensitivity_C9 <- c()
  Working_Specificity_C9 <- c()
  # Find list of samples with true positive expansion
  # leave this as EHV3 it will be right for all
  Positive_Samples<- pmid_vs_cc_status$Sample_well[pmid_vs_cc_status$c9orf72_status_clean==1 & !is.na(pmid_vs_cc_status$c9orf72_status_clean) & (pmid_vs_cc_status$Sample_well %in% unique(EHV3$Sample[EHV3$Dataset=="ALS_WGS_PCR_FREE"]))]
  # Find list of samples without expansion
  Negative_Samples <- pmid_vs_cc_status$Sample_well[!(pmid_vs_cc_status$Sample_well %in% Positive_Samples)]
  positive_count<-length(Positive_Samples) 
  negative_count<-total_sample_count - positive_count
  # Number of positive samples that are genotyped (either positively or negatively)
  positives_genotyped_count <- nrow(Working_Dataset[Working_Dataset$Sample %in% Positive_Samples & !is.na(Working_Dataset$A2) & Working_Dataset$CC_C9=="C9_Positive_Patient",])
  # Number of negative samples that are genotyped (either positively or negatively)
  negatives_genotyped_count <- nrow(Working_Dataset[Working_Dataset$Sample %in% Negative_Samples & !is.na(Working_Dataset$A2),])
  positive_allele_calls <- sort(Working_Dataset$A2[Working_Dataset$CC_C9=="C9_Positive_Patient"])
  # Cycle through each value of predicted length in true positive samples, treating each as a cutoff
  for (min_call in unique(positive_allele_calls)){
    # Number of C9 negatives which are above or equal to the current minimum
    false_positives <- length(Working_Dataset$Sample[Working_Dataset$CC_C9!="C9_Positive_Patient" & Working_Dataset$A2 >= min_call])
    # Number of C9 positives which are below the current minimum 
    false_negatives <- length(positive_allele_calls[positive_allele_calls<min_call & !is.na(positive_allele_calls)])
    # Calculate Sens and Spec
    Sensitivity <- (positives_genotyped_count-false_negatives)*100/positives_genotyped_count
    Specificity <- (negatives_genotyped_count-false_positives)*100/negatives_genotyped_count
    # Add to previous results
    Working_Sensitivity_C9 <- c(Sensitivity,Working_Sensitivity_C9)
    Working_Specificity_C9 <- c(Specificity,Working_Specificity_C9)
  }
  # Add a 0,0 for plotting
  Working_Sensitivity_C9 <- c(0,Working_Sensitivity_C9)
  Working_Specificity_C9 <- c(0,Working_Specificity_C9)
  # Calculate overall point values for sens and spec
  Sensitivity_Point <- (positives_genotyped_count-length(positive_allele_calls[positive_allele_calls<30]))*100/positives_genotyped_count
  Specificity_Point <- (negatives_genotyped_count-length(Working_Dataset$Sample[Working_Dataset$CC_C9!="C9_Positive_Patient" & Working_Dataset$A2 >= 30]))*100/negatives_genotyped_count
  if (primary_dataset=="exSTRa"){
    Sensitivity_Point <- (positives_genotyped_count-nrow(Working_Dataset[Working_Dataset$P_Value>0.05/34 & Working_Dataset$CC_C9=="C9_Positive_Patient",]))*100/positives_genotyped_count
    Specificity_Point <- (negatives_genotyped_count-nrow(Working_Dataset[Working_Dataset$P_Value<0.05/34 & Working_Dataset$CC_C9!="C9_Positive_Patient",]))*100/negatives_genotyped_count
  }
  to_return=list(Sensitivity=Working_Sensitivity_C9,Specificity=Working_Specificity_C9,Sensitivity_Point=Sensitivity_Point,Specificity_Point=Specificity_Point)
  return(to_return)
}

# Cycle through each dataset
for (primary_dataset in vector_of_datasets2){
  # Call Dataset 
  Working_Dataset <-get(primary_dataset)[get(primary_dataset)$Gene=="C9orf72" & get(primary_dataset)$Dataset=="ALS_WGS_PCR_FREE",]
  # Merge 
  Working_Dataset <- merge(Working_Dataset,ALS_CC_STATUS,all.x=T,all.y=F)
  Working_Dataset$C9orf72_Status <- ifelse(is.na(Working_Dataset$C9orf72_Status)==T,0,Working_Dataset$C9orf72_Status)
  Working_Dataset$CC_C9 <- ifelse(Working_Dataset$C9orf72_Status==1,"C9_Positive_Patient",ifelse(Working_Dataset$CC==1,"C9_Negative_Patient","Control"))

  if(primary_dataset != "STRetch"){
    # total 408 cases and controls 
    total_sample_count=408
  } else {
    # STRetch only gives output for cases 
    total_sample_count=272
  }
  # Doing exSTRa based on p-value so easiest to create an artificial A2 column
  if (primary_dataset=="exSTRa"){
    Working_Dataset$A2 <- 1-Working_Dataset$P_Value
  }
  assign(paste(primary_dataset,"_Sensitivity_Point",sep=""),sens_spec_calc_c9(Working_Dataset,primary_dataset)$Sensitivity_Point)
  assign(paste(primary_dataset,"_Specificity_Point",sep=""),sens_spec_calc_c9(Working_Dataset,primary_dataset)$Specificity_Point)
  assign(paste(primary_dataset,"_Sensitivity",sep=""),sens_spec_calc_c9(Working_Dataset,primary_dataset)$Sensitivity)
  assign(paste(primary_dataset,"_Specificity",sep=""),sens_spec_calc_c9(Working_Dataset,primary_dataset)$Specificity)
}

# PDF = 210x297 - margins are 35,20,20,20 mm 
# PDF = 8.3x11.7 - margins are 1.38,0.79,0.79,0.79 inches
pdf(paste(output_directory,"Chapter3_C9orf72_Sens_Spec_Gradient.pdf",sep=""),width=(8.3-1.38-0.8),height=(11.7-0.8-0.8))
par(mfrow=c(4,2),mar=c(2.5,4.5,3.5,0.5))
###
# Create line plot of sens / spec
###
blank_plot(0,100,0,100)
for (primary_dataset in vector_of_datasets2){
  lines(x=get(paste(primary_dataset,"_Sensitivity",sep="")),y=get(paste(primary_dataset,"_Specificity",sep="")),col=vector_of_software_colours[match(primary_dataset,vector_of_datasets2)],lwd=1.2)
}
axis(1,at=c(0,25,50,75,100),lab=c(0,25,50,75,100),las=1,lwd=1,col=nicegrey,col.axis=nicegrey,cex.axis=1)
axis(2,at=c(0,25,50,75,100),lab=c(0,25,50,75,100),las=1,lwd=1,col=nicegrey,col.axis=nicegrey,cex.axis=1)
# Title for y-axis 1 
title(xlab="Sensitivity (%)",ylab="Specificity (%)", cex.lab=1,col.lab=nicegrey,line=2.5,font=1.1,xpd=NA)
title(main=expression(paste("Software Accuracy at ",italic("C9orf72")," Locus",sep="")), cex.lab=1,col.lab=nicegrey,line=1.5,font=1.1,xpd=NA)
###
# Legend Plot 
###
blank_plot(0,100,0,100)
height=90
for (primary_dataset in vector_of_datasets2){
  points(x=0,y=height,cex=1.5,pch=21,col="white",bg=vector_of_software_colours[match(primary_dataset,vector_of_datasets2)],xpd=NA)
  text(x=5,y=height,cex=1.3,vector_of_software_names[match(primary_dataset,vector_of_datasets2)],col=nicegrey,adj=0,xpd=NA)
  height=height-10
}
dev.off()
print("Plot of Chapter3_C9orf72_Sens_Spec_Gradient.pdf complete")

################################################################################
################################################################################

# Figure 3.8, Supplementary FigureSXXX- : Comparison of Gold Standard PCR Genotyping with in silico Predictions

################################################################################
################################################################################

#########################################
# Compare to gold standard PCR genotypes
#########################################

for (primary_dataset in vector_of_datasets){
  # This software name is used for the plot title and the name the output PDF 
  software=vector_of_software_names[match(primary_dataset,vector_of_datasets)]
  # Filter to PCR FREE genotypes
  Working_Dataset <- get(primary_dataset)[get(primary_dataset)$Dataset=="ALS_WGS_PCR_FREE",]
  Working_Dataset <- merge(Working_Dataset,PCR_Genotypes_Working,all.x=T,all.y=F)
  Working_Dataset <- Working_Dataset[!is.na(Working_Dataset$PCR_A1) & !is.na(Working_Dataset$PCR_A2),]
  # Currently a problem with the SCA1 & DIP2B PCR genotypes and C9orf72 has been looked at elsewhere
  Working_Dataset <- Working_Dataset[Working_Dataset$Gene!="ATXN1" & Working_Dataset$Gene!="C9orf72" & Working_Dataset$Gene!="DIP2B",]

  # Different softwares simply start counting the alleles at certain loci from different a different point to PCR, using a correction to bring them in concurrence with the PCR 
  #
  EHV2_PCR_Gene_Correction_Gene <-c(  "AR", "ATN1", "ATXN3","ATXN7",  "GIPC1","HTT","NOP56","PPP2R2B","TBP")
  EHV2_PCR_Gene_Correction_Value <-c(  -1,    -4,     3,     -5,       -2,      -3,     2,    -2,        -1)
  EHV3_PCR_Gene_Correction_Gene <-c(  "AR", "ATN1", "ATXN3","ATXN7",  "GIPC1","HTT","NOP56","PPP2R2B","TBP")
  EHV3_PCR_Gene_Correction_Value <-c(  -1,    -4,     3,     -5,             -2,      -3,     2,    -2,        -1)
  GangSTR_NonTargeted_PCR_Gene_Correction_Gene <-c(  "AR", "ATN1", "ATXN3","ATXN7",  "GIPC1","HTT","NOP56","PPP2R2B","TBP")
  GangSTR_NonTargeted_PCR_Gene_Correction_Value <-c(  0,    0,     6,     -5,         -2,      -3,     4,    -2,        -1)
  GangSTR_Targeted_PCR_Gene_Correction_Gene <-c(  "AR", "ATN1", "ATXN3","ATXN7",  "GIPC1","HTT","NOP56","PPP2R2B","TBP")
  GangSTR_Targeted_PCR_Gene_Correction_Value <-c(  0,    0,     6,     -5,            -2,      -3,     4,    -2,        -1)
  hipstr_PCR_Gene_Correction_Gene <-c(  "AR", "ATN1", "ATXN3","ATXN7",  "GIPC1","HTT","NOP56","PPP2R2B","TBP","JPH3")
  hipstr_PCR_Gene_Correction_Value <-c(  3,    -4,     3,     -5,            -6,  0,     1,    -3,        -1,-2)
  tredparse_PCR_Gene_Correction_Gene <- c("AR", "ATN1", "ATXN3","ATXN7",  "GIPC1","HTT","NOP56","PPP2R2B","TBP")
  tredparse_PCR_Gene_Correction_Value <- c(  1,    0,     6,     -5,       -2,      -3,     1,    -2, -2     )
  STRetch_PCR_Gene_Correction_Gene <-c(  "AR", "ATN1", "ATXN3","ATXN7",  "GIPC1","HTT","NOP56","PPP2R2B","TBP")
  STRetch_PCR_Gene_Correction_Value <-c(  -1,    -4,     3,     -5,             -2,      -3,     2,    -2,        -1)
  repeatseq_PCR_Gene_Correction_Gene <-c(  "AR", "ATN1", "ATXN3","ATXN7", "GIPC1","HTT","NOP56","PPP2R2B","TBP")
  repeatseq_PCR_Gene_Correction_Value <-c(  0,    -4,     3,     0,      -2,      -7,     2,    -2,        -1)
  # Compare WGS to WES calls 
  # PDF = 210x297 - margins are 35,20,20,20 mm 
  # PDF = 8.3x11.7 - margins are 1.38,0.79,0.79,0.79 inches
  pdf(paste(output_directory,"Chapter3_WGS_PCR_Comparison",primary_dataset,".pdf",sep=""),width=(11.7-0.8-0.8),height=(8.3-1.38-0.8))
  par(mfrow=c(3,4),mar=c(3.5,3.5,4.5,2))
  i=1
  for (gene in sort(unique(Working_Dataset$Gene))){
    Working_Dataset_1<-Working_Dataset[Working_Dataset$Gene==gene,]
    # randomly shuffle to avoid colour clumps 
    set.seed(42)
    rows <- sample(nrow(Working_Dataset_1))
    Working_Dataset_1 <- Working_Dataset_1[rows,]
    if (gene %in% get(paste(primary_dataset,"_PCR_Gene_Correction_Gene",sep=""))){
      Working_Dataset_1$PCR_A1 <- as.numeric(Working_Dataset_1$PCR_A1)-as.numeric(get(paste(primary_dataset,"_PCR_Gene_Correction_Value",sep=""))[match(gene,get(paste(primary_dataset,"_PCR_Gene_Correction_Gene",sep="")))])
      Working_Dataset_1$PCR_A2 <- as.numeric(Working_Dataset_1$PCR_A2)-as.numeric(get(paste(primary_dataset,"_PCR_Gene_Correction_Value",sep=""))[match(gene,get(paste(primary_dataset,"_PCR_Gene_Correction_Gene",sep="")))])
    }
    if (software=="STRetch"){
      EH_alleles <-c(Working_Dataset_1$A2)
      PCR_alleles <- c(Working_Dataset_1$PCR_A2)
    } else {
      EH_alleles <-c(Working_Dataset_1$A1,Working_Dataset_1$A2)
      PCR_alleles <- c(Working_Dataset_1$PCR_A1,Working_Dataset_1$PCR_A2)
    }
    difference <- c(PCR_alleles-EH_alleles)
    rmsd <- (sum(difference^2)/length(difference))^(1/2)
    point_colour <- c(Working_Dataset_1$CC,Working_Dataset_1$CC)
    point_colour <- ifelse(point_colour==1,see_through_colour(hexcode=nicered,fade_proportion=0.4),see_through_colour(hexcode=nicedarkblue,fade_proportion=0.4))
    plot_min=min(c(EH_alleles,PCR_alleles))
    plot_max=max(c(EH_alleles,PCR_alleles))
    blank_plot(
      xmin=plot_min
      ,xmax=plot_max
      ,ymin=plot_min
      ,ymax=plot_max
    )
    # Plot 45 degree line across plot 
    lines(x=c(-20,1.5*plot_max),y=c(-20,1.5*plot_max),col=nicegrey,lwd=1,lty=1)
    # Add individual plot titles 
    title(main=bquote(italic(.(gene))), cex.main=1.1,cex.lab=1,col.lab=nicegrey,xpd=NA,line=0.7,font=3)
    # Add axis labels 
    title(ylab="Software Allele Prediction", xlab="PCR Allele Measurement", cex.main=1.1,cex.lab=1,col.lab=nicegrey,xpd=NA,line=2.2)
    # Draw external box 
    box(lwd=2,col=nicegrey)
    # plot x axis ticks 
    axis(side=1,lwd=1,cex.axis=1,col.axis=nicegrey,col=nicegrey)
    # plot y axis ticks 
    axis(side=2,lwd=1,las=1,cex.axis=1,col.axis=nicegrey,col=nicegrey)
    # plot points giving each a unique colour to correlate with other plots 
    points(x=PCR_alleles,y=EH_alleles,col=point_colour,bg=point_colour,pch=21,cex=1.2,xpd=NA,lwd=0.3)
    # Plot RMSD 
    text(x=plot_max-0.35*(plot_max-plot_min),y=plot_min+0.1*(plot_max-plot_min),paste("rmsd= ",round(rmsd,2),sep=""),cex=1,adj=0,col=nicegrey)
    # Plot Legend on first plot 
    if (i==1){
      # Plot title above first plot 
      text(x=plot_min,y=(plot_max+0.4*(plot_max-plot_min)),paste(software,": Comparison of Gold Standard PCR Genotyping with Software Allele Prediction"),cex=1.5,adj=0,col=nicegrey,xpd=NA,font=2)

      # Plot Coloured Points
      points(x=c(plot_min,plot_min),y=c(plot_max-0.05*(plot_max-plot_min),plot_max-0.15*(plot_max-plot_min)),col=c(nicered,nicedarkblue),bg=c(nicered,nicedarkblue),pch=21,cex=1.4,lwd=0.6)
      # Plot Case/Control Labels
      text(x=c(plot_min+0.1*(plot_max-plot_min),plot_min+0.1*(plot_max-plot_min)),y=c(plot_max-0.05*(plot_max-plot_min),plot_max-0.15*(plot_max-plot_min)),c("Patient","Control"),col=nicegrey,adj=0,cex=1)
    }
    i=1+1
  }
  dev.off()
}
print("Plot Chapter3_WGS_PCR_Comparison Complete")

################################################################################
################################################################################

# Figure 3.9 : Comparison of Genotype Calls from Samples Sequenced with WES and WGS

################################################################################
################################################################################

for (primary_dataset in vector_of_datasets){
  if (primary_dataset != "STRetch"){
    total_vector_WES<-c()
    total_vector_WGS<-c()
    # This software name is used for the plot title and the name the output PDF 
    software=vector_of_software_names[match(primary_dataset,vector_of_datasets)]
    # Compare WGS to WES calls 
    # PDF = 210x297 - margins are 35,20,20,20 mm 
    # PDF = 8.3x11.7 - margins are 1.38,0.79,0.79,0.79 inches
    pdf(paste(output_directory,"Chapter3_WES_WGS_rmsd_",software,".pdf",sep=""),width=(11.7-0.8-0.8),height=(8.3-1.38-0.8))
    par(mfrow=c(3,4),mar=c(3.5,3.5,4.5,2))
    i=1
    # Cycle through each gene 
    for (gene in sort(unique(get(paste(primary_dataset,"_doubles",sep=""))$Gene[get(paste(primary_dataset,"_doubles",sep=""))$Dataset=="RCSI_WES"]))) {
      # Create working dataset 
      Working_Dataset=get(paste(primary_dataset,"_doubles",sep=""))
      # Filter to gene 
      Working_Dataset <- Working_Dataset[Working_Dataset$Gene==gene,]
      # Find samples that appear twice (not all samples will have a call in both exome and genome for many reasons)
      retain_samples <- names(table(Working_Dataset$Sample)[table(Working_Dataset$Sample)>1])
      if (length(retain_samples)>=1){
        # Filter to samples with both WES and WGS result 
        Working_Dataset <- Working_Dataset[Working_Dataset$Sample %in% retain_samples,]
        # Fill Lists to get overall RMSD
        total_vector_WES<-c(total_vector_WES,as.numeric(Working_Dataset$A1[Working_Dataset$Dataset=="RCSI_WES"]),as.numeric(Working_Dataset$A2[Working_Dataset$Dataset=="RCSI_WES"]))
        total_vector_WGS<-c(total_vector_WGS,as.numeric(Working_Dataset$A1[Working_Dataset$Dataset=="RCSI_WGS_PCR_FREE"]),as.numeric(Working_Dataset$A2[Working_Dataset$Dataset=="RCSI_WGS_PCR_FREE"]))
        # Convert dataframe format 
        Working_Dataset_1 <- data.frame(Sample=c(paste(retain_samples,"_A1",sep=""),paste(retain_samples,"_A2",sep="")),WES=NA,WGS=NA)
        for (sample in retain_samples){
          Working_Dataset_1$WES[Working_Dataset_1$Sample==paste(sample,"_A1",sep="")] <- Working_Dataset$A1[Working_Dataset$Sample==sample & Working_Dataset$Dataset=="RCSI_WES"]
          Working_Dataset_1$WES[Working_Dataset_1$Sample==paste(sample,"_A2",sep="")] <- Working_Dataset$A2[Working_Dataset$Sample==sample & Working_Dataset$Dataset=="RCSI_WES"]
          Working_Dataset_1$WGS[Working_Dataset_1$Sample==paste(sample,"_A1",sep="")] <- Working_Dataset$A1[Working_Dataset$Sample==sample & Working_Dataset$Dataset=="RCSI_WGS_PCR_FREE"]
          Working_Dataset_1$WGS[Working_Dataset_1$Sample==paste(sample,"_A2",sep="")] <- Working_Dataset$A2[Working_Dataset$Sample==sample & Working_Dataset$Dataset=="RCSI_WGS_PCR_FREE"]
        }
        ###
        # Calculate rmsd
        ###
        Working_Dataset_1$WES <- as.numeric(Working_Dataset_1$WES)
        Working_Dataset_1$WGS <- as.numeric(Working_Dataset_1$WGS)
        Working_Dataset_1$diff <- ((Working_Dataset_1$WES-Working_Dataset_1$WGS)^2)^0.5
        rmsd <- (sum(Working_Dataset_1$diff^2)/nrow(Working_Dataset_1))^(1/2)
        # Some gene only have a single allele and this plots funny so adding one point either side of this 
        unique_alleles <- c(sort(as.numeric(unique(c(Working_Dataset_1$WGS,Working_Dataset_1$WES)))))
        # Get locations of axis ticks 
        ifelse(length(unique_alleles) >= 2,axis_tick_locations <- unique_alleles ,axis_tick_locations <- c(unique_alleles -1,unique_alleles,unique_alleles +1))
        # Create blank plot with the right dimensions
        blank_plot(
          xmin=min(axis_tick_locations)
          ,xmax=max(axis_tick_locations)
          ,ymin=min(axis_tick_locations)
          ,ymax=max(axis_tick_locations)
          )
        # Plot 45 degree line across plot 
        lines(x=c(0.5*min(axis_tick_locations),1.5*max(axis_tick_locations)),y=c(0.5*min(axis_tick_locations),1.5*max(axis_tick_locations)),col=nicegrey,lwd=1,lty=1)
        # Add individual plot titles 
        title(main=bquote(italic(.(gene))), cex.main=1.1,cex.lab=1,col.lab=nicegrey,xpd=NA,line=0.7,font=3)
        # Add axis labels 
        title(ylab="WES Allele Call", xlab="WGS Allele Call", cex.main=1.1,cex.lab=1,col.lab=nicegrey,xpd=NA,line=2.2)
        # Draw external box 
        box(lwd=2,col=nicegrey)
        # plot x axis ticks 
        axis(side=1,lwd=1,cex.axis=1,col.axis=nicegrey,col=nicegrey,at=axis_tick_locations)
        # plot y axis ticks 
        axis(side=2,lwd=1,las=1,cex.axis=1,col.axis=nicegrey,col=nicegrey,at=axis_tick_locations)
        # plot points giving each a unique colour to correlate with other plots 
        points(x=Working_Dataset_1$WGS,y=Working_Dataset_1$WES,col=see_through_colour(hexcode=colorRampPalette(vangogh_palette)(length(gene_list_quotes))[match(gene,gene_list_quotes)],fade_proportion=0.9),bg=see_through_colour(hexcode=colorRampPalette(vangogh_palette)(length(gene_list_quotes))[match(gene,gene_list_quotes)],fade_proportion=0.2),pch=21,cex=1.2,xpd=NA,lwd=0.3)
        # Want red font if rmsd > 1 
        if (rmsd <1){
          text(x=min(axis_tick_locations),y=max(axis_tick_locations)-0.1*(max(axis_tick_locations)-min(axis_tick_locations)),paste("rmsd= ",round(rmsd,2),sep=""),cex=1,adj=0,col=nicegrey)
        } else {
          text(x=min(axis_tick_locations),y=max(axis_tick_locations)-0.1*(max(axis_tick_locations)-min(axis_tick_locations)),paste("rmsd= ",round(rmsd,2),sep=""),cex=1,adj=0,col=nicered)
        }
        # Plot title  above first plot 
        if(i==1){
          text(x=(min(axis_tick_locations)),y=(max(axis_tick_locations)+0.4*(max(axis_tick_locations)-min(axis_tick_locations))),paste(software,": Comparison of WGS and WES Allele Calls in the Same Samples",sep=""),cex=1.5,adj=0,col=nicegrey,xpd=NA)
        }
        # Plot letter 
        plot_letter_label_b(x_min=min(axis_tick_locations),x_max=max(axis_tick_locations),y_min=min(axis_tick_locations),y_max=max(axis_tick_locations),letter=LETTERS[i])
        i=i+1
      }
    }
    dev.off()
    # Calculate overall RMSD 
    overall_diff <- ((total_vector_WES-total_vector_WGS)^2)^0.5
    overall_rmsd <- (sum(overall_diff^2)/length(overall_diff))^(1/2)
    print(paste(software,": Overall WES WGS RMSD : ",overall_rmsd,sep=""))
    print(paste(software,": Provides a total of ", length(total_vector_WES)/2," calls. The max WES call is ",max(total_vector_WES),". The max WGS call is ", max(total_vector_WGS),". The mean WES call is ",round(mean(total_vector_WES),2)," (SD: ",round(sd(total_vector_WES),2),"). The mean WGS call is ",round(mean(total_vector_WGS),2)," SD: ",round(sd(total_vector_WGS),2),").",sep=""))
  }
}

print("Plot Chapter3_WES_WGS_rmsd_ Complete")


################################################################################
################################################################################

# Figure 3.10, supplementary figure SXXX : Prediction of STR Lengths in Epilepsy Patients

################################################################################
################################################################################

for (primary_dataset in vector_of_datasets){
  # This software name is used for the plot title and the name the output PDF 
  software=vector_of_software_names[match(primary_dataset,vector_of_datasets)]
  Working_Dataset <- get(primary_dataset)
  # Filter to RCSI and ALS Controls and remove parents 
  Working_Dataset <- Working_Dataset[(is.na(Working_Dataset$CC)==T | Working_Dataset$CC==0)&Working_Dataset$Dataset!="ALS_WES"&Working_Dataset$Parent!="Yes",]
  # Filter to genes that are called in epilepsy data 
  epilepsy_genes <- unique(Working_Dataset$Gene[Working_Dataset$Dataset_Master=="RCSI"])
  Working_Dataset <- Working_Dataset[Working_Dataset$Gene %in% epilepsy_genes,]
  # For TREDPARSE don't include WGS PCR 
  if (primary_dataset=="tredparse"){
    Working_Dataset<-Working_Dataset[Working_Dataset$Dataset !="RCSI_WGS_PCR",]
  }
  # Want all plots to have the max number of columns and to align with plot above so need to find the max numbero of columns 
  max_columns <- c()
  for (gene in sort(unique(Working_Dataset$Gene))){
    max_columns <- c(max_columns,length(unique(Working_Dataset$A2[Working_Dataset$Gene==gene])))
  }
  max_columns <- max(max_columns)
  # Add case control status 
  Working_Dataset$Status<- gsub("_.*","",Working_Dataset$Dataset)
  Working_Dataset$Status<- gsub("RCSI","Patient",Working_Dataset$Status)
  Working_Dataset$Status<- gsub("ALS","Control",Working_Dataset$Status)
  # For the barplots we want to precreate the barplots to get the positions we need 
  # Create DF for each gene 
  for (gene in sort(unique(Working_Dataset$Gene))){
    assign(paste(gene,"_working_dataframe",sep=""),create_gene_table(temp_df=Working_Dataset,gene=gene))
  }
  # Precreate each barplot for first type of plots - need this to get x positions 
  for (gene in sort(unique(Working_Dataset$Gene))){
    assign(paste(gene,"_bp1",sep=""),create_bp(gene))
  }
  # Assign new dataframes including blank spaces 
  for (gene in sort(unique(Working_Dataset$Gene))){
    assign(paste(gene,"_working_dataframe",sep=""),bring_dfs_to_max_columns(gene=gene,max_columns=max_columns))
  }

  ##########################################
  # Plot 
  ##########################################

  e<-2.718281828459045

  # Cycle through each gene
  # PDF = 210x297 - margins are 35,20,20,20 mm 
  # PDF = 8.3x11.7 - margins are 1.38,0.79,0.79,0.79 inches
  pdf(paste(output_directory,"Chapter3_REs_BP_OR_",software,".pdf",sep=""),width=(11.7-0.8-0.8),height=(8.3-1.38-0.8))
  par(mfrow=c(4,1))
  i=1
  # Cycle through each gene 
  for (gene in sort(unique(Working_Dataset$Gene))) {
    par(mar=c(2.3,8,3.5,1))
    # Filter to gene of interest
    Working_Dataset_1 <- Working_Dataset[Working_Dataset$Gene==gene,]

    #############
    # First plot - barplot 
    ############

    # Fetch the repeat motif 
    repeat_motif=repeat_motifs[match(gene,gene_list_quotes)]
    # Plot barplot
    # Remove border on bars so that there is no zero line at the bottom of the plot 
    par(lty = 0)
    # For most genes setting the ylimit at 75 is visually best but some require higher - have found these through trial and error 
    ymax=100
    # Calculate barplot 
    bp <- barplot(get(paste(gene,"_working_dataframe",sep="")),
      beside=T,
      ylim=c(0,ymax),         
      col=c(nicedarkblue,nicered),
      #lwd=2,
      yaxt="n",
      xaxt="n",
      space=c(0,0.35),
      cex.names=0.6,
      col.axis=nicegrey,
      xlab="",
      cex.lab=0.1,
      col.lab=nicegrey,
    )
    # Title for y-axis 1 
    title(ylab=bquote(italic(.(gene)) ~ " allele carrier frequency,"),cex.lab=0.9,col.lab=nicegrey,line=4.5,font=3,xpd=NA)
    # Title for y-axis 2
    title(ylab="longer allele (%)",cex.lab=0.9,col.lab=nicegrey,line=3.5,font=3,xpd=NA)
    # Plot Main Title
    text(x=min(bp),y=ymax*1.1,bquote(bold(italic(.(gene)) ~ " " ~ .(repeat_motif) ~ " repeats   (" ~ .(software) ~ ")")),adj=0,cex=1,col=nicegrey,xpd=TRUE,font=2)
    # Grid Lines:  only want to plot these over the original distance of the plot, not the fake distance that keeps all the bars uniform 
    lines(x=c(-5,max(get(paste(gene,"_bp1",sep="")))),y=c(25,25),lty=2,lwd=0.8,col="grey")
    lines(x=c(-5,max(get(paste(gene,"_bp1",sep="")))),y=c(50,50),lty=2,lwd=0.8,col="grey")
    lines(x=c(-5,max(get(paste(gene,"_bp1",sep="")))),y=c(75,75),lty=2,lwd=0.8,col="grey")
    lines(x=c(-5,max(get(paste(gene,"_bp1",sep="")))),y=c(100,100),lty=2,lwd=0.8,col="grey")
    # Replot bars over grid lines
    barplot(get(paste(gene,"_working_dataframe",sep="")),
      beside=T,
      ylim=c(0,ymax),         
      col=c(nicedarkblue,nicered),
      #lwd=2,
      yaxt="n",
      space=c(0,0.35),
      cex.names=0.6,
      col.axis=nicegrey,
      xlab="",
      cex.lab=0.1,
      col.lab=nicegrey,
      #border=c(nicedarkblue,nicered),
      #border=c("white","white"),
      add=T
    )
    # white line at 0 to hide 0 values (that otherwise plot as tiny thin lines)
    abline(h=0,lty=1,lwd=0.1,col="white")
    # Plot axis
    axis(2,at=c(0,25,50,75,100),lab=c(0,25,50,75,100),las=1,lwd=1,col=nicegrey,col.axis=nicegrey,cex.axis=0.9)
    # Plot "repeats:"
    axis(1,at=bp[1],lab="repeats: ",las=1,col=nicegrey,col.axis=nicegrey,cex.axis=0.9,hadj=1,line=NA,tick=F)
    # Plot legend on the first plot on each page 
    # Want to plot this in a different place for the plot with the most columns 
    if (length(unique(Working_Dataset$A2[Working_Dataset$Gene==gene]))>2){
      # Plot Coloured Points
      points(x=c(max(get(paste(gene,"_bp1",sep="")))*0.9,max(get(paste(gene,"_bp1",sep="")))*0.9),y=c(83,67),col=c(nicered,nicedarkblue),bg=c(nicered,nicedarkblue),pch=21,cex=1.4,lwd=0.6)
      # Plot Case/Control Labels
      text(x=c(max(get(paste(gene,"_bp1",sep="")))*0.9+max(bp)*0.01,max(get(paste(gene,"_bp1",sep="")))*0.9+max(bp)*0.01),y=c(83,67),c(paste("Patients (n= ",length(unique(Working_Dataset_1$Sample[Working_Dataset_1$Status=="Patient"])),")",sep=""),paste("Controls (n= ",length(unique(Working_Dataset_1$Sample[Working_Dataset_1$Status=="Control"])),")",sep="")),col=nicegrey,adj=0,cex=0.9)
    } else {
      # Plot Coloured Points
      points(x=c(max(get(paste(gene,"_bp1",sep="")))*2,max(get(paste(gene,"_bp1",sep="")))*2),y=c(83,67),col=c(nicered,nicedarkblue),bg=c(nicered,nicedarkblue),pch=21,cex=1.4,lwd=0.6)
      # Plot Case/Control Labels
      text(x=c(max(get(paste(gene,"_bp1",sep="")))*2+max(bp)*0.01,max(get(paste(gene,"_bp1",sep="")))*2+max(bp)*0.01),y=c(83,67),c(paste("Patients (n= ",length(unique(Working_Dataset_1$Sample[Working_Dataset_1$Status=="Patient"])),")",sep=""),paste("Controls (n= ",length(unique(Working_Dataset_1$Sample[Working_Dataset_1$Status=="Control"])),")",sep="")),col=nicegrey,adj=0,cex=0.9)
    }
    # Plot Letter 
    if (i <= 26){
      plot_letter=LETTERS[i]
    } else {
      j=i
      j=j-26
      plot_letter=paste(LETTERS[j],".2",sep="")
    }
    plot_letter_label(min(bp),max(bp),0,ymax,plot_letter)


    #############
    # Second plot - odds ratio  
    ############
    par(mar=c(4.5,8,0.8,1))
    # Create Dataframe to store odds ratio values 
    Odds_Ratio_DF <- data.frame(rep_count=c(sort(as.numeric(names(table(Working_Dataset_1$A2))))),OR=NA,ORuci=NA,ORlci=NA,fisherman=NA)
    for (rep_count in Odds_Ratio_DF$rep_count){
      ort<-table(Working_Dataset_1$A2>=as.numeric(rep_count),Working_Dataset_1$Status)
      if(nrow(ort)==1) { 
        ort<-rbind(ort,c(0,0)) 
      } # FALSE counts don't make it in if they're all 0 
      ort<-ort+(corr<-0.5) # Haldane-Anscombe correction
      Odds_Ratio_DF$OR[Odds_Ratio_DF$rep_count==rep_count]<-(ort[1]*ort[4])/(ort[2]*ort[3])
      #Odds_Ratio_DF$OR[Odds_Ratio_DF$rep_count==rep_count]<-(ort[2]*ort[3])/(ort[4]*ort[1])
      Odds_Ratio_DF$ORuci[Odds_Ratio_DF$rep_count==rep_count]<-e^(log(Odds_Ratio_DF$OR[Odds_Ratio_DF$rep_count==rep_count])+1.96*sqrt(1/ort[2]+1/ort[3]+1/ort[4]+1/ort[1]))
      #Odds_Ratio_DF$ORuci[Odds_Ratio_DF$rep_count==rep_count]<-e^(log(Odds_Ratio_DF$OR[Odds_Ratio_DF$rep_count==rep_count])+1.96*sqrt(1/ort[2]+1/ort[3]+1/ort[4]+1/ort[1]))
      Odds_Ratio_DF$ORlci[Odds_Ratio_DF$rep_count==rep_count]<-e^(log(Odds_Ratio_DF$OR[Odds_Ratio_DF$rep_count==rep_count])-1.96*sqrt(1/ort[2]+1/ort[3]+1/ort[4]+1/ort[1]))
      Odds_Ratio_DF$OR[Odds_Ratio_DF$OR==Inf] <- 0
      ort<-ort-corr # Un-correct counts for Fisher test etc
      if(as.numeric(rep_count)>min(as.numeric(names(table(Working_Dataset_1$A2))))) { 
        Odds_Ratio_DF$fisherman[Odds_Ratio_DF$rep_count==rep_count]<-fisher.test(ort)$p 
      }
    }
    # Apply bonferroni correction to p-values 
    Odds_Ratio_DF$fisherman <- Odds_Ratio_DF$fisherman < 0.05/(nrow(Odds_Ratio_DF)*length(unique(Working_Dataset$Gene)))
    # Convert odds ratios to log10
    Odds_Ratio_DF$OR <- log10(Odds_Ratio_DF$OR)
    Odds_Ratio_DF$ORuci <- log10(Odds_Ratio_DF$ORuci)
    Odds_Ratio_DF$ORlci <- log10(Odds_Ratio_DF$ORlci)
    blank_plot(0,max_columns,-3,3)
    # Grid Lines:  only want to plot these over the original distance of the plot, not the fake distance that keeps all the bars uniform 
    for (g in c(-3,-2,-1,1,2,3)){
      lines(x=c(-2,nrow(Odds_Ratio_DF)),y=c(g,g),lty=2,lwd=0.8,col="grey")
    }
    lines(x=c(-2,nrow(Odds_Ratio_DF)),y=c(0,0),lty=1,lwd=0.8,col="grey")
    j=0
    for (rep_count in Odds_Ratio_DF$rep_count){
      # plot vertical CI lines
      lines(x=c(j+0.4,j+0.4),y=c(Odds_Ratio_DF$ORuci[Odds_Ratio_DF$rep_count==rep_count],Odds_Ratio_DF$ORlci[Odds_Ratio_DF$rep_count==rep_count]),lwd=0.5,col=nicegrey,lty=1)
      # plot CI caps 
      lines(x=c(j+0.3,j+0.5),y=c(Odds_Ratio_DF$ORuci[Odds_Ratio_DF$rep_count==rep_count],Odds_Ratio_DF$ORuci[Odds_Ratio_DF$rep_count==rep_count]),lwd=0.5,col=nicegrey,lty=1)
      lines(x=c(j+0.3,j+0.5),y=c(Odds_Ratio_DF$ORlci[Odds_Ratio_DF$rep_count==rep_count],Odds_Ratio_DF$ORlci[Odds_Ratio_DF$rep_count==rep_count]),lwd=0.5,col=nicegrey,lty=1)
      points(x=j+0.4,y=Odds_Ratio_DF$OR[Odds_Ratio_DF$rep_count==rep_count],pch=21,cex=1.5,lwd=0.5,bg="grey",col=nicegrey)
      # x axis labels 
      #text(x=j,y=-3.5,paste("\u2265",rep_count,sep=""),xpd=NA,cex=0.6,col=nicegrey)
      text(x=j+0.4,y=-3.9,paste(">",rep_count-1,sep=""),xpd=NA,cex=0.6,col=nicegrey)
      # Add asterix if any are significant 
      if (isTRUE(Odds_Ratio_DF$fisherman[Odds_Ratio_DF$rep_count==rep_count& !is.na(Odds_Ratio_DF$fisherman)])){
        text(x=j+0.05*length(Odds_Ratio_DF$rep_count),y=1.7,"*",pos=3,cex=2,xpd=NA,col=nicegrey)
      }
      j=j+1
    }
    # Plot axis
    axis(2,at=c(-3:3),lab=c(-3:3),las=1,lwd=1,col=nicegrey,col.axis=nicegrey,cex.axis=0.9)
    title(ylab=expression(paste(log[10](OR))),line=4.5,cex.lab=0.9,col.lab="grey20",xpd=NA)
    title(ylab="(95% CI)",line=3.5,cex.lab=0.9,col.lab="grey20",xpd=NA)
    # Plot "repeats:"
    text(x=-1.7,y=-4,"repeats: ",las=1,col=nicegrey,cex=0.9,xpd=NA)
    #axis(1,at=-1,lab="repeats >= : ",las=1,col=nicegrey,col.axis=nicegrey,cex.axis=0.9,hadj=1,line=NA,tick=F,xpd=NA)
    i=i+1
  }
  dev.off()
}
print("Plot Chapter3_REs_BP_OR Complete")

################################################################################
################################################################################

# Figure 3.11 : Exploration of TREDPARSE Predicted TCF4 Expansions

################################################################################
################################################################################

# PDF = 210x297 - margins are 35,20,20,20 mm 
# PDF = 8.3x11.7 - margins are 1.38,0.79,0.79,0.79 inches
pdf(paste(output_directory,"Chapter3_TREDPARSE_TCF4_Exploration.pdf",sep=""),width=(8.3-1.38-0.8),height=(11.7-0.8-0.8))
par(mfrow=c(4,2),mar=c(4,5,3.5,1))

#######
# Compare to EHV3 Results 
#######
# Want to check if coverage is sufficient to rely on WES results 
# Find the samples which tredparse calls above 60
tredparse_long_samples <- tredparse$Sample[tredparse$A2>60 & tredparse$Gene=="TCF4" & tredparse$Dataset=="RCSI_WGS_PCR_FREE"]
#
EHV3_allele_size <- c(EHV3$A1[EHV3$Sample %in% tredparse_long_samples & EHV3$Gene=="TCF4"],EHV3$A2[EHV3$Sample %in% tredparse_long_samples & EHV3$Gene=="TCF4"])
tredparse_allele_size <- c(tredparse$A1[tredparse$Sample %in% tredparse_long_samples & tredparse$Gene=="TCF4"],tredparse$A2[tredparse$Sample %in% tredparse_long_samples & tredparse$Gene=="TCF4"])
# Create blank plot 
blank_plot(0,200,0,100)
# Grid Lines 
lines(x=c(-5,250),y=c(0,0),lty=2,lwd=0.8,col="grey")
lines(x=c(-5,250),y=c(50,50),lty=2,lwd=0.8,col="grey")
lines(x=c(-5,250),y=c(100,100),lty=2,lwd=0.8,col="grey")
# 45 degree line
abline(0,1,lwd=1,col=nicegrey)
# plot points 
points(x=tredparse_allele_size,y=EHV3_allele_size,pch=21,cex=1.2,bg=c(rep("grey",length(EHV3_allele_size)/2),rep(niceorange,length(EHV3_allele_size)/2)),col=NA,lwd=0.01)
# Axes
axis(1,at=c(0,50,100,150,200),lab=c(0,50,100,150,200),las=1,lwd=1,col=nicegrey,col.axis=nicegrey,cex.axis=1)
axis(2,at=c(0,25,50,75,100),lab=c(0,25,50,75,100),las=1,lwd=1,col=nicegrey,col.axis=nicegrey,cex.axis=1)
box(lwd=1.2,col=nicegrey)
# Title for y-axis 1 
title(xlab=expression(paste("TREDPARSE ",italic("TCF4")," Allele",sep="")),ylab=expression(paste("ExpansionHunter v3 ",italic("TCF4")," Allele",sep="")), cex.lab=1,col.lab=nicegrey,line=2.5,font=1.1,xpd=NA)
title(main=expression(paste("ExpansionHunter / TREDPARSE ",italic("TCF4"),sep="")),cex.main=1, cex.lab=1,col.main=nicegrey,line=1.5,font=1.1,xpd=NA)
# Legend 
points(x=0,y=90,pch=21,bg="grey",col="white",lwd=0.1,cex=1)
points(x=0,y=80,pch=21,bg=niceorange,col="white",lwd=0.1,cex=1)
text(x=5,y=90,"Shorter Allele",cex=0.7,col=nicegrey,adj=0)
text(x=5,y=80,"Longer Allele",cex=0.7,col=nicegrey,adj=0)
plot_letter_label(0,200,0,100,"A")
#######
# TCF4 Exome Coverage 
#######

# Set col names
colnames(TCF4_coverage) <- c("TCF4_coverage","WES_coverage","Sample")
# Set colour 
TCF4_coverage$Colour <- "grey"
TCF4_coverage$Colour[TCF4_coverage$Sample %in% tredparse_long_samples] <- niceorange
# Create blank plot 
blank_plot(0,125,0,75)
# Grid Lines 
lines(x=c(-5,150),y=c(25,25),lty=2,lwd=0.8,col="grey")
lines(x=c(-5,150),y=c(50,50),lty=2,lwd=0.8,col="grey")
lines(x=c(-5,150),y=c(75,75),lty=2,lwd=0.8,col="grey")
abline(0,1,lwd=1,col=nicegrey)
# Plot points 
points(x=TCF4_coverage$WES_coverage,y=TCF4_coverage$TCF4_coverage,pch=21,cex=1.2,bg=TCF4_coverage$Colour,col=NA,lwd=0.01)
# Want highlighted points on top 
points(x=TCF4_coverage$WES_coverage[TCF4_coverage$Sample %in% tredparse_long_samples],y=TCF4_coverage$TCF4_coverage[TCF4_coverage$Sample %in% tredparse_long_samples],pch=21,cex=1.2,bg=TCF4_coverage$Colour[TCF4_coverage$Sample %in% tredparse_long_samples],col=NA,lwd=0.1)
# Axes
axis(1,at=c(0,25,50,75,100,125),lab=c(0,25,50,75,100,125),las=1,lwd=1,col=nicegrey,col.axis=nicegrey,cex.axis=1)
axis(2,at=c(0,25,50,75,100,125),lab=c(0,25,50,75,100,125),las=1,lwd=1,col=nicegrey,col.axis=nicegrey,cex.axis=1)
box(lwd=1.2,col=nicegrey)
# Title for y-axis 1 
title(xlab="Sample Mean Exome Coverage (X)",ylab=expression(paste("Sample Mean ",italic("TCF4")," Coverage (X)",sep="")), cex.lab=1,col.lab=nicegrey,line=2.5,font=1.1,xpd=NA)
title(main=expression(paste("Epilepsy Coverage at ",italic("TCF4")," Locus",sep="")),cex.main=1, cex.lab=1,col.lab=nicegrey,line=1.5,font=1.1,xpd=NA)
# Legend 
points(x=0,y=70,pch=21,bg="grey",col="white",lwd=0.1,cex=1)
points(x=0,y=65,pch=21,bg=niceorange,col="white",lwd=0.1,cex=1)
text(x=5,y=70,"All Epilepsy Samples with WES",cex=0.7,col=nicegrey,adj=0)
text(x=5,y=65,expression(paste("TREDPARSE predicted significants",sep="")),,cex=0.7,col=nicegrey,adj=0)
plot_letter_label(0,125,0,75,"B")



#######
# exSTRa ecdf
#######
# Create a plot of the ecdf values of the exSTRa output for TCF4 for WES samples
# Set Colours
TCF4_exSTRa_WES$Colour <- "grey"
TCF4_exSTRa_WES$Colour[TCF4_exSTRa_WES$Sample %in% tredparse_long_samples] <- niceorange

# Set blank plot 
blank_plot(
  xmin=0
  ,xmax=150
  ,ymin=0
  ,ymax=1
)
# Plot grey samples 
for (sample in sample(unique(TCF4_exSTRa_WES$Sample[TCF4_exSTRa_WES$Colour=="grey"]))){
  colour_to_plot <- TCF4_exSTRa_WES$Colour[TCF4_exSTRa_WES$Sample==sample]
  sample_ecdf <- ecdf(sort(TCF4_exSTRa_WES$rep[TCF4_exSTRa_WES$Sample==sample]))
  # Plot gridlines
  lines(x=c(-17,150),y=c(25,25),lty=2,lwd=0.8,col="grey",xpd=NA)
  lines(x=c(-17,150),y=c(50,50),lty=2,lwd=0.8,col="grey",xpd=NA)
  lines(x=c(-17,150),y=c(75,75),lty=2,lwd=0.8,col="grey",xpd=NA)
  lines(x=c(-17,150),y=c(100,100),lty=2,lwd=0.8,col="grey",xpd=NA)
  plot(sample_ecdf, verticals=TRUE, do.points=FALSE,add=TRUE,col=colour_to_plot,lwd=0.4,col.01line = NULL)
}
# Plot orange samples
for (sample in sample(unique(TCF4_exSTRa_WES$Sample[TCF4_exSTRa_WES$Colour==niceorange]))){
  colour_to_plot <- TCF4_exSTRa_WES$Colour[TCF4_exSTRa_WES$Sample==sample]
  sample_ecdf <- ecdf(sort(TCF4_exSTRa_WES$rep[TCF4_exSTRa_WES$Sample==sample]))
  # Plot gridlines
  lines(x=c(-17,150),y=c(25,25),lty=2,lwd=0.8,col="grey",xpd=NA)
  lines(x=c(-17,150),y=c(50,50),lty=2,lwd=0.8,col="grey",xpd=NA)
  lines(x=c(-17,150),y=c(75,75),lty=2,lwd=0.8,col="grey",xpd=NA)
  lines(x=c(-17,150),y=c(100,100),lty=2,lwd=0.8,col="grey",xpd=NA)
  plot(sample_ecdf, verticals=TRUE, do.points=FALSE,add=TRUE,col=colour_to_plot,lwd=0.4,col.01line = NULL)
}
# plot white rectangle to cover edges
rect(xleft=155,xright=200,ytop=1.1,ybottom=0.9,col="white",border=NA,lwd=0.8,xpd=NA)
rect(xleft=-25,xright=-5,ytop=0.1,ybottom=-0.1,col="white",border=NA,lwd=0.8,xpd=NA)
#lines(x=c(-17,-17),y=c(0,1),lwd=1,col=nicegrey,xpd=NA)
#lines(x=c(0,150),y=c(-0.1,-0.1),lwd=1,col=nicegrey,xpd=NA)
# Axes
axis(1,at=c(0,150),lab=c(0,150),las=1,lwd=1,col=nicegrey,col.axis=nicegrey,cex.axis=1)
axis(2,at=c(0,0.25,0.5,0.75,1),lab=c(0,0.25,0.5,0.75,1),las=1,lwd=1,col=nicegrey,col.axis=nicegrey,cex.axis=1)
box(lwd=1.2,col=nicegrey)
# Title for y-axis 1 
title(xlab="Repeated Bases", cex.lab=1,col.lab=nicegrey,line=2.5,font=1.1,xpd=NA)
title(ylab="ecdf(Repeated Bases)", cex.lab=1,col.lab=nicegrey,line=3.5,font=1.1,xpd=NA)
title(main=expression(paste("exSTRa: ",italic("TCF4"),sep="")), cex.main=1,cex.lab=1,col.main=nicegrey,line=1.5,font=1.1,xpd=NA)
# Plot Legend
points(x=60,y=0.25,pch=21,bg="grey",col="white",lwd=0.1,cex=1)
points(x=60,y=0.20,pch=21,bg=niceorange,col="white",lwd=0.1,cex=1)
text(x=65,y=0.25,"All Epilepsy Samples with WES",cex=0.7,col=nicegrey,adj=0)
text(x=65,y=0.20,"TREDPARSE Significant",cex=0.7,col=nicegrey,adj=0,xpd=NA)
plot_letter_label(0,150,0,1,"C")

#######
# exSTRa P Values
#######
# Create a QQ Plot of TCF4 in exSTRa Exome Samples 
# Filter to Data of Interest
exSTRa_TCF4<-subset(exSTRa[exSTRa$Gene=="TCF4" & exSTRa$Dataset=="RCSI_WES",],select=c(Sample,P_Value))
# Sort by p-value (v important of colours will be mislabeled)
exSTRa_TCF4<-exSTRa_TCF4[order(exSTRa_TCF4$P_Value),]
# Assign Colours 
exSTRa_TCF4$Colour <- "grey"
exSTRa_TCF4$Colour[exSTRa_TCF4$Sample %in% tredparse_long_samples] <- niceorange
# Our Observed Values 
o = -log10(sort(exSTRa_TCF4$P_Value,decreasing=FALSE))
# Our Expected Values 
e = -log10( ppoints(length(exSTRa_TCF4$P_Value) ))
# Create blank plot 
blank_plot(0,4,0,4)
# 45 degree abline 
abline(a=0,b=1,col="grey20",lty=1,lwd=1.2)
# Plot points 
points(e,o,pch=21,lwd=0.01,col=NA,bg=exSTRa_TCF4$Colour,cex=1.2)
# p value threshold line
lines(x=c(-1,5),y=c(-log10(0.05/34),-log10(0.05/34)),lty=2,lwd=0.8,col=nicered)
text(x=0,y=3,"p-value threshold",cex=0.8,col=nicegrey,adj=0,xpd=NA)
# Axes
axis(1,at=c(0,-log10(0.1),-log10(0.01),-log10(0.001),-log10(0.0001)),lab=c(0,1,2,3,4),las=1,lwd=1,col=nicegrey,col.axis=nicegrey,cex.axis=1)
axis(2,at=c(0,-log10(0.1),-log10(0.01),-log10(0.001),-log10(0.0001)),lab=c(0,1,2,3,4),las=1,lwd=1,col=nicegrey,col.axis=nicegrey,cex.axis=1)
box(lwd=1.2,col=nicegrey)
# Title for y-axis 1 
title(xlab=expression('Expected -log'[10]*' (p-value)'),ylab=expression('Observed -log'[10]*' (p-value)'), cex.lab=1,col.lab=nicegrey,line=2.5,font=1.1,xpd=NA)
title(main=expression(paste("exSTRa ",italic("TCF4")," QQ Plot",sep="")), cex.main=1,cex.lab=1,col.main=nicegrey,line=1.5,font=1.1,xpd=NA)
# Legend 
points(x=0,y=3.9,pch=21,bg="grey",col="white",lwd=0.1,cex=1)
points(x=0,y=3.6,pch=21,bg=niceorange,col="white",lwd=0.1,cex=1)
text(x=0.1,y=3.9,"All Epilepsy Samples with WES",cex=0.7,col=nicegrey,adj=0)
text(x=0.1,y=3.6,"TREDPARSE Significant",cex=0.7,col=nicegrey,adj=0)
plot_letter_label(0,4,0,4,"D")
dev.off()


################################################################################
################################################################################

# Figure 3.12 : exploration of loci predicted significant by STRetch

################################################################################
################################################################################
#########################################
# STRetch WGS PCR v EHv3 Controls 
#########################################
# A plot to compare the 'significant' loci identified by STRetch to Irish WGS PCR FREE controls 
# NOTE: no significant loci were identified from the WES data and the whole genome PCR free data was not tested as there were no controls 
STR_EH_comp <- subset(STRetch[STRetch$Dataset=="RCSI_WGS_PCR",],select=c(Dataset,A2,Gene))
EH_STR_comp <- subset(EHV3[EHV3$Dataset=="ALS_WGS_PCR_FREE" & EHV3$CC==0 & !is.na(EHV3$CC) & EHV3$Gene %in% unique(STR_EH_comp$Gene),],select=c(Dataset,A2,Gene))
Working_Dataset <- rbind(STR_EH_comp,EH_STR_comp)
max_columns <- c()
for (gene in sort(unique(Working_Dataset$Gene))){
  max_columns <- c(max_columns,length(unique(Working_Dataset$A2[Working_Dataset$Gene==gene])))
  max_columns <- c(max_columns,length(unique(Working_Dataset$A2[Working_Dataset$Gene==gene])))
}
max_columns <- max(max_columns)
# Add case control status 
Working_Dataset$Status<- gsub("_.*","",Working_Dataset$Dataset)
Working_Dataset$Status<- gsub("RCSI","Patient",Working_Dataset$Status)
Working_Dataset$Status<- gsub("ALS","Control",Working_Dataset$Status)


# For the barplots we want to precreate the barplots to get the positions we need 
# Create DF for each gene 
for (gene in sort(unique(Working_Dataset$Gene))){
  assign(paste(gene,"_working_dataframe",sep=""),create_gene_table(temp_df=Working_Dataset,gene=gene))
}
# Precreate each barplot for first type of plots - need this to get x positions 
for (gene in sort(unique(Working_Dataset$Gene))){
  assign(paste(gene,"_bp1",sep=""),create_bp(gene))
}

# Assign new dataframes including blank spaces 
for (gene in sort(unique(Working_Dataset$Gene))){
  assign(paste(gene,"_working_dataframe",sep=""),bring_dfs_to_max_columns(gene=gene,max_columns=max_columns))
}

# PDF = 210x297 - margins are 35,20,20,20 mm 
# PDF = 8.3x11.7 - margins are 1.38,0.79,0.79,0.79 inches
pdf(paste(output_directory,"Chapter3_STR_RCSI_WGS_PCR_SIGNIFICANTS.pdf",sep=""),width=(11.7-0.8-0.8),height=(8.3-1.38-0.8))
par(mfrow=c(4,1))
i=1
# Cycle through each gene 
for (gene in sort(unique(Working_Dataset$Gene))) {
  par(mar=c(2.3,8,3.5,1))
  # Filter to gene of interest
  Working_Dataset_1 <- Working_Dataset[Working_Dataset$Gene==gene,]

  #############
  # First plot - barplot 
  ############

  # Fetch the repeat motif 
  #repeat_motif=repeat_motifs[match(gene,gene_list_quotes)]
  repeat_motif=repeat_motifs[match(gene,gene_list_quotes)]
  # Plot barplot
  # Remove border on bars so that there is no zero line at the bottom of the plot 
  par(lty = 0)
  # For most genes setting the ylimit at 75 is visually best but some require higher - have found these through trial and error 
  ymax=100
  # Calculate barplot 
  bp <- barplot(get(paste(gene,"_working_dataframe",sep="")),
    beside=T,
    ylim=c(0,ymax),         
    col=c(nicedarkblue,nicered),
    #lwd=2,
    yaxt="n",
    xaxt="n",
    space=c(0,0.35),
    cex.names=0.6,
    col.axis=nicegrey,
    xlab="",
    cex.lab=0.1,
    col.lab=nicegrey,
    #border=c("white","white"),
    #border=c(nicedarkblue,nicered)
  )
  # Title for y-axis 1 
  title(ylab=bquote(italic(.(gene)) ~ " allele carrier frequency,"),cex.lab=0.9,col.lab=nicegrey,line=4.5,font=3,xpd=NA)
  # Title for y-axis 2
  title(ylab="longer allele (%)",cex.lab=0.9,col.lab=nicegrey,line=3.5,font=3,xpd=NA)
  # Plot Main Title
  text(x=min(bp),y=ymax*1.1,bquote(bold(italic(.(gene)) ~ " " ~ .(repeat_motif) ~ " repeats   (" ~ .("STRetch Cases and ExpansionHunterv3 Controls") ~ ")")),adj=0,cex=1,col=nicegrey,xpd=TRUE,font=2)
  # Grid Lines:  only want to plot these over the original distance of the plot, not the fake distance that keeps all the bars uniform 
  lines(x=c(-5,max(get(paste(gene,"_bp1",sep="")))),y=c(25,25),lty=2,lwd=0.8,col="grey")
  lines(x=c(-5,max(get(paste(gene,"_bp1",sep="")))),y=c(50,50),lty=2,lwd=0.8,col="grey")
  lines(x=c(-5,max(get(paste(gene,"_bp1",sep="")))),y=c(75,75),lty=2,lwd=0.8,col="grey")
  lines(x=c(-5,max(get(paste(gene,"_bp1",sep="")))),y=c(100,100),lty=2,lwd=0.8,col="grey")
  # Replot bars over grid lines
  barplot(get(paste(gene,"_working_dataframe",sep="")),
    beside=T,
    ylim=c(0,ymax),         
    col=c(nicedarkblue,nicered),
    #lwd=2,
    yaxt="n",
    space=c(0,0.35),
    cex.names=0.6,
    col.axis=nicegrey,
    xlab="",
    cex.lab=0.1,
    col.lab=nicegrey,
    #border=c(nicedarkblue,nicered),
    #border=c("white","white"),
    add=T
  )
  # white line at 0 to hide 0 values (that otherwise plot as tiny thin lines)
  abline(h=0,lty=1,lwd=0.1,col="white")
  # Plot axis
  axis(2,at=c(0,25,50,75,100),lab=c(0,25,50,75,100),las=1,lwd=1,col=nicegrey,col.axis=nicegrey,cex.axis=0.9)
  # Plot "repeats:"
  axis(1,at=bp[1],lab="repeats: ",las=1,col=nicegrey,col.axis=nicegrey,cex.axis=0.9,hadj=1,line=NA,tick=F)
  # Plot legend on the first plot on each page 
  # Plot Coloured Points
  points(x=c(max(get(paste(gene,"_bp1",sep="")))*0.9,max(get(paste(gene,"_bp1",sep="")))*0.9),y=c(83,67),col=c(nicered,nicedarkblue),bg=c(nicered,nicedarkblue),pch=21,cex=1.4,lwd=0.6)
  # Plot Case/Control Labels
  text(x=c(max(get(paste(gene,"_bp1",sep="")))*0.9+max(bp)*0.01,max(get(paste(gene,"_bp1",sep="")))*0.9+max(bp)*0.01),y=c(83,67),c("Patients","Controls"),col=nicegrey,adj=0,cex=0.9)
  # Plot Letter 
  if (i <= 26){
    plot_letter=LETTERS[i]
  } else {
    j=i
    j=j-26
    plot_letter=paste(LETTERS[j],".2",sep="")
  }
  plot_letter_label(min(bp),max(bp),0,ymax,plot_letter)
  i=i+1
}
dev.off()
print("Plot Chapter3_STR_RCSI_WGS_PCR_SIGNIFICANTS.pdf Complete")


################################################################################
################################################################################

# Figure 3.13 , SXXX-XXX : Potential de novo REs in Epilepsy Patients

################################################################################
################################################################################

# This makes a plot comparing the longest index alleles to their parent alleles to find any de novo REs in target genes 
# It separately checks the WGS PCR data (9 trios) and the WES data (XXX trios)

plot_parent_proband_comparisons=function(proband_A2s,parents_A2s,gene,i,j,software,dataset,fail_genes){
  if (dataset=='RCSI_WGS_PCR'){
    dataset_literal <- "WGS PCR"
  } else {
    dataset_literal <- "WES"
    # Create blank plots so WES appears on new line 
    if (i==1 & j %% 4 == 0){
        blank_plot(0,1,0,1)
    } else if (i==1 & j %% 4 == 3){
      blank_plot(0,1,0,1)
      blank_plot(0,1,0,1)
    } else if (i==1 & j %% 4 == 2){
      blank_plot(0,1,0,1)
      blank_plot(0,1,0,1)
      blank_plot(0,1,0,1)
    }
  }
  label_locations <- sort(unique(c(proband_A2s,parent_A2s)))
  if (length(label_locations)<2){
    label_locations <- c(label_locations-1,label_locations,label_locations+1)
  } 
  blank_plot(xmin=min(label_locations),xmax=max(label_locations),ymin=min(label_locations),ymax=max(label_locations))
  # Plot 45 degree line across plot 
  lines(x=c(0.5*min(label_locations),1.5*max(label_locations)),y=c(0.5*min(label_locations),1.5*max(label_locations)),col=nicegrey,lwd=1,lty=1)
  # Add individual plot titles 
  title(main=bquote(italic(.(gene))), cex.main=1.1,cex.lab=1,col.lab=nicegrey,xpd=NA,line=0.7,font=3)
  # Add axis labels 
  title(ylab="Proband Longest Allele", xlab="Parent Corresponding Allele", cex.main=1.1,cex.lab=1,col.lab=nicegrey,xpd=NA,line=2.2)
  # Draw external box 
  box(lwd=2,col=nicegrey)
  # plot x axis ticks 
  #axis(side=1,lwd=1,cex.axis=1,col.axis=nicegrey,col=nicegrey,at=axis_tick_locations)
  axis(side=1,lwd=1,cex.axis=1,col.axis=nicegrey,col=nicegrey,at=label_locations,labels=label_locations)
  # plot y axis ticks 
  #axis(side=2,lwd=1,las=1,cex.axis=1,col.axis=nicegrey,col=nicegrey,at=axis_tick_locations)
  axis(side=2,lwd=1,las=1,cex.axis=1,col.axis=nicegrey,col=nicegrey,at=label_locations,labels=label_locations)
  # plot points giving each a unique colour to correlate with other plots 
  points(x=proband_A2s,y=parent_A2s,col=see_through_colour(hexcode=colorRampPalette(vangogh_palette)(length(gene_list_quotes))[match(gene,gene_list_quotes)],fade_proportion=0.9),bg=see_through_colour(hexcode=colorRampPalette(vangogh_palette)(length(gene_list_quotes))[match(gene,gene_list_quotes)],fade_proportion=0.4),pch=21,cex=1.2,xpd=NA,lwd=0.3)
  # Plot title  above first plot 
  if(i==1){
      text(x=min(label_locations),y=(max(label_locations)+0.4*(max(label_locations)-min(label_locations))),paste(software,": Comparison of Longest Allele in Proband & Corresponding Parent Allele from ",dataset_literal," Data",sep=""),cex=1.5,adj=0,col=nicegrey,xpd=NA,font=2)
  }
  # Add an asterix if gene is found to have failed the RMSD check earlier
  if (gene %in% fail_genes & dataset_literal=="WES"){
    text(x=(max(label_locations)),y=max(label_locations)+0.15*(max(label_locations)-min(label_locations)),"*",cex=3.5,adj=0,col=nicered,xpd=NA,font=2)
  } 
}


for (primary_dataset in vector_of_datasets){
  if (primary_dataset!='STRetch'){
    # This software name is used for the plot title and the name the output PDF 
    software=vector_of_software_names[match(primary_dataset,vector_of_datasets)]

    # PDF = 210x297 - margins are 35,20,20,20 mm 
    # PDF = 8.3x11.7 - margins are 1.38,0.79,0.79,0.79 inches
    pdf(paste(output_directory,"Chapter3_Parent_Proband_Comparison_",software,".pdf",sep=""),width=(11.7-0.8-0.8),height=(8.3-1.38-0.8))
    par(mfrow=c(3,4),mar=c(3.5,3.5,4.5,2))
    j=1
    for (dataset in c("RCSI_WGS_PCR","RCSI_WES")){
      i=1
      # Filter dataset to the WGS PCR data - this contains 11 probands of which 9 are trios 
      primary_dataset_trios <- get(paste(primary_dataset,"_pre_exclusion",sep=""))[get(paste(primary_dataset,"_pre_exclusion",sep=""))$Dataset==dataset,]
      # Filter to the 9 trios 
      # Find samples which are trios 
      samples_present <- gsub("A|B|C","",unique(primary_dataset_trios$Sample))
      #trio_samples <- c(paste(names(table(samples_present )[table(samples_present )==3]),"A",sep=""),paste(names(table(samples_present )[table(samples_present )==3]),"B",sep=""),paste(names(table(samples_present )[table(samples_present)==3]),"C",sep=""))
      trio_samples <- names(table(samples_present)[table(samples_present )==3])
      # Retain just these 
      primary_dataset_trios <- primary_dataset_trios[primary_dataset_trios$Sample %in% c(paste(trio_samples,"A",sep=""),paste(trio_samples,"B",sep=""),paste(trio_samples,"C",sep="")),]
      for (gene in sort(unique(primary_dataset_trios$Gene))) {
        proband_A2s <- c() 
        parent_A2s <- c()
        for (sample in trio_samples){
          # The proband allele is the longest allele - this may not be the longest allele in the parents.
          # Checking to see if the parents have this allele and if not, then taking the longest allele in the parents and assuming expanded from this allele 
          # If all alleles are called 
          if (length(primary_dataset_trios$A2[primary_dataset_trios$Gene==gene & primary_dataset_trios$Sample==paste(sample,"A",sep="")]) >0 & length(primary_dataset_trios$A1[primary_dataset_trios$Gene==gene & primary_dataset_trios$Sample==paste(sample,"B",sep="")]) >0 & length(primary_dataset_trios$A1[primary_dataset_trios$Gene==gene & primary_dataset_trios$Sample==paste(sample,"C",sep="")]) >0){
            # Get longest proband allele 
            proband_A2 <- primary_dataset_trios$A2[primary_dataset_trios$Gene==gene & primary_dataset_trios$Sample==paste(sample,"A",sep="")]
            proband_A2s <- c(proband_A2s,proband_A2) 
            # Get all parent alleles 
            parent_alleles <- c(primary_dataset_trios$A1[primary_dataset_trios$Gene==gene & primary_dataset_trios$Sample==paste(sample,"B",sep="")],primary_dataset_trios$A2[primary_dataset_trios$Gene==gene & primary_dataset_trios$Sample==paste(sample,"B",sep="")],primary_dataset_trios$A1[primary_dataset_trios$Gene==gene & primary_dataset_trios$Sample==paste(sample,"C",sep="")],primary_dataset_trios$A2[primary_dataset_trios$Gene==gene & primary_dataset_trios$Sample==paste(sample,"C",sep="")])
            # If one of hte parents carry the longest allele this is straightforward 
            if (proband_A2 %in% parent_alleles){
              parent_allele <- proband_A2 
            } else {
              # If neither parent carriers the longest allele we want to exclude A1 and find the nearest allele (i.e. the most likely to be the parent allele)
              proband_A1 <- primary_dataset_trios$A1[primary_dataset_trios$Gene==gene & primary_dataset_trios$Sample==paste(sample,"A",sep="")]
              # If A1 is in parent alleles one time we need to remove it 
              if (length(parent_alleles[parent_alleles==proband_A1])==1){
                parent_alleles<- parent_alleles[parent_alleles!=proband_A1]
              }
              parent_allele <- parent_alleles[which.min(abs(parent_alleles - proband_A2))]
            }
            parent_A2s <- c(parent_A2s,parent_allele)
          }
        }
        if (length(proband_A2s)>0){
          plot_parent_proband_comparisons(proband_A2s=proband_A2s,parents_A2s=parents_A2s,gene=gene,i=i,j=j,software=software,dataset=dataset,fail_genes=get(paste(primary_dataset,"_fail_genes",sep="")))
        i=i+1
        j=j+1
        }
      }
    }
  dev.off()
  }
}

print("Plot Chapter3_Parent_Proband_Comparison Complete")




################################################################################
################################################################################

# Supplementary Figure S3.24 : Coverage of Samples with exSTRa Predicted NOTCH2 REs

################################################################################
################################################################################

#############
# NOTCH2 Coverage Plot 
#############
# PDF = 210x297 - margins are 35,20,20,20 mm 
# PDF = 8.3x11.7 - margins are 1.38,0.79,0.79,0.79 inches
pdf(paste(output_directory,"Chapter3_NOTCH2_WES_COVERAGE_Exploration.pdf",sep=""),width=(8.3-1.38-0.8),height=(11.7-0.8-0.8))
par(mfrow=c(4,2),mar=c(4,5,3.5,1))

# Create blank plot 
blank_plot(0,120,0,1200)
# 45 degree abline 
abline(lm(WES_NOTCH2_Coverage$V2 ~ WES_NOTCH2_Coverage$V1),col=nicered,lty=1,lwd=1.2)
points(x=WES_NOTCH2_Coverage$V1,y=WES_NOTCH2_Coverage$V2,pch=21,bg="grey",col=NA,lwd=0.1,cex=1)
# Have manually checked that these are the correct coverages 
points(x=c(110.399,89.8823,73.7135),y=c(623,525,361),pch=21,bg="orange",col=NA,lwd=0.1,cex=1)
# Axes
axis(1,at=c(0,40,80,120),lab=c(0,40,80,120),las=1,lwd=1,col=nicegrey,col.axis=nicegrey,cex.axis=1)
axis(2,at=c(0,400,800,1200),lab=c(0,400,800,1200),las=1,lwd=1,col=nicegrey,col.axis=nicegrey,cex.axis=1)
box(lwd=1.2,col=nicegrey)
# Title for y-axis 1 
title(xlab='Mean WES Coverage (X)',ylab=expression(paste(italic("NOTCH2")," Coverage (X)",sep="")), cex.lab=1,col.lab=nicegrey,line=2.5,font=1.1,xpd=NA)
title(main=expression(paste("Coverage of ",italic("NOTCH2")," Repeat in WES Samples",sep="")), cex.main=1,cex.lab=1,col.main=nicegrey,line=1.5,font=1.1,xpd=NA)
# Legend 
points(x=0,y=1150,pch=21,bg="grey",col="white",lwd=0.1,cex=1)
points(x=0,y=1080,pch=21,bg=niceorange,col="white",lwd=0.1,cex=1)
text(x=5,y=1150,"All WES Samples",cex=0.7,col=nicegrey,adj=0)
text(x=5,y=1080,"exSTRa Predicted Significant",cex=0.7,col=nicegrey,adj=0)
lines(x=c(0,3),y=c(1010,1010),col=nicered,lty=1,lwd=1.2)
text(x=5,y=1010,"Regression",cex=0.7,col=nicegrey,adj=0)
dev.off()

