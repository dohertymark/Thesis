################################################################################
################################################################################
# A script to plot the Irish, global and European freuqencies for ALS and FTD variants
# email: dohertm7@tcd.ie
################################################################################
################################################################################

rm(list=ls())

################################################################################
################################################################################

# Define Directories 

################################################################################
################################################################################

data_directory=""
output_directory=""

################################################################################
################################################################################

# Load packages

################################################################################
################################################################################

library(data.table)
library(tidyr)
library(binom)
library(berryFunctions)

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
#vangogh_palette<<-c(nicered,hull,peach,sunset,niceorange,lightgreen,darkgreen,reallylightblue,nicelightblue,nicedarkblue,reallydarkblue,pink,green,lime,medium_green,brown,bright_yellow,lavendar,algae,yellow1,purple,light_purple,pink_purple)
vangogh_palette<<-c(nicered,hull,peach,sunset,niceorange,lightgreen,darkgreen,reallylightblue,nicelightblue,nicedarkblue,reallydarkblue,pink,green,lime,brown,bright_yellow,lavendar,algae,yellow1,purple,light_purple,pink_purple)

#####
# FTD Individuals Data 
#####

ftd_population_studies      <<- fread(paste(data_directory"/ftd_population_studies.txt"),header=TRUE,sep='\t',na.strings=c("","NA"))
ftd_population_people_dataset <<- fread(paste(data_directory"/ftd_population_people_dataset.txt"),header=TRUE,sep='\t',na.strings=c("","NA"))

#####
# ALS Individuals Data 
#####

als_population_studies      <<- fread(paste(data_directory"/als_population_studies.txt"),header=TRUE,sep='\t',na.strings=c("","NA"))
als_population_people_dataset <<- fread(paste(data_directory"/als_population_people_dataset.txt"),header=TRUE,sep='\t',na.strings=c("","NA"))

#####
# journALS Data 
#####

# Read in journALS Population Studies 
population_studies      <<- fread(paste(data_directory"/journALS_population_studies_long_format.tsv"),header=TRUE,sep='\t',na.strings=c("","NA"))

# Read in journALS Data 
lit_review_pen<- fread(paste(data_directory"/journALS_post_analysis.tsv.gz"),header=TRUE,sep='\t',na.strings=c("","NA"))

# Age P-Value Cutoff
age_p_cutoff <<-  data.frame(fread(paste(data_directory"/rev14_p_cutoff_df.tsv"),header=TRUE,sep='\t',na.strings=c("","NA")))

age_p_cutoff <<- age_p_cutoff$p_cutoff
# All ages df
#all_ages_df <<- data.frame(fread(paste(data_directory"/journALS_journALS_all_ages_df.tsv"),header=TRUE,sep='\t',na.strings=c("","NA")))
all_ages_df <<- data.frame(fread(paste(data_directory"/journALS/Datafiles/rev14_journALS_all_ages_df.tsv"),header=TRUE,sep='\t',na.strings=c("","NA")))

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
#variants_in_p_or_lp_genes<<-unique(c(lit_review_pen$HGVS[lit_review_pen$population_carriers_count>0 & !is.na(lit_review_pen$population_carriers_count) & lit_review_pen$gene %in% p_or_lp_genes & !is.na(lit_review_pen$gene)],ftd_population_people_dataset$HGVS,als_population_people_dataset$HGVS))
variants_in_p_or_lp_genes<<-unique(c(lit_review_pen$HGVS[lit_review_pen$population_carriers_count>0 & !is.na(lit_review_pen$population_carriers_count) & lit_review_pen$gene %in% p_or_lp_genes & !is.na(lit_review_pen$gene)],als_population_people_dataset$HGVS[als_population_people_dataset$gene %in% p_or_lp_genes]))
list_of_pathogenic_variants<<-list(pathogenic_variants,p_or_lp_variants,variants_in_p_or_lp_genes)

#####
# Define useful functions
#####

# A function to create blank plots, can then plot exactly what is wanted on top
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
  text(c(left,left+0.25*(right-left),left+0.5*(right-left),left+0.75*(right-left),right),5,c("0","0.1","1","10","100%") ,col=nicegrey,adj=0.5,cex=0.8,xpd=NA)
}

#####
# Parse Data 
#####


population_people_dataset <<- data.frame(population_people_dataset_function(lit_review_pen))
# The current existing Cuba data in journALS is a subset of this data 
population_people_dataset <<- population_people_dataset[population_people_dataset$population_carriers_nationality!="Ireland",]
population_studies <<- population_studies[population_studies$Country != "Ireland",]
# Get gene for Cuban Data 
ftd_population_people_dataset$gene <- gsub(":.*","",ftd_population_people_dataset$HGVS)
# Merge Back in new Cuban Data 
population_studies <<- rbind(population_studies,ftd_population_studies)
ftd_population_people_dataset <- subset(ftd_population_people_dataset,select=colnames(population_people_dataset))
population_people_dataset <<- rbind(population_people_dataset,ftd_population_people_dataset)
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
#NAM_countries <- unique(population_studies$Country[population_studies$Continent=="NAM"])[unique(population_studies$Country[population_studies$Continent=="NAM"]) !="Cuba"]
#SAM_countries <- unique(population_studies$Country[population_studies$Continent=="SAM"])[unique(population_studies$Country[population_studies$Continent=="SAM"]) !="Cuba"]

Overall_Proportion_Plot = function(population_people_use_fam,population_people_use_spor,population_studies_use_fam,population_studies_use_spor,title,fam_prop_slider){
  # A function to create publication style plots when using the overall proportion (rather than Familial or Sporadic)
  # Collpase files 
  population_studies_use_collapse_fam <- aggregate(population_studies_use_fam$Count,by=list(gene=population_studies_use_fam$Gene),FUN=sum)
  # Rename columns 
  colnames(population_studies_use_collapse_fam) <- c("gene","x_fam")
  # Calculate familial carrier carrier 
  population_studies_use_collapse_fam$Carrier_count_fam<-0 
  for (Gene in population_studies_use_collapse_fam$gene){
    population_studies_use_collapse_fam$Carrier_count_fam <- ifelse(population_studies_use_collapse_fam$gene==Gene,nrow(population_people_use_fam[population_people_use_fam$gene==Gene,]),population_studies_use_collapse_fam$Carrier_count_fam)
  }
  # Calculate binomial proportion, upper and lower bound for familial
  population_studies_use_collapse_fam$proportion_fam <-binom.confint(x=population_studies_use_collapse_fam$Carrier_count_fam,n=population_studies_use_collapse_fam$x_fam,method='wilson')$mean
  population_studies_use_collapse_fam$proportion.lower_fam <-binom.confint(x=population_studies_use_collapse_fam$Carrier_count_fam,n=population_studies_use_collapse_fam$x_fam,method='wilson')$lower
  population_studies_use_collapse_fam$proportion.upper_fam <-binom.confint(x=population_studies_use_collapse_fam$Carrier_count_fam,n=population_studies_use_collapse_fam$x_fam,method='wilson')$upper
  population_studies_use_collapse_spor <- aggregate(population_studies_use_spor$Count,by=list(gene=population_studies_use_spor$Gene),FUN=sum)
  colnames(population_studies_use_collapse_spor) <- c("gene","x_spor")
  # Create blank column
  population_studies_use_collapse_spor$Carrier_count_spor<-0
  for (Gene in population_studies_use_collapse_spor$gene){
    population_studies_use_collapse_spor$Carrier_count_spor <- ifelse(population_studies_use_collapse_spor$gene==Gene,nrow(population_people_use_spor[population_people_use_spor$gene==Gene,]),population_studies_use_collapse_spor$Carrier_count_spor)
  }
  # Calculate binomial proportion, upper and lower bound for sporadic
  population_studies_use_collapse_spor$proportion_spor <-binom.confint(x=population_studies_use_collapse_spor$Carrier_count_spor,n=population_studies_use_collapse_spor$x_spor,method='wilson')$mean
  population_studies_use_collapse_spor$proportion.lower_spor <-binom.confint(x=population_studies_use_collapse_spor$Carrier_count_spor,n=population_studies_use_collapse_spor$x_spor,method='wilson')$lower
  population_studies_use_collapse_spor$proportion.upper_spor <-binom.confint(x=population_studies_use_collapse_spor$Carrier_count_spor,n=population_studies_use_collapse_spor$x_spor,method='wilson')$upper
  # Merge familial and sporadic 
  population_studies_use_collapse <- merge(population_studies_use_collapse_fam,population_studies_use_collapse_spor,all.x=T,all.y=T)
  # Need to make some corrections for NA values after merging: set NA to 0 for proportion explained and lower bound and 100 for upper bound (i.e. completely unknown)
  population_studies_use_collapse$proportion_fam<-ifelse(is.na(population_studies_use_collapse$proportion_fam)==TRUE,0,population_studies_use_collapse$proportion_fam)
  population_studies_use_collapse$proportion.lower_fam<-ifelse(is.na(population_studies_use_collapse$proportion.lower_fam)==TRUE,0,population_studies_use_collapse$proportion.lower_fam)
  population_studies_use_collapse$proportion.lower_fam<-ifelse(population_studies_use_collapse$proportion.lower_fam<0,0,population_studies_use_collapse$proportion.lower_fam)
  population_studies_use_collapse$proportion.upper_fam<-ifelse(is.na(population_studies_use_collapse$proportion.upper_fam)==TRUE,100,population_studies_use_collapse$proportion.upper_fam)    
  population_studies_use_collapse$proportion_spor<-ifelse(is.na(population_studies_use_collapse$x_spor)==TRUE,0,population_studies_use_collapse$proportion_spor)
  population_studies_use_collapse$proportion.lower_spor<-ifelse(is.na(population_studies_use_collapse$proportion.lower_spor)==TRUE,0,population_studies_use_collapse$proportion.lower_spor)
  population_studies_use_collapse$proportion.lower_spor<-ifelse(population_studies_use_collapse$proportion.lower_spor<0,0,population_studies_use_collapse$proportion.lower_spor)
  population_studies_use_collapse$proportion.upper_spor<-ifelse(is.na(population_studies_use_collapse$proportion.upper_spor)==TRUE,100,population_studies_use_collapse$proportion.upper_spor)    
  # Weight the binomial proportions by familial proportions slider 
  population_studies_use_collapse$proportion <-(fam_prop_slider*population_studies_use_collapse$proportion_fam) + ((1-fam_prop_slider)*population_studies_use_collapse$proportion_spor)
  # Weight the lower CIs 
  population_studies_use_collapse$proportion.lower <-(fam_prop_slider*population_studies_use_collapse$proportion.lower_fam) + ((1-fam_prop_slider)*population_studies_use_collapse$proportion.lower_spor)
  # Weight the upper CIs
  population_studies_use_collapse$proportion.upper <-(fam_prop_slider*population_studies_use_collapse$proportion.upper_fam) + ((1-fam_prop_slider)*population_studies_use_collapse$proportion.upper_spor)
  population_studies_use_collapse$proportion.upper[population_studies_use_collapse$proportion.upper>1 & !is.na(population_studies_use_collapse$proportion.upper)] <- 1
  population_studies_use_collapse$proportion.lower[population_studies_use_collapse$proportion.lower<0 & !is.na(population_studies_use_collapse$proportion.lower)] <- 0
  population_studies_use_collapse<- population_studies_use_collapse[order(-population_studies_use_collapse$proportion),]
  population_studies_use_collapse$proportion<- 100*population_studies_use_collapse$proportion
  population_studies_use_collapse$proportion.lower<-100*population_studies_use_collapse$proportion.lower
  population_studies_use_collapse$proportion.upper<- 100*population_studies_use_collapse$proportion.upper
  # Create Plot 
  blank_plot(0,100,0,125)
  text(50,120,title,col=nicegrey,adj=0.5,cex=1.1,xpd=NA)
  scale_bar_plot_log(0,100)
  # position of first box 
  height=13
  # Need to cycle through each gene 
  for (gene in rev(population_studies_use_collapse$gene)) {
    colour<-vangogh_palette[match(gene,p_or_lp_genes)]
    mean_position=population_studies_use_collapse$proportion[population_studies_use_collapse$gene==gene]
    lower_position=population_studies_use_collapse$proportion.lower[population_studies_use_collapse$gene==gene]
    upper_position=population_studies_use_collapse$proportion.upper[population_studies_use_collapse$gene==gene]
    # if mean_position is zero set to 0.01 which will plot at zero 
    if (mean_position == 0 & !is.na(mean_position)){
      mean_position=0.01
    }
    # if lower_position is zero set to 0.01 which will plot at zero 
    if ((lower_position <= 0 | log_position(4,lower_position) < 0) & !is.na(lower_position)) {
      lower_position=0.01
    }
    # if upper_position is zero set to 0.01 which will plot at zero 
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
    roundedRect(left_pos+log_position(4,mean_position)-box_dimension/2.2,height-box_dimension/2,left_pos+log_position(4,mean_position)+box_dimension/2.2,height+box_dimension/2,border=NA,col=colour, bothsame=TRUE, aspcorrect=TRUE,corfactor=1)
    # Plot gene name on first pass
    text(left_pos-5,height,gene,col=nicegrey,adj=1,cex=0.8,font=3,xpd=NA)
    # Move up to next box position 
    height=height+5
  }
}

#####
# Create Plot
#####

pdf(paste(output_directory,"/Chapter5_Proportion_of_FTD_Cases_Explained_Results.pdf"),width=(11.7-0.8-0.8),height=(8.3-1.38-0.8))

filter="Ireland"
# Unlike other plots want to create a 3 across 2 down plot (not a 2x2)

par(mfrow=c(2,3), mai = c(0.7,1,0.2,0.1))

####
# Create 6 plots
####

###########################################################
# Plot 1 : Ireland FTD Pathogenic or LP Variants 
###########################################################

regions_pathogenicity_selection=2 #pathogenic or likely pathogenic variants 
regions_pathogenicity_selection1=list_of_pathogenic_variants[[as.numeric(regions_pathogenicity_selection)]]
# Create DF of mean proportions 
row_count<-length(unique(population_studies$Gene[population_studies$Gene %in% p_or_lp_genes & population_studies$Phenotype=="FTD"]))
mean_gene_dataframe<-data.frame(FTD_overall=rep(NA,row_count))
rownames(mean_gene_dataframe) <- unique(population_studies$Gene[population_studies$Gene %in% unique(population_studies$Gene[population_studies$Gene %in% p_or_lp_genes & population_studies$Phenotype=="FTD"])])
for(gene in rownames(mean_gene_dataframe)){
  total=sum(population_studies$Count[population_studies$Gene == gene & population_studies$Phenotype=="FTD" & !is.na(population_studies$Phenotype) & !is.na(population_studies$Gene) & population_studies$Country=="Ireland"])
  cases=nrow(population_people_dataset[!is.na(population_people_dataset$gene) & population_people_dataset$gene==gene & !is.na(population_people_dataset$HGVS) & population_people_dataset$HGVS %in% regions_pathogenicity_selection1 & !is.na(population_people_dataset$population_carriers_primary_phenotype) & population_people_dataset$population_carriers_primary_phenotype == "FTD" & (population_people_dataset$population_carriers_nationality=="Ireland"),])
  mean_gene_dataframe[gene,] <- 100*binom.confint(x=cases,n=total,method='wilson')$mean
}
left_pos <- 0 
blank_plot(0,100,0,125)
text(50,120,"Ireland",col=nicegrey,adj=0.5,cex=1.1,xpd=NA)
scale_bar_plot_log(0,100)
# Find order to plot genes y reordering dataframe 
# Add column to aid sorting 
mean_gene_dataframe$test <- NA 
mean_gene_dataframe<- mean_gene_dataframe[order(-mean_gene_dataframe$FTD_overall),]
# position of first box 
height=13
# Need to cycle through each gene 
for (gene in rev(rownames(mean_gene_dataframe))){
  colour<-vangogh_palette[match(gene,p_or_lp_genes)]
  total=sum(population_studies$Count[population_studies$Gene == gene & population_studies$Phenotype=="FTD" & !is.na(population_studies$Phenotype) & !is.na(population_studies$Gene) & population_studies$Country=="Ireland"])
  cases=nrow(population_people_dataset[!is.na(population_people_dataset$gene) & population_people_dataset$gene==gene & !is.na(population_people_dataset$HGVS) & population_people_dataset$HGVS %in% regions_pathogenicity_selection1 & !is.na(population_people_dataset$population_carriers_primary_phenotype) & population_people_dataset$population_carriers_primary_phenotype == "FTD" & (population_people_dataset$population_carriers_nationality=="Ireland"),])
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
  roundedRect(left_pos+log_position(4,mean_position)-box_dimension/2.2,height-box_dimension/2,left_pos+log_position(4,mean_position)+box_dimension/2.2,height+box_dimension/2,border=NA,col=colour, bothsame=TRUE, aspcorrect=TRUE,corfactor=1)
  # Plot gene name on first pass
  text(left_pos-5,height,gene,col=nicegrey,adj=1,cex=0.8,font=3,xpd=NA)
  # Move up to next box position 
  height=height+5
}


###########################################################
# Plot 2 : Europe FTD Pathogenic or LP Variants 
###########################################################

# This plot is complicated because we are separately calculating the familial proportion and sporadic proportion 
regions_pathogenicity_selection=2 #pathogenic or likely pathogenic 
regions_pathogenicity_selection1=list_of_pathogenic_variants[[as.numeric(regions_pathogenicity_selection)]]
regions_phenotype_selection="FTD"
# Filter population individuals: phenotype of interest; history of interest ; variant matches 
population_people_use_fam    <- population_people_dataset[population_people_dataset$gene %in% p_or_lp_genes & !is.na(population_people_dataset$gene) & population_people_dataset$population_carriers_primary_phenotype==regions_phenotype_selection & !is.na(population_people_dataset$population_carriers_primary_phenotype) & population_people_dataset$population_carriers_family_history=="Familial" & !is.na(population_people_dataset$population_carriers_family_history) & (population_people_dataset$HGVS %in% regions_pathogenicity_selection1),]
population_people_use_spor   <- population_people_dataset[population_people_dataset$gene %in% p_or_lp_genes & !is.na(population_people_dataset$gene) & population_people_dataset$population_carriers_primary_phenotype==regions_phenotype_selection & !is.na(population_people_dataset$population_carriers_primary_phenotype) & population_people_dataset$population_carriers_family_history=="Sporadic" & !is.na(population_people_dataset$population_carriers_family_history) & (population_people_dataset$HGVS %in% regions_pathogenicity_selection1),]
# Filter population studies: phenotype of interest; history of interest 
population_studies_use_fam   <- population_studies[population_studies$Gene %in% p_or_lp_genes & !is.na(population_studies$Gene) & population_studies$Phenotype==regions_phenotype_selection & !is.na(population_studies$Phenotype) & population_studies$History=="Familial" & !is.na(population_studies$History),]
population_studies_use_spor   <- population_studies[population_studies$Gene %in% p_or_lp_genes & !is.na(population_studies$Gene) & population_studies$Phenotype==regions_phenotype_selection & !is.na(population_studies$Phenotype) & population_studies$History=="Sporadic" & !is.na(population_studies$History),]
# Filter to Europe 
population_people_use_fam  <- population_people_use_fam[population_people_use_fam$population_carriers_continent=="EUR" & !is.na(population_people_use_fam$population_carriers_continent) & population_people_use_fam$population_carriers_nationality != "Ireland" & !is.na(population_people_use_fam$population_carriers_nationality),]
population_people_use_spor  <- population_people_use_spor[population_people_use_spor$population_carriers_continent=="EUR" & !is.na(population_people_use_spor$population_carriers_continent) & population_people_use_spor$population_carriers_nationality != "Ireland" & !is.na(population_people_use_spor$population_carriers_nationality),]
population_studies_use_fam <- population_studies_use_fam[population_studies_use_fam$Continent=="EUR" & !is.na(population_studies_use_fam$Continent) & population_studies_use_fam$Country != "Ireland" & !is.na(population_studies_use_fam$Country),]           
population_studies_use_spor <- population_studies_use_spor[population_studies_use_spor$Continent=="EUR" & !is.na(population_studies_use_spor$Continent) & population_studies_use_spor$Country != "Ireland" & !is.na(population_studies_use_spor$Country),]     
Overall_Proportion_Plot(population_people_use_fam,population_people_use_spor,population_studies_use_fam,population_studies_use_spor,"Europe",0.4)
text(x=log_position(4,1),y=130,"FTD: Pathogenic or Likely Pathogenic Variants",cex=1.2,col=nicegrey,adj=0.5,xpd=NA)
###########################################################
# Plot 3 : Global FTD Pathogenic or LP Variants 
###########################################################

# This plot is complicated because we are separately calculating the familial proportion and sporadic proportion 
regions_pathogenicity_selection=2 #pathogenic or likely pathogenic 
regions_pathogenicity_selection1=list_of_pathogenic_variants[[as.numeric(regions_pathogenicity_selection)]]
regions_phenotype_selection="FTD"
# Filter population individuals: phenotype of interest; history of interest ; variant matches 
population_people_use_fam    <- population_people_dataset[population_people_dataset$gene %in% p_or_lp_genes & !is.na(population_people_dataset$gene) & population_people_dataset$population_carriers_primary_phenotype==regions_phenotype_selection & !is.na(population_people_dataset$population_carriers_primary_phenotype) & population_people_dataset$population_carriers_family_history=="Familial" & !is.na(population_people_dataset$population_carriers_family_history) & (population_people_dataset$HGVS %in% regions_pathogenicity_selection1),]
population_people_use_spor   <- population_people_dataset[population_people_dataset$gene %in% p_or_lp_genes & !is.na(population_people_dataset$gene) & population_people_dataset$population_carriers_primary_phenotype==regions_phenotype_selection & !is.na(population_people_dataset$population_carriers_primary_phenotype) & population_people_dataset$population_carriers_family_history=="Sporadic" & !is.na(population_people_dataset$population_carriers_family_history) & (population_people_dataset$HGVS %in% regions_pathogenicity_selection1),]
# Filter population studies: phenotype of interest; history of interest 
population_studies_use_fam   <- population_studies[population_studies$Gene %in% p_or_lp_genes & !is.na(population_studies$Gene) & population_studies$Phenotype==regions_phenotype_selection & !is.na(population_studies$Phenotype) & population_studies$History=="Familial" & !is.na(population_studies$History),]
population_studies_use_spor   <- population_studies[population_studies$Gene %in% p_or_lp_genes & !is.na(population_studies$Gene) & population_studies$Phenotype==regions_phenotype_selection & !is.na(population_studies$Phenotype) & population_studies$History=="Sporadic" & !is.na(population_studies$History),]
# Filter to not Ireland 
population_people_use_fam  <- population_people_use_fam[population_people_use_fam$population_carriers_nationality != "Ireland" & !is.na(population_people_use_fam$population_carriers_nationality),]
population_people_use_spor  <- population_people_use_spor[population_people_use_spor$population_carriers_nationality != "Ireland" & !is.na(population_people_use_spor$population_carriers_nationality),]
population_studies_use_fam <- population_studies_use_fam[population_studies_use_fam$Country != "Ireland" & !is.na(population_studies_use_fam$Country),]           
population_studies_use_spor <- population_studies_use_spor[population_studies_use_spor$Country != "Ireland" & !is.na(population_studies_use_spor$Country),]     
Overall_Proportion_Plot(population_people_use_fam,population_people_use_spor,population_studies_use_fam,population_studies_use_spor,"Global",0.4)

###########################################################
# Plot 4 : Ireland FTD Reported Variants in P or LP Genes 
###########################################################

regions_pathogenicity_selection=3 #pathogenic or likely pathogenic variants 
regions_pathogenicity_selection1=list_of_pathogenic_variants[[as.numeric(regions_pathogenicity_selection)]]
# Create DF of mean proportions 
row_count<-length(unique(population_studies$Gene[population_studies$Gene %in% p_or_lp_genes & population_studies$Phenotype=="FTD"]))
mean_gene_dataframe<-data.frame(FTD_overall=rep(NA,row_count))
rownames(mean_gene_dataframe) <- unique(population_studies$Gene[population_studies$Gene %in% unique(population_studies$Gene[population_studies$Gene %in% p_or_lp_genes & population_studies$Phenotype=="FTD"])])
for(gene in rownames(mean_gene_dataframe)){
  total=sum(population_studies$Count[population_studies$Gene == gene & population_studies$Phenotype=="FTD" & !is.na(population_studies$Phenotype) & !is.na(population_studies$Gene) & population_studies$Country=="Ireland"])
  cases=nrow(population_people_dataset[!is.na(population_people_dataset$gene) & population_people_dataset$gene==gene & !is.na(population_people_dataset$HGVS) & population_people_dataset$HGVS %in% regions_pathogenicity_selection1 & !is.na(population_people_dataset$population_carriers_primary_phenotype) & population_people_dataset$population_carriers_primary_phenotype == "FTD" & (population_people_dataset$population_carriers_nationality=="Ireland"),])
  mean_gene_dataframe[gene,] <- 100*binom.confint(x=cases,n=total,method='wilson')$mean
}
left_pos <- 0 
blank_plot(0,100,0,125)
text(50,120,"Ireland",col=nicegrey,adj=0.5,cex=1.1,xpd=NA)
scale_bar_plot_log(0,100)
# Find order to plot genes y reordering dataframe 
# Add column to aid sorting 
mean_gene_dataframe$test <- NA 
mean_gene_dataframe<- mean_gene_dataframe[order(-mean_gene_dataframe$FTD_overall),]
# position of first box 
height=13
# Need to cycle through each gene 
for (gene in rev(rownames(mean_gene_dataframe))){
  colour<-vangogh_palette[match(gene,p_or_lp_genes)]
  total=sum(population_studies$Count[population_studies$Gene == gene & population_studies$Phenotype=="FTD" & !is.na(population_studies$Phenotype) & !is.na(population_studies$Gene) & population_studies$Country=="Ireland"])
  cases=nrow(population_people_dataset[!is.na(population_people_dataset$gene) & population_people_dataset$gene==gene & !is.na(population_people_dataset$HGVS) & population_people_dataset$HGVS %in% regions_pathogenicity_selection1 & !is.na(population_people_dataset$population_carriers_primary_phenotype) & population_people_dataset$population_carriers_primary_phenotype == "FTD" & (population_people_dataset$population_carriers_nationality=="Ireland"),])
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
  roundedRect(left_pos+log_position(4,mean_position)-box_dimension/2.2,height-box_dimension/2,left_pos+log_position(4,mean_position)+box_dimension/2.2,height+box_dimension/2,border=NA,col=colour, bothsame=TRUE, aspcorrect=TRUE,corfactor=1)
  # Plot gene name on first pass
  text(left_pos-5,height,gene,col=nicegrey,adj=1,cex=0.8,font=3,xpd=NA)
  # Move up to next box position 
  height=height+5
}

###########################################################
# Plot 5 : Europe FTD Reported Variants in P or LP Genes  
###########################################################

# This plot is complicated because we are separately calculating the familial proportion and sporadic proportion 
regions_pathogenicity_selection=3 #pathogenic or likely pathogenic 
regions_pathogenicity_selection1=list_of_pathogenic_variants[[as.numeric(regions_pathogenicity_selection)]]
regions_phenotype_selection="FTD"
# Filter population individuals: phenotype of interest; history of interest ; variant matches 
population_people_use_fam    <- population_people_dataset[population_people_dataset$gene %in% p_or_lp_genes & !is.na(population_people_dataset$gene) & population_people_dataset$population_carriers_primary_phenotype==regions_phenotype_selection & !is.na(population_people_dataset$population_carriers_primary_phenotype) & population_people_dataset$population_carriers_family_history=="Familial" & !is.na(population_people_dataset$population_carriers_family_history) & (population_people_dataset$HGVS %in% regions_pathogenicity_selection1),]
population_people_use_spor   <- population_people_dataset[population_people_dataset$gene %in% p_or_lp_genes & !is.na(population_people_dataset$gene) & population_people_dataset$population_carriers_primary_phenotype==regions_phenotype_selection & !is.na(population_people_dataset$population_carriers_primary_phenotype) & population_people_dataset$population_carriers_family_history=="Sporadic" & !is.na(population_people_dataset$population_carriers_family_history) & (population_people_dataset$HGVS %in% regions_pathogenicity_selection1),]
# Filter population studies: phenotype of interest; history of interest 
population_studies_use_fam   <- population_studies[population_studies$Gene %in% p_or_lp_genes & !is.na(population_studies$Gene) & population_studies$Phenotype==regions_phenotype_selection & !is.na(population_studies$Phenotype) & population_studies$History=="Familial" & !is.na(population_studies$History),]
population_studies_use_spor   <- population_studies[population_studies$Gene %in% p_or_lp_genes & !is.na(population_studies$Gene) & population_studies$Phenotype==regions_phenotype_selection & !is.na(population_studies$Phenotype) & population_studies$History=="Sporadic" & !is.na(population_studies$History),]
# Filter to Europe 
population_people_use_fam  <- population_people_use_fam[population_people_use_fam$population_carriers_continent=="EUR" & !is.na(population_people_use_fam$population_carriers_continent) & population_people_use_fam$population_carriers_nationality != "Ireland" & !is.na(population_people_use_fam$population_carriers_nationality),]
population_people_use_spor  <- population_people_use_spor[population_people_use_spor$population_carriers_continent=="EUR" & !is.na(population_people_use_spor$population_carriers_continent) & population_people_use_spor$population_carriers_nationality != "Ireland" & !is.na(population_people_use_spor$population_carriers_nationality),]
population_studies_use_fam <- population_studies_use_fam[population_studies_use_fam$Continent=="EUR" & !is.na(population_studies_use_fam$Continent) & population_studies_use_fam$Country != "Ireland" & !is.na(population_studies_use_fam$Country),]           
population_studies_use_spor <- population_studies_use_spor[population_studies_use_spor$Continent=="EUR" & !is.na(population_studies_use_spor$Continent) & population_studies_use_spor$Country != "Ireland" & !is.na(population_studies_use_spor$Country),]     
Overall_Proportion_Plot(population_people_use_fam,population_people_use_spor,population_studies_use_fam,population_studies_use_spor,"Europe",0.4)
text(x=log_position(4,1),y=130,"FTD: Reported Variants in Pathogenic or Likely Pathogenic Genes",cex=1.2,col=nicegrey,adj=0.5,xpd=NA)

###########################################################
# Plot 6 : Global FTD Reported Variants in P or LP Genes
###########################################################

# This plot is complicated because we are separately calculating the familial proportion and sporadic proportion 
regions_pathogenicity_selection=3 #pathogenic or likely pathogenic 
regions_pathogenicity_selection1=list_of_pathogenic_variants[[as.numeric(regions_pathogenicity_selection)]]
regions_phenotype_selection="FTD"
# Filter population individuals: phenotype of interest; history of interest ; variant matches 
population_people_use_fam    <- population_people_dataset[population_people_dataset$gene %in% p_or_lp_genes & !is.na(population_people_dataset$gene) & population_people_dataset$population_carriers_primary_phenotype==regions_phenotype_selection & !is.na(population_people_dataset$population_carriers_primary_phenotype) & population_people_dataset$population_carriers_family_history=="Familial" & !is.na(population_people_dataset$population_carriers_family_history) & (population_people_dataset$HGVS %in% regions_pathogenicity_selection1),]
population_people_use_spor   <- population_people_dataset[population_people_dataset$gene %in% p_or_lp_genes & !is.na(population_people_dataset$gene) & population_people_dataset$population_carriers_primary_phenotype==regions_phenotype_selection & !is.na(population_people_dataset$population_carriers_primary_phenotype) & population_people_dataset$population_carriers_family_history=="Sporadic" & !is.na(population_people_dataset$population_carriers_family_history) & (population_people_dataset$HGVS %in% regions_pathogenicity_selection1),]
# Filter population studies: phenotype of interest; history of interest 
population_studies_use_fam   <- population_studies[population_studies$Gene %in% p_or_lp_genes & !is.na(population_studies$Gene) & population_studies$Phenotype==regions_phenotype_selection & !is.na(population_studies$Phenotype) & population_studies$History=="Familial" & !is.na(population_studies$History),]
population_studies_use_spor   <- population_studies[population_studies$Gene %in% p_or_lp_genes & !is.na(population_studies$Gene) & population_studies$Phenotype==regions_phenotype_selection & !is.na(population_studies$Phenotype) & population_studies$History=="Sporadic" & !is.na(population_studies$History),]
# Filter to not Ireland 
population_people_use_fam  <- population_people_use_fam[population_people_use_fam$population_carriers_nationality != "Ireland" & !is.na(population_people_use_fam$population_carriers_nationality),]
population_people_use_spor  <- population_people_use_spor[population_people_use_spor$population_carriers_nationality != "Ireland" & !is.na(population_people_use_spor$population_carriers_nationality),]
population_studies_use_fam <- population_studies_use_fam[population_studies_use_fam$Country != "Ireland" & !is.na(population_studies_use_fam$Country),]           
population_studies_use_spor <- population_studies_use_spor[population_studies_use_spor$Country != "Ireland" & !is.na(population_studies_use_spor$Country),]     
Overall_Proportion_Plot(population_people_use_fam,population_people_use_spor,population_studies_use_fam,population_studies_use_spor,"Global",0.4)


dev.off()


################################################################################
################################################################################

# ALS - Proportion Explained - Proportions explained in Ireland and Europe 

################################################################################
################################################################################


#####
# Parse Data 
#####


population_people_dataset <<- data.frame(population_people_dataset_function(lit_review_pen))
# The current existing Cuba data in journALS is a subset of this data 
population_people_dataset <<- population_people_dataset[population_people_dataset$population_carriers_nationality!="Ireland",]
population_studies <<- population_studies[population_studies$Country != "Ireland",]
# Get gene for Cuban Data 
als_population_people_dataset$gene <- gsub(":.*","",als_population_people_dataset$HGVS)
# Merge Back in new Cuban Data 
population_studies <<- rbind(population_studies,als_population_studies)
als_population_people_dataset <- subset(als_population_people_dataset,select=colnames(population_people_dataset))
population_people_dataset <<- rbind(population_people_dataset,als_population_people_dataset)
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
EUR_countries <- unique(population_studies$Country[population_studies$Continent=="EUR" & population_studies$Country != "Ireland"])
#####
# Create Plot
#####
pdf(paste(output_directory,"/Chapter5_Proportion_of_ALS_Cases_Explained_Results.pdf",height=11.7,width=(8.3),onefile=T)

for (filter in list("Ireland",EUR_countries)){


  ####
  # OVERVIEW PLOT 
  ###
  par(mfrow=c(4,1), mai = c(0.1,0.1,0.1,0.1))

  ####
  # Will Create 3 2x2 plots - Pathogenic, Pathogenic and Likely Pathogenic, Reported Variants
  ###
  for (regions_pathogenicity_selection in c(2,3)){ # 1 is P variants, 2 is P or LP variants, 3 is all variants 
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
          if (regions_pathogenicity_selection==2){
            if (filter == "Ireland"){
              text((left.1+right.2)/2,120,"Ireland: Pathogenic or Likely Pathogenic Variants",col=nicegrey,adj=0.5,cex=1.1,xpd=NA)
              text(left.1-(0.12*(right.2-left.1)),120,col=nicegrey,"A",adj=0,cex=1.2,xpd=NA,font=2)
            }
            if (filter == EUR_countries){
              text((left.1+right.2)/2,120,"Europe: Pathogenic or Likely Pathogenic Variants",col=nicegrey,adj=0.5,cex=1.1,xpd=NA)
              text(left.1-(0.12*(right.2-left.1)),120,"B",col=nicegrey,adj=0,cex=1.2,xpd=NA,font=2)
            }
          }
          if (regions_pathogenicity_selection==3){
            if (filter == "Ireland"){
              text((left.1+right.2)/2,120,"Ireland: Reported Variants in Genes with Pathogenic or Likely Pathogenic Variants",col=nicegrey,adj=0.5,cex=1.1,xpd=NA)
            }
            if (filter == EUR_countries){
              text((left.1+right.2)/2,120,"Europe: Reported Variants in Genes with Pathogenic or Likely Pathogenic Variants",col=nicegrey,adj=0.5,cex=1.1,xpd=NA)
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
              text(left.1-30,70,regions_history_selection,col=nicegrey,adj=0.5,cex=0.8,srt=90,xpd=NA)
            }
            else {
              text(left.2-30,70,regions_history_selection,col=nicegrey,adj=0.5,cex=0.8,srt=90,xpd=NA)
            }
          }
      }
  }
  dev.off()

  
