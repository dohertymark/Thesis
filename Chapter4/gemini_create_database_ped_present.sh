#!/bin/bash
###########################################################################################################
###########################################################################################################

# A script to create a gemini database from a vcf and ped file 
# Requires a user to input locations of reference files 
# Script should be placed in $PATH and run in the directory you want files to be created
# email: dohertm7@tcd.ie

###########################################################################################################
###########################################################################################################

echo "FILE SHOULD BE A .vcf , IF FILE IS .vcf.gz  , ENTER CRTL+C FOLLOWED BY gunzip file, IF FILE IS .vcf HIT ENTER"
read
echo "PLEASE ENTER FILE NAME"
read FILE
echo "PLEASE ENTER REFERENCE NAME AND PATH"
read REFERENCE
echo "ENTER NAME FOR FINAL DATABASE e.b. database.db"
read DATABASE
echo "ENTER FULL NAME PF .ped FILE"
read PED
echo "HOW MANY CORES"
read NUM_CORES

###########################################################################################################
###########################################################################################################

vt decompose -s $FILE | vt normalize -r $REFERENCE - > POST_VT_${FILE}
java -Xmx16G -jar snpEff.jar -c snpEff.config GRCh37.75 POST_VT_${FILE} > POST_SNPEFF_POST_VT_${FILE}
gemini load -v POST_SNPEFF_POST_VT_${FILE} -p ${PED} -t snpEff ${DATABASE} --core ${NUM_CORES}

###########################################################################################################
###########################################################################################################

# Remove unnecessary files 
rm POST_VT_${FILE}*
rm POST_SNPEFF_POST_VT_${FILE}*
rm snpEff*
###########################################################################################################
###########################################################################################################

exit 0
