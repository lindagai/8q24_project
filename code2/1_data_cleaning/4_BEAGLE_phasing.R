#Description: 
#Sections

#Author: Linda Gai
#Last update: 12/12/18

#TO DO: Clean up code to get rid of cruft

##########################################################

#Sign into cluster
ssh -X lgai@jhpce01.jhsph.edu 

#Allocate memory
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G
module load R
R

##########################################################

################ 1. Download BEAGLEv4.0 ##############

#Un-comment this out if you don't already have BEAGLE4.0
#NOTE: Do not use BEAGLE4.1 or 5.0: these do not use pedigree information.

#Download vcf2beagle
# filepath.beagle4<-"/users/lgai/beagle.r1399.jar"
# dl.beagle4<-paste0("wget -O ",filepath.beagle4," https://faculty.washington.edu/browning/beagle/beagle.r1399.jar")
# system(dl.beagle4)

#Phase vcf
filepath.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/8q24.cleaned.vcf"
filepath.ped<-"/users/lgai/8q24_project/data/processed_data/gmkf_euro_completetrios_ids_match_vcf_mend_errors_removed.txt"
filepath.phased.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/8q24.cleaned.phased"


phase.command<-paste0("java -Xmx10000m -jar ", filepath.beagle4,
                      " gt=",filepath.vcf,
                      " ped=",filepath.ped,
                      " out=",filepath.phased.vcf)
phase.command

system(phase.command)

# TODO: Check to make sure missing calls are dealt with appropriately... is it '.' or './.'?

#NOTE: BEAGLE only needs the family structure (1st 4 columns)
# So you don't need to recode the affected/unaffected status, etc.

##########################################################

#Old files

# filepath.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/8q24.recode.15.mend.error.individuals.removed.recode.vcf"
# filepath.ped<-"/users/lgai/8q24_project/data/processed_data/gmkf_euro_completetrios_ids_match_vcf_mend_errors_removed.txt"
# filepath.phased.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/8q24.recode.15.mend.error.individuals.removed.recode.BEAGLE"



#java -Xmx10000m -jar beagle.r1399.jar gt=/users/lgai/8q24_project/data/processed_data/vcfs/8q24.recode.15.mend.error.individuals.removed.recode.vcf ped=/users/lgai/8q24_project/data/processed_data/gmkf_euro_completetrios_ids_match_vcf_mend_errors_removed.txt out=/users/lgai/8q24_project/data/processed_data/vcfs/8q24.recode.15.mend.error.individuals.removed.recode.BEAGLE

#java -Xmx10000m -jar beagle.r1399.jar gt=/users/lgai/8q24_project/data/processed_data/vcfs/8q24.recode.15.mend.error.individuals.removed.recode.vcf ped=/users/lgai/8q24_project/data/processed_data/gmkf_euro_completetrios_ids_match_vcf_mend_errors_removed.txt out=/users/lgai/8q24_project/data/processed_data/vcfs/8q24.recode.15.mend.error.individuals.removed.recode.BEAGLE
