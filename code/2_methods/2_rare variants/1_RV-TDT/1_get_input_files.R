################################################################################
ssh -Y lgai@jhpce01.jhsph.edu

#Allocate memory
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G
module load R
R

library(VariantAnnotation)
library(dplyr)

################################################################################

#Before running the commands in this file, run functions.R.

######################## 1. Get tped ############################################

#CADD filtered version
#filepath.tped<-"~/cadd_10/cadd_10.tped"
#get.tped(vcf.geno,filepath.tped)

#Checks
#tped<-read.table(filepath.tped)
#tped[1:5,1:10]
#dim(tped) #810 1963

filepath.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/rare_var_only/8q24.cleaned.phased.filtered.annotation.rarevar.monomorphs.removed.recode.vcf"

vcf<-readVcf(filepath.vcf)
vcf.geno<-geno(vcf)$GT

filepath.tped<-"/users/lgai/8q24_project/data/processed_data/RV-TDT/8q24.cleaned.phased.rarevar.monomorphs.removed.1000G.euro.annotation.tped"

get.tped(vcf.geno,filepath.tped)

#Checks
tped<-read.table(filepath.tped)
tped[1:5,1:10]
dim(tped) #368 1963

############################## 2. Get map file ###################################################

# filepath_tped<-"~/cadd_10/cadd_10.tped"
# filepath_rvtdt_map<-"~/cadd_10/cadd_10.map"
# 
# get.map(filepath_tped,filepath_rvtdt_map)
# 
# map<-read.table(filepath_rvtdt_map)
# #head(map)
# #dim(map)
# #1283    3

filepath_tped<-"/users/lgai/8q24_project/data/processed_data/RV-TDT/8q24.cleaned.phased.rarevar.monomorphs.removed.1000G.euro.annotation.tped"
filepath_rvtdt_map<-"/users/lgai/8q24_project/data/processed_data/RV-TDT/8q24.cleaned.phased.rarevar.monomorphs.removed.1000G.euro.annotation.map"

get.map(filepath_tped,filepath_rvtdt_map)

map<-read.table(filepath_rvtdt_map)
head(map)
dim(map)
#368   3

############################ 3. Get phen file #############################################

#TRIO columns
#famid                pid              fatid              motid    sex  affected

filepath.trio.ped<-"/users/lgai/8q24_project/data/processed_data/gmkf_euro_completetrios_ids_match_vcf_mend_errors_removed.txt"
filepath.phen<-"/users/lgai/8q24_project/data/processed_data/RV-TDT/gmkf_euro_completetrios_ids_match_vcf_mend_errors_removed.phen"

ped<-read.table(filepath.trio.ped,header=TRUE)
head(ped)
dim(ped)

trio.ped.to.phen(filepath.trio.ped,vcf.geno,filepath.phen)