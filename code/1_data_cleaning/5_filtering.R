#Sign into cluster
ssh -X lgai@jhpce01.jhsph.edu 

#Allocate memory
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G
module load R
R

########################

library(dplyr)
library(trio)
library(VariantAnnotation)

###########################################


get.subsetted.vcf<-function(filepath.annovar,filepath.vcf,filepath.filtered.vcf,hg.assembly){
  filepath.annovar
  annovar<-read.table(filepath.annovar,sep="\t",header=TRUE, quote ="")
  pos<-annovar$StartPosition
    
  vcf <- readVcf(filepath.vcf,hg.assembly)
  geno(vcf)$GT<-apply(geno(vcf)$GT,2,function(x)gsub('\\|','/',x))
  geno(vcf)$GT<-apply(geno(vcf)$GT,2,function(x)gsub('1/0','0/1',x))
  
  vcf<-vcf[which(start(rowRanges(vcf)) %in% pos),]
  
  #TODO: remove the monomorphic SNPs from the vcf #########
  #Filter to MAF <0.01
  geno(vcf)$GT
  
  writeVcf(vcf,filepath.filtered.vcf)
}

get.geno<-function(filepath.ped,filepath.filtered.vcf,filepath.geno,hg.assembly){
  filepath.ped
  filepath.filtered.vcf
  vcf <- readVcf(filepath.filtered.vcf,hg.assembly)
  ped<-read.table(filepath.ped,header=TRUE)
  geno <- vcf2geno(vcf,ped)
  saveRDS(geno, filepath.geno) 
}

###########################################



################# Choose filtering #########################

#TODO: Add in more functions

#A. Subset Annovar report to find positions with allele frequency information for rvTDT

filepath.sm.annovar<-"/users/lgai/8q24_project/data/processed_data/annotation/8q24_annovar_report_useful_variables.txt"

annovar<-read.table(filepath.sm.annovar,sep="\t",header=TRUE, quote ="")

annovar.filt<-annovar %>%
  filter(!is.na(Afalt_1000g2015aug_eur)) %>%
  filter(TotalDepth>20)

dim(annovar.filt)

#B. Subset ANNOVAR based on functional annotation/TDT peak location/both.
annovar.filt.annotation <- annovar.filt %>%
  filter(CADDgt20 > 10 | WG_GWAVA_score > 0.4 | WG_EIGEN_score > 4)
filepath.annovar.filt.annotation<-"/users/lgai/8q24_project/data/processed_data/annotation/filtered_annovar/8q24_annovar_report_functional_annnotation_and_1000G_euro_af.txt"
write.table(annovar.filt.annotation, filepath.annovar.filt.annotation, sep="\t",row.names = FALSE,quote = FALSE)
#This one is actually exactly the same as the one you made before--you already filtered to 1000G
  
annovar.filt.filepaths<-c(filepath.annovar.filt.annotation)

#C. Filter vcf to these positions

########## Cleaned phased vcf#########
#TODO: rewrite this in R 
#TODO: Change this so it doesn't have to be unzipped

#Remove monomorphic SNPs

vcftools --vcf /users/lgai/8q24_project/data/processed_data/vcfs/8q24.cleaned.phased.vcf --max-maf 0.01 --non-ref-ac-any 1 --recode --out /users/lgai/8q24_project/data/processed_data/vcfs/rare_var_only/8q24.cleaned.phased.0.01RareVarOnly.NoMonomorph

filepath.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/rare_var_only/8q24.cleaned.phased.0.01RareVarOnly.NoMonomorph.recode.vcf"


#####

#OK, if you already have the filtered vcf, you can just filter to remove monomorphic SNPs/MAF>0.01
#TODO: fix this for R later

vcftools --vcf /users/lgai/8q24_project/data/processed_data/vcfs/filtered_by_annotation/8q24.cleaned.phased.filtered.annotation.vcf --max-maf 0.01 --non-ref-ac-any 1 --recode --out /users/lgai/8q24_project/data/processed_data/vcfs/rare_var_only/8q24.cleaned.phased.filtered.annotation.rarevar.monomorphs.removed

################################### USE THIS ONE ###################################################
filepath.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/rare_var_only/8q24.cleaned.phased.filtered.annotation.rarevar.monomorphs.removed.recode.vcf"


#################

#Filtered vcf

#   filepath.annovar<-annovar.filt.filepaths[i]
#   filepath.filtered.vcf<-filtered.vcf.filepaths[i]
#   filepath.filtered.geno<-filtered.geno.filepaths[i]
#   
#   #get.subsetted.vcf(filepath.annovar,filepath.vcf,filepath.filtered.vcf,hg.assembly)
#   get.geno(filepath.ped,filepath.filtered.vcf,filepath.filtered.geno,hg.assembly)

######################################## SCRATCH #############################################


# 368 SNPs after filtering with annotation information, monomorphs removed
filepath.filtered.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/rare_var_only/8q24.cleaned.phased.filtered.annotation.rarevar.monomorphs.removed.recode.vcf"
hg.assembly<-"hg19"
filepath.ped<-"/users/lgai/8q24_project/data/processed_data/gmkf_euro_completetrios_ids_match_vcf_mend_errors_removed.txt"
filepath.filtered.geno <- "/users/lgai/8q24_project/data/processed_data/geno_matrix/filtered/8q24.cleaned.phased.filtered.annotation.rarevar.monomorphs.removed.rds"

#get.subsetted.vcf(filepath.annovar,filepath.vcf,filepath.filtered.vcf,hg.assembly)
get.geno(filepath.ped,filepath.filtered.vcf,filepath.filtered.geno,hg.assembly)

# #C. Convert filtered, phased vcf to geno matrix
# filepath.ped<-"/users/lgai/8q24_project/data/processed_data/gmkf_euro_completetrios_ids_match_vcf_mend_errors_removed.txt"
# 
# filepath.geno.annotation <- "/users/lgai/8q24_project/data/processed_data/geno_matrix/filtered/geno_phased_annotation.rds"
# filtered.geno.filepaths<-c(filepath.geno.annotation)
# filtered.geno.filepaths
# 
# #Obtain the filtered vcfs and geno matrices
# #Filtered vcf
# filepath.filtered.vcf.annotation<-"/users/lgai/8q24_project/data/processed_data/vcfs/filtered_by_annotation/8q24.cleaned.phased.filtered.annotation.1000G.euro.af.vcf"
# hg.assembly<-"hg19"
# 
# filtered.vcf.filepaths<-c(filepath.filtered.vcf.annotation)
# filtered.vcf.filepaths
# 
# for (i in 1:1){
#   filepath.annovar<-annovar.filt.filepaths[i]
#   filepath.filtered.vcf<-filtered.vcf.filepaths[i]
#   filepath.filtered.geno<-filtered.geno.filepaths[i]
#   
#   #get.subsetted.vcf(filepath.annovar,filepath.vcf,filepath.filtered.vcf,hg.assembly)
#   get.geno(filepath.ped,filepath.filtered.vcf,filepath.filtered.geno,hg.assembly)
# }