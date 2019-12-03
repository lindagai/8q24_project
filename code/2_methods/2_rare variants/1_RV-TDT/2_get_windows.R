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

########################### 4. Break files into windows ##############################################

#TODO: must add system command for creating directory to put these files in

#25 SNP windows, whole 8q24 region, no overlap
# window.size<-25
# overlap<-0
window.size<-25
overlap<-24
window.type<-"snps"

filepath.tped<-"/users/lgai/8q24_project/data/processed_data/RV-TDT/8q24.cleaned.phased.rarevar.monomorphs.removed.1000G.euro.annotation.tped"
filepath.map<-"/users/lgai/8q24_project/data/processed_data/RV-TDT/8q24.cleaned.phased.rarevar.monomorphs.removed.1000G.euro.annotation.map"
results.dir<-"/users/lgai/8q24_project/data/processed_data/RV-TDT/"

get.rv_tdt.files.in.windows(window.size,overlap,filepath.tped,filepath.map,results.dir)