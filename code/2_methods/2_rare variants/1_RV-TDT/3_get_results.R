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

########################### 5. Run RV TDT ##############################################

n.snps<-368
window.size<-25
n.windows<-ceiling(n.snps/window.size)
n.windows
25*n.windows

#Whole 8q24 region, no overlap, 25 SNPs
#TODO: add code to change the window directory automatically

# i<-1
# n.window.sizes<-1
# window.size=25
# overlap=0
# window.type<-"snps"

filepath.tped<-"/users/lgai/8q24_project/data/processed_data/RV-TDT/8q24.cleaned.phased.rarevar.monomorphs.removed.1000G.euro.annotation.tped"
filepath.map<-"/users/lgai/8q24_project/data/processed_data/RV-TDT/8q24.cleaned.phased.rarevar.monomorphs.removed.1000G.euro.annotation.map"

#results.dir<-"/users/lgai/8q24_project/data/processed_data/RV-TDT/"
results.dir<-"./8q24_project/data/processed_data/RV-TDT/"
filepath.phen <- "./8q24_project/data/processed_data/RV-TDT/gmkf_euro_completetrios_ids_match_vcf_mend_errors_removed.phen"

get.rvtdt.results(window.size,overlap,filepath.phen,filepath.tped,filepath.map,results.dir)