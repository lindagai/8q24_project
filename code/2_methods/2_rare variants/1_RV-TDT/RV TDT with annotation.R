#######NOTE: Before running this code, run functions.R to define functions used here.

#Description:
#Filters 8q24 region to SNPs with CADD>10 and runs RV TDT

#1. Read in phased vcf file

#2. Subset vcf files depending on annotation information (CADD, eigen, gwava)
#i.e., subset to list of positions in AnnoVar that satisfy requirements
#TODO:
#i. Obtain additional scores for 8q24 region: e.g. SIFT, PolyPhen and plot
#       SIFT, PolyPhen scores not available for this regions
#ii. Use existing VariantAnnotation functions e.g. locateVariant and predictCoding to filter

#3. Get vcf file into RV TDT format
#The phen file does not need to be recreated each time.
#i. trio ped format -> phen file
#TODO: change this so it converts raw plink -> phen file directly

#The tped and map files need to be recreated each time you run this
#ii. vcf -> tped file
#iii. vcf -> map file

#4. Create n-SNP windows
#(max. RV TDT can handle at a time is ~1000)

#5. Run RV TDT

#6. TODO: Add visualization

#NOTE: ~400 SNPS max seems to be what RV-TDT can handle.

#Author: Linda Gai
#Last update: 3/23/19

#TODO: Put each section into its own file

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

#IMPORTANT NOTE: RV-TDT is unable to handle absolute directories -- u must give relative directories

# for (i in 1:n.windows){
#         rv.tdt.dir<-"./rv-tdt/rvTDT "
#         work.dir <- "cd /users/lgai/"
#         
#         data.dir <- "./8q24_project/data/processed_data/RV-TDT/window25_no_overlap/"
#         data.file <- paste0("rarevar.1000G.euro.annotation.window-size=25_snps_",i)
#         results.dir <- paste0("./8q24_project/data/processed_data/RV-TDT/window25_no_overlap/window-size=25_snps_",i)
#         g.file <- paste0(data.file,".tped")
#         p.file <- "./8q24_project/data/processed_data/RV-TDT/gmkf_euro_completetrios_ids_match_vcf_mend_errors_removed.phen"
#         
#         m.file <- paste0(data.file,".map")
#         
#         #This works
#         command<-paste0(rv.tdt.dir, results.dir, 
#                         " -G ", data.dir, g.file,
#                         " -P ", p.file,
#                         " -M ", data.dir, m.file,
#                         " --adapt 100 --alpha 0.00001 --permut 1000",
#                         " -u 0.01"
#         )
#         
#         command
#         
#         system(work.dir)
#         system(command)     
# }

#################################################################################