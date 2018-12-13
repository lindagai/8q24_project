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

################# Filter geno matrix #########################

#A. Subset Annovar report to find positions with allele frequency information

filepath.sm.annovar<-"/users/lgai/8q24_project/data/processed_data/annotation/8q24_annovar_report_useful_variables.txt"

annovar<-read.table(filepath.sm.annovar,sep="\t",header=TRUE, quote ="")

#TODO: Choose filters
# filterrules:

# Window size
# Window overlap size

window.size<-100

head(annovar)
dim(annovar)
annovar.filt<-annovar %>%
  filter(!is.na(Afalt_1000g2015aug_eur)) %>%
  filter(TotalDepth>20)

#B. Subset ANNOVAR based on functional annotation/TDT peak location/both/divide into windows.

#TODO: This could be coded more coded more efficiently
annovar.filt.peak <- annovar.filt %>%
  filter(StartPosition > 129875000 & StartPosition < 130000000)

pos<-annovar.filt$StartPosition

annovar.filt.functional <- annovar.filt %>%
  filter(CADDgt20 > 10 | WG_GWAVA_score > 0.4 | WG_EIGEN_score > 4)

annovar.filt.both <- annovar.filt 

#C. Filter vcf to these positions
#TODO: Change this so it doesn't have to be unzipped

filepath.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/8q24.cleaned.phased.vcf"
filepath.filtered.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/filtered_by_annotation/8q24.cleaned.phased.annotation.filtered.vcf"
hg.assembly<-"hg19"

#C. Convert filtered, phased vcf to geno matrix

filepath.filtered.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/filtered_by_annotation/8q24.cleaned.phased.annotation.filtered.vcf"
filepath.ped<-"/users/lgai/8q24_project/data/processed_data/gmkf_euro_completetrios_ids_match_vcf_mend_errors_removed.txt"
filepath.geno<-"/users/lgai/8q24_project/data/processed_data/geno_matrix/geno.phased.rds"
hg.assembly<-"hg19"

###########################################

get.subsetted.vcf<-function(filepath.annovar,filepath.vcf,filepath.filtered.vcf,hg.assembly){
  vcf <- readVcf(filepath.filtered.vcf,hg.assembly)
  geno(vcf)$GT<-apply(geno(vcf)$GT,2,function(x)gsub('\\|','/',x))
  geno(vcf)$GT<-apply(geno(vcf)$GT,2,function(x)gsub('1/0','0/1',x))
  vcf<-vcf[which(start(rowRanges(vcf)) %in% pos),]
  writeVcf(vcf,filepath.filtered.vcf)
}

get.geno<-function(filepath.ped,filepath.filtered.vcf,filepath.geno,hg.assembly){
  vcf <- readVcf(filepath.filtered.vcf,hg.assembly)
  ped<-read.table(filepath.ped,header=TRUE)
  geno <- vcf2geno(vcf,ped)
  saveRDS(geno, filepath.geno) 
}

###########################################













# Scratch

###########################################

vcf <- readVcf(filepath.vcf,hg.assembly)

#Change BEAGLE4.0's output to vcf style
geno(vcf)$GT<-apply(geno(vcf)$GT,2,function(x)gsub('\\|','/',x))
geno(vcf)$GT<-apply(geno(vcf)$GT,2,function(x)gsub('1/0','0/1',x))

length(which(start(rowRanges(vcf)) %in% pos)) #869

#Filter vcf to relevant positions
vcf2<-vcf[which(start(rowRanges(vcf)) %in% pos),]
vcf2
geno(vcf2)$GT

#Save vcf so you can transform it to BEAGLE v3 format for ScanTrio
writeVcf(vcf2,filepath.filtered.vcf)


###


vcf <- readVcf(filepath.filtered.vcf,hg.assembly)
ped<-read.table(filepath.ped,header=TRUE)
geno <- vcf2geno(vcf,ped)
saveRDS(geno, filepath.geno) 


###########################################

get.geno<-function(filepath.ped,filepath.filtered.vcf,filepath.geno,hg.assembly){
  vcf <- readVcf(filepath.filtered.vcf,hg.assembly)
  ped<-read.table(filepath.ped,header=TRUE)
  geno <- vcf2geno(vcf,ped)
  saveRDS(geno, filepath.geno) 
}

directory.name<-"filtered_by_annotation"
typeofdata<-"vcfs"
filename<-"8q24.cleaned.phased.annotation.filtered.vcf"
command<-paste0("mkdir /users/lgai/8q24_project/data/processed_data/",typeofdata,"/",directory.name,"/",filename)
#paste0("mkdir /users/lgai/8q24_project/data/processed_data/",typeofdata,"/",directory.name,"/")

system(command)

###########################################


#TODO: Add filtering to peak TDT region

############################# 0.5. Filter to windows ##########################

#i. window containing all variants around peak in gTDT
#Estimate: 
129750000
129875000-130000000

#Using positions from 90th percentile CADD
#peak.positions<-positions[positions > 129875000 && positions > 130000000]
positions<-caddscore_positions_90th$Position
length(positions)

peak.positions<-positions[positions > 129875000 & positions < 130000000]
12500
length(peak.positions)
151

100, 100 (overlapping windows)
75, 76 (non-overlapping windows) (of a different size, though)

#Filter to CADD score > 10 and within peak
vcf.geno<-filter.vcf.to.selected.pos(peak.positions,filepath.vcf)

#ii. windows of 100 variants around peak in gTDT, non-overlapping
#iii. windows of 100 variants of all variants in 8q24 region, non-overlapping
#iv. windows of 100 variants of 

#################

#TODO: Add in part about filtering to peak:

############################# 0.5. Filter to windows ##########################

#i. window containing all variants around peak in gTDT
#Estimate: 
129750000
129875000-130000000

#Using positions from 90th percentile CADD
#peak.positions<-positions[positions > 129875000 && positions > 130000000]
positions<-caddscore_positions_90th$Position
length(positions)

peak.positions<-positions[positions > 129875000 & positions < 130000000]
12500
length(peak.positions)
151

100, 100 (overlapping windows)
75, 76 (non-overlapping windows) (of a different size, though)

#Filter to CADD score > 10 and within peak
vcf.geno<-filter.vcf.to.selected.pos(peak.positions,filepath.vcf)

#ii. windows of 100 variants around peak in gTDT, non-overlapping
#iii. windows of 100 variants of all variants in 8q24 region, non-overlapping
#iv. windows of 100 variants of 

#################
#################
#################













####### Scratch
################# Unfiltered geno matrix  #########################

#Read in cleaned vcf
filepath.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/8q24.cleaned.vcf"
#filepath.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/8q24.cleaned.phased.vcf"
hg.assembly<-"hg19"
vcf <- readVcf(filepath.vcf,hg.assembly)

table(geno(vcf)$GT)
geno(vcf)$GT[1:10,1:10]

#Read in cleaned ped file
filepath.ped<-"/users/lgai/8q24_project/data/processed_data/gmkf_euro_completetrios_ids_match_vcf_mend_errors_removed.txt"
ped<-read.table(filepath.ped,header=TRUE)
geno(vcf)$GT<-apply(geno(vcf)$GT,2,function(x)gsub('\\|','/',x))
geno <- vcf2geno(vcf,ped)

#Save it
filepath.geno<-"/users/lgai/8q24_project/data/processed_data/geno_matrix/geno.cleaned.rds"
filepath.geno<-"/users/lgai/8q24_project/data/processed_data/geno_matrix/geno.phased.rds"
saveRDS(geno, filepath.geno)