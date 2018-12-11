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

################# Unfiltered geno matrix  #########################

#Read in cleaned vcf
filepath.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/8q24.cleaned.vcf"
filepath.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/8q24.cleaned.phased.vcf"
hg.assembly<-"hg19"
vcf <- readVcf(filepath.vcf,hg.assembly)

#Read in cleaned ped file
filepath.ped<-"/users/lgai/8q24_project/data/processed_data/gmkf_euro_completetrios_ids_match_vcf_mend_errors_removed.txt"
ped<-read.table(filepath.ped,header=TRUE)

geno <- vcf2geno(vcf,ped)

#Save it
filepath.geno<-"/users/lgai/8q24_project/data/processed_data/geno_matrix/geno.cleaned.rds"
filepath.geno<-"/users/lgai/8q24_project/data/processed_data/geno_matrix/geno.phased.rds"
saveRDS(geno, filepath.geno)

################# Filter geno matrix #########################

#A. Subset Annovar report to find the appropriate positions

filepath.sm.annovar<-"/users/lgai/8q24_project/data/processed_data/annotation/8q24_annovar_report_useful_variables.txt"

annovar_report<-read.table(filepath.sm.annovar,sep="\t",header=TRUE, quote ="")

annovar.filt<-annovar_report %>%
  filter(!is.na(Afalt_1000g2015aug_eur)) %>%
  filter(TotalDepth>20) %>%
  filter(CADDgt20 > 10 | WG_GWAVA_score > 0.4 | WG_EIGEN_score > 4)
dim(annovar.filt)
head(annovar.filt)

pos<-annovar.filt$StartPosition
length(pos) #1048

#B. Filter vcf to these positions
#TODO: Change this so it doesn't have to be unzipped
filepath.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/8q24.cleaned.phased.vcf"
hg.assembly<-"hg19"
vcf <- readVcf(filepath.vcf,hg.assembly)

length(which(start(rowRanges(vcf)) %in% pos)) #869

#Filter vcf to relevant positions
vcf2<-vcf[which(start(rowRanges(vcf)) %in% pos),]
vcf2
geno(vcf2)$GT

#Save vcf
filepath.filtered.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/filtered_by_annotation/8q24.cleaned.phased.annotation.filtered.vcf"
writeVcf(vcf2,filepath.filtered.vcf)

#C. Convert filtered, phased vcf to geno matrix
filepath.filtered.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/filtered_by_annotation/8q24.cleaned.phased.annotation.filtered.vcf"
hg.assembly<-"hg19"
vcf <- readVcf(filepath.filtered.vcf,hg.assembly)

table(geno(vcf)$GT)
# 0|0    0|1    1|0    1|1
# 730917  45856  44878  30838

#Change BEAGLE4.0's coding back to /, not |
geno(vcf)$GT<-apply(geno(vcf)$GT,2,function(x)gsub('\\|','/',x))

#table(geno(vcf)$GT)

filepath.ped<-"/users/lgai/8q24_project/data/processed_data/gmkf_euro_completetrios_ids_match_vcf_mend_errors_removed.txt"
ped<-read.table(filepath.ped,header=TRUE)

#TODO: Why are so many SNPs removed? Unsure why this did not work
geno <- vcf2geno(vcf,ped)
dim(geno)
geno(vcf)$GT[1:5,1:5]

#Save geno matrix
filepath.geno<-"/users/lgai/8q24_project/data/processed_data/geno_matrix/geno.phased.rds"
saveRDS(geno, filepath.geno)

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


#Scratch

none="0/0"
geno(vcf)$GT[1:5,1]!=none

, one=c("0/1"), both="1/1"

c("0/1")



table(geno(vcf)$GT)


if(removeNonBiallelic){
  ped
  geno <- matrix(-1, nrow=nrow(vcf), ncol=ncol(vcf), dimnames=dimnames(vcf))
  geno[1:5,1:5]
  
  ids.kid1 <- ped$fatid != 0
  mat.kid <- as.matrix(ped[ids.kid1, c("fatid", "motid", "pid")])
  
  vec.ids <- as.vector(t(mat.kid))
  
  
  mat.trio <- t(geno[,vec.ids])
  mat.trio[1:5,1:5]
  
  as.vector(colSums(mat.trio == -1, na.rm=TRUE))
  
  idsMore <- colSums(mat.trio == -1, na.rm=TRUE) > 0
  as.vector(idsMore)
  unique(idsMore)
  table(idsMore)
  if(any(idsMore)){
    mat.trio <- mat.trio[,!idsMore]
    warning("Since ", sum(idsMore), " of the SNVs show other/addtional genotypes than/to\n",
            "the ones specified by none, one, and both, these SNVs are removed.")
  }
  
  Warning messages:
    1: In vcf2geno(vcf2, ped) :
    Since 805 of the SNVs show other/addtional genotypes than/to
  the ones specified by none, one, and both, these SNVs are removed.
  2: In vcf2geno(vcf2, ped) :
    Since 21 of the SNVs were monomorphic, these SNVs were removed.
  > dim(geno)


###########################################

directory.name<-"filtered_by_annotation"
typeofdata<-"vcfs"
filename<-"8q24.cleaned.phased.annotation.filtered.vcf"
command<-paste0("mkdir /users/lgai/8q24_project/data/processed_data/",typeofdata,"/",directory.name,"/",filename)
#paste0("mkdir /users/lgai/8q24_project/data/processed_data/",typeofdata,"/",directory.name,"/")

system(command)