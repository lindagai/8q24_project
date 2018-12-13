#Description: 
#Sections

#1. Read in vcf file

#2. Delete:
#  tri-allelic calls
#  duplicated SNPs

#Inputs:
#Outputs:

#Author: Linda Gai
#Last update: 12/12/18

################################################################################

#Sign into cluster
ssh -X lgai@jhpce01.jhsph.edu 

#Allocate memory
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G
module load R
R

################################################################################

library(VariantAnnotation)

# Read in raw vcf
filepath.vcf<-"/dcl01/beaty/data/gmkf/euro/vcfs/filtered/8q24.recode.vcf"
hg.assembly<-"hg19"
vcf <- readVcf(filepath.vcf,hg.assembly)

table(geno(vcf)$GT)

#Convert missing values into correct format for BEAGLE
geno(vcf)$GT[geno(vcf)$GT == "."]<-"./."

# Remove tri-allelic SNPS/duplicate sites
vcf2<-vcf[-(which(start(rowRanges(vcf)) %in% start(rowRanges(vcf))[duplicated(start(rowRanges(vcf)))])),]

#Remove individuals with extra Mendelian errors
#TODO: Fix this txt file, it shouldn't have an x in it
filepath_ids_with_mend_errors<-"/users/lgai/8q24_project/data/raw_data/vcfs/ids_with_mend_errors.txt"
pids.with.mendelian.errors<-unlist(read.table(filepath_ids_with_mend_errors,skip=1,stringsAsFactors = FALSE),use.names = FALSE)
pids.with.mendelian.errors

pids<-colnames(geno(vcf)$PID)

vcf3<-vcf2[,-which(colnames(geno(vcf)$PID) %in% pids.with.mendelian.errors)]
vcf3

filepath.filtered.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/8q24.cleaned.vcf"
writeVcf(vcf3,filepath.filtered.vcf)

###################### create table of IDs and positions in the vcf  ##########################

pos<-start(rowRanges(chr8_vcf))
snp <- names(rowRanges(chr8_vcf))
snp.pos<-cbind(snp,pos)

#Remove duplicated rows
snp.pos<-snp.pos[!duplicated(snp.pos),]

#Save
filepath.vcf.snp.pos<-"/users/lgai/8q24_project/data/processed_data/vcf.snp.pos.txt"
write.table(snp.pos, filepath.vcf.snp.pos, row.names=FALSE,quote=FALSE)

################################################################################

#Scratch is below














#For Trio, you will need to filter to MAF > 0.01 or MAF > 0.05
chr8_geno_filt<-removeSNPs(chr8_geno, maf=0.01, perc.na =0.1)

################################################################################










#Other considerations?

#You should not have missing values like this if you phase it
# #Delete "half" missing values, where one of the alleles is called but other isn't
# geno(vcf)$GT[geno(vcf)$GT == "0/."]<-"."
# geno(vcf)$GT[geno(vcf)$GT == "1/."]<-"."

## For RVTDT, remove entries where all genotype calls are either "." or "0/0" 
#Gets the list of genotypes
# tmp<-geno(vcf)$GT
#dim(tmp)#24170   981

#Code MT used to remove the duplicate sites, sites where all calls were either "." or "0/0"
# vcfSub<-vcf[!(rowSums(tmp == "./.") + rowSums(tmp == "0/0") == ncol(vcf)),]
# vcfSub2<-vcfSub[-(which(start(rowRanges(vcfSub)) %in% start(rowRanges(vcfSub))[duplicated(start(rowRanges(vcfSub)))])),]

#TODO: Make sure all half-calls (0/., ./1) are set to missing
#vcfSub<-vcf[!(rowSums(tmp == "./.") + rowSums(tmp == "0/0") == ncol(vcf)),]

#Check: Why does the ped file not work?

#The pedfile you should have been using is clean_pedfile.csv?

#Check: Are there any weird values?
table(geno(vcf)$GT)
# .      0/0      0/1      1/1
# 3326561 19438019   903362   405378

#This could also be a function?

#######

#Improve this code with VariantAnnotation functions?

## Return all 'fixed' fields, "LAF" from 'info' and "GT" from 'geno'
> svp <- ScanVcfParam(info="LDAF", geno="GT")
> vcf1 <- readVcf(fl, "hg19", svp)
> names(geno(vcf1))

########################################

exons <- exons(txdb)
exons22 <- exons[seqnames(exons) == "chr22"]
seqlevelsStyle(exons22) <- "NCBI" ## match chrom names in VCF file
## Range-based filter:
withinRange <- function(rng)
  function(x) x
## The first filter identifies SNVs and the second applies the
## range restriction.
filters <- FilterRules(list(
  isSNV = isSNV,
  withinRange = withinRange(exons22)))
## Apply
## Not run:
filt1 <- filterVcf(fl, "hg19", tempfile(), filters=filters, verbose=TRUE)

########################################