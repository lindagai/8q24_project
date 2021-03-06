ssh -Y lgai@jhpce01.jhsph.edu

#Allocate memory
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G
module load R
R

################

library(GenomicRanges)
library(VariantAnnotation)

################ 1. Use vcf2beagle to convert BEAGLEv4 vcf into a BEAGLE v3 file. ##############

#Download vcf2beagle
#Un-comment this out if you don't already have vcf2beagle.jar

# #Download vcf2beagle
# filepath.vcf2beagle<-"/users/lgai/vcf2beagle.jar"
# dl.vcf2beagle<-paste0("wget -O ",filepath.vcf2beagle," https://faculty.washington.edu/browning/beagle_utilities/vcf2beagle.jar")
# system(dl.vcf2beagle)
# 
#Convert BEAGLEv4.gz to BEAGLEv3
# filepath.BEAGLE.v4.gz <- "/users/lgai/scratch/data/8q24.recode.FILTERED.plink.nohalfcalls.DuplicatesRemoved.BEAGLE.vcf"
# filepath.BEAGLE.v3<-"/users/lgai/scratch/data/8q24.recode.FILTERED.plink.nohalfcalls.DuplicatesRemoved.BEAGLE.BEAGLEv3"

#VCF must be gzipped for this to work
filepath.vcf.gz<-"/users/lgai/8q24_project/data/processed_data/vcfs/rare_var_only/8q24.cleaned.phased.filtered.annotation.rarevar.monomorphs.removed.recode.vcf"
filepath.BEAGLE.v3<-"/users/lgai/8q24_project/data/processed_data/Scan-Trio/8q24.cleaned.phased.filtered.annotation.rarevar.monomorphs.removed.BEAGLEv3"

convert2beaglev3<-paste("zcat",filepath.vcf.gz, "|","java -jar /users/lgai/vcf2beagle.jar .", filepath.BEAGLE.v3)
system(convert2beaglev3)

################ 2. Create markerDF, which gives the corresponding position for each SNP name/rsID ##############

#Don't need to comment this back in unless you want to remake the markerDF

# TODO: Unsure that this is supposed to be from the vcf, or from something like annovar?
# #This can possibly be sped up by readGeno or something like that
# #Fix this later
# TODO: why are some of the alt alleles missing? do we use them?

#filepath_vcf<-"/users/lgai/scratch/data/8q24.recode.FILTERED.plink.nohalfcalls.DuplicatesRemoved.BEAGLE.vcf"
filepath.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/rare_var_only/8q24.cleaned.phased.filtered.annotation.rarevar.monomorphs.removed.recode.vcf"
vcf <- readVcf(filepath.vcf,"hg19")

vcf.rowranges<-rowRanges(vcf)
head(vcf.rowranges)
rsids<-names(rowRanges(vcf))
pos<-start(rowRanges(vcf))
ref<-as.character(ref(vcf))
alt<-as.character(unlist(alt(vcf)))

markerDF<-cbind(rsids,pos,ref,alt)
head(markerDF)

filepath.markerDF<-"/users/lgai/8q24_project/data/processed_data/Scan-Trio/rarevar.monomorphs.removed.markerDF.txt"

write.table(markerDF, filepath.markerDF, sep=" ",row.names = FALSE,quote = TRUE)
