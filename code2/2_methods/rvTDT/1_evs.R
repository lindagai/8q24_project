ssh -Y lgai@jhpce01.jhsph.edu

#Allocate memory
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G
module load R
R

################################

library(dplyr)

################################

get.evs<-function(filepath.annovar,filepath.evs,filepath.geno){
  filepath.annovar
  
  annovar<-read.table(filepath.annovar,sep="\t",header=TRUE, quote ="")
  annovar
  
  evs.raw<-annovar %>%
    mutate(geno0 = (1-Afalt_1000g2015aug_eur)^2 + 0.0000001) %>%
    mutate(geno1 = (1-Afalt_1000g2015aug_eur)*Afalt_1000g2015aug_eur + 0.0000001)  %>%
    mutate(geno2 = (Afalt_1000g2015aug_eur)^2*10000 + 0.0000001)  %>%
    dplyr::select("StartPosition",
           "geno2",
           "geno1",
           "geno0",
           "TotalDepth"
    )
  evs.raw[1:5,1:5]
  
  ################
  
  #Get positions in geno matrix
  # #TODO: Do this more efficiently :/
  geno<-readRDS(filepath.geno)
  snps<-colnames(geno) 
  
  filepath.vcf.snp.pos<-"/users/lgai/8q24_project/data/processed_data/vcf_snp_pos.txt"
  vcf.snp.pos<-read.table(filepath.vcf.snp.pos,sep="\t",quote ="",header=TRUE)
  head(vcf.snp.pos)
  pos<-vcf.snp.pos[vcf.snp.pos$snp %in% snps,"pos"]
  
  ################

  evs<-evs.raw[evs.raw$StartPosition %in% pos,]
  write.table(evs, filepath.evs, sep="\t",row.names = FALSE,quote = FALSE)
}

##############################

filepath.annovar.filt.annotation<-"/users/lgai/8q24_project/data/processed_data/annotation/filtered_annovar/8q24_annovar_report_functional_annnotation.txt"
filepath.annovar.filt.peak<-"/users/lgai/8q24_project/data/processed_data/annotation/filtered_annovar/8q24_annovar_report_TDT_peak.txt"
filepath.annovar.filt.both<-"/users/lgai/8q24_project/data/processed_data/annotation/filtered_annovar/8q24_annovar_report_filtered_by_both.txt"
annovar.filt.filepaths<-c(filepath.annovar.filt.annotation,filepath.annovar.filt.peak,filepath.annovar.filt.both)

filepath.evs.annotation<-"/users/lgai/8q24_project/data/processed_data/rvTDT/evs/evs.annotation.txt"
filepath.evs.peak<-"/users/lgai/8q24_project/data/processed_data/rvTDT/evs/evs.peak.txt"
filepath.evs.both<-"/users/lgai/8q24_project/data/processed_data/rvTDT/evs/evs.both.txt"
evs.filepaths<-c(filepath.evs.annotation,filepath.evs.peak,filepath.evs.both)

filepath.geno.annotation <- "/users/lgai/8q24_project/data/processed_data/geno_matrix/filtered/geno_phased_annotation.rds"
filepath.geno.peak <-"/users/lgai/8q24_project/data/processed_data/geno_matrix/filtered/geno_phased_peak.rds"
filepath.geno.both <-"/users/lgai/8q24_project/data/processed_data/geno_matrix/filtered/geno_phased_both.rds"

filtered.geno.filepaths<-c(filepath.geno.annotation,filepath.geno.peak,filepath.geno.both)

for (i in 1:3){
  filepath.annovar<-annovar.filt.filepaths[i]
  filepath.evs<-evs.filepaths[i]
  filepath.geno<-filtered.geno.filepaths[i]
  get.evs(filepath.annovar,filepath.evs,filepath.geno)
}

##############################















#Scratch





# vcf<-readVcf(filepath.vcf,"hg19")
# pos<- start(rowRanges(vcf))
# pos







> ## Return all 'fixed' fields, "LAF" from 'info' and "GT" from 'geno'
  filepath<-"/users/lgai/8q24_project/data/processed_data/vcfs/filtered_by_annotation/8q24.cleaned.phased.filtered.annotation.vcf"
fl <- system.file(file=filepath, package="VariantAnnotation")
svp <- ScanVcfParam(rowRanges="start")
vcf1 <- readVcf(filepath, "hg19", svp)
names(geno(vcf1))  

> ## Return all 'fixed' fields, "LAF" from 'info' and "GT" from 'geno'
  > svp <- ScanVcfParam(info="LDAF", geno="GT")
> vcf1 <- readVcf(fl, "hg19", svp)
> names(geno(vcf1))





filepath.sm.annovar<-"/users/lgai/8q24_project/data/processed_data/rvTDT/rvTDT_ped.txt"
#file.path("/users","lgai","8q24_project","data","processed_data","annotation","8q24_annovar_report.txt")

annovar_report<-read.table(filepath.sm.annovar,sep="\t",header=TRUE, quote ="")

annovar.filt<-annovar_report %>%
  filter(!is.na(Afalt_1000g2015aug_eur)) %>%
  filter(TotalDepth>20) %>%
  filter(CADDgt20 > 10 | WG_GWAVA_score > 0.4 | WG_EIGEN_score > 4)
dim(annovar.filt)

#TODO: This could be a function
evs.raw<-annovar.filt %>%
  mutate(geno0 = (1-Afalt_1000g2015aug_eur)^2) %>%
  mutate(geno1 = (1-Afalt_1000g2015aug_eur)*Afalt_1000g2015aug_eur)  %>%
  mutate(geno2 = (Afalt_1000g2015aug_eur)^2*10000)  %>%
  select("StartPosition",
         "geno2",
         "geno1",
         "geno0",
         "TotalDepth"
  )

head(evs.raw)

filepath.evs.raw<-"/users/lgai/rvTDT_test/evs_raw.txt"
write.table(evs.raw, filepath.evs.raw, sep="\t",row.names = FALSE,quote = FALSE)

filepath.cluster<-"/users/lgai/rvTDT_test/evs_raw.txt"
filepath.destination<-" '/Users/lindagai 1/Documents/classes/4th year/2nd term/Research/testing/' "

#On your laptop's R console, if you want to transfer it
scp.command<-paste0("scp lgai@jhpce01.jhsph.edu:", filepath.cluster, " ", filepath.destination)
scp.command
system(scp.command)

#################

#B. Subset to the part contained in your vcf
#TODO: eventually, you should do all this on the cluster
#TODO: add filtering

#Load the raw evs file
filepat.evs.raw<-"/Users/lindagai 1/Documents/classes/4th year/2nd term/Research/testing/evs_raw.txt"
evs.raw<-read.table(filepat.evs.raw,sep="\t",quote ="",header=TRUE,stringsAsFactors=FALSE)
head(evs.raw)

#Get the positions/rsIDs in the vcf
filepath.vcf.snp.pos<-"/Users/lindagai 1/Documents/classes/4th year/2nd term/Research/testing/vcf_snp_positions.txt"
vcf.snp.pos<-read.table(filepath.vcf.snp.pos, header=TRUE)
head(vcf.snp.pos)

#Subset the raw evs to the positions/rsIDs in the vcf and remove position #
evs.processed<-left_join(vcf.snp.pos,evs.raw, by = c("pos" = "StartPosition"))[,-2]
evs.processed$geno0<-evs.processed2$geno0 + 0.0000001
evs.processed$geno1<-evs.processed2$geno1 + 0.0000001
evs.processed$geno2<-evs.processed2$geno1 + 0.0000001
dim(evs.processed)

#filepath.rvTDT.ped<-"/Users/lindagai 1/Documents/classes/4th year/2nd term/Research/testing/rvTDT_ped.txt"
filepath.evs.processed<-"/users/lgai/8q24_project/data/processed_data/rvTDT/evs.processed.txt"
write.table(evs.processed,filepath.evs.processed,sep="\t",quote=FALSE,row.name=FALSE)