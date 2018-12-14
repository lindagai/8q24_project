#Sign into cluster
ssh -X lgai@jhpce01.jhsph.edu 

#Allocate memory
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G
module load R
R

##################################################################

get.rvTDT.ped<-function(filepath.geno,filepath.rvTDT.ped){
  geno<-readRDS(filepath.geno)
  dim(geno)  
  nrows<-nrow(geno)
  f.index<-seq(1,nrows,by=3)
  m.index<-seq(2,nrows,by=3)
  c.index<-seq(3,nrows,by=3) 
  rvTDT.ped<-rbind(geno[c.index,], geno[m.index,],geno[f.index,])
  dim(rvTDT.ped)
    
  write.table(rvTDT.ped,filepath.rvTDT.ped,sep="\t",quote=FALSE,row.name=FALSE)
}

##################################################################

#1. Convert trio geno matrix to rvTDT ped format

#Additional geno matrices
#filepath.geno<-"/users/lgai/8q24_project/data/processed_data/geno_matrix/geno.phased.rds"
filepath.geno.annotation <- "/users/lgai/8q24_project/data/processed_data/geno_matrix/filtered/geno_phased_annotation.rds"
filepath.geno.peak <-"/users/lgai/8q24_project/data/processed_data/geno_matrix/filtered/geno_phased_peak.rds"
filepath.geno.both <-"/users/lgai/8q24_project/data/processed_data/geno_matrix/filtered/geno_phased_both.rds"
geno.filepaths<-c(filepath.geno.annotation,filepath.geno.peak,filepath.geno.both)
geno.filepaths

#filepath.rvTDT.ped<-"/users/lgai/8q24_project/data/processed_data/rvTDT/ped/rvTDT_ped.txt"
filepath.rvTDT.ped.annotation<-"/users/lgai/8q24_project/data/processed_data/rvTDT/ped/rvTDT_ped_annotation.txt"
filepath.rvTDT.ped.peak<-"/users/lgai/8q24_project/data/processed_data/rvTDT/ped/rvTDT_ped_peak.txt"
filepath.rvTDT.ped.both<-"/users/lgai/8q24_project/data/processed_data/rvTDT/ped/rvTDT_ped_both.txt"
rvTDT.ped.filepaths<-c(filepath.rvTDT.ped.annotation,filepath.rvTDT.ped.peak,filepath.rvTDT.ped.both)

rvTDT.ped.filepaths

i<-1

for (i in 1:3){
  filepath.geno<-geno.filepaths[i]
  filepath.geno
  filepath.rvTDT.ped<-rvTDT.ped.filepaths[i]
  filepath.rvTDT.ped
  get.rvTDT.ped(filepath.geno,filepath.rvTDT.ped)
}

##############################















#Scratch

#1. Convert trio geno matrix to rvTDT ped format

filepath.geno<-"/users/lgai/8q24_project/data/processed_data/geno_matrix/geno.phased.rds"

geno<-readRDS(filepath.geno)

rvTDT.ped<-get.rvTDT.ped.from.trio.geno(geno)
rm(geno)

filepath.rvTDT.ped<-"/users/lgai/8q24_project/data/processed_data/rvTDT/rvTDT_ped.txt"
write.table(rvTDT.ped,filepath.rvTDT.ped,sep="\t",quote=FALSE,row.name=FALSE)

#Additional geno matrices
filepath.geno.annotation <- "/users/lgai/8q24_project/data/processed_data/geno_matrix/filtered/geno.phased.annotation.rds"
filepath.geno.peak <-"/users/lgai/8q24_project/data/processed_data/geno_matrix/filtered/geno.phased.peak.rds"
filepath.geno.both <-"/users/lgai/8q24_project/data/processed_data/geno_matrix/filtered/geno.phased.annotation.peak.rds"
filtered.geno.filepaths<-c(filepath.geno.annotation,filepath.geno.peak,filepath.geno.both)
