#Sign into cluster
ssh -X lgai@jhpce01.jhsph.edu 

#Allocate memory
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G
module load R
R

##################################################################

library(rvTDT)

##################################################################

get.rvTDT.results<-function(filepath.rvTDT.ped,filepath.evs,filepath.rvTDT.results){
  evs<-read.table(filepath.evs,sep="\t",quote ="",header=TRUE,stringsAsFactors=FALSE)
  rvTDT.ped<-read.table(filepath.rvTDT.ped,sep="\t",quote ="",header=TRUE,stringsAsFactors=FALSE)
  dim(rvTDT.ped) #981 4402...should be 1048
  dim(evs) #1048 is correct
  apply(evs, 2, function(x) any(is.na(x)))
  apply(rvTDT.ped, 2, function(x) any(is.na(x)))
  
  evs[1:5,1:4]
  rvTDT.ped[1:5,1:5]
  
  rvTDT.results <- rvTDT(rvTDT.ped,evs,maf.threshold=1)
  saveRDS(rvTDT.results, file=filepath.rvTDT.results)
}

##############################

filepath.rvTDT.ped.annotation<-"/users/lgai/8q24_project/data/processed_data/rvTDT/ped/rvTDT_ped_annotation.txt"
filepath.rvTDT.ped.peak<-"/users/lgai/8q24_project/data/processed_data/rvTDT/ped/rvTDT_ped_peak.txt"
filepath.rvTDT.ped.both<-"/users/lgai/8q24_project/data/processed_data/rvTDT/ped/rvTDT_ped_both.txt"
rvTDT.ped.filepaths<-c(filepath.rvTDT.ped.annotation,filepath.rvTDT.ped.peak,filepath.rvTDT.ped.both)

filepath.evs.annotation<-"/users/lgai/8q24_project/data/processed_data/rvTDT/evs/evs.annotation.txt"
filepath.evs.peak<-"/users/lgai/8q24_project/data/processed_data/rvTDT/evs/evs.peak.txt"
filepath.evs.both<-"/users/lgai/8q24_project/data/processed_data/rvTDT/evs/evs.both.txt"
evs.filepaths<-c(filepath.evs.annotation,filepath.evs.peak,filepath.evs.both)

#filepath.rvTDT.results<-"/users/lgai/8q24_project/data/processed_data/rvTDT/rvTDT_results.rds"
filepath.rvTDT.annotation<-"/users/lgai/8q24_project/data/processed_data/rvTDT/rvTDT_results_annotation.rds"
filepath.rvTDT.peak<-"/users/lgai/8q24_project/data/processed_data/rvTDT/rvTDT_results_peak.rds"
filepath.rvTDT.both<-"/users/lgai/8q24_project/data/processed_data/rvTDT/rvTDT_results_both.rds"
rvTDT.results<-c(filepath.rvTDT.annotation,filepath.rvTDT.peak,filepath.rvTDT.both)

i<-1

for (i in 1:3){
  filepath.rvTDT.ped<-rvTDT.ped.filepaths[i]
  filepath.rvTDT.ped
  filepath.evs<-evs.filepaths[i]
  filepath.rvTDT.results<-rvTDT.results[i]
  
  get.rvTDT.results(filepath.rvTDT.ped,filepath.evs,filepath.rvTDT.results)

}

readRDS(filepath.rvTDT.results)

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
