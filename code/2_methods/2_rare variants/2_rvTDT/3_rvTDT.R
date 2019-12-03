#Sign into cluster
ssh -X lgai@jhpce01.jhsph.edu 

#Allocate memory
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G
module load R
R

##################################################################

library(rvTDT)

#################### 1. Break the rvTDT files into windows ##########
n.window.sizes<-1
window.size=25
overlap=24
window.type<-"snps"

results.dir<-"/users/lgai/8q24_project/data/processed_data/rvTDT/"
filepath.evs <- "/users/lgai/8q24_project/data/processed_data/rvTDT/evs/8q24.cleaned.phased.filtered.annotation.rarevar.monomorphs.removed.evs"
filepath.ped <- "/users/lgai/8q24_project/data/processed_data/rvTDT/ped/8q24.cleaned.phased.filtered.annotation.rarevar.monomorphs.removed.ped"

get.rvTDT.files.in.windows(window.size,overlap,filepath.ped,filepath.evs,results.dir)

#3. Read the files in and run the test; save the results in one df
#This function and the above function could just be put together?
get.rvTDT.results.in.df(window.size,overlap,filepath.ped,filepath.evs,results.dir)
  
#4. Graph the files
filepath.rvTDT.results<-"/users/lgai/8q24_project/data/processed_data/rvTDT/window-size=25snps.overlap=24/rvTDT.results.txt"
results<-read.table(filepath.rvTDT.results,header=TRUE)
head(results)

scp lgai@jhpce01.jhsph.edu:/users/lgai/8q24_project/data/processed_data/rvTDT/window-size=25snps.overlap=24/rvTDT.results.txt '/Users/lindagai 1/Documents/classes/4th year/4th term'

results<-read.table(filepath.rvTDT.results,header=TRUE)












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
