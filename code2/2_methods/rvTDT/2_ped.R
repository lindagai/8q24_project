get.rvTDT.ped.from.trio.geno <-function(geno){
  nrows<-nrow(geno)
  f.index<-seq(1,nrows,by=3)
  m.index<-seq(2,nrows,by=3)
  c.index<-seq(3,nrows,by=3) 
  rbind(geno[c.index,], geno[m.index,],geno[f.index,])
}

##################################################################

#1. Convert trio geno matrix to rvTDT ped format

#should be:
filepath.geno<-"/users/lgai/8q24_project/data/processed_data/geno_matrix/geno.phased.rds"
geno<-readRDS(filepath.geno)
geno[1:5,1:5]

rvTDT.ped<-get.rvTDT.ped.from.trio.geno(geno)
rm(geno)

#filepath.rvTDT.ped<-"/Users/lindagai 1/Documents/classes/4th year/2nd term/Research/testing/rvTDT_ped.txt"
filepath.rvTDT.ped<-"/users/lgai/8q24_project/data/processed_data/rvTDT/rvTDT_ped.txt"
write.table(rvTDT.ped,filepath.rvTDT.ped,sep="\t",quote=FALSE,row.name=FALSE)