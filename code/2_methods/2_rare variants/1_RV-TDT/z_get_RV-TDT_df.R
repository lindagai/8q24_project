#Description: 

#Creates datasets for ggplot, GViz
#Adds position information
#adds indicator variables for all filtering info

#Sections

#TODO: This should probably be broken down further into different files


#Inputs:
#Outputs:

#Author: Linda Gai

################################################################################

library(dplyr)
library(ggplot2)

################################################################################


#1. Read in common variants dataset
  
#2. Read in RV-TDT results from every directory (graph RV-TDT results.R)
#TODO: Fix the RV-TDT code so it doesn't output 10s of little files? or it cleans up the files automatically

n.snps<-368
window.size<-25
n.windows<-ceiling(n.snps/window.size)
n.windows

data.dir<-"/users/lgai/8q24_project/data/processed_data/RV-TDT/"
df<-as.data.frame(matrix(nrow=n.windows,ncol=7))

#TODO: The gene name (8q24) is coerced to a 1 when you put it in the data frame--fix this!

#Read in all RV-TDT output to get the pvals
for (i in 1:n.windows){
  filepath.results<-paste0(data.dir, 
                           "scratch/window", window.size,"_no_overlap/window-size=",
                           window.size,"_snps_",i,"_pval/8q24.pval")
  # print(filepath.results)
  pval.df<-read.table(filepath.results,comment.char = "",header=TRUE)
  
  if (i==1){
    colnames(df)<-names(pval.df)
  }
  
  #print(head(pval.df))
  df[i,]<-pval.df[1,]
}

head(df)
df

#Your filepath
/users/lgai/8q24_project/data/processed_data/RV-TDT/scratch/window25_no_overlap/window-size=25_snps_1_pval/8q24.pval

#ii. Write out RV-TDT results in one .txt file
filepath.results<-"/users/lgai/8q24_project/data/processed_data/RV-TDT/scratch/window25_no_overlap/all_results.txt"

write.table(df,filepath.results, sep="\t", row.names=FALSE,quote=FALSE)

df2<-read.table(filepath.results,sep="\t",header=TRUE, quote ="")
head(df2)
df2

################################################################################

#2. Get the position of the middle of the window 
#TODO: Read in 
#a. get the "middle" positions of all the windows

#For each file, read in the .map

  #Creates map and tped files using the first n SNPs
  map<-read.table(filepath.map)
  tped<-read.table(filepath.tped)
  
  #Split the map and tped files into windows of size 100
  n.windows<-ceiling(nrow(map) / window.size)
  nrow(map)
  n.windows
  
  #TODO: make this more efficient/clean it up somehow
  filepath.vcf.snp.pos<-"/users/lgai/8q24_project/data/processed_data/vcf_snp_pos.txt"
  vcf.snp.pos<-read.table(filepath.vcf.snp.pos,sep="\t",quote ="",header=TRUE)
  head(vcf.snp.pos)
  
  i<-1
  
  df.window.pos<-as.data.frame(matrix(nrow=n.windows,ncol=1))
  df.window.pos

    
  for (i in 1:(n.windows)){
    #Get all positions in the window
    filepath.window.map.with.pos<-paste0(filepath.test,"window-size=",window.size,"_snps","_",i,".map.with.pos.txt")
    map.with.pos<-as.data.frame(read.table(filepath.window.map.with.pos,sep="\t",quote ="",header=FALSE))
    head(map.with.pos)
    window.pos<-sort(as.numeric(map.with.pos$V4))
    window.pos
    
    #Get the largest and smallest, and then set the window's position as halfway between them
    start.pos<-window.pos[1]
    end.pos<-tail(window.pos, 1)
    mid.pos<- start.pos + (end.pos - start.pos)/2
    mid.pos
    df.window.pos[i,]<-mid.pos
  }
  
colnames(df.window.pos)<-c("Position")
filepath.window.mid.pos<-"/users/lgai/8q24_project/data/processed_data/RV-TDT/scratch/window25_no_overlap/window.midpoints.txt"
write.table(df.window.pos,filepath.window.mid.pos,sep="\t", row.names=FALSE,quote=FALSE)

filepath.window.mid.pos<-"/users/lgai/8q24_project/data/processed_data/RV-TDT/scratch/window25_no_overlap/window.midpoints.txt"
df.window.pos  <-read.table(filepath.window.mid.pos,sep="\t",header=TRUE, quote ="")
head(df.window.pos)

df2<-read.table(filepath.results,sep="\t",header=TRUE, quote ="")
head(df2)
df2

#Create new df with results and proper window position
df.results.with.pos<-cbind(df2,df.window.pos)
head(df.results.with.pos)
filepath.df.results.with.pos<-"/users/lgai/8q24_project/data/processed_data/RV-TDT/scratch/window25_no_overlap/results.with.pos.txt"

# Save
write.table(df.results.with.pos,filepath.df.results.with.pos,sep="\t", row.names=FALSE,quote=FALSE)

#Copy over to laptop
#This one works
#????
scp lgai@jhpce01.jhsph.edu:"/users/lgai/8q24_project/data/processed_data/RV-TDT/scratch/window25_no_overlap/results.with.pos.txt" '/Users/lindagai 1/Documents/classes/4th year/4th term'

################################################################################

