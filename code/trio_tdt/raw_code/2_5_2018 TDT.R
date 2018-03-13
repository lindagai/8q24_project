#Description: 
# Does TDTs using trio on 8q24 region on GMKF trio data.

# Sections:
# 1. Read in raw ped and vcf files and select relevant columns
# 2. Edit IDs in ped file to match format of IDs in vcf file
# 3. Create genotype matrix for trio 
# 4. Do TDTs using trio

#Author: Linda Gai
#Last update: 2/5/18
#TO DO: 
# Run through all of it to make sure it works & send to Jackie
# Clean up code to get rid of cruft

################################################################################

#Sign into cluster
ssh -X lgai@jhpce01.jhsph.edu 

#Allocate memory
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G
module load R
R

######### 1. Read in raw ped and vcf files and select relevant columns ################
#Load R packages
library(VariantAnnotation)
library(trio)
library(dplyr)

#Read in vcf file; takes a minute
filepath_vcf<-"/dcl01/beaty/data/gmkf/euro/vcfs/filtered/8q24.recode.vcf"
chr8_vcf <- readVcf(filepath_vcf,"hg19")

#Read in pedfile
filepath_ped<-"/dcl01/beaty/data/gmkf/euro/peds/gmkf_euro_completetrios.csv"

#Create dataframe
chr8_ped <- read.csv(filepath_ped)

chr8_ped<-chr8_ped %>%
  rename("famid"="Family.ID",
         "pid"="Individual.ID",
         "fatid"="Father.ID",
         "motid"="Mother.ID",
         "sex" = "Gender",
         "affected"="Clinical.Status")

chr8_ped<-chr8_ped %>%
  select( "famid", "pid", "fatid", "motid", "sex","affected")

head(chr8_ped$pid)
head(chr8_ped)
#Order is:
#(family)
#child
#father
#mother

#In the trio documentation, it's
#father
#mother
#chlid

names(chr8_ped)

#Checks:

#Is there a differene between the id column and dbGaP.submitted.sample_id column?
#Using completetrios.csv
#setdiff(chr8_ped$id, chr8_ped$dbGaP.submitted.sample_id)
#NULL

#PA1985B is an id in the vcf file. Does it appear in the pedfile?
#match("PA1985B",chr8_ped$pid)

#################### 2. Edit IDs in ped file to match vcf file ########################

edit_pids<-function(pid){
  new_pid <-paste("H_TZ-",pid,"-",pid,sep="")
  new_pid
}

#Before editing the PIDs to the "vcf" format
#match("PA1985",chr8_ped$pid)
##[1] 493
#match("PA1985B",chr8_ped$pid)
#[1] NA

#Can you run vcf2geno if you change the format?
#No, Looks like some of the vcf file IDs have a B after them--why?
#Jing doesn't know, so you'll have to add them manually

pids_vcf_style<-sapply(chr8_ped$pid,edit_pids)
head(pids_vcf_style)

#Now check which of the vcf IDs have B's in them

#In Terminal:
#vcftools --vcf 8q24.recode.vcf --out 8q24.recode.GT --extract-FORMAT-info GT

#Get list of vcf IDs
#To create the list again, just comment this chunk back in

# vcf_gt_path<-"/users/lgai/8q24.recode.GT.GT.FORMAT"
# 
# vcf_gt<-as.data.frame(read.table(vcf_gt_path))
# vcf_ids<-vcf_gt[1,3:998]
# vcf_ids<-t(vcf_ids)
# dim(vcf_ids)
# head(vcf_ids)
# 
# vcf_ids_path<-file.path("/users","lgai","8q24.recode.vcfids.txt")
# 
# write.table(vcf_ids,vcf_ids_path, sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)

###########

vcf_ids_path<-file.path("/users","lgai","8q24.recode.vcfids.txt")

vcf_ids<-as.data.frame(read.table(vcf_ids_path))

############

rows<-grep("*B$",vcf_ids$V1)
rows
vcf_B<-vcf_ids[rows,] #26 such IDs
head(vcf_B)
length(vcf_B)
#26 for old sample.txt file
#23 for new sample.txt file! :D
#typeof(vcf_B[1])

#Check: are they in the same order in complete trios?
#chr8_ped$pid[rows]
#Nope

#Find every ID in test that appears with a B in the vcf and paste a B on it
#NOTE: no "duplicates", i.e. id either exists with B or not

for(i in 1:length(vcf_B)){
  curr_id <-toString(vcf_B[i])
  print("curr")
  print(curr_id)
  
  search <-substr(curr_id,1,nchar(curr_id)-1)
  #print("search")
  #print(search)
  
  row<-match(search,pids_vcf_style)
  #print(row)
  #print(pids_vcf_style[row])
  
  pids_vcf_style[row]<-curr_id
  #print("    *****    ")
}

#NOTE: All of the VCF IDs in the sample1.txt are in the ped file

#Which rows have Bs in the IDs?
rows<-grep("*B$",chr8_ped$pid)
rows
chr8_ped[rows,]

#Check to see if vcf_ids$V1 and test/chr8_ped$pid matches
setdiff(pids_vcf_style, vcf_ids$V1)
#character(0)

#Put pids_vcf_style into chr8_ped
chr8_ped$pid<-pids_vcf_style

########

#Rename father and mother ID columns
#Does the famid column need to be changed?

names(chr8_ped)
head(chr8_ped)

#Edit parent IDs in pedfile to match vcf file
edit_parent_ids<-function(pid){
  if (pid==0){
    "0"
  } else {
    new_pid <-paste("H_TZ-",pid,"-",pid,sep="")
    new_pid
  }
}

fatids_vcf_style<-sapply(chr8_ped$fatid,edit_parent_ids)
head(fatids_vcf_style)
head(chr8_ped$fatid)
chr8_ped$fatid<-fatids_vcf_style

motids_vcf_style<-sapply(chr8_ped$motid,edit_parent_ids)
head(motids_vcf_style)
head(chr8_ped$motid)
chr8_ped$motid<-motids_vcf_style

#Deal with ids with B first
#Which rows have Bs in the IDs?
rows<-grep("*B$",chr8_ped$pid)
rows
chr8_ped[rows,1:4]

#Get the list of IDs that end in B
vcf_B

#Which rows have Bs in the IDs?
rows<-grep("*B$",chr8_ped$pid)
rows
chr8_ped[rows,1:5]

#Get the list of rows that end in B, AND are not children
#(i.e., we'd need to edit the fatid and motids for those rows)

#7 IDs
#Would be good to write some code to automate this, but you'll just have to get to it later

v1<-gsub("PA1578$","PA1578B",chr8_ped$motid)
v2<-gsub("PA1641$","PA1641B",v1)
v3<-gsub("PA2067$","PA2067B",v2)

chr8_ped$motid<-v3

w1<-gsub("PA1499$","PA1499B",chr8_ped$fatid)
w2<-gsub("PA1703$","PA1703B",w1)
w3<-gsub("PA2084$","PA2084B",w2)
w4<-gsub("PA2093$","PA2093B",w3)

chr8_ped$fatid<-w4

write.csv(chr8_ped, "clean_pedfile_chr8.csv")

########## Remove 15 individuals from families with Mendelian errors ################

#In Terminal:

#module load vcftools
#vcftools --vcf 8q24.recode.vcf --out 8q24.recode.FILTERED --remove ids_with_mend_errors.txt --recode

#R

#Get rid of the IDs in the ped file
chr8_ped<-read.csv("clean_pedfile_chr8.csv")
chr8_ped<-as.data.frame(chr8_ped)

ids_with_mend_errors<-read.csv("ids_with_mend_errors.txt",header=FALSE)
print(ids_with_mend_errors)

ids_with_mend_errors<-trimws(ids_with_mend_errors$V1)
#ids_with_mend_errors
write.csv(ids_with_mend_errors, "ids_with_mend_errors.csv")

#head(chr8_ped)
dim(chr8_ped)
#996   7

chr8_filt<-chr8_ped[!chr8_ped$pid %in% ids_with_mend_errors,]
#dim(chr8_filt)
#981   7

write.csv(chr8_filt, "clean_pedfile_chr8_FILTERED.csv")

################# 4. Create genotype matrix for trio  ########################

#Read in the current version of the pedfile
chr8_ped<-as.matrix(read.table("clean_pedfile_chr8.csv"))
head(chr8_ped)

chr8_geno <- vcf2geno(chr8_vcf,chr8_ped)
write.csv(chr8_geno, "chr8_geno_19_11_2017.csv")

#################### 5. Do TDTs using trio  ########################

filepath_geno<-"/Users/lindagai/Documents/classes/3rd year/3rd term/Margaret:Ingo/2:2 update/chr8_geno_filt_12_06_2017.rda"

load(filepath_geno)

#A. Allelic TDT
allelicTDT(chr8_geno_filt)

#########

#B. Test the 1st SNP with a genotypic TDT
tdt(chr8_geno_filt[,1])

#########

#C. Test interaction between 1st and 2nd SNP with a genotypic TDT
tdtGxG(chr8_geno_filt[,1], chr8_geno_filt[,2])

#########

#D. Testing all SNPs in a Matrix in Genotype Format with a Genotypic TDT
colTDT(chr8_geno_filt)