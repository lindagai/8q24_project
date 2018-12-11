######### 1. Read in raw ped and vcf files and select relevant columns ################

#Read in vcf file; takes a minute

filepath_vcf<-"/dcl01/beaty/data/gmkf/euro/vcfs/filtered/8q24.recode.vcf"
chr8_vcf <- readVcf(filepath_vcf,"hg19")

#Read in pedfile
filepath_ped<-"/dcl01/beaty/data/gmkf/euro/peds/gmkf_euro_completetrios.csv"

#Create dataframe
chr8_ped <- read.csv(filepath_ped)
head(chr8_ped)

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


#Checks:

#Is there a difference between the id column and dbGaP.submitted.sample_id column?
#No

#Using completetrios.csv
#setdiff(chr8_ped$id, chr8_ped$dbGaP.submitted.sample_id)
#NULL

#PA1985B is an id in the vcf file. Does it appear in the pedfile?
#match("PA1985B",chr8_ped$pid)

#################### 2. Edit IDs in ped file to match vcf file ########################

######## Edit all PIDs to match VCF ID format ############

edit_pids<-function(pid){
  new_pid <-paste("H_TZ-",pid,"-",pid,sep="")
  new_pid
}

#Check a random PID: did it work?
#Before editing the PIDs to the "vcf" format
#match("PA1985",chr8_ped$pid)
##[1] 493

pids_vcf_style<-sapply(chr8_ped$pid,edit_pids)
head(pids_vcf_style)

#After editing the PIDs to the "vcf" format
#match("PA1985B",chr8_ped$pid)
#[1] NA

#NOTE:
#Some of the vcf file IDs have a B after them (Jing and Margaret don't know why)
#Need to add the B's manually

######## Get list of IDs from the VCF ############

#You only need to comment this back if if you need to re-create list of all vcf ids
# filepath_vcf<-"/dcl01/beaty/data/gmkf/euro/vcfs/filtered/8q24.recode.vcf"
# chr8_vcf <- readVcf(filepath_vcf,"hg19")
# vcf_ids<-colnames(geno(chr8_vcf)$GQ)
# head(vcf_ids)
# 
# vcf_ids_path<-"/users/lgai/8q24_project/data/processed_data/vcfs/8q24_vcf_ids.txt"
# write.table(vcf_ids,vcf_ids_path, sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)

######## Check which of the vcf IDs have B's in them ############

vcf_ids_path<-"/users/lgai/8q24_project/data/processed_data/vcfs/8q24_vcf_ids.txt"

vcf_ids<-as.data.frame(read.table(vcf_ids_path))
head(vcf_ids)

rows<-grep("*B$",vcf_ids$V1)
rows
vcf_B<-vcf_ids[rows,]
head(vcf_B)
length(vcf_B)
#23 IDs with Bs

#Check: are they in the same order in complete trios?

#chr8_ped$pid[rows]
#No

#Find every ID in test that appears with a B in the vcf and paste a B on it
#NOTE: no "duplicates", i.e. id either exists with B or not

#TODO: rewrite this to use sapply instead of a for loop

for(i in 1:length(vcf_B)){
  curr_id <-toString(vcf_B[i])
  #print("curr")
  #print(curr_id)
  
  search <-substr(curr_id,1,nchar(curr_id)-1)
  #print("search")
  #print(search)
  
  row<-match(search,pids_vcf_style)
  #print(row)
  #print(pids_vcf_style[row])
  
  pids_vcf_style[row]<-curr_id
  #print("    *****    ")
}

#Which rows have Bs in the IDs?
rows<-grep("*B$",chr8_ped$pid)
rows
#integer(0)
#chr8_ped[rows,]

#Check to see if vcf_ids$V1 and test/chr8_ped$pid matches
setdiff(pids_vcf_style, vcf_ids$V1)
#character(0)

#TODO: more error checking?

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
#chr8_ped[rows,1:5]

#Get the list of rows that end in B, AND are not children
#(i.e., we'd need to edit the fatid and motids for those rows)

#7 IDs
#TODO: Would be good to write some code to automate this, but you'll just have to get to it later

v1<-gsub("PA1578$","PA1578B",chr8_ped$motid)
v2<-gsub("PA1641$","PA1641B",v1)
v3<-gsub("PA2067$","PA2067B",v2)

chr8_ped$motid<-v3

w1<-gsub("PA1499$","PA1499B",chr8_ped$fatid)
w2<-gsub("PA1703$","PA1703B",w1)
w3<-gsub("PA2084$","PA2084B",w2)
w4<-gsub("PA2093$","PA2093B",w3)

chr8_ped$fatid<-w4

filepath_ped<-"/users/lgai/8q24_project/data/processed_data/gmkf_euro_completetrios_ids_match_vcf.txt"

write.table(chr8_ped, filepath_ped, sep=" ",row.names = FALSE,quote = FALSE)

#Checks
chr8_ped<-read.table(filepath_ped,header=TRUE)
head(chr8_ped)

########## 3. Remove 15 individuals from families with Mendelian errors  ################

#In Terminal:
#TODO: Fix this part to only use VariantAnnotation

#module load vcftools

#vcftools --vcf 8q24.recode.vcf --out /users/lgai/8q24_project/data/processed_data/vcfs/8q24.recode.15.mend.error.individuals.removed 

#Back in R

#Get rid of the IDs in the ped file
filepath_ped<-"/users/lgai/8q24_project/data/processed_data/gmkf_euro_completetrios_ids_match_vcf.txt"

chr8_ped<-read.table(filepath_ped,header=TRUE)
head(chr8_ped)
dim(chr8_ped)

#You only need to run this code if you are cleaning the file from Jing again
# 
# filepath_ids_with_mend_errors<-"/users/lgai/8q24_project/data/raw_data/vcfs/ids_with_mend_errors.txt"
# ids_with_mend_errors<-read.csv(filepath_ids_with_mend_errors,header=FALSE)
# 
# ids_with_mend_errors<-trimws(ids_with_mend_errors$V1)
# ids_with_mend_errors
# write.table(ids_with_mend_errors, filepath_ids_with_mend_errors, row.names=FALSE, col.names=FALSE)
#

filepath_ids_with_mend_errors<-"/users/lgai/8q24_project/data/raw_data/vcfs/ids_with_mend_errors.txt"

ids_with_mend_errors<-read.table(filepath_ids_with_mend_errors,skip=1)
print(ids_with_mend_errors)
#There should be 15 IDs with mendelian errors

head(chr8_ped)
dim(chr8_ped)
#996   6

chr8_ped_filt<-chr8_ped[!chr8_ped$pid %in% ids_with_mend_errors$V1,]
dim(chr8_ped_filt)
#981   6

filepath_ped<-"/users/lgai/8q24_project/data/processed_data/gmkf_euro_completetrios_ids_match_vcf_mend_errors_removed.txt"

write.table(chr8_ped_filt, filepath_ped, sep=" ",row.names = FALSE,quote = FALSE)
