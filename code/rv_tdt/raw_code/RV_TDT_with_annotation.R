#Description:
# Subsets 8q24 based on annotation information and runs RV TDT.
#
# Sections:
# 1. Filter 8q24 vcf
#         A. Positions with at least 1 score > 90th percentile
#
# 2. Process files into RV TDT format
#         A. tped
#         B. map
#         B. phen
# TODO: write step 2 as a sequence of function calls 
# 
# 3. Do RV TDT based on combinations of filters
#
# Input:
#
#
# Output:
#
#

#Author: Linda Gai
#Last update: 3/12/18
#TO DO: Clean up code to get rid of cruft

################################################################################

ssh -Y lgai@jhpce01.jhsph.edu

#Allocate memory
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G
module load R
R

################################################################################

############################## 1. Filter 8q24 vcf ##############################

############# A. Get file of chromosomes, relevant positions #################

filepath_annovar_report_updated_CADD<-file.path("/users","lgai","8q24_project","data","processed_data","annotation","8q24_annovar_reports_CADD13.txt")

#filepath_annovar_report_updated_CADD

annovar_report<-read.table(filepath_annovar_report_updated_CADD,sep="\t",quote ="",header=TRUE)

colnames(annovar_report)

#CADD:
#  90th percentile score: 7.7207
#filepath_cadd_8q24<-file.path("/users","lgai","2_23_18_filtering","cadd_8q24.tsv")
#cadd_8q24<-read.table(filepath_cadd_8q24,sep="\t",quote ="",header=TRUE)

library(dplyr)

cadd_8q24<-annovar_report[which(!is.na(annovar_report$CADD13)),] %>% 
  select("StartPosition",
         "CADD13")

head(cadd_8q24)
colnames(cadd_8q24)

cadd_90_quantile <- quantile(cadd_8q24$CADD13, probs = 0.9)[[1]] #90th quantile = 7.7207
cadd_90_quantile

caddscore_mat_90th<-cadd_8q24[cadd_8q24$CADD13>cadd_90_quantile,]
head(caddscore_mat_90th)
dim(caddscore_mat_90th) #1284    2
#1285    2 now, which makes sense

chr<-rep(8,nrow(caddscore_mat_90th))
length(chr)
cadd_for_vcf<-cbind(chr,caddscore_mat_90th$StartPosition)
colnames(cadd_for_vcf)<-c("Chr","Position")
head(cadd_for_vcf)
dim(cadd_for_vcf) #1284    2

filepath_caddscore_positions_90th<-file.path("/users","lgai","8q24_project","data","processed_data","annotation","caddscore_90th_percentile_positions_for_vcf.txt")

write.table(cadd_for_vcf, filepath_caddscore_positions_90th, sep="\t",row.names = FALSE,quote = FALSE)

#filepath_caddscore_positions_90th

############# B. Filter vcf using vcftools #################
#TODO: rewrite this part in R

#In Terminal, run vcf to filter

#vcftools --gzvcf 8q24.recode.FILTERED.FILTERED2.BEAGLE.vcf.gz --chr 8 --to-bp 129642137 --recode --out 8q24.recode.FILTERED.FILTERED2.BEAGLE.to.bp129642137

#Format:
#vcftools --vcf input_file.vcf --snps mySNPs.txt --recode --recode-INFO-all --out SNPs_only

#TODO: move vcf file to new directory, re-clean and re-filter

#Run in directory /users/lgai/2_23_18_filtering/filtered_data
vcftools --gzvcf /users/lgai/2_5_2018_RVTDT/8q24.recode.FILTERED.FILTERED2.BEAGLE.vcf.gz --positions /users/lgai/8q24_project/data/processed_data/annotation/caddscore_90th_percentile_positions_for_vcf.txt --recode --recode-INFO-all --out /users/lgai/8q24_project/data/processed_data/vcfs/8q24.caddscore_90th_percentile_positions_old

#read to see if it worked
library(VariantAnnotation)
filepath_vcf_cadd_filtered<-file.path("/users","lgai","8q24_project","data","processed_data","vcfs","8q24.caddscore_90th_percentile_positions_old.recode.vcf")

cadd_filtered_vcf<-readVcf(filepath_vcf_cadd_filtered)
#starts<-start(rowRanges(cadd_filtered_vcf))
#length(starts) #1284
#1283 unsure why -- is it because of the duplicated SNP at the beginning?

#################### 2. Process files into RV TDT format ############################################################

#TODO: The below code has not been edited.

module load vcftools

cd 2_5_2018_RVTDT

#Convert vcf files that have been phased by BEAGLE files into RV TDT format

#Convert vcf files into ped/map and tped
#Run from 2_5_2018_RVTDT directory
#vcftools --gzvcf 8q24.recode.FILTERED.FILTERED2.BEAGLE.vcf.gz  --plink --out phased_data/8q24.recode.8q24.recode.FILTERED.FILTERED2.BEAGLE
#vcftools --gzvcf 8q24.recode.FILTERED.FILTERED2.BEAGLE.vcf.gz --plink-tped --out phased_data/8q24.recode.8q24.recode.FILTERED.FILTERED2.BEAGLE

vcftools --vcf /users/lgai/2_23_18_filtering/filtered_data/8q24.caddscore_90th_percentile_positions.recode.vcf  --plink --out /users/lgai/2_23_18_filtering/filtered_data/cadd_rv_tdt/8q24.caddscore_90th_percentile_positions.recode
vcftools --vcf /users/lgai/2_23_18_filtering/filtered_data/8q24.caddscore_90th_percentile_positions.recode.vcf --plink-tped --out /users/lgai/2_23_18_filtering/filtered_data/cadd_rv_tdt/8q24.caddscore_90th_percentile_positions.recode

#Get the files in the correct format

#### tped ####
#RV TDT columns:
# 1) SNP/variant id
# 2-n) genotype on every individual (for the rest of the columns)

#PLINK tped files
#Columns:
# 2) rs or snp identifier
# 5-n) allele calls

#TO DO: would possibly be better to do it via R/Unix code directly

library(VariantAnnotation)

vcf<-readVcf("/users/lgai/2_23_18_filtering/filtered_data/8q24.caddscore_90th_percentile_positions.recode.vcf")

#List of genotypes/SNPs in the format n/n
vcfGeno<-geno(vcf)$GT # genotypes in "0|0"

dim(vcfGeno) #1284  981
#vcfGeno[1:5,1:5]

tped_ids<-colnames(vcfGeno)
writeLines(tped_ids, "/users/lgai/2_23_18_filtering/filtered_data/tpedFromMargaret_ids.txt")

#vcfGenoParents<-vcfGeno[,setdiff(1:ncol(vcfGeno), seq(1, ncol(vcfGeno), by=3))]
#dim(vcfGenoParents)
#vcfGenoParents[1:5,1:5]

#for every row of vcfGeno, replace | with tabs, separate every row with tab
vcfGenoBool<-apply(vcfGeno, 1, function(x) paste(gsub("\\|", "\t", x), collapse="\t"))

toWrite<-paste(rownames(vcfGeno), vcfGenoBool, sep="\t")
writeLines(toWrite, "/users/lgai/2_23_18_filtering/filtered_data/cadd_rv_tdt/8q24.caddscore_90th_percentile_positions.recode.rvtdt.tped")

#### map ####
# RV TDT columns:
#1.Gene ID (8q24)
#2. Variant id. The variant id must matches with the variant id in tped file
#3. MAF (?) #Can be obtained in PLINK, most similar to .frq file

#PLINK .frq file
#Columns
#1. CHR
#2. SNP
#3. A1
#4. A2
#5. MAF
#6. NM

#Select 2, 5

#Have to cd to plink folder first

#Get PLINK .frq file
./plink --file /users/lgai/2_23_18_filtering/filtered_data/cadd_rv_tdt/8q24.caddscore_90th_percentile_positions.recode --freq --out /users/lgai/2_23_18_filtering/filtered_data/cadd_rv_tdt/8q24.caddscore_90th_percentile_positions.recode

#Read in plink frq file
filepath_plink_frq<-"/users/lgai/2_23_18_filtering/filtered_data/cadd_rv_tdt/8q24.caddscore_90th_percentile_positions.recode.frq"

plink_frq<-read.table(filepath_plink_frq)

head(plink_frq) #Looks right

#Select the relevant columns
rvtdt_map<-plink_frq[,c(2,5)]
head(rvtdt_map) #Looks right

#Put ‘8q24’ as the 1st column.
v1<-rep('8q24',nrow(rvtdt_map))
rvtdt_map2<-cbind(v1,rvtdt_map)
head(rvtdt_map2)

#remove column names in 1st row
rvtdt_map2<-rvtdt_map2[-1,]
head(rvtdt_map2)

#Save
filepath_rvtdt_map<-"/users/lgai/2_23_18_filtering/filtered_data/cadd_rv_tdt/8q24.caddscore_90th_percentile_positions.recode.rvtdt.map"

#write.table(rvtdt_map2,filepath_rvtdt_map, sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)
write.table(rvtdt_map2,filepath_rvtdt_map, col.names=FALSE,row.names = FALSE,quote = FALSE)


#Check
rvtdt_map2<-read.table(filepath_rvtdt_map)
head(rvtdt_map2)

##### phen ####
#RV TDT Columns:
#1) sample ID
#2) family ID
#3) father ID
#4) mother ID
#5) sex (1 for male, 0 for female)
#6) case(1) / control(0).
#The order of the sample ID must match the order of individuals in tped file.

filepath_rvtdt_ped<-file.path("/users","lgai","data","8q24.recode.FILTERED.recode.BEAGLE.rvtdt.phen")

#Check the file
rvtdt_ped<-read.table(filepath_rvtdt_ped)
head(rvtdt_ped)
dim(rvtdt_ped) #981   6

####### 1) Make sure the order of the IDs in phen file match the order of IDs in tped file ########

#We can get this from the tfam file

#Read in the phen file and get the list of IDs
filepath_rvtdt_ped<-"/users/lgai/2_23_18_filtering/filtered_data/cadd_rv_tdt/8q24.caddscore_90th_percentile_positions.recode.rvtdt.phen"

rvtdt_ped<-read.table(filepath_rvtdt_ped)
head(rvtdt_ped)[1:5]
#phen_ids <-rvtdt_ped[,1]

#Check: are the IDs in the phen unique?
length(unique(rvtdt_ped$V1)) #981, so yes
nrow(rvtdt_ped) #981

#Read in the tfam file and get the list of IDs in the order of the tped file

filepath_plink_tfam<-"/users/lgai/2_23_18_filtering/filtered_data/cadd_rv_tdt/8q24.caddscore_90th_percentile_positions.recode.tfam"

tfam<-read.table(filepath_plink_tfam)
head(tfam)[1:5] #Definitely in a different order from ped/phen file

tped_ids <-as.data.frame(tfam[,1])
colnames(tped_ids) <- "V1"
head(tped_ids)

#check the length
nrow(tped_ids)
length(unique(tped_ids$V1))

library(dplyr)

#Get all the rows from phen that match the list of IDs from tped, in that order
reduced_phen <-left_join(tped_ids,rvtdt_ped)
head(reduced_phen)

head(tped_ids)

#Save
filepath_rvtdt_phen<-"/users/lgai/2_23_18_filtering/filtered_data/cadd_rv_tdt/8q24.caddscore_90th_percentile_positions.recode.rvtdt.phen"

#write.table(reduced_phen,filepath_rvtdt_phen, sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)
write.table(rvtdt_phen,filepath_rvtdt_phen,col.names=FALSE,row.names = FALSE,quote = FALSE)


#Check the file
rvtdt_phen<-read.table(filepath_rvtdt_phen)
head(rvtdt_phen)


###################### 3. Do RV TDT ###########################################################



##############

#filepath_plink_ped<-file.path("/users","lgai","2_5_2018_RVTDT","phased_data", "8q24.recode.8q24.recode.FILTERED.FILTERED2.BEAGLE.ped")

#Run RV TDT
/users/lgai/rv-tdt/rvTDT 8q24cadd.filtered -G ./2_23_18_filtering/filtered_data/cadd_rv_tdt/8q24.caddscore_90th_percentile_positions.recode.rvtdt.tped -P ./2_23_18_filtering/filtered_data/cadd_rv_tdt/8q24.caddscore_90th_percentile_positions.recode.rvtdt.phen \
-M ./2_23_18_filtering/filtered_data/cadd_rv_tdt/8q24.caddscore_90th_percentile_positions.recode.rvtdt.map \
--adapt 500 --alpha 0.00001 --permut 2000 \
--lower_cutoff 0 --upper_cutoff 100 \
--minVariants 3 \
--maxMissRatio 1


################################################################################