#Description: 
#Sections

#1. Read in vcf file

#2. Delete:
#  half-missing calls
#  tri-allelic calls
#  duplicated SNPs

#3. Use BEAGLE to phase

#4. Get the files into RV TDT format

#Inputs:
#Outputs:

#Author: Linda Gai
#Last update: 2/5/18

#TO DO: Clean up code to get rid of cruft

##########################################################

#Sign into cluster
ssh -X lgai@jhpce01.jhsph.edu 

#Allocate memory
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G
module load R
R

##########################################################

library(VariantAnnotation)

#1. Read in vcf file and delete:
#  half-missing calls
#  tri-allelic calls
#  duplicated SNPs

vcf<-readVcf("/users/lgai/8q24.recode.FILTERED.recode.vcf")

#Make sure missing sites are coded as "./.", not "."
geno(vcf)$GT[geno(vcf)$GT == "."]<-"./."

## remove entries where all genotype calls are either "." or "0/0" 
#Gets the list of genotypes
tmp<-geno(vcf)$GT
#dim(tmp)#24170   981

#Code MT used to remove the duplicate sites, sites where all calls were either "." or "0/0"
vcfSub<-vcf[!(rowSums(tmp == "./.") + rowSums(tmp == "0/0") == ncol(vcf)),]
vcfSub2<-vcfSub[-(which(start(rowRanges(vcfSub)) %in% start(rowRanges(vcfSub))[duplicated(start(rowRanges(vcfSub)))])),]

#TODO: Make sure all half-calls (0/., ./1) are set to missing
#vcfSub<-vcf[!(rowSums(tmp == "./.") + rowSums(tmp == "0/0") == ncol(vcf)),]

writeVcf(vcfSub2,"/users/lgai/8q24.recode.FILTERED.FILTERED2.vcf")

#Check: Why does the ped file not work?

vcf<-readVcf("/users/lgai/8q24.recode.FILTERED.FILTERED2.vcf")

#The pedfile you should have been using is clean_pedfile.csv?

############# SCRATCH #######################

#TEST: Does BEAGLE like it if you just save the same vcf without changing anything?
writeVcf(vcf, "/users/lgai/8q24.recode.FILTERED.recode.writeVcf2.vcf")
  
#Try reading in the old vcf file--how is it coding the missing values?

period_table<-tmp[grep("./.", rownames(tmp)), ]
dim(period_table)

#rowSums(tmp == ".")[1:2] #no. of rows s.t. all entries are missing

#rowSums(tmp == "0/0")[1:2] #no. of rows s.t. all calls are 0/0

#gets the rows where the genotype = "0/1"?
rowSums(tmp == "0/0")[1:2]
tmp[967]

dot<-which(tmp[,] == ".")
dim(dot)
NULL

dot<-which(tmp == "0/0")

#tmp[grep(".", rownames(tmp)), ]

###################################

#3. Use BEAGLE to phase

ped=[file] specifies a Linkage-format pedigree file for specifying family relationships.
The pedigree file has one line per individual. The first 4 white-space delimited fields of
each line are 1) pedigree ID, 2) individual ID, 3) father’s ID, and 4) mother’s ID. A “0” is
used in column 3 or 4 if the father or mother is unknown. The individual IDs are required
to be unique. Beagle uses the data in columns 2-4 to identify parent-offspring duos and
trios in the input data. Any or all columns of the pedigree file after column 4 may be
omitted. See also the duoscale and trioscale parameters.

filepath_rvtdt_ped<-file.path("/users","lgai","data","8q24.recode.FILTERED.recode.BEAGLE.rvtdt.phen")

java -Xmx10000m -jar beagle.r1399.jar gt=/users/lgai/8q24.recode.FILTERED.FILTERED2.vcf ped=/users/lgai/data/8q24.recode.FILTERED.recode.BEAGLE.rvtdt.phen out=/users/lgai/2_2_2017_RVTDT/8q24.recode.FILTERED.FILTERED2.BEAGLE
#Exception in thread "main" java.lang.IllegalArgumentException: missing genotype separator (.): 

8q24.recode.FILTERED.recode.vcf

java -Xmx10000m -jar beagle.r1399.jar gt=8q24.recode.FILTERED.recode.vcf ped=/users/lgai/data/8q24.recode.FILTERED.recode.BEAGLE.rvtdt.phen out=/users/lgai/2_2_2017_RVTDT/8q24.recode.FILTERED.FILTERED2.BEAGLE
#This worked, so it must be something in the previous command

java -Xmx10000m -jar beagle.r1399.jar gt=/users/lgai/8q24.recode.FILTERED.recode.writeVcf.vcf ped=/users/lgai/data/8q24.recode.FILTERED.recode.BEAGLE.rvtdt.phen out=/users/lgai/2_2_2017_RVTDT/8q24.recode.FILTERED.FILTERED2.BEAGLE
#This did NOT work, so it has to do with writeVcf?

java -Xmx10000m -jar beagle.r1399.jar gt=/users/lgai/8q24.recode.FILTERED.recode.writeVcf2.vcf ped=/users/lgai/data/8q24.recode.FILTERED.recode.BEAGLE.rvtdt.phen out=/users/lgai/2_2_2017_RVTDT/8q24.recode.FILTERED.FILTERED2.BEAGLE
#this worked. so just make sure you change the missing values before writing in

8q24.recode.FILTERED.FILTERED2.recodeMiss.vcf 

java -Xmx10000m -jar beagle.r1399.jar gt=/users/lgai/8q24.recode.FILTERED.FILTERED2.recodeMiss.vcf ped=/users/lgai/data/8q24.recode.FILTERED.recode.BEAGLE.rvtdt.phen out=/users/lgai/2_2_2017_RVTDT/8q24.recode.FILTERED.FILTERED2.BEAGLE


cd /users/lgai/data/8q24.recode.FILTERED.recode.BEAGLE.rvtdt.phen



############# TEST: How to get BEAGLE to recognize family structure #####

java -Xmx10000m -jar beagle.r1399.jar gt=/users/lgai/8q24.recode.FILTERED.FILTERED2.vcf ped=/users/lgai/data/8q24.recode.FILTERED.recode.BEAGLE.rvtdt.phen out=/users/lgai/2_2_2017_RVTDT/8q24.recode.FILTERED.FILTERED2.BEAGLE
#Did not work

##############################

java -Xmx10000m -jar beagle.r1399.jar gt=/users/lgai/8q24.recode.FILTERED.FILTERED2.vcf ped=clean_pedfile_chr8_FILTERED.txt out=/users/lgai/2_2_2017_RVTDT/8q24.recode.FILTERED.FILTERED2.BEAGLE
#This worked

#Crashed before it finished :/


############## Currently running #################
java -Xmx10000m -jar beagle.r1399.jar gt=/users/lgai/8q24.recode.FILTERED.FILTERED2.vcf ped=clean_pedfile_chr8_FILTERED.txt out=/users/lgai/2_5_2018_RVTDT/8q24.recode.FILTERED.FILTERED2.BEAGLE

##############################

#4. Get the files into RV TDT format and run RV TDT

##################### RV TDT #####################
mkdir phased_data
module load vcftools

cd 2_5_2018_RVTDT

#Convert vcf files that have been phased by BEAGLE files into RV TDT format

#Convert vcf files into ped/map and tped
#Run from 2_5_2018_RVTDT directory
vcftools --gzvcf 8q24.recode.FILTERED.FILTERED2.BEAGLE.vcf.gz  --plink --out phased_data/8q24.recode.8q24.recode.FILTERED.FILTERED2.BEAGLE
vcftools --gzvcf 8q24.recode.FILTERED.FILTERED2.BEAGLE.vcf.gz --plink-tped --out phased_data/8q24.recode.8q24.recode.FILTERED.FILTERED2.BEAGLE

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

#module load R
#R

#Read in file
filepath_plink_tped<-file.path("/users","lgai","2_5_2018_RVTDT","phased_data", "8q24.recode.8q24.recode.FILTERED.FILTERED2.BEAGLE.tped")

plink_tped<-read.table(filepath_plink_tped)
head(plink_tped)[1:5] #Looks right

#Select the relevant columns
rvtdt_tped<-plink_tped[,c(2, 5:ncol(plink_tped))]
head(rvtdt_tped)[1:5] #Looks right

#Save
filepath_rvtdt_tped<-file.path("/users","lgai","2_5_2018_RVTDT","phased_data","8q24.recode.8q24.recode.FILTERED.FILTERED2.BEAGLE.rvtdt.tped")

write.table(rvtdt_tped,filepath_rvtdt_tped, sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)

#Read and file to check
rvtdt_tped<-read.table(filepath_rvtdt_tped)
head(rvtdt_tped)[1:8]
ncol(rvtdt_tped) #1963

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
./plink --file /users/lgai/2_5_2018_RVTDT/phased_data/8q24.recode.8q24.recode.FILTERED.FILTERED2.BEAGLE --freq --out /users/lgai/2_5_2018_RVTDT/phased_data/8q24.recode.8q24.recode.FILTERED.FILTERED2.BEAGLE

#Read in plink frq file
filepath_plink_frq<-file.path("/users","lgai","2_5_2018_RVTDT","phased_data", "8q24.recode.8q24.recode.FILTERED.FILTERED2.BEAGLE.frq")

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
filepath_rvtdt_map<-file.path("/users","lgai","2_5_2018_RVTDT","phased_data", "8q24.recode.8q24.recode.FILTERED.FILTERED2.BEAGLE.rvtdt.map")

write.table(rvtdt_map2,filepath_rvtdt_map, sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)

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

#PLINK ped columns
The PED file is a white-space (space or tab) delimited file: the first six columns are mandatory:
#Family ID
#Individual ID
#Paternal ID
#Maternal ID
#Sex (1=male; 2=female; other=unknown)
#Phenotype [1] 0 0 0 0 0

plink_ped[1:5,4]
[1] 0 0 0 0 0

> plink_ped[1:5,3]
[1] 0 0 0 0 0

> plink_ped[1:5,1]
[1] 0 0 0 0 0

> plink_ped[1:5,6]
[1] 0 0 0 0 0

#TEST: use the ped file instead of Use clean_pedfile_chr8_FILTERED.csv?
#NO it doesn't work

#This won't work b/c they are all 0's...
#filepath_plink_ped<-file.path("/users","lgai","2_5_2018_RVTDT","phased_data", "8q24.recode.8q24.recode.FILTERED.FILTERED2.BEAGLE.ped")

#Use clean_pedfile_chr8_FILTERED.csv

#########################

#Read in file
filepath_plink_ped<-file.path("/users","lgai","clean_pedfile_chr8_FILTERED.txt")
plink_ped<-as.data.frame(read.table(filepath_plink_ped))
names(plink_ped)
head(plink_ped)

#Switch the order of famid and pid
#        "famid","pid","fatid","motid","sex","affected"
rvtdt_ped<-cbind(plink_ped[,2,drop=FALSE],plink_ped[,1,drop=FALSE],plink_ped[,3:6])
head(rvtdt_ped)

#Change coding of...

#sex (1 for male, 0 for female)
v1<-gsub("female","0",rvtdt_ped$V5)
v2<-gsub("male","1",v1)
head(v2)
rvtdt_ped$V5<-v2
head(rvtdt_ped)

#case(1) / control(0)
v1<-gsub("Unaffected","0",rvtdt_ped$V6)
head(v1)
v2<-gsub("Affected","1",v1)
head(v2)
rvtdt_ped$V6<-v2

#Save
filepath_rvtdt_ped<-file.path("/users","lgai","data","8q24.recode.FILTERED.recode.BEAGLE.rvtdt.phen")
#filepath_rvtdt_ped

write.table(rvtdt_ped,filepath_rvtdt_ped, sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)

#Check the file
#rvtdt_ped<-read.table(filepath_rvtdt_ped)
#head(rvtdt_ped)

####### Make sure the order of the IDs in phen file match the order of IDs in tped file ########

#We can get this from the tfam file

#Read in the phen file and get the list of IDs
filepath_rvtdt_ped<-file.path("/users","lgai","data","8q24.recode.FILTERED.recode.BEAGLE.rvtdt.phen")

rvtdt_ped<-read.table(filepath_rvtdt_ped)
head(rvtdt_ped)[1:5]
#phen_ids <-rvtdt_ped[,1]

#Check: are the IDs in the phen unique?
length(unique(rvtdt_ped$V1)) #981, so yes
nrow(rvtdt_ped) #981

#Read in the tfam file and get the list of IDs in the order of the tped file

filepath_plink_tfam<-file.path("/users","lgai","2_5_2018_RVTDT","phased_data", "8q24.recode.8q24.recode.FILTERED.FILTERED2.BEAGLE.tfam")

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
filepath_rvtdt_phen<-file.path("/users","lgai","2_5_2018_RVTDT","phased_data","8q24.recode.8q24.recode.FILTERED.FILTERED2.BEAGLE.phen")

write.table(reduced_phen,filepath_rvtdt_phen, sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)

#Check the file
rvtdt_phen<-read.table(filepath_rvtdt_phen)
head(rvtdt_phen)

##############

#filepath_plink_ped<-file.path("/users","lgai","2_5_2018_RVTDT","phased_data", "8q24.recode.8q24.recode.FILTERED.FILTERED2.BEAGLE.ped")

#Run RV TDT
/users/lgai/rv-tdt/rvTDT 8q24test6 -G ./2_5_2018_RVTDT/phased_data/8q24.recode.8q24.recode.FILTERED.FILTERED2.BEAGLE.rvtdt.tped -P ./2_5_2018_RVTDT/phased_data/8q24.recode.8q24.recode.FILTERED.FILTERED2.BEAGLE.phen \
-M ./2_5_2018_RVTDT/phased_data/8q24.recode.8q24.recode.FILTERED.FILTERED2.BEAGLE.map \
--adapt 500 --alpha 0.00001 --permut 2000 \
--lower_cutoff 0 --upper_cutoff 100 \
--minVariants 3 \
--maxMissRatio 1

#TO DO:

#The ped file you are using
#/2_5_2018_RVTDT/phased_data/8q24.recode.8q24.recode.FILTERED.FILTERED2.BEAGLE.phen
#
#Will flip-scan work with this one? NO
  
  #Use plink flip-scan to get the MAF to see if you can get the RV TDT to work
  
  #In Terminal:
  
  ./plink --vcf /users/lgai/8q24.recode.FILTERED.recode.vcf --vcf-half-call m --recode --out /users/lgai/data/8q24.recode.FILTERED.plink.nohalfcalls

#Get the PLINK version of the vcf file you want
./plink --gzvcf /users/lgai/8q24.recode.FILTERED.FILTERED2.BEAGLE.vcf.gz --recode --out /users/lgai/data/8q24.recode.FILTERED.plink.nohalfcalls

#Convert it to binary format

#Use PLINK's flipscan to get the major/minor alleles
plink --bfile mydata --flip-scan --out flip

./plink --file /users/lgai/2_5_2018_RVTDT/phased_data/8q24.recode.8q24.recode.FILTERED.FILTERED2.BEAGLE --flip-scan --out /users/lgai/2_5_2018_RVTDT/phased_data/flip
#Did not work, using the ped file generated by PLINK
#It just has 1000s of columns, only 1st 2 are non-zero

#Try Ped file with binary coding
./plink --file /2_5_2018_RVTDT/phased_data/8q24.recode.8q24.recode.FILTERED.FILTERED2.BEAGLE.phen
 --flip-scan --out /users/lgai/2_5_2018_RVTDT/phased_data/flip

#It needs the map file
#Try renaming it to ped

filepath_rvtdt_ped<-file.path("/users","lgai","2_5_2018_RVTDT", "phased_data","8q24.recode.8q24.recode.FILTERED.FILTERED2.BEAGLE.ped")
rvtdt_ped<-read.table(filepath_rvtdt_ped)
dim(rvtdt_ped)
head(rvtdt_ped)[1:5,1:6]

#So the ped file isn't correct

##### phen ####
#RV TDT Columns:
#1) sample ID
#2) family ID
#3) father ID
#4) mother ID
#5) sex (1 for male, 0 for female)
#6) case(1) / control(0).
#The order of the sample ID must match the order of individuals in tped file.

#PLINK ped columns
The PED file is a white-space (space or tab) delimited file: the first six columns are mandatory:
  #Family ID
  #Individual ID
  #Paternal ID
  #Maternal ID
  #Sex (1=male; 2=female; other=unknown)
  #Phenotype

#Use clean_pedfile_chr8_FILTERED.csv
#filepath_plink_ped<-file.path("/users","lgai","clean_pedfile_chr8_FILTERED.txt")
#Did not work -- PLINK flip-scan expects more columns than the 1st 6
  #now does it work?
#No, pedfile still has to have the extra colums
#filepath_rvtdt_ped<-file.path("/users","lgai","2_5_2018_RVTDT","phased_data","8q24.recode.8q24.recode.FILTERED.FILTERED2.BEAGLE.ped")
#write.table(rvtdt_ped,filepath_rvtdt_ped, sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)



cp thefile thecopy
cp 8q24.recode.8q24.recode.FILTERED.FILTERED2.BEAGLE.phen 8q24.recode.8q24.recode.FILTERED.FILTERED2.BEAGLE.ped

#########THIS WORKS

#bim file 

#Try to get get A1 (minor) and A2 (major) via bim file

#./plink --file mydata --out mydata --make-bed
./plink --file /users/lgai/2_5_2018_RVTDT/phased_data/8q24.recode.8q24.recode.FILTERED.FILTERED2.BEAGLE --out /users/lgai/2_5_2018_RVTDT/phased_data/binary --make-bed


#1  Chromosome code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or name
#2  Variant identifier
#3  Position in morgans or centimorgans (safe to use dummy value of '0')
#4  Base-pair coordinate (normally 1-based, but 0 ok; limited to 231-2)
#5  Allele 1 (corresponding to clear bits in .bed; usually minor)
#6  Allele 2 (corresponding to set bits in .bed; usually major)

###
#Read in file
filepath_plink_bim<-file.path("/users","lgai","2_5_2018_RVTDT","phased_data", "binary.bim")
plink_bim<-read.table(filepath_plink_bim)

head(plink_bim)[1:5] #Looks right

#take the SNP (2), A1-minor (5), A2-major (6) fields
snp_mat <-plink_bim[,c(2,5,6)]
head(snp_mat)
dim(snp_mat)

#The first column is SNP/variant id, and followed
#by the genotype on every individual.
#Since RV-TDT only takes phased data,
#every two columns present the haplotypes of one individual

#irisSubset <- iris[grep("osa", iris$Species), ]
#Subset <- snp_mat[grep(currSNP, snp_mat), ]

filepath_tped<-file.path("/users","lgai","2_5_2018_RVTDT","phased_data", "8q24.recode.8q24.recode.FILTERED.FILTERED2.BEAGLE.rvtdt.tped")
tped<-read.table(filepath_tped)
tped[1:5,1:5]
i=1

new_tped_mat<-matrix(, nrow = nrow(tped), ncol = ncol(tped))
new_tped_mat[1:5,1:5]

for (i in 1:nrow(tped)){
  currSNP = toString(tped[i,1])
  #currSNP
  new_tped_mat[i,1]<-currSNP
  #tped[i,1:5]
  currSNP_info = snp_mat[grep(currSNP, snp_mat[,1]), ]
  #currSNP_info
  minor=toString(currSNP_info[[2]])
  #minor
  major=toString(currSNP_info[[3]])
  #major
#  test<-tped[i,1:10]
#  test
  currRowSNPs<-unlist(tped[i,2:length(tped[i,])])
  #currRowSNPs[1:5]
  tped_row1<-gsub(minor,1,currRowSNPs)
  #tped_row1[1:5]
  tped_row2<-gsub(major,0,tped_row1)
  #tped_row2[1:5]
  #tped_mat[i,2:length(tped[i,])]<-tped_row2
  new_tped_mat[i,2:length(tped[i,])]<-tped_row2
  #new_tped_mat[i,1:5]
}

#toString(minor)

#Change this to lapply later
#lapply(function, matrix)

#rvtdt_bool_tped

#Save
filepath_rvtdt_bool_tped<-file.path("/users","lgai","2_5_2018_RVTDT","phased_data","8q24.recode.8q24.recode.FILTERED.FILTERED2.BEAGLE.rvtdt.bool.tped")

write.table(new_tped_mat,filepath_rvtdt_bool_tped, sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)

#Read and file to check
#rvtdt_bool_tped<-read.table(filepath_rvtdt_bool_tped)
#head(rvtdt_bool_tped)[1:8]
#ncol(rvtdt_bool_tped) #1963

#Run RV TDT
/users/lgai/rv-tdt/rvTDT 8q24test7 -G ./2_5_2018_RVTDT/phased_data/tpedFromMargaret.tped -P ./2_5_2018_RVTDT/phased_data/8q24.recode.8q24.recode.FILTERED.FILTERED2.BEAGLE.phen \
-M ./2_5_2018_RVTDT/phased_data/8q24.recode.8q24.recode.FILTERED.FILTERED2.BEAGLE.map \
--adapt 500 --alpha 0.00001 --permut 2000 \
--lower_cutoff 0 --upper_cutoff 100 \
--minVariants 3 \
--maxMissRatio 1

#Ran but still didn't work
