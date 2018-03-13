#Description:
# Cleans Annovar report file and obtain CADD1.3 scores for 8q24 region
# Exploratory data analysis of annotation information
#
#
# Sections:
# 1. Clean raw files
#         A. Annovar Report
# 2. Exploratory data analysis of Annovar Report
#         A. Degree of missingness
#         B. Graphs of CADDgt20, EIGEN, GWAVA scores
# 3. Additional CADD scores from CADD1.3 browser
#         A. Downloading scores from CADD browser
#         B. Add CADD1.3 scores to Annovar Report
# 4. Exploratory data analysis of Annovar Report/CADD1.3
#         A. Graphs of CADD1.3 scores
#         B. Comparison of EIGEN, GWAVA, CADD13 scores
#
# Input:
#
#
# Output:
#
#

#Author: Linda Gai
#Last update: 3/11/18
#TO DO: Clean up code to get rid of cruft

################################################################################

ssh -Y lgai@jhpce01.jhsph.edu

#Allocate memory
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G
module load R
R

################# 1. Clean raw files #########################

######## Cleaning Annovar report #########

######## Subset chr8 annotation file to 8q24 region #########

#Original location
#filepath_annovar_report<-"/dcl01/beaty/data/gmkf/euro/anno/ANNOVAR_Uploaded/fullannot.gmkf.snvs.c8.T3.ss_ANNOVAR_REPORT.txt"

#Copied location
filepath_annovar_report<-"/users/lgai/8q24_project/data/raw_data/fullannot.gmkf.snvs.c8.T3.ss_ANNOVAR_REPORT.txt"

annovar_report<-read.table(filepath_annovar_report, quote = "", sep = "\t",header=TRUE)

#annovar_report_columns<-colnames(annovar_report)
#saved in annovar contents.txt
#annovar_report_columns[c(1:3,5,85,86,88,90)]

#dim(annovar_report)# 3705599     102
#annovar_report[1:5,1:5]

library(VariantAnnotation)

#TODO: finish moving raw VCF to 8q24_project/data/raw_data and repeat cleaning steps
vcf<-readVcf("/users/lgai/2_5_2018_RVTDT/8q24.recode.FILTERED.FILTERED2.BEAGLE.vcf.gz","hg19")
vcf_geno<-geno(vcf)

#dim(vcf)
#length(vcf_geno)
#vcf_geno[1:3]
#vcf[1:5,1:5]

#Get the part of annoVar report that is relevant to the 8q24 region
#gwava_90_quantile <- quantile(GWAVA_score_mat$WG_GWAVA_score, probs = 0.9)[[1]] #90th quantile = 0.44

starts<-start(rowRanges(vcf))
head(starts)

ends<-end(rowRanges(vcf))

start <-min(starts)
end <-max(ends)

#NOTE: the start and end rows are the same here.
#setdiff(starts,ends)

######## Obtain quantiles of start/end positions #########

#Get quantile of BP positions to split vcf into ~2MB files.
#You will upload the 2 MB to CADD browsers to get CADD scores.
#TO DO: Write how to download full CADD dataset and subset to your region(s)

quantile(starts, probs = 0.33) #129642137
quantile(starts, probs = 0.66) #129970267

quantile(ends, probs = 0.33) #129642137
quantile(ends, probs = 0.66) #129970267

#################

#Get the positions between start and end in AnnoVar
annovar_start_rows<-annovar_report[which(annovar_report$StartPosition > start),]
#dim(annovar_start_rows)
#annovar_start_rows$StartPosition[1:5]

annovar_8q24<-annovar_start_rows[which(annovar_start_rows$EndPosition < end),]
#dim(annovar_8q24)

filepath_annovar_8q24<-file.path("/users","lgai","8q24_project","data","processed_data","annotation","fullannot.gmkf.snvs.8q24.T3.ss_ANNOVAR_REPORT.txt")

write.table(annovar_8q24, filepath_annovar_8q24, sep="\t",row.names = FALSE,quote = FALSE)

######## Subset 8q24 annotation to relevant columns #########

filepath_annovar_8q24<-file.path("/users","lgai","8q24_project","data","processed_data","annotation","fullannot.gmkf.snvs.8q24.T3.ss_ANNOVAR_REPORT.txt")

annovar_8q24<-read.table(filepath_annovar_8q24,sep="\t",quote ="",header=TRUE)

#Check
#head(annovar_8q24)
#dim(annovar_8q24) #24166   102

#Choose the relevant columns
#cols<-colnames(annovar_8q24)
#cols[c(2:5,56:58,45,77,93:94)]

library(dplyr)

#sm_annovar_8q24<-annovar_8q24[,c(2:5,56:58,45,77,93:94)]

sm_annovar_8q24<-annovar_8q24 %>%
  select("StartPosition",
         "EndPosition",
         "ReferenceAllele",
          "AlternativeAllele",
         "Score_Ljb_phylop",
         "Score_Ljb_pp2hdiv",
        "Score_Ljb_pp2hvar",
         "SiftScore",
         "CADDgt20",
        "WG_GWAVA_score",
        "WG_EIGEN_score")

sm_filepath_annovar_8q24<-file.path("/users","lgai","8q24_project","data","processed_data","annotation","8q24_annovar_report.txt")
sm_filepath_annovar_8q24

write.table(sm_annovar_8q24, sm_filepath_annovar_8q24, sep="\t",row.names = FALSE,quote = FALSE)

#Check
annovar_report<-read.table(sm_filepath_annovar_8q24,sep="\t",header=TRUE, quote ="")
head(annovar_report)
################# # Part B. Exploratory Data Analysis

################# Degree of missingness in each variable ########################

sm_filepath_annovar_8q24<-file.path("/users","lgai","8q24_project","data","processed_data","annotation","8q24_annovar_report.txt")

annovar_report<-read.table(sm_filepath_annovar_8q24,sep="\t",header=TRUE, quote ="")

dim(annovar_report)
head(annovar_report)

############# Non-missing variables ###############

#CADD score - initial filtering

#Select rows with non-missing CADD
caddscore_mat<-annovar_report[which(!is.na(annovar_report$CADDgt20)),]
head(caddscore_mat)
dim(caddscore_mat) #124  11

#Select rows with CADD >5
caddscore_mat_5<-caddscore_mat[caddscore_mat$CADDgt20>5,]
head(caddscore_mat_5)
dim(caddscore_mat_5) # 6 11

#No rows with CADD > 10 

################### GWAVA_score #######################

GWAVA_score_mat<-annovar_report[which(!is.na(annovar_report$WG_GWAVA_score)),]
head(GWAVA_score_mat)
dim(GWAVA_score_mat) #8899   11

range(GWAVA_score_mat$WG_GWAVA_score) #0.09 0.76

#Filter by all with GWAVA_score>0.5, e.g.
GWAVA_score_mat_0.5<-GWAVA_score_mat[GWAVA_score_mat$WG_GWAVA_score>0.5,]
head(GWAVA_score_mat_0.5)
dim(GWAVA_score_mat_0.5) #279  11

################### EIGEN_score #######################

EIGEN_score_mat<-annovar_report[which(!is.na(annovar_report$WG_EIGEN_score)),]
head(EIGEN_score_mat)
dim(EIGEN_score_mat) #24165    11

range(EIGEN_score_mat$WG_EIGEN_score) #-1.4515  2.3631

#Filter by all with EIGEN_score>0
EIGEN_score_mat_0<-EIGEN_score_mat[EIGEN_score_mat$WG_EIGEN_score>0,]
head(EIGEN_score_mat_0)
dim(EIGEN_score_mat_0) #6048   11

############# Missing variables ###############

#Sift score
siftscore_mat<-annovar_report[which(!is.na(annovar_report$SiftScore)),]
dim(siftscore_mat)
#all NA

#Phylop_score
phylop_score_mat<-annovar_report[which(!is.na(annovar_report$Score_Ljb_phylop)),]
head(phylop_score_mat)
dim(phylop_score_mat)
#all NA

#Pp2hdiv_score
pp2hdiv_score_mat<-annovar_report[which(!is.na(annovar_report$Score_Ljb_pp2hdiv)),]
head(pp2hdiv_score_mat)
#all NA

#Pp2hvar_score for now
pp2hvar_score_mat<-annovar_report[which(!is.na(annovar_report$Score_Ljb_pp2hvar)),]
head(pp2hvar_score_mat)
#all NA

################### Graphs #######################

#Graph the GWAVA

hist_gwava_filepath<-file.path("/users","lgai","8q24_project","figures","exploratory_figures","gwava_hist.pdf"))

pdf(file=hist_gwava_filepath)
gwava_90_quantile <- quantile(GWAVA_score_mat$WG_GWAVA_score, probs = 0.9)[[1]] #90th quantile = 0.44

hist_gwava<-hist(as.numeric(GWAVA_score_mat$WG_GWAVA_score),main="WG_GWAVA_score frequency",
                 xlab="GWAVA_score",ylab="Frequency",xlim=c(0,0.8))
abline(v=gwava_90_quantile,col="red")

#-1.4515  2.3631

dev.off()


#Graph the EIGEN

hist_eigen_filepath<-file.path("/users","lgai","8q24_project","figures","exploratory_figures","eigen_hist.pdf")

pdf(file=hist_eigen_filepath)
eigen_90_quantile <- quantile(EIGEN_score_mat$WG_EIGEN_score, probs = 0.9)[[1]] #90th quantile = 0.23

hist_eigen<-hist(as.numeric(EIGEN_score_mat$WG_EIGEN_score),main="WG_EIGEN_score frequency",
                 xlab="EIGEN_score",ylab="Frequency",xlim=c(-1.5,2.5))
abline(v=eigen_90_quantile,col="red")

dev.off()

################### 3. Additional CADD scores from CADD browser ########################

############ A. Downloading scores from CADD browser ###########

#Break vcf into thirds for cadd scores to make them <2 MB for uploading

#TODO: download CADD score file for whole genome and use that method instead

#TODO: move raw vcf files into /users/lgai/8q24_project/figures and re-clean/re-name
module load vcftools

vcftools --gzvcf 8q24.recode.FILTERED.FILTERED2.BEAGLE.vcf.gz --chr 8 --to-bp 129642137 --recode --out 8q24.recode.FILTERED.FILTERED2.BEAGLE.to.bp129642137
vcftools --gzvcf 8q24.recode.FILTERED.FILTERED2.BEAGLE.vcf.gz --chr 8 --from-bp 129642137 --recode --to-bp 129970267 --out 8q24.recode.FILTERED.FILTERED2.BEAGLE.from.bp129642137.to.bp129970267
vcftools --gzvcf 8q24.recode.FILTERED.FILTERED2.BEAGLE.vcf.gz --chr 8 --from-bp 129970267 --recode --out 8q24.recode.FILTERED.FILTERED2.BEAGLE.from.bp129970267

#Upload the vcf files to CADD browser to get the CADD scores

#Read in the files

filepath_cadd_scores1<-file.path("/users","lgai","8q24_project", "data","raw_data","cadd.8q24.from.beginning.to.bp129642137.tsv")

cadd_scores1<-read.table(filepath_cadd_scores1,sep="\t",quote ="")
head(cadd_scores1)
dim(cadd_scores1) #4236    6

#from.to.tsv
filepath_cadd_scores2<-file.path("/users","lgai","8q24_project", "data","raw_data","cadd.8q24.from.bp.129642137.to.bp.129970267.tsv")

cadd_scores2<-read.table(filepath_cadd_scores2,sep="\t",quote ="")

head(cadd_scores2)
dim(cadd_scores2) #4365    6

#to.bp129642137.tsv
filepath_cadd_scores3<-file.path("/users","lgai","8q24_project", "data","raw_data","cadd.8q24.from.bp.12977026.to.end.tsv")

cadd_scores3<-read.table(filepath_cadd_scores3,sep="\t",quote ="")
head(cadd_scores3)
tail(cadd_scores3)
dim(cadd_scores3) #4236    6

#Get the list of CADD scores for every position

cadd_8q24<-rbind(cadd_scores1,cadd_scores2,cadd_scores3)
dim(cadd_8q24) #12837     6

#Get base pair position and scaled cadd score
cadd_8q24_reduced<-cadd_8q24[,c(2,6)]
head(cadd_8q24_reduced)

colnames(cadd_8q24_reduced)<-(c("Position","CADD13"))

#Save
filepath_cadd_8q24<-file.path("/users","lgai","8q24_project", "data","processed_data","annotation", "cadd_scores_8q24.tsv")
                              
write.table(cadd_8q24_reduced, filepath_cadd_8q24, sep="\t",row.names = FALSE,quote = FALSE)



################ B. Add CADD1.3 scores to Annovar Report ###########################

#Read in

filepath_cadd_8q24<-file.path("/users","lgai","8q24_project", "data","processed_data","annotation", "cadd_scores_8q24.tsv")

cadd_8q24<-read.table(filepath_cadd_8q24,sep="\t",quote ="",header=TRUE)
head(cadd_8q24)
#colnames(cadd_8q24)

#Read 
filepath_annovar_8q24<-file.path("/users","lgai","8q24_project","data","processed_data","annotation","8q24_annovar_report.txt")
annovar_report<-read.table(filepath_annovar_8q24,sep="\t",header=TRUE, quote ="")
head(annovar_report)

dim(annovar_report) #24164    11

#In the cleaned version, the first 2 SNP positions are duplicated but have diff. rows...
#24165    11

colnames(annovar_report)

#Add CADD to appropriate rows
library(dplyr)
annovar_report_updated_CADD<-left_join(annovar_report,cadd_8q24,by=c("StartPosition"="Position"))
head(annovar_report_updated_CADD)

#Save again
filepath_annovar_report_updated_CADD<-file.path("/users","lgai","8q24_project","data","processed_data","annotation","8q24_annovar_reports_CADD13.txt")

write.table(annovar_report_updated_CADD, filepath_annovar_report_updated_CADD, sep="\t",row.names = FALSE,quote = FALSE)

################## A. CADD Exploratory Data Analysis #########################

#Read in the annovar report with updated CADD

filepath_annovar_report_updated_CADD<-file.path("/users","lgai","8q24_project","data","processed_data","annotation","8q24_annovar_reports_CADD13.txt")

annovar_report<-read.table(filepath_annovar_report_updated_CADD,sep="\t",quote ="",header=TRUE)

head(annovar_report)
dim(annovar_report) #24165    12

#Remove missing rows
caddscore_mat<-annovar_report[which(!is.na(annovar_report$CADD13)),]
head(caddscore_mat)
dim(caddscore_mat) #12842    12

#Filter by all with caddscore>5
caddscore_mat_5<-caddscore_mat[caddscore_mat$CADD13>5,]
head(caddscore_mat_5)
dim(caddscore_mat_5) # 2613   12

#Filter by all with caddscore>10
caddscore_mat_10<-caddscore_mat[caddscore_mat$CADD13>10,]
head(caddscore_mat_10)
dim(caddscore_mat_10) # 812  12

#Filter by all with caddscore>15
caddscore_mat_15<-caddscore_mat[caddscore_mat$CADD13>15,]
head(caddscore_mat_15)
dim(caddscore_mat_15) #243  12

#################### 4. Exploratory data analysis of Annovar Report/CADD1.3 ########################

############ A. Graphs of CADD1.3 scores ###########

#Graph the cadd

hist_cadd_filepath<-file.path("/users","lgai","8q24_project","figures","exploratory_figures","cadd_hist.pdf")

pdf(file=hist_cadd_filepath)
cadd_90_quantile <- quantile(caddscore_mat$CADD13, probs = 0.9)[[1]] #90th quantile = 7.7207

hist_cadd<-hist(as.numeric(caddscore_mat$CADD13),main="CADDv1.3 score frequency",
                xlab="cadd_score",ylab="Frequency",xlim=c(0,20))
abline(v=cadd_90_quantile,col="red")

dev.off()

############ B. Comparison of EIGEN, GWAVA, CADD13 scoress ###########

#Missingness

#Score cut-offs: Literature

#Score cut-offs: 90th percentile