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
#Last update: 3/12/18
#TO DO: Clean up code to get rid of cruft

################################################################################



############ B. Comparison of EIGEN, GWAVA, CADD13 scores ###########

#Read in the annovar report with updated CADD

#filepath_annovar_report_updated_CADD<-file.path("/users","lgai","8q24_project","data","processed_data","annotation","8q24_annovar_reports_CADD13.txt")
#filepath_annovar_report_updated_CADD<-file.path("/users","lgai","8q24_project","data","processed_data","annotation","8q24_annovar_reports_CADD13.txt")

filepath_annovar_report_updated_CADD<-"/Users/lindagai/Documents/classes/3rd year/3rd term/Margaret:Ingo/3:13 meeting/8q24_annovar_reports_CADD13.txt"

annovar_report<-read.table(filepath_annovar_report_updated_CADD,sep="\t",quote ="",header=TRUE)
colnames(annovar_report)
head(annovar_report)
dim(annovar_report) #24165    12

#Missingness
devtools::install_github('njtierney/neato')
library(neato)
#dependency ‘ggalt’ is not available for package ‘neato’
#install.packages("ggalt", repos='http://cran.us.r-project.org')

#install.packages("Amelia", repos='http://cran.us.r-project.org',dependencies=TRUE)
#library(Amelia)
#library(ggalt)

ggplot_missing(annovar_report[,10:12])

###############################################

#Filter to rows WITH missing values in any row


################ Filter to rows WITHOUT missing values in any row ##################

annovar_report_no_missing<-annovar_report[(which(!is.na(annovar_report$CADD13) & !is.na(annovar_report$WG_GWAVA_score)) ),]
dim(annovar_report_no_missing) #5074   12
dim(annovar_report) #24165    12

annovar_report[which(is.na(annovar_report$WG_EIGEN_score))]

#Don't need to check EIGEN score bc EIGEN available for all rows

################

#Score cut-offs: Literature

#How much overlap is there when you cut off by the 

#Score cut-offs: 90th percentile
#Score cut-offs: 90th percentile
eigen_90_quantile <- quantile(annovar_report_no_missing$WG_EIGEN_score, probs = 0.9)[[1]]
gwava_90_quantile <- quantile(annovar_report_no_missing$WG_GWAVA_score, probs = 0.9)[[1]]
cadd_90_quantile <- quantile(annovar_report_no_missing$CADD13, probs = 0.9)[[1]]

#90th percentile CADD => 90th percentile GWAVA?
GWAVA_90th<-rep(NA,nrow(annovar_report_no_missing))

for (i in 9:nrow(annovar_report_no_missing)){
  if (annovar_report_no_missing$WG_GWAVA_score[i]> gwava_90_quantile){
    GWAVA_90th[i] = 1
  }# else {
  #  GWAVA_90th[i] = 0
  #}
}

GWAVA_90th

EIGEN_90th<-rep(NA,nrow(annovar_report_no_missing))

for (i in 9:nrow(annovar_report_no_missing)){
  if (annovar_report_no_missing$WG_EIGEN_score[i]> eigen_90_quantile){
    EIGEN_90th[i] = 1
  }
}

CADD13_90th<-rep(NA,nrow(annovar_report_no_missing))

for (i in 9:nrow(annovar_report_no_missing)){
  if (annovar_report_no_missing$CADD13[i]> cadd_90_quantile){
    CADD13_90th[i] = 1
  }
}

annovar_90th_percentile_indicator<-cbind(CADD13_90th,EIGEN_90th,GWAVA_90th)
title="Is the score greater than the 90th percentile in this dataset?"
ggplot_missing(annovar_90th_percentile_indicator)

