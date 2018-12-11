################# 0. Subset chr8 annovar to 8q24 #########################

library(dplyr)

################# 1. Subset 8q24 Annovar to relevant columns #########################

filepath.annovar<-"/users/lgai/8q24_project/data/processed_data/annotation/fullannot.gmkf.snvs.8q24.T3.ss_ANNOVAR_REPORT.txt"

annovar<-read.table(filepath.annovar,sep="\t",quote ="",header=TRUE,stringsAsFactors=FALSE)

#TODO: This could be a function?
sm.annovar<-annovar %>%
  select("StartPosition",
         "EndPosition",
         "ReferenceAllele",
         "AlternativeAllele",
         "Genotype",
         "Quality",
         "TotalDepth",
         "Afalt_1000g2015aug_all",
         "Afalt_1000g2015aug_afr",
         "Afalt_1000g2015aug_amr",
         "Afalt_1000g2015aug_eur",
         "Afalt_1000g2015aug_eas",
         "Afalt_1000g2015aug_sas",
         "gnomAD_genome_freq",
         "Score_Ljb_phylop",
         "Score_Ljb_pp2hdiv",
         "Score_Ljb_pp2hvar",
         "SiftScore",
         "CADDgt20",
         "WG_GWAVA_score",
         "WG_EIGEN_score"
  )

#Trim whitespace
sm.annovar<-as.data.frame(apply(sm.annovar,2,function(x)gsub('\\s+', '',x)))

filepath.sm.annovar<-file.path("/users","lgai","8q24_project","data","processed_data","annotation","8q24_annovar_report_useful_variables.txt")
write.table(sm.annovar, filepath.sm.annovar, sep="\t",row.names = FALSE,quote = FALSE)

# Check
# annovar<-read.table(filepath.sm.annovar,sep="\t",header=TRUE, quote ="")
# head(annovar)
