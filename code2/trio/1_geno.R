################# 4. Create genotype matrix for trio  ########################

#Read in ped file
filepath_ped<-"/users/lgai/8q24_project/data/processed_data/gmkf_euro_completetrios_ids_match_vcf_mend_errors_removed.txt"

chr8_ped<-read.table(filepath_ped,header=TRUE)
head(chr8_ped)
dim(chr8_ped)

#Remove half calls from the vcf

#Read in vcf file; takes a minute
filepath_vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/8q24.recode.15.mend.error.individuals.removed.recode.vcf"
chr8_vcf <- readVcf(filepath_vcf,"hg19")

#Delete "half" missing values, where one of the alleles is called but other isn't
geno(chr8_vcf)$GT[geno(chr8_vcf)$GT == "0/."]<-"."
geno(chr8_vcf)$GT[geno(chr8_vcf)$GT == "1/."]<-"."

#Create geno matrix
chr8_geno <- vcf2geno(chr8_vcf,chr8_ped)

#Filter to remove
#perc.na: percentage of missing values that a SNP or a subject is allowed to have

chr8_geno_filt<-removeSNPs(chr8_geno, maf=0.01, perc.na =0.1)
dim(chr8_geno_filt)
#981 3393

filepath_geno_matrix_rds<-"/users/lgai/8q24_project/data/processed_data/trio_tdt_format/chr8_geno_matrix_15_individuals_removed_half_calls_removed_maf_0.01.RDS"

saveRDS(chr8_geno_filt,filepath_geno_matrix_rds)
