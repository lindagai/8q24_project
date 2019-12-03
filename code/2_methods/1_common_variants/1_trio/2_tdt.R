#################### 5. Do TDTs using trio  ########################

filepath_geno_matrix_rds<-"/users/lgai/8q24_project/data/processed_data/trio_tdt_format/chr8_geno_matrix_15_individuals_removed_half_calls_removed_maf_0.01.RDS"

chr8_geno_filt<-readRDS(filepath_geno_matrix_rds)

dim(chr8_geno_filt)
#981 3393

#A. Allelic TDT
allelic_results<-allelicTDT(chr8_geno_filt)

#########

#B. Test the 1st SNP with a genotypic TDT
tdt(chr8_geno_filt[,1])

#########

#C. Test interaction between 1st and 2nd SNP with a genotypic TDT
tdtGxG(chr8_geno_filt[,1], chr8_geno_filt[,2])

#########

#D. Testing all SNPs in a Matrix in Genotype Format with a Genotypic TDT
colTDT(chr8_geno_filt)

