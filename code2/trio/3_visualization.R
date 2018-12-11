#################### 6. Create graphs of -10logp vs position  ########################

#A. Get dataframe of SNP names in the genotype matrix and position

#If you need to re-create DF of snp_name and position, comment this code chunk back in

# #Read in list of SNPs and positions
# filepath_vcf<-"/users/lgai/8q24_project/data/raw_data/vcfs/8q24.recode.vcf"
# chr8_vcf <- readVcf(filepath_vcf,"hg19")
# 
# #Get snp positions
# position<-start(rowRanges(chr8_vcf))
# 
# # #Checks
# #vcf_starts<-start(rowRanges(chr8_vcf))
# #vcf_ends<-end(rowRanges(chr8_vcf))
# # head(vcf_starts)
# # head(vcf_ends)
# # setdiff(vcf_starts, vcf_ends)
# #integer(0)
# #No difference between start and end positions
# 
# #Get snp names
# snp_name <- names(rowRanges(chr8_vcf))
# #head(vcf_snps)
# 
# #Create dataframe of vcf_snp_and_pos, save it
# vcf_snp_and_pos<-cbind(snp_name,position)
# #head(vcf_snp_and_pos)
# #dim(vcf_snp_and_pos) #24170
# 
# #Remove duplicated rows
# vcf_snp_and_pos_no_duplicates<-vcf_snp_and_pos[!duplicated(vcf_snp_and_pos),]
# 
# #dim(vcf_snp_and_pos_no_duplicates) #24136
# 
# #Save
# filepath_vcf_snp_and_pos<-"/users/lgai/8q24_project/data/processed_data/vcf_snp_names_and_positions.txt"
# write.table(vcf_snp_and_pos_no_duplicates, filepath_vcf_snp_and_pos, row.names=FALSE,quote=FALSE)

#Read in 
filepath_vcf_snp_and_pos<-"/users/lgai/8q24_project/data/processed_data/vcf_snp_names_and_positions.txt"
vcf_snp_and_pos<-as.data.frame(read.table(filepath_vcf_snp_and_pos,header=TRUE))
head(vcf_snp_and_pos)

#Read in genotype matrix and get list of SNPs
filepath_geno_matrix_rds<-"/users/lgai/8q24_project/data/processed_data/trio_tdt_format/chr8_geno_matrix_15_individuals_removed_half_calls_removed_maf_0.01.RDS"
chr8_geno_filt<-readRDS(filepath_geno_matrix_rds)

dim(chr8_geno_filt)
#981 3393

geno_snp_names<-as.data.frame(colnames(chr8_geno_filt))
colnames(geno_snp_names)<-"snp_name"

head(geno_snp_names)
nrow(geno_snp_names) #3393

#Use left_join to get the positions for every SNP in the genotype matrix
geno_snp_and_pos<-left_join(geno_snp_names,vcf_snp_and_pos)
filepath_geno_snp_and_pos<-"/users/lgai/8q24_project/data/processed_data/geno_snp_names_and_positions.txt"

#dim(geno_snp_and_pos) #3393    2
#head(geno_snp_and_pos)

write.table(geno_snp_and_pos, filepath_geno_snp_and_pos, row.names=FALSE,quote=FALSE)

#B. Run the TDTs and create a dataframe of snp names, position,
#tstatistic, p-val, and -log10p

#Read in
filepath_geno_snp_and_pos<-"/users/lgai/8q24_project/data/processed_data/geno_snp_names_and_positions.txt"
geno_snp_and_pos<-read.table(filepath_geno_snp_and_pos,header=TRUE)
head(geno_snp_and_pos)
dim(geno_snp_and_pos) #3393    2

#i. Allelic TDT
allelic_results<-allelicTDT(chr8_geno_filt)

#Create dataframe of allelic TDT p-values and snp_names
tdt_results<-as.data.frame(cbind(allelic_results$stat,allelic_results$pval))
colnames(tdt_results)<-c("stat","pval")

#head(tdt_results)
#dim(tdt_results)
#3393    2

allelic_tdt_df<-cbind(geno_snp_and_pos,tdt_results)
#head(allelic_tdt_df)

#Sort rows by p-value (ascending)
allelic_tdt_df <- allelic_tdt_df[order(allelic_tdt_df$pval),] 
allelic_tdt_df<-allelic_tdt_df %>%
  mutate(neglogp = -log(pval))
head(allelic_tdt_df)
#range(allelic_tdt_df$pval[!is.na(allelic_tdt_df$pval)])

filepath_allelic_tdt_results<-"/users/lgai/8q24_project/results/8q24_allelic_tdt_results.txt"
write.table(allelic_tdt_df, filepath_allelic_tdt_results, row.names=FALSE,quote=FALSE)

#ii. Genotypic TDT
genotypic_results<-colTDT(chr8_geno_filt)

#Create dataframe of genotypic TDT p-values and snp_names
tdt_results<-as.data.frame(cbind(genotypic_results$stat,genotypic_results$pval))
colnames(tdt_results)<-c("stat","pval")

#head(tdt_results)
#dim(tdt_results)
#13403     2

genotypic_tdt_df<-cbind(geno_snp_and_pos,tdt_results)

#Sort rows by p-value (ascending)
genotypic_tdt_df <- genotypic_tdt_df[order(genotypic_tdt_df$pval),] 
genotypic_tdt_df<-genotypic_tdt_df %>%
  mutate(neglogp = -log(pval))
#head(genotypic_tdt_df)
#range(genotypic_tdt_df$pval[!is.na(genotypic_tdt_df$pval)])

filepath_genotypic_tdt_results<-"/users/lgai/8q24_project/results/8q24_genotypic_tdt_results.txt"
write.table(genotypic_tdt_df, filepath_genotypic_tdt_results, row.names=FALSE,quote=FALSE)

#C. Make graphs of -10logp vs position

#i. Allelic TDT
filepath_allelic_tdt_results<-"/users/lgai/8q24_project/results/8q24_allelic_tdt_results.txt"
allelic_tdt_results<-read.table(filepath_allelic_tdt_results,header=TRUE)
head(allelic_tdt_results)

allelic_tdt_plot_filepath<-file.path("/users","lgai","8q24_project","figures","exploratory_figures","8q24_allelic_tdt_neg10logp_position.pdf")

pdf(file=allelic_tdt_plot_filepath)

ggplot(allelic_tdt_results)+
  geom_point(aes(x=position,y=neglogp)) + 
  labs(title="Allelic TDT results",
       x="SNP position", y = "-logp")
dev.off()

#ii. Genotypic TDT
filepath_genotypic_tdt_results<-"/users/lgai/8q24_project/results/8q24_genotypic_tdt_results.txt"
genotypic_tdt_results<-read.table(filepath_genotypic_tdt_results,header=TRUE)
head(genotypic_tdt_results)

genotypic_tdt_plot_filepath<-file.path("/users","lgai","8q24_project","figures","exploratory_figures","8q24_genotypic_tdt_neg10logp_position.pdf")

pdf(file=genotypic_tdt_plot_filepath)

ggplot(genotypic_tdt_results)+
  geom_point(aes(x=position,y=neglogp)) + 
  labs(title="Genotypic TDT results",
       x="SNP position", y = "-logp")

#WARNING: Removed 18 rows containing missing values (geom_point).
#TODO: Figure out why there are missing values

which(is.na(genotypic_tdt_results$neglogp))
# [1] 3376 3377 3378 3379 3380 3381 3382 3383 3384 3385 3386 3387 3388 3389 3390
#[16] 3391 3392 3393


dev.off()
