#Get transmitted RV counts

#1. Get the SNPs where the allele frequency <0.05
filepath_geno_matrix_rds<-"/users/lgai/8q24_project/data/processed_data/trio_tdt_format/chr8_geno_matrix_15_individuals_removed_half_calls_removed_maf_0.01.RDS"
geno<-readRDS(filepath_geno_matrix_rds)
geno[1:5,1:5]
table(geno)

# 0       1       2
# 2077149  802860  383410

dim(geno)

#TODO: Fix NAs
mafs<-colMeans(geno,na.rm=TRUE)/2 #divide by 2 bc 2 alleles per person
range(mafs)
length(mafs)

# rare.snps.geno<-geno[,which(mafs<0.05)]
# dim(rare.snps.geno)
# #981 952

rare.snps.geno<-geno[,which(mafs<0.01)]
mafs2<-colMeans(rare.snps.geno,na.rm=TRUE)/2
range(mafs2) #ok
dim(rare.snps.geno)
#981  40

#Every rare SNP was in at least 1 child

#2. Find the SNPs that are present in at least one child
child.seq<-seq(3,981,by=3)
child.seq
#For maf <0.05
# child.geno<-rare.snps.geno[child.seq,]
# dim(child.geno) # 327 952
# ct<-colSums(child.geno,na.rm=TRUE)
# range(ct) #2 35

#For maf<0.01
child.geno<-rare.snps.geno[child.seq,]
dim(child.geno) # 327 40
rare.snps<-colnames(child.geno)
ct<-colSums(child.geno,na.rm=TRUE)
range(ct) #2 6
"rs76300445" 

#Every rare SNP was in at least 1 child.
#For the SNPs of MAF<0.01, each child had at least 2 rare SNPs, up to 6.

#3. For every SNP, check to see how many times it was transmitted
#(i.e. it was in a child and in one parent)

rare.snp.table<-as.data.frame(rep(NA,length(rare.snps)))
rare.snp.table
rownames(rare.snp.table)<-rare.snps
colnames(rare.snp.table)<-"trans.ct"
#rare.snp.table[,1]<-as.character(rare.snps)
head(rare.snp.table)

i<-rare.snps[1]
i
rare.snp.table[i,2]
as.factor(i)

for (i in rare.snps){
  child.with.rv<-names(which(child.geno[,i]>0))
  #get row.names of
  rows.child.with.rv<-which(rownames(rare.snps.geno) %in% child.with.rv)
  child.ct<-rare.snps.geno[rows.child.with.rv,i]
  
  #identify parents and check if they have it
  parent.rows<-cbind(rows.child.with.rv-2,rows.child.with.rv-1)
  trans.rv.counts<-0
  j<-1
  for (j in (1:length(child.ct))){
    curr.parents<-parent.rows[j,]
    rare.snps.geno[curr.parents,i]
    parent.ct<-sum(na.omit(rare.snps.geno[curr.parents,i]))
    print(parent.ct)
    trans.ct.trio<-min(parent.ct,child.ct[j])
    trans.ct.trio
    trans.rv.counts<-trans.rv.counts + trans.ct.trio
  }
  trans.rv.counts
  rare.snp.table[i,]<-trans.rv.counts
}

# "rs76300445"
# 
# head(rare.snp.table)
# 
# range(rare.snp.table) #1 6
# 
# rare.snp.table

filepath.trans.rv.ct<-"/users/lgai/GViz_test/8q24_trans_rv_ct.txt"
write.table(rare.snp.table,filepath.trans.rv.ct,sep="\t",row.names = TRUE,quote = FALSE)

#yaaas it worked
