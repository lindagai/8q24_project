#ii. Genotypic TDT
scp lgai@jhpce01.jhsph.edu:"/users/lgai/8q24_project/results/8q24_genotypic_tdt_results.txt" "/Users/lindagai 1/Documents/classes/4th year/Research/8q24_project/Update 12:13"

library(dplyr)
library(ggplot2)

filepath_genotypic_tdt_results<-"/Users/lindagai 1/Documents/classes/4th year/Research/8q24_project/Update 12:13/8q24_genotypic_tdt_results.txt"
genotypic_tdt_results<-read.table(filepath_genotypic_tdt_results,header=TRUE)
head(genotypic_tdt_results)

filepath.annovar<-"/Users/lindagai 1/Documents/classes/4th year/Research/8q24_project/Update 12:13/8q24_annovar_report_functional_annnotation.txt"
annovar<-read.table(filepath.annovar,sep="\t",header=TRUE, quote ="")
pos.annotation<-annovar$StartPosition
head(annovar)

genotypic_tdt_results <- genotypic_tdt_results %>%
  mutate(annot.pos = (position %in% pos.annotation)) %>%
  mutate(peak.pos = (position > 129875000 & position < 130000000)) %>%
  mutate(both.pos = (position > 129875000 & position < 130000000) & (position %in% pos.annotation) )

  head(genotypic_tdt_results)




#

genotypic_tdt_plot_filepath_annotation<-"/Users/lindagai 1/Documents/classes/4th year/Research/8q24_project/Update 12:13/8q24_gTDT_highlight_annotation.pdf"

pdf(genotypic_tdt_plot_filepath_annotation)

ggplot(genotypic_tdt_results)+
  geom_point(aes(x=position,y=neglogp,color=annot.pos)) +
  labs(title="Genotypic TDT results, high annotation score positions labeled green",
       x="SNP position", y = "-logp")

dev.off()

#

genotypic_tdt_plot_filepath_peak<-"/Users/lindagai 1/Documents/classes/4th year/Research/8q24_project/Update 12:13/8q24_gTDT_highlight_peak.pdf"

pdf(genotypic_tdt_plot_filepath_peak)

ggplot(genotypic_tdt_results)+
  geom_point(aes(x=position,y=neglogp,color=peak.pos)) +
  labs(title="Genotypic TDT results, peak positions labeled green",
       x="SNP position", y = "-logp")

dev.off()

#

genotypic_tdt_plot_filepath_both<-"/Users/lindagai 1/Documents/classes/4th year/Research/8q24_project/Update 12:13/8q24_gTDT_highlight_both.pdf"

pdf(genotypic_tdt_plot_filepath_both)

ggplot(genotypic_tdt_results)+
  geom_point(aes(x=position,y=neglogp,color=both.pos)) +
  labs(title="Genotypic TDT results, high annotation score positions in peak labeled green",
       x="SNP position", y = "-logp")

dev.off()

#####################

filepath.no.overlap.map.6<-"/Users/lindagai 1/Documents/classes/4th year/Research/8q24_project/Update 12:13/cadd_10/window100_10_overlap/cadd_10_window-size=100_snps_6.map"
filepath.no.overlap.map.7<-"/Users/lindagai 1/Documents/classes/4th year/Research/8q24_project/Update 12:13/cadd_10/window100_10_overlap/cadd_10_window-size=100_snps_6.map"


#####################

filepath.10.overlap.map.6<-"/Users/lindagai 1/Documents/classes/4th year/Research/8q24_project/Update 12:13/cadd_10/window100_10_overlap/cadd_10_window-size=100_snps_6.map"
overlap.map.6<-read.table(filepath.10.overlap.map.6,header=FALSE,stringsAsFactors = FALSE)
head(overlap.map.6)
filepath.10.overlap.map.7<-"/Users/lindagai 1/Documents/classes/4th year/Research/8q24_project/Update 12:13/cadd_10/window100_10_overlap/cadd_10_window-size=100_snps_7.map"
overlap.map.7<-read.table(filepath.10.overlap.map.7,header=FALSE,stringsAsFactors = FALSE)
overlap.map.7$V2
overlap.pos<-unique(c(overlap.map.6$V2,overlap.map.7$V2))
overlap.pos

genotypic_tdt_results <- genotypic_tdt_results %>%
  mutate(overlap.pos = (snp_name %in% overlap.pos))
unique(genotypic_tdt_results$overlap.pos)

ggplot(genotypic_tdt_results)+
  geom_point(aes(x=position,y=neglogp,color=overlap.pos)) +
  labs(title="Genotypic TDT results, positions in low-p-value 10 overlap data in green",
       x="SNP position", y = "-logp")

filepath.no.overlap.map.6<-"/Users/lindagai 1/Documents/classes/4th year/Research/8q24_project/Update 12:13/cadd_10/window100_no_overlap/cadd_10_window-size=100_snps_6.map"
no.overlap.map.6<-read.table(filepath.no.overlap.map.6,header=FALSE,stringsAsFactors = FALSE)
head(no.overlap.map.6)
filepath.no.overlap.map.7<-"/Users/lindagai 1/Documents/classes/4th year/Research/8q24_project/Update 12:13/cadd_10/window100_no_overlap/cadd_10_window-size=100_snps_7.map"
no.overlap.map.7<-read.table(filepath.no.overlap.map.7,header=FALSE,stringsAsFactors = FALSE)
no.overlap.map.7$V2
no.overlap.pos<-as.data.frame(unique(c(no.overlap.map.6$V2,no.overlap.map.7$V2)))
colnames(no.overlap.pos)<-"snp"
no.overlap.pos

genotypic_tdt_results <- genotypic_tdt_results %>%
  mutate(no.overlap.pos = (snp_name %in% no.overlap.pos))
unique(genotypic_tdt_results$no.overlap.pos)

ggplot(genotypic_tdt_results)+
  geom_point(aes(x=position,y=neglogp,color=overlap.pos)) +
  labs(title="Genotypic TDT results, positions in low-p-value no overlap data in green",
       x="SNP position", y = "-logp")
unique(genotypic_tdt_results$no.overlap.pos)


filepath.10.overlap.map.6<-"/Users/lindagai 1/Documents/classes/4th year/Research/8q24_project/Update
12:13/cadd_10/window100_10_overlap/cadd_10_window-size=100_snps_6.map"
overlap.map.6<-read.table(filepath.10.overlap.map.6,header=FALSE,stringsAsFactors = FALSE)

















filepath.10.overlap.map.6<-"/Users/lindagai 1/Documents/classes/4th year/Research/8q24_project/Update 12:13/cadd_10/window100_10_overlap/cadd_10_window-size=100_snps_6.map"
#overlap.map.6<-read.table(filepath.10.overlap.map.6,header=FALSE,stringsAsFactors = FALSE)
overlap.map.6<-read.table(filepath.10.overlap.map.6,header=FALSE)

head(overlap.map.6)
filepath.10.overlap.map.7<-"/Users/lindagai 1/Documents/classes/4th year/Research/8q24_project/Update 12:13/cadd_10/window100_10_overlap/cadd_10_window-size=100_snps_7.map"
#overlap.map.7<-read.table(filepath.10.overlap.map.7,header=FALSE,stringsAsFactors = FALSE)
overlap.map.7<-read.table(filepath.10.overlap.map.7,header=FALSE)
#overlap.map.7$V2
overlap.pos<-unique(c(overlap.map.6$V2,overlap.map.7$V2))
overlap.pos

genotypic_tdt_results <- genotypic_tdt_results %>%
  mutate(overlap.pos = (snp_name %in% overlap.pos))
unique(genotypic_tdt_results$overlap.pos)

head(genotypic_tdt_results)
typeof(genotypic_tdt_results$snp_name)

ggplot(genotypic_tdt_results)+
  geom_point(aes(x=position,y=neglogp,color=overlap.pos)) +
  labs(title="Genotypic TDT results, positions in low-p-value 10 overlap data in green",
       x="SNP position", y = "-logp")



filepath.10.overlap.map.6<-"/Users/lindagai 1/Documents/classes/4th year/Research/8q24_project/Update 12:13/cadd_10/window100_10_overlap/cadd_10_window-size=100_snps_6.map"
overlap.map.6<-read.table(filepath.10.overlap.map.6,header=FALSE,stringsAsFactors = FALSE)

head(overlap.map.6)
filepath.10.overlap.map.7<-"/Users/lindagai 1/Documents/classes/4th year/Research/8q24_project/Update 12:13/cadd_10/window100_10_overlap/cadd_10_window-size=100_snps_7.map"
overlap.map.7<-read.table(filepath.10.overlap.map.7,header=FALSE,stringsAsFactors = FALSE)
#overlap.map.7$V2
overlap.pos<-unique(c(overlap.map.6$V2,overlap.map.7$V2))
overlap.pos
list<-as.character(overlap.pos$snp)

genotypic_tdt_results <- genotypic_tdt_results %>%
  mutate(overlap.pos = (snp_name %in% overlap.pos))

ggplot(genotypic_tdt_results)+
  geom_point(aes(x=position,y=neglogp,color=overlap.pos)) +
  labs(title="Genotypic TDT results, positions in low-p-value 10 overlap data in green",
       x="SNP position", y = "-logp")

filepath.no.overlap.map.6<-"/Users/lindagai 1/Documents/classes/4th year/Research/8q24_project/Update 12:13/cadd_10/window100_no_overlap/cadd_10_window-size=100_snps_6.map"
no.overlap.map.6<-read.table(filepath.no.overlap.map.6,header=FALSE,stringsAsFactors = FALSE)
head(no.overlap.map.6)
filepath.no.overlap.map.7<-"/Users/lindagai 1/Documents/classes/4th year/Research/8q24_project/Update 12:13/cadd_10/window100_no_overlap/cadd_10_window-size=100_snps_7.map"
no.overlap.map.7<-read.table(filepath.no.overlap.map.7,header=FALSE,stringsAsFactors = FALSE)
no.overlap.map.7$V2
no.overlap.pos<-unique(c(no.overlap.map.6$V2,no.overlap.map.7$V2))
no.overlap.pos

genotypic_tdt_results <- genotypic_tdt_results %>%
  mutate(no.overlap.pos = (snp_name %in% no.overlap.pos))

ggplot(genotypic_tdt_results)+
  geom_point(aes(x=position,y=neglogp,color=no.overlap.pos)) +
  labs(title="Genotypic TDT results, positions in low-p-value no overlap data in green",
       x="SNP position", y = "-logp")