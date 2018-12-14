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

