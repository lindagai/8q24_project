################################################################################

#Sign into cluster
ssh -X lgai@jhpce01.jhsph.edu 

#Allocate memory
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G
module load R
R

################################################################################

#TODO: 
# Resize things so it looks better
# Create LR graph for Scan Trio

#load package
install.packages("GViz")
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("Gviz")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")

library(Gviz)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(dplyr)
library(data.table)

################### 1. Basic information and ideogram and gene region tracks

chr <- "chr8"
#TODO: Positions are taken from VCF, may need editing
from <- 129295902
to <- 130354740
regCode <- "8q24"
gen <- "hg19"

#Set up initial tracks
itrack <- IdeogramTrack(genome = gen, chromosome = chr, showBandId = TRUE,
                        cex.bands=0.5#, 
                        #showID=TRUE
)
gtrack <- GenomeAxisTrack()
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

#Get the gene region track for the database results in the relevant region
txTr <- GeneRegionTrack(txdb, genome = gen, chromosome = chr,
                        start = from, end = to, 
                        name = "Genes", #showId=TRUE,
                        geneSymbols=TRUE, collapseTranscripts="meta", 
                        col.title="black", cex.title=1.5, cex.group=1.5)
#col.title="black", cex.title=0.5, cex.group=1.5)

################### 2. SymbolMapping track
data(genesymbol, package = "biovizBase")
granges.8q24 <- subsetByOverlaps(genesymbol, GRanges(seqnames="chr8", ranges=IRanges(start=129295902, end=130354740)))

track.8q24<-GeneRegionTrack(granges.8q24,
                            name = "Genes",
                            genome = gen, chromosome = chr,
                            start = from, end = to, showId=TRUE,
                            geneSymbols=TRUE, 
                            collapseTranscripts="meta", 
                            col.title="black", cex.title=1.5, cex.group=1.5)

################### 3. gTDT track

#A. gTDT track
#filepath_genotypic_tdt_results<-"/users/lgai/8q24_project/results/8q24_genotypic_tdt_results.txt"
filepath_genotypic_tdt_results<-"/Users/lindagai 1/Documents/classes/4th year/Research/8q24_project/Update 12:13/8q24_genotypic_tdt_results.txt"
genotypic_tdt_results<-read.table(filepath_genotypic_tdt_results,header=TRUE)
head(genotypic_tdt_results)
log10p.gTDT <- -log10(genotypic_tdt_results$pval)

dtrack.gTDT <- DataTrack(data=log10p.gTDT, start=genotypic_tdt_results$position-1, end=genotypic_tdt_results$position,
                         genome="hg19", chromosome=chr, name= "-log10p (gTDT)",
                         ylim=c(-0.5,8),baseline=0,v=0,col.line="grey92", cex=2,
                         cex.title=1.5, col.title="black", 
                         col.axis="black")

################### 4. aTDT track
#scp lgai@jhpce01.jhsph.edu:"/users/lgai/8q24_project/results/8q24_allelic_tdt_results.txt" "/Users/lindagai 1/Documents/classes/4th year/Research"

filepath_allelic_tdt_results<-"/Users/lindagai 1/Documents/classes/4th year/Research/8q24_allelic_tdt_results.txt"
allelic_tdt_results<-read.table(filepath_allelic_tdt_results,header=TRUE)
head(allelic_tdt_results)
log10p.aTDT <- -log10(allelic_tdt_results$pval)

dtrack.aTDT <- DataTrack(data=log10p.gTDT, start=allelic_tdt_results$position-1, end=allelic_tdt_results$position,
                         genome="hg19", chromosome=chr, name= "-log10p (aTDT)",
                         ylim=c(-0.5,8.5),baseline=0,v=0,col.line="grey92", cex=2,
                         cex.title=1.5, col.title="black", 
                         col.axis="black")

################### 5. ScanTrio p-val track

#filepath.scantrio<-"/users/lgai/GViz_test/scantrio.pos.and.lr.txt"
#filepath.scantrio<-"/users/lgai/8q24_project/data/processed_data/Scan-Trio/scantrio.pos.and.lr.txt"

#Filepath for windows = 100, common and rare
filepath.scantrio<-'/Users/lindagai 1/Documents/classes/4th year/4th term/Research/scantrio.pos.p.and.lr.txt'

df<-read.table(filepath.scantrio,header=TRUE)
head(df)

df2<-left_join(df,genotypic_tdt_results,by=c("avg.pos"="position"))
head(df2)
range(df$p)
which(df$p==-Inf)
# 
# ScanTrioLRtrack<-DataTrack(data = -log10(df2$p),
#                            start=as.numeric(df2$avg.pos-1),
#                            end=as.numeric(df2$avg.pos-1),
#                            genome="hg19",
#                            chromosome=chr, name="-log10p (ScanTrio)",
#                            type="l",
#                            ylim=c(0, 8),
#                            col=c("blue"),
#                            baseline=-log10(bonferroni.sig.level),
#                            col.baseline = "red",
#                            # baseline=0, v=0, col.line="blue", cex=2, cex.title=1.5,
#                            #baseline=0,
#                            v=0, col.line="blue", cex=2, cex.title=1,
#                            col.title="black", col.axis="black")

n.windows.common<-16
bonferroni.sig.level.common<-0.05/n.windows.common

range(-log10(df$p))
-log10(df$p)[which(-log10(df$p) > 100)]

#NOTE: GViz automatically removes infinite values, so you must replace them
# with a real number

test<--log10(df$p)

test[which(-log10(df$p) > 100)]<-rep(100,length(which(-log10(df$p) > 100)))

test[which(-log10(df$p) > 99)]
test[which(-log10(df$p) > 2)]

#id="name", name = "Visel", rot.title=0, cex.title=1, 

ScanTrioLRtrack2<-DataTrack(data = test,
                            start=as.numeric(df$avg.pos-1),
                            end=as.numeric(df$avg.pos-1),
                            genome="hg19",
                            chromosome=chr, name="-log10p (ScanTrio)",
                            type="l",
                            ylim=c(0, 8),
                            col=c("blue"),
                            baseline=-log10(bonferroni.sig.level.common),
                            col.baseline = "red",
                            # baseline=0, v=0, col.line="blue", cex=2, cex.title=1.5,
                            #baseline=0,
                            v=0, col.line="blue", cex=2, cex.title=1,
                            col.title="black", col.axis="black")

################### Create graph

filepath.gviz<-"/Users/lindagai 1/Documents/classes/4th year/4th term/Research/GViz_commonv2_fixed.pdf"

pdf(filepath.gviz,height=16,width=16)

plotTracks(list(itrack,
                gtrack,
                #track.8q24,
                dtrack.gTDT,
                dtrack.aTDT,
                ScanTrioLRtrack2#,
                #txTr
),sizes=c(1,1,3,3,3),
background.title="darkgray",from= from, to=to)

dev.off()
