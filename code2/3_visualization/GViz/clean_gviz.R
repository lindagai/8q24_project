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
                        start = from, end = to, name = "Genes", #showId=TRUE,
                        geneSymbols=TRUE, collapseTranscripts="meta", 
                        col.title="black", cex.title=1.5, cex.group=1.5)
                        #col.title="black", cex.title=0.5, cex.group=1.5)

################### 2. SymbolMapping track
data(genesymbol, package = "biovizBase")
granges.8q24 <- subsetByOverlaps(genesymbol, GRanges(seqnames="chr8", ranges=IRanges(start=129295902, end=130354740)))

track.8q24<-GeneRegionTrack(granges.8q24,
                #             genome = gen, chromosome = chr,
                # start = from, end = to, name = "Genes", showId=TRUE,
                geneSymbols=TRUE, collapseTranscripts="meta", 
                col.title="black", cex.title=1.5, cex.group=1.5)

################### 3.

#A. gTDT track
filepath_genotypic_tdt_results<-"/users/lgai/8q24_project/results/8q24_genotypic_tdt_results.txt"
genotypic_tdt_results<-read.table(filepath_genotypic_tdt_results,header=TRUE)
head(genotypic_tdt_results)
log10p.gTDT <- -log10(genotypic_tdt_results$pval)

dtrack.gTDT <- DataTrack(data=log10p.gTDT, start=genotypic_tdt_results$position-1, end=genotypic_tdt_results$position,
                        genome="hg19", chromosome=chr, name= "gTDT -log10p",
                        ylim=c(-0.5,8.5),baseline=0,v=0,col.line="grey92", cex=2,
                        cex.title=1.5, col.title="black", 
                        col.axis="black")

#B. aTDT track

filepath_allelic_tdt_results<-"/users/lgai/8q24_project/results/8q24_allelic_tdt_results.txt"
allelic_tdt_results<-read.table(filepath_allelic_tdt_results,header=TRUE)
head(allelic_tdt_results)
log10p.aTDT <- -log10(allelic_tdt_results$pval)

dtrack.aTDT <- DataTrack(data=log10p.gTDT, start=allelic_tdt_results$position-1, end=allelic_tdt_results$position,
                        genome="hg19", chromosome=chr, name= "aTDT -log10p",
                        ylim=c(-0.5,8.5),baseline=0,v=0,col.line="grey92", cex=2,
                        cex.title=1.5, col.title="black", 
                        col.axis="black")

################### 4. Transmitted RV counts track

filepath.trans.rv.ct<-"/users/lgai/GViz_test/8q24_trans_rv_ct.txt"

#Read in list of counts
trans.rv.ct<-read.table(filepath.trans.rv.ct,header=TRUE)
head(trans.rv.ct)
setDT(trans.rv.ct,keep.rownames=TRUE)
head(trans.rv.ct)
dim(trans.rv.ct)

#Join the trans.rv.ct with the position, and the SNPs in the data track of previous parts
filepath_genotypic_tdt_results<-"/users/lgai/8q24_project/results/8q24_genotypic_tdt_results.txt"
genotypic_tdt_results<-read.table(filepath_genotypic_tdt_results,header=TRUE)
head(genotypic_tdt_results)

df<-left_join(trans.rv.ct,genotypic_tdt_results,by=c("rn"="snp_name"))
head(df)
# df

RVtrack<-DataTrack(data=df$trans.ct,
                   start=df$position-1,
                   end=df$position,
                   genome="hg19",
                   chromosome=chr, name="Transmitted RV Counts\nMAF<0.01 ",
                   type=c("p", "g"),
                   ylim=c(0, 6),  col="blue",
                   baseline=0, v=0, col.line="grey92",
                   cex=2, cex.title=1,
                   col.title="black", col.axis="black")

# RVtrack<-DataTrack(data = df$trans.ct,
#                    start=df$position-1,
#                    end=df$position,
#                    genome="hg19",
#                    chromosome=chr, name="Transmitted RV Counts",
#                    type=c("p", "g"),
#                    ylim=c(0, 6), col=c("blue"),
#                    baseline=0, v=0, col.line="grey92", cex=2, cex.title=1.5,
#                    col.title="black", col.axis="black")

################### 5. ScanTrio LR track

filepath.scantrio<-"/users/lgai/GViz_test/scantrio.pos.and.lr.txt"
df<-read.table(filepath.scantrio,header=TRUE)
head(df)

df2<-left_join(df,genotypic_tdt_results,by=c("avg.pos"="position"))
head(df2)

ScanTrioLRtrack<-DataTrack(data = df2$loglr,
                           start=df2$avg.pos-1,
                           end=as.numeric(df2$avg.pos-1),
                           genome="hg19",
                           chromosome=chr, name="logLR ScanTrio, \n win=100M",
                           type="l",
                           ylim=c(0, 35), col=c("blue"),
                           # baseline=0, v=0, col.line="blue", cex=2, cex.title=1.5,
                           baseline=0, v=0, col.line="blue", cex=2, cex.title=1,
                           col.title="black", col.axis="black")

#id="name", name = "Visel", rot.title=0, cex.title=1, 

################### Create graph

pdf("/users/lgai/GViz_test/my_gviz_scanTrio.pdf",height=9,width=16)

plotTracks(list(itrack,
                gtrack,
                track.8q24,
                dtrack.gTDT,
                dtrack.aTDT,
                RVtrack,
                ScanTrioLRtrack#,
                #txTr
),
background.title="darkgray",from= from, to=to)

dev.off()

