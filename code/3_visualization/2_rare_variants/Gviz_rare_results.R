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

################### 4. Transmitted RV counts track

filepath.trans.rv.ct<-'/Users/lindagai 1/Documents/classes/4th year/4th term/Research/8q24_trans_rv_ct.txt'

#Read in list of counts
trans.rv.ct<-read.table(filepath.trans.rv.ct,header=TRUE)
head(trans.rv.ct)
setDT(trans.rv.ct,keep.rownames=TRUE)
head(trans.rv.ct)
dim(trans.rv.ct)

#Join the trans.rv.ct with the position, and the SNPs in the data track of previous parts
#filepath_genotypic_tdt_results<-"/users/lgai/8q24_project/results/8q24_genotypic_tdt_results.txt"
#genotypic_tdt_results<-read.table(filepath_genotypic_tdt_results,header=TRUE)
head(genotypic_tdt_results)
nrow(genotypic_tdt_results)

df<-left_join(trans.rv.ct,genotypic_tdt_results,by=c("rn"="snp_name"))
head(df)
# df

RVtrack<-DataTrack(data=df$trans.ct,
                   start=df$position-1,
                   end=df$position,
                   genome="hg19",
                   chromosome=chr, name="# of transmitted rare variants ",
                   type=c("p", "g"),
                   ylim=c(0, 6),  col="blue",
                   baseline=0, v=0, col.line="grey92",
                   cex=2, cex.title=1,
                   col.title="black", col.axis="black")

################### 5. RV-TDT track
#filepath.RV_TDT<-"/Users/lindagai 1/Documents/classes/4th year/4th term/Research/results.with.pos.txt"

filepath.RV_TDT<-"/Users/lindagai 1/Documents/classes/4th year/4th term/results.from.all.windows.txt"
df.RV_TDT<-read.table(filepath.RV_TDT,header=TRUE)
head(df.RV_TDT)
dim(df.RV_TDT)

n.windows<-nrow(df.RV_TDT)
n.windows

bonferroni.sig.level<-0.05/n.windows
bonferroni.sig.level

df.RV_TDT<-left_join(df.RV_TDT,genotypic_tdt_results,by=c("Position"="position"))
head(df.RV_TDT)
df.RV_TDT$CMC.Haplo

-log10(bonferroni.sig.level)

RV_TDTtrack<-DataTrack(data = -log10(df.RV_TDT$CMC.Haplo),
                           start=df.RV_TDT$Position-1,
                           end=df.RV_TDT$Position,
                           genome="hg19",
                           chromosome=chr, name="-log10p (RV-TDT)",
                           type="l",
                           ylim=c(0, 8),
                       baseline=-log10(bonferroni.sig.level),
                       col.baseline = "red",
                       col=c("blue"),
                           # baseline=0, v=0, col.line="blue", cex=2, cex.title=1.5,
                           #baseline=0, 
                       v=0, col.line="blue", cex=2, cex.title=1,
                           col.title="black", col.axis="black")

################### 6. rvTDT track

filepath.rvTDT.results<-'/Users/lindagai 1/Documents/classes/4th year/4th term/rvTDT.results.txt'

results<-read.table(filepath.rvTDT.results,header=TRUE)
results<-cbind(df.RV_TDT$Position,results)
head(results)
colnames(results)[1]<-"Position"

rvTDTtrack<-DataTrack(data = -log10(results$p_k_maf),
                       start=results$Position-1,
                       end=results$Position,
                       genome="hg19",
                       chromosome=chr, name="-log10p (rvTDT)",
                       type="l",
                       ylim=c(0, 8),
                       baseline=-log10(bonferroni.sig.level),
                       col.baseline = "red",
                       col=c("blue"),
                       # baseline=0, v=0, col.line="blue", cex=2, cex.title=1.5,
                       #baseline=0,
                       v=0, col.line="blue", cex=2, cex.title=1,
                       col.title="black", col.axis="black")


# #Long version -- fix this to allow all graphs
# 
# results.long<-gather(results, key="test", value="pval",
#                      p_lc_1, p_lc_maf,p_lc_pc,p_k_1,  p_k_maf,p_k_pc)
# head(results.long)
# head(results)
# length(results.long$test)
# length(results.long$test)
# length(results.long$pval)
# 
# rvTDTtrack<-DataTrack(data = as.data.frame(-log10(results.long$pval)),
#                       #data = -log10(results.long$pval),
#                        start=results.long$Position-1,
#                        end=results.long$Position,
#                        genome="hg19",
#                        chromosome=chr, name="-log10p (RV-TDT), \n window size=25M",
#                        type="l",
#                        ylim=c(0, 3),
#                        baseline=-log10(bonferroni.sig.level),
#                        col.baseline = "red",
#                        col=c("blue"),
#                        # baseline=0, v=0, col.line="blue", cex=2, cex.title=1.5,
#                        #baseline=0,
#                        v=0, col.line="blue", cex=2, cex.title=1,
#                        col.title="black", col.axis="black")
# # 
# # plotTracks(dTrack, groups=rep(c("control", "treated"), each=3), type=c("a", "p", "confint"))
# plotTracks(rvTDTtrack, groups=results.long$test)
# 
# nrow(as.data.frame(results.long$pval))

#Notes
# data(twoGroups)
# head(twoGroups)
# dTrack <- DataTrack(twoGroups, name="uniform")
# plotTracks(dTrack)
# plotTracks(dTrack, groups=rep(c("control", "treated"), each=3), type=c("a", "p"), legend=TRUE)


################### 7. ScanTrio p-val track

#filepath.scantrio<-"/users/lgai/GViz_test/scantrio.pos.and.lr.txt"
#filepath.scantrio<-"/users/lgai/8q24_project/data/processed_data/Scan-Trio/scantrio.pos.and.lr.txt"

#Filepath for windows = 100, common and rare
#filepath.scantrio<-'/Users/lindagai 1/Documents/classes/4th year/4th term/Research/scantrio.pos.p.and.lr.txt'

filepath.scantrio<-"/Users/lindagai 1/Documents/classes/4th year/4th term/Research/8q24-euro-stPerm-10k-from-BEAGLEv3.filtered.annotation.rarevar.monomorphs.removed.windowsize=25.overlap24.scantrio.pos.and.lr.txt"

df<-read.table(filepath.scantrio,header=TRUE)
head(df)

df2<-left_join(df,genotypic_tdt_results,by=c("avg.pos"="position"))
head(df2)

ScanTrioLRtrack<-DataTrack(data = -log10(df2$p),
                           start=as.numeric(df2$avg.pos-1),
                           end=as.numeric(df2$avg.pos-1),
                           genome="hg19",
                           chromosome=chr, name="-log10p (ScanTrio)",
                           type="l",
                           ylim=c(0, 8),
                           col=c("blue"),
                           baseline=-log10(bonferroni.sig.level),
                           col.baseline = "red",
                           # baseline=0, v=0, col.line="blue", cex=2, cex.title=1.5,
                           #baseline=0,
                           v=0, col.line="blue", cex=2, cex.title=1,
                           col.title="black", col.axis="black")

#id="name", name = "Visel", rot.title=0, cex.title=1, 

################### Create graph

#filepath.gviz<-"/Users/lindagai 1/Documents/classes/4th year/4th term/Research/GViz_rarev2.pdf"
filepath.gviz<-"/Users/lindagai 1/Documents/classes/4th year/4th term/Research/GViz_rarev2.pdf"

pdf(filepath.gviz,height=16,width=16)

plotTracks(list(itrack,
                gtrack,
                #track.8q24,
                dtrack.gTDT,
                RVtrack,
                RV_TDTtrack,
                rvTDTtrack,
                ScanTrioLRtrack#,
                #txTr
),sizes=c(1,1,5,3,3,3,3),
background.title="darkgray",from= from, to=to)

dev.off()

