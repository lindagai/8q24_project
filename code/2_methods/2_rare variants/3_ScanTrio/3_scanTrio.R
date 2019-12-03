#Sign into cluster
ssh -X lgai@jhpce01.jhsph.edu 

#Allocate memory
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G
module load R
R

################################################

library(GenomicRanges)
library(VariantAnnotation)

################################ 4. Run Scan Trio ################################################

nPerms<-100
nCores<-25
permIdx<-split(seq(1,nPerms), rep(1:nCores, each=nPerms/nCores))

currChr<-8
regCode<-"8q24"
pop<-"euro"
windowType<-"M"
#windowSz<-25
windowSz<-5
nPerms<-1000
nCores<-25
permIdx<-split(seq(1,nPerms), rep(1:nCores, each=nPerms/nCores))

#load("/users/mtaub/Research/Trio/targetedSeq/data/8q24/8q24-euro-stPerm-1K.rda")
filepath.perm<-"/users/lgai/8q24_project/data/processed_data/Scan-Trio/8q24-euro-stPerm-10k-from-BEAGLEv3.filtered.annotation.rarevar.monomorphs.removed.rda"
load(filepath.perm)

#Get IDs of rare variants
rareIdx<-which(mafs <= 0.01)
length(rareIdx)

# window specified in kb
if (windowType == "K"){
  window<-GRanges(seqnames=currChr, ranges=IRanges(start=start(tGR)[rareIdx], width=windowSz))
} 

windowType

# window specified in number of markers
if (windowType == "M"){
  #edited this line to get the chr8 to work
  window<-GRanges(seqnames=currChr, ranges=IRanges(start=start(tGR)[rareIdx[1:(length(rareIdx)-windowSz+1)]], end=start(tGR[rareIdx[windowSz:length(rareIdx)]])))
}
window

stOutEither<-myScanTrioEither(window, tGR[rareIdx])

#stPermOutTmpEither<-mclapply(permIdx, function(x) myScanTrioPermEither(window, tGR[rareIdx], permCtMin[rareIdx, x], permCtMaj[rareIdx, x]), mc.cores=nCores)
stPermOutTmpEither<-lapply(permIdx, function(x) myScanTrioPermEither(window, tGR[rareIdx], permCtMin[rareIdx, x], permCtMaj[rareIdx, x]))

stPermOutEither<-do.call(cbind, stPermOutTmpEither)
stPEither<-rowMeans(stOutEither$lr < stPermOutEither)
range(stPEither)

range(stOutEither$lr)

new.filepath<-paste0("/users/lgai/8q24_project/data/processed_data/Scan-Trio/8q24-euro-stPerm-10k-from-BEAGLEv3.filtered.annotation.rarevar.monomorphs.removed.RESULTS.","windowSz",".rda")

save(stOutEither, stPermOutEither, stPEither, mafs, minorAlleles, window, file=new.filepath)
