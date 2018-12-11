

################################ 4. Run Scan Trio ################################################

#load("/users/mtaub/Research/Trio/targetedSeq/data/8q24/8q24-euro-stPerm-1K.rda")
load("/users/lgai/scantrio_tests/8q24-euro-stPerm-1K-from-BEAGLEv3.rda")

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


stOutEither<-myScanTrioEither(window, tGR[rareIdx])

#stPermOutTmpEither<-mclapply(permIdx, function(x) myScanTrioPermEither(window, tGR[rareIdx], permCtMin[rareIdx, x], permCtMaj[rareIdx, x]), mc.cores=nCores)
stPermOutTmpEither<-lapply(permIdx, function(x) myScanTrioPermEither(window, tGR[rareIdx], permCtMin[rareIdx, x], permCtMaj[rareIdx, x]))

stPermOutEither<-do.call(cbind, stPermOutTmpEither)
stPEither<-rowMeans(stOutEither$lr < stPermOutEither)
new.filepath<-"/users/lgai/scantrio_tests/testrun_myScanTrio.Rda"
save(stOutEither, stPermOutEither, stPEither, mafs, minorAlleles, window, file=new.filepath)
