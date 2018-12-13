############ 3. Create permuted data to obtain p-values for the ScanTrio statistic. ##############
library(GenomicRanges)
library(VariantAnnotation)
#library(parallel)

nPerms<-100
nCores<-25
permIdx<-split(seq(1,nPerms), rep(1:nCores, each=nPerms/nCores))

currChr<-8
regCode<-"8q24"
pop<-"euro"
windowType<-"M"
windowSz<-100
nPerms<-1000
nCores<-25
permIdx<-split(seq(1,nPerms), rep(1:nCores, each=nPerms/nCores))

# create permuted data sets and calculate counts
# these can be reused for different window configurations

filepath.beagle.v3<-"/users/lgai/scratch/data/8q24.recode.FILTERED.plink.nohalfcalls.DuplicatesRemoved.BEAGLE.BEAGLEv3.bgl.gz"
colCt<-length(scan(gzfile(filepath.beagle.v3), what="character", nlines=1))
phaseMat<-matrix(scan(gzfile(filepath.beagle.v3), what="character"),ncol=colCt, byrow=TRUE)

markerDF<-as.data.frame(read.table(filepath.markerDF,header=TRUE))

tmpRow<-phaseMat[,2][-1]
head(tmpRow)
posInfo<-as.numeric(markerDF[match(tmpRow, markerDF[,1]),2])
head(posInfo)
sum(!is.na(posInfo))

tmpCol<-phaseMat[1, -c(1:2)]
phaseMat<-phaseMat[-1,-c(1,2)]
rownames(phaseMat)<-posInfo
colnames(phaseMat)<-tmpCol

mafs<-apply(phaseMat, 1, function(x){min(prop.table(table(x)))	})
head(mafs)

minorAlleles<-apply(phaseMat, 1, function(x){names(which.min(prop.table(table(x))))})

minBase<-as.numeric(posInfo[1])
maxBase<-as.numeric(posInfo[length(posInfo)])

#TODO: Why are there NAs in posInfo?
#get rid of NAs, figure out why they're there later
varGR<-GRanges(seqnames=currChr, ranges=IRanges(start=as.numeric(na.omit(posInfo)), width=1))
#varGR<-GRanges(seqnames=currChr, ranges=IRanges(start=as.numeric(posInfo), width=1))

tIdx<-seq(1, ncol(phaseMat), by=2)
uIdx<-seq(2, ncol(phaseMat), by=2)

tGR<-myTransCount(phaseMat, varGR, minorAlleles, tIdx, uIdx, TRUE)

## perform permutations to allow significance assessment for scanTrio
nFams<-ncol(phaseMat)/2
idxMat<-matrix(1:ncol(phaseMat), ncol=2, byrow=TRUE)
randIdxMat<-matrix(sample(c(1,2), nFams*nPerms, replace=TRUE), nrow=nFams, ncol=nPerms)
tPermIdx<-matrix(0, nrow=nFams, ncol=nPerms)
uPermIdx<-matrix(0, nrow=nFams, ncol=nPerms)
i<-1

for (i in 1:nPerms){
  tPermIdx[,i]<-idxMat[cbind(1:nrow(idxMat), randIdxMat[,i])]
  uPermIdx[,i]<-idxMat[cbind(1:nrow(idxMat), 3-randIdxMat[,i])]	
}

doPermCt<-function(idx, phaseMat, varGR, minorAlleles, tPermIdx, uPermIdx){
  permCtMin<-matrix(0, nrow=nrow(phaseMat), ncol=length(idx))
  permCtMaj<-matrix(0, nrow=nrow(phaseMat), ncol=length(idx))
  for (i in 1:length(idx)){
    if (i%%10 == 0) print(i)
    tmp<-myTransCount(phaseMat, varGR, minorAlleles, tPermIdx[,idx[i]], uPermIdx[,idx[i]], FALSE)
    permCtMin[,i]<-tmp$tMinor
    permCtMaj[,i]<-tmp$tMajor			
  }
  return(list(permCtMin=permCtMin, permCtMaj=permCtMaj))
}

#permOut<-mclapply(permIdx, doPermCt, phaseMat, varGR, minorAlleles, tPermIdx, uPermIdx, mc.cores=nCores)
permOut<-lapply(permIdx, doPermCt, phaseMat, varGR, minorAlleles, tPermIdx, uPermIdx)

permCtMin<-do.call(cbind, lapply(permOut, "[[", "permCtMin"))
permCtMaj<-do.call(cbind, lapply(permOut, "[[", "permCtMaj"))
permCtMaj[1:5]

filepath.perm<-"/users/lgai/scantrio_tests/8q24-euro-stPerm-1K-from-BEAGLEv3.rda"
save(tGR, tPermIdx, uPermIdx, permCtMin, permCtMaj, mafs, minorAlleles, file=filepath.perm)


##########################################################################

myTransCount<-function(phaseMat, varGR, minorAlleles, tIdx, uIdx, returnFull=TRUE){
  #transmitted allele counts?
  isHet<-phaseMat[,tIdx] != phaseMat[,uIdx]
  tmp<-phaseMat[,tIdx] == minorAlleles #get rows in phaseMat[,tIdx] with minor alleles
  tmp[!isHet]<-NA #set
  
  tMinor<-rowSums(tmp, na.rm=TRUE) # no. of minor alleles in phaseMat's tIdx col.
  # tmp<-phaseMat[,tIdx] != minorAlleles
  # tmp[!isHet]<-NA
  tMajor<-rowSums(!tmp, na.rm=TRUE) #remove non-het individuals, and 
  
  if (returnFull){
    elementMetadata(varGR)<-cbind(tMinor=tMinor, tMajor=tMajor)
    return(varGR)
  } else {
    return(list(tMinor=tMinor, tMajor=tMajor))
  }
}
