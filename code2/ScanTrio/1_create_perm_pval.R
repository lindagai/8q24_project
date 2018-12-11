
#Un-comment this out if you don't already have BEAGLE4.0
#NOTE: Do not use BEAGLE4.1 or 5.0: these do not use pedigree information.

#Download vcf2beagle
filepath.beagle4<-"/users/lgai/beagle.r1399.jar"
dl.beagle4<-paste0("wget -O ",filepath.beagle4," https://faculty.washington.edu/browning/beagle/beagle.r1399.jar")
system(dl.beagle4)

#Phase vcf
filepath.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/8q24.recode.15.mend.error.individuals.removed.recode.vcf"
filepath.ped<-"/users/lgai/8q24_project/data/processed_data/gmkf_euro_completetrios_ids_match_vcf_mend_errors_removed.txt"
filepath.phased.vcf<-"users/lgai/8q24_project/data/processed_data/vcfs/8q24.recode.15.mend.error.individuals.removed.recode.BEAGLE"

phase.command<-paste0("java -Xmx10000m -jar ", filepath.beagle4,
                      " gt=",filepath.vcf,
                      " ped=",filepath.ped,
                      " out=",filepath.phased.vcf)

system(phase.command)

################ 1. Use vcf2beagle to convert BEAGLEv4 vcf into a BEAGLE v3 file. ##############

#Download vcf2beagle
#Un-comment this out if you don't already have vcf2beagle.jar

# #Download vcf2beagle
# filepath.vcf2beagle<-"/users/lgai/vcf2beagle.jar"
# dl.vcf2beagle<-paste0("wget -O ",filepath.vcf2beagle," https://faculty.washington.edu/browning/beagle_utilities/vcf2beagle.jar")
# system(dl.vcf2beagle)
# 
# #Convert BEAGLEv4.gz to BEAGLEv3
# filepath.BEAGLE.v4.gz <- "/users/lgai/scratch/data/8q24.recode.FILTERED.plink.nohalfcalls.DuplicatesRemoved.BEAGLE.vcf"
# filepath.BEAGLE.v3<-"/users/lgai/scratch/data/8q24.recode.FILTERED.plink.nohalfcalls.DuplicatesRemoved.BEAGLE.BEAGLEv3"
# 
# # TODO: Check to make sure missing calls are dealt with appropriately... is it '.' or './.'?
# convert2beaglev3<-paste("zcat",filepath.BEAGLE.v4.gz, "|","java -jar /users/lgai/vcf2beagle.jar .", filepath.BEAGLE.v3)
# system(convert2beaglev3)

################ 2. Create markerDF, which gives the corresponding position for each SNP name/rsID ##############

#Don't need to comment this back in unless you want to remake the markerDF

# TODO: Unsure that this is supposed to be from the vcf, or from something like annovar?
# #This can possibly be sped up by readGeno or something like that
# #Fix this later
# TODO: why are some of the alt alleles missing? do we use them?

# filepath_vcf<-"/users/lgai/scratch/data/8q24.recode.FILTERED.plink.nohalfcalls.DuplicatesRemoved.BEAGLE.vcf"
# vcf <- readVcf(filepath_vcf,"hg19")
# 
# vcf.rowranges<-rowRanges(vcf)
# head(vcf.rowranges)
# rsids<-names(rowRanges(vcf))
# pos<-start(rowRanges(vcf))
# ref<-as.character(ref(vcf))
# alt<-as.character(unlist(alt(vcf)))
# 
# markerDF<-cbind(rsids,pos,ref,alt)
# head(markerDF)
# 
# filepath.markerDF<-"/users/lgai/scantrio_tests/markerDF.txt"

#write.table(markerDF, filepath.markerDF, sep=" ",row.names = FALSE,quote = TRUE)

############ 3. Create permuted data to obtain p-values for the ScanTrio statistic. ##############o

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
