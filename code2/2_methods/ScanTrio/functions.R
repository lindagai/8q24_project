#ScanTrio functions, should go in a separate file

myScanTrioEither<-function(window, tGR){
  tOL<-findOverlaps(tGR, window)
  tOL
  
  minor.in<-sapply(split(queryHits(tOL), subjectHits(tOL)), function(x) sum(tGR[x]$tMinor))
  major.in<-sapply(split(queryHits(tOL), subjectHits(tOL)), function(x) sum(tGR[x]$tMajor))
  
  minor.out<-sum(tGR$tMinor) - minor.in
  major.out<-sum(tGR$tMajor) - major.in
  
  n <- (minor.in + minor.out + 1)/(minor.in + minor.out + 
                                     major.in + major.out + 2)
  y.in <- minor.in
  y.out <- minor.out
  n.in <- minor.in + major.in
  n.out <- minor.out + major.out
  p.in <- (minor.in + 1)/(n.in + 2)
  p.out <- (minor.out + 1)/(n.out + 2)
  
  #Log-likelihood
  loglr <- y.in * log(p.in/n) + y.out * 
    log(p.out/n) + (n.in - y.in) * log((1 - p.in)/(1 - 
                                                     n)) + (n.out - y.out) * log((1 - p.out)/(1 - n))
  lr <- exp(loglr)
  
  windowSub<-window[unique(subjectHits(tOL))]
  meta <- values(windowSub)
  meta2 <- DataFrame(meta, lr = lr, minor.in = as.integer(minor.in), 
                     major.in = as.integer(major.in), minor.out = as.integer(minor.out), 
                     major.out = as.integer(major.out), mendel.in = NA, 
                     mendel.out = NA)
  values(windowSub) <- meta2
  return(windowSub)
}

myScanTrioPermEither<-function(window, tGR, permCtMin, permCtMaj){
  tOL<-findOverlaps(tGR, window)
  minor.in<-sapply(split(queryHits(tOL), subjectHits(tOL)), function(x) colSums(permCtMin[x,, drop=FALSE]))
  major.in<-sapply(split(queryHits(tOL), subjectHits(tOL)), function(x) colSums(permCtMaj[x,, drop=FALSE]))
  minor.out<-colSums(permCtMin) - minor.in
  major.out<-colSums(permCtMaj) - major.in
  
  n <- (minor.in + minor.out + 1)/(minor.in + minor.out + 
                                     major.in + major.out + 2)
  y.in <- minor.in
  y.out <- minor.out
  n.in <- minor.in + major.in
  n.out <- minor.out + major.out
  p.in <- (minor.in + 1)/(n.in + 2)
  p.out <- (minor.out + 1)/(n.out + 2)
  loglr <- y.in * log(p.in/n) + y.out * 
    log(p.out/n) + (n.in - y.in) * log((1 - p.in)/(1 - 
                                                     n)) + (n.out - y.out) * log((1 - p.out)/(1 - n))
  lr <- exp(loglr)
  return(t(lr))
}

#########################