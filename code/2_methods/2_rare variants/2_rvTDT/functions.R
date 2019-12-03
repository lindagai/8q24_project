#get evs files

get.evs.raw<-function(filepath.annovar,filepath.evs.raw){
  
  annovar<-read.table(filepath.annovar,sep="\t",header=TRUE, quote ="")
  
  # evs.raw<-annovar %>%
  #   filter(!is.na(Afalt_1000g2015aug_eur)) %>%
  #   filter(TotalDepth>20) %>%
  #   mutate(geno0 = (1-Afalt_1000g2015aug_eur)^2 + 0.0000001) %>%
  #   mutate(geno1 = (1-Afalt_1000g2015aug_eur)*Afalt_1000g2015aug_eur + 0.0000001)  %>%
  #   mutate(geno2 = (Afalt_1000g2015aug_eur)^2*10000 + 0.0000001)  %>%
  #   dplyr::select("StartPosition",
  #                 "geno2",
  #                 "geno1",
  #                 "geno0",
  #                 "TotalDepth"
  #   )
  
  evs.raw<-annovar %>%
    filter(!is.na(Afalt_1000g2015aug_eur)) %>%
    filter(TotalDepth>20) %>%
    mutate(geno0 = (1-Afalt_1000g2015aug_eur)^2) %>%
    mutate(geno1 = (1-Afalt_1000g2015aug_eur)*Afalt_1000g2015aug_eur)  %>%
    mutate(geno2 = (Afalt_1000g2015aug_eur)^2)  %>%
    dplyr::select("StartPosition",
                  "geno2",
                  "geno1",
                  "geno0",
                  "TotalDepth"
    )
  
  write.table(evs.raw, filepath.evs.raw, sep="\t",row.names = FALSE,quote = FALSE)
}

get.evs.filtered<-function(filepath.evs.raw,filepath.vcf,filepath.evs.filtered){
  evs.raw<-read.table(filepath.evs.raw,sep="\t",header=TRUE, quote ="")
  
  #TODO: MAKE THIS MORE EFFICIENT
  vcf<-readVcf(filepath.vcf)
  pos<-start(rowRanges(vcf)) 
  
  evs.filtered<-evs.raw[evs.raw$StartPosition %in% pos,]
  write.table(evs.filtered, filepath.evs.filtered, sep="\t",row.names = FALSE,quote = FALSE)
}

##########################################################################################

#get rvTDT results

get.rvTDT.files.in.windows <- function(window.size,overlap,filepath.ped,filepath.evs,results.dir){

  #Creates map and tped files using the first n SNPs
  evs<-read.table(filepath.evs,header=TRUE)
  dim(evs)
  ped<-read.table(filepath.ped,header=TRUE)
  dim(ped)
  
  #i. Create directory for results
  #TODO: This could probably be a new function
  window.param<-paste0(
    "window-size=",window.size,window.type,
    ".overlap=",overlap)
  
  window.param
  
  data.dir <-paste0(results.dir, window.param,"/")
  data.dir
  
  if(!dir.exists(data.dir)) {
    dir.create(data.dir)
  }
  
  #ii. Read in files and get # of windows
  
  #####################
  
  n.snps<-nrow(map)
  
  if (overlap==0){
    n.windows<-ceiling(nrow(map)/window.size)  
  } else{
    #TODO: This only works when you shift over by 1 (i.e. overlap = window.size-1) -- fix this to be more general
    n.windows<- n.snps - window.size + 1
  }
  
  #####################
  
  #iii. Split the map and tped files into windows of specified size
  for (i in 1:(n.windows)){
    if (overlap==0){
      start.index<-(i-1)*window.size+1
      end.index<-min(i*window.size, nrow(evs))    
    } else{
      start.index<-(i-1)*window.size-(i-1)*overlap+1
      end.index<-min(i*window.size-(i-1)*overlap,nrow(evs))
    }
    
    file.param<-paste0("window",i,".",
                       start.index,"-",end.index,window.type)
    
    filepath.evs.new<-paste0(data.dir,file.param,".evs.txt")
    filepath.ped.new<-paste0(data.dir,file.param,".ped.txt")
    filepath.evs.new    
    
    sm.evs<-evs[start.index:end.index,]
    sm.ped<-ped[,start.index:end.index]
    dim(sm.ped)
    dim(ped)
    
    write.table(sm.evs,filepath.evs.new, sep="\t",col.names=TRUE,row.names = FALSE,quote = FALSE)
    write.table(sm.ped,filepath.ped.new, sep="\t",col.names=TRUE,row.names = FALSE,quote = FALSE)
    
  }
  
}

get.rvTDT.results.in.df <- function(window.size,overlap,filepath.ped,filepath.evs,results.dir){
  #i. Create directory for results
  #TODO: This could probably be a new function
  window.param<-paste0(
    "window-size=",window.size,window.type,
    ".overlap=",overlap)
  
  window.param
  
  data.dir <-paste0(results.dir, window.param,"/")
  data.dir
  
  if(!dir.exists(data.dir)) {
    dir.create(data.dir)
  }
  
  #ii. Read in files and get # of windows
  evs<-read.table(filepath.evs,header=TRUE)
  ped<-read.table(filepath.ped,header=TRUE)
  
  #####################
  
  n.snps<-nrow(map)
  
  if (overlap==0){
    n.windows<-ceiling(nrow(map)/window.size)  
  } else{
    #TODO: This only works when you shift over by 1 (i.e. overlap = window.size-1) -- fix this to be more general
    n.windows<- n.snps - window.size + 1
  }
  
  #####################
  
  all.results<-as.data.frame(matrix(nrow=n.windows,ncol=9))
  all.results
  
  for (i in 1:(n.windows)){
    if (overlap==0){
      start.index<-(i-1)*window.size+1
      end.index<-min(i*window.size, nrow(evs))    
    } else{
      start.index<-(i-1)*window.size-(i-1)*overlap+1
      end.index<-min(i*window.size-(i-1)*overlap,nrow(evs))
    }
    
    file.param<-paste0("window",i,".",
                       start.index,"-",end.index,window.type)
    
    filepath.evs.new<-paste0(data.dir,file.param,".evs.txt")
    filepath.ped.new<-paste0(data.dir,file.param,".ped.txt")
    
    sm.evs<-read.table(filepath.evs.new,header=TRUE)
    sm.ped<-read.table(filepath.ped.new,header=TRUE)
    head(sm.evs)
    head(sm.ped)
    
    rvTDT.results <- rvTDT(sm.ped,sm.evs,maf.threshold=1)
    
    if(i==1) {
      colnames(all.results)<-names(rvTDT.results)
    }
    
    all.results[i,]<-unlist(rvTDT.results)
    all.results
  }
  
  all.results
  filepath.rvTDT.results<-paste0(data.dir,"rvTDT.results.txt")
  filepath.rvTDT.results
  write.table(all.results,filepath.rvTDT.results, sep="\t",col.names=TRUE,row.names = FALSE,quote = FALSE)
}
