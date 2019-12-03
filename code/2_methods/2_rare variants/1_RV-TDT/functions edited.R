######################################################

get.rv_tdt.files.in.windows <- function(window.size,overlap,filepath.tped,filepath.map,results.dir){
  #Creates map and tped files using the first n SNPs
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
  map<-read.table(filepath.map,header=FALSE)
  tped<-read.table(filepath.tped,header=FALSE)

  #####################
  n.snps<-nrow(map)
  
  if (overlap==0){
    n.windows<-ceiling(nrow(map)/window.size)  
  } else{
    #TODO: This only works when you shift over by 1 (i.e. overlap = window.size-1) -- fix this to be more general
    n.windows<- n.snps - window.size + 1
  }
  
  #####################
  
  #i<-1
  
  #iii. Split the map and tped files into windows of specified size
  for (i in 1:(n.windows)){
    if (overlap==0){
      start.index<-(i-1)*window.size+1
      end.index<-min(i*window.size, nrow(map))    
    } else{
      start.index<-(i-1)*window.size-(i-1)*overlap+1
      end.index<-min(i*window.size-(i-1)*overlap,n.snps)
    }
    
    file.param<-paste0("window",i,".",
                       start.index,"-",end.index,window.type)
    
    filepath.map.new<-paste0(data.dir,file.param,".map")
    filepath.tped.new<-paste0(data.dir,file.param,".tped")
    filepath.tped.new
    
    sm.map<-map[start.index:end.index,]
    sm.tped<-tped[start.index:end.index,]
    
    write.table(sm.map,filepath.map.new, sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)
    write.table(sm.tped,filepath.tped.new, sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)
    
    #iii. Get a table of snps and positions of each window
    #TODO: make this more efficient/clean it up somehow
    filepath.vcf.snp.pos<-"/users/lgai/8q24_project/data/processed_data/vcf_snp_pos.txt"
    vcf.snp.pos<-read.table(filepath.vcf.snp.pos,sep="\t",quote ="",header=TRUE)

    window.map.with.pos<-left_join(sm.map,vcf.snp.pos,by=c("V2"="snp"))
    
    filepath.window.map.with.pos<-paste0(data.dir,file.param,".map.with.pos.txt")
    
    write.table(window.map.with.pos,filepath.window.map.with.pos, sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)
    
  }
  
}

######################################################

get.rvtdt.results <- function(window.size,overlap,filepath.phen,filepath.tped,filepath.map,results.dir){
  
  rv.tdt.dir<-"./rv-tdt/rvTDT "
  work.dir <- "cd /users/lgai/"
  
  #Creates map and tped files using the first n SNPs
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
  
  #TODO: you need the # of SNPs to calculate the # of SNPs in the window--make this more efficient
  map<-read.table(filepath.map,header=FALSE)
  
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
      end.index<-min(i*window.size, nrow(map))    
    } else{
      start.index<-(i-1)*window.size-(i-1)*overlap+1
      end.index<-min(i*window.size-(i-1)*overlap,nrow(map))
    }
    
    file.param<-paste0("window",i,".",
                       start.index,"-",end.index,window.type)
    
    filepath.map<-paste0(data.dir,file.param,".map")
    filepath.tped<-paste0(data.dir,file.param,".tped")
    
    rv.tdt.results.dir<-paste0(data.dir,file.param)
    rv.tdt.results.dir
    
    command<-paste0(rv.tdt.dir, rv.tdt.results.dir, 
                    " -G ", filepath.tped,
                    " -P ", filepath.phen,
                    " -M ", filepath.map,
                    " --adapt 100 --alpha 0.00001 --permut 1000",
                    " -u 0.01"
    )
    
    command
    
    system(work.dir)
    system(command) 	
    
  }
  
}

######################################################





















################################# All below functions are unedited ##############################

filter.snps.with.annotation <- function(filepath.annotation,cadd=NA,gwava=NA,eigen=NA){
        
        annovar.report<-read.table(filepath.annotation,sep="\t",quote ="",header=TRUE)
        
        subset.command<-paste0("subset.annovar.report<-annovar.report %>%",NA,
                               "filter(WG_GWAVA_score>gwava)")
        
        subset.annovar.report <- annovar.report   
        
        if (!is.na(gwava)){
                subset.annovar.report<-subset.annovar.report %>% filter(WG_GWAVA_score>gwava)
        }
        if (!is.na(cadd)){
                subset.annovar.report<-subset.annovar.report %>% filter(CADD13>cadd) 
        }
        if (!is.na(eigen)){
                subset.annovar.report<-subset.annovar.report %>% filter(WG_EIGEN_score>eigen)  
        }
        
        #dim(subset.annovar.report)        
        #head(subset.annovar.report)
        #colnames(subset.annovar.report)
        
        return(subset.annovar.report$StartPosition)
}

filter.vcf.to.selected.pos <- function(selected.pos,filepath.vcf){
        #filepath.positions<-filepath_caddscore_positions_90th
        #selected.pos<-positions.cadd.filtered.10
        #length(selected.pos)
        vcf<-readVcf(filepath.vcf)
        vcf.geno<-geno(vcf)$GT
        
        #selected.pos<-read.table(filepath.positions,header=TRUE)
        #head(selected.pos)
        
        vcf.snp.pos<-as.data.frame(cbind(names(rowRanges(vcf)), start(rowRanges(vcf))))
        colnames(vcf.snp.pos)<-c("SNP","Position")
        dim(vcf.snp.pos)
        head(vcf.snp.pos)
        
        vcf.snps.filtered<- subset(vcf.snp.pos, Position %in% selected.pos)
        dim(vcf.snps.filtered)
        
        vcf.geno.filtered <- subset(vcf.geno, rownames(vcf.geno) %in% vcf.snps.filtered$SNP)
        #dim(vcf.geno.filtered)
        #vcf.geno.filtered[1:5,1:5]
        #return(vcf.geno.filtered)
}

get.tped <- function(vcf.geno,filepath.tped){
        #RV TDT tped columns:
        # 1) SNP/variant id
        # 2-n) genotype on every individual (for the rest of the columns)
        #vcf.geno.bool<-apply(vcf.geno, 1, function(x) paste(gsub("\\|", "\t", x), collapse="\t"))
        vcf.geno.bool<-apply(vcf.geno, 1, function(x) paste(gsub("\\/", "\t", x), collapse="\t"))
        toWrite<-paste(rownames(vcf.geno), vcf.geno.bool, sep=" ")
        writeLines(toWrite, filepath.tped)
}

get.map <- function(filepath.tped,filepath.map){
        # RV TDT map columns:
        #1.Gene ID (8q24)
        #2. Variant id. The variant id must match with the variant ids in tped file
        #3. MAF
        tped<-read.table(filepath.tped,stringsAsFactors = FALSE)  
        num.snps<-nrow(tped)
        gene.id.vec<-rep('8q24',num.snps)
        snps<-tped[,1]
        mafs<-rowMeans(tped[,-1])
        map<-cbind(gene.id.vec,snps,mafs)
        write.table(map,filepath.map,col.names=FALSE,row.names = FALSE,quote = FALSE)
}

trio.ped.to.phen <- function(filepath.ped,vcf.geno,filepath.phen){
        # RV TDT phen columns:
        #1) sample ID
        #2) family ID
        #3) father ID
        #4) mother ID
        #5) sex (1 for male, 0 for female)
        #6) case(1) / control(0).
        #The order of the sample ID must match the order of individuals in tped file.
        ped<-read.table(filepath.ped,header = TRUE)
        
        rvtdt.phen<-ped %>%
                select(pid,famid,fatid,motid,sex,affected)
        
        rvtdt.phen$sex<-rvtdt.phen$sex %>%
                gsub("female","0", .) %>%
                gsub("male","1",.)
        
        rvtdt.phen$affected<-rvtdt.phen$affected%>%
                gsub("Unaffected","0",.) %>%
                gsub("Affected","1",.)
        
        #IDs must be in same order as in tped
        ids.vcf.order<-as.data.frame(colnames(vcf.geno))
        colnames(ids.vcf.order)<-"pid"
        rvtdt.phen <-left_join(ids.vcf.order,rvtdt.phen)
        write.table(rvtdt.phen,filepath.phen, sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)
}

get.test.case <- function(n,filepath.tped,filepath.map,filepath.test){
        #Creates map and tped files using the first n SNPs
        map<-read.table(filepath.map)
        tped<-read.table(filepath.tped)
        filepath.map.new<-paste0(filepath.test,n,"snps.map")
        filepath.tped.new<-paste0(filepath.test,n,"snps.tped")
        write.table(map[1:n,],filepath.map.new, sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)
        write.table(tped[1:n,],filepath.tped.new, sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)
}