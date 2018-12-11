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
        vcf.geno.bool<-apply(vcf.geno, 1, function(x) paste(gsub("\\|", "\t", x), collapse="\t"))
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


get.rvtdt.files.in.windows <- function(window.size=100,overlap=0,filepath.tped,filepath.map,filepath.test){
        #Creates map and tped files using the first n SNPs
        map<-read.table(filepath.map)
        tped<-read.table(filepath.tped)
        
        #Split the map and tped files into windows of size 100
        n.windows<-nrow(map) %/% window.size
        #nrow(map)
        #n.windows
        
        for (i in 1:(n.windows+1)){
                if (overlap==0){
                        start.index<-(i-1)*window.size+1
                        end.index<-min(i*window.size, nrow(map))    
                } else{
                        start.index<-(i-1)*window.size-(i-1)*overlap+1
                        end.index<-min(i*window.size-(i-1)*overlap,nrow(map))
                }
                #print(start.index)
                #print(end.index)
                filepath.map.new<-paste0(filepath.test,"window-size=",window.size,"_snps","_",i,".map")
                filepath.tped.new<-paste0(filepath.test,"window-size=",window.size,"_snps","_",i,".tped")
                #print(filepath.map.new)
                #print(filepath.tped.new)
                
                sm.map<-map[start.index:end.index,]
                sm.tped<-tped[start.index:end.index,]
                
                write.table(sm.map,filepath.map.new, sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)
                write.table(sm.tped,filepath.tped.new, sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)
        }
        
}