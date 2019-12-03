ssh -Y lgai@jhpce01.jhsph.edu

#Allocate memory
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G
module load R
R

################################################################################

library(VariantAnnotation)
library(tidyverse)

############### 1. Extract RV-TDT information from all windows ###############################

# i<-1
# n.window.sizes<-1
# window.size=25
# overlap=0
# window.type<-"snps"

window.size<-25
overlap<-24
window.type<-"snps"

filepath.tped<-"/users/lgai/8q24_project/data/processed_data/RV-TDT/8q24.cleaned.phased.rarevar.monomorphs.removed.1000G.euro.annotation.tped"
filepath.map<-"/users/lgai/8q24_project/data/processed_data/RV-TDT/8q24.cleaned.phased.rarevar.monomorphs.removed.1000G.euro.annotation.map"

#results.dir<-"/users/lgai/8q24_project/data/processed_data/RV-TDT/"
results.dir<-"./8q24_project/data/processed_data/RV-TDT/"
filepath.phen <- "./8q24_project/data/processed_data/RV-TDT/gmkf_euro_completetrios_ids_match_vcf_mend_errors_removed.phen"

get.rvtdt.results.in.df<-function(window.size,overlap,filepath.phen,filepath.tped,filepath.map,results.dir){
  map<-read.table(filepath.map)
  
  #####################
  n.snps<-nrow(map)
  
  if (overlap==0){
    n.windows<-ceiling(nrow(map)/window.size)  
  } else{
    #TODO: This only works when you shift over by 1 (i.e. overlap = window.size-1) -- fix this to be more general
    n.windows<- n.snps - window.size + 1
  }
  
  #####################
  
  df<-as.data.frame(matrix(nrow=n.windows,ncol=8))
  # head(df)
  
  #i. Read in the results and save them as a df
  for (i in 1:(n.windows)){
    if (overlap==0){
      start.index<-(i-1)*window.size+1
      end.index<-min(i*window.size, nrow(map))    
    } else{
      start.index<-(i-1)*window.size-(i-1)*overlap+1
      end.index<-min(i*window.size-(i-1)*overlap,nrow(map))
    }
    
    window.param<-paste0(
      "window-size=",window.size,window.type,
      ".overlap=",overlap)
    
    # window.param
    
    data.dir <-paste0(results.dir, window.param,"/")
    data.dir
    
    file.param<-paste0("window",i,".",
                       start.index,"-",end.index,window.type)
    
    rv.tdt.results.dir<-paste0(data.dir,file.param,"_pval/8q24.pval")
    # rv.tdt.results.dir
    
    pval.df<-read.table(rv.tdt.results.dir,comment.char = "",header=TRUE,stringsAsFactors=FALSE)
    
    if (i==1){
      colnames(df)<-c(names(pval.df),"Position")
    }
    
    #print(head(pval.df))
    df[i,]<-pval.df[1,]
    # head(df)
    
    #ii. Get the middle position of each window	
    filepath.window.map.with.pos<-paste0(data.dir,file.param,".map.with.pos.txt")
    # filepath.window.map.with.pos
    
    map.with.pos<-read.table(filepath.window.map.with.pos,comment.char = "",header=FALSE)
    # head(map.with.pos)
    
    window.pos<-sort(as.numeric(map.with.pos$V4))
    # window.pos
    
    #Get the largest and smallest, and then set the window's position as halfway between them
    start.pos<-window.pos[1]
    end.pos<-tail(window.pos, 1)
    mid.pos<- start.pos + (end.pos - start.pos)/2
    # typeof(mid.pos)
    # typeof(df$CMC.Analytical)
    df[i,8]<-mid.pos
    # head(df)
  }
  data.dir
  filepath.df.results.with.pos<-paste0(data.dir,"results.from.all.windows.txt")
  filepath.df.results.with.pos
  df
  write.table(df,filepath.df.results.with.pos,sep="\t", row.names=FALSE,quote=FALSE)
}

get.rvtdt.results.in.df(window.size,overlap,filepath.phen,filepath.tped,filepath.map,results.dir)
  

###################### 2. Plot all RV-TDT results on one graph ###################################

scp lgai@jhpce01.jhsph.edu:"/users/lgai/8q24_project/data/processed_data/RV-TDT/window-size=25snps.overlap=24/results.from.all.windows.txt" '/Users/lindagai 1/Documents/classes/4th year/4th term'

filepath.results <-"/Users/lindagai 1/Documents/classes/4th year/4th term/results.from.all.windows.txt"

df.results.with.pos<-read.table(filepath.results,comment.char = "",header=TRUE)

df.results.with.pos.long<-gather(df.results.with.pos, key="test", value="pval",
                                 CMC.Analytical,BRV.Haplo,CMC.Haplo,VT.BRV.Haplo,VT.CMC.Haplo,WSS.Haplo )
head(df.results.with.pos.long)

#i. Plot all RV-TDT results
#Bonferroni corrected
n.windows<-nrow(df.results.with.pos)
n.windows

bonferroni.sig.level<-0.05/n.windows

ggplot() +
  geom_line(data = df.results.with.pos.long, aes(group=test, color = test,
                                                 x = Position, y = pval))+
  geom_hline(yintercept=bonferroni.sig.level, linetype=2, color = "red", size=2) +
  labs(title='RV-TDT (rare var only, annot. filtered, \n window size =25 SNPs, 24 overlap)',
       x ='Position (hg19)', y = 'p-value at center of window')+
  guides(color=guide_legend("RV-TDT test type")) + 
  scale_linetype_manual(name = "Bonferroni-corrected significance", values = 2, 
                        guide = guide_legend(override.aes = list(color = c("red"))))

#ii. Overlay on gTDT results
filepath_genotypic_tdt_results<-"/Users/lindagai 1/Documents/classes/4th year/Research/8q24_project/Update 12:13/8q24_genotypic_tdt_results.txt"
genotypic_tdt_results<-read.table(filepath_genotypic_tdt_results,header=TRUE)
head(genotypic_tdt_results)

#Plot the annot.pos TRUE 2nd
ggplot(genotypic_tdt_results)+
  geom_point(data = genotypic_tdt_results, 
             aes(x=position,y=neglogp),alpha = 0.1) +
  #geom_vline(data = df2, aes(xintercept = as.numeric(pos)))
  labs(title="RV-TDT",
       x="Position (hg19)", y = "-log10(p)") +
  geom_line(data = df.results.with.pos.long, aes(group=test, color = test,
                                                 x = Position, y = -log10(pval)))+
  geom_hline(yintercept=-log10(bonferroni.sig.level), linetype=2, color = "red", size=1) +
  guides(color=guide_legend("Test type"))+
  theme(legend.position = c(.2,.72), plot.title = element_text(hjust = 0.5),
        legend.background = element_rect(color = "black"))

###################### 3. Plot all rvTDT results ###################################
filepath.rvTDT.results<-'/Users/lindagai 1/Documents/classes/4th year/4th term/rvTDT.results.txt'

results<-read.table(filepath.rvTDT.results,header=TRUE)
head(results)
results<-cbind(df.results.with.pos$Position,results)
colnames(results)[1]<-"Position"
head(results)
names(results)

results.long<-gather(results, key="test", value="pval",
                     p_lc_1, p_lc_maf,p_lc_pc,p_k_1,  p_k_maf,p_k_pc)
head(results.long)

n.windows<-nrow(results)
n.windows
bonferroni.sig.level<-0.05/n.windows

################################################################################

#i. Plot all rvTDT results

ggplot() +
  geom_line(data = results.long, aes(group=test, color = test,alpha=0.5,
                                                 x = Position, y = pval))+
  geom_hline(yintercept=bonferroni.sig.level, linetype=2, color = "red", size=2) +
  labs(title='rvTDT (rare var only, annot. filtered, \n window size = 25 SNPs,  24 overlap)',
       x ='Position (hg19)', y = 'p-value at center of window')+
  guides(color=guide_legend("RV-TDT test type")) + 
  scale_linetype_manual(name = "Bonferroni-corrected significance", values = 2, 
                        guide = guide_legend(override.aes = list(color = c("red"))))

################################################################################

#ii. Compare RV-TDT and rvTDT

df.results.with.pos.long
temp1<-df.results.with.pos.long[,2:4]
temp1

temp2<-results.long[,c("Position","test","pval")]
temp2

df<-rbind(temp1,temp2)
head(df)

#TODO: change the colors of the RV and rv TDTs so it's easier to tell what's what
ggplot() +
  geom_line(data = df, aes(group=test, color = test,
                                     x = Position, y = pval))+
  geom_hline(yintercept=bonferroni.sig.level, linetype=2, color = "red", size=2) +
  labs(title='RV and rv TDT (rare var only, annot. filtered, \n window size = 25 SNPs,  24 overlap)',
       x ='Position (hg19)', y = 'p-value at center of window')+
  guides(color=guide_legend("RV-TDT test type")) + 
  scale_linetype_manual(name = "Bonferroni-corrected significance", values = 2, 
                        guide = guide_legend(override.aes = list(color = c("red"))))

################################################################################

#iii. Compare rvTDT to gTDT
filepath_genotypic_tdt_results<-"/Users/lindagai 1/Documents/classes/4th year/Research/8q24_project/Update 12:13/8q24_genotypic_tdt_results.txt"
genotypic_tdt_results<-read.table(filepath_genotypic_tdt_results,header=TRUE)
head(genotypic_tdt_results)

#Plot the annot.pos TRUE 2nd
ggplot(genotypic_tdt_results)+
  geom_point(data = genotypic_tdt_results, 
             aes(x=position,y=neglogp),alpha = 0.1) +
  #geom_vline(data = df2, aes(xintercept = as.numeric(pos)))
  labs(title="rvTDT",
       x="Position (hg19)", y = "-log10(p)") +
  geom_line(data = results.long, aes(group=test, color = test,
                                                 x = Position, y = -log10(pval)))+
  geom_hline(yintercept=-log10(bonferroni.sig.level), linetype=2, color = "red", size=1) +
  guides(color=guide_legend("Test type"))+
  scale_color_hue(labels = c("Kernel", "Kernel-MAF","Kernel-PC","Linear \ncombination","Linear \ncombination-MAF","Linear \ncombination-PC"))+
  theme(legend.position = c(.2,.72), plot.title = element_text(hjust = 0.5),
        legend.background = element_rect(color = "black"))


#Actual clean graph
ggplot(genotypic_tdt_results)+
  geom_point(data = genotypic_tdt_results, 
             aes(x=position,y=neglogp),alpha = 0.1) +
  #geom_vline(data = df2, aes(xintercept = as.numeric(pos)))
  labs(title="RV-TDT",
       x="Position (hg19)", y = "-log10(p)") +
  geom_line(data = df.results.with.pos.long, aes(group=test, color = test,
                                                 x = Position, y = -log10(pval)))+
  geom_hline(yintercept=-log10(bonferroni.sig.level), linetype=2, color = "red", size=1) +
  guides(color=guide_legend("Test type"))+
  theme(legend.position = c(.2,.72), plot.title = element_text(hjust = 0.5),
        legend.background = element_rect(color = "black"))

###################### 3. Plot Scan-Trio results ###################################

#i. gTDT and Scan-Trio
#ii. RV-TDT and Scan-Trio
#iii. rvTDT and Scan-Trio



