#Load in raw scantrio results
#filepath.scantrio<-"/users/lgai/scantrio_tests/testrun_myScanTrio.Rda"
#filepath.scantrio<-"/users/lgai/8q24_project/data/processed_data/Scan-Trio/ScanTrio.results.annotation.filtered.window=25snps.overlap=0.Rda"
#filepath.scantrio<-"/users/lgai/8q24_project/data/processed_data/Scan-Trio/ScanTrio.results.annotation.filtered.window=5snps.overlap=4.Rda"
filepath.scantrio<-"/users/lgai/8q24_project/data/processed_data/Scan-Trio/8q24-euro-stPerm-10k-from-BEAGLEv3.filtered.annotation.rarevar.monomorphs.removed.RESULTS.rda"

filepath.scantrio
load(filepath.scantrio)

head(stOutEither)
library(dplyr)
library(ggplot2)

#get the LR
newdf<-as.data.frame(cbind(start(ranges(stOutEither)), end(ranges(stOutEither)), stPEither, stOutEither$lr))
colnames(newdf)<-c("start","end","p","lr")
head(newdf)
newdf2<- newdf %>% mutate(avg.pos = round(start + (end-start)/2)) %>% mutate(neglog10p = -log10(p))
head(newdf2)

head(newdf2)
range(newdf2$lr)
#0.3759068 0.4384633

filepath.scantrio<-"/users/lgai/8q24_project/data/processed_data/Scan-Trio/8q24-euro-stPerm-10k-from-BEAGLEv3.filtered.annotation.rarevar.monomorphs.removed.windowsize=25.overlap24.scantrio.pos.and.lr.txt"
write.table(newdf2,filepath.scantrio,sep="\t",row.names = FALSE,quote = FALSE)

scp lgai@jhpce01.jhsph.edu:"/users/lgai/8q24_project/data/processed_data/Scan-Trio/8q24-euro-stPerm-10k-from-BEAGLEv3.filtered.annotation.rarevar.monomorphs.removed.windowsize=25.overlap24.scantrio.pos.and.lr.txt"  '/Users/lindagai 1/Documents/classes/4th year/4th term/Research/'
