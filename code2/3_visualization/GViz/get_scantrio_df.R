#Load in raw scantrio results
filepath.scantrio<-"/users/lgai/scantrio_tests/testrun_myScanTrio.Rda"
filepath.scantrio
load(filepath.scantrio)

head(stOutEither)
library(dplyr)
library(ggplot2)

#get the LR
newdf<-as.data.frame(cbind(start(ranges(stOutEither)), end(ranges(stOutEither)), stOutEither$lr))
colnames(newdf)<-c("start","end","lr")
head(newdf)
newdf2<- newdf %>% mutate(avg.pos = round(start + (end-start)/2)) %>% mutate(loglr = log(lr))
head(newdf2)

filepath.scantrio<-"/users/lgai/GViz_test/scantrio.pos.and.lr.txt"
write.table(newdf2,filepath.scantrio,sep="\t",row.names = FALSE,quote = FALSE)
