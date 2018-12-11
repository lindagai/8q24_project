#testthat.R

######################

###### Checks for get.rvTDT.ped.from.trio.geno -- this should be automated

clean_ped<-read.csv("/Users/lindagai 1/Documents/classes/4th year/2nd term/Research/testing/clean_pedfile_chr8.csv")
library(dplyr)

#Are there nonempty values in both motid and fatid?
kids <-rownames(geno[c.index,])
kids_test <- clean_ped %>%
  filter(pid %in% kids)  %>%
  filter(motid !=0) %>%
  filter(fatid != 0)
dim(kids_test)

#Are there empty values in both motid and fatid, and sex is female?
moms <-rownames(geno[m.index,])
moms_test <- clean_ped %>%
  filter(pid %in% moms)  %>%
  filter(motid ==0) %>%
  filter(fatid ==0) %>%
  filter(sex == "female")
dim(moms_test)

#Are there empty values in both motid and fatid, and sex is male?
dads <-rownames(geno[f.index,])
dads_test <- clean_ped %>%
  filter(pid %in% dads)  %>%
  filter(motid==0) %>%
  filter(fatid== 0) %>%
  filter(sex== "male")
dim(dads_test)
