#######NOTE: Before running this code, run functions.R to define functions used here.

#Description:
#Filters 8q24 region to SNPs with CADD>10 and runs RV TDT

#1. Read in phased vcf file

#2. Subset vcf files depending on annotation information (CADD, eigen, gwava)
#i.e., subset to list of positions in AnnoVar that satisfy requirements
#TODO:
#i. Obtain additional scores for 8q24 region: e.g. SIFT, PolyPhen and plot
#       SIFT, PolyPhen scores not available for this regions
#ii. Use existing VariantAnnotation functions e.g. locateVariant and predictCoding to filter

#3. Get vcf file into RV TDT format
#The phen file does not need to be recreated each time.
#i. trio ped format -> phen file
#TODO: change this so it converts raw plink -> phen file directly

#The tped and map files need to be recreated each time you run this
#ii. vcf -> tped file
#iii. vcf -> map file

#4. Create n-SNP test files if needed
#(max. RV TDT can handle at a time is ~1000)

#5. Run RV TDT

#6. TODO: Add visualization

#Author: Linda Gai
#Last update: 7/14/18

################################################################################
ssh -Y lgai@jhpce01.jhsph.edu

#Allocate memory
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G
module load R
R

library(VariantAnnotation)
library(dplyr)

################################################################################

#Before running the commands in this file, run functions.R.

############################# 0. Filter vcf ##########################

#i. Find positions you want
#Get function for selecting the criterion you want to filter by
#Writes out list of positions which can then be used to filter the vcf

#Get ANNOVAR report for 8q24 region with updated CADD scores
filepath.annovar<-"/users/lgai/8q24_project/data/processed_data/annotation/8q24_annovar_reports_CADD13.txt"

#90th percentile values from exploratory analysis of 8q24
#cadd<-7.09
#gwava<-0.45
#eigen<-0.22685

#CADD: suggested minimum cut-off
cadd<-10

#Obtain SNP positions with cadd, eigen, and gwava scores above 90th percentile of 8q24
positions.cadd.filtered.10<-filter.snps.with.annotation(filepath.annotation=filepath.annovar,cadd)
length(positions.cadd.filtered.10)
positions.cadd.filtered.10[1:10]

#GWAVA: 90th percentile
gwava<-0.45

positions.gwava.filtered.10<-filter.snps.with.annotation(filepath.annotation=filepath.annovar,gwava)
length(positions.gwava.filtered.10)
positions.gwava.filtered.10[1:10]

#eigen: 90th percentile
eigen<-0.22685

positions.eigen.filtered.10<-filter.snps.with.annotation(filepath.annotation=filepath.annovar,eigen)
length(positions.eigen.filtered.10)
positions.eigen.filtered.10[1:10]

#ii. Filter the vcf to the selected positions
#Currently using:
#CADD filtered
filepath.vcf<-"/users/lgai/2_5_2018_RVTDT/8q24.recode.FILTERED.FILTERED2.BEAGLE.vcf.gz"

#Read in selected positions
#To recreate the results using 90th percentile cut-off for CADD scores:
filepath_caddscore_positions_90th<-"/users/lgai/8q24_project/data/processed_data/annotation/caddscore_90th_percentile_positions_for_vcf.txt"
caddscore_positions_90th<-read.table(filepath_caddscore_positions_90th,header=TRUE)
head(caddscore_positions_90th)
vcf.geno<-filter.vcf.to.selected.pos(caddscore_positions_90th$Position,filepath.vcf)

#Not using 
vcf.geno<-filter.vcf.to.selected.pos(positions.cadd.filtered.10,filepath.vcf)
#vcf.geno[1:5,1:5]
dim(vcf.geno)
#810 981

#########

#GWAVA 90th percentile filtered
vcf.geno<-filter.vcf.to.selected.pos(positions.gwava.filtered.10,filepath.vcf)
#vcf.geno[1:5,1:5]
dim(vcf.geno)
#11380   981
positions<-rownames(vcf.geno)

#########

#EIGEN 90th percentile filtered

vcf.geno<-filter.vcf.to.selected.pos(positions.eigen.filtered.10,filepath.vcf)
#vcf.geno[1:5,1:5]
dim(vcf.geno)
#12011   981

############################# 0.5. Filter to windows ##########################

#i. window containing all variants around peak in gTDT
#Estimate: 
129750000
129875000-130000000

#Using positions from 90th percentile CADD
#peak.positions<-positions[positions > 129875000 && positions > 130000000]
positions<-caddscore_positions_90th$Position
length(positions)

peak.positions<-positions[positions > 129875000 & positions < 130000000]
12500
length(peak.positions)
151

100, 100 (overlapping windows)
75, 76 (non-overlapping windows) (of a different size, though)

#Filter to CADD score > 10 and within peak
vcf.geno<-filter.vcf.to.selected.pos(peak.positions,filepath.vcf)

#ii. windows of 100 variants around peak in gTDT, non-overlapping
#iii. windows of 100 variants of all variants in 8q24 region, non-overlapping
#iv. windows of 100 variants of 

######################## 1. Get tped ############################################

#CADD filtered version
#filepath.tped<-"~/cadd_10/cadd_10.tped"
#get.tped(vcf.geno,filepath.tped)

#Checks
#tped<-read.table(filepath.tped)
#tped[1:5,1:10]
#dim(tped) #810 1963

filepath.tped<-"~/cadd_10/cadd_10_peak.tped"
get.tped(vcf.geno,filepath.tped)

#Checks
tped<-read.table(filepath.tped)
tped[1:5,1:10]
dim(tped) #151 1963

############################## 2. Get map file ###################################################

# filepath_tped<-"~/cadd_10/cadd_10.tped"
# filepath_rvtdt_map<-"~/cadd_10/cadd_10.map"
# 
# get.map(filepath_tped,filepath_rvtdt_map)
# 
# map<-read.table(filepath_rvtdt_map)
# #head(map)
# #dim(map)
# #1283    3

filepath_tped<-"~/cadd_10/cadd_10_peak.tped"
filepath_rvtdt_map<-"~/cadd_10/cadd_10_peak.map"

get.map(filepath_tped,filepath_rvtdt_map)

map<-read.table(filepath_rvtdt_map)
head(map)
dim(map)
#151   3

############################ 3. Get phen file #############################################

#TRIO columns
#famid                pid              fatid              motid    sex  affected

filepath.trio.ped<-"/users/lgai/8q24_project/data/processed_data/gmkf_euro_completetrios_ids_match_vcf_mend_errors_removed.txt"
filepath.phen<-"~/cadd_10/cadd_10.phen"

chr8_ped<-read.table(filepath.trio.ped,header=TRUE)
head(chr8_ped)
dim(chr8_ped)



trio.ped.to.phen(filepath.trio.ped,vcf.geno,filepath.phen)

########################### 4. Get test case ##############################################

#Unnecessary if using <1000 SNPs
test.filepath<-"~/cadd_10/cadd_10."
filepath_tped<-"~/cadd_10/cadd_10.tped"
filepath_rvtdt_map<-"~/cadd_10/cadd_10.map"
#100-950 SNPs works

get.test.case(500,filepath_tped,filepath_rvtdt_map,test.filepath)

###########################

#Break files into windows

#100 SNP windows, whole 8q24 region, no overlap
window.size=100
filepath.tped<-"~/cadd_10/cadd_10.tped"
filepath.map<-"~/cadd_10/cadd_10.map"
filepath.test<-"~/cadd_10/window100_no_overlap/cadd_10_"
#overlap=0

get.rvtdt.files.in.windows(window.size=100,overlap=0,filepath.tped,filepath.map,filepath.test)

#######

#100 SNP windows, whole 8q24 region, 10 SNP overlap
window.size=100
filepath.tped<-"~/cadd_10/cadd_10.tped"
filepath.map<-"~/cadd_10/cadd_10.map"
filepath.test<-"~/cadd_10/window100_10_overlap/cadd_10_"
#overlap=0
overlap=10

get.rvtdt.files.in.windows(window.size=100,overlap=10,filepath.tped,filepath.map,filepath.test)

#######

#100 SNP windows, peak only, no overlap

mkdir peak_10_overlap
mkdir peak_no_overlap

window.size=100
filepath.tped<-"~/cadd_10/cadd_10_peak.tped"
filepath.map<-"~/cadd_10/cadd_10_peak.map"
filepath.test<-"~/cadd_10/peak_10_overlap/cadd_10_peak_"
#overlap=0
overlap=10

get.rvtdt.files.in.windows(window.size=100,overlap=10,filepath.tped,filepath.map,filepath.test)

#######

#76 SNP windows, peak only, no overlap
window.size=76
filepath.tped<-"~/cadd_10/cadd_10_peak.tped"
filepath.map<-"~/cadd_10/cadd_10_peak.map"
filepath.test<-"~/cadd_10/peak_no_overlap/cadd_10_peak_"
#overlap=0

get.rvtdt.files.in.windows(window.size,overlap=0,filepath.tped,filepath.map,filepath.test)

#######

########################### 5. Run RV TDT ##############################################

#Whole 8q24 region, 10 overlap, 100 SNPs
for (i in 1:9){
        rv.tdt.dir<-"./rv-tdt/rvTDT "
        work.dir <- "cd /users/lgai/"
        data.dir <- "./cadd_10/window100_no_overlap/"
        data.file <- paste0("cadd_10_window-size=100_snps_",i)
        results.dir <- paste0("./cadd_10/window100_no_overlap/cadd_10_window-size=100_snps_",i)
        g.file <- paste0(data.file,".tped")
        p.file <- "./cadd_10/cadd_10.phen"
        #p.file <- paste0(data.file,".phen")
        m.file <- paste0(data.file,".map")
        
        #This works
        command<-paste0(rv.tdt.dir, results.dir, 
                        " -G ", data.dir, g.file,
                        " -P ", p.file,
                        " -M ", data.dir, m.file,
                        " --adapt 100 --alpha 0.00001 --permut 1000",
                        " -u 0.01"
        )
        
        command
        
        system(work.dir)
        system(command)     
}

#Whole 8q24 region, 10 overlap, 100 SNPs
for (i in 1:9){
        rv.tdt.dir<-"./rv-tdt/rvTDT "
        work.dir <- "cd /users/lgai/"
        data.dir <- "./cadd_10/window100_10_overlap/"
        data.file <- paste0("cadd_10_window-size=100_snps_",i)
        results.dir <- paste0("./cadd_10/window100_10_overlap/cadd_10_window-size=100_snps_",i)
        g.file <- paste0(data.file,".tped")
        p.file <- "./cadd_10/cadd_10.phen"
        #p.file <- paste0(data.file,".phen")
        m.file <- paste0(data.file,".map")
        
        #This works
        command<-paste0(rv.tdt.dir, results.dir, 
                        " -G ", data.dir, g.file,
                        " -P ", p.file,
                        " -M ", data.dir, m.file,
                        " --adapt 100 --alpha 0.00001 --permut 1000",
                        " -u 0.01"
        )
        
        command
        
        system(work.dir)
        system(command)     
}

#Peaks, 10 overlap, 100 SNPs
for (i in 1:2){
        rv.tdt.dir<-"./rv-tdt/rvTDT "
        work.dir <- "cd /users/lgai/"
        data.dir <- "./cadd_10/peak_10_overlap/"
        data.file <- paste0("cadd_10_peak_window-size=100_snps_",i)
        results.dir <- paste0("./cadd_10/peak_10_overlap/cadd_10_peak_window-size=100_snps_",i)
        g.file <- paste0(data.file,".tped")
        p.file <- "./cadd_10/cadd_10.phen"
        #p.file <- paste0(data.file,".phen")
        m.file <- paste0(data.file,".map")
        
        #This works
        command<-paste0(rv.tdt.dir, results.dir, 
                        " -G ", data.dir, g.file,
                        " -P ", p.file,
                        " -M ", data.dir, m.file,
                        " --adapt 100 --alpha 0.00001 --permut 1000",
                        " -u 0.01"
        )
        
        command
        
        system(work.dir)
        system(command)     
}


#Peaks, no overlap, 76 SNPs
for (i in 1:2){
        rv.tdt.dir<-"./rv-tdt/rvTDT "
        work.dir <- "cd /users/lgai/"
        data.dir <- "./cadd_10/peak_no_overlap/"
        data.file <- paste0("cadd_10_peak_window-size=76_snps_",i)
        results.dir <- paste0("./cadd_10/peak_no_overlap/cadd_10_peak_window-size=76_snps_",i)
        g.file <- paste0(data.file,".tped")
        p.file <- "./cadd_10/cadd_10.phen"
        #p.file <- paste0(data.file,".phen")
        m.file <- paste0(data.file,".map")
        
        #This works
        command<-paste0(rv.tdt.dir, results.dir, 
                        " -G ", data.dir, g.file,
                        " -P ", p.file,
                        " -M ", data.dir, m.file,
                        " --adapt 100 --alpha 0.00001 --permut 1000",
                        " -u 0.01"
        )
        
        command
        
        system(work.dir)
        system(command)     
}




#################################################################################







#This does not work.

#CADD filtered version
rv.tdt.dir<-"./rv-tdt/rvTDT "
work.dir <- "cd /users/lgai/"
data.dir <- "./cadd_10/"
data.file <- "cadd_10"
results.dir <- "cadd_10"
g.file <- paste0(data.file,".tped")
p.file <- paste0(data.file,".phen")
m.file <- paste0(data.file,".map")

command<-paste0(rv.tdt.dir, results.dir, 
                " -G ", data.dir, g.file,
                " -P ", data.dir, p.file,
                " -M ", data.dir, m.file,
                " --adapt 500 --alpha 0.00001 --permut 2000",
                " --lower_cutoff 0 --upper_cutoff 100",
                " --minVariants 3",
                " --maxMissRatio 1")
command

system(work.dir)
system(command)

################################################################################

#CADD filtered version to 400 SNPs, which does work
rv.tdt.dir<-"./rv-tdt/rvTDT "
work.dir <- "cd /users/lgai/"
data.dir <- "./cadd_10/"
data.file <- "cadd_10.400snps"
results.dir <- "cadd_10.400snps"
g.file <- paste0(data.file,".tped")
p.file <- "./cadd_10/cadd_10.phen"
#p.file <- paste0(data.file,".phen")
m.file <- paste0(data.file,".map")

#This works
command<-paste0(rv.tdt.dir, results.dir, 
                " -G ", data.dir, g.file,
                " -P ", p.file,
                " -M ", data.dir, m.file,
                " --adapt 100 --alpha 0.00001 --permut 1000",
                " -u 0.01"
)

command

system(work.dir)
system(command)

#This works as well
# command<-paste0(rv.tdt.dir, results.dir, 
#                 " -G ", data.dir, g.file,
#                 #               " -P ", data.dir, p.file,
#                 " -P ", "./cadd_10/cadd_10.phen",
#                 " -M ", data.dir, m.file,
#                 " --adapt 500 --alpha 0.00001 --permut 2000",
#                 " --lower_cutoff 0 --upper_cutoff 100",
#                 " --minVariants 3",
#                 " --maxMissRatio 1")
# 
# system(work.dir)
# system(command)

################################################################################

#CADD filtered version to 100 SNPs
rv.tdt.dir<-"./rv-tdt/rvTDT "
work.dir <- "cd /users/lgai/"
data.dir <- "./cadd_10/"
data.file <- "cadd_10.100snps"
results.dir <- "cadd_10.100snps"
g.file <- paste0(data.file,".tped")
p.file <- "./cadd_10/cadd_10.100snps.phen"
#p.file <- paste0(data.file,".phen")
m.file <- paste0(data.file,".map")

command<-paste0(rv.tdt.dir, results.dir, 
                " -G ", data.dir, g.file,
                #               " -P ", data.dir, p.file,
                " -P ", "./cadd_10/cadd_10.phen",
                " -M ", data.dir, m.file,
                " --adapt 500 --alpha 0.00001 --permut 2000",
                " --lower_cutoff 0 --upper_cutoff 100",
                " --minVariants 3",
                " --maxMissRatio 1")
command

system(work.dir)
system(command)

################################################################################

#CADD filtered version to 500 SNPs, which does not work
rv.tdt.dir<-"./rv-tdt/rvTDT "
work.dir <- "cd /users/lgai/"
data.dir <- "./cadd_10/"
data.file <- "cadd_10.500snps"
results.dir <- "cadd_10.500snps"
g.file <- paste0(data.file,".tped")
p.file <- "./cadd_10/cadd_10.phen"
#p.file <- paste0(data.file,".phen")
m.file <- paste0(data.file,".map")

command<-paste0(rv.tdt.dir, results.dir, 
                " -G ", data.dir, g.file,
                #                " -P ", p.file,
                " -P ", "./cadd_10/cadd_10.phen",
                " -M ", data.dir, m.file,
                " --adapt 100 --alpha 0.00001 --permut 1000",
                " -u 0.01"
)

command

system(work.dir)
system(command)

################################################################################
#Scratch

#Neither of the below commands work
command<-paste0(rv.tdt.dir, results.dir, 
                " -G ", data.dir, g.file,
                #               " -P ", data.dir, p.file,
                " -P ", "./cadd_10/cadd_10.phen",
                " -M ", data.dir, m.file,
                " --adapt 500 --alpha 0.00001 --permut 2000",
                " --lower_cutoff 0 --upper_cutoff 100",
                " --minVariants 3",
                " --maxMissRatio 1")

command
system(work.dir)
system(command)

command<-paste0(rv.tdt.dir, results.dir, 
                " -G ", data.dir, g.file,
                " -P ", p.file,
                " -M ", data.dir, m.file,
                " --adapt 100 --alpha 0.00001 --permut 1000",
                " -u 0.01"
)

system(paste0("./rvTDT ", reg, "-", k, "-real", currCadd, 
              " -G ", outPath, "/", reg, "_", k,"_tped_forRVTDTSim.txt -P ",
              outPath, "/", reg, "_", k, "_pheno_forRVTDTSim.txt -M ", outPath,
              "/", reg, "_", k, 
              "_map_forRVTDTSim.txt --adapt 100 --alpha 0.00001 --permut 1000 -u 0.01"))

################################################################################