This folder contains files for data cleaning. None of the files need to be run more than once, except 5_filtering_by_annotation.R if you want to change the filtering paramters.

####################################

1_clean_vcf.R

Converts missing values into the correct format for BEAGLE
Removes tri-allelic SNPs/duplicate sites
Removes individuals with Mendelian errors

DOES NOT: remove half-calls (1/., ./0). These will be phased.

Creates table of SNP names and positions.

####################################

2_clean_ped.R

Edits pedfile IDs to match vcf file IDs.

####################################

3_clean_annovar.R

Reads in 8q24 ANNOVAR and subsets it to "useful" variables.

####################################

4_BEAGLE_phasing.R

####################################

5_filtering.R

Filters the phased file in 3 ways:

1. Filtered by functional annotation information (CADD/GWAVA/EIGEN) using cut-offs (1048 SNPs)
2. Filtered to peak TDT signal (679 SNPs)
3. Both 1 and 2 (119 SNPs)

and outputs a trio geno txt file (for RV-TDT, rvTDT) and a vcf file (for ScanTrio).