### AMOVA in R ###

### set working directory containing all files needed ###
setwd("YOUR/DIRECTORY")

### load libraries (after installation of vcfR and poppr) ###
library(vcfR)
library(poppr)

### read in vcf used for calculating MDS plot ###
SPECIES_vcf <- read.vcfR("SPECIESmergeDPmm22homosnpsstar.vcf", verbose = FALSE)

### convert vcf to genlight format ###
SPECIES_genlight <- vcfR2genlight(SPECIES_vcf)

### make table of strata in working directory like this ###
### population	lineage ###
### Hainich	1 ###
### Hainich	1 ###
### Hainich	2 ###
### Kwald	2 ###
### ...	... ###

### read in table of strata ###
SPECIES.pop.data <- read.table("TABLE", sep = "\t", header = TRUE)

### define strata in the genlight according to the table and set population ###
strata(SPECIES_genlight) <- SPECIES.pop.data
setPop(SPECIES_genlight) <- ~population

### do amova using population as hierarchy (multiple nested levels are possible e.g. ~population/subpopulation etc.) ###
SPECIES_amova <- poppr.amova(SPECIES_genlight, ~population)

### show amova results ###
SPECIES_amova

### check for significance of variances at different hierarchical levels ###
randtest(SPECIES_amova, nrepet = 10000)

### see https://grunwaldlab.github.io/Population_Genetics_in_R/AMOVA.html for more ###
