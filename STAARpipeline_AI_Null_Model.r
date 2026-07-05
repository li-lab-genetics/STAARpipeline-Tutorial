##################################################################
# fit ancestry-informed STAAR null model for single-trait analysis
# Wenbo Wang, Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 07/05/2026
##################################################################

rm(list=ls())
gc()

library(GENESIS)
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)

###########################################################
#           User Input
###########################################################
## Phenotype file
phenotype <- read.csv("/path_to_the_file/pheno.csv")
## (sparse) GRM file
sgrm <- get(load("/path_to_the_file/sGRM.Rdata"))
## file directory for the output file 
output_path <- "/path_to_the_output_file/"
## output file name
output_name <- "obj_nullmodel_AI.Rdata"

###########################################################
#           Main Function 
###########################################################

obj_nullmodel <- fit_nullmodel(trait.norm~age+age2+male+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+as.factor(study_ancestry),
                               data=phenotype,kins=sgrm,use_sparse=TRUE,kins_cutoff=0.022,id="sample.id",
                               groups="ancestry",pop.groups=phenotype$ancestry,family=gaussian(link = "identity"),verbose=TRUE,
                               B = 500)

save(obj_nullmodel,file=paste0(output_path,output_name))

