###########################################################
# Ancestry-Informed Individual analysis using STAARpipeline
# Wenbo Wang, Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 07/03/2026
###########################################################
rm(list=ls())
gc()

## load required packages
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)

###########################################################
#           User Input
###########################################################
jobs_num <- get(load("/path_to_the_file/jobs_num.Rdata"))
## aGDS directory
agds_dir <- get(load("/path_to_the_file/agds_dir.Rdata"))
## Null model
obj_nullmodel <- get(load("/path_to_the_file/obj_nullmodel.Rdata"))

## QC_label
QC_label <- "annotation/filter"
## variant_type
variant_type <- "variant"
## geno_missing_imputation
geno_missing_imputation <- "mean"

## output path
output_path <- "/path_to_the_output_file/"
## output file name
output_file_name <- "TOPMed_F8_Individual_Analysis_AI_GWA"
## input array id from batch file
arrayid <- as.numeric(commandArgs(TRUE)[1])

#individual analysis results
individual_results <- get(load("/path_to_the_file/results_individual_analysis.Rdata"))
#or, use alternative functions to read data in other formats

###########################################################
#           Main Function 
###########################################################

chr <- which.max(arrayid <= cumsum(jobs_num$individual_analysis_num))
group.num <- jobs_num$individual_analysis_num[chr]

if (chr == 1){
  groupid <- arrayid
}else{
  groupid <- arrayid - cumsum(jobs_num$individual_analysis_num)[chr-1]
}

## aGDS file
agds.path <- agds_dir[chr]
genofile <- seqOpen(agds.path)

start_loc <- (groupid-1)*(10e6) + jobs_num$start_loc[chr]
end_loc <- start_loc + (10e6) - 1
end_loc <- min(end_loc,jobs_num$end_loc[chr])

start_time <- Sys.time()
results_individual_analysis <- Individual_Analysis(chr=chr, individual_results = individual_results, 
                                                   start_loc=start_loc,end_loc=end_loc,
                                                   genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                   QC_label=QC_label,variant_type=variant_type,
                                                   geno_missing_imputation=geno_missing_imputation,
                                                   use_ancestry_informed = TRUE, find_weight = TRUE)


end_time <- Sys.time()
end_time - start_time

save(results_individual_analysis,file=paste0(output_path,output_file_name,"_",chr,"_",groupid,".Rdata"))

seqClose(genofile)