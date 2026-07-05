#####################################################################
# Ancestry-Informed Sliding window analysis using STAARpipeline
# Wenbo Wang, Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 07/05/2026
#####################################################################
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
## aGDS directory
agds_dir <- get(load("/path_to_the_file/agds_dir.Rdata"))
## Null model
obj_nullmodel <- get(load("/path_to_the_file/obj_nullmodel.Rdata"))

## QC_label
QC_label <- "annotation/filter"
## variant_type
variant_type <- "SNV"
## geno_missing_imputation
geno_missing_imputation <- "mean"
## sliding_window_length
sliding_window_length <- 2000

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("/path_to_the_file/Annotation_name_catalog.Rdata"))
# Or equivalently
# Annotation_name_catalog <- read.csv("/path_to_the_file/Annotation_name_catalog.csv")
## Use_annotation_weights
Use_annotation_weights <- TRUE
## Annotation name
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")

## output path
output_path <- "/path_to_the_output_file/"
## output file name
output_file_name <- "TOPMed_F8_Sliding_Window_AI_STAAR"
## input chromosome from batch file
chr <- as.numeric(commandArgs(TRUE)[1])

#sliding window dataset - with columns for "Chr", "Start Loc", "End Loc"
sliding_window_results <- get(load("/path_to_the_file/sliding_window_results.Rdata"))
#or, use alternative functions to read data in other formats

###########################################################
#           Main Function 
###########################################################
## aGDS file
agds.path <- agds_dir[chr]
genofile <- seqOpen(agds.path)
sliding_window_results_chr <- sliding_window_results[unlist(sliding_window_results[,"Chr"]) == chr,]

AI_results_sliding_window <- c()

if(!is.null(nrow(sliding_window_results_chr))){
  for(i in 1:nrow(sliding_window_results_chr)){
    results <- c()
    
    start_loc <- unlist(sliding_window_results_chr[,"Start Loc"])
    end_loc <- unlist(sliding_window_results_chr[,"End Loc"])
    start_loc_sub <- start_loc[i]
    end_loc_sub <- end_loc[i]
    
    results <- try(Sliding_Window(chr=chr,start_loc=start_loc_sub,end_loc=end_loc_sub,
                                  sliding_window_length=sliding_window_length,type="multiple",
                                  genofile=genofile,obj_nullmodel=obj_nullmodel,
                                  rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                  QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                  Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                  Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name))
    
    if(class(results)[1]!="try-error")
    {
      AI_results_sliding_window <- c(AI_results_sliding_window,results)
    }
  }
}else{
  start_loc <- unlist(sliding_window_results_chr["Start Loc"])
  end_loc <- unlist(sliding_window_results_chr["End Loc"])
  AI_results_sliding_window <- try(Sliding_Window(chr=chr,start_loc=start_loc,end_loc=end_loc,
                                               sliding_window_length=sliding_window_length,type="multiple",
                                               genofile=genofile,obj_nullmodel=obj_nullmodel,
                                               rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                               QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                               Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                               Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name))
}

save(AI_results_sliding_window,file=paste0(output_path,output_file_name,"_",chr,".Rdata"))

seqClose(genofile)
