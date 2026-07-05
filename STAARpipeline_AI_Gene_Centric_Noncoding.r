########################################################################
# Ancestry-Informed Gene-centric Noncoding analysis using STAARpipeline
# Wenbo Wang, Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 07/03/2026
########################################################################
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

## Annotation_dir
Annotation_dir <- "annotation/info/TOPMedAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("/path_to_the_file/Annotation_name_catalog.Rdata"))
# Or equivalently
# Annotation_name_catalog <- read.csv("/path_to_the_file/Annotation_name_catalog.csv")
## Use_annotation_weights
Use_annotation_weights <- TRUE
## Annotation name
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")

#set gene chromosome
chr <- 1
#set gene name
gene_name <- "GeneName"

agds.path <- agds_dir[chr]
genofile <- seqOpen(agds.path)

## output path
output_path <- "/path_to_the_output_file/"
## output file name
output_file_name <- "TOPMed_F8_GC_Noncoding_AI_STAAR"

#set coding functional category
category <- "all_categories"

###########################################################
#           Main Function 
###########################################################

start_time <- Sys.time()
results_noncoding <- Gene_Centric_Noncoding(chr, gene_name, genofile, obj_nullmodel, category = category,
                                            rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                            QC_label = QC_label, variant_type = variant_type,
                                            use_ancestry_informed=TRUE,find_weight=TRUE,
                                            geno_missing_imputation = geno_missing_imputation,
                                            Annotation_dir = Annotation_dir, Annotation_name_catalog = Annotation_name_catalog,
                                            Use_annotation_weights = Use_annotation_weights,
                                            Annotation_name = Annotation_name)
end_time <- Sys.time()
end_time - start_time

save(results_noncoding,file=paste0(output_path, output_file_name, "_", chr, "_", gene_name, "_", category,".Rdata"))

seqClose(genofile)

