#####################################################################
# Ancestry-Informed Gene-centric ncRNA analysis using STAARpipeline
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
output_file_name <- "TOPMed_F8_ncRNA_AI_STAAR"

#ncRNA dataset - with columns for "Gene name", "Chr"
ncRNA_results <- get(load("/path_to_the_file/results_ncRNA.Rdata"))
#or, use alternative functions to read data in other formats

###########################################################
#           Main Function 
###########################################################

results_ncRNA <- c()

if(!is.null(nrow(ncRNA_results))){ #multiple ncRNAs
  gene_name <- as.character(unlist(ncRNA_results[,"Gene name"]))
  chr <- as.integer(unlist(ncRNA_results[,"Chr"]))
  
  for(i in 1:nrow(ncRNA_results)){
    gene_name_sub <- gene_name[i]
    chr_sub <- chr[i]
    
    agds.path <- agds_dir[chr_sub]
    genofile <- seqOpen(agds.path)
    
    results <- c()
    results <- try(ncRNA(chr=chr_sub,gene_name=gene_name_sub,genofile=genofile,obj_nullmodel=obj_nullmodel,
                         rare_maf_cutoff=0.01,rv_num_cutoff=2,
                         QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                         Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                         use_ancestry_informed=TRUE,find_weight=TRUE))
    seqClose(genofile)
    results_ncRNA <- c(results_ncRNA,results)
  }
}else if(is.null(length(ncRNA_results))){
  results_ncRNA <- c()
}else{
  gene_name <- as.character(unlist(ncRNA_results["Gene name"]))
  chr <- as.integer(unlist(ncRNA_results["Chr"]))
  
  agds.path <- agds_dir[chr]
  genofile <- seqOpen(agds.path)
  
  results_ncRNA <- try(ncRNA(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                             rare_maf_cutoff=0.01,rv_num_cutoff=2,
                             QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                             Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                             Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                             use_ancestry_informed=TRUE,find_weight=TRUE))
  
}

save(results_ncRNA,file=paste0(output_path,output_file_name,".Rdata"))
