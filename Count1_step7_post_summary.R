##################################################################
## (9)                                                          ##
## Code for running at the cluster (should be similar to        ##
## Count1_step7_post_summary.ipynb)                             ##
## Get p values for each PheCode                                ##
## should be similar to Count1_step6_post_summary_test.ipynb    ##
##################################################################


# PheWAS results in: gs://fc-secure-fa05fe0c-78d6-4ab2-843a-e14e9400d98a/result/
# Save summary to: gs://fc-secure-fa05fe0c-78d6-4ab2-843a-e14e9400d98a/result/PheWAS_summary
# Groupfile:  
# paste0("gs://fc-secure-fa05fe0c-78d6-4ab2-843a-e14e9400d98a/data/annotation/vep109/grouping/hclof_noflag_missense0.8_7tools_POPMAX0.001/AoU_250K_exome_annotation_vep109_chr",
# 1:22,".vcf.gz.hclof_noflag_missense0.8_7tools_POPMAX0.001.RData")


args <- (commandArgs(TRUE))
pheno_i     <- as.character(args[1])
study_pathi <- as.character(args[2])   # PheWAS results file
chr_i       <- as.numeric(args[3])
group_i     <- as.character(args[4])
outfolder_i <- as.character(args[5])


library(GENESIS)
library(CompQuadForm)
library(GMMAT)
library(survey)


print(pheno_i)
print(paste0("Calculating for chr ",chr_i))


if(!dir.exists(outfolder_i) ){ 
  dir.create( outfolder_i, recursive=TRUE)
}


    
## summarizing 
   source("transcript_meta_analysis_v2.R")
   pvals_i_ver2 <- transcript_meta_analysis_ver2(study_path_vector=study_pathi, 
                                                 grouping_path_vector=group_i,
                                                 test=c("Burden", "SKAT", "SMMAT")) 
   save(pvals_i_ver2, file = paste0(outfolder_i, pheno_i, "_pvals_chr", chr_i, "_ver2.rdata")  )  

  
  
   
