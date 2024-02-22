
#######################################################################
### Matching Intogen driver genes with CN state and RNA expression 
######################################################################
# This script plots the correlations between copy number states of Intogen driver genes and RNA expression 
# Bulk RNA expression of cancer driver genes is adjusted for the tumour cell percentage estimated by a board-certified pathologist

library(data.table)
library(dplyr)
library(rstatix)
require(tidyverse)

# prepare output dirs:
fig_dir = "~/analysis/multi_omic/Intogen_driver/plots"
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)

# My theme 
theme_mypub <- function(base_size = 14,
                        base_family = ""){
  theme_bw(base_size = base_size, 
           base_family = base_family) %+replace%
    theme(
      panel.grid.major = element_blank(),   
      #panel.grid.minor = element_blank(),   
      
      complete = TRUE
    )
}

my_theme <- theme(
  plot.title = element_text(color="black", size=16, face="bold"),
  axis.title.x = element_text(color="black", size=16, face="bold"),
  axis.text.x = element_text(color="black", size=14, face="bold"),
  axis.title.y = element_text(color="black", size=16, face="bold"),
  axis.text.y = element_text(color="black", size=14, face="bold"),
  legend.title = element_text(color = "black", size = 14, face="bold"),
  legend.text = element_text(color = "black", size = 14, face="bold"),
  
)


## Read in Overview of matching RNA and WES data
matching_overview <- read.table("~/analysis/General_files/MEMORI_specific/Matching_RNA_DNA_samples.txt", header = T, stringsAsFactors = F, sep = "\t")
matching_overview <- matching_overview %>% filter_at(vars(c("RNA_file", "WES_file")), all_vars(!is.na(.)))  

# Read in Intogen list with driver genes form Esophageal Cancer and gastric cancer 
# https://www.intogen.org/ 
Intogen_EC_GC_list  <- read.table("~/analysis/General_files/external_files/Intogen_drivers.txt", stringsAsFactors = F, header = T)

# ENSEMBLE/ HUGO Gene ID and gene location list
ENSEMBLE_list<- read.table("~/analysis/General_files/external_files/compiledGeneInfo.txt", stringsAsFactors = F, header = T)
chr_filt <- paste0("chr", c(1:23, "X", "Y"))
ENSEMBLE_list <- ENSEMBLE_list[which(ENSEMBLE_list$Chr %in% chr_filt),]

# Get ENSEMBL Gene Names and gene location for Oncogenes/TSG from Intogen list
Intogen_loc <- merge(Intogen_EC_GC_list, ENSEMBLE_list, by.x = "Symbol", by.y = "Name")
Intogen_loc_reduced <- Intogen_loc[ ,c("Symbol", "GeneID")]

# Create Mastertable of CNA per driver 
  files_CNA_Intogen <- list.files(path = "~/analysis/WES/CNA/datasets/CNA_Intogen_driver/", pattern = "_CNA_Intogen.csv", full.names = T) 
                                  
  for (fn in files_CNA_Intogen) { 
    
    genes <- read.csv(files_CNA_Intogen[1], header=TRUE,stringsAsFactors = F)[,1]
    samples <- gsub("\\S+(MEMORI_\\S+D\\d)\\S+","\\1", files_CNA_Intogen)
    df    <- do.call(cbind,lapply(files_CNA_Intogen,function(fn)read.csv(fn, header = TRUE,stringsAsFactors = F)[,2]))
    MasterTable_detailed    <- as.data.frame(cbind(genes,df),stringsAsFactors=F)
    colnames(MasterTable_detailed)  <- c("Symbol",samples)
  }
  
# Select only samples, that have matching RNA
  MasterTable_detailed <- merge(MasterTable_detailed, Intogen_loc_reduced, by = "Symbol")
  columns_to_select <- c("GeneID", matching_overview$WES_file)
  MasterTable_detailed <- MasterTable_detailed[, c(names(MasterTable_detailed) %in% columns_to_select)]
  rownames(MasterTable_detailed) <- MasterTable_detailed$GeneID; MasterTable_detailed <- MasterTable_detailed[ , -which(names(MasterTable_detailed) %in% c("GeneID"))]
  MasterTable_detailed <- as.data.frame(t(MasterTable_detailed))
  
  # Read in Gene expression table 
  columns_RNA_to_select <- c("Gene_ENSEMBLE", matching_overview$RNA_file)
  GeneIndex_MasterTable1.5M_norm <- read.table("~/analysis/RNA_Seq/input_files/GeneExpression_Mastertables/GeneIndex_MasterTable1.5M_norm_batch_1_2_3.txt_all_genes.txt", header = T, stringsAsFactors = F)
  GeneIndex_MasterTable1.5M_norm_only_drivers <- GeneIndex_MasterTable1.5M_norm[which(GeneIndex_MasterTable1.5M_norm$Gene_ENSEMBLE %in% unique(Intogen_loc_reduced$GeneID)), ]
  GeneIndex_MasterTable1.5M_norm_only_drivers <- GeneIndex_MasterTable1.5M_norm_only_drivers[,names(GeneIndex_MasterTable1.5M_norm) %in%  columns_RNA_to_select]
  GeneIndex_MasterTable1.5M_norm_only_drivers<- as.data.frame(t(GeneIndex_MasterTable1.5M_norm_only_drivers))
  colnames(GeneIndex_MasterTable1.5M_norm_only_drivers) <- GeneIndex_MasterTable1.5M_norm_only_drivers[1,]; GeneIndex_MasterTable1.5M_norm_only_drivers <- GeneIndex_MasterTable1.5M_norm_only_drivers[-1,]
  
  driver_genes <- unique(Intogen_loc_reduced$GeneID)
  list_corr_coef_RNA_cell_adj <- list() 
  list_corr_coef_RNA_NO_adj <- list() 
  
  ############################################
  # Plot correlation plots
  
  for (g in driver_genes){
    gene_name <- Intogen_loc_reduced[which(Intogen_loc_reduced$GeneID == g), 1]
    single_driver_RNA_df <-  subset(GeneIndex_MasterTable1.5M_norm_only_drivers, select=g)
    colnames(single_driver_RNA_df) <- "RNA_expression"
    single_driver_RNA_df$RNA_expression <- as.numeric(single_driver_RNA_df$RNA_expression)
    
    single_driver_CNA_df <-  subset(MasterTable_detailed, select=g)
    colnames(single_driver_CNA_df) <- "CNt"
    single_driver_CNA_df$CNt <- as.numeric(single_driver_CNA_df$CNt)
   
     Merged_RNA_CNt <- cbind(single_driver_RNA_df, single_driver_CNA_df)
     Merged_RNA_CNt$SampleID <- rownames(Merged_RNA_CNt )
     
     # Adjust RNA expression of Intogen driver genes for tumour cell percentage (estimated by pathologist)
     RNA_pathology_cellularity <- read.table("~/analysis/General_files/MEMORI_specific/RNA_samples_pathologist_review.txt", sep = "\t", header = T, stringsAsFactors = F)
     RNA_pathology_cellularity$Patho.percentage <- ifelse((RNA_pathology_cellularity$Patho.percentage == "<0.1" | RNA_pathology_cellularity$Patho.percentage == ""), NA, RNA_pathology_cellularity$Patho.percentage)
     Merged_RNA_CNt_patho <- merge(Merged_RNA_CNt , RNA_pathology_cellularity, by.x = "SampleID", by.y = "SampleID")
     Merged_RNA_CNt_patho$Patho.percentage <- as.numeric(Merged_RNA_CNt_patho$Patho.percentage); Merged_RNA_CNt_patho$CNt <- as.numeric(Merged_RNA_CNt_patho$CNt)
     Merged_RNA_CNt_patho$Cellul_adapted_RNA_expression <- Merged_RNA_CNt_patho$RNA_expression/Merged_RNA_CNt_patho$Patho.percentage
     Merged_RNA_CNt_patho$CNt
     
     pdf(paste0(fig_dir, "/", gene_name, "_RNA_CNt_corr_plot.pdf"), height = 5, width = 5)
     p <- ggplot(Merged_RNA_CNt_patho, aes(x=CNt, y=Cellul_adapted_RNA_expression)) +
       geom_point(size=2, shape=23) + geom_smooth() + stat_cor(method = "pearson", size = 8) +
       labs(x="CNS", y="Adjusted RNA expression [cpm]") + ggtitle(paste0(gene_name, " correlation CNS and RNA"))  + theme_mypub() + my_theme
     plot(p)
     dev.off()
  
    
     }

  
