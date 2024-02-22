
######################################################################
### Expression of genes coding for neoantigens during treatment
######################################################################
# This script 
# 1. merges neoantigenic SNVs of each sample with gene expression data of the respective gene that harbours the neoantigen
# 2. calculates the mean expression of genes coding for neoantigens
# 3. Plots the mean expression of genes coding for neoantigens during treatment

library(data.table)
library(dplyr)
library(rstatix)
library(ggplot2)

# prepare output dirs:
data_dir = "~/analysis/multi_omic/neoantigen_expression/tables/mean_RNA_expression_of_neoantigens"
dir.create(data_dir, showWarnings=FALSE, recursive=TRUE)
fig_dir = "~/analysis/multi_omic/neoantigen_expression/plots"
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)

# My theme 
my_theme <- theme(
  plot.title = element_text(color="black", size=16, face="bold"),
  axis.title.x = element_text(color="black", size=16, face="bold"),
  axis.text.x = element_text(color="black", size=16, face="bold"),
  axis.title.y = element_text(color="black", size=16, face="bold"),
  axis.text.y = element_text(color="black", size=16, face="bold"),
  legend.title = element_text(color = "black", size = 14, face="bold"),
  legend.text = element_text(color = "black", size = 14, face="bold"),
  
)

## Read in file indicating file names of matching RNA-Seq and WES data
matching_overview <- read.table("~/analysis/General_files/MEMORI_specific/Matching_RNA_DNA_samples.txt", header = T, stringsAsFactors = F, sep = "\t")
matching_overview <- matching_overview %>% filter_at(vars(c("RNA_file", "WES_file")), all_vars(!is.na(.)))  

list_neoant_expr <- list()

for(i in 1:nrow(matching_overview)){
  
## Prepare filtered neoantigen file  
DNA_SAM <- matching_overview[i, 3] 
single_sample_neonatigen_df <- read.table(paste0("~/analysis/WES/neoantigens/filtering_steps/tables/neoSNV_tables/",DNA_SAM, ".weakly_filtered_neoSNVs"), header = T, stringsAsFactors = F, sep ="\t")
single_sample_neonatigen_df$start <- single_sample_neonatigen_df$Start
single_sample_neonatigen_df$end <- single_sample_neonatigen_df$Start
single_sample_neonatigen_df$Sample <- sapply(single_sample_neonatigen_df$Sample, function(x) substr(x,1, 23))
single_sample_neonatigen_df <- single_sample_neonatigen_df[,c("Sample", "LineID", "Ref", "Alt", "peptide", "Score", "Affinity", "BindLevel", "Chr", "start", "end")]
single_sample_neonatigen_df <- as.data.table(single_sample_neonatigen_df)

## Prepare normalized RNA expression file
RNA_SAM <- matching_overview[i, 2]
GeneIndex_MasterTable1.5M_norm <- read.table("~/analysis/RNA_Seq/input_files/GeneExpression_Mastertables/GeneIndex_MasterTable1.5M_norm_batch_1_2_3.txt_all_genes.txt", header = T, stringsAsFactors = F)
Gene_info_table <- read.table("~/analysis/General_files/external_files/compiledGeneInfo.txt", header = T, stringsAsFactors = F)
Gene_info_table <- Gene_info_table[,c("GeneID", "Chr", "Start", "End")]
df <- merge(GeneIndex_MasterTable1.5M_norm, Gene_info_table, by.x = "Gene_ENSEMBLE",  by.y = "GeneID")
colnames(df)[colnames(df) == "Start"] <- "start"
colnames(df)[colnames(df) == "End"] <- "end"

## Select requested columns and reoder them in correct order for foverlaps function
single_sample_RNA_df <-  as.data.table(df[,c(RNA_SAM, "Gene_ENSEMBLE", "Chr", "start", "end")])

# Matching neoantigenic mutations with the respective gene they occur in
setkey(single_sample_RNA_df, Chr, start, end)
z <- foverlaps(single_sample_neonatigen_df, single_sample_RNA_df, type="any")

# Calculation of mean expression of strong binding (SB), weak binding (WB) and all neoantigens
z_Total_neoantigens <- z[which(z$BindLevel == "SB" | z$BindLevel == "WB"), ]
z_SB_neoantigens <- z[which(z$BindLevel == "SB"), ]
z_WB_neoantigens <- z[which(z$BindLevel == "WB"), ]
mean_exp_WB_neoantigens <- colMeans(z_WB_neoantigens %>% dplyr::select(starts_with("MEME")), na.rm = TRUE)
mean_exp_SB_neoantigens <- colMeans(z_SB_neoantigens %>% dplyr::select(starts_with("MEMEX")), na.rm = TRUE)
mean_exp_Total_neoantigens <- colMeans(z_Total_neoantigens  %>% dplyr::select(starts_with("MEMEX")), na.rm = TRUE)
Mean_exp_neoantigen_df <- data.frame("SampleID" = DNA_SAM, "mean_exp_Total_neoantigens" = mean_exp_Total_neoantigens, "mean_exp_SB_neoantigens"=  mean_exp_SB_neoantigens, "mean_exp_WB_neoantigens" = mean_exp_WB_neoantigens)
list_neoant_expr[[DNA_SAM]] <- Mean_exp_neoantigen_df

# Export files 
write.table(Mean_exp_neoantigen_df, file = paste0(data_dir, "/",DNA_SAM, ".mean.exp.neoant.txt"), row.names = F, quote = F)
} 

# Create Mastertable with mean expression of all neoantigens per sample
Neo_SNV_MasterTable <- rbindlist(list_neoant_expr)

# Add clinical data
Neo_SNV_MasterTable$patient <- sapply(Neo_SNV_MasterTable$SampleID, function(x) substr(x,8,13))
Neo_SNV_MasterTable$timepoint <- sapply(Neo_SNV_MasterTable$SampleID, function(x) substr(x,14,14)) 
Neo_SNV_MasterTable$Responsiveness <- sapply(Neo_SNV_MasterTable$SampleID, function(x) substr(x,8,9)) 
Neo_SNV_MasterTable$mean_exp_Total_neoantigens <- as.numeric(Neo_SNV_MasterTable$mean_exp_Total_neoantigens)


## Plot neoantigen expression during treatment
my_comparisons <- list(c("A", "B"), c("A", "C"), c("B", "C"))
pdf(paste0(fig_dir, "/Total_neoantigen_expr_per_timepoint.pdf"), height = 5, width = 5)
p <- ggplot2.violinplot(data=Neo_SNV_MasterTable, xName="timepoint", yName="mean_exp_Total_neoantigens",
                        groupName="timepoint",
                        groupColors= alpha(c("khaki2", '#cac3df', '#02818a'), c(0.9, 1, 0.8)), showLegend=FALSE,
                        backgroundColor="white", xtitle="", ytitle="Neoantigen mean expression [cpm]", 
                        mainTitle="Expression of Neoantigens during treatment",
                        addDot=T, dotSize=0.5,
                        removePanelGrid=F,removePanelBorder=F,
                        axisLine=c(0.5, "solid", "black")) + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size = 5, face="bold") + ylim(c(-0.01,4.7))
q <- p + my_theme + theme(plot.title = element_text(color = "black", size = 13, face="bold"))
plot(q)
dev.off()
