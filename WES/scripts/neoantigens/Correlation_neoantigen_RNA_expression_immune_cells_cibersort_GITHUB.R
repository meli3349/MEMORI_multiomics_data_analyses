
####################################################################################################################
# Correlation of RNA expression of neoantigens and infiltration of immune cells deconvoluted via CIBERSORT
####################################################################################################################
# This script...
# 1. classifies samples into high immune infiltrated or low immune infiltrated using the mean of respective immune cell score
# 2. plots the mean neoantigen expression in high and low infiltrated groups for several immune cell types 

library(dplyr)
library(corrplot)
library(ggplot2)
library(easyGgplot2)
library(ggpubr)
library(rstatix)

# prepare output dirs:
fig_dir = "~/analysis/multi_omic/neoantigen_expression/plots/Correlation_Cibersort_neoantigen_expression"
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)

# My theme
my_theme <- theme(
  plot.title = element_text(color="black", size=16, face="bold"),
  axis.title.x = element_text(color="black", size=16, face="bold"),
  axis.text.x = element_text(color="black", size=14, face="bold"),
  axis.title.y = element_text(color="black", size=16, face="bold"),
  axis.text.y = element_text(color="black", size=14, face="bold"),
  legend.title = element_text(color = "black", size = 14, face="bold"),
  legend.text = element_text(color = "black", size = 14, face="bold"),
  
)

## Read in file indicating file names of matching RNA-Seq and WES data
matching_overview <- read.table("~/analysis/General_files/MEMORI_specific/Matching_RNA_DNA_samples.txt", header = T, stringsAsFactors = F, sep = "\t")
matching_overview <- matching_overview %>% filter_at(vars(c("RNA_file", "WES_file")), all_vars(!is.na(.)))  
samples <- matching_overview$WES_file

# Create Mastertable with mean expression of all neoantigens per sample
Neo_SNV_MasterTable <- do.call(rbind,lapply(samples,function(sam)read.table(paste0("~/analysis/multi_omic/neoantigen_expression/tables/mean_RNA_expression_of_neoantigens/", sam, ".mean.exp.neoant.txt"), header = TRUE,stringsAsFactors = F)))

# Read in Cibersort table with deconvoluted immune cell types
CibersortTable <- read.delim("~/analysis/RNA_Seq/input_files/Cibersort_tables/CIBERSORT.Output_Abs_Job6.txt")

# Adding absolute cell scores of different cell subtype a score for main cell types (e.g. Macrophages M0, M1 and M2 to Macrophages)
CibersortTable$TotalMacrophages <- rowSums(CibersortTable[,c(15:17)])
CibersortTable$TotalNKcells <- rowSums(CibersortTable[,c(12:13)])
CibersortTable$TotalMonocytes <- CibersortTable[,14]
CibersortTable$TotalDendricCells <- rowSums(CibersortTable[,c(18:19)])
CibersortTable$TotalMastCells <- rowSums(CibersortTable[,c(20:21)])
CibersortTable$TotalEosinophils <- CibersortTable[,22]
CibersortTable$TotalNeutrophils <- CibersortTable[,23]
CibersortTable$OtherBcells <- rowSums(CibersortTable[,c(2:4)])
CibersortTable$OtherTcells <- rowSums(CibersortTable[,c(9:11)])
CibersortTable$TotalCD4cells <- rowSums(CibersortTable[,c(6:8)])
CibersortTable$TotalCD8cells <- CibersortTable[,5]

# Selecting immune cell types of interest 
CibersortTablePlot <- CibersortTable[,c(1,27:37)]
colnames(CibersortTablePlot) <- c("patient", "Macrophages", "NK cells", "Monocytes", "Dendric cells", "Mast cells", "Eosinphils", "Neutrophils", "Other B-cells", "Other T-cells", "CD4 T-cells", "CD8 T-cells")

# Merging Cibersort data with neoantigen expression data
y <- merge(matching_overview, CibersortTablePlot, by.x = "RNA_file", by.y = "patient")
merged_Neo_SNV_Cibersort <- merge(y, Neo_SNV_MasterTable, by.x = "WES_file", by.y = "SampleID")

# Classification of samples into low and high infiltrated samples using the mean of individual immune cell types
merged_Neo_SNV_Cibersort$CD8_score <- if_else(merged_Neo_SNV_Cibersort$`CD8 T-cells` > mean(merged_Neo_SNV_Cibersort$`CD8 T-cells`), "high", "low")
merged_Neo_SNV_Cibersort$CD8_score <- factor(merged_Neo_SNV_Cibersort$CD8_score, levels = c("low", "high"))
merged_Neo_SNV_Cibersort$CD4_score <- if_else(merged_Neo_SNV_Cibersort$`CD4 T-cells` > mean(merged_Neo_SNV_Cibersort$`CD4 T-cells`), "high", "low")
merged_Neo_SNV_Cibersort$CD4_score <- factor(merged_Neo_SNV_Cibersort$CD4_score, levels = c("low", "high"))
merged_Neo_SNV_Cibersort$NK_score <- if_else(merged_Neo_SNV_Cibersort$`NK cells` > mean(merged_Neo_SNV_Cibersort$`NK cells`), "high", "low")
merged_Neo_SNV_Cibersort$NK_score <- factor(merged_Neo_SNV_Cibersort$NK_score, levels = c("low", "high"))
merged_Neo_SNV_Cibersort$Macrophages_score <- if_else(merged_Neo_SNV_Cibersort$Macrophages > mean(merged_Neo_SNV_Cibersort$Macrophages), "high", "low")
merged_Neo_SNV_Cibersort$Macrophages_score <- factor(merged_Neo_SNV_Cibersort$Macrophages_score, levels = c("low", "high"))


# Plotting neoantigen mean expression in low and high infiltrated CD8 cell groups
pdf(paste0(fig_dir, "/Total_Neoantigen_RNA_expr_CD8_Cibersort_score.pdf"), width = 5, height = 5)
compar_CD8 <- list(c("high", "low"))
p <- ggplot2.violinplot(data=merged_Neo_SNV_Cibersort, xName="CD8_score",yName="mean_exp_Total_neoantigens", #fill = "Clonality",
                        groupName="CD8_score",
                        groupColors= alpha(c("darkseagreen3","#6c8f6c"), c(0.8,1)), showLegend=FALSE,
                        backgroundColor="white", xtitle="CD8 T-cell score", ytitle="Neoantigen Mean expression", 
                        mainTitle="",
                        addDot=T, dotSize=0.5,
                        removePanelGrid=F,removePanelBorder=F,
                        axisLine=c(0.5, "solid", "black"),
                        legendPosition="bottom") + ylim(c(-0.1,4.2)) +
  stat_compare_means(comparisons = compar_CD8, method = "t.test", size = 6, face="bold") + my_theme 
plot(p)
dev.off()


# Plotting neoantigen mean expression in low and high infiltrated CD4 cell groups
pdf(paste0(fig_dir, "/Total_Neoantigen_RNA_expr_CD4_Cibersort_score.pdf"), width = 5, height = 5)
compar_CD4 <- list(c("high", "low"))
q <- ggplot2.violinplot(data=merged_Neo_SNV_Cibersort, xName="CD4_score",yName="mean_exp_Total_neoantigens", #fill = "Clonality",
                        groupName="CD4_score",
                        groupColors= alpha(c("darkseagreen3","#6c8f6c"), c(0.8,1)), showLegend=FALSE,
                        backgroundColor="white", xtitle="CD4 T-cell score", ytitle="Neoantigen Mean expression", 
                        mainTitle="",
                        addDot=T, dotSize=0.5,
                        removePanelGrid=F,removePanelBorder=F,
                        axisLine=c(0.5, "solid", "black"),
                        legendPosition="bottom") + ylim(c(-0.1,4.2)) +
  stat_compare_means(comparisons = compar_CD4, method = "t.test", size = 6, face="bold") + my_theme
plot(q) 
dev.off()

# Plotting neoantigen mean expression in low and high infiltrated NK cell groups
pdf(paste0(fig_dir, "/Total_Neoantigen_RNA_expr_NKcells_Cibersort_score.pdf"), width = 5, height = 5)
compar_NK<- list(c("high", "low"))
r <- ggplot2.violinplot(data=merged_Neo_SNV_Cibersort, xName="NK_score",yName="mean_exp_Total_neoantigens", #fill = "Clonality",
                        groupName="NK_score",
                        groupColors= alpha(c("#2166AC","#B2182B"), .8), showLegend=FALSE,
                        backgroundColor="white", xtitle="NK cell score", ytitle="Neoantigen Mean expression", 
                        mainTitle="",
                        addDot=T, dotSize=0.5,
                        legendPosition="bottom") + ylim(c(-0.1,4.2)) +
  stat_compare_means(comparisons = compar_NK, method = "wilcox.test", size = 6, face="bold") + my_theme 
plot(r) 
dev.off()

# Plotting neoantigen mean expression in low and high infiltrated macrophage groups
pdf(paste0(fig_dir, "/Total_Neoantigen_RNA_expr_Macrophages_Cibersort_score.pdf"), width = 5, height = 5)
compar_Macrophages<- list(c("high", "low"))
s <- ggplot2.violinplot(data=merged_Neo_SNV_Cibersort, xName="Macrophages_score",yName="mean_exp_Total_neoantigens", #fill = "Clonality",
                        groupName="NK_score",
                        groupColors= alpha(c("#2166AC","#B2182B"), .8), showLegend=FALSE,
                        backgroundColor="white", xtitle="Macrophage score", ytitle="Neoantigen Mean expression", 
                        mainTitle="",
                        addDot=T, dotSize=0.5,
                        legendPosition="bottom") + ylim(c(-0.1,4.2)) +
  stat_compare_means(comparisons = compar_Macrophages, method = "wilcox.test", size = 6, face="bold") + my_theme 
plot(s) 
dev.off()


