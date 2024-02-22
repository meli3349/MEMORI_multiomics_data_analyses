
##########################################################################
## Plotting number of Intogen drivers with genetic alterations
##########################################################################
# This script calls copy number states, presence of SNVs or indels for Intogen driver for each sample
# and plots... 
# 1. the number of Intogen drivers with genetic alterations during treatment in REP and NRP
# 2. the number of Intogen drivers with CNA during treatment in REP and NRP 
# 3. the number of Intogen drivers with SNVs or indels during treatment in REP and NRP

library(data.table)
library(dplyr)

# prepare output dirs:
fig_dir = "~/analysis/WES/Intogen_driver/plots"
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

# Read in Intogene list with driver genes form esophageal Cancer and gastric cancer 
# https://www.intogen.org/ 
Intogen_EC_GC_list  <- read.table("~/analysis/General_files/external_files/Intogen_drivers.txt", stringsAsFactors = F, header = T)

# ENSEMBLE/ HUGO Gene ID and gene location list
ENSEMBLE_list<- read.table("~/analysis/General_files/external_files/compiledGeneInfo.txt", stringsAsFactors = F, header = T)
chr_filt <- paste0("chr", c(1:23, "X", "Y"))
ENSEMBLE_list <- ENSEMBLE_list[which(ENSEMBLE_list$Chr %in% chr_filt),]

# Get ENSEMBL Gene Names and gene location for drivers in Intogen list
Intogen_loc <- merge(Intogen_list, ENSEMBLE_list, by.x = "Symbol", by.y = "Name")
Intogen_loc  <- Intogen_loc[ ,c("Symbol", "GeneID", "Chr", "Strand", "Start", "End")]

# Select samples included for genetic analyses
selected_samples <- read.table("~/analysis/General_files/MEMORI_specific/WES_samples_include_for_analysis.txt", sep = "\t", header = T)
selected_samples_Seq <- selected_samples[which(selected_samples$PASS_Sequenza_best_fitting_ploidy == "y"), ]
samples <- selected_samples_Seq$SampleID


# Table per sample with Copy number state for each driver 
List_driver_CNS <- list()
for (sam in samples) { 
  df <- read.table(paste0("~/analysis/WES/input_files/Sequenza/best_fitting_ploidies_per_sample/", sam, "_segments.txt"), header = T, stringsAsFactors = F)
  colnames(df)[colnames(df) == "chromosome"] <- "Chr"
  colnames(df)[colnames(df) == "start.pos"] <- "Start"
  colnames(df)[colnames(df) == "end.pos"] <- "End"
  df$copy_state <- paste0(df$A, "_", df$B)
  df <- df[,c("Chr", "Start", "End", "copy_state", "CNt")]
  
  Intogen_loc  <- as.data.table(Intogen_loc)
  df <- as.data.table(df)
  
  setkey(df, Chr, Start, End)
  merged_df <- foverlaps(Intogen_loc, df, type="within")
  
  export_df_1 <- merged_df[,c("copy_state")]
  List_driver_CNS[[sam]] <- export_df_1 
  genes <- merged_df[,c("Symbol")]
  
  
}

# Create Mastertable with copy number state for each Intogen driver
df <- list.cbind(List_driver_CNS)
MasterTable_detailed <- as.data.frame(cbind(genes,df),stringsAsFactors=F)
colnames(MasterTable_detailed)  <- c("GeneID", samples)


melt_MasterTable_detailed<- melt(MasterTable_detailed, id="GeneID")
melt_MasterTable_detailed$colID <- melt_MasterTable_detailed$value


# Combining copy number states into copy number categories
loss <- "1_0"
high_LOH <- c("20_0", "19_0", "18_0", "17_0", "16_0", "15_0", "14_0", "13_0", "12_0", "11_0", "10_0", "9_0", "8_0", "7_0", "6_0", "5_0", "4_0", "3_0")
normal_LOH <- "2_0"
triploid_ampl <- "2_1"
normal_diploid <- "1_1"

melt_MasterTable_detailed$colID[melt_MasterTable_detailed$colID %in% loss] <- "Loss 1:0"
melt_MasterTable_detailed$colID[melt_MasterTable_detailed$colID %in% normal_LOH] <- "LOH 2:0"
melt_MasterTable_detailed$colID[melt_MasterTable_detailed$colID %in% high_LOH] <- "LOH >2:0"
melt_MasterTable_detailed$colID[melt_MasterTable_detailed$colID %in% triploid_ampl] <- "Amplification = 3"
melt_MasterTable_detailed$colID[melt_MasterTable_detailed$colID %in% normal_diploid] <- "diploid"
melt_MasterTable_detailed$colID[which(grepl("_", melt_MasterTable_detailed$colID, fixed = TRUE))] <- "Amplification > 3"
melt_MasterTable_detailed$colID[is.na(melt_MasterTable_detailed$colID)] <- "unknown"


## Add indels in driver genes to Mastertable 
# MasterTable_indels shows cancer cell fraction (CCF) harbouring indels in driver genes. NA indicates no indel in cancer driver.
# Only indels with CCF>0.1 and <1.4 are included, as these are defined as biologically plausible and relevant
MasterTable_indels <- read.table("~/analysis/WES/input_files/mutation_results/Intogen_driver/Mastertable_CCF_indels_Intogen_drivers.txt", stringsAsFactors = F, header = T)
melt_MasterTable_indels <- melt(MasterTable_indels, id="GeneID")
melt_MasterTable_indels$colID_indels <- ifelse(melt_MasterTable_indels$value> 0.1 & melt_MasterTable_indels$value <1.4, "Indel", NA)
melt_MasterTable_detailed <- merge(melt_MasterTable_detailed, melt_MasterTable_indels, by = c("GeneID", "variable"), all.x = T)


## Add SNVs in driver genes to Mastertable 
# MasterTable_SNVs shows cancer cell fraction (CCF) harbouring SNVs in driver genes. NA indicates no indel in cancer driver.
# Only SNVs with CCF>0.1 and <1.4 are included, as these are defined as biologically plausible and relevant
MasterTable_SNVs <- read.table("~/analysis/WES/input_files/mutation_results/Intogen_driver/Mastertable_CCF_SNVs_Intogen_drivers.txt", stringsAsFactors = F, header = T)
melt_MasterTable_SNVs <- melt(MasterTable_SNVs, id="GeneID")
melt_MasterTable_SNVs$colID_SNV <- ifelse((melt_MasterTable_SNVs$value> 0.1 & melt_MasterTable_SNVs$value <1.4), "SNV", NA)
melt_MasterTable_detailed <- merge(melt_MasterTable_detailed, melt_MasterTable_SNVs, by = c("GeneID", "variable"), all.x = T)


# Define 1.) altered driver genes (if CNS non-dipoloid, SNV or Indel in driver gene)
# Define 2.) CNA in driver gene (if CNS non-dipoloid )
# Define 3.) SNV in driver genes (if  SNV in driver gene)
# Define 4.) Indel in driver genes (if Indel in driver gene)
melt_MasterTable_detailed$Altered_driver_gene <- if_else((melt_MasterTable_detailed$colID == "diploid" & is.na(melt_MasterTable_detailed$colID_indels) &  is.na(melt_MasterTable_detailed$colID_SNV)), "No", "Yes")
melt_MasterTable_detailed$CNA_driver_gene <- if_else((melt_MasterTable_detailed$colID == "diploid"), "No", "Yes")
melt_MasterTable_detailed$SNV_or_Indel_driver_gene <- if_else(is.na(melt_MasterTable_detailed$colID_SNV) & is.na(melt_MasterTable_detailed$colID_indels), "No", "Yes")


# 1) Calculate number of altered driver genes for each sample
melt_MasterTable_detailed_sum_alterations <- melt_MasterTable_detailed[, c("GeneID", "variable", "Altered_driver_gene")]
MasterTable_detailed_driver_genes_sum_alterations <- dcast(melt_MasterTable_detailed_sum_alterations, variable ~ GeneID)
MasterTable_detailed_driver_genes_sum_alterations$Sum_altered_genes <- rowSums(MasterTable_detailed_driver_genes_sum_alterations[,c(2:ncol(MasterTable_detailed_driver_genes_sum_alterations)), ] == "Yes")
colnames(MasterTable_detailed_driver_genes_sum_alterations)[colnames(MasterTable_detailed_driver_genes_sum_alterations) == "variable"] <- "SampleID"

# 2.) Calculate number of driver genes with CNA for each sample
melt_MasterTable_detailed_only_CNA<- melt_MasterTable_detailed[, c("GeneID", "variable", "CNA_driver_gene")]
MasterTable_detailed_driver_genes_only_CNA <- dcast(melt_MasterTable_detailed_only_CNA, variable ~ GeneID)
MasterTable_detailed_driver_genes_only_CNA$Sum_driver_genes_with_CNA <- rowSums(MasterTable_detailed_driver_genes_only_CNA[,c(2:ncol(MasterTable_detailed_driver_genes_only_CNA)), ] == "Yes")
colnames(MasterTable_detailed_driver_genes_only_CNA)[colnames(MasterTable_detailed_driver_genes_only_CNA) == "variable"] <- "SampleID"

# 3.) Calculate number of driver genes with SNV or Indel for each sample
melt_MasterTable_detailed_only_Indel_SNV<- melt_MasterTable_detailed[, c("GeneID", "variable", "SNV_or_Indel_driver_gene")]
MasterTable_detailed_driver_genes_only_Indel_SNV <- dcast(melt_MasterTable_detailed_only_Indel_SNV, variable ~ GeneID)
MasterTable_detailed_driver_genes_only_Indel_SNV$Sum_driver_genes_with_SNV_or_Indel <- rowSums(MasterTable_detailed_driver_genes_only_Indel_SNV[,c(2:ncol(MasterTable_detailed_driver_genes_only_Indel_SNV)), ] == "Yes")
colnames(MasterTable_detailed_driver_genes_only_Indel_SNV)[colnames(MasterTable_detailed_driver_genes_only_Indel_SNV) == "variable"] <- "SampleID"



####### 1.)  Plot number of driver alterations 
## 1a) for REP and NRP together
MasterTable_detailed_driver_genes_sum_alterations2 <- MasterTable_detailed_driver_genes_sum_alterations
MasterTable_detailed_driver_genes_sum_alterations2$Timepoint <- sapply(MasterTable_detailed_driver_genes_sum_alterations2$SampleID, function(x) substr(x,14, 14))
MasterTable_detailed_driver_genes_sum_alterations2$Responsiveness <- sapply(MasterTable_detailed_driver_genes_sum_alterations2$SampleID, function(x) substr(x,8, 9))

pdf(paste0(fig_dir, "/Total_driver_altertions_over_time.pdf"), width = 5, height = 5)
p <- ggplot2.violinplot(data=MasterTable_detailed_driver_genes_sum_alterations2, xName="Timepoint",yName="Sum_altered_genes",
                        groupName="Timepoint",
                        groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue')), .8), showLegend=FALSE,
                        backgroundColor="white", xtitle="       Timepoint", ytitle="Number of altered driver genes", 
                        mainTitle="Genetic alterations in driver genes",
                        addDot=T, dotSize=0.5,
                        removePanelGrid=F,removePanelBorder=F,
                        axisLine=c(0.5, "solid", "black")) + 
  stat_compare_means(comparisons=list(c('A','B'),c('A','C'), c('B','C')), p.adjust.methods = "bonferroni", method = "wilcox.test", size = 4.5, face="bold")
q <- p + my_theme
plot(q)
dev.off()


## 1b) for REP and NRP separately
pdf(paste0(fig_dir, "/Total_driver_altertions_RE_NE_over_time.pdf"), width = 4.5, height = 4.5)
p <- ggplot2.violinplot(data=MasterTable_detailed_driver_genes_sum_alterations2, xName="Timepoint",yName="Sum_altered_genes",
                        groupName="Timepoint",
                        groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue')), .8), showLegend=FALSE,
                        backgroundColor="white", xtitle="       Timepoint", ytitle="Number of altered driver genes", 
                        mainTitle="Genetic alterations in driver genes",
                        addDot=T, dotSize=0.5,
                        legendPosition="bottom") + 
  facet_grid(~Responsiveness , scales = "free_x", space = "free_x",  labeller = as_labeller(c(`NE` = "NRP",`RE` = "REP"))) +
  theme(panel.spacing = unit(0.5, "cm"),
        strip.text = element_text(size = 12, face = "bold"))+
  stat_compare_means(comparisons=list(c('A','B'),c('A','C'), c('B','C')), p.adjust.methods = "bonferroni", method = "wilcox.test")
q <- p + my_theme + theme(plot.title = element_text(color="black", size=14, face="bold")) + theme(axis.title.x = element_text(color="black", size=14, face="bold"))+ theme(axis.title.y = element_text(color="black", size=14, face="bold")) 
plot(q)
dev.off()

# Wilcox test for sum of altered driver genes for NR and RE 
stat.test_RE_sum_alterations <- MasterTable_detailed_driver_genes_sum_alterations2[which(MasterTable_detailed_driver_genes_sum_alterations2$Responsiveness == "RE"), ] %>%
  wilcox_test(Sum_altered_genes ~ Timepoint) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

stat.mean_RE_sum_alterations <- MasterTable_detailed_driver_genes_sum_alterations2[which(MasterTable_detailed_driver_genes_sum_alterations2$Responsiveness == "RE"), ] %>%
  group_by(Timepoint) %>%
  get_summary_stats(Sum_altered_genes, show = c("mean", "sd", "median", "iqr"))


stat.test_NR_sum_alterations <- MasterTable_detailed_driver_genes_sum_alterations2[which(MasterTable_detailed_driver_genes_sum_alterations2$Responsiveness == "NE"), ] %>%
  wilcox_test(Sum_altered_genes ~ Timepoint) %>%
  get_summary_stats(Sum_altered_genes)
  adjust_pvalue(method = "BH") %>%
  add_significance()

  stat.mean_NR_sum_alterations <- MasterTable_detailed_driver_genes_sum_alterations2[which(MasterTable_detailed_driver_genes_sum_alterations2$Responsiveness == "NE"), ] %>%
    group_by(Timepoint) %>%
    get_summary_stats(Sum_altered_genes, show = c("mean", "sd", "median", "iqr"))

####### 2.)  Plot number of driver genes with CNA
## 2a) for REP and NRP together
MasterTable_detailed_driver_genes_only_CNA2<- MasterTable_detailed_driver_genes_only_CNA
MasterTable_detailed_driver_genes_only_CNA2$Timepoint <- sapply(MasterTable_detailed_driver_genes_only_CNA2$SampleID, function(x) substr(x,14, 14))
MasterTable_detailed_driver_genes_only_CNA2$Responsiveness <- sapply(MasterTable_detailed_driver_genes_only_CNA2$SampleID, function(x) substr(x,8, 9))

pdf(paste0(fig_dir, "/Total_driver_CNA_over_time.pdf"), height = 5, width = 5)
p <- ggplot2.violinplot(data=MasterTable_detailed_driver_genes_only_CNA2, xName="Timepoint",yName="Sum_driver_genes_with_CNA",
                        groupName="Timepoint",
                        groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue')), .8), showLegend=FALSE,
                        backgroundColor="white", xtitle="       Timepoint", ytitle="Number of driver genes with CNAs", 
                        mainTitle="Driver genes with CNAs",
                        addDot=T, dotSize=0.5,
                        removePanelGrid=F,removePanelBorder=F,
                        axisLine=c(0.5, "solid", "black")) + 
  stat_compare_means(comparisons=list(c('A','B'),c('A','C'), c('B','C')), method = "wilcox.test", size = 4.5, face="bold")
q <- p + my_theme
plot(q)
dev.off()

## 2b) for REP and NRP separately
pdf(paste0(fig_dir, "/Total_driver_CNA_RE_NE_over_time.pdf"), width = 4.5, height = 4.5)
p <- ggplot2.violinplot(data=MasterTable_detailed_driver_genes_only_CNA2, xName="Timepoint",yName="Sum_driver_genes_with_CNA",
                        groupName="Timepoint",
                        groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue')), .8), showLegend=FALSE,
                        backgroundColor="white", xtitle="       Timepoint", ytitle="Number of driver genes with CNAs", 
                        mainTitle="Driver genes with CNAs",
                        addDot=T, dotSize=0.5,
                        legendPosition="bottom") + 
  facet_grid(~Responsiveness , scales = "free_x", space = "free_x",  labeller = as_labeller(c(`NE` = "NRP",`RE` = "REP"))) +
  theme(panel.spacing = unit(0.5, "cm"),
        strip.text = element_text(size = 12, face = "bold"))+
  stat_compare_means(comparisons=list(c('A','B'),c('A','C'), c('B','C')), method = "wilcox.test")
q <- p + my_theme + theme(plot.title = element_text(color="black", size=14, face="bold")) + theme(axis.title.x = element_text(color="black", size=14, face="bold"))+ theme(axis.title.y = element_text(color="black", size=14, face="bold")) 
plot(q)
dev.off()


# Wilcox test for sum of altered driver genes for NR and RE 
stat.test_RE_sum_CNA <- MasterTable_detailed_driver_genes_only_CNA2[which(MasterTable_detailed_driver_genes_only_CNA2$Responsiveness == "RE"), ] %>%
  wilcox_test(Sum_driver_genes_with_CNA ~ Timepoint) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

stat.mean_RE_sum_CNA <- MasterTable_detailed_driver_genes_only_CNA2[which(MasterTable_detailed_driver_genes_only_CNA2$Responsiveness == "RE"), ] %>%
  group_by(Timepoint) %>%
  get_summary_stats(Sum_driver_genes_with_CNA, show = c("mean", "sd", "median", "iqr"))


stat.test_NR_sum_CNA <- MasterTable_detailed_driver_genes_only_CNA2[which(MasterTable_detailed_driver_genes_only_CNA2$Responsiveness == "NE"), ] %>%
  wilcox_test(Sum_driver_genes_with_CNA ~ Timepoint) %>%
adjust_pvalue(method = "BH") %>%
  add_significance()

stat.mean_NR_sum_CNA <- MasterTable_detailed_driver_genes_only_CNA2[which(MasterTable_detailed_driver_genes_only_CNA2$Responsiveness == "NE"), ] %>%
  group_by(Timepoint) %>%
  get_summary_stats(Sum_driver_genes_with_CNA, show = c("mean", "sd", "median", "iqr"))


####### 3.)  Plot number of driver genes with SNVs or indels
## 3a) for REP and NRP together
MasterTable_detailed_driver_genes_only_Indel_SNV2 <- MasterTable_detailed_driver_genes_only_Indel_SNV
MasterTable_detailed_driver_genes_only_Indel_SNV2$Timepoint <- sapply(MasterTable_detailed_driver_genes_only_Indel_SNV2$SampleID, function(x) substr(x,14, 14))
MasterTable_detailed_driver_genes_only_Indel_SNV2$Responsiveness <- sapply(MasterTable_detailed_driver_genes_only_Indel_SNV2$SampleID, function(x) substr(x,8, 9))

pdf(paste0(fig_dir, "/Total_driver_SNV_Indel_over_time.pdf"), width = 5, height = 5)
p <- ggplot2.violinplot(data=MasterTable_detailed_driver_genes_only_Indel_SNV2, xName="Timepoint",yName="Sum_driver_genes_with_SNV_or_Indel",
                        groupName="Timepoint",
                        groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue')), .8), showLegend=FALSE,
                        backgroundColor="white", xtitle="       Timepoint", ytitle="Number of driver genes with SNVs or indels", 
                        mainTitle="Driver genes with SNVs or indels",
                        addDot=T, dotSize=0.5,
                        removePanelGrid=F,removePanelBorder=F,
                        axisLine=c(0.5, "solid", "black")) + 
  stat_compare_means(comparisons=list(c('A','B'),c('A','C'), c('B','C')), method = "wilcox.test", size = 4.5, face="bold")
q <- p + my_theme
plot(q)
dev.off()

## 3b) for REP and NRP separately
pdf(paste0(fig_dir, "/Total_driver_SNV_Indel_RE_NE_over_time.pdf"), width = 4.5, height = 4.5)
p <- ggplot2.violinplot(data=MasterTable_detailed_driver_genes_only_Indel_SNV2, xName="Timepoint",yName="Sum_driver_genes_with_SNV_or_Indel",
                        groupName="Timepoint",
                        groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue')), .8), showLegend=FALSE,
                        backgroundColor="white", xtitle="       Timepoint", ytitle="Driver genes with SNVs or indels", 
                        mainTitle="Driver genes with SNVs or indels",
                        addDot=T, dotSize=0.5,
                        legendPosition="bottom") + 
  facet_grid(~Responsiveness , scales = "free_x", space = "free_x",  labeller = as_labeller(c(`NE` = "NRP",`RE` = "REP"))) +
  theme(panel.spacing = unit(0.5, "cm"),
        strip.text = element_text(size = 12, face = "bold"))+
  stat_compare_means(comparisons=list(c('A','B'),c('A','C'), c('B','C')), method = "wilcox.test")
q <- p + my_theme + theme(plot.title = element_text(color="black", size=14, face="bold")) + theme(axis.title.x = element_text(color="black", size=14, face="bold"))+ theme(axis.title.y = element_text(color="black", size=14, face="bold")) 
plot(q)
dev.off()

# Wilcox test for sum of altered driver genes for NRP and REP 
stat.test_RE_sum_CNA <- MasterTable_detailed_driver_genes_only_Indel_SNV2[which(MasterTable_detailed_driver_genes_only_Indel_SNV2$Responsiveness == "RE"), ] %>%
  wilcox_test(Sum_driver_genes_with_SNV_or_Indel ~ Timepoint) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

stat.mean_RE_sum_CNA <- MasterTable_detailed_driver_genes_only_Indel_SNV2[which(MasterTable_detailed_driver_genes_only_Indel_SNV2$Responsiveness == "RE"), ] %>%
  group_by(Timepoint) %>%
  get_summary_stats(Sum_driver_genes_with_SNV_or_Indel, show = c("mean", "sd", "median", "iqr"))


stat.test_NR_sum_CNA <- MasterTable_detailed_driver_genes_only_Indel_SNV2[which(MasterTable_detailed_driver_genes_only_Indel_SNV2$Responsiveness == "NE"), ] %>%
  wilcox_test(Sum_driver_genes_with_SNV_or_Indel ~ Timepoint) %>%
  get_summary_stats(Sum_driver_genes_with_SNV_or_Indel)
adjust_pvalue(method = "BH") %>%
  add_significance()

stat.mean_NR_sum_CNA <- MasterTable_detailed_driver_genes_only_Indel_SNV2[which(MasterTable_detailed_driver_genes_only_Indel_SNV2$Responsiveness == "NE"), ] %>%
  group_by(Timepoint) %>%
  get_summary_stats(Sum_driver_genes_with_SNV_or_Indel, show = c("mean", "sd", "median", "iqr"))

