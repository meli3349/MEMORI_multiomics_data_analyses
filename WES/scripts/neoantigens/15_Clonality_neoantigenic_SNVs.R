
########################################################
## Plotting neoantigenic SNV burden during treatment ##
#######################################################

# This script plots...
# 1a) Number of neoantigenic SNVs during treatment
# 1b) Number of neoantigenic SNVs during treatment for REP and NRP separately
# 1c) Number of neoantigenic SNVs broken down by treatment (for A, B and C broken down by treatment)
# 2)  Number of neoantigenic SNVs in different Becker remission grades
# 3a) Number of neoantigenic clonal and subclonal SNVs over time 
# 3b) Number of neoantigenic clonal and subclonal SNVs for REP and NRP separately
# 3c) Number of neoantigenic clonal and subclonal SNVs broken down by treatment (for A, B and C broken down by treatment)

library(ggplot2)
library(ggpubr)
library(easyGgplot2)
library(reshape2)
library(rstatix)
library("scales")

# prepare output dirs:
fig_dir = "~/analysis/WES/neoantigens/plots/NeoSNV_burden"
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)


## My theme 
my_theme <- theme(
  plot.title = element_text(color="black", size=16, face="bold"),
  axis.title.x = element_text(color="black", size=16, face="bold"),
  axis.text.x = element_text(color="black", size=16, face="bold"),
  axis.title.y = element_text(color="black", size=16, face="bold"),
  axis.text.y = element_text(color="black", size=16, face="bold"),
  legend.title = element_text(color = "black", size = 14, face="bold"),
  legend.text = element_text(color = "black", size = 14, face="bold"),
  
)


# Select samples included for genetic analyses
selected_samples <- read.table("~/analysis/General_files/MEMORI_specific/WES_samples_for_analysis.txt", sep = "\t", header = T)
selected_samples_Seq <- selected_samples[which(selected_samples$PASS_Sequenza_best_fitting_ploidy == "y"), ]
samples <- selected_samples_Seq$SampleID

## Create Mastertable of NeoSNV summaries and add clinical annotations
Neo_SNV_MasterTable <- do.call(rbind,lapply(samples,function(sam)read.table(paste0("~/analysis/WES/neoantigens/filtering_steps/tables/summary_tables_neoSNVs/", sam, ".weakly_filtered_summary_neoSNVs"), header = TRUE,stringsAsFactors = F)))
Neo_SNV_MasterTable$patient <- sapply(Neo_SNV_MasterTable$SampleID, function(x) substr(x,8,13))
Neo_SNV_MasterTable$timepoint <- sapply(Neo_SNV_MasterTable$SampleID, function(x) substr(x,14,14)) 
Neo_SNV_MasterTable$Responsiveness <- sapply(Neo_SNV_MasterTable$SampleID, function(x) substr(x,8,9)) 
Neo_SNV_MasterTable$Treatment_broken_down <- ifelse(Neo_SNV_MasterTable$timepoint == "A", "A", 
                                                     ifelse(Neo_SNV_MasterTable$timepoint == "B", "B",
                                                            ifelse ((Neo_SNV_MasterTable$timepoint == "C" & Neo_SNV_MasterTable$Responsiveness == "RE"), "C_RE",
                                                                    ifelse ((Neo_SNV_MasterTable$timepoint == "C" & Neo_SNV_MasterTable$Responsiveness == "NE"), "C_NE", "XX"))))


# Calculation of NeoSNV per Mb 
# Covered region: 35718732 bases (result from from Qualimap Report)
Neo_SNV_MasterTable$Total_NeoSNVs_per_Mb <-  Neo_SNV_MasterTable$Total_NeoSNVs/35718732*10^6


## 1a) Neoantigens during treatment
my_comparisons <- list(c("A", "B"), c("B", "C"), c("A", "C"))
pdf(paste0(fig_dir, "/Total_neoSNV_burden_all_samples_over_time.pdf"), width = 5, height = 5)
p <- ggplot2.violinplot(data=Neo_SNV_MasterTable, xName="timepoint", yName="Total_NeoSNVs_per_Mb",
                        groupName="timepoint",
                        groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue')), .8), showLegend=FALSE,
                        backgroundColor="white", xtitle="", ytitle="Neoantigenic mutations [SNVs/Mb]", 
                        mainTitle="Neoantigen burden during treatment",
                        addDot=T, dotSize=0.5,
                        removePanelGrid=F,removePanelBorder=F,
                        axisLine=c(0.5, "solid", "black"))+ 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size = 4.5, face="bold")  +
  scale_y_continuous(trans='log10', limits = c(1,400)) # + ylim(c(0, 260))
q <- p + my_theme
plot(q)
dev.off()


## 1b) Neoantigens for REP and NRP separately
pdf(paste0(fig_dir, "/Total_neoSNV_burden_RE_NE_over_time.pdf"), width = 5, height = 5)
p <- ggplot2.violinplot(data=Neo_SNV_MasterTable, xName="timepoint",yName="Total_NeoSNVs_per_Mb", #fill = "Clonality",
                        groupName="timepoint",
                        groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue')), .8), showLegend=FALSE,
                        backgroundColor="white", xtitle="", ytitle="Neoantigenic mutations [SNVs/Mb]", 
                        mainTitle="Neoantigenic mutation burden in NRP and REP",
                        addDot=T, dotSize=0.5,
                        legendPosition="bottom") + 
  facet_grid(~Responsiveness, scales = "free_x", space = "free_x", labeller = as_labeller(c(`NE` = "NRP",`RE` = "REP"))) +
  theme(panel.spacing = unit(0.5, "cm"),
        strip.text = element_text(size = 12, face = "bold")) + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size = 4.5, face="bold") + 
  scale_y_continuous(trans='log10', limits = c(1,400))
q <- p + my_theme + theme(plot.title = element_text(color="black", size=14, face="bold")) + theme(axis.title.x = element_text(color="black", size=14, face="bold"))+ theme(axis.title.y = element_text(color="black", size=14, face="bold")) 
plot(q)
dev.off()

## 1c) for A, B and C broken down by treatment
pdf(paste0(fig_dir, "/Total_neoSNV_burden_broken_down_by_treatment.pdf"), width = 5, height = 5)
p <- ggplot2.violinplot(data=Neo_SNV_MasterTable, xName="Treatment_broken_down",yName="Total_NeoSNVs_per_Mb",
                        groupName="Treatment_broken_down",
                        groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue'), 'khaki1'), .8), showLegend=FALSE,
                        backgroundColor="white", xtitle="", ytitle="Neoantigenic mutations [SNVs/Mb]", 
                        mainTitle="Neoantigenic mutations by treatment",
                        addDot=T, dotSize=0.5,
                        removePanelGrid=F,removePanelBorder=F,
                        axisLine=c(0.5, "solid", "black")) + 
  scale_y_continuous(trans='log10', limits = c(1,400)) +
  stat_compare_means(comparisons=list(c('B','C_NE'),c('B','C_RE'), c('A','C_NE'),c('A','C_RE')), method = "wilcox.test",  p.adjust.method= "bonferroni", size = 4.5, face="bold", 
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.1), symbols = c("****", "***", "**", "*", "ns"))) +
  scale_x_discrete(labels=c('A', 'B', 'C \n (NRP)', 'C \n (REP)'))
q <- p + my_theme 
plot(q)
dev.off()



# 2.) Neoantigens in different Becker remission grades
# Remission grade
regression_df<- read.table("~/analysis/General_files/MEMORI_specific/WES_samples_include_for_analysis.txt", sep = "\t", header = T)
regression_df <- regression_df[,c("SampleID", "Becker_remission_grade")]
Neo_SNV_MasterTable_remission <- merge(Neo_SNV_MasterTable, regression_df, by.x = "SampleID", by.y ="SampleID", all.x = T)


#.a) Neoantigens per remission grade 
pdf(paste0(fig_dir, "/Total_neoantigen_burden_Remission_grades_over_time.pdf"), width = 5, height = 4)
p <- ggplot2.violinplot(data=Neo_SNV_MasterTable_remission, xName="timepoint",yName="Total_NeoSNVs_per_Mb", #fill = "Clonality",
                        groupName="timepoint",
                        groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue')), .8), showLegend=FALSE,
                        backgroundColor="white", xtitle="", ytitle="Neoantigenic mutations [SNVs/Mb]", 
                        mainTitle="Neoantigen burden in remission grades",
                        addDot=T, dotSize=0.5,
                        legendPosition="bottom") + 
  facet_grid(~Becker_remission_grade, scales = "free_x", space = "free_x") +
  theme(panel.spacing = unit(0.5, "cm"),
        strip.text = element_text(size = 12, face = "bold")) +
  scale_y_continuous(trans='log10', limits = c(1,400)) + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size = 4.5, face="bold") #+ ylim(c(0, 8.4))
q <- p + my_theme + theme(plot.title = element_text(color="black", size=14, face="bold")) + theme(axis.title.x = element_text(color="black", size=14, face="bold"))+ theme(axis.title.y = element_text(color="black", size=14, face="bold")) 
plot(q)
dev.off() 



# 3.) Plot clonal and non-clonal neoSNVs

list_neo_mob <- list()
for (sam in samples){
  
  # Read in output of MOBSTER analyses, where SNVs are classified into tail SNVs (=subclonal SNVs) and nontail (in dataframe C1-Cx) (= clonal)
  df <- read.table(paste0("~/analysis/WES/mutations/clonality_analyses/MOBSTER/output_tables/", sam, ".Mobster.CCF.txt"), header = T, stringsAsFactors = F, sep="\t")
  df$SampleID <- sam
  new_df <-  df [,c("SampleID", "chr", "from", "to","ref", "alt", "Mobster_adj_CCF_VAF", "cluster")]
  new_df$SNV <- paste0(new_df$chr, "_", new_df$from, "_", new_df$ref, "_", new_df$alt)
  new_df <- new_df[,c("SNV", "Mobster_adj_CCF_VAF", "cluster")] 

  # Import Neoantigen file
  neoantigen_file <-  read.table(paste0("~/analysis/WES/neoantigens/filtering_steps/tables/neoSNV_tables/", sam, ".weakly_filtered_neoSNVs"), stringsAsFactors = F, header = T)
  neoantigen_file <- neoantigen_file[,c("SNV", "BindLevel", "Novelty")]
  neoantigen_file <- neoantigen_file[!duplicated(neoantigen_file$SNV), ]
  
  # Merge Mobster and NeoSNV files
  merged_neoSNVs_Mobster <- merge(neoantigen_file, new_df, by="SNV")

  ## Calculate tail (=subclonal) and non-Tail (= clonal) SNVs per sample
  NeoSNV_non_tail <- merged_neoSNVs_Mobster[which(merged_neoSNVs_Mobster$cluster != "Tail"), ]
  NeoSNV_total_non_tail <- as.numeric(nrow(NeoSNV_non_tail))
  
  NeoSNV_tail <-  merged_neoSNVs_Mobster[which(merged_neoSNVs_Mobster$cluster == "Tail"), ]
  NeoSNV_total_tail <- as.numeric(nrow(NeoSNV_tail))
  
  headers <- c("SampleID", "NeoSNV_total_non_tail", "NeoSNV_total_tail")
  first_row <- c(sam, NeoSNV_total_non_tail, NeoSNV_total_tail)
  
  NeoSNV_df <- as.data.frame(rbind(headers, first_row))
  names(NeoSNV_df) <- as.matrix(NeoSNV_df[1, ])
  NeoSNV_df <- NeoSNV_df[-1, ]
  rownames(NeoSNV_df) <- NULL
  list_neo_mob[[sam]] <- NeoSNV_df
  }


# Create Mastertable with sum of tail and non-tail neoSNVs
Sum_NeoSNV_MasterTable <- rbindlist(list_neo_mob)

# Calculate sum of tail and non-tail neoSNVs per Mb (covered region 35718732 bases is taken from Qualimap Rerport)
Sum_NeoSNV_MasterTable$Total_NeoSNV_burden <- Sum_NeoSNV_MasterTable$NeoSNV_total_non_tail + Sum_NeoSNV_MasterTable$NeoSNV_total_tail  
Sum_NeoSNV_MasterTable$NeoSNV_total_non_tail_per_Mb <-  Sum_NeoSNV_MasterTable$NeoSNV_total_non_tail/35718732*10^6
Sum_NeoSNV_MasterTable$NeoSNV_total_tail_per_Mb <- Sum_NeoSNV_MasterTable$NeoSNV_total_tail/35718732*10^6
Sum_NeoSNV_MasterTable$Total_NeoSNV_burden_per_Mb <- Sum_NeoSNV_MasterTable$Total_NeoSNV_burden/35718732*10^6

# Select columns of interest
Sum_NeoSNV_MasterTable2 <- Sum_NeoSNV_MasterTable[,c("SampleID", "NeoSNV_total_non_tail_per_Mb", "NeoSNV_total_tail_per_Mb")]
reshaped_Sum_NeoSNV_MasterTable2 <- melt(Sum_NeoSNV_MasterTable2, id="SampleID")

# Add clinical information
reshaped_Sum_NeoSNV_MasterTable2$Timepoint <- sapply(reshaped_Sum_NeoSNV_MasterTable2$SampleID, function(x) substr(x,14, 14))
reshaped_Sum_NeoSNV_MasterTable2$Responsiveness <- sapply(reshaped_Sum_NeoSNV_MasterTable2$SampleID, function(x) substr(x,8, 9))
reshaped_Sum_NeoSNV_MasterTable2$Treatment_broken_down <- ifelse(reshaped_Sum_NeoSNV_MasterTable2$Timepoint == "A", "A", 
                                                    ifelse(reshaped_Sum_NeoSNV_MasterTable2$Timepoint == "B", "B",
                                                           ifelse ((reshaped_Sum_NeoSNV_MasterTable2$Timepoint == "C" & reshaped_Sum_NeoSNV_MasterTable2$Responsiveness == "RE"), "C_RE",
                                                                   ifelse ((reshaped_Sum_NeoSNV_MasterTable2$Timepoint == "C" & reshaped_Sum_NeoSNV_MasterTable2$Responsiveness == "NE"), "C_NE", "XX"))))

# Rename "non-tail NeoSNVs" as "clonal NeoSNVs" and "tail NeoSNVs" as "subclonal NeoSNVs"
reshaped_Sum_NeoSNV_MasterTable2$Clonality <- ifelse(reshaped_Sum_NeoSNV_MasterTable2$variable == "NeoSNV_total_non_tail_per_Mb", "Clonal NeoSNVs", "Subclonal NeoSNVs")

## 3a.) Plot clonal and subclonal NeoSNVs during treatment 
pdf(paste0(fig_dir, "/NeoSNV_clonality_over_time.pdf"), width = 5, height = 5)
p <- ggplot2.violinplot(data=reshaped_Sum_NeoSNV_MasterTable2, xName="Timepoint",yName="value",fill = "Clonality",
                        groupName="Timepoint",
                        groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue')), .8), showLegend=FALSE,
                        backgroundColor="white", xtitle="", ytitle="NeoSNV burden [NeoSNVs/Mb]", 
                        mainTitle="Clonal and subclonal NeoSNV burden",
                        addDot=T, dotSize=0.5,
                        legendPosition="bottom") + 
  facet_grid(~Clonality, scales = "free_x", space = "free_x") +
  theme(panel.spacing = unit(0.5, "cm"),
        strip.text = element_text(size = 12, face = "bold")) + 
  scale_y_continuous(trans='log2', limits = c(0.1,3000)) + 
  stat_compare_means(comparisons=list(c('A','B'),c('A','C'), c('B','C')), method = "wilcox.test")
q <- p + my_theme + theme(plot.title = element_text(color="black", size=14, face="bold")) + theme(axis.title.x = element_text(color="black", size=16, face="bold"))+ theme(axis.title.y = element_text(color="black", size=14, face="bold")) 
plot(q)
dev.off()


## 3b.) Plot subclonal SNVs for REs and NRs separated
reshaped_Sum_NeoSNV_MasterTable2_only_subclonal <- reshaped_Sum_NeoSNV_MasterTable2[which(reshaped_Sum_NeoSNV_MasterTable2$Clonality == "Subclonal NeoSNVs"), ]
  
pdf(paste0(fig_dir, "/NeoSNV_subclonal_RE_NR_over_time.pdf"), width = 5, height = 5)
p <- ggplot2.violinplot(data=reshaped_Sum_NeoSNV_MasterTable2_only_subclonal, xName="Timepoint",yName="value",fill = "Clonality",
                        groupName="Timepoint",
                        groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue')), .8), showLegend=FALSE,
                        backgroundColor="white", xtitle="", ytitle="NeoSNV burden [NeoSNVs/Mb]", 
                        mainTitle="Subclonal NeoSNV burden in REPs and NRPs",
                        addDot=T, dotSize=0.5,
                        legendPosition="bottom") + 
  facet_grid(~Responsiveness, scales = "free_x", space = "free_x", labeller = as_labeller(c(`NE` = "NRP",`RE` = "REP"))) +
  theme(panel.spacing = unit(0.5, "cm"),
        strip.text = element_text(size = 12, face = "bold")) + 
 # scale_y_continuous(trans='log2', limits = c(0.01,3000)) + 
  stat_compare_means(comparisons=list(c('A','B'),c('A','C'), c('B','C')), method = "wilcox.test")
q <- p + my_theme + theme(plot.title = element_text(color="black", size=14, face="bold")) + theme(axis.title.x = element_text(color="black", size=16, face="bold"))+ theme(axis.title.y = element_text(color="black", size=14, face="bold")) 
plot(q)
dev.off()



## 3c) Plot clonal and subclonal SNVs for A, B and C broken down by treatment
pdf(paste0(fig_dir, "/NeoSNV_clonality_broken_down_by_treatment.pdf"), width = 6, height = 5)
p <- ggplot2.violinplot(data=reshaped_Sum_NeoSNV_MasterTable2, xName="Treatment_broken_down",yName="value", fill = "Clonality",
                        groupName="Treatment_broken_down",
                        groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue'), 'khaki1'), .8), showLegend=FALSE,
                        backgroundColor="white", xtitle="", ytitle="NeoSNV burden [SNVs/Mb]", 
                        mainTitle="Clonal and subclonal NeoSNV burden by treatment",
                        addDot=T, dotSize=0.5,
                        removePanelGrid=F,removePanelBorder=F,
                        axisLine=c(0.5, "solid", "black")) + 
  facet_grid(~Clonality, scales = "free_x", space = "free_x") +
  theme(panel.spacing = unit(0.5, "cm"),
        strip.text = element_text(size = 12, face = "bold")) + 
  scale_y_continuous(trans='log2', limits = c(0.1,3000)) +
  stat_compare_means(comparisons=list(c('B','C_NE'),c('B','C_RE'), c('A','C_NE'),c('A','C_RE')), method = "wilcox.test",  p.adjust.method= "bonferroni", size = 4.5, face="bold", 
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.1), symbols = c("****", "***", "**", "*", "ns"))) +
  scale_x_discrete(labels=c('A', 'B', 'C \n (NRP)', 'C \n (REP)'))
q <- p + my_theme 
plot(q)
dev.off()
