
########################################################
##  SNV burden during treatment 
########################################################

# This script plots...
# 1a) Number of SNVs during treatment
# 1b) Number of SNVs during treatment for REP and NRP separately
# 1c) Number of SNVs broken down by treatment (for A, B and C broken down by treatment)
# 2a) Number of clonal and subclonal SNVs during treatment
# 2b) Number of clonal and subclonal SNVs for REP and NRP separately
# 2c) Number of clonal and subclonal SNVs for REP and NRP separately

library(reshape2)
library(ggplot2)
library(easyGgplot2)
library(rstatix)
library("scales")
library(data.table)

# prepare output dirs:
fig_dir = "~/analysis/WES/mutations/plots/SNV_burden"
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

# Select samples included for genetic analyses
selected_samples <- read.table("~/analysis/General_files/MEMORI_specific/WES_samples_for_analysis.txt", sep = "\t", header = T)
selected_samples_Seq <- selected_samples[which(selected_samples$PASS_Sequenza_best_fitting_ploidy == "y"), ]
samples <- selected_samples_Seq$SampleID

list_clonality <- list()
for (sam in samples){
  
  df <- read.table(paste0("~/analysis/WES/mutations/MOBSTER/output_tables/", sam, ".Mobster.CCF.txt"), header = T, stringsAsFactors = F, sep="\t")
  df$SampleID <- sam
  new_df <-  df [,c("SampleID", "chr", "from", "to","ref", "alt", "Mobster_adj_CCF_VAF", "cluster")]
  # write.table(new_df, file = paste0("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/mutation_results/SNV_clonality_over_time/Tables/Tables_from_Mobster_CCF_best_fit_ploidy/reduced_tables/", sam, ".reduced.Mobster.CCF.txt"), sep = "\t", row.names = F, quote = F)
  
  ## Calculate Tail and Non-Tail SNV per sample
  SNV_non_tail <- new_df[which(new_df$cluster != "Tail"), ]
  SNV_total_non_tail <- as.numeric(nrow(SNV_non_tail))
  
  SNV_tail <- new_df[which(new_df$cluster == "Tail"), ]
  SNV_total_tail <- as.numeric(nrow(SNV_tail))
  
  headers <- c("SampleID", "SNV_total_non_tail", "SNV_total_tail")
  first_row <- c(sam, SNV_total_non_tail, SNV_total_tail)
  
  SNV_df <- as.data.frame(rbind(headers, first_row))
  names(SNV_df) <- as.matrix(SNV_df[1, ])
  SNV_df <- SNV_df[-1, ]
  rownames(SNV_df) <- NULL
  list_clonality[[sam]] <- SNV_df
  
}

## Create Mastertable of SNV clonality summary
Sum_SNV_MasterTable <- rbindlist(list_clonality)

# Calculation of SNVs per Mb 
# Covered region: 35718732 bases (result from from Qualimap Report)
Sum_SNV_MasterTable$Total_Mutation_burden <- Sum_SNV_MasterTable$SNV_total_non_tail + Sum_SNV_MasterTable$SNV_total_tail  
Sum_SNV_MasterTable$SNV_total_non_tail_per_Mb <-  Sum_SNV_MasterTable$SNV_total_non_tail/35718732*10^6
Sum_SNV_MasterTable$SNV_total_tail_per_Mb <- Sum_SNV_MasterTable$SNV_total_tail/35718732*10^6
Sum_SNV_MasterTable$Total_Mutation_burden_per_Mb <- Sum_SNV_MasterTable$Total_Mutation_burden/35718732*10^6

# Add clinical information
Sum_SNV_MasterTable1 <- Sum_SNV_MasterTable
Sum_SNV_MasterTable1$Timepoint <- sapply(Sum_SNV_MasterTable1$SampleID, function(x) substr(x,14, 14))
Sum_SNV_MasterTable1$Responsiveness <- sapply(Sum_SNV_MasterTable1$SampleID, function(x) substr(x,8, 9))
Sum_SNV_MasterTable1$Treatment_broken_down <- ifelse(Sum_SNV_MasterTable1$Timepoint == "A", "A", 
                                                     ifelse(Sum_SNV_MasterTable1$Timepoint == "B", "B",
                                                            ifelse ((Sum_SNV_MasterTable1$Timepoint == "C" & Sum_SNV_MasterTable1$Responsiveness == "RE"), "C_RE",
                                                                    ifelse ((Sum_SNV_MasterTable1$Timepoint == "C" & Sum_SNV_MasterTable1$Responsiveness == "NE"), "C_NE", "XX"))))
########
## 1.)  Plot Total SNV burden during treatment 
## 1a) for overall cohort
pdf(paste0(fig_dir, "/Total_Mutation_burden_from_Mobster_over_time.pdf"), width = 5, height = 5)
p <- ggplot2.violinplot(data=Sum_SNV_MasterTable1, xName="Timepoint",yName="Total_Mutation_burden_per_Mb",
                   groupName="Timepoint",
                   groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue')), .8), showLegend=FALSE,
                   backgroundColor="white", xtitle="", ytitle="Mutation burden [SNVs/Mb]", 
                   mainTitle="Mutation burden during treatment",
                   addDot=T, dotSize=0.5,
                   removePanelGrid=F,removePanelBorder=F,
                   axisLine=c(0.5, "solid", "black")) + 
  scale_y_continuous(trans='log10', limits = c(4,4000)) +
  stat_compare_means(comparisons=list(c('A','B'),c('A','C'), c('B','C')), method = "wilcox.test", size = 4.5, face="bold")
q <- p + my_theme
plot(q)
dev.off()


## 1b) for REPs and NRPs separately
pdf(paste0(fig_dir, "/Total_Mutation_burden_from_Mobster_RE_NE_over_time.pdf"), width = 5, height = 5)
p <- ggplot2.violinplot(data=Sum_SNV_MasterTable1, xName="Timepoint",yName="Total_Mutation_burden_per_Mb", fill = "Clonality",
                        groupName="Timepoint",
                        groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue')), .8), showLegend=FALSE,
                        backgroundColor="white", xtitle="", ytitle="Mutation burden [SNVs/Mb]", 
                        mainTitle="Mutation burden in NR and RE",
                        addDot=T, dotSize=0.5,
                        legendPosition="bottom") + 
  facet_grid(~Responsiveness , scales = "free_x", space = "free_x",  labeller = as_labeller(c(`NE` = "NR",`RE` = "RE"))) +
  theme(panel.spacing = unit(0.5, "cm"),
        strip.text = element_text(size = 12, face = "bold"))+
  stat_compare_means(comparisons=list(c('A','B'),c('A','C'), c('B','C')), method = "wilcox.test") + 
  scale_y_continuous(trans='log10', limits = c(4,4000)) 
q <- p + my_theme + theme(plot.title = element_text(color="black", size=14, face="bold")) + theme(axis.title.x = element_text(color="black", size=14, face="bold"))+ theme(axis.title.y = element_text(color="black", size=14, face="bold")) 
plot(q)
dev.off()


## 1c) for A, B and C broken down by treatment
pdf(paste0(fig_dir, "/Total_Mutation_burden_from_Mobster_broken_down_by_treatment.pdf"), width = 5, height = 5)
p <- ggplot2.violinplot(data=Sum_SNV_MasterTable1, xName="Treatment_broken_down",yName="Total_Mutation_burden_per_Mb",
                        groupName="Treatment_broken_down",
                        groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue'), 'khaki1'), .8), showLegend=FALSE,
                        backgroundColor="white", xtitle="", ytitle="Mutation burden [SNVs/Mb]", 
                        mainTitle="Mutation burden by treatment",
                        addDot=T, dotSize=0.5,
                        removePanelGrid=F,removePanelBorder=F,
                        axisLine=c(0.5, "solid", "black")) + 
  scale_y_continuous(trans='log10', limits = c(1,10000)) +
  stat_compare_means(comparisons=list(c('B','C_NE'),c('B','C_RE'), c('A','C_NE'),c('A','C_RE')), method = "wilcox.test",  p.adjust.method= "bonferroni", size = 4.5, face="bold", 
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.1), symbols = c("****", "***", "**", "*", "ns"))) +
  scale_x_discrete(labels=c('A', 'B', 'C \n (NRP)', 'C \n (REP)'))
q <- p + my_theme
plot(q)
q <- p + my_theme 
plot(q)
dev.off()

### 2.) PLot Violiplots of clonal and subclonal SNVs during treatment 
## 2a.) for overall cohort 
Sum_SNV_MasterTable2 <- Sum_SNV_MasterTable[,c("SampleID", "SNV_total_non_tail_per_Mb", "SNV_total_tail_per_Mb")]
reshaped_Sum_SNV_MasterTable <- melt(Sum_SNV_MasterTable2, id="SampleID")
reshaped_Sum_SNV_MasterTable$Timepoint <- sapply(reshaped_Sum_SNV_MasterTable$SampleID, function(x) substr(x,14, 14))
reshaped_Sum_SNV_MasterTable$Responsiveness <- sapply(reshaped_Sum_SNV_MasterTable$SampleID, function(x) substr(x,8, 9))
reshaped_Sum_SNV_MasterTable$Clonality <- ifelse(reshaped_Sum_SNV_MasterTable$variable == "SNV_total_non_tail_per_Mb", "Clonal SNVs", "Subclonal SNVs")

pdf(paste0(fig_dir, "/Mutation_burden_from_Mobster_with_clonality_over_time.pdf"), width = 6, height = 6)
p <- ggplot2.violinplot(data=reshaped_Sum_SNV_MasterTable, xName="Timepoint",yName="value",fill = "Clonality",
                   groupName="Timepoint",
                   groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue')), .8), showLegend=FALSE,
                   backgroundColor="white", xtitle="", ytitle="Mutation burden [SNVs/Mb]", 
                   mainTitle="Clonal and subclonal mutation burden",
                   addDot=T, dotSize=0.5,
                   legendPosition="bottom") + 
  facet_grid(~Clonality, scales = "free_x", space = "free_x") +
  theme(panel.spacing = unit(0.5, "cm"),
        strip.text = element_text(size = 12, face = "bold")) + 
  scale_y_continuous(trans='log10', limits = c(0.1,9000)) + 
  stat_compare_means(comparisons=list(c('A','B'),c('A','C'), c('B','C')), method = "wilcox.test")
q <- p + my_theme + theme(plot.title = element_text(color="black", size=14, face="bold")) + theme(axis.title.x = element_text(color="black", size=16, face="bold"))+ theme(axis.title.y = element_text(color="black", size=14, face="bold")) 
plot(q)
dev.off()

## 2b.) for REPs and NRPs separately
RE_reshaped_Sum_SNV_MasterTable <- reshaped_Sum_SNV_MasterTable[(reshaped_Sum_SNV_MasterTable$Responsiveness == "RE"), ]
NE_reshaped_Sum_SNV_MasterTable <- reshaped_Sum_SNV_MasterTable[(reshaped_Sum_SNV_MasterTable$Responsiveness == "NE"), ]


# for REPs
pdf(paste0(fig_dir, "/Graphs/Mutation_burden_from_Mobster_with_clonality_over_time_only_RE.pdf"), width = 6, height = 6)
p_RE <- ggplot2.violinplot(data=RE_reshaped_Sum_SNV_MasterTable, xName="Timepoint",yName="value",fill = "Clonality",
                        groupName="Timepoint",
                        groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue')), .8), showLegend=FALSE,
                        backgroundColor="white", xtitle="Timepoint", ytitle="Mutation burden [SNVs/Mb]", 
                        mainTitle="Clonal and subclonal mutation burden in RE",
                        addDot=T, dotSize=0.5,
                        legendPosition="bottom") + 
  facet_grid(~Clonality, scales = "free_x", space = "free_x") +
  theme(panel.spacing = unit(0.5, "cm"),
        strip.text = element_text(size = 12, face = "bold")) +
  stat_compare_means(comparisons=list(c('A','B'),c('A','C'), c('B','C')), method = "wilcox.test")
q <- p_RE + my_theme + theme(plot.title = element_text(color="black", size=14, face="bold")) + theme(axis.title.x = element_text(color="black", size=14, face="bold"))+ theme(axis.title.y = element_text(color="black", size=14, face="bold")) 
plot(q)
dev.off()


# for NRPs
pdf(paste0(fig_dir, "/Mutation_burden_from_Mobster_with_clonality_over_time_only_NR.pdf"), width = 6, height = 6)
p_NE <- ggplot2.violinplot(data=NE_reshaped_Sum_SNV_MasterTable, xName="Timepoint",yName="value",fill = "Clonality",
                           groupName="Timepoint",
                           groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue')), .8), showLegend=FALSE,
                           backgroundColor="white", xtitle="Timepoint", ytitle="Mutation burden [SNVs/Mb]", 
                           mainTitle="Clonal and subclonal mutation burden in NR",
                           addDot=T, dotSize=0.5,
                           legendPosition="bottom")  + 
  facet_grid(~Clonality, scales = "free_x", space = "free_x") +
  theme(panel.spacing = unit(0.5, "cm"),
        strip.text = element_text(size = 12, face = "bold"))+
  stat_compare_means(comparisons=list(c('A','B'),c('A','C'), c('B','C')), method = "wilcox.test")
q <- p_NE + my_theme + theme(plot.title = element_text(color="black", size=14, face="bold")) + theme(axis.title.x = element_text(color="black", size=14, face="bold"))+ theme(axis.title.y = element_text(color="black", size=14, face="bold")) 
plot(q)
dev.off()

## 2c.) for A, B and C broken down by treatment
pdf(paste0(fig_dir, "/Mutation_burden_from_Mobster_with_clonality_broken_down_by_treatment.pdf"), width = 7, height = 6)
p <- ggplot2.violinplot(data=reshaped_Sum_SNV_MasterTable, xName="Treatment_broken_down",yName="value",fill = "Clonality",
                        groupName="Treatment_broken_down",
                        groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue'), 'khaki1'), .8), showLegend=FALSE,
                        backgroundColor="white", xtitle="", ytitle="Mutation burden [SNVs/Mb]", 
                        mainTitle="Clonal and subclonal mutation burden",
                        addDot=T, dotSize=0.5,
                        legendPosition="bottom") + 
  facet_grid(~Clonality, scales = "free_x", space = "free_x") +
  theme(panel.spacing = unit(0.5, "cm"),
        strip.text = element_text(size = 12, face = "bold")) + 
  scale_y_continuous(trans='log10', limits = c(0.1,20000)) + 
  stat_compare_means(comparisons=list(c('B','C_NE'),c('B','C_RE'), c('A','C_NE'),c('A','C_RE')), label = "p.adj",  method = "wilcox.test",  p.adjust.method= "bonferroni", size = 4.5, face="bold", 
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.1), symbols = c("****", "***", "**", "*", "ns"))) +
  scale_x_discrete(labels=c('A', 'B', 'C \n NRP', 'C \n REP')) 
q <- p + my_theme + theme(plot.title = element_text(color="black", size=14, face="bold")) + theme(axis.title.x = element_text(color="black", size=16, face="bold"))+ theme(axis.title.y = element_text(color="black", size=14, face="bold")) 
plot(q)
dev.off()
