

# This scripts uses subsampled CD4 and CD8 cells and subtypes and plots: 
# 1.) Cell counts of CD4 or CD8 subtypes
# 2.) Ratio of GzmB+/PD1+ and GzmB+/Eomes+ CD4 and CD8 cells
# 3.) Fraction of CD4 or CD8 subtypes of total CD4 or CD8 cells

library(ggplot2); 
require(ggpubr); 
require(scales); 
require(reshape2); 
library(plyr); 
library(easyGgplot2); 
library(stringr)

# prepare output dirs:
fig_dir1a = "~/analysis/IMC/CD4_CD8_subsampled_cells/plots/cell_counts/Per_Responsiveness"
dir.create(fig_dir1a, showWarnings=FALSE, recursive=TRUE)
fig_dir1b = "~/analysis/IMC/CD4_CD8_subsampled_cells/plots/cell_counts/Per_Timepoint"
dir.create(fig_dir1b, showWarnings=FALSE, recursive=TRUE)
fig_dir1c = "~/analysis/IMC/CD4_CD8_subsampled_cells/plots/plots/cell_counts/Per_Treatment_broken_down"
dir.create(fig_dir1c, showWarnings=FALSE, recursive=TRUE)
fig_dir1d = "~/analysis/IMC/CD4_CD8_subsampled_cells/plots/plots/cell_counts/Per_Regression_grade"
dir.create(fig_dir1d, showWarnings=FALSE, recursive=TRUE)

fig_dir3a= "~/analysis/IMC/CD4_CD8_subsampled_cells/plots/cell_fractions/Per_Responsiveness"
dir.create(fig_dir3a, showWarnings=FALSE, recursive=TRUE)
fig_dir3b = "~/analysis/IMC/CD4_CD8_subsampled_cells/plots/cell_fractions/Per_Timepoint"
dir.create(fig_dir3b, showWarnings=FALSE, recursive=TRUE)
fig_dir3c = "~/analysis/IMC/CD4_CD8_subsampled_cells/plots/plots/cell_fractions/Per_Treatment_broken_down"
dir.create(fig_dir3c, showWarnings=FALSE, recursive=TRUE)
fig_dir3d = "~/analysis/IMC/CD4_CD8_subsampled_cells/plots/plots/cell_fractions/Per_Regression_grade"
dir.create(fig_dir3d, showWarnings=FALSE, recursive=TRUE)

# My theme 
theme_mypub <- function(base_size = 14,
                        base_family = ""){
  theme_bw(base_size = base_size, 
           base_family = base_family) %+replace%
    theme(
      panel.grid.major = element_blank(),   
      panel.grid.minor = element_blank(),   
      
      complete = TRUE
    )
}

my_theme_adapted <- theme(
  plot.title = element_text(color="black", size=20, face="bold"),
  axis.title.x = element_text(color="black", size=22, face="bold"),
  axis.text.x = element_text(color="black", size=18, face="bold"),
  axis.title.y = element_text(color="black", size=18, face="bold"),
  axis.text.y = element_text(color="black", size=18, face="bold"),
  legend.title = element_text(color = "black", size = 14, face="bold"),
  legend.text = element_text(color = "black", size = 16, face="bold"),
  
)

########################################
# 1) Read in subsampled CD4 or CD8 cells (subsampling performed in OMIQ) 
subsampled_immune_cells_df <- read.table("~/analysis/IMC/input_files/CD8_CD4_subsets_cell_counts_exploratory_set.csv", header = T, stringsAsFactors = F, sep = ",")

subsampled_immune_cells_df$Treatment_broken_down <-  ifelse(subsampled_immune_cells_df$Timepoint == "A", "A", 
                                                                 ifelse(subsampled_immune_cells_df$Timepoint == "B", "B",
                                                                        ifelse ((subsampled_immune_cells_df$Timepoint == "C" & subsampled_immune_cells_df$Responsiveness == "RE"), "C_RE",
                                                                                ifelse ((subsampled_immune_cells_df$Timepoint == "C" & subsampled_immune_cells_df$Responsiveness == "NE"), "C_NE", "XX"))))

  
# Read in area of ROIs calculated in Fiji
area_df <- read.csv("~/analysis/IMC/input_files/Area_exploratory_set_NucCyto.csv")
area_df <- area_df[,c("SampleID", "Final_area")]

subsampled_immune_cells_df <- merge(subsampled_immune_cells_df, area_df, by.x = "MEMORI_Sample_ID", by.y = "SampleID")
cell_columns <- c("Unfiltered.count", colnames((subsampled_immune_cells_df)[ , grepl( "CD" , names( subsampled_immune_cells_df ) ) ]))

# Calculation of immune cells per mm2
for(c in cell_columns) {subsampled_immune_cells_df[,c] <- as.numeric(subsampled_immune_cells_df[,c])}; subsampled_immune_cells_df$Final_area <- as.numeric(subsampled_immune_cells_df$Final_area)
for(c in cell_columns) {subsampled_immune_cells_df[,c] <- subsampled_immune_cells_df[,c]/subsampled_immune_cells_df$Final_area*10^6}
for(c in cell_columns) {subsampled_immune_cells_df[,c] <- as.numeric(subsampled_immune_cells_df[,c])}
subsampled_immune_cells_df[,c(3:24)] <- subsampled_immune_cells_df[,c(3:24)] +1


# 1.) Cell count plots 
my_comparisons <- list(c("A", "B"), c("B", "C"), c("A", "C"))

for (i in cell_columns){
  immune_cell <- str_sub(i, end=-8)
 Other_info_list = strsplit(immune_cell, ".", fixed = TRUE)
 extractinfo <- function(lst,n) { sapply(lst,'[', n)}
 immune_marker1 <-  extractinfo(Other_info_list ,1); immune_marker2 <-  extractinfo(Other_info_list ,2)

 # 1a.) Plot CD4 and CD8 subsets per Responsiveness and timepoint
  pdf(paste0(fig_dir1a, "/", immune_cell, "_cell_counts_per_responsiveness_TP.pdf"), height = 5, width = 5)
  p <- ggplot2.violinplot(data=subsampled_immune_cells_df, xName="Timepoint",yName=i,
                          groupName="Timepoint",
                          groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue')), .8), showLegend=FALSE,
                          backgroundColor="white", xtitle="", ytitle="cell count", 
                          mainTitle=bquote(bold(paste(.(immune_marker1)^"+"~.(immune_marker2)^"+" ~"T-cells"))),
                          addDot=T, dotSize=0.5,
                          legendPosition="bottom") + 
    facet_grid(~Responsiveness, scales = "free_x", space = "free_x", labeller = as_labeller(c(`NE` = "NRP",`RE` = "REP"))) +
    theme(panel.spacing = unit(0.5, "cm"),
          strip.text = element_text(size = 16, face = "bold")) + 
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", p.adjust.method= "bonferroni", size = 6, face="bold", label = "p.signif",
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.1), 
                                          symbols = c("****", "***", "**", "*", "ns"))) + 
    scale_y_continuous(trans='log10', limits = c(1,10000), labels = trans_format("log10", math_format(10^.x))) + 
    xlab("") + ylab(expression("cell count/"~mm^2))
  q <- p + my_theme_adapted 
  plot(q)
  dev.off()
  
  
  # 1b.) Plot per Timepoint and compare Responsiveness  
  RE_NR_comparison <- list(c("NE", "RE"))
  pdf(paste0(fig_dir1b, "/", immune_cell, "_cell_counts_per_TP.pdf"), height = 5, width = 5)
  p <- ggplot2.violinplot(subsampled_immune_cells_df, xName="Responsiveness",yName=i,
                          groupName="Responsiveness",
                          groupColors= alpha(c('#b9adde','#26A69A'), .8), showLegend=FALSE,
                          backgroundColor="white", xtitle="", ytitle="cell count", 
                          mainTitle=bquote(bold(paste(.(immune_marker1)^"+"~.(immune_marker2)^"+" ~"T-cells"))),
                          addDot=T, dotSize=0.5,
                          legendPosition="bottom") + 
    facet_grid(~Timepoint, scales = "free_x", space = "free_x") +
    theme(panel.spacing = unit(0.5, "cm"),
          strip.text = element_text(size = 16, face = "bold")) + 
    stat_compare_means(comparisons = RE_NR_comparison, method = "wilcox.test", size = 6, face="bold",label = "p.signif",
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.1), 
                                          symbols = c("****", "***", "**", "*", "ns"))) + 
    scale_y_continuous(trans='log10', limits = c(1,10000), labels = trans_format("log10", math_format(10^.x))) +
    xlab("") + ylab(expression("cell count/"~mm^2)) + scale_x_discrete(labels=c("NE"= "NRP", "RE" = "REP"))
  q <- p + my_theme_adapted + theme(axis.text.x = element_text(color="black", size=16, face="bold"))
  plot(q)
  dev.off()
  
  
  # 1c.) Plot for A, B and C broken down by treatment
  pdf(paste0(fig_dir1c, "/", immune_cell, "_cell_counts_by_treatment.pdf"), width = 5, height = 5)
  p <- ggplot2.violinplot(data=subsampled_immune_cells_df, xName="Treatment_broken_down",yName=i,
                          groupName="Treatment_broken_down",
                          groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue'), 'khaki1'), .8), showLegend=FALSE,
                          backgroundColor="white", xtitle="", ytitle="cell count", 
                          mainTitle=bquote(bold(paste(.(immune_marker1)^"+"~.(immune_marker2)^"+" ~"T-cells"))),
                          addDot=T, dotSize=0.5,
                          legendPosition="bottom") + 
    stat_compare_means(comparisons = list(c('B','C_NE'),c('B','C_RE'), c('A','C_NE'),c('A','C_RE')), method = "wilcox.test",  p.adjust.method= "bonferroni", size = 6, face="bold", label = "p.signif",
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.1), 
                                          symbols = c("****", "***", "**", "*", "ns"))) + 
    scale_y_continuous(trans='log10', limits = c(1,10000), labels = trans_format("log10", math_format(10^.x))) + 
    scale_x_discrete(labels=c('A', 'B', 'C \n (NRP)', 'C \n (REP)')) +
    xlab("") + ylab(expression("cell count/"~mm^2))
  q <- p + my_theme_adapted 
  plot(q)
  dev.off()
  
  
  
  
  # 1d.) Plot per regression grade and TP
  pdf(paste0(fig_dir1d, "/", immune_cell, "_cell_counts_per_remission_TP.pdf"), height = 5, width = 5)
  p <- ggplot2.violinplot(subsampled_immune_cells_df, xName="Timepoint",yName=i,
                          groupName="Timepoint",
                          groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue')), .8), showLegend=FALSE,
                          backgroundColor="white", xtitle="", ytitle="cell count", 
                          mainTitle=bquote(bold(paste(.(immune_marker1)^"+"~.(immune_marker2)^"+" ~"T-cells"))),
                          addDot=T, dotSize=0.5,
                          legendPosition="bottom") + 
    facet_grid(~simplified_regression_score, scales = "free_x", space = "free_x") +
    theme(panel.spacing = unit(0.5, "cm"),
          strip.text = element_text(size = 12, face = "bold")) + 
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size = 4, face="bold", label = "p.signif",
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.1), 
                                          symbols = c("****", "***", "**", "*", "ns"))) + 
    scale_y_continuous(trans='log10', limits = c(1,10000), labels = trans_format("log10", math_format(10^.x)))
  q <- p + my_theme_adapted + theme(plot.title = element_text(color="black", size=14, face="bold")) + theme(axis.title.x = element_text(color="black", size=14, face="bold"))+ theme(axis.title.y = element_text(color="black", size=14, face="bold")) +
    xlab("") + ylab(expression("cell count/"~mm^2))
  plot(q)
  dev.off()
  
  
}

# 2.) Calculate the ratio of GzmB+/PD1+ and GzmB+/Eomes+ CD4 and CD8 cells
subsampled_immune_cells_df$ratio_CD4GzmB_CD4PD1 <- subsampled_immune_cells_df$CD4.GzmB..count/subsampled_immune_cells_df$CD4.PD1..count
subsampled_immune_cells_df$ratio_CD8GzmB_CD8PD1 <- subsampled_immune_cells_df$CD8.GzmB..count/subsampled_immune_cells_df$CD8.PD1..count
subsampled_immune_cells_df$ratio_CD4GzmB_CD4Eomes  <- subsampled_immune_cells_df$CD4.GzmB..count/subsampled_immune_cells_df$CD4.Eomes..count
subsampled_immune_cells_df$ratio_CD8GzmB_CD8Eomes  <- subsampled_immune_cells_df$CD8.GzmB..count/subsampled_immune_cells_df$CD8.Eomes..count

for (i in names(subsampled_immune_cells_df)[c(33:36)]){
  cell_type <- sapply(i, function(x) substr(x, 7,9))
  
  # 2a.) Plot per Responsiveness and compare timepoints
  pdf(paste0(fig_dir1a, "/", i, "_cell_counts_per_responsiveness_TP.pdf"), height = 5, width = 5)
  p <- ggplot2.violinplot(data=subsampled_immune_cells_df, xName="Timepoint",yName=i,
                          groupName="Timepoint",
                          groupColors= c('#e63e00','darkseagreen3',('steelblue')), showLegend=FALSE,
                          backgroundColor="white", xtitle="", ytitle="ratio", 
                          mainTitle=i,
                          addDot=T, dotSize=0.5,
                          legendPosition="bottom") + 
    facet_grid(~Responsiveness, scales = "free_x", space = "free_x", labeller = as_labeller(c(`NE` = "NRP",`RE` = "REP"))) +
    theme(panel.spacing = unit(0.5, "cm"),
          strip.text = element_text(size = 16, face = "bold")) + 
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size = 6, face="bold", vjust = 0, label = "p.signif", 
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.1), 
                                          symbols = c("****", "***", "**", "*", "ns"))) + 
    scale_y_continuous(trans='log10', limits = c(0.1,180), labels = trans_format("log10", math_format(10^.x))) + 
    xlab("") + ylab("ratio")
  q <- p + my_theme_adapted 
  plot(q)
  dev.off()
  
  
  # 2b.) Plot per Timepoint and compare Responsiveness  
  RE_NR_comparison <- list(c("NE", "RE"))
  pdf(paste0(fig_dir1b, "/", i, "_cell_counts_per_TP.pdf"), height = 5, width = 5)
  p <- ggplot2.violinplot(subsampled_immune_cells_df, xName="Responsiveness",yName=i,
                        groupName="Responsiveness",
                        groupColors= c('#b9adde','#26A69A'), showLegend=FALSE,
                        backgroundColor="white", xtitle="", ytitle="ratio", 
                        mainTitle=i,
                        addDot=T, dotSize=0.5,
                        legendPosition="bottom") + 
  facet_grid(~Timepoint, scales = "free_x", space = "free_x") +
  theme(panel.spacing = unit(0.5, "cm"),
        strip.text = element_text(size = 16, face = "bold")) + 
  stat_compare_means(comparisons = RE_NR_comparison, method = "wilcox.test", size = 6, face="bold", vjust = 0, label = "p.signif", 
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.1), 
                                        symbols = c("****", "***", "**", "*", "ns"))) + 
  scale_y_continuous(trans='log10', limits = c(0.1,180)) + scale_x_discrete(labels=c("NE"= "NRP", "RE" = "REP"))
  q <- p + my_theme_adapted + theme(axis.text.x = element_text(color="black", size=16, face="bold"))
  plot(q)
  dev.off()

  # 2d.) Plot per Timepoint and compare Remission score  
  pdf(paste0(fig_dir1d, "/", i, "_cell_counts_per_remission_TP.pdf"), height = 5, width = 5)
  p <- ggplot2.violinplot(subsampled_immune_cells_df, xName="Timepoint",yName=i,
                        groupName="Timepoint",
                        groupColors= c('#e63e00','darkseagreen3','steelblue'), showLegend=FALSE,
                        backgroundColor="white", xtitle="", ytitle="ratio", 
                        mainTitle=i,
                        addDot=T, dotSize=0.5,
                        legendPosition="bottom") + 
  facet_grid(~simplified_regression_score, scales = "free_x", space = "free_x") +
  theme(panel.spacing = unit(0.5, "cm"),
        strip.text = element_text(size = 16, face = "bold")) + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size = 6, face="bold", vjust = 0, label = "p.signif", 
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.1), 
                                        symbols = c("****", "***", "**", "*", "ns"))) + 
  scale_y_continuous(trans='log10', limits = c(0.1,800))
  q <- p + my_theme_adapted 
  plot(q)
  dev.off()



}


################################################
# 3.) Fraction of subsets per CD4 or CD8 parent

# Read in percentage of CD4 or CD8 subtypes from total CD4 or CD8 cell counts (subsampling performed in OMIQ)  
fraction_subsampled_immune_cells_df <- read.table("~/analysis/IMC/input_files/CD8_CD4_subsets_fraction_of_parent_exploratory_set.csv", header = T, stringsAsFactors = F, sep = ",")
fraction_subsampled_immune_cells_df$Treatment_broken_down <-  ifelse(fraction_subsampled_immune_cells_df$Timepoint == "A", "A", 
                                                            ifelse(fraction_subsampled_immune_cells_df$Timepoint == "B", "B",
                                                                   ifelse ((fraction_subsampled_immune_cells_df$Timepoint == "C" & fraction_subsampled_immune_cells_df$Responsiveness == "RE"), "C_RE",
                                                                           ifelse ((fraction_subsampled_immune_cells_df$Timepoint == "C" & fraction_subsampled_immune_cells_df$Responsiveness == "NE"), "C_NE", "XX"))))
fraction_subsampled_immune_cells_df[,c(3:24)] <- fraction_subsampled_immune_cells_df[,c(3:24)] +1


for (i in cell_columns){
  immune_cell <- str_sub(i, end=-17)
  Other_info_list = strsplit(immune_cell, ".", fixed = TRUE)
  extractinfo <- function(lst,n) { sapply(lst,'[', n)}
  immune_marker1 <-  extractinfo(Other_info_list ,1); immune_marker2 <-  extractinfo(Other_info_list ,2)
  
 
  # 3a.) Plot per Responsiveness and compare TP  
  my_comparisons <- list(c("A", "B"), c("B", "C"), c("A", "C"))
  pdf(paste0(fig_dir3a, "/", immune_cell, "_fractions_per_responsiveness_TP.pdf"), height = 5, width = 5)
  p <- ggplot2.violinplot(fraction_subsampled_immune_cells_df, xName="Timepoint",yName=i,
                          groupName="Timepoint",
                          groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue')), .8), showLegend=FALSE,
                          backgroundColor="white", xtitle="", ytitle= bquote(bold(paste("Fraction of"~ .(immune_marker1)^"+" ~ "cells [%]"))),
                          mainTitle=bquote(bold(paste(.(immune_marker1)^"+"~.(immune_marker2)^"+" ~"T-cells"))),
                          addDot=T, dotSize=0.5,
                          legendPosition="bottom") + 
    facet_grid(~Responsiveness, scales = "free_x", space = "free_x", labeller = as_labeller(c(`NE` = "NRP",`RE` = "REP"))) +
    theme(panel.spacing = unit(0.5, "cm"),
          strip.text = element_text(size = 16, face = "bold")) + 
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size = 6, face="bold", label = "p.signif",
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.1), 
                                          symbols = c("****", "***", "**", "*", "ns"))) + 
    scale_y_continuous(trans='log10', limits = c(-1,500))
  q <- p + my_theme_adapted 
  plot(q)
  dev.off()
  
  # 3b.) Plot per Timepoint and compare Responsiveness  
  RE_NR_comparison <- list(c("NE", "RE"))
  pdf(paste0(fig_dir3b, "/", immune_cell, "_fractions_per_TP.pdf"), height = 5, width = 5)
  p <- ggplot2.violinplot(fraction_subsampled_immune_cells_df, xName="Responsiveness",yName=i,
                          groupName="Responsiveness",
                          groupColors= alpha(c('orange','green'), .8), showLegend=FALSE,
                          backgroundColor="white", xtitle="", ytitle=bquote(bold(paste("Fraction of"~ .(immune_marker1)^"+" ~ "cells [%]"))), 
                          mainTitle=bquote(bold(paste(.(immune_marker1)^"+"~.(immune_marker2)^"+" ~"T-cells"))),
                          addDot=T, dotSize=0.5,
                          legendPosition="bottom") + 
    facet_grid(~Timepoint, scales = "free_x", space = "free_x") +
    theme(panel.spacing = unit(0.5, "cm"),
          strip.text = element_text(size = 16, face = "bold")) + 
    stat_compare_means(comparisons = RE_NR_comparison, method = "wilcox.test", size = 6, face="bold", label = "p.signif",
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.1), 
                                          symbols = c("****", "***", "**", "*", "ns"))) + 
    scale_y_continuous(trans='log10', limits = c(-1,500))
  q <- p + my_theme_adapted  
  plot(q)
  dev.off()
  
  # 3c.) Plot for A, B and C broken down by treatment
  pdf(paste0(fig_dir3c, "/", immune_cell, "_cell_counts_by_treatment.pdf"), width = 5, height = 5)
  p <- ggplot2.violinplot(data=fraction_subsampled_immune_cells_df, xName="Treatment_broken_down",yName=i,
                          groupName="Treatment_broken_down",
                          groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue'), 'khaki1'), .8), showLegend=FALSE,
                          backgroundColor="white", xtitle="", ytitle=bquote(bold(paste("Fraction of"~ .(immune_marker1)^"+" ~ "cells [%]"))), 
                          mainTitle=bquote(bold(paste(.(immune_marker1)^"+"~.(immune_marker2)^"+" ~"T-cells"))),
                          addDot=T, dotSize=0.5,
                          legendPosition="bottom") + 
    stat_compare_means(comparisons = list(c('B','C_NE'),c('B','C_RE'), c('A','C_NE'),c('A','C_RE')), method = "wilcox.test",  p.adjust.method= "bonferroni", size = 6, face="bold", label = "p.signif",
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.1), 
                                          symbols = c("****", "***", "**", "*", "ns"))) + 
    scale_y_continuous(trans='log10', limits = c(-1,500)) + 
    scale_x_discrete(labels=c('A', 'B', 'C \n (NRP)', 'C \n (REP)')) +
    xlab("")
  q <- p + my_theme_adapted 
  plot(q)
  dev.off()
  
 
  # 3d.) Plot per Remission grade and TP
  pdf(paste0(fig_dir3d, "/", immune_cell, "_fractions_per_remission_TP.pdf"), height = 5, width = 5)
  p <- ggplot2.violinplot(fraction_subsampled_immune_cells_df, xName="Timepoint",yName=i,
                          groupName="Timepoint",
                          groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue')), .8), showLegend=FALSE,
                          backgroundColor="white", xtitle="", ytitle=bquote(bold(paste("Fraction of"~ .(immune_marker1)^"+" ~ "cells [%]"))), 
                          mainTitle= bquote(bold(paste(.(immune_marker1)^"+"~.(immune_marker2)^"+" ~"T-cells"))),
                          addDot=T, dotSize=0.5,
                          legendPosition="bottom") + 
    facet_grid(~simplified_regression_score, scales = "free_x", space = "free_x") +
    theme(panel.spacing = unit(0.5, "cm"),
          strip.text = element_text(size = 16, face = "bold")) + 
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size = 6, face="bold", label = "p.signif",
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.1), 
                                          symbols = c("****", "***", "**", "*", "ns"))) + 
    scale_y_continuous(trans='log10', limits = c(-1,500))
  q <- p + my_theme_adapted 
  plot(q)
  dev.off()
  
}