

# This scripts uses cell count data for immune metacluster IC 1-24 and makes following plots: 
# 1.) Individual plots for each immune metaclusters, plotting cell counts of respective metacluster 
# 2.) Joint plots for immune metaclusters subtypes plotting cell counts
library(ggplot2); 
require(ggpubr); 
require(scales); 
require(reshape2); 
library(plyr); 
library(easyGgplot2); 
library(ggpubr); 
library(stringr)

# prepare output dirs:
fig_dir1a = "~/analysis/IMC/metacluster_analyses/CD45_cells/plots/cell_counts/Per_Responsiveness"
dir.create(fig_dir1a, showWarnings=FALSE, recursive=TRUE)
fig_dir1b = "~/analysis/IMC/metacluster_analyses/CD45_cells/plots/cell_counts/Per_Timepoint"
dir.create(fig_dir1b, showWarnings=FALSE, recursive=TRUE)
fig_dir2 = "~/analysis/IMC/metacluster_analyses/CD45_cells/plots/cell_counts/Cluster_subtypes"
dir.create(fig_dir2, showWarnings=FALSE, recursive=TRUE)

# My themes
theme_mypub <- function(base_size = 14,
                        base_family = ""){
  theme_bw(base_size = base_size, 
           base_family = base_family) %+replace%
    theme(
      panel.grid.minor = element_blank(),   
      
      complete = TRUE
    )
}

my_theme_adapted <- theme(
  plot.title = element_text(color="black", size=18, face="bold"),
  axis.title.x = element_text(color="black", size=18, face="bold"),
  axis.text.x = element_text(color="black", size=18, face="bold"),
  axis.title.y = element_text(color="black", size=18, face="bold"),
  axis.text.y = element_text(color="black", size=18, face="bold"),
  legend.title = element_text(color = "black", size = 14, face="bold"),
  legend.text = element_text(color = "black", size = 14, face="bold"),
  
)

# Read in CD45 metaclusters (phenograph clustering performed in OMIQ)
Metacluster_df <- read.csv("~/analysis/IMC/input_files/CD45_metacluster_MEMORI.csv", header = T, stringsAsFactors = F, sep = ",")
colnames(Metacluster_df)[3:27] <- c(paste0("IC",c(1:24)), "CD45")
clusters <- c(paste0("IC", c(1:24)), "CD45")

# Read in area of ROIs calculated in Fiji
area_df <- read.csv("~/analysis/IMC/input_files/Area_exploratory_set_NucCyto.csv")
area_df <- area_df[,c("SampleID", "Final_area")]

# Calculate cell count per mm2
Metacluster_df <- merge(Metacluster_df, area_df, by.x = "MEMORI_Sample_ID", by.y = "SampleID")
for(c in clusters) {Metacluster_df[,c] <- as.numeric(Metacluster_df[,c])}; Metacluster_df$Final_area <- as.numeric(Metacluster_df$Final_area)
for(c in clusters) {Metacluster_df[,c] <- (Metacluster_df[,c]/Metacluster_df$Final_area*10^6) + 1}


## 1a.) Metacluster cell counts per Responsiveness and timepoint
my_comparisons <- list(c("A", "B"), c("B", "C"), c("A", "C"))
for (i in names(Metacluster_df)[4:28]){
  cluster <- i

pdf(paste0(fig_dir1a, "/" , cluster, "_cell_counts_per_responsiveness_TP.pdf"), height = 5, width = 5)
p <- ggplot2.violinplot(data=Metacluster_df, xName="Timepoint",yName=i, #fill = "Clonality",
                        groupName="Timepoint",
                        groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue')), .8), showLegend=FALSE,
                        backgroundColor="white", xtitle="", ytitle="cell count/mm2", 
                        mainTitle=paste0("Metacluster ", cluster),
                        addDot=T, dotSize=0.5,
                        legendPosition="bottom") + 
  facet_grid(~Responsiveness, scales = "free_x", space = "free_x", labeller = as_labeller(c(`NE` = "NRP",`RE` = "REP"))) +
  theme(panel.spacing = unit(0.5, "cm"),
        strip.text = element_text(size = 16, face = "bold")) + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size = 5, face="bold") + 
  scale_y_continuous(trans='log10', limits = c(1,1000000), labels = trans_format("log10", math_format(10^.x))) +
  ylab(expression("cell count/"~mm^2))
q <- p + my_theme_adapted 
plot(q)
dev.off()


}

## 1b.) Metacluster cell counts per TP 
my_comparisons_TP <- list(c("RE", "NE"))
for (i in names(Metacluster_df)[4:28]){
  cluster <- i
  
  pdf(paste0(fig_dir1b, "/" , cluster, "_cell_counts_per_TP.pdf"), height = 5, width = 5)
  p <- ggplot2.violinplot(data=Metacluster_df, xName="Responsiveness",yName=i, #fill = "Clonality",
                          groupName="Responsiveness",
                          groupColors= alpha(c("#b9adde","#26A69A"), .8), showLegend=FALSE,
                          backgroundColor="white", xtitle="", ytitle="cell count/mm2", 
                          mainTitle=paste0("Metacluster ", cluster),
                          addDot=T, dotSize=0.5,
                          legendPosition="bottom") + 
    facet_grid(~Timepoint, scales = "free_x", space = "free_x") +
    theme(panel.spacing = unit(0.5, "cm"),
          strip.text = element_text(size = 16, face = "bold")) + 
    stat_compare_means(comparisons = my_comparisons_TP, method = "wilcox.test", size = 5, face="bold") + 
    scale_y_continuous(trans='log10', limits = c(1,1000000), labels = trans_format("log10", math_format(10^.x))) +
    scale_x_discrete(labels= c("NRP", "REP")) +
    ylab(expression("cell count/"~mm^2))
  q <- p + my_theme_adapted 
  plot(q)
  dev.off()
  
  
  }
  
  
##################################################################
# 2.) Joint plots for metacluster subtypes

# Read in data frame where clusters are assigned to cell types and cluster subtypes based on marker expression (see heatmap of marker expression in C1-C28 from OMIQ)
cluster_assignment_df <- read.table("~/analysis/IMC/input_files/CD45_metacluster_assignments.txt", header = T, stringsAsFactors = F, sep = "\t")
cols <- c("A" = '#e63e00', "B" = "darkseagreen3", "C" = 'steelblue')
  
# CD8 cell cluster 
  CD8_cluster <- cluster_assignment_df[which(cluster_assignment_df$Biology == "activated_CD8" | cluster_assignment_df$Biology == "exhausted_CD8" | cluster_assignment_df$Biology == "CD8"), c("cluster_name")]
  Metacluster_CD8 <- Metacluster_df[, colnames(Metacluster_df) %in% CD8_cluster]
  Metacluster_CD8$Responsiveness <- Metacluster_df$Responsiveness
  Metacluster_CD8$Timepoint <- Metacluster_df$Timepoint
  Metacluster_CD8_2 <- melt(Metacluster_CD8)
  

  pdf(paste0(fig_dir2, "/CD8_clusters.pdf"), height = 3, width = 13)
  q <- ggplot(Metacluster_CD8_2, aes(x = variable, y = value)) +
    geom_violin(
      aes(color = Timepoint), trim = FALSE,
      position = position_dodge(0.6) 
    ) +
    geom_boxplot(
      aes(color = Timepoint), width = 0.15,
      position = position_dodge(0.6)
    ) +
    scale_color_manual(values = cols) + 
    xlab("Cluster") +
    ylab(expression("cell count/"~mm^2)) +
    facet_grid(~Responsiveness, scales = "free_x", space = "free_x", labeller = as_labeller(c(`NE` = "NRP",`RE` = "REP"))) +
    theme(panel.spacing = unit(0.5, "cm"),
          strip.text.x = element_text(size = 30, face = "bold")) + 
    stat_compare_means(aes(comparison = Timepoint), method = "wilcox", label = "p.signif", label.y = 6, size = 6, face = "bold",  hide.ns = T) +  
    scale_y_continuous(trans='log10', limits = c(1,5000000), labels = trans_format("log10", math_format(10^.x)))+ 
    theme_mypub() + my_theme_adapted + theme(axis.text.x = element_text(color="black", size=16, face="bold"))
  plot(q) + theme(strip.text.x = element_text(size = 30, face = "bold"))
  dev.off()
  # CD8 cell cluster
  
  
# Other T-cell cluster  
  Other_t_cell_cluster <- cluster_assignment_df[which(cluster_assignment_df$Biology == "other_t_cells"), c("cluster_name")]
  Metacluster_t_cells <- Metacluster_df[, colnames(Metacluster_df) %in% Other_t_cell_cluster ]
  Metacluster_t_cells$Responsiveness <- Metacluster_df$Responsiveness
  Metacluster_t_cells$Timepoint <- Metacluster_df$Timepoint
  Metacluster_t_cells2 <- melt(Metacluster_t_cells)

  pdf(paste0(fig_dir2, "/Other_t_cell_clusters.pdf"), height = 3, width = 13)
  q <- ggplot(Metacluster_t_cells2, aes(x = variable, y = value)) +
    geom_violin(
      aes(color = Timepoint), trim = FALSE,
      position = position_dodge(0.6) 
    ) +
    geom_boxplot(
      aes(color = Timepoint), width = 0.15,
      position = position_dodge(0.6)
    ) +
    scale_color_manual(values = cols) + 
    xlab("Cluster") +
    ylab(expression("cell count/"~mm^2)) +
    facet_grid(~Responsiveness, scales = "free_x", space = "free_x", labeller = as_labeller(c(`NE` = "NRP",`RE` = "REP"))) +
    theme(panel.spacing = unit(0.5, "cm"),
          strip.text.x = element_text(size = 30, face = "bold")) + 
    stat_compare_means(aes(comparison = Timepoint), method = "wilcox", label = "p.signif", label.y = 6, size = 6, face = "bold",  hide.ns = T) +  
    scale_y_continuous(trans='log10', limits = c(1,5000000), labels = trans_format("log10", math_format(10^.x)))+ 
    theme_mypub() + my_theme_adapted + theme(axis.text.x = element_text(color="black", size=16, face="bold"))
  plot(q) + theme(strip.text.x = element_text(size = 30, face = "bold"))
  dev.off()
  
# Mixed cluster with granulocytes, macrophages and myeloid cells    
  granulo_macro_myeloid_cluster <- cluster_assignment_df[which(cluster_assignment_df$Biology == "granulocytes" | cluster_assignment_df$Biology == "macrophages" | cluster_assignment_df$Biology == "myleoid_cells"), c("cluster_name")]
  Metacluster_granulo_macro_myeloid_cluster <- Metacluster_df[, colnames(Metacluster_df) %in% granulo_macro_myeloid_cluster ]
  Metacluster_granulo_macro_myeloid_cluster$Responsiveness <- Metacluster_df$Responsiveness
  Metacluster_granulo_macro_myeloid_cluster$Timepoint <- Metacluster_df$Timepoint
  Metacluster_granulo_macro_myeloid_cluster2 <- melt(Metacluster_granulo_macro_myeloid_cluster)

  pdf(paste0(fig_dir2, "/Granulo_macroph_myeloid_clusters.pdf"), height = 3, width = 17)
  q <- ggplot(Metacluster_granulo_macro_myeloid_cluster2, aes(x = variable, y = value)) +
    geom_violin(
      aes(color = Timepoint), trim = FALSE,
      position = position_dodge(0.6) 
    ) +
    geom_boxplot(
      aes(color = Timepoint), width = 0.15,
      position = position_dodge(0.6)
    ) +
    scale_color_manual(values = cols) + 
    xlab("Cluster") +
    ylab(expression("cell count/"~mm^2)) +
    facet_grid(~Responsiveness, scales = "free_x", space = "free_x", labeller = as_labeller(c(`NE` = "NRP",`RE` = "REP"))) +
    theme(panel.spacing = unit(0.5, "cm"),
          strip.text.x = element_text(size = 30, face = "bold")) + 
    stat_compare_means(aes(comparison = Timepoint), method = "wilcox", label = "p.signif", label.y = 6, size = 6, face = "bold",  hide.ns = T) +  
    scale_y_continuous(trans='log10', limits = c(1,5000000), labels = trans_format("log10", math_format(10^.x)))+ 
    theme_mypub() + my_theme_adapted + theme(axis.text.x = element_text(color="black", size=16, face="bold"))
  plot(q) + theme(strip.text.x = element_text(size = 30, face = "bold"))
  dev.off()
  
  
  
