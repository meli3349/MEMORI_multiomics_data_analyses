

# This scripts uses cell count data for metacluster C1-28 from OMIQ analyses and makes following plots: 
# 1.) Individual plots for each metaclusters, plotting cell counts of respective metacluster 
# 2.) Joint plots for metaclusters subtypes (immune cell, tumour cell, other cell type) plotting cell counts

library(ggplot2); 
require(ggpubr); 
require(scales); 
require(reshape2); 
library(plyr); 
library(easyGgplot2); 
library(ggpubr); 
library(stringr)

# prepare output dirs:
fig_dir1a = "~/analysis/IMC/metacluster_analyses/all_cells/plots/cell_counts/Per_Responsiveness"
dir.create(fig_dir1a, showWarnings=FALSE, recursive=TRUE)
fig_dir1b = "~/analysis/IMC/metacluster_analyses/all_cells/plots/cell_counts/Per_Timepoint"
dir.create(fig_dir1b, showWarnings=FALSE, recursive=TRUE)
fig_dir1c = "~/analysis/IMC/metacluster_analyses/all_cells/plots/cell_counts/Per_Treatment_broken_down"
dir.create(fig_dir1c, showWarnings=FALSE, recursive=TRUE)
fig_dir2 = "~/analysis/IMC/metacluster_analyses/all_cells/plots/cell_counts/Cluster_subtypes"
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

# Read in metaclusters (phenograph clustering performed in OMIQ)
Metacluster_df <- read.csv("~/analysis/IMC/input_files/All_cell_metcluster_MEMORI.csv", header = T, stringsAsFactors = F, sep = ",")
colnames(Metacluster_df)[3:30] <- paste0("C",c(1:28))
clusters <- paste0("C", c(1:28))
Metacluster_df$Treatment_broken_down <- ifelse(Metacluster_df$Timepoint == "A", "A", 
                                               ifelse(Metacluster_df$Timepoint == "B", "B",
                                                      ifelse ((Metacluster_df$Timepoint == "C" & Metacluster_df$Responsiveness == "RE"), "C_RE",
                                                              ifelse ((Metacluster_df$Timepoint == "C" & Metacluster_df$Responsiveness == "NE"), "C_NE", "XX"))))


# Read in area of ROIs calculated in Fiji
area_df <- read.csv("~/analysis/IMC/input_files/Area_exploratory_set_NucCyto.csv")
area_df <- area_df[,c("SampleID", "Final_area")]

# Calculate cell count per mm2
Metacluster_df <- merge(Metacluster_df, area_df, by.x = "MEMORI_Sample_ID", by.y = "SampleID")
for(c in clusters) {Metacluster_df[,c] <- as.numeric(Metacluster_df[,c])}; Metacluster_df$Final_area <- as.numeric(Metacluster_df$Final_area)
for(c in clusters) {Metacluster_df[,c] <- (Metacluster_df[,c]+1)/Metacluster_df$Final_area*10^6}



## 1a.) Plot cell counts for each metacluster per Responsiveness and timepoint
my_comparisons <- list(c("A", "B"), c("B", "C"), c("A", "C"))
for (i in names(Metacluster_df)[3:31]){
  cluster <- i

pdf(paste0(fig_dir1a, "/" , cluster, "_cell_counts_per_responsiveness_TP.pdf"), height = 5, width = 5)
p <- ggplot2.violinplot(data=Metacluster_df, xName="Timepoint",yName=i, 
                        groupName="Timepoint",
                        groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue')), .8), showLegend=FALSE,
                        backgroundColor="white", xtitle="", ytitle="cell count/mm2", 
                        mainTitle=paste0("Metacluster ", cluster),
                        addDot=T, dotSize=0.5,
                        legendPosition="bottom") + 
  facet_grid(~Responsiveness, scales = "free_x", space = "free_x", labeller = as_labeller(c(`NE` = "NRP",`RE` = "REP"))) +
  theme(panel.spacing = unit(0.5, "cm"),
        strip.text = element_text(size = 16, face = "bold")) + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", p.adjust.method= "bonferroni",size = 5, face="bold") + 
  scale_y_continuous(trans='log10', limits = c(0.09,1000000), labels = trans_format("log10", math_format(10^.x))) +
  ylab(expression("cell count/"~mm^2))
q <- p + my_theme_adapted 
plot(q)
dev.off()


}

## 1b.) Plot metacluster cell counts per Timepoint 
my_comparisons_TP <- list(c("RE", "NE"))
for (i in names(Metacluster_df)[3:31]){
  cluster <- i
  
  pdf(paste0(fig_dir1b, "/" , cluster, "_cell_counts_per_TP.pdf"), height = 5, width = 5)
  p <- ggplot2.violinplot(data=Metacluster_df, xName="Responsiveness",yName=i,
                          groupName="Responsiveness",
                          groupColors= alpha(c("#b9adde","#26A69A"), .8), showLegend=FALSE,
                          backgroundColor="white", xtitle="", ytitle="cell count/mm2", 
                          mainTitle=paste0("Metacluster ", cluster),
                          addDot=T, dotSize=0.5,
                          legendPosition="bottom") + 
    facet_grid(~Timepoint, scales = "free_x", space = "free_x") +
    theme(panel.spacing = unit(0.5, "cm"),
          strip.text = element_text(size = 16, face = "bold")) + 
    stat_compare_means(comparisons = my_comparisons_TP, method = "wilcox.test", p.adjust.method= "bonferroni", size = 5, face="bold") + 
    scale_y_continuous(trans='log10', limits = c(0.09,1000000), labels = trans_format("log10", math_format(10^.x))) +
    scale_x_discrete(labels= c("NRP", "REP")) +
    ylab(expression("cell count/"~mm^2))
  q <- p + my_theme_adapted 
  plot(q)
  dev.off()
  
  ## 1c) Plot metacluster cell counts for Timepoints A, B and C broken down by treatment
  pdf(paste0(fig_dir1c, "/" , cluster, "_cell_counts_per_treatment.pdf"), height = 5, width = 5)
  #pdf(paste0("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/CYTOF/Metacluster_analysis/2022.05.18_Phenograph_CA-CS/statistics/cell_counts/Treatment_broken_down/", cluster, "_cell_counts_per_treatment.pdf"), width = 5, height = 5)
  p <- ggplot2.violinplot(data=Metacluster_df, xName="Treatment_broken_down",yName=i,
                          groupName="Treatment_broken_down",
                          groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue'), 'khaki1'), .8), showLegend=FALSE,
                          backgroundColor="white", xtitle="", ytitle="cell count/mm2", 
                          mainTitle=paste0("Metacluster ", cluster),
                          addDot=T, dotSize=0.5) + 
    scale_y_continuous(trans='log10', limits = c(0.09,1000000), labels = trans_format("log10", math_format(10^.x))) +
    stat_compare_means(comparisons=list(c('B','C_NE'),c('B','C_RE'), c('A','C_NE'),c('A','C_RE')), method = "wilcox.test",  p.adjust.method= "bonferroni", size = 4.5, face="bold", 
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.1), symbols = c("****", "***", "**", "*", "ns"))) +
    scale_x_discrete(labels=c('A', 'B', 'C \n (NRP)', 'C \n (REP)'))
  q <- p +  my_theme_adapted  
  plot(q)
  dev.off()
  
  
}



## 2.) Joint plots for metacluster subtypes
# Read in data frame where clusters are assigned to cell types and cluster subtypes based on marker expression (see heatmap of marker expression in C1-C28 from OMIQ)
cluster_assignment_df <- read.table("~/analysis/IMC/input_files/All_cell_metacluster_assignments.txt", header = T, stringsAsFactors = F, sep = "\t")

tumour_cluster <- cluster_assignment_df[which(cluster_assignment_df$Cluster_assignment == "tumour_cluster"), c("cluster_name")]
immune_cluster <- cluster_assignment_df[which(cluster_assignment_df$Cluster_assignment == "immune_cells"), c("cluster_name")]
Other_cell_cluster <- cluster_assignment_df[which(cluster_assignment_df$Cluster_assignment == "other_cells"), c("cluster_name")]



######################################  
# Immune Cluster

Metacluster_immune <- Metacluster_df[, colnames(Metacluster_df) %in% immune_cluster]
Metacluster_immune$Responsiveness <- Metacluster_df$Responsiveness
Metacluster_immune$Timepoint <- Metacluster_df$Timepoint
Metacluster_immune2 <- melt(Metacluster_immune)
cols <- c("A" = '#e63e00', "B" = "darkseagreen3", "C" = 'steelblue')

pdf(paste0(fig_dir2, "/Immune_clusters.pdf"), height = 3, width = 15)
q <- ggplot(Metacluster_immune2, aes(x = variable, y = value)) +
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
        strip.text = element_text(size = 23, face = "bold")) + 
  stat_compare_means(aes(group = Timepoint), method = "anova", label = "p.signif", label.y = 6, size = 6, face = "bold") +  
  scale_y_continuous(trans='log10', limits = c(0.1,50000000), labels = trans_format("log10", math_format(10^.x)))+ 
  theme_mypub() + my_theme_adapted
plot(q)
  dev.off()
  

 
######################################  
### Tumour Cluster
  
  Metacluster_tumour <- Metacluster_df[, colnames(Metacluster_df) %in% tumour_cluster]
  Metacluster_tumour$Responsiveness <- Metacluster_df$Responsiveness
  Metacluster_tumour$Timepoint <- Metacluster_df$Timepoint
  Metacluster_tumour2 <- melt(Metacluster_tumour)
  cols <- c("A" = '#e63e00', "B" = "darkseagreen3", "C" = 'steelblue')
 
  pdf(paste0(fig_dir2, "/Tumour_clusters.pdf"), height = 3, width = 16)
  p <- ggplot(Metacluster_tumour2, aes(x = variable, y = value)) +
    geom_violin(
      aes(color = Timepoint), trim = FALSE,
      position = position_dodge(0.9) 
    ) +
    geom_boxplot(
      aes(color = Timepoint), width = 0.15,
      position = position_dodge(0.9)
    ) +
    scale_color_manual(values = cols) + 
    xlab("Cluster") +
    ylab(expression("cell count/"~mm^2)) +
    facet_grid(~Responsiveness, scales = "free_x", space = "free_x", labeller = as_labeller(c(`NE` = "NRP",`RE` = "REP"))) +
    theme(panel.spacing = unit(0.5, "cm")) + 
    theme(strip.text = element_text(size = 16, face = "bold", colour = "orange")) +
    stat_compare_means(aes(group = Timepoint), method = "anova", label = "p.signif", label.y = 6, size = 5, face = "bold") +  
    scale_y_continuous(trans='log10', limits = c(0.1,5000000), labels = trans_format("log10", math_format(10^.x)))+ 
    theme_mypub() + my_theme_adapted
  plot(p)
  dev.off()
  

  
  ######################################  
  ### Other cell cluster
  
  Metacluster_ECM <- Metacluster_df[, colnames(Metacluster_df) %in% Other_cell_cluster]
  Metacluster_ECM$Responsiveness <- Metacluster_df$Responsiveness
  Metacluster_ECM$Timepoint <- Metacluster_df$Timepoint
  Metacluster_ECM2 <- melt(Metacluster_ECM)
  cols <- c("A" = '#e63e00', "B" = "darkseagreen3", "C" = 'steelblue')

  pdf(paste0(fig_dir2, "/ECM_clusters.pdf"), height = 3, width = 12)
  ggplot(Metacluster_ECM2, aes(x = variable, y = value)) +
    geom_violin(
      aes(color = Timepoint), trim = FALSE,
      position = position_dodge(0.9) 
    ) +
    geom_boxplot(
      aes(color = Timepoint), width = 0.15,
      position = position_dodge(0.9)
    ) +
    scale_color_manual(values = cols) + 
    xlab("") +
    ylab(expression("cell count/"~mm^2)) +
    facet_grid(~Responsiveness, scales = "free_x", space = "free_x", labeller = as_labeller(c(`NE` = "NRP",`RE` = "REP"))) +
    theme(panel.spacing = unit(0.5, "cm"),
          strip.text = element_text(size = 20, face = "bold")) + 
    stat_compare_means(aes(group = Timepoint), method = "anova", label = "p.signif", label.y = 6, size = 5, face = "bold") +  
    scale_y_continuous(trans='log10', limits = c(0.1,5000000), labels = trans_format("log10", math_format(10^.x)))+ 
    theme_mypub() + my_theme_adapted
  plot(p) +
    theme(panel.spacing = unit(0.5, "cm"),
          strip.text = element_text(size = 20, face = "bold"))
  dev.off()
  
  
  
