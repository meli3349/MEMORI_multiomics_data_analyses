

## This script plots Cibersort results (absolute values) as bar graphs
library(ggplot2); library(reshape2); library(tidyselect); library(dplyr)

# prepare output dirs:
fig_dir = "~/analysis/RNA_Seq/immune_cell_deconvolution/Cibersort/plots"
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)

# themes for plots
theme_mypub_grid <- function(base_size = 14,
                             base_family = ""){
  theme_bw(base_size = base_size, 
           base_family = base_family) %+replace%
    theme(
      #panel.grid.major = element_blank(),   
      panel.grid.minor = element_blank(),   
      
      complete = TRUE
    )
}

my_theme <- theme(
  plot.title = element_text(color="black", size=16, face="bold"),
  axis.title.x = element_text(color="black", size=16, face="bold"),
  axis.text.x = element_text(color="black", size=16, face="bold"),
  axis.title.y = element_text(color="black", size=16, face="bold"),
  axis.text.y = element_text(color="black", size=16, face="bold"),
  legend.title = element_text(color = "black", size = 14, face="bold"),
  legend.text = element_text(color = "black", size = 14, face="bold"),
  
)

# Read in Cibersort table with absolute immune cell values generated from https://cibersortx.stanford.edu
CibersortTable <- read.table("~/analysis/RNA_Seq/input_files/Cibersort_tables/CIBERSORT.Output_Abs_Job6.txt", header = T, sep = "\t")

# Sum up absolute Cibersort values from different cell subtypes to main cell type (e.g. Macrophages M0, Macrophages M1 and Macrophages M2 are combined to "Total Macrophages")
CibersortTable$TotalMacrophages <- rowSums(CibersortTable %>% dplyr::select(starts_with("Macrophages"))) 
CibersortTable$TotalNKcells <- rowSums(CibersortTable %>% dplyr::select(starts_with("NK.cells")))
CibersortTable$TotalMonocytes <- rowSums(CibersortTable %>% dplyr::select(starts_with("Monocytes")))
CibersortTable$TotalDendricCells <- rowSums(CibersortTable %>% dplyr::select(starts_with("Dendric.cells")))
CibersortTable$TotalMastCells <- rowSums(CibersortTable %>% dplyr::select(starts_with("Mast.cells")))
CibersortTable$TotalEosinophils <- rowSums(CibersortTable %>% dplyr::select(starts_with("Eosinophils")))
CibersortTable$TotalNeutrophils <- rowSums(CibersortTable %>% dplyr::select(starts_with("Neutrophils")))
CibersortTable$OtherBcells <- rowSums(CibersortTable %>% dplyr::select(starts_with("B.cells") | starts_with("Plasma.cells")))
CibersortTable$OtherTcells <- rowSums(CibersortTable[,c(9:11)])
CibersortTable$TotalCD4cells <- rowSums(CibersortTable %>% dplyr::select(starts_with("T.cells.CD4")))
CibersortTable$TotalCD8cells <- rowSums(CibersortTable %>% dplyr::select(starts_with("T.cells.CD8")))


# Analyses of main cell types (combined scores from same cell types) 
CibersortTablePlot <- CibersortTable[,c(1,27:37)]
colnames(CibersortTablePlot) <- c("sampleID", "Macrophages", "NK cells", "Monocytes", "Dendric cells", "Mast cells", "Eosinphils", "Neutrophils", "Other B-cells", "Other T-cells", "CD4 T-cells", "CD8 T-cells")
NewCibersortTablePlot <- melt(CibersortTablePlot)


# Create bargraphs for each sample 
NewCibersortTablePlot_1 <- NewCibersortTablePlot 
NewCibersortTablePlot_1$sampleID <- str_replace_all(NewCibersortTablePlot_1$sampleID, "NE", "NRP")
NewCibersortTablePlot_1$sampleID <- str_replace_all(NewCibersortTablePlot_1$sampleID, "RE", "REP")

pdf(paste0(fig_dir, "/Cibersort_absolute_score_1.5M_samples.pdf"), width = 18, height = 10)
p <- ggplot(NewCibersortTablePlot_1, aes(sampleID, value, fill= variable)) + geom_bar(stat= "Identity") +  
  theme(axis.title.y = element_text(size = 14, face = "bold"))+ theme(axis.text.y = element_text(size = 14, face = "bold", color ="black")) + theme(legend.text=element_text(size= 14, face = "bold"))+ theme(legend.title = element_text(size= 16, face = "bold")) + 
  scale_fill_manual(values = c("firebrick2", "orange","grey", "chocolate4", "pink","darkmagenta","khaki2", "lightgreen", "lightblue", "steelblue", "darkblue")) +
  ylab("Cibersort absolute score") + xlab("") +
  scale_x_discrete(label=sapply(NewCibersortTablePlot$sampleID, function(x) substr(x, 7,17))) + labs(fill = "Immune cell types") + theme_mypub_grid()  + my_theme + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 14, face = "bold", color = "black"))
plot(p)
dev.off()


# Create Bargraph for REP and NRP per timepoint for summmarized immune cells
NewCibersortTablePlot_2 <- NewCibersortTablePlot
NewCibersortTablePlot_2$Responsiveness <- sapply(NewCibersortTablePlot_2$sampleID, function(x) substr(x,7, 8))
NewCibersortTablePlot_2$Timepoint <- sapply(NewCibersortTablePlot_2$sampleID, function(x) substr(x,13, 13))

NewCibersortTablePlot_2 <- aggregate(NewCibersortTablePlot_2$value, by = list(NewCibersortTablePlot_2$variable, NewCibersortTablePlot_2$Responsiveness, NewCibersortTablePlot_2$Timepoint), FUN = mean)
colnames(NewCibersortTablePlot_2) <- c("variable", "Responsiveness", "Timepoint", "value")

pdf(paste0(fig_dir, "/Cibersort_absolute_score_1.5M_per_group.pdf"), width = 6.5, height = 5) 
q <- ggplot(NewCibersortTablePlot_2, aes(Timepoint, value, fill= variable)) + geom_bar(stat= "Identity") + theme(axis.text.x = element_text(size = 16, face = "bold", color = "black")) + 
  theme(axis.title.y = element_text(size = 16, face = "bold"))+ theme(axis.text.y = element_text(size = 16, face = "bold", color ="black")) + theme(legend.text=element_text(size= 14, face = "bold"))+ theme(legend.title = element_text(size= 16, face = "bold")) + 
  scale_fill_manual(values = c("firebrick2", "orange","grey", "chocolate4", "pink","darkmagenta","khaki2", "lightgreen", "lightblue", "steelblue", "darkblue"))+
  theme_mypub_grid() + my_theme + 
  facet_grid(~Responsiveness, scales = "free_x", space = "free_x", labeller = as_labeller(c("RE" = "REP", "NE" = "NRP"))) +
  theme(panel.spacing = unit(0.5, "cm"),
        strip.text = element_text(size = 14, face = "bold")) + 
  ylab("Cibersort absolute score") + xlab("") + labs(fill = "Immune cell types")
plot(q)
dev.off()


# Create Bargraph for REP and NRP per timepoint for all Cibersort immune cells
Complete_CibersortTablePlot <- CibersortTable[,c(1:23)]
Complete_CibersortTablePlot$Mast_cells <- Complete_CibersortTablePlot$Mast.cells.resting + Complete_CibersortTablePlot$Mast.cells.activated
Complete_CibersortTablePlot$Others <- Complete_CibersortTablePlot$T.cells.follicular.helper + Complete_CibersortTablePlot$T.cells.gamma.delta
Complete_CibersortTablePlot <- Complete_CibersortTablePlot %>% select(-T.cells.follicular.helper, -T.cells.gamma.delta, -Mast.cells.resting, -Mast.cells.activated)
colnames(Complete_CibersortTablePlot) <- c("sampleID", "B-cells naive", "B-cells memory", "Plasma cells", "CD8 T-cells", "CD4 T-cells naive", "CD4 T-cells rest", "CD4 T-cells active", "T-cell reg", 
                                           "NK cells rest", "NK cells active", "Monocytes", "Macrophages M0", "Macrophages M1", "Macrophages M2", "Dendritic cells rest", "Dendritic cells active",   
                                            "Eosinophils", "Neutrophils", "Mast cells", "Others")

NewCibersortTablePlot_3 <- melt(Complete_CibersortTablePlot)
NewCibersortTablePlot_3$Responsiveness <- sapply(NewCibersortTablePlot_3$sampleID, function(x) substr(x,7, 8))
NewCibersortTablePlot_3$Timepoint <- sapply(NewCibersortTablePlot_3$sampleID, function(x) substr(x,13, 13))

NewCibersortTablePlot_3 <- aggregate(NewCibersortTablePlot_3$value, by = list(NewCibersortTablePlot_3$variable, NewCibersortTablePlot_3$Responsiveness, NewCibersortTablePlot_3$Timepoint), FUN = mean)
colnames(NewCibersortTablePlot_3) <- c("variable", "Responsiveness", "Timepoint", "value")
NewCibersortTablePlot_3$variable <- factor(NewCibersortTablePlot_3$variable, levels = c("Others", "Macrophages M2", "Macrophages M1", "Macrophages M0", "NK cells active", "NK cells rest", "Monocytes", "Dendritic cells active", "Dendritic cells rest", "Mast cells", "Eosinophils", "Neutrophils", "Plasma cells", "B-cells memory", "B-cells naive", "T-cell reg", "CD4 T-cells active", "CD4 T-cells rest", "CD4 T-cells naive", "CD8 T-cells"))

pdf(paste0(fig_dir, "/Cibersort_absolute_score_1.5M_per_group_extended_imme_cells.pdf"), width = 6.5, height = 6) 
q <- ggplot(NewCibersortTablePlot_3, aes(Timepoint, value, fill= variable)) + geom_bar(stat= "Identity") + theme(axis.text.x = element_text(size = 16, face = "bold", color = "black")) + 
  theme(axis.title.y = element_text(size = 16, face = "bold"))+ theme(axis.text.y = element_text(size = 16, face = "bold", color ="black")) + theme(legend.text=element_text(size= 12, face = "bold"))+ theme(legend.title = element_text(size= 13, face = "bold")) + 
  scale_fill_manual(values = c("grey33", "lightsalmon", "coral2", "firebrick", "orange" , "orange2", "grey", "goldenrod3","chocolate4", "pink", "darkmagenta", "khaki2", "lightgreen", "aquamarine3", "forestgreen", "floralwhite", "slateblue1", "lightskyblue3", "steelblue", "darkblue"))+
  theme_mypub_grid() + my_theme + 
  facet_grid(~Responsiveness, scales = "free_x", space = "free_x", labeller = as_labeller(c("RE" = "REP", "NE" = "NRP"))) +
  theme(panel.spacing = unit(0.5, "cm"),
        strip.text = element_text(size = 16, face = "bold")) +
  ylab("Cibersort absoulte score") + xlab("") + labs(fill = "Immune cell types")
plot(q)
dev.off()


