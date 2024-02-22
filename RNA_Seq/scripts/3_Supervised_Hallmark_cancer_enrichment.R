##########################################################################################
##### Hierarchical clustering based on expression of Hallmark of cancer pathways #########
##########################################################################################

# This scripts creates: 
# a.) unsupervised hierarchical clustering based on the expression of all OAC-relevant hallmark of cancer pathways
# b.) supervised hierarchical clustering based on the expression of significantly differentially expressed Hallmark of cancer pathways between samples from Timepoint A/B and samples from Timepoint C

# prepare output dirs:
fig_dir = "~/analysis/RNA_Seq/pathway_analyses/Hallmark_cancer_pathway/Plots"
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)
##########################################################################################################
# a.) Unsupervised hierarchical clustering based on all OAC-relevant Hallmark of cancer pathways
##########################################################################################################
library(pheatmap)
library(dendextend)
library(easyGgplot2)
library(rstatix)
library(ggpubr)

# Hallmark of cancer pathway expression in study sample
Hallmark_df <- read.table("~/analysis/RNA_Seq/pathway_analyses/Hallmark_cancer_pathway/Data_frames/Hallmark_cancer_scores.txt", sep = "\t", header = T, stringsAsFactors = F)

# Exclusion of pathways, that are not relevant to oesophageal adenocarcinoma
Hallmark_df$HALLMARK_SPERMATOGENESIS <- NULL
Hallmark_df$HALLMARK_MYOGENESIS <- NULL
Hallmark_df$HALLMARK_ANDROGEN_RESPONSE <- NULL

Hallmark_df <- data.frame(t(Hallmark_df))
rownames(Hallmark_df) <- sapply(rownames(Hallmark_df), function(x) substr(x, 10, nchar(x)))

# Hierarchical cluster analysis 
hclust_samples <- hclust(dist(t(Hallmark_df)), method = "complete")
my_sample_cluster <- cutree(hclust_samples, k = 2)

# Annotation of Hallmark of cancer category and short name of pathways
Hallmark_cancer_annotation <- read.delim("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/General_lists/GSEA_cancer_hallmarks/Hallmark_cancer_pathway.txt", header = T, stringsAsFactors = F)
Hallmark_cancer_annotation <- Hallmark_cancer_annotation[,c("Hallmark_Cancer_Pathway", "Pathway_category")]
rownames(Hallmark_cancer_annotation) <- Hallmark_cancer_annotation$Hallmark_Cancer_Pathway
rownames(Hallmark_cancer_annotation) <- sapply(rownames(Hallmark_cancer_annotation), function(x) substr(x, 10, nchar(x)))

# Include row annotation with pathway names and pathway category
annotation_df <- data.frame(Hallmark_cancer_annotation[, c("Pathway_category")]); rownames(annotation_df) <- rownames(Hallmark_cancer_annotation)
colnames(annotation_df) <- "Pathway_category"
annotation_df$Pathway_category <- factor(annotation_df$Pathway_category , level = c("oncogenic", "cellular stress", "immune", "stromal", "other"))

# Include column annotation with information on cluster, Responsiveness, Timepoint and sampleID
sample_cluster_df <- data.frame(Cluster = ifelse(test = my_sample_cluster == 1, "cluster1", ifelse(test = my_sample_cluster == 2, "cluster2", "")))
sample_cluster_df$Responsiveness <- sapply(rownames(sample_cluster_df), function(x) substr(x,7, 8))
sample_cluster_df$Responsiveness <- ifelse(sample_cluster_df$Responsiveness == "RE", "REP", "NRP")
sample_cluster_df$Responsiveness <- factor(sample_cluster_df$Responsiveness, level = c("NRP", "REP"))
sample_cluster_df$Timepoint <- sapply(rownames(sample_cluster_df), function(x) substr(x,13, 13))


new_sam_names <- paste0(sample_cluster_df$Responsiveness, sapply(rownames(sample_cluster_df), function(x) substr(x,11, 16)))
Hallmark_df_plot <- Hallmark_df
colnames(Hallmark_df_plot) <- new_sam_names
rownames(sample_cluster_df) <- new_sam_names

my_colour = list(
  Responsiveness = c(REP = "#26A69A", NRP ="#b9adde"),
  Timepoint = c(A = '#e63e00', B = "darkseagreen3", C ="steelblue"),
  Cluster = c(cluster1 = "lightblue", cluster2 = "khaki2"),
  Pathway_category = c('cellular stress' = "orange", immune = "forestgreen", oncogenic = "firebrick3", other = "grey", stromal = "darkblue")
)


# Plot heatmap using all Hallmark of cancer pathways
pdf(paste0(fig_dir, "/Unsupervised_Hallmark_cancer_heatmap.pdf"), height = 9, width = 15)
pheatmap(Hallmark_df_plot,
         annotation_colors = my_colour,
         annotation_col = sample_cluster_df,
         annotation_row = annotation_df,
         cutree_cols = 2)
dev.off()


#############################################################################################################################################################################################################
# b.) Supervised hierarchical clustering based on significantly differentially expressed Hallmark of cancer pathways between samples from Timepoint A/B and samples from Timepoint C
#############################################################################################################################################################################################################


# Caclulate wilcox test for each pathway between samples from TPA/B and samples from TP C
Wilcox_pathway_df <- data.frame(t(Hallmark_df))
Wilcox_pathway_df$Timepoint <- sapply(rownames(Wilcox_pathway_df), function(x) substr(x,13, 13))
Wilcox_pathway_df <- Wilcox_pathway_df[,c(48, 1:47)]

list <- list()
for (i in names(Wilcox_pathway_df)[2:ncol(Wilcox_pathway_df)]){
  single_pathway_df <- Wilcox_pathway_df[,c("Timepoint",i)]
  pathway <- colnames(single_pathway_df[2])
  single_pathway_df$pre_post_treatment <- ifelse(single_pathway_df$Timepoint == "A" | single_pathway_df$Timepoint == "B", "pre", "post")
  single_pathway_df <-single_pathway_df[,c(3, 2)]
  melt_df <- melt(single_pathway_df)
  
  wilcox <- melt_df %>%
    wilcox_test(value ~ pre_post_treatment)
  res <- data.frame(pathway = wilcox$p)
  res <- data.frame(t(res))
  res$pathway <- pathway
  rownames(res) <- NULL; res <- res[, c(2,1)]
  colnames(res) <- c("Pathway", "p_val_compar_TPAB_versus_C")
  list[[i]] <- res
}

Master_wilcox_df <- rbindlist(list)
Master_wilcox_df$p_adj <- p.adjust(Master_wilcox_df$p_val_compar_TPAB_versus_C, method = "fdr")

# Select significant pathways for heatmap
sign_pathways <- Master_wilcox_df[which(Master_wilcox_df$p_adj <0.05), ]; sign_pathways <- sign_pathways$Pathway

Hallmark_df_selected_sign <- Hallmark_df[which((rownames(Hallmark_df)) %in% sign_pathways), ]

# Hierarchical cluster analysis 
hclust_samples_sign <- hclust(dist(t(Hallmark_df_selected_sign)), method = "complete")
my_sample_cluster_sign <- cutree(hclust_samples_sign, k = 2)

# Include column annotation with information on Responsiveness, Timepoint and sampleID
sample_cluster_df_sig <- data.frame(Cluster = ifelse(test = my_sample_cluster_sign == 1, "cluster1", ifelse(test = my_sample_cluster_sign == 2, "cluster2", "")))
sample_cluster_df_sig$Responsiveness <- sapply(rownames(sample_cluster_df_sig), function(x) substr(x,7, 8))
sample_cluster_df_sig$Responsiveness <- ifelse(sample_cluster_df_sig$Responsiveness == "RE", "REP", "NRP")
sample_cluster_df_sig$Responsiveness <- factor(sample_cluster_df_sig$Responsiveness, level = c("NRP", "REP"))
sample_cluster_df_sig$Timepoint <- sapply(rownames(sample_cluster_df_sig), function(x) substr(x,13, 13))
sample_cluster_df_sig <- sample_cluster_df_sig[,-1]

new_sam_names <- paste0(sample_cluster_df_sig$Responsiveness, sapply(rownames(sample_cluster_df_sig), function(x) substr(x,11, 16)))
colnames(Hallmark_df_selected_sign) <- new_sam_names
rownames(sample_cluster_df_sig) <- new_sam_names

pdf(paste0(fig_dir, "/Supervised_Hallmark_cancer_heatmap_sign_pathways.pdf"), height = 7, width = 15)
pheatmap(Hallmark_df_selected_sign,
         annotation_colors = my_colour,
         annotation_col = sample_cluster_df_sig,
         annotation_row = annotation_df,
         cutree_cols = 2)
dev.off()

