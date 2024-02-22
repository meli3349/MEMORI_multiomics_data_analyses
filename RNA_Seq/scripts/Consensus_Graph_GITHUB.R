
##########################################
# Plot ConsensusTME immune cell scores   # 
##########################################
# Tumour microenvironment cell estimation with ConsensusTME
# https://github.com/cansysbio/ConsensusTME

library(ConsensusTME)
library(biomaRt)
library(reshape2)

# prepare output dirs:
fig_dir = "~/analysis/RNA_Seq/immune_cell_deconvolution/ConsensusTME/Plots"
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)
data_dir = "~/analysis/RNA_Seq/immune_cell_deconvolution/ConsensusTME/Data_frames"
dir.create(data_dir, showWarnings=FALSE, recursive=TRUE)

# Read in non-filtered nomralized gene expression Mastertable 
logcpm_nonfiltered <- read.table("~/analysis/RNA_Seq/input_files/GeneExpression_Mastertables/GeneIndex_MasterTable1.5M_norm_batch_1_2_3.txt_all_genes.txt", header = T, stringsAsFactors = F, sep = "\t")

# Get the full ensembl dataset
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

# assign hgnc IDs to logcpm expression table to run Consensus
geneMap = getBM(attributes=c("ensembl_gene_id","entrezgene_id","hgnc_symbol"), mart = ensembl,filters="ensembl_gene_id",values = logcpm_nonfiltered$Gene_ENSEMBLE)
logcpm_with_hgnc <- merge(geneMap, logcpm_nonfiltered, by.x = "ensembl_gene_id", by.y = "Gene_ENSEMBLE", all.y = T)

# exclude genes without HUGO ID. These genes can not be included in ConsensusTME analysis 
logcpm_with_hgnc <- logcpm_with_hgnc[!(is.na(logcpm_with_hgnc$hgnc_symbol) | logcpm_with_hgnc$hgnc_symbol==""), ]

# Identify the duplicated HUGO IDs
duplicates <-logcpm_with_hgnc[which(duplicated(logcpm_with_hgnc$hgnc_symbol)), 'hgnc_symbol']

# Subset the dataframe for all HUGO IDs which aren't duplicated and remove duplicated ENSEMBLE IDs
`%ni%` <- Negate(`%in%`)
logcpm_with_hgnc <- logcpm_with_hgnc[which(logcpm_with_hgnc$hgnc_symbol %ni% duplicates),]
logcpm_with_hgnc <- logcpm_with_hgnc[which(!duplicated(logcpm_with_hgnc$ensembl_gene_id)),]
logcpm_with_hgnc_unique <- logcpm_with_hgnc[order(logcpm_with_hgnc$hgnc_symbol),];row.names(logcpm_with_hgnc_unique) <- logcpm_with_hgnc_unique$hgnc_symbol

# Make input data frame for ConsensusTME numerical
df <- as.data.frame(logcpm_with_hgnc_unique[,c(4:ncol(logcpm_with_hgnc_unique))], stringsAsFactors=F)
samples <- gsub("\\S+(MEMEX_\\S+Q\\d)\\S+","\\1",colnames(df))
bulkExpMatrix    <- do.call(cbind,lapply(samples,function(sam)as.numeric(df[,sam])))
colnames(bulkExpMatrix) <- samples; rownames(bulkExpMatrix) <- rownames(df)


# Cell Type Estimation with ConsensusTME
Consensus_df  <- ConsensusTME::consensusTMEAnalysis(bulkExpMatrix, cancer = "ESCA", statMethod = "ssgsea")
write.table(Consensus_df, file = paste0(fig_dir, "/Consensus_immune_cell_scores_per_sample.txt"), quote = F, row.names = T, sep = "\t")

# Eliminate cell types with no interest
Consensus_df <- Consensus_df[which(rownames(Consensus_df)!= "Immune_Score" & rownames(Consensus_df)!= "Macrophages_M1" & rownames(Consensus_df)!= "Macrophages_M2" & rownames(Consensus_df)!= "Endothelial" & rownames(Consensus_df)!= "Fibroblasts"), ]

Consensus_df <- t(Consensus_df)
colnames(Consensus_df) <- c("B-cells", "Cytotoxic cells", "Dendritic cells", "Eosinophils", "Macrophages", "Mast cells", "NK cells", "Neutrophils", 
                            "CD4 T-cells", "CD8 T-cells", "T-cells gamma delta", "T-cell reg", "Monocytes", "Plasma cells")


# Prepare data frame for plots
NewConsensus_df <- melt(Consensus_df) 
colnames(NewConsensus_df) <- c("Sample", "variable", "value")
NewConsensus_df$variable <- factor(NewConsensus_df$variable, levels = c("Macrophages", "NK cells", "Monocytes", "Dendritic cells", "Mast cells", "Eosinophils", "Neutrophils", "Plasma cells", 
                                                                        "B-cells", "Cytotoxic cells", "T-cell reg", "CD4 T-cells", "CD8 T-cells", "T-cells gamma delta"))

# Create Bargraph for each sample 
NewConsensus_df$Sample <- str_replace_all(NewConsensus_df$Sample, "NE", "NRP")
NewConsensus_df$Sample <- str_replace_all(NewConsensus_df$Sample, "RE", "REP")

pdf(paste0(fig_dir, "/Consensus_absolute_score_1.5M_samples.pdf"), width = 18, height = 10)
p <- ggplot(NewConsensus_df, aes(Sample, value, fill= variable)) + geom_bar(stat= "Identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 14, face = "bold", color = "black")) + 
  theme(axis.title.y = element_text(size = 14, face = "bold"))+ theme(axis.text.y = element_text(size = 14, face = "bold", color ="black")) + theme(legend.text=element_text(size= 14, face = "bold"))+ theme(legend.title = element_text(size= 16, face = "bold")) + 
  scale_fill_manual(values = c("firebrick2", "orange","grey", "chocolate4", "pink","darkmagenta","khaki2", "forestgreen", "lightgreen", "coral2", "lightblue", "steelblue", "darkblue", "slateblue1")) +
  ylab("Consensus score") + xlab("") + labs(fill = "Immune cell types")
plot(p)
dev.off()



# Create Bargraph for REP and NRP for each timepoint
NewConsensus_df_2 <- melt(Consensus_df)
colnames(NewConsensus_df_2) <- c("Sample", "variable", "value")
NewConsensus_df_2$variable <- factor(NewConsensus_df_2$variable, levels = c("Macrophages", "NK cells", "Monocytes", "Dendritic cells", "Mast cells", "Eosinophils", "Neutrophils", "Plasma cells", 
                                                                            "B-cells", "Cytotoxic cells", "T-cell reg", "CD4 T-cells", "CD8 T-cells", "T-cells gamma delta"))

NewConsensus_df_2$Responsiveness_TP <- paste0(sapply(NewConsensus_df_2$Sample, function(x) substr(x,7, 8)), "_", sapply(NewConsensus_df_2$Sample, function(x) substr(x,13, 13)))
NewConsensus_df_2.5 <- aggregate(NewConsensus_df_2, by = list(NewConsensus_df_2$Responsiveness_TP, NewConsensus_df_2$variable), FUN=mean)
NewConsensus_df_2.5$Responsiveness <- sapply(NewConsensus_df_2.5$Group.1, function(x) substr(x,1, 2))
NewConsensus_df_2.5$Timepoint <- sapply(NewConsensus_df_2.5$Group.1, function(x) substr(x,4, 4))

pdf(paste0(fig_dir, "/Consensus_absolute_score_1.5M_per_responsiveness.pdf"), width = 7, height = 5) 
q <- ggplot(NewConsensus_df_2.5, aes(Timepoint, value, fill= Group.2)) + geom_bar(stat= "Identity") + theme(axis.text.x = element_text(size = 16, face = "bold", color = "black")) + 
  theme(axis.title.y = element_text(size = 16, face = "bold"))+ theme(axis.text.y = element_text(size = 16, face = "bold", color ="black")) + theme(legend.text=element_text(size= 14, face = "bold"))+ theme(legend.title = element_text(size= 16, face = "bold")) + 
  scale_fill_manual(values = c("firebrick2", "orange","grey", "chocolate4", "pink","darkmagenta","khaki2", "forestgreen", "lightgreen", "coral2", "lightblue", "steelblue", "darkblue", "slateblue1"))+
  facet_grid(~Responsiveness, scales = "free_x", space = "free_x", labeller = as_labeller(c("NE" = "NRP", 'RE' = "REP"))) +
  theme(panel.spacing = unit(1, "cm"),
        strip.text = element_text(size = 14, face = "bold")) +
  ylab("Consensus score") + xlab("") + labs(fill = "Immune cell types")
plot(q)
dev.off()

