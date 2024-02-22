###############################################################################
## KEGG and GO analyses for samples from different timepoints and treatments ##
###############################################################################
library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# prepare output dirs:
data_dirGO = "~/analysis/RNA_Seq/pathway_analyses/GO_analysis/Tables"
dir.create(data_dirGO, showWarnings=FALSE, recursive=TRUE)

data_dirKEGG = "~/analysis/RNA_Seq/pathway_analyses/KEGG_analysis/Tables"
dir.create(data_dirKEGG, showWarnings=FALSE, recursive=TRUE)

# Read in Mastertable with raw gene counts 
MasterTable <- read.table("~/analysis/RNA_Seq/input_files/GeneExpression_Mastertables/GeneIndex_MasterTable_raw_gene_counts_batch_1_2_3.txt", header = T, sep = "\t")

# Droping samples with < 1.5Mio gene counts
morethan1.5 <- colnames(MasterTable[,c(10:ncol(MasterTable))])[which(colSums(MasterTable[,c(10:ncol(MasterTable))])>=1.5e6)]
MasterTable_1.5M <- MasterTable[,c("GeneID",morethan1.5)]
MasterTable_1.5M <- as.data.frame(MasterTable_1.5M)

# Build a df with clinical information on responsiveness, timepoint, patient, treatment for subsequent differential expression analyses between clinical subgroups
clin.df.MergedMasterTable1.5M = data.frame(patient=names(MergedMasterTable1.5M))
clin.df.MergedMasterTable1.5M$Responsiveness = sapply(clin.df.MergedMasterTable1.5M$patient, function(x) substr(x,7,8))
clin.df.MergedMasterTable1.5M$PatientNumber = sapply(clin.df.MergedMasterTable1.5M$patient, function(x) substr(x,7,12))
clin.df.MergedMasterTable1.5M$Timepoint = sapply(clin.df.MergedMasterTable1.5M$patient, function(x) substr(x,13, 13))
clin.df.MergedMasterTable1.5M$Responsiveness_Timepoint <- paste0(clin.df.MergedMasterTable1.5M$Responsiveness, "_", clin.df.MergedMasterTable1.5M$Timepoint)

# add column with timepoint C split up for Non-Responder and Responder, in order to later split groups by treatment (chemo versus radiochemo)
clin.df.MergedMasterTable1.5M$Treatment_broken_down <-  ifelse(clin.df.MergedMasterTable1.5M$Timepoint == "A", "A", 
                                                                                      ifelse(clin.df.MergedMasterTable1.5M$Timepoint == "B", "B",
                                                                                             ifelse ((clin.df.MergedMasterTable1.5M$Timepoint == "C" & clin.df.MergedMasterTable1.5M$Responsiveness == "RE"), "C_RE",
                                                                                                     ifelse ((clin.df.MergedMasterTable1.5M$Timepoint == "C" & clin.df.MergedMasterTable1.5M$Responsiveness == "NE"), "C_NE", "XX"))))
# Normalization with DESeq
dds <- DESeqDataSetFromMatrix(countData = MergedMasterTable1.5M,
                              colData = clin.df.MergedMasterTable1.5M,
                              design = ~ Responsiveness_Timepoint)

dds_TP <- DESeqDataSetFromMatrix(countData = MergedMasterTable1.5M,
                              colData = clin.df.MergedMasterTable1.5M,
                              design = ~ Timepoint)

dds_by_treatment <- DESeqDataSetFromMatrix(countData = MergedMasterTable1.5M,
                              colData = clin.df.MergedMasterTable1.5M,
                              design = ~ Treatment_broken_down)


dds <- DESeq(dds)
dds_TP <- DESeq(dds_TP)
dds_by_treatment <- DESeq(dds_by_treatment)



#### Calculate differential expression between different clinical subgroups
# subgroups from different timepoints
resA_B <- results(dds_TP, contrast= c("Timepoint", "A", "B"))
resB_C <- results(dds_TP, contrast= c("Timepoint", "B", "C"))
resA_C <- results(dds_TP, contrast= c("Timepoint", "A", "C"))

# subgroups broken down by treatments
res_A_ReC  <- results(dds_by_treatment, contrast= c("Treatment_broken_down", "A", "C_RE"))
res_B_ReC  <- results(dds_by_treatment, contrast= c("Treatment_broken_down", "B", "C_RE"))
res_A_NeC  <- results(dds_by_treatment, contrast= c("Treatment_broken_down", "A", "C_NE"))
res_B_NeC  <- results(dds_by_treatment, contrast= c("Treatment_broken_down", "B", "C_NE"))

# subgroups with different responsiveness and from different timepoints
res_ReA_ReB  <- results(dds, contrast= c("Responsiveness_Timepoint", "RE_A", "RE_B"))
res_ReB_ReC  <- results(dds, contrast= c("Responsiveness_Timepoint", "RE_B", "RE_C"))
res_ReA_ReC  <- results(dds, contrast= c("Responsiveness_Timepoint", "RE_A", "RE_C"))

res_ReA_NeA  <- results(dds, contrast= c("Responsiveness_Timepoint", "RE_A", "NE_A"))
res_ReB_NeB  <- results(dds, contrast= c("Responsiveness_Timepoint", "RE_B", "NE_B"))
res_ReC_NeC  <- results(dds, contrast= c("Responsiveness_Timepoint", "RE_C", "NE_C"))

res_NeA_NeB  <- results(dds, contrast= c("Responsiveness_Timepoint", "NE_A", "NE_B"))
res_NeB_NeC  <- results(dds, contrast= c("Responsiveness_Timepoint", "NE_B", "NE_C"))
res_NeA_NeC  <- results(dds, contrast= c("Responsiveness_Timepoint", "NE_A", "NE_C"))

# check how many genes are significantly differentially expressed
sum(resA_C$padj < 0.05, na.rm = TRUE)
sum(resB_C$padj < 0.05, na.rm = TRUE)
sum(resA_B$padj < 0.05, na.rm = TRUE)

sum(res_A_ReC$padj < 0.05, na.rm = TRUE)
sum(res_B_ReC$padj < 0.05, na.rm = TRUE)
sum(res_A_NeC$padj < 0.05, na.rm = TRUE)
sum(res_B_NeC$padj < 0.05, na.rm = TRUE)

sum(res_ReA_ReB$padj < 0.05, na.rm = TRUE)
sum(res_ReB_ReC$padj < 0.05, na.rm = TRUE)
sum(res_ReA_ReC$padj < 0.05, na.rm = TRUE)

sum(res_ReA_NeA$padj < 0.05, na.rm = TRUE)
sum(res_ReB_NeB$padj < 0.05, na.rm = TRUE)
sum(res_ReC_NeC$padj < 0.05, na.rm = TRUE)

sum(res_NeA_NeB$padj < 0.05, na.rm = TRUE)
sum(res_NeB_NeC$padj < 0.05, na.rm = TRUE)
sum(res_NeA_NeC$padj < 0.05, na.rm = TRUE)

#######################################################################################
### Run Gene Set Enrichment Analysis with ClusterProfiler for GO and KEGG pathways ### 
#######################################################################################

organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

# a.) GO enrichment analysis
################################

comparisons <- c("resA_B", "resB_C", "resA_C", "res_A_ReC", "res_A_NeC", "res_B_ReC", "res_B_NeC", "res_ReA_ReB", "res_ReB_ReC", "res_ReA_ReC", "res_ReA_NeA", "res_ReB_NeB", "res_ReC_NeC", "res_NeA_NeB", "res_NeB_NeC", "res_NeA_NeC")

for (i in comparisons[13]){
# reading in data from Deseq2
  comparison_name <- i
df = get(i)

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- rownames(df)

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 40, 
             maxGSSize = 500, 
             pvalueCutoff = 0.5, 
             verbose = T, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "fdr")

GO_df <- as.data.frame(gse)
write.table(GO_df, file = paste0(data_dirGO, "/GO_table_", comparison_name, ".p.adj.txt"), quote = F, sep = "\t", row.names = F)



# b.) KEGG Gene Set Enrichment Analysis
#########################################

# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
# remove duplicate IDS
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped
df2 = df[rownames(df) %in% dedup_ids$ENSEMBL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Entrez_ID = dedup_ids$ENTREZID

# Create a vector of the gene universe
kegg_gene_list <- df2$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Entrez_ID

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kegg_organism = "hsa"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 1000,
               minGSSize    = 40,
               maxGSSize    = 500,
               pvalueCutoff = 0.9,
               pAdjustMethod = "fdr",
               keyType       = "ncbi-geneid")


kk2_df <- as.data.frame(kk2)
write.table(kk2_df, file = paste0(data_dirKEGG, "/KEGG_table_", comparison_name, ".p.adj.txt"), quote = F, sep = "\t", row.names = F)

}
