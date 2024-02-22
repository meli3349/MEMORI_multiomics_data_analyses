

### Normalisation of raw count gene expression in samples >1.5M reads


dataset_dir = "~/analysis/RNA_Seq/input_files/GeneExpression_Mastertables"

# prepare output dirs:
dir.create(dataset_dir, showWarnings=FALSE, recursive=TRUE)

# Read in raw count Mastertable
GeneIndex_MasterTable_raw_counts <- read.table(paste0(dataset_dir, "/GeneIndex_MasterTable_raw_gene_counts_batch_1_2_3.txt"), header = T, stringsAsFactors = F)

# Filter samples with > 1.5 million reads 
morethan1.5 <- colnames(GeneIndex_MasterTable_raw_counts[,c(10:ncol(GeneIndex_MasterTable_raw_counts))])[which(colSums(GeneIndex_MasterTable_raw_counts[,c(10:ncol(GeneIndex_MasterTable_raw_counts))])>=1.5e6)]
GeneIndex_MasterTable_raw_counts_1.5M <- GeneIndex_MasterTable_raw_counts[,c("GeneID", "Name",morethan1.5)]


# Read in clin.df and drop samples < 1.5 Mio reads
clin.df <- read.table("~/analysis/General_files/MEMORI_specific/Clinical_df_RNA_batch_1_2_3.txt", header = T, stringsAsFactors = T)
rownames(clin.df) <- clin.df[,1]
clin.df_1.5M <- clin.df[morethan1.5, ]
clin.df_1.5M$batch[clin.df_1.5M$Batch == 1] <- "a"
clin.df_1.5M$batch[clin.df_1.5M$Batch == 2] <- "b"
clin.df_1.5M$batch[clin.df_1.5M$Batch == 3] <- "c"

all(rownames(clin.df_1.5M) %in% colnames(GeneIndex_MasterTable_raw_counts_1.5M))

### Build DeSeq data frame
library(edgeR)
rownames(GeneIndex_MasterTable_raw_counts_1.5M) <- GeneIndex_MasterTable_raw_counts_1.5M$GeneID
y <- DGEList(counts = GeneIndex_MasterTable_raw_counts_1.5M[,c(3:ncol(GeneIndex_MasterTable_raw_counts_1.5M))])
keep <- filterByExpr(y)
y_filtered <- y[keep, keep.lib.sizes=FALSE]


# Transform genecounts into counts per million (cpm) and filter for expressed genes
logcpm_filtered <- cpm(y_filtered, log=TRUE)
Gene_ENSEMBLE <- rownames(logcpm_filtered)
logcpm_filtered <- as.data.frame(cbind(Gene_ENSEMBLE,logcpm_filtered))


# Transform genecounts into counts per million (cpm)
logcpm_nonfiltered <- cpm(y, log=TRUE)
Gene_ENSEMBLE <- rownames(logcpm_nonfiltered)
logcpm_nonfiltered <- as.data.frame(cbind(Gene_ENSEMBLE,logcpm_nonfiltered))

write.table(logcpm_filtered, file=paste0(dataset_dir, "/GeneIndex_MasterTable1.5M_norm_batch_1_2_3_filtered_expressed_genes.txt"), sep = "\t", row.names = F)
write.table(logcpm_nonfiltered, file=paste0(dataset_dir, "/GeneIndex_MasterTable1.5M_norm_batch_1_2_3.txt_all_genes.txt"), sep = "\t", row.names = F)
