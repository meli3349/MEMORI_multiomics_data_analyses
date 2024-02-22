##############################################################################################################
### Immune escape through PDL1 overexpression 
### This scripts... 
# 1. plots PDL1 expression (RNA-Seq data) during treatment 
# 2. classifies samples according to their PDL1 expression into "immune-escaped through PDL1 overexpression" and "non-immune-escaped through PDL1 overexpression"
# 3. creates a Mastertable with different immune escape mechanisms (HLA LOH, mutations in HLA or B2M, PDL1 overexpression)

require(ggplot2); 
require(ggpubr); 
require(scales); 
require(reshape2); 
library(edgeR)

# prepare output dirs:
fig_dir = "~/analysis/multi_omic/immune_escape/plots/PDL1_expression"
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)
data_dir = "~/analysis/multi_omic/immune_escape/tables"
dir.create(data_dir, showWarnings=FALSE, recursive=TRUE)

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


# Read in non-normalized, non filtered GeneIndex_MasterTable
GeneIndex_MasterTable_raw_counts <- read.table("~/analysis/RNA_Seq/input_files/GeneExpression_Mastertables/merged_non_normalized_Mastertable/GeneIndex_MasterTable_raw_gene_counts_batch_1_2_3.txt", header = T, stringsAsFactors = F)

# Filter for samples with > 1.5 million readcounts
morethan1.5 <- colnames(GeneIndex_MasterTable_raw_counts[,c(10:ncol(GeneIndex_MasterTable_raw_counts))])[which(colSums(GeneIndex_MasterTable_raw_counts[,c(10:ncol(GeneIndex_MasterTable_raw_counts))])>=1.5e6)]
GeneIndex_MasterTable_raw_counts_1.5M <- GeneIndex_MasterTable_raw_counts[,c("GeneID", "Name",morethan1.5)]

### Normalization with edgeR
rownames(GeneIndex_MasterTable_raw_counts_1.5M) <- GeneIndex_MasterTable_raw_counts_1.5M$GeneID
y <- DGEList(counts = GeneIndex_MasterTable_raw_counts_1.5M[,c(3:ncol(GeneIndex_MasterTable_raw_counts_1.5M))])
logcpm_nonfiltered <- cpm(y, log=TRUE)
Gene_ENSEMBLE <- rownames(logcpm_nonfiltered)

# Calculate PDL1 expression in logcpm
logcpm_nonfiltered <- as.data.frame(cbind(Gene_ENSEMBLE,logcpm_nonfiltered))

# Calculate PDL1 expression in cpm
cpm_nonfiltered <- cpm(y, log=F)
Gene_ENSEMBLE <- rownames(cpm_nonfiltered)
cpm_nonfiltered <- as.data.frame(cbind(Gene_ENSEMBLE,cpm_nonfiltered))

# Get PDL1 gene expression (= ENSG00000120217)
# for logcpm
PDL1_df_logcpm <- logcpm_nonfiltered[which(logcpm_nonfiltered$Gene_ENSEMBLE == "ENSG00000120217"), ]
PDL1_df_logcpm <- as.data.frame(t(PDL1_df_logcpm)); 
PDL1_df_logcpm <- as.data.frame(cbind(rownames(PDL1_df_logcpm), PDL1_df_logcpm$ENSG00000120217))
colnames(PDL1_df_logcpm ) <- c("Sample", "PDL1_logcpm"); PDL1_df_logcpm  <- PDL1_df_logcpm[-1,]

# for cpm
PDL1_df_cpm <- cpm_nonfiltered[which(cpm_nonfiltered$Gene_ENSEMBLE == "ENSG00000120217"), ]
PDL1_df_cpm <- as.data.frame(t(PDL1_df_cpm))
PDL1_df_cpm <-  as.data.frame(cbind(rownames(PDL1_df_cpm), PDL1_df_cpm$ENSG00000120217))
colnames(PDL1_df_cpm ) <- c("Sample", "PDL1_cpm"); PDL1_df_cpm  <- PDL1_df_cpm[-1,]

# logcpm and cpm
PDL1_df <-  merge(PDL1_df_logcpm, PDL1_df_cpm, by = "Sample")
PDL1_df$Timepoint <- sapply(PDL1_df$Sample, function(x) substr(x, 13,13))
PDL1_df$PDL1_logcpm <- as.numeric(PDL1_df$PDL1_logcpm); PDL1_df$PDL1_cpm <- as.numeric(PDL1_df$PDL1_cpm)

# 1a. Plot PDL1 expression for TP A-C in cpm 
pdf(paste0(fig_dir, "/PDL1_cpm_Expression.pdf"), height = 5, width = 5)
p <- ggplot2.violinplot(data=PDL1_df, xName="Timepoint", yName="PDL1_cpm",
                        groupName="Timepoint",
                        groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue')), .8), showLegend=FALSE,
                        backgroundColor="white", xtitle="", ytitle="PDL1 expression [cpm]", 
                        mainTitle="PDL1 expression during treatment",
                        addDot=T, dotSize=0.5,
                        removePanelGrid=F,removePanelBorder=F,
                        axisLine=c(0.5, "solid", "black")) + 
  stat_compare_means(comparisons = list(c('A','B'),c('A','C')), method = "wilcox.test", size = 5, face="bold") +
  scale_y_continuous(trans='log10', limits=c(0.1,200)) # + ylim(c(0, 260))
q <- p + my_theme + theme(plot.title = element_text(color = "black", size = 16, face="bold"))
plot(q)
dev.off()


# 1b. Plot PDL1 expression for TP A-C in logcpm 
pdf(paste0(fig_dir, "/PDL1_Expression.pdf"), height = 5, width = 5)
p <- ggplot2.violinplot(data=PDL1_df, xName="Timepoint", yName="PDL1_logcpm",
                        groupName="Timepoint",
                        groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue')), .8), showLegend=FALSE,
                        backgroundColor="white", xtitle="", ytitle="PDL1 expression [logcpm]", 
                        mainTitle="PDL1 expression during treatment",
                        addDot=T, dotSize=0.5,
                        removePanelGrid=F,removePanelBorder=F,
                        axisLine=c(0.5, "solid", "black")) + 
  stat_compare_means(comparisons = list(c('A','B'),c('A','C')), method = "wilcox.test", size = 5, face="bold") +
  scale_y_continuous(breaks=c(0,2,4,6), limits=c(-.5,7))
q <- p + my_theme + theme(plot.title = element_text(color = "black", size = 16, face="bold"))
plot(q)
dev.off()


###########################################################################################################################
# 2. Classification of samples as immune escaped through PDL1 overexpression 
# PDL1 immune escape is defined as overexpression of PDL1 > 2 standard deviations from PDL1 mean expression at timepoint A
df_A <- PDL1_df[which(PDL1_df$Timepoint == "A"), ]
PDL1_escape_limit_cpm <- mean(df_A$PDL1_cpm) + 2*sd(df_A$PDL1_cpm)
PDL1_df$PDL1_immune_escape <- ifelse(PDL1_df$PDL1_cpm > PDL1_escape_limit_cpm, "PDL1_high", "FALSE")

###########################################################################################################################
# 3. Adding PDL1 overexpression as immune escape mechanism to immune escape Mastertable

# Read in table with genetic immune escape data 
escape_df <- read.table("~/analysis/multi_omic/input_files/MEMORI_genetic_escape_master_file.txt", stringsAsFactors = F, header = T, sep = "\t")
escape_df <- as.data.frame(escape_df)

# Read in file with overvie of samples with matching RNA-Seq and WES data 
matching_overview <- read.table("~/analysis/General_files/MEMORI_specific/Matching_RNA_DNA_samples.txt", header = T, stringsAsFactors = F, sep = "\t")
matching_overview <- matching_overview[, c("RNA_file", "WES_file")]

PDL1_to_merge <- merge(matching_overview, PDL1_df, by.x = "RNA_file", by.y = "Sample", all = T)
PDL1_to_merge <- as.data.frame(PDL1_to_merge[,c("WES_file", "RNA_file", "PDL1_immune_escape")])

final_immune_escape_df <- merge(PDL1_to_merge, escape_df, by.x = "WES_file", by.y = "Sample",  all = T)
final_immune_escape_df$Patient <- ifelse(is.na(final_immune_escape_df$Patient), sapply(final_immune_escape_df$RNA_file, function(x) substr(x, 7,12)), sapply(final_immune_escape_df$Patient, function(x) substr(x, 8,13)))
final_immune_escape_df$PDL1_immune_escape <- ifelse(final_immune_escape_df$PDL1_immune_escape == "y", "PDL1_high", ifelse(final_immune_escape_df$PDL1_immune_escape == "n", "FALSE", final_immune_escape_df$PDL1_immune_escape))
final_immune_escape_df$Sample <- ifelse(is.na(final_immune_escape_df$WES_file), paste0("MEMORI_", sapply(final_immune_escape_df$RNA_file, function(x) substr(x, 7,22))), final_immune_escape_df$WES_file)

# Classification of final HLA LOH result if HLA LOH detected either in Sequenza (-> HLAregionCNLOH == "LOH (?)") or with LOHHLA software (-> HLAloh == "LOH" | "LOH (?)")
final_immune_escape_df$final_HLA_LOH <-  ifelse(final_immune_escape_df$HLAloh == "LOH" | final_immune_escape_df$HLAloh == "LOH (?)" |
                                     final_immune_escape_df$HLAregionCNLOH == "LOH (?)", "LOH", "FALSE")
final_immune_escape_df$genetic_escape <- ifelse(final_immune_escape_df$HLAmut ==  "mutated" | final_immune_escape_df$HLAloh == "LOH" | final_immune_escape_df$HLAloh == "LOH (?)" |
                                                  final_immune_escape_df$B2Mmut == "mutated" | final_immune_escape_df$HLAregionCNLOH == "LOH (?)", "y", "n")
final_immune_escape_df$final_escape_result <- ifelse(final_immune_escape_df$PDL1_immune_escape == "PDL1_high" | final_immune_escape_df$HLAmut ==  "mutated" | final_immune_escape_df$HLAloh == "LOH" | final_immune_escape_df$HLAloh == "LOH (?)" |
                                          final_immune_escape_df$B2Mmut == "mutated" | final_immune_escape_df$HLAregionCNLOH == "LOH (?)", "y", "n")
# Export table
final_immune_escape_df <- final_immune_escape_df[c("Sample", "Patient", "WES_file", "RNA_file", "Patient",  "PDL1_immune_escape", "final_HLA_LOH", "HLAloh", "HLAregionCNLOH", "HLAmut", "B2Mmut", "genetic_escape", "final_escape_result")]
write.table(final_immune_escape_df, file = paste0(data_dir, "/MEMORI_escape_master_file_incl_PDL1.txt"), sep = "\t", quote = F, row.names = F)

