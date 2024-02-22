

#################################################################
### Plotting expression of Intogen drivers during treatment   ### 
#################################################################
# This script plots the RNA expression of selected driver genes during treatment

require(ggplot2); 
require(ggpubr); 
require(scales); 
require(reshape2); 
library(edgeR)

# prepare output dirs:
fig_dir = "~/analysis/RNA_Seq/Intogen_driver/plots"
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)

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
my_theme <- theme(
  plot.title = element_text(color="black", size=16, face="bold"),
  axis.title.x = element_text(color="black", size=16, face="bold"),
  axis.text.x = element_text(color="black", size=14, face="bold"),
  axis.title.y = element_text(color="black", size=16, face="bold"),
  axis.text.y = element_text(color="black", size=14, face="bold"),
  legend.title = element_text(color = "black", size = 14, face="bold"),
  legend.text = element_text(color = "black", size = 14, face="bold"),
  
)

# Read in Mastertable with raw gene counts 
MasterTable <- read.table("~/analysis/RNA_Seq/input_files/GeneExpression_Mastertables/GeneIndex_MasterTable_raw_gene_counts_batch_1_2_3.txt", header = T, sep = "\t")

# Droping samples with < 1.5Mio gene counts
morethan1.5 <- colnames(MasterTable[,c(10:ncol(MasterTable))])[which(colSums(MasterTable[,c(10:ncol(MasterTable))])>=1.5e6)]
MasterTable_1.5M <- MasterTable[,c("GeneID", "Name", morethan1.5)]
MasterTable_1.5M <- as.data.frame(MasterTable_1.5M)


# Normalization with EdgeR
rownames(MasterTable_1.5M) <- MasterTable_1.5M$GeneID
y <- DGEList(counts = MasterTable_1.5M[,c(3:ncol(MasterTable_1.5M))])
cpm_nonfiltered <- cpm(y, log=F)
Gene_ENSEMBLE <- rownames(cpm_nonfiltered)
cpm_nonfiltered <- as.data.frame(cbind(Gene_ENSEMBLE,cpm_nonfiltered))


# Read in Intogen list with driver genes form Esophageal Cancer and gastric cancer 
# https://www.intogen.org/ 
Intogen_EC_GC_list  <- read.table("~/analysis/General_files/external_files/Intogen_drivers.txt", stringsAsFactors = F, header = T)

# ENSEMBLE/ HUGO Gene ID and gene location list
ENSEMBLE_list<- read.table("~/analysis/General_files/external_files/compiledGeneInfo.txt", stringsAsFactors = F, header = T)
chr_filt <- paste0("chr", c(1:23, "X", "Y"))
ENSEMBLE_list <- ENSEMBLE_list[which(ENSEMBLE_list$Chr %in% chr_filt),]

# Get ENSEMBL Gene Names and gene location for Oncogenes/TSG from Intogen list
Intogen_loc <- merge(Intogen_EC_GC_list, ENSEMBLE_list, by.x = "Symbol", by.y = "Name")
Intogen_loc  <- Intogen_loc[ ,c("Symbol", "GeneID", "Chr", "Strand", "Start", "End")]
genes_interest <- as.factor(c("TP53", "KRAS", "CDKN2A", "MYC", "ERBB2", "ERBB3", "GNAS", "TOP2A", "PIK3CA", "AXIN1", "KMT2C", "KMT2D", "SMARCA4", "SMAD4", "APC", "ARID1A"))

# Select drivers of interest
Intogen_loc2 <- Intogen_loc[which(Intogen_loc$Symbol %in% genes_interest),]
genes_interest_ENSG <- Intogen_loc2$GeneID

# RNA expresion of selected Intogen drivers
for (d in genes_interest_ENSG){
   gene_name_df <- Intogen_loc2[which(Intogen_loc2$GeneID == d), ]; gene_name <- gene_name_df[1,1]
   single_gene_cpm <- cpm_nonfiltered[which(cpm_nonfiltered$Gene_ENSEMBLE == d), ]
   single_gene_cpm <- as.data.frame(t(single_gene_cpm))
   single_gene_cpm <-  as.data.frame(cbind(rownames(single_gene_cpm), single_gene_cpm[,1]))
   colnames(single_gene_cpm) <- c("Sample", paste0(gene_name, "_cpm")); single_gene_cpm  <- single_gene_cpm[-1,]

# Add clinical data
single_gene_cpm$Timepoint <- sapply(single_gene_cpm$Sample, function(x) substr(x, 13,13))
single_gene_cpm$Responsiveness <- sapply(single_gene_cpm$Sample, function(x) substr(x, 7,8))
single_gene_cpm[,2] <- as.numeric(single_gene_cpm[,2])

# Plot driver gene expression for timepoint A-C
pdf(paste0(fig_dir, "/", gene_name ,"_cpm_Expression.pdf"), height = 5, width = 5)
p <- ggplot2.violinplot(data=single_gene_cpm, xName="Timepoint", yName=paste0(gene_name, "_cpm"),
                        groupName="Timepoint",
                        groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue')), .8), showLegend=FALSE,
                        backgroundColor="white", xtitle="", ytitle=paste0(gene_name, " expression [cpm]"), 
                        mainTitle=paste0(gene_name, " expression during treatment"),
                        addDot=T, dotSize=0.5,
                        removePanelGrid=F,removePanelBorder=F,
                        axisLine=c(0.5, "solid", "black")) + 
  stat_compare_means(comparisons = list(c('A','B'),c('A','C')), method = "wilcox.test", size = 5, face="bold") +
  scale_y_continuous(trans='log10')
q <- p + my_theme + theme(plot.title = element_text(color = "black", size = 16, face="bold"))
plot(q)
dev.off()

# Plot driver gene expression for timepoint A-C for REP and NRPS seperatelty
pdf(paste0(fig_dir, "/", gene_name ,"_NR_RE_cpm_Expression.pdf"), height = 5, width = 5)
p <- ggplot2.violinplot(data=single_gene_cpm, xName="Timepoint", yName=paste0(gene_name, "_cpm"),
                        groupName="Timepoint",
                        groupColors= alpha(c('#e63e00','darkseagreen3',('steelblue')), .8), showLegend=FALSE,
                        backgroundColor="white", xtitle="", ytitle=paste0(gene_name, " expression [cpm]"), 
                        mainTitle=paste0(gene_name, " expression during treatment"),
                        addDot=T, dotSize=0.5,
                        removePanelGrid=F,removePanelBorder=F,
                        axisLine=c(0.5, "solid", "black")) + 
                        facet_grid(~Responsiveness , scales = "free_x", space = "free_x",  labeller = as_labeller(c(`NE` = "NRP",`RE` = "REP"))) +
                        theme(panel.spacing = unit(0.5, "cm"),
                        strip.text = element_text(size = 12, face = "bold"))+
                        stat_compare_means(comparisons = list(c('A','B'),c('A','C')), method = "wilcox.test", size = 5, face="bold") +
                        scale_y_continuous(trans='log10')
q <- p + my_theme + theme(plot.title = element_text(color = "black", size = 16, face="bold"))
plot(q)
dev.off()


} 

 