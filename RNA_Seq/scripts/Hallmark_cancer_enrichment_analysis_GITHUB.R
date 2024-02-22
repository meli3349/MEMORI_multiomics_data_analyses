

#########################################################################################
## Expression analyses using the Hallmark of Cancer gene set ############################
#########################################################################################
library(GSVAdata)
library(GSVA)
data(c2BroadSets)
data(commonPickrellHuang)
library(msigdbr)
library(biomaRt)
library(ggplot2)

# prepare output dirs:
fig_dir = "~/analysis/RNA_Seq/pathway_analyses/Hallmark_cancer_pathway/Plots"
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)
data_dir = "~/analysis/RNA_Seq/pathway_analyses/Hallmark_cancer_pathway/Data_frames"
dir.create(data_dir, showWarnings=FALSE, recursive=TRUE)

# Read in normalized GeneIndex_MasterTables
logcpm <- read.table("~/analysis/RNA_Seq/input_files/GeneExpression_Mastertables/GeneIndex_MasterTable1.5M_norm_batch_1_2_3_filtered_expressed_genes.txt", header = T, stringsAsFactors = F, sep = "\t")
logcpm_nonfiltered <- read.table("~/analysis/RNA_Seq/input_files/GeneExpression_Mastertables/GeneIndex_MasterTable1.5M_norm_batch_1_2_3.txt_all_genes.txt", header = T, stringsAsFactors = F, sep = "\t")


# Load gene set for enrichment analyses
# H: Cancer hallmark gene sets from GSEA-msigdb
# https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

m_df = msigdbr(species = "Homo sapiens", category = "H")
Hallmark_cancer_list = split(x = m_df$entrez_gene, f = m_df$gs_name)


# Get the full ensembl dataset
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

# Dataframe of matches between ensembl, symbol and entrez
geneMap = getBM(attributes=c("ensembl_gene_id","entrezgene_id","hgnc_symbol"), mart = ensembl,filters="ensembl_gene_id",values = logcpm$Gene_ENSEMBLE)

# Get the Entrez gene ID for normalized gene expression df
logcpm_with_entrez <- merge(geneMap, logcpm, by.x = "ensembl_gene_id", by.y = "Gene_ENSEMBLE", all = T)

# exclude genes without entrez, as these can not be included in GSEA enrichment analysis 
logcpm_with_entrez <- logcpm_with_entrez[!is.na(logcpm_with_entrez$entrezgene_id), ]

# Identify duplicated entrez IDs and exclude then
duplicates <- logcpm_with_entrez[which(duplicated(logcpm_with_entrez$entrezgene_id)), 'entrezgene_id']

logcpm_with_entrez <- logcpm_with_entrez[which(!duplicated(logcpm_with_entrez$entrezgene_id)),];row.names(logcpm_with_entrez) <- c(1:nrow(logcpm_with_entrez))
logcpm_with_entrez_unique <- logcpm_with_entrez[order(logcpm_with_entrez$entrezgene_id),];row.names(logcpm_with_entrez_unique) <- logcpm_with_entrez_unique$entrezgene_id


# Calculate Geneset enrichment score per sample 
# with Cancer hallmark gene sets from GSEA-msigdb
Hallmark_cancer_per_sample <- gsva(data.matrix(logcpm_with_entrez_unique[, c(4:ncol(logcpm_with_entrez_unique))]),
                              Hallmark_cancer_list, min.sz=5, max.sz=500, mx.diff=T,
                              verbose=FALSE, parallel.sz=1)

Hallmark_df <- t(Hallmark_cancer_per_sample)
write.table(Hallmark_df, file = paste0(data_dir, "/Hallmark_cancer_scores.txt"), sep = "\t", row.names = T, quote = F)


##### Plot PCA using Hallmark cancer genes

# theme for plots 
my_theme <- theme(
  plot.title = element_text(color="black", size=16, face="bold"),
  axis.title.x = element_text(color="black", size=16, face="bold"),
  axis.text.x = element_text(color="black", size=16, face="bold"),
  axis.title.y = element_text(color="black", size=16, face="bold"),
  axis.text.y = element_text(color="black", size=16, face="bold"),
  legend.title = element_text(color = "black", size = 14, face="bold"),
  legend.text = element_text(color = "black", size = 14, face="bold"),
  
)

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

# Principal component analysis
res_hallmark_cancer <- prcomp(t(Hallmark_cancer_per_sample));sum_hallmark_cancer <- summary(res_hallmark_cancer)

# Vector with patients
patients <- gsub('MEMEX_(\\S+\\d+)\\S+_\\S\\d_\\S\\d_\\S\\d','\\1',row.names(res_hallmark_cancer$x))

# PCA with PC1 and PC2
pca_df <- data.frame(PC1 = res_hallmark_cancer$x[, 1], PC2 = res_hallmark_cancer$x[, 2], SampleID = rownames(res_hallmark_cancer$x))
merged_pca_df <- merge(pca_df,clin.df_1.5M, by.x = "row.names", by.y = "row.names")
colnames(merged_pca_df)[which(colnames(merged_pca_df) == "responsiveness")] <- "Responsiveness"

# PCA axis labeling: Use % of explained variance to label x and y axis
labx_hall_cancer <- paste0('PC1 (',signif(sum_hallmark_cancer$importance[2,1]*100,3),'% explained variance)')
laby_hall_cancer <- paste0('PC2 (',signif(sum_hallmark_cancer$importance[2,2]*100,3),'% explained variance)')
color_NE_RE <- c("#b9adde", "#26A69A")


# Plot PCA using Hallmark cancer genes for all samples
# a.) color per patient
pdf(paste0(fig_dir, "/PCA_Hallmark_cancer_patient_colored.pdf"), height = 9, width = 12)
p <- ggplot(data = merged_pca_df, aes_string(x = "PC1", y = "PC2", color="patient", shape="Timepoint", size = 4)) + geom_point() +
  xlab(labx_hall_cancer) + ylab(laby_hall_cancer) + my_theme + theme(legend.text = element_text(size = 8, face="bold"))
plot(p)
dev.off()

# b.) color per responsiveness group
pdf(paste0(fig_dir, "/PCA_Hallmark_cancer_responsiveness_colored.pdf"), height = 5, width = 7)
r <- ggplot(data = merged_pca_df, aes_string(x = "PC1", y = "PC2", color="Responsiveness", shape="Timepoint")) + geom_point(size = 5) +
  xlab(labx_hall_cancer) + ylab(laby_hall_cancer) +theme_mypub_grid()+ my_theme + scale_color_manual(values = color_NE_RE, labels = c("NRP", "REP"))
plot(r)
dev.off()

# c.) color per sequencing batch
pdf(paste0(fig_dir, "/PCA_Hallmark_cancer_batch_colored.pdf"), height = 5, width = 7)
b <- ggplot(data = merged_pca_df, aes_string(x = "PC1", y = "PC2", color="Batch", shape="Timepoint", size = 4)) + geom_point() +
  xlab(labx_hall_cancer) + ylab(laby_hall_cancer)  +theme_bw() + my_theme
plot(b)
dev.off()

##### Plot Loadings for PCAs
Hallmark_cancer_df <- read.delim("~/analysis/General_files/external_files/Hallmark_cancer_pathway.txt", header = T, stringsAsFactors = F)

Hallmark_cancer_df$Label_name <- ifelse(Hallmark_cancer_df$Pathway_category == "other", NA, Hallmark_cancer_df$Pathway_short_name)
res_hallmark_cancer_for_loading <- as.data.frame(res_hallmark_cancer$rotation)
res_hallmark_cancer_for_loading$Hallmark_Cancer_Pathway <- rownames(res_hallmark_cancer_for_loading)
df_loading <- merge(res_hallmark_cancer_for_loading, Hallmark_cancer_df, by.x = "Hallmark_Cancer_Pathway", by.y = "Hallmark_Cancer_Pathway")
rownames(df_loading) <- df_loading[,1]
df_loading <- df_loading[,-1]
df_loading$color <- ifelse(df_loading$Pathway_category == "oncogenic", "firebrick3", ifelse(df_loading$Pathway_category == "immune", "forestgreen", ifelse(df_loading$Pathway_category == "cellular stress", "orange", ifelse(df_loading$Pathway_category == "stromal", "steelblue", "grey"))))

#Plot
pdf(paste0(fig_dir, "/Loading_Hallmark_cancer_unlabeled.pdf"), height = 5, width = 7)
plot(c(-0.35,0.35),c(-0.1,0.35),axes=F, xlab='PC1 loadings',ylab='PC2 loadings', cex.lab=1.5, cex.names = 5.5, cex.axis=7.5, cex.main=1.5, cex.sub=2.5, font.lab =2, font.axis =2)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="white",border=T)
grid(nx = NULL, ny = NULL, lty=1,col='grey91');par(new=T)
axis(side=1,font=2);axis(side=2,font=2)
arrows(x0=0,y0=0,x1=df_loading[,1],y1=df_loading[,2], col = df_loading$color, lwd = 2)
text(x=df_loading[,1],y=df_loading[,2], labels = "", cex=1.2, col = "black")  + theme_bw() + my_theme
dev.off()

pdf(paste0(fig_dir, "/Loading_Hallmark_cancer_labeled.pdf"), height = 15, width = 15)
plot(c(-0.35,0.35),c(-0.1,0.35),axes=F, xlab='PC1 loadings',ylab='PC2 loadings', cex.lab=1.5, cex.names = 5.5, cex.axis=7.5, cex.main=1.5, cex.sub=2.5, font.lab =2, font.axis =2)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="white",border=T)
grid(nx = NULL, ny = NULL, lty=1,col='grey91');par(new=T)
axis(side=1,font=2);axis(side=2,font=2)
arrows(x0=0,y0=0,x1=df_loading[,1],y1=df_loading[,2], col = df_loading$color, lwd = 2)
text(x=df_loading[,1],y=df_loading[,2], labels = df_loading$Label_name, cex=1.2, col = "black")  + theme_bw() + my_theme
dev.off()



##### Plot PCA using Hallmark cancer genes for treatment-naive biospies (only Timepoint A (TPA))

# select samples TPA
clinical_df <- data.frame(SampleID = colnames(logcpm_with_entrez_unique)[4:ncol(logcpm_with_entrez_unique)], Timepoint = "x")
clinical_df$Timepoint <- sapply(clinical_df$SampleID, function(x) substr (x, 13,13))
TPA_df <- clinical_df[which(clinical_df$Timepoint == "A"), ]
TPA_samples <- TPA_df$SampleID
TPA_df_subset <- logcpm_with_entrez_unique[, TPA_samples]

# Gene enrichment analyses with H: Cancer hallmark gene sets from GSEA-msigdb
Hallmark_cancer_TPA_sample <- gsva(data.matrix(TPA_df_subset),
                                      Hallmark_cancer_list, min.sz=5, max.sz=500, mx.diff=TRUE,
                                      verbose=FALSE, parallel.sz=1)



# PCA using Hallmark cancer genes for treatment-naive biospies (only Timepoint A)
res_hallmark_cancer_TPA <- prcomp(t(Hallmark_cancer_TPA_sample));sum_hallmark_cancer_TPA <- summary(res_hallmark_cancer_TPA)

pca_df_TPA <- data.frame(PC1 = res_hallmark_cancer_TPA$x[, 1], PC2 = res_hallmark_cancer_TPA$x[, 2], SampleID = rownames(res_hallmark_cancer_TPA$x))
pca_df_TPA$Responsiveness <- sapply(pca_df_TPA$SampleID, function(x) substr(x,7,8))
pca_df_TPA$Timepoint <- sapply(pca_df_TPA$SampleID, function(x) substr(x,13,13))
pca_df_TPA$PatientNumber <- sapply(pca_df_TPA$SampleID, function(x) substr(x,7,12))

# PCA axis labeling: Use % of explained variance to label x and y axis
labx_hall_cancer_TPA <- paste0('PC1 (',signif(sum_hallmark_cancer_TPA$importance[2,1]*100,3),'% explained variance)')
laby_hall_cancer_TPA <- paste0('PC2 (',signif(sum_hallmark_cancer_TPA$importance[2,2]*100,3),'% explained variance)')


# Plot PCA using only Hallmark cancer genes
# b.) color per responsiveness
pdf(paste0(fig_data, "/PCA_Hallmark_cancer_only_TPA_responsiveness.pdf"), height = 5, width = 7)
p_TPA <- ggplot(data =pca_df_TPA, aes_string(x = "PC1", y = "PC2", color="Responsiveness", shape="Timepoint")) + geom_point(size = 5) +
  xlab(labx_hall_cancer_TPA) + ylab(laby_hall_cancer_TPA) +theme_mypub_grid()+ my_theme + scale_color_manual(values = color_NE_RE, labels = c("NRP", "REP"))
plot(p_TPA)
dev.off()

# d.) color per SUVmax value in screening PET CT
PET_SUV_df <- read.csv("~/analysis/General_files/MEMORI_specific/PET_CT_SUV_MEMORI.csv", header = T)
PET_SUVmax_df <- PET_SUV_df[,c("Patient.ID", "SUV.max.in.Screening.PET", "SUV.max.of.2nd.PET")]
PET_SUVmax_df$SUV.max.in.Screening.PET <- as.numeric(PET_SUVmax_df$SUV.max.in.Screening.PET); PET_SUVmax_df$SUV.max.of.2nd.PET <- as.numeric(PET_SUVmax_df$SUV.max.of.2nd.PET)
 
PETmax_quartile_df <- as.data.frame(quantile(PET_SUVmax_df$SUV.max.in.Screening.PET, na.rm= TRUE)); colnames(PETmax_quartile_df)[1] <- "SUVmax_ScreeningPET_quartiles"
Q0_SUVmax_screening <- PETmax_quartile_df[which(rownames(PETmax_quartile_df) == "0%"), 1]
Q25_SUVmax_screening <- PETmax_quartile_df[which(rownames(PETmax_quartile_df) == "25%"), 1]
Q50_SUVmax_screening <- PETmax_quartile_df[which(rownames(PETmax_quartile_df) == "50%"), 1]
Q75_SUVmax_screening <- PETmax_quartile_df[which(rownames(PETmax_quartile_df) == "75%"), 1]
Q100_SUVmax_screening <- PETmax_quartile_df[which(rownames(PETmax_quartile_df) == "100%"), 1]

PET_SUVmax_df$Quantile_SUVmax <- ifelse((PET_SUVmax_df$SUV.max.in.Screening.PET >= Q0_SUVmax_screening & PET_SUVmax_df$SUV.max.in.Screening.PET <= Q25_SUVmax_screening), "Q1", 
                                          ifelse((PET_SUVmax_df$SUV.max.in.Screening.PET > Q25_SUVmax_screening & PET_SUVmax_df$SUV.max.in.Screening.PET <= Q50_SUVmax_screening), "Q2", 
                                                 ifelse((PET_SUVmax_df$SUV.max.in.Screening.PET > Q50_SUVmax_screening & PET_SUVmax_df$SUV.max.in.Screening.PET <= Q75_SUVmax_screening), "Q3", 
                                                        ifelse((PET_SUVmax_df$SUV.max.in.Screening.PET > Q75_SUVmax_screening & PET_SUVmax_df$SUV.max.in.Screening.PET <= Q100_SUVmax_screening), "Q4", NA))))
                                          
pca_df_TPA2 <- merge(pca_df_TPA, PET_SUVmax_df, by.x = "PatientNumber", by.y = "Patient.ID")
color_SUVmax <- c("#FFEB3B", "#FFC107", "#FF9800", "#FF5722")

pdf(paste0(fig_data, "/PCA_Hallmark_cancer_only_TPA_SUVmaxQuantiles_ScreeningPET.pdf"), height = 5, width = 7)
p_TPA <- ggplot(data = pca_df_TPA2, aes_string(x = "PC1", y = "PC2", color="Quantile_SUVmax", shape="Responsiveness")) + geom_point(size = 5) +
  xlab(labx_hall_cancer_TPA) + ylab(laby_hall_cancer_TPA) +theme_mypub_grid()+ my_theme + scale_color_manual(values = color_SUVmax, labels = c("Q1", "Q2", "Q3", "Q4")) + scale_shape_manual(values = c(16,17), labels = c("NRP", "REP"))
plot(p_TPA)
dev.off()

# e.) color per absolute delta SUVmax between screening PET CT and 2nd PET
PET_SUVmax_df$Absolute_delta_PET_CTs <- PET_SUVmax_df$SUV.max.in.Screening.PET - PET_SUVmax_df$SUV.max.of.2nd.PET

Absolute_delta_PETmax_quartile_df <- as.data.frame(quantile((PET_SUVmax_df$SUV.max.in.Screening.PET-PET_SUVmax_df$SUV.max.of.2nd.PET), na.rm= TRUE)); colnames(PETmax_quartile_df)[1] <- "Quantiles_of_absolute_Delta_SUVmax_between_PET_CTs"
Q0_absolute_delta_SUVmax <- Absolute_delta_PETmax_quartile_df[which(rownames(Absolute_delta_PETmax_quartile_df) == "0%"), 1]
Q25_absolute_delta_SUVmax  <- Absolute_delta_PETmax_quartile_df[which(rownames(Absolute_delta_PETmax_quartile_df) == "25%"), 1]
Q50_absolute_delta_SUVmax  <- Absolute_delta_PETmax_quartile_df[which(rownames(Absolute_delta_PETmax_quartile_df) == "50%"), 1]
Q75_absolute_delta_SUVmax  <- Absolute_delta_PETmax_quartile_df[which(rownames(Absolute_delta_PETmax_quartile_df) == "75%"), 1]
Q100_absolute_delta_SUVmax  <- Absolute_delta_PETmax_quartile_df[which(rownames(Absolute_delta_PETmax_quartile_df) == "100%"), 1]

PET_SUVmax_df$Quantile_absolute_delta_SUVmax <- ifelse((PET_SUVmax_df$Absolute_delta_PET_CTs >= Q0_absolute_delta_SUVmax & PET_SUVmax_df$Absolute_delta_PET_CTs <= Q25_absolute_delta_SUVmax), "Q1", 
                                        ifelse((PET_SUVmax_df$Absolute_delta_PET_CTs > Q25_absolute_delta_SUVmax & PET_SUVmax_df$Absolute_delta_PET_CTs <= Q50_absolute_delta_SUVmax), "Q2", 
                                               ifelse((PET_SUVmax_df$Absolute_delta_PET_CTs > Q50_absolute_delta_SUVmax & PET_SUVmax_df$Absolute_delta_PET_CTs <= Q75_absolute_delta_SUVmax), "Q3", 
                                                      ifelse((PET_SUVmax_df$Absolute_delta_PET_CTs > Q75_absolute_delta_SUVmax & PET_SUVmax_df$Absolute_delta_PET_CTs <= Q100_absolute_delta_SUVmax), "Q4", NA))))

merged_pca_df_TPA2 <- merge(pca_df_TPA, PET_SUVmax_df, by.x = "PatientNumber", by.y = "Patient.ID")
color_absolute_delta_SUVmax <- c("#FFEB3B", "#CDDC39", "#BBC34A", "#4CAF50")

pdf(paste0(fig_data, "/PCA_Hallmark_cancer_only_TPA_Quantiles_absolute_SUVmax_delta_between_PETs.pdf"), height = 5, width = 10)
p_TPA <- ggplot(data = merged_pca_df_TPA2, aes_string(x = "PC1", y = "PC2", color="Quantile_absolute_delta_SUVmax", shape="Responsiveness")) + geom_point(size = 5) +
  xlab(labx_hall_cancer_TPA) + ylab(laby_hall_cancer_TPA) +theme_mypub_grid()+ my_theme + scale_color_manual(values = color_absolute_delta_SUVmax, labels = c("Q1", "Q2", "Q3", "Q4")) + scale_shape_manual(values = c(16,17), labels = c("NRP", "REP"))
plot(p_TPA)
dev.off()
