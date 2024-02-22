##################################################################################
## Correlation of CD4 and CD8 cells and their activation status (from IMC data) with immune escape status (from WES and RNA-Seq data)

# This scripts correlates...
# 1a) Number of CD8 cells and immune ecape
# 1b) Number of CD4 cell score and immune ecape
# 2a) CD8 activation status and immune ecape
# 2b) CD4 activation status and immune ecape
# 3a) Number of CD8 and HLA LOH
# 3b) Number of CD4 and HLA LOH
# 4a) CD8 activation status and HLA LOH
# 4b) CD4 activation status and HLA LOH

# prepare output dirs:
fig_dir = "~/analysis/IMC/CD4_CD8_subsampled_cells/plots/correlation_immune_escape/"
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)


# My themes
theme_mypub_grid <- function(base_size = 14,
                             base_family = ""){
  theme_bw(base_size = base_size, 
           base_family = base_family) %+replace%
    theme(
      panel.grid.major = element_blank(),   
      panel.grid.minor = element_blank(), 
      panel.background = element_blank(),
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

# Read in CD4 and CD8 cell counts from IMC analyses for each region of interest (ROI)
CD4_8_df <- read.table("~/analysis/IMC/input_files/CD8_CD4_subsets_cell_counts_exploratory_set.csv", header = T, stringsAsFactors = F, sep = ",")
CD4_8_df <- CD4_8_df[,c("file", "CD8..count", "CD4..count", "CD8.GzmB..count", "CD8.PD1..count", "CD4.GzmB..count", "CD4.PD1..count")]
CD4_8_df$SampleID <- sapply(CD4_8_df$file, function(x) substr(x,1, 26))

# Read in area of ROIs calculated in Fiji
area_df <- read.csv("~/analysis/IMC/input_files/Area_exploratory_set_NucCyto.csv")
area_df <- area_df[,c("SampleID", "Final_area")]

# Calculation of immune cells per mm2
CD4_8_df <- merge(CD4_8_df, area_df, by.x = "SampleID", by.y = "SampleID")
cell_columns <- c("CD8..count", "CD4..count", "CD8.GzmB..count", "CD8.PD1..count", "CD4.GzmB..count", "CD4.PD1..count")
for(c in cell_columns) {CD4_8_df[,c] <- as.numeric(CD4_8_df[,c])};CD4_8_df$Final_area <- as.numeric(CD4_8_df$Final_area)
for(c in cell_columns) {CD4_8_df[,c] <- CD4_8_df[,c]/CD4_8_df$Final_area*10^6}
CD4_8_df$patient_TP <- sapply(CD4_8_df$file, function(x) substr(x,8,17)) 

# Calculate mean cell counts for multi-ROI samples
mean_CD4_8_df = aggregate(CD4_8_df,by = list(CD4_8_df$patient_TP), FUN = mean)

# Caculate activation status of CD4 and CD8 cells (ratio of GzmB+/PD1+ CD4 and CD8 cells)
mean_CD4_8_df$ratio_CD8GzmB_CD8PDL1 <- as.numeric((mean_CD4_8_df$CD8.GzmB..count+1)/(mean_CD4_8_df$CD8.PD1..count+1))
mean_CD4_8_df$ratio_CD4GzmB_CD4PDL1 <- as.numeric((mean_CD4_8_df$CD4.GzmB..count+1)/(mean_CD4_8_df$CD4.PD1..count+1))
colnames(mean_CD4_8_df)[1] <- "patient_TP"
mean_CD4_8_df <- mean_CD4_8_df[,c("patient_TP", "CD8..count", "CD4..count", "CD8.GzmB..count", "CD8.PD1..count", "CD4.GzmB..count", "CD4.PD1..count", "ratio_CD8GzmB_CD8PDL1", "ratio_CD4GzmB_CD4PDL1")]


# Read in table with immune escape results 
escape_df <- read.table("~/analysis/multi_omic/immune_escape/tables/MEMORI_escape_master_file_incl_PDL1.txt", stringsAsFactors = F, header = T, sep = "\t")
escape_df <- escape_df[complete.cases(escape_df[ ,c("WES_file", "RNA_file")]),]
escape_df$patient_TP <- sapply(escape_df$Sample, function(x) substr(x,8,17)) 


# Merge immune escape results with with IMC data 
merged_df <- merge(escape_df, mean_CD4_8_df, by.x = "patient_TP", by.y = "patient_TP")
merged_df$Timepoint <- sapply(merged_df$Sample, function(x) substr(x,14,14))


## 1-2.) Correlation between immune escape and CD4 and CD8 cell counts/ their activation status
# Categorization of samples into high and low immune infiltration/ activation  
y <- merged_df 
y$CD8..count <- as.numeric(y$CD8..count); y$CD4..count <- as.numeric(y$CD4..count)
y$CD8..count_category <- ifelse(y$CD8..count > median(y$CD8..count, na.rm = T), "high", "low")
y$CD4..count_category <- ifelse(y$CD4..count > median(y$CD4..count, na.rm = T), "high", "low")
y$ratio_CD8GzmB_CD8PDL1_category <- ifelse(y$ratio_CD8GzmB_CD8PDL1 > median(y$ratio_CD8GzmB_CD8PDL1, na.rm = T), "high", "low")
y$ratio_CD4GzmB_CD4PDL1_category <- ifelse(y$ratio_CD4GzmB_CD4PDL1 > median(y$ratio_CD4GzmB_CD4PDL1, na.rm = T), "high", "low")



# 1a) Plotting CD8 cell scores and immune escape 
df_CD8 <- y[,c("final_escape_result", "CD8..count_category")]

Sum_df <- t(table(df_CD8))
chisq_CD8 <- chisq.test(Sum_df)
Sum_df <- scale(Sum_df, FALSE, colSums(Sum_df)) * 100
melt_Sum_df <- melt(Sum_df)
melt_Sum_df$CD8..count_category  = factor(melt_Sum_df$CD8..count_category, levels=c("low", "high"))



pdf(paste0(fig_dir, "/IMC_quantification_CD8cells_immune_escape.pdf"), width = 6, height = 4)
q <- ggplot(melt_Sum_df, aes(final_escape_result, value, fill= CD8..count_category)) + geom_bar(stat= "Identity") + theme(axis.text.x = element_text(size = 14, face = "bold", color = "black")) + 
  theme(axis.title.y = element_text(size = 14, face = "bold"))+ theme(axis.text.y = element_text(size = 14, face = "bold", color ="black")) + theme(legend.text=element_text(size= 14, face = "bold"))+ theme(legend.title = element_text(size= 16, face = "bold")) + 
  scale_fill_manual(values = c("grey", "lightblue"))  + 
  ylab("percentage of samples") + xlab("Immune escaped") + labs(fill ="IMC CD8 cell score") + theme_mypub_grid() + my_theme +
  scale_y_continuous(limits = c(0,120))
plot(q)
dev.off()


# 1b) Plotting CD4 cell scores and immune escape 
df_CD4 <- y[,c("final_escape_result", "CD4..count_category")]

Sum_df <- t(table(df_CD4))
chisq_CD4 <- chisq.test(Sum_df)
Sum_df <- scale(Sum_df, FALSE, colSums(Sum_df)) * 100
melt_Sum_df <- melt(Sum_df)
melt_Sum_df$CD4..count_category  = factor(melt_Sum_df$CD4..count_category, levels=c("low", "high"))


pdf(paste0(fig_dir, "/IMC_quantification_CD4cells_immune_escape.pdf"), width = 6, height = 4)
q <- ggplot(melt_Sum_df, aes(final_escape_result, value, fill= CD4..count_category)) + geom_bar(stat= "Identity") + theme(axis.text.x = element_text(size = 14, face = "bold", color = "black")) + 
  theme(axis.title.y = element_text(size = 14, face = "bold"))+ theme(axis.text.y = element_text(size = 14, face = "bold", color ="black")) + theme(legend.text=element_text(size= 14, face = "bold"))+ theme(legend.title = element_text(size= 16, face = "bold")) + 
  scale_fill_manual(values = c("grey", "lightblue"))  +
  ylab("percentage of samples") + xlab("Immune escaped") + labs(fill = "IMC CD4 cell score") + theme_mypub_grid() + my_theme 
plot(q)
dev.off()

# 2a) Plotting CD8 activation status and immune escape
df_ratio_CD8GzmB_CD8PDL1 <- y[,c("final_escape_result", "ratio_CD8GzmB_CD8PDL1_category")]
Sum_df <- t(table(df_ratio_CD8GzmB_CD8PDL1))
chisq_activation_status_CD8 <- chisq.test(Sum_df)
Sum_df <- scale(Sum_df, FALSE, colSums(Sum_df)) * 100
melt_Sum_df <- melt(Sum_df)

melt_Sum_df$ratio_CD8GzmB_CD8PDL1_category = factor(melt_Sum_df$ratio_CD8GzmB_CD8PDL1_category, levels=c("low", "high"))

pdf(paste0(fig_dir, "/IMC_activation_status_CD8cells_immune_escape.pdf"), width = 6, height = 4)
q <- ggplot(melt_Sum_df, aes(final_escape_result, value, fill= ratio_CD8GzmB_CD8PDL1_category)) + geom_bar(stat= "Identity") + theme(axis.text.x = element_text(size = 14, face = "bold", color = "black")) + 
  theme(axis.title.y = element_text(size = 14, face = "bold"))+ theme(axis.text.y = element_text(size = 14, face = "bold", color ="black")) + theme(legend.text=element_text(size= 14, face = "bold"))+ theme(legend.title = element_text(size= 16, face = "bold")) + 
  scale_fill_manual(values = c("grey", "lightblue", "aquamarine3"))  +
  ylab("percentage of samples") + xlab("Immune escape") + labs(fill = bquote(bold("ratio CD8"^"+"~"GzmB"^"+"~"/ CD8"^"+"~"PD1"^"+"))) + theme_mypub_grid() + my_theme 
plot(q)
dev.off()


# 2b) Plotting CD4 activation status and immune escape
df_ratio_CD4GzmB_CD4PDL1 <- y[,c("final_escape_result", "ratio_CD4GzmB_CD4PDL1_category")]
Sum_df <- t(table(df_ratio_CD4GzmB_CD4PDL1))
chisq_activation_status_CD4 <- chisq.test(Sum_df)
Sum_df <- scale(Sum_df, FALSE, colSums(Sum_df)) * 100
melt_Sum_df <- melt(Sum_df)

melt_Sum_df$ratio_CD4GzmB_CD4PDL1_category = factor(melt_Sum_df$ratio_CD4GzmB_CD4PDL1_category, levels=c("low", "high"))

pdf(paste0(fig_dir, "/IMC_activation_status_CD4cells_immune_escape.pdf"),  width = 6, height = 4)
q <- ggplot(melt_Sum_df, aes(final_escape_result, value, fill= ratio_CD4GzmB_CD4PDL1_category)) + geom_bar(stat= "Identity") + 
  scale_fill_manual(values = c("grey", "lightblue", "aquamarine3"))  +
  ylab("percentage of samples") + xlab("Immune escape") + labs(fill = bquote(bold("ratio CD4"^"+"~"GzmB"^"+"~"/ CD4"^"+"~"PD1"^"+"))) + theme_mypub_grid() + my_theme +
  scale_y_continuous(limits = c(0,120))
plot(q)
dev.off()


## 3-4.) Correlation between HLA LOH state and CD4 and CD8 cell counts/ their activation status
# Create dataframe with samples with/ without HLA LOH, but no other immune escape mechanism
# and categorization of samples into high and low immune infiltration/ activation
y_HLA_LOH_without_other_escape <- merged_df[which(y$PDL1_immune_escape  == "FALSE" & y$HLAmut == "FALSE" & y$B2Mmut == "FALSE"), ]
y_HLA_LOH_without_other_escape$CD8..count <- as.numeric(y_HLA_LOH_without_other_escape$CD8..count); y_HLA_LOH_without_other_escape$CD4..count <- as.numeric(y_HLA_LOH_without_other_escape$CD4..count)
y_HLA_LOH_without_other_escape$CD8..count_category <- ifelse(y_HLA_LOH_without_other_escape$CD8..count > median(y_HLA_LOH_without_other_escape$CD8..count, na.rm = T), "high", "low")
y_HLA_LOH_without_other_escape$CD4..count_category <- ifelse(y_HLA_LOH_without_other_escape$CD4..count > median(y_HLA_LOH_without_other_escape$CD4..count, na.rm = T), "high", "low")
y_HLA_LOH_without_other_escape$ratio_CD8GzmB_CD8PDL1_category <- ifelse(y_HLA_LOH_without_other_escape$ratio_CD8GzmB_CD8PDL1 > median(y_HLA_LOH_without_other_escape$ratio_CD8GzmB_CD8PDL1, na.rm = T), "high", "low")
y_HLA_LOH_without_other_escape$ratio_CD4GzmB_CD4PDL1_category <- ifelse(y_HLA_LOH_without_other_escape$ratio_CD4GzmB_CD4PDL1 > median(y_HLA_LOH_without_other_escape$ratio_CD4GzmB_CD4PDL1, na.rm = T), "high", "low")

# 3a) Plotting CD8 cell scores and HLA LOH
df_CD8 <- y_HLA_LOH_without_other_escape[,c("final_HLA_LOH", "CD8..count_category")]

Sum_df <- t(table(df_CD8))
chisq_CD8 <- chisq.test(Sum_df)
Sum_df <- scale(Sum_df, FALSE, colSums(Sum_df)) * 100
melt_Sum_df <- melt(Sum_df)
melt_Sum_df$CD8..count_category  = factor(melt_Sum_df$CD8..count_category, levels=c("low", "high"))


pdf(paste0(fig_dir, "/IMC_quantification_CD8cells_HLALOH.pdf"),  width = 6, height = 4)
q <- ggplot(melt_Sum_df, aes(final_HLA_LOH, value, fill= CD8..count_category)) + geom_bar(stat= "Identity") + theme(axis.text.x = element_text(size = 14, face = "bold", color = "black")) + 
  theme(axis.title.y = element_text(size = 14, face = "bold"))+ theme(axis.text.y = element_text(size = 14, face = "bold", color ="black")) + theme(legend.text=element_text(size= 14, face = "bold"))+ theme(legend.title = element_text(size= 16, face = "bold")) + 
  scale_fill_manual(values = c("grey", "lightblue", "aquamarine3"))  +
  ylab("percentage of samples") + xlab("HLA LOH") + labs(fill = "IMC CD8 cell score") + theme_mypub_grid() + my_theme 
plot(q)
dev.off()

# 3b) Plotting CD4 cell scores and HLA LOH
df_CD4 <-  y_HLA_LOH_without_other_escape[,c("final_HLA_LOH", "CD4..count_category")]
Sum_df <- t(table(df_CD4))
chisq_CD4 <- chisq.test(Sum_df)
Sum_df <- scale(Sum_df, FALSE, colSums(Sum_df)) * 100
melt_Sum_df <- melt(Sum_df)
melt_Sum_df$CD4..count_category  = factor(melt_Sum_df$CD4..count_category, levels=c("low", "high"))


pdf(paste0(fig_dir, "/IMC_quantification_CD4cells_HLALOH.pdf"),  width = 6, height = 4)
q <- ggplot(melt_Sum_df, aes(final_HLA_LOH, value, fill= CD4..count_category)) + geom_bar(stat= "Identity") + theme(axis.text.x = element_text(size = 14, face = "bold", color = "black")) + 
  theme(axis.title.y = element_text(size = 14, face = "bold"))+ theme(axis.text.y = element_text(size = 14, face = "bold", color ="black")) + theme(legend.text=element_text(size= 14, face = "bold"))+ theme(legend.title = element_text(size= 16, face = "bold")) + 
  scale_fill_manual(values = c("grey", "lightblue", "aquamarine3"))  +
  ylab("percentage of samples") + xlab("HLA LOH") + labs(fill = "IMC CD4 cell score") + theme_mypub_grid() + my_theme 
plot(q)
dev.off()


# 4a) Plotting CD8 activation status and HLA LOH
df_ratio_CD8GzmB_CD8PDL1 <- y_HLA_LOH_without_other_escape[,c("final_HLA_LOH", "ratio_CD8GzmB_CD8PDL1_category")]
Sum_df <- t(table(df_ratio_CD8GzmB_CD8PDL1))
chisq_activation_status_CD8 <- chisq.test(Sum_df)
Sum_df <- scale(Sum_df, FALSE, colSums(Sum_df)) * 100
melt_Sum_df <- melt(Sum_df)
melt_Sum_df$ratio_CD8GzmB_CD8PDL1_category = factor(melt_Sum_df$ratio_CD8GzmB_CD8PDL1_category, levels=c("low", "high"))

pdf(paste0(fig_dir, "/IMC_activation_status_CD8cells_HLALOH.pdf"), width = 6, height = 4)
q <- ggplot(melt_Sum_df, aes(final_HLA_LOH, value, fill= ratio_CD8GzmB_CD8PDL1_category)) + geom_bar(stat= "Identity") + 
  scale_fill_manual(values = c("grey", "lightblue", "aquamarine3"))  +
  ylab("percentage of samples") + xlab("HLA LOH") + labs(fill = bquote(bold("Ratio CD8"^"+"~"GzmB"^"+"~"/ CD8"^"+"~"PD1"^"+"))) + theme_mypub_grid() + my_theme +
  scale_y_continuous(limits = c(0,120))
plot(q)
dev.off()


# 4b) Plotting CD4 activation status and HLA LOH
df_ratio_CD4GzmB_CD4PDL1 <- y_HLA_LOH_without_other_escape[,c("final_HLA_LOH", "ratio_CD4GzmB_CD4PDL1_category")]
Sum_df <- t(table(df_ratio_CD4GzmB_CD4PDL1))
chisq_activation_status_CD4 <- chisq.test(Sum_df)
Sum_df <- scale(Sum_df, FALSE, colSums(Sum_df)) * 100
melt_Sum_df <- melt(Sum_df)
melt_Sum_df$ratio_CD4GzmB_CD4PDL1_category = factor(melt_Sum_df$ratio_CD4GzmB_CD4PDL1_category, levels=c("low", "high"))

pdf(paste0(fig_dir, "/IMC_activation_status_CD4cells_HLALOH.pdf"), width = 6, height = 4)
q <- ggplot(melt_Sum_df, aes(final_HLA_LOH, value, fill= ratio_CD4GzmB_CD4PDL1_category)) + geom_bar(stat= "Identity") + theme(axis.text.x = element_text(size = 14, face = "bold", color = "black")) + 
  theme(axis.title.y = element_text(size = 14, face = "bold"))+ theme(axis.text.y = element_text(size = 14, face = "bold", color ="black")) + theme(legend.text=element_text(size= 14, face = "bold"))+ theme(legend.title = element_text(size= 16, face = "bold")) + 
  scale_fill_manual(values = c("grey", "lightblue", "aquamarine3"))  + 
  ylab("percentage of samples") + xlab("HLA LOH") + labs(fill =bquote(bold("Ratio CD4"^"+"~"GzmB"^"+"~"/ CD4"^"+"~"PD1"^"+"))) + theme_mypub_grid() + my_theme +
  scale_y_continuous(limits = c(0,120))
plot(q)
dev.off()

