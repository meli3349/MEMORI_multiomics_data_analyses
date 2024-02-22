
# This scripts uses subsampled CD4 and CD8 cells and subtypes and plots: 
# 1) Heterogeneity  of CD4 and CD8 cell counts between multiregion ROIs of the same sample
# 2) Heterogeneity  of GzmB+/PD1+ ratio ofCD4 and CD8 cell between multiregion ROIs of the same sample

library(ggplot2); 
library(dplyr); 
library(ggpubr); 
library(reshape2)

# prepare output dirs:
fig_dir = "~/analysis/IMC/CD4_CD8_subsampled_cells/plots/MultiROI_CD4_CD8_heterogeneity"
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)

# My theme 
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

my_theme_adapted <- theme(
  plot.title = element_text(color="black", size=20, face="bold"),
  axis.title.x = element_text(color="black", size=22, face="bold"),
  axis.text.x = element_text(color="black", size=18, face="bold"),
  axis.title.y = element_text(color="black", size=18, face="bold"),
  axis.text.y = element_text(color="black", size=18, face="bold"),
  legend.title = element_text(color = "black", size = 14, face="bold"),
  legend.text = element_text(color = "black", size = 16, face="bold"),
  
)



# 1) Read in subsampled CD4 or CD8 cells (subsampling performed in OMIQ) 
CD4_8_df <- read.table("~/analysis/IMC/input_files/CD8_CD4_subsets_cell_counts_including_revision_dataset.csv", stringsAsFactors = F, header = T, sep = ",")
CD4_8_df <- CD4_8_df[,c("file", "MEMORI_Sample_ID", "CD8..count", "CD4..count", "CD8.PD1..count", "CD8.GzmB..count", "CD4.PD1..count", "CD4.GzmB..count")]


# Read in area of ROIs calculated in Fiji
area_df_exploratory  <- read.csv("~/analysis/IMC/input_files/Area_exploratory_set_NucCyto.csv")
area_df_revision <- read.csv("~/analysis/IMC/input_files/Area_revision_set_NucCyto.csv")
area_df <- rbind (area_df_exploratory, area_df_revision)
area_df <- area_df[,c("SampleID", "Final_area")]

CD4_8_df <- merge(CD4_8_df, area_df, by.x = "MEMORI_Sample_ID", by.y = "SampleID", all = T)

# Calculation of immune cells per mm2
cell_columns <- c("CD8..count", "CD4..count", "CD8.PD1..count", "CD8.GzmB..count", "CD4.PD1..count", "CD4.GzmB..count")
for(c in cell_columns) {CD4_8_df[,c] <- as.numeric(CD4_8_df[,c])}; CD4_8_df$Final_area <- as.numeric(CD4_8_df$Final_area)
for(c in cell_columns) {CD4_8_df[,c] <- ((CD4_8_df[,c]/CD4_8_df$Final_area*10^6) + 1)}

# Caculate ratio of GzmB+/PD+ CD4 and CD8 cells
CD4_8_df$Ratio_CD8GzmB_CD8PD1 <-(CD4_8_df$CD8.GzmB..count)/(CD4_8_df$CD8.PD1..count)
CD4_8_df$Ratio_CD4GzmB_CD4PD1 <- (CD4_8_df$CD4.GzmB..count)/(CD4_8_df$CD4.PD1..count)

# Select patients with multiregion ROIs
CD4_8_df$patient_TP <- sapply(CD4_8_df$file, function(x) substr(x,8,17)) 
CD4_8_df$patient_TP <- str_replace_all(CD4_8_df$patient_TP, 'NE', 'NRP'); CD4_8_df$patient_TP <- str_replace_all(CD4_8_df$patient_TP, 'RE', 'REP')
duplicated_patient_TP <- CD4_8_df$patient_TP[duplicated(CD4_8_df$patient_TP)]
CD4_8_df_duplicated <- CD4_8_df[which(CD4_8_df$patient_TP %in% duplicated_patient_TP), ]


# Plot heterogeneity of CD4 and CD8 cell count in multi-ROI samples
CD4_8_df_duplicated_to_melt <- CD4_8_df_duplicated[,c("patient_TP", "CD8..count", "CD4..count")]
colnames(CD4_8_df_duplicated_to_melt) <- c("patient_TP", "CD8 cells", "CD4 cells")
df_1 <- melt(CD4_8_df_duplicated_to_melt)  

pdf(paste0(fig_dir, "/MultiROI_CD4_CD8_heterogeneity_including_revision_dataset.pdf"), height = 6, width = 10)
p1<-ggplot(df_1, aes(x=patient_TP, y=as.numeric(value), group=variable)) + geom_point(size=4,alpha=0.7,aes(color=variable,fill=variable)) + theme_mypub_grid()  + 
  theme(axis.text.x = element_text(angle =90)) + 
  scale_y_continuous(trans='log2', limits = c(0.1,2000), labels = trans_format("log2", math_format(2^.x)))+
  scale_color_manual(values = c("CD8 cells" = "red", "CD4 cells" = "dodgerblue2")) +
  xlab("Sample") + ylab(expression(bold("cell count/"~mm^2))) + ggtitle("Heterogeneity CD4 and CD8 cell counts") + my_theme_adapted + theme(legend.position = "bottom", legend.title = element_blank())
p1
dev.off()

# Plot heterogeneity of ratios CD4GzmB+/CD4PD1+ and CD8GzmB+/CD8PD1+ in multi-ROI samples 
CD4_8_phenotype_df_duplicated_to_melt <- CD4_8_df_duplicated[,c("patient_TP", "Ratio_CD8GzmB_CD8PD1", "Ratio_CD4GzmB_CD4PD1")]
CD4_8_phenotype_df_duplicated_to_melt$patient_TP <- str_replace_all(CD4_8_phenotype_df_duplicated_to_melt$patient_TP, 'NE', 'NR')

colnames(CD4_8_phenotype_df_duplicated_to_melt) <- c("patient_TP", "Ratio_CD8GzmB_CD8PD1", "Ratio_CD4GzmB_CD4PD1")
df2 <- melt(CD4_8_activation_df_duplicated_to_melt)  
df2_CD4 <- df2[which(df2$variable == "Ratio_CD4GzmB_CD4PD1"), ]
df2_CD8 <- df2[which(df2$variable == "Ratio_CD8GzmB_CD8PD1"), ]

# CD4
pdf(paste0(fig_dir,"/MultiROI_CD4_activation_status_heterogeneity_including_revision_dataset.pdf"), height = 5, width = 10)
p1<-ggplot(df2_CD4, aes(x=patient_TP, y=as.numeric(value), group=variable)) + geom_point(size=3,alpha=0.5,aes(color=variable,fill=variable), fill = "dodgerblue2", color = "dodgerblue2", pch=21,position=position_jitter(h=0.3, w=0.3)) + theme_mypub_grid() + 
  theme(axis.text.x = element_text(angle =90)) + scale_y_continuous(trans='log2', limits = c(0.1,100), labels = trans_format("log2", math_format(2^.x)))
p1<-p1  + xlab("Sample") + ylab(expression(bold("ratio"))) + ggtitle("Ratio CD4+GzmB+/CD4+PD1+") + my_theme_adapted
p1
dev.off()

# CD8
pdf(paste0(fig_dir, "/MultiROI_CD8_activation_status_heterogeneity_including_revision_dataset.pdf"), height = 5, width = 10)
p1<-ggplot(df2_CD8, aes(x=patient_TP, y=as.numeric(value), group=variable)) + geom_point(size=3,alpha=0.5, aes(color=variable,fill=variable), fill = "red", color = "red",pch=21,position=position_jitter(h=0.3, w=0.3)) + theme_mypub_grid()+ 
  theme(axis.text.x = element_text(angle =90)) + scale_y_continuous(trans='log2', limits = c(0.1,100), labels =trans_format("log2", math_format(2^.x)))
p1<-p1  + xlab("Sample") + ylab(expression(bold("ratio"))) + ggtitle("Ratio CD8+GzmB+/CD8+PD1+") + my_theme_adapted
p1 
dev.off()

