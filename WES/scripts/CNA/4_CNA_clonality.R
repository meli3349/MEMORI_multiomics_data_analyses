

###################################################################################################
###      Clonality analyses of CNA in patients with multi-timepoint, multi-region samples       ###
###################################################################################################
library(dplyr)

# prepare output dirs:
fig_dir = "~/analysis/WES/CNA/plots"
dataset_dir = "~/analysis/WES/CNA/datasets/CNA_clonality_files"
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(dataset_dir, showWarnings=FALSE, recursive=TRUE)

selected_samples_df <- read.table("/analysis/General_files/MEMORI_specific/WES_samples_include_for_analysis.txt", sep = "\t", header = T)
patient_for_CNA_clon <- selected_samples_df[which(selected_samples_df$PASS_Sequenza_best_fitting_ploidy == "y"),]
samples <- patient_for_CNA_clon$SampleID
patients <- unique(patient_for_CNA_clon$patients)

# Rename all binned 100kb files with SampleID in colnames and join Minor allele and Major allele CNS to an overview CNS
list_CNA_files <- list()
for (sam in samples){
  CNA_df <- read.table(paste0("/analysis/WES/data_sets/Sequenza/best_fitting_ploidy_binned_100kb/", sam, ".100kb_chrom_bin.txt"), header = T, stringsAsFactors = F)
  CNA_df$CNA <- paste0(CNA_df$A, "_", CNA_df$B)
  colnames(CNA_df) = c('chr', 'from', 'to', 'CNt', 'Major', 'minor', "CNA")
  colnames(CNA_df)[5:7] = paste0(sam, "_", colnames(CNA_df)[5:7])
  assign(paste0(sam, "_100kb_CNA"), CNA_df)
  list_CNA_files[[sam]] <- CNA_df
  
}


# Create a table with CNS of samples per patient 
bindf <- list_CNA_files[[1]] 
bindf <- bindf[,c('chr','from','to')]


for (pn in patients){
  df_patient <- selected_samples_df[which(selected_samples_df$patients == pn & selected_samples_df$PASS_Sequenza_best_fitting_ploidy == "y"),]
  samples_for_CNA_analysis <- df_patient$SampleID
  Patient_CNA_MasterTable <- bindf
  
  for(sam in samples_for_CNA_analysis){
    dfsam <- list_CNA_files[[sam]]
    Patient_CNA_MasterTable[[sam]] <- dfsam[,paste0(sam,'_CNA')] #
    write.table(Patient_CNA_MasterTable, file = paste0(dataset_dir, "/CNA_100kb_per_patient/", pn, ".100kb_chrom_bin.txt"), sep = "\t", quote = F, row.names = F)
  }
}

#############################
## Assign clonality of CNA ##
#############################


### RULES FOR CLONALITY DEFINITIONS ######
## 1.) for 2 samples per patient
# Clonal CNA: present in all samples (length(numdiff)==1)
# private/subclonal CNA: present in 1 of 2 samples (length(numdiff)==2)

## 2.) for =>3 samples per patient
## Clonal CNA: present in all samples (length(numdiff)==1)
## subclonal CNA: present in > 1/3 of samples, but not all samples
## private CNA: present in <= 1/3 of samples


require(dplyr)

Two_sample_patients_df <- selected_samples_df[which(selected_samples_df$PASS_Sequenza_best_fitting_ploidy == "y" & selected_samples_df$Samples_for_CNA_clonality == 2),]
Two_sample_patients <- unique(Two_sample_patients_df$patients)


##### 1.) for 2 samples per patient

for (tp in Two_sample_patients){
  cna <- read.table(paste0(dataset_dir, "/CNA_100kb_per_patient/", tp, ".100kb_chrom_bin.txt"), header = T, stringsAsFactors = F, sep ="\t")
  stability_df <- select(cna, starts_with('MEMORI'))
  
  # Assign clonality of CNA according to RULES FOR CLONALITY
  clonality<-c();for(i in c(1:nrow(stability_df))) { numdiff <- table(as.matrix(stability_df[i,])); clonality <- c(clonality,ifelse(length(numdiff)==1,'clonal',ifelse(length(numdiff)==2,'subclonal/private', "NA")))}
  stability_df$CNA_over_time <- clonality
  
  CNA_stability_df <- as.data.frame(cbind(cna$chr, cna$from, cna$to, stability_df))
  colnames(CNA_stability_df)[colnames(CNA_stability_df) == "cna$chr"] <- "chr"
  colnames(CNA_stability_df)[colnames(CNA_stability_df) == "cna$from"] <- "from"
  colnames(CNA_stability_df)[colnames(CNA_stability_df) == "cna$to"] <- "to"
  
  
  CNA_stability_df$CNA_sequence <- apply(CNA_stability_df[, startsWith(colnames(CNA_stability_df), "MEMORI")] , 1 , paste , collapse = "_" )
  
  CNA_stability_df$CNA_changes_over_time <- if_else(CNA_stability_df$CNA_over_time == "clonal", paste0("clonal_", CNA_stability_df[,4]), 
                                                    if_else(CNA_stability_df$CNA_over_time == "subclonal/private", paste0("subclonal/private_", CNA_stability_df$CNA_sequence), ""))
  
  # exclude clonal diploid regions and regions with NA Copy status
  CNA_stability_df <- CNA_stability_df[-grep('NA',CNA_stability_df$CNA_changes_over_time), ]
  CNA_stability_df <- CNA_stability_df[which(CNA_stability_df$CNA_changes_over_time != "clonal_1_1"), ]
  
  ## Calculate the fragment size of clonal and subclonal/private fragments with copy number alterations
  y <- rle(CNA_stability_df$CNA_changes_over_time)
  yy <- as.data.frame(cbind(y$values, y$lengths))
  colnames(yy) <- c("CNA_states", "number_of_100kb_bins")
  
  yy$CNA_clonality <- if_else(startsWith(yy$CNA_states, "clonal"), "clonal",  
                              if_else(startsWith(yy$CNA_states, "subclonal/private"), "subclonal/private", ""))
  yy$number_of_100kb_bins <- as.numeric(yy$number_of_100kb_bins)
  yy$Patient <- tp
  
  ## Export CNA clonality table for single patients 
  write.table(yy, file = paste0(dataset_dir, "/", tp, "_CNA_clonality.txt"), sep = "\t", quote = F, row.names = F)
  
}



#### 2.) for >=3 samples per patients 

Multisample_patients_df <- selected_samples_df[which(selected_samples_df$PASS_Sequenza_best_fitting_ploidy == "y" & selected_samples_df$Samples_for_CNA_clonality == ">=3"),]
Multisample_patients <- unique(Multisample_patients_df$patients)


for (mp in Multisample_patients){
  cna <- read.table(paste0(dataset_dir, "/CNA_100kb_per_patient/", mp, ".100kb_chrom_bin.txt"), header = T, stringsAsFactors = F, sep ="\t")
  stability_df <- select(cna, starts_with('MEMORI'))
  
  # Assign clonality of CNA according to RULES FOR CLONALITY
  x <- ncol(select(cna, starts_with('MEMORI')))
  clonality<-c();for(i in c(1:nrow(stability_df))) {numdiff <- table(as.matrix(stability_df[i,])); clonality <- c(clonality,ifelse(max(numdiff) == x,'clonal',ifelse(max(numdiff) >(1/3*x) & max(numdiff) < x,'subclonal', ifelse(max(numdiff) <=(1/3*x) ,'private', "NA"))))}
  stability_df$CNA_over_time <- clonality
  
  CNA_stability_df <- as.data.frame(cbind(cna$chr, cna$from, cna$to, stability_df))
  colnames(CNA_stability_df)[colnames(CNA_stability_df) == "cna$chr"] <- "chr"
  colnames(CNA_stability_df)[colnames(CNA_stability_df) == "cna$from"] <- "from"
  colnames(CNA_stability_df)[colnames(CNA_stability_df) == "cna$to"] <- "to"
  
  CNA_stability_df$CNA_sequence <- apply(CNA_stability_df[, startsWith(colnames(CNA_stability_df), "MEMORI")] , 1 , paste , collapse = "_" )
  
  CNA_stability_df$CNA_changes_over_time <- if_else(CNA_stability_df$CNA_over_time == "clonal", paste0("clonal_", CNA_stability_df[,4]), 
                                                    if_else(CNA_stability_df$CNA_over_time == "subclonal", paste0("subclonal_", CNA_stability_df$CNA_sequence), 
                                                            if_else(CNA_stability_df$CNA_over_time == "private", paste0("private_", CNA_stability_df$CNA_sequence), "")))
  
  # exclude clonal diploid regions and regions with NA Copy status
  CNA_stability_df <- CNA_stability_df[-grep('NA',CNA_stability_df$CNA_changes_over_time), ]
  CNA_stability_df <- CNA_stability_df[which(CNA_stability_df$CNA_changes_over_time != "clonal_1_1"), ]
  
  ## Calculate the fragment size of clonal, subclonal and private fragments
  y <- rle(CNA_stability_df$CNA_changes_over_time)
  yy <- as.data.frame(cbind(y$values, y$lengths))
  colnames(yy) <- c("CNA_states", "number_of_100kb_bins")
  
  yy$CNA_clonality <- if_else(startsWith(yy$CNA_states, "clonal"), "clonal",  
                              if_else(startsWith(yy$CNA_states, "subclonal"), "subclonal", 
                                      if_else(startsWith(yy$CNA_states, "private"), "private", "")))
  
  yy$number_of_100kb_bins <- as.numeric(yy$number_of_100kb_bins)
  yy$Patient <- mp
  ## Export CNA clonality table for individual patients 
  write.table(yy, file = paste0(dataset_dir, "/", mp, "_CNA_clonality.txt"), quote = F, row.names = F,  sep = "\t")
  
}



#############################
## Plot Results 
#############################
library(ggplot2)
library(ggpubr)
library(rstatix)

my_theme <- theme(
  plot.title = element_text(color="black", size=20, face="bold"),
  axis.title.x = element_text(color="black", size=20, face="bold"),
  axis.text.x = element_text(color="black", size=26, face="bold"),
  axis.title.y = element_text(color="black", size=18, face="bold"),
  axis.text.y = element_text(color="black", size=20, face="bold"),
  legend.title = element_text(color = "black", size = 20, face="bold"),
  legend.text = element_text(color = "black", size = 18, face="bold"),
  panel.border = element_rect(colour = "black", fill=NA, size=1.5)
  
)

clonality_files <- list.files(path = dataset_dir, pattern = "*_CNA_clonality.txt", full.names = T)

for (fn in clonality_files){
  CNA_MasterTable <- do.call(rbind,lapply(clonality_files, function(fn)read.table(fn, header = TRUE, stringsAsFactors = F, sep ="\t")))
}
CNA_MasterTable$CNA_clonality_ordered <- factor(CNA_MasterTable$CNA_clonality, levels = c("clonal", "subclonal", "private", "subclonal/private"))
CNA_MasterTable$log_number_of_100kb_bins <- log2(CNA_MasterTable$number_of_100kb_bins)


# Boxplot for fragment sizes of clonal, subclonal and private CNA in REP and NRP
CNA_MasterTable$Responsiveness <- sapply(CNA_MasterTable$Patient, function(x) substr(x,1,2))

pdf(paste0(fig_dir, "/CNA_clonality_fragmentsize.pdf"), width = 7, height = 5)
p <- ggplot(data = CNA_MasterTable, aes(x = CNA_clonality_ordered, y = log_number_of_100kb_bins, fill=Responsiveness)) + ylim(c(0,13)) +
  geom_boxplot() + scale_fill_manual(labels = c("NRP", "REP"), values = c("#b9adde","#26A69A")) + 
  xlab("") + ylab("Fragment size [log2(100kb)]") + scale_x_discrete(labels= c("clonal", "subclonal", "private", "subcl./priv.")) 
q <- p + theme_bw() + my_theme + theme(legend.title = element_blank()) + theme(axis.text.x = element_text(color="black", vjust = -0.00004, size=18, face="bold"))
plot(q)
dev.off()

# Statistics
stat.test2_between_clonalities <- CNA_MasterTable %>%
  t_test(log_number_of_100kb_bins ~ CNA_clonality_ordered) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")

stat.test2_clonal <- CNA_MasterTable[which(CNA_MasterTable$CNA_clonality == "clonal"), ] %>%
  t_test(log_number_of_100kb_bins ~ Responsiveness) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")

stat.test2_subclonal <- CNA_MasterTable[which(CNA_MasterTable$CNA_clonality == "subclonal"), ] %>%
  t_test(log_number_of_100kb_bins ~ Responsiveness) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")


stat.test2_private <- CNA_MasterTable[which(CNA_MasterTable$CNA_clonality == "private"), ] %>%
  t_test(log_number_of_100kb_bins ~ Responsiveness) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")


stat.test2_subclonal_private <- CNA_MasterTable[which(CNA_MasterTable$CNA_clonality == "subclonal/private"), ] %>%
  t_test(log_number_of_100kb_bins ~ Responsiveness) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")
