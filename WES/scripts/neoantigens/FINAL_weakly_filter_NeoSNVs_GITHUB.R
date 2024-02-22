
##################################################################
# Filtering steps for neoantigenic SNVs 
##################################################################
# filtering criteria for neoantigenic SNVs are: 
#  * PASS in Annovar output file
#  * reads at SNV position in the tumour and in germline control > 5 reads
#  * > 3 alternative reads at SNV position in the tumour and < 2 alternative reads in the germline control

library(stringr)

# prepare output dirs:
data_dir = "~/analysis/WES/neoantigens/filtering_steps/tables"
dir.create(data_dir, showWarnings=FALSE, recursive=TRUE)

# Select samples included for genetic analyses
selected_samples <- read.table("~/analysis/WES/input_files/WES_samples_for_analysis.txt", sep = "\t", header = T)
selected_samples_Seq <- selected_samples[which(selected_samples$PASS_Sequenza_best_fitting_ploidy == "y"), ]
samples <- selected_samples_Seq$SampleID

### filter the Annovar output files including all SNVs equally as for all other SNV downstream analyses 
for (sam in samples){

# read in and process Annovar file 
  df <- read.csv(paste0("~/analysis/WES/input_files/mutation_results/Annovar_output/", sam, ".hg38_multianno.csv"), stringsAsFactors = F, header=T)
  Other_info_list <- str_split(df$Otherinfo,"\t")
  extractinfo <- function(lst,n) { sapply(lst,'[', n)}
  df$Filter <- extractinfo(Other_info_list ,10)
  
  ##### Filter for mutations that passed the Filter criterias in Annovar
  df <- df[which(df$Filter=="PASS" | df$Filter=="."),]
  
  ###### filter to get only the SNVs and filter out the bigger inserts or deletions
  nucs <- c("A","T","C","G")
  df <- df[which(df$Ref %in% nucs & df$Alt %in% nucs),]
  
  
  ############################################################################################
  ####### filter for total reads in the tumour > 5 and alternative (variant reads) > 3
  ############################################################################################
  
  ## AD in vcf file has 2 numbers separated by a ",": The 1st is the sum of referece reads, the 2nd the sum of alternative (variant) reads
  library(stringr)
  Other_info_list <- str_split(df$Otherinfo,"\t")
  extractinfo <- function(lst,n) { sapply(lst,'[', n)}
  df$vcf_tumor <-  extractinfo(Other_info_list ,13)
  
  vcf_tum_list <- str_split(df$vcf_tumor,":")
  extractinfo <- function(lst,n) { sapply(lst,'[', n)}
  df$AD_tumour <- extractinfo(vcf_tum_list,2)
  AD_tumour_split <- str_split(df$AD_tumour,",")
  
  
  df$sum_refer_tumour <- as.numeric(extractinfo(AD_tumour_split,1))
  df$sum_alter_tumour <- as.numeric(extractinfo(AD_tumour_split,2))
  df$totalsum_tumour <- df$sum_refer_tumour + df$sum_alter_tumour
  
  ### filter for total reads in tumour > 5
  bin1 <- df[(which(df$totalsum_tumour <= 5)), ]
  df <- df[(which(df$totalsum_tumour > 5)), ]
  
  ### filter for alternative reads in tumour > 3
  bin2 <- df[(which(df$sum_alter_tumour < 3)), ]
  df <- df[(which(df$sum_alter_tumour >= 3)), ]
  
  
  ############################################################################################
  ####### filter for total reads in the normal control > 5 reads and alternative (variant reads) <= 2
  ############################################################################################
  
  ## AD in vcf file has 2 numbers separated by a ",": The 1st is the sum of referece reads, the 2nd the sum of alternative (variant) reads
  Other_info_list <- str_split(df$Otherinfo,"\t")
  extractinfo <- function(lst,n) { sapply(lst,'[', n)}
  df$vcf_normal<-  extractinfo(Other_info_list ,14)
  
  vcf_norm_list <- str_split(df$vcf_normal,":")
  extractinfo <- function(lst,n) { sapply(lst,'[', n)}
  df$AD_normal <- extractinfo(vcf_norm_list,2)
  
  AD_normal_split <- str_split(df$AD_normal,",")
  
  df$sum_refer_normal <- as.numeric(extractinfo(AD_normal_split,1))
  df$sum_alter_normal <- as.numeric(extractinfo(AD_normal_split,2))
  df$totalsum_normal <- df$sum_refer_normal + df$sum_alter_normal
  
  ### filter for total reads in normal > 5 
  bin3 <- df[(which(df$totalsum_normal <= 5)), ]
  df <- df[(which(df$totalsum_normal > 5)), ]
  
  ### filter for alternative reads in normal <= 2
  bin4 <- df[(which(df$sum_alter_normal > 2)), ]
  df <- df[(which(df$sum_alter_normal <= 2)), ]
  
  df$VAF_annovar <- df$sum_alter_tumour/(df$sum_refer_tumour + df$sum_alter_tumour)
  df$SNV <- paste0(df$Chr, "_", df$Start, "_", df$Ref, "_", df$Alt)
  reduced_df <- df[,c("SNV", "VAF_annovar")]
  
  
  ## Import neo_SNV file for matching patient
  neo_SNV <- read.table(paste0("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Neoantigens/Neoepitopes/results/neopred_output_Optitype/", sam, ".neoantigens.txt"), header = T, stringsAsFactors = F)
  colnames(neo_SNV) <- c('Sample', 'LineID', 'Chr', 'Start',
                         'Ref', 'Alt', 'Gene', 'pos', 'hla', 'peptide', 'core', 'Of', 'Gp',
                         'Gl', 'Ip', 'Il', 'Icore', 'Identity', 'Score','Affinity', 'Rank', 'Cand', 'BindLevel', 'Novelty')
  neo_SNV$SNV <- paste0(neo_SNV$Chr, "_", neo_SNV$Start, "_", neo_SNV$Ref, "_", neo_SNV$Alt)
  
  
  # merge the filtered Annovar SNV file with neoSNV file and keep ONLY matching SNVs from both dataframes and neoantigenic SNVs with Novelty = 1
  filtered_neo_SNVs <- merge(reduced_df, neo_SNV, by = "SNV")
  filtered_neo_SNVs <- filtered_neo_SNVs[which(filtered_neo_SNVs$Novelty == 1), ]
  write.table(filtered_neo_SNVs, file = paste0(data_dir, "/neoSNV_tables/", sam, ".weakly_filtered_neoSNVs"), quote= F, row.names=F, sep="\t")
  
  # Create summary table with number of total, strong binder (SB) and weak binder (WB) neopeptides
  summary_neopeptides <- as.data.frame(cbind(sam, nrow(filtered_neo_SNVs), nrow(filtered_neo_SNVs[which(filtered_neo_SNVs$BindLevel == "SB"), ]), nrow(filtered_neo_SNVs[which(filtered_neo_SNVs$BindLevel == "WB"), ])))
  colnames(summary_neopeptides) <- c("SampleID", "Total", "Total_SB", "Total_WB")
  write.table(summary_neopeptides, file = paste0(data_dir, "/summary_tables_neopeptides/", sam, ".weakly_filtered_summary_neopeptides"), quote= F, row.names=F, sep="\t")
  
  
  # Create summary table with number of neoantigenic mutations 
  summary_neoSNVs <- filtered_neo_SNVs[!duplicated(filtered_neo_SNVs$SNV), ]
  summary_neoSNVs <- as.data.frame(cbind(sam, nrow(summary_neoSNVs)))
  colnames(summary_neoSNVs) <- c("SampleID", "Total_NeoSNVs")
  write.table(summary_neoSNVs,  file = paste0(data_dir, "/summary_tables_neoSNVs/", sam, ".weakly_filtered_summary_neoSNVs"), quote= F, row.names=F, sep="\t")
  
}

