

##################################################################
###     Filtering called SNVs from Annovar output files        ###
##################################################################
# This script filters SNV with the follwing criteria: 
#  * PASS in Annovar output file
#  * reads at SNV position in the tumour and in germline control > 5 reads
#  * > 3 alternative reads at SNV position in the tumour and < 2 alternative reads in the germline control

library(stringr)

selected_samples <- read.table("~/analysis/WES/input_files/WES_samples_for_analysis.txt", sep = "\t", header = T)
selected_samples_Seq <- selected_samples[which(selected_samples$PASS_Sequenza_best_fitting_ploidy == "y"), ]
samples <- selected_samples_Seq$SampleID


for (sam in samples){
  
  df <- read.csv(paste0("~/analysis/WES/input_files/mutation_results/Annovar_output/", sam, ".hg38_multianno.csv"), stringsAsFactors = F, header=T)
  Other_info_list <- str_split(df$Otherinfo,"\t")
  extractinfo <- function(lst,n) { sapply(lst,'[', n)}
  df$Filter <- extractinfo(Other_info_list ,10)
  
  
  ##### Filter for mutations that passed the filter criterias
  df <- df[which(df$Filter=="PASS" | df$Filter=="."),]
  
  
  ###### filter to get only the SNVs and filter out indels
  nucs <- c("A","T","C","G")
  df <- df[which(df$Ref %in% nucs & df$Alt %in% nucs),]
  
  
  ############################################################################################
  ####### filter for total reads in the tumour > 5 and alternative (variant reads) >= 3
  ############################################################################################
  
  ## AD in vcf file has 2 numbers separated by ",": The 1st is the sum of referece reads, the 2nd the sum of alternative (variant) reads
  
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
  
  ## AD in vcf file has 2 numbers separated by ",": The 1st is the sum of referece reads, the 2nd the sum of alternative (variant) reads
  
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
  
  
  #### calculate the VAF = alternative read/ (reference + alternative read) ####
  df$VAF_manually <- df$sum_alter_tumour/(df$sum_refer_tumour + df$sum_alter_tumour)

  write.table(df, file = paste0("~/analysis/WES/input_files/mutation_results/mutations_filtered/", sam, ".filtered.SNVs.txt"), row.names=F, sep="\t")
  
}

