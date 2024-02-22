

#################################################################
# Calling COSMIC mutation signatures with SigDeconstruct
#################################################################
# This script
# 1. runs SigDeconstruct for each sample including all COSMIC signatures
# 2. runs SigDeconstruct for each sample limiting called COSMIC signatures to signatures that are present at >5% in NRP or REP at any timepoint
# 3. runs SigDeconstruct for treatment-naive SNVs, chemo-induced SNVs and radiochemo induced SNVs

library(deconstructSigs)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyr)

# prepare output dirs:
data_dir1 = "~/analysis/WES/mutations/COSMIC_signatures/output_tables/Per_sample/signature_weight/all_signatures"
dir.create(data_dir1, showWarnings=FALSE, recursive=TRUE)
data_dir2 = "~/analysis/WES/mutations/COSMIC_signatures/output_tables/Per_sample/trinucleotide_output/"
dir.create(data_dir2, showWarnings=FALSE, recursive=TRUE)
data_dir3 = "~/analysis/WES/mutations/COSMIC_signatures/output_tables/Per_sample/base_output/"
dir.create(data_dir3, showWarnings=FALSE, recursive=TRUE)
data_dir4 = "~/analysis/WES/mutations/COSMIC_signatures/output_tables/Per_sample/signature_weight/more_5percent_signatures"
dir.create(data_dir4, showWarnings=FALSE, recursive=TRUE)

data_dir5 = "~/analysis/WES/mutations/COSMIC_signatures/output_tables/Per_SNV_inducer/signature_weight/"
dir.create(data_dir5, showWarnings=FALSE, recursive=TRUE)
data_dir6 = "~/analysis/WES/mutations/COSMIC_signatures/output_tables/Per_SNV_inducer/trinucleotide_output/"
dir.create(data_dir6, showWarnings=FALSE, recursive=TRUE)
data_dir7 = "~/analysis/WES/mutations/COSMIC_signatures/output_tables/Per_SNV_inducer/base_output/"
dir.create(data_dir7, showWarnings=FALSE, recursive=TRUE)


# Select samples included for genetic analyses
selected_samples <- read.table("~/analysis/General_files/MEMORI_specific/WES_samples_for_analysis.txt", sep = "\t", header = T)
selected_samples_Seq <- selected_samples[which(selected_samples$PASS_Sequenza_best_fitting_ploidy == "y"), ]
samples <- selected_samples_Seq$SampleID

# Create mastertable with filtered SNVs as input for SigDeconstruct
list_weak_filtered_SNV_files <- list()
for (sam in samples){
  df <- read.table(paste0("~/analysis/WES/input_files/mutation_results/mutations_filtered/", sam, ".weak.filtered.SNVs.txt"), header = TRUE,stringsAsFactors = F)
  selected_chr <- paste0("chr", c(1:22, "X", "Y"))
  df$Sample <- sam 
  
  ## filter for strange chromosomes
  norm_chrom <-paste0("chr", c(1:22, "X", "X"))
  df <- df[which(df$Chr %in% norm_chrom),]
  
  df <- df[, c("Sample", "Chr", "Start", "Ref", "Alt")]
  colnames(df) <-  c("Sample", "chr", "pos", "ref", "alt")
  list_weak_filtered_SNV_files[[sam]] <- df
}

input_deconsigs <- rbindlist(list_weak_filtered_SNV_files) 

 

# 1. Run SigDeconstruct for each sample including all COSMIC signatures
 sigs.input_per_sample <- mut.to.sigs.input(mut.ref = input_deconsigs, 
                                 sample.id = "Sample", 
                                 chr = "chr", 
                                 pos = "pos", 
                                 ref = "ref", 
                                 alt = "alt", 
                                 bsg=BSgenome.Hsapiens.UCSC.hg38)
 
 
 for (sam in samples) {
   
 y  = whichSignatures(tumor.ref = sigs.input_per_sample, 
                          signatures.ref = signatures.cosmic, 
                          sample.id = sam, 
                          contexts.needed = TRUE,
                          tri.counts.method = 'default')
 
 # export table with COSMIC weight signatures
 sig_table_1 <- as.data.frame(t(y$weights))
 sig_table_1$COSMIC_signatures <- rownames(sig_table_1)
 sig_table_1 <- sig_table_1[,c(2,1)]
 write.table(sig_table_1, file = paste0(data_dir1, "/", sam, "_signature_weights.txt"), row.names = F, sep = "\t", quote = F)

 # export table with base changes in different trinucleotides
 sig_table_2 <- as.data.frame(t(y$tumor))
 sig_table_2$Signature_bases <- rownames(sig_table_2)
 sig_table_2 <- sig_table_2[,c(2,1)]
 write.table(sig_table_2, file = paste0(data_dir2, "/", sam, "_signature_tumours.txt"), row.names = F, sep = "\t", quote = F)
                                                                                                                                                                   mutation_signatures/Tables/Per_sample/weight_output/non_mobster2/
 # summarize base changes in different trinucleotides and export table
 sig_table_3 <- sig_table_2
 sig_table_3$presummarize_Signature_bases <- sapply(sig_table_3$Signature_bases, function(x) substr(x,3, 5))
 base_changes <- unique(sig_table_3$presummarize_Signature_bases)
 
 for (b in base_changes){ 
 df <- sig_table_3[which(sig_table_3$presummarize_Signature_bases == b), ]
 
 sums_df <- data.frame(t(colSums(df[2]))); rownames(sums_df) <- paste0("base_change_", b)
 assign(paste0("base_change_", b), sums_df)
 }
 
 summarized_base_change_df <- rbind(`base_change_C>A`, `base_change_C>G`, `base_change_C>T`, `base_change_T>A`, `base_change_T>C`, `base_change_T>G`)
 summarized_base_change_df$Base_change_signature <- rownames(summarized_base_change_df)
 summarized_base_change_df <- summarized_base_change_df[,c(2,1)]
 write.table(summarized_base_change_df,file = paste0(data_dir3, "/", sam, "_signature_tumours_summarized.txt"), row.names = F, sep = "\t", quote = F)

 
 }
 
 ## 2. Run SigDeconstruct  for each sample limiting called signatures to signatures that are present at >5% in NRP or REP at any timepoint
signatures_more5percent <- read.table("~/analysis/WES/mutations/COSMIC_signatures/output_tables/Per_sample/signature_weight/More_5_percent_in_Responsiveness_Timepoint_groups.txt", sep = "\t", header = T)
signatures_more5percent <- signatures_more5percent[which(signatures_more5percent$x != "unknown"), ]
seletced_sig <- as.integer(stringr::str_extract(signatures_more5percent, "\\d+"))
seletced_sig <- paste0("Signature.", seletced_sig)

 for (sam in samples) {
   
   y  = whichSignatures(tumor.ref = sigs.input_per_sample, 
                        signatures.ref = signatures.cosmic, 
                        sample.id = sam, 
                        associated = seletced_sig,
                        contexts.needed = TRUE,
                        tri.counts.method = 'default')
   
   # export table on base COSMIC weight sigantures
   sig_table_1 <- as.data.frame(t(y$weights))
   sig_table_1$COSMIC_signatures <- rownames(sig_table_1)
   sig_table_1 <- sig_table_1[,c(2,1)]
   write.table(sig_table_1, file = paste0(data_dir4, "/", sam, "_signature_weights_more5percent_resp_TP_group.txt"), row.names = F, sep = "\t", quote = F)
    }

 
 
 # 3. Run SigDeconstruct for treatment-naive SNVs, SNVs newly occuring under chemotherapy and SNVs newly occuring under radiochemotherapy

 input_per_patient_treatment_naive_RE <- list.files("~/analysis/WES/phylogenetic_analyses/phyl_tree_branch_tables/pre_treatment_SNVs/",  pattern = "MEMORI_RE*", full.names = T)
 input_per_patient_treatment_naive_NE <- list.files("~/analysis/WES/phylogenetic_analyses/phyl_tree_branch_tables/pre_treatment_SNVs/",  pattern = "MEMORI_NE*", full.names = T)
 input_per_patient_chemo_induced <- list.files("~/analysis/WES/phylogenetic_analyses/phyl_tree_branch_tables/pre_treatment_SNVs/chemo_induced_SNVs/",  pattern = "*.chemo_induced_SNVs.txt", full.names = T)
 input_per_patient_radio_induced <- list.files("~/analysis/WES/phylogenetic_analyses/phyl_tree_branch_tables/pre_treatment_SNVs/radiochemo_induced_SNVs/",  pattern = "*.radio_induced_SNVs.txt", full.names = T)
 
 for (x in input_per_patient_treatment_naive_RE){
   input_deconsigs_treat_naive_RE <- do.call(rbind,lapply(input_per_patient_treatment_naive_RE, function(x)read.table(x, header = T ,stringsAsFactors = F)))}
    input_deconsigs_treat_naive_RE$SNV_origin <- "treatment_naive_RE"
 
  for (x in input_per_patient_treatment_naive_NE){
      input_deconsigs_treat_naive_NE <- do.call(rbind,lapply(input_per_patient_treatment_naive_NE, function(x)read.table(x, header = T ,stringsAsFactors = F)))}
    input_deconsigs_treat_naive_NE$SNV_origin <- "treatment_naive_NE"
    
  for (x in input_per_patient_chemo_induced){
     input_deconsigs_chemo_ind <- do.call(rbind,lapply(input_per_patient_chemo_induced, function(x)read.table(x, header = T ,stringsAsFactors = F)))}
      input_deconsigs_chemo_ind$SNV_origin <- "chemo_induced"
   
 for (x in input_per_patient_radio_induced){
   input_deconsigs_radio_ind <- do.call(rbind,lapply(input_per_patient_radio_induced, function(x)read.table(x, header = T ,stringsAsFactors = F)))}
    input_deconsigs_radio_ind$SNV_origin <- "radio_induced"
    
    input_deconsigs <- rbind(input_deconsigs_treat_naive_RE, input_deconsigs_treat_naive_NE, input_deconsigs_chemo_ind, input_deconsigs_radio_ind)
    colnames(input_deconsigs) <- c("Sample", "chr", "pos", "ref", "alt", "SNV_origin")
 
    input_deconsigs$SNV_origin_3_cat <- ifelse(input_deconsigs$SNV_origin == "treatment_naive_NE" | input_deconsigs$SNV_origin == "treatment_naive_RE", "treatment_naive", input_deconsigs$SNV_origin)
 
    # Run SigDeconstruct for treatment naive, CTx induced and RCTx induced SNVs
    sigs.input <- mut.to.sigs.input(mut.ref = input_deconsigs, 
                                    sample.id = "SNV_origin_3_cat", 
                                    chr = "chr", 
                                    pos = "pos", 
                                    ref = "ref", 
                                    alt = "alt", 
                                    bsg=BSgenome.Hsapiens.UCSC.hg38)
    
    SNV_inducer <- c("treatment_naive", "chemo_induced", "radio_induced")
      
    
    for (fn in SNV_inducer){ 
    
    y = whichSignatures(tumor.ref = sigs.input, 
                           signatures.ref = signatures.cosmic, 
                           sample.id = fn, 
                           contexts.needed = TRUE,
                           tri.counts.method = 'default')
    
    
    # export table with COSMIC weight signatures
    sig_table_1 <- as.data.frame(t(y$weights))
    sig_table_1$COSMIC_signatures <- rownames(sig_table_1)
    sig_table_1 <- sig_table_1[,c(2,1)]
    write.table(sig_table_1, file = paste0(data_dir5, "/", fn, "_SNV_signature_weights.txt"), row.names = F, sep = "\t", quote = F)
    
    # export table with base changes in different trinucleotides
    sig_table_2 <- as.data.frame(t(y$tumor))
    sig_table_2$Signature_bases <- rownames(sig_table_2)
    sig_table_2 <- sig_table_2[,c(2,1)]
    write.table(sig_table_2, file = paste0(data_dir6, "/", fn, "_SNV_signature_tumours.txt"), row.names = F, sep = "\t", quote = F)
    
    # summarize same base change in different trinucleotides and export table
    sig_table_3 <- sig_table_2
    sig_table_3$presummarize_Signature_bases <- sapply(sig_table_3$Signature_bases, function(x) substr(x,3, 5))
    base_changes <- unique(sig_table_3$presummarize_Signature_bases)
    
    for (b in base_changes){ 
      df <- sig_table_3[which(sig_table_3$presummarize_Signature_bases == b), ]
      
      sums_df <- data.frame(t(colSums(df[2]))); rownames(sums_df) <- paste0("base_change_", b)
      assign(paste0("base_change_", b), sums_df)
    }
    
    summarized_base_change_df <- rbind(`base_change_C>A`, `base_change_C>G`, `base_change_C>T`, `base_change_T>A`, `base_change_T>C`, `base_change_T>G`)
    summarized_base_change_df$Base_change_signature <- rownames(summarized_base_change_df)
    summarized_base_change_df <- summarized_base_change_df[,c(2,1)]
    write.table(summarized_base_change_df,file = paste0(data_dir7, "/", fn, "_SNV_signature_tumours_summarized.txt"), row.names = F, sep = "\t", quote = F)
    
    
    }