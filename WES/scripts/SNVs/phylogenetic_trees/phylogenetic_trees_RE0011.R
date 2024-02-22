
##############################################################################
# Phylogenetic trees for patients with 3 or more samples

library(dplyr)
library(phangorn)
library(stringr)

# prepare output dirs:
fig_dir = "~/analysis/WES/mutations/plots/phylogenetic_trees"
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)
data_dir = "~/analysis/WES/mutations/phylogenetic_analyses/phyl_tree_branch_tables"
dir.create(data_dir, showWarnings=FALSE, recursive=TRUE)

# Select patient and respective samples to run phylogenetic analyses
patient_name <- "MEMORI_RE0011"

sample_A1 <- "MEMORI_RE0011A_T1_B1_D1"  
sample_A2 <-"MEMORI_RE0011A_T3_B1_D1"  

sample_B1 <- "MEMORI_RE0011B_T1_C1_D1"
sample_B2 <-"MEMORI_RE0011B_T2_C1_D1" 
sample_B3 <-"MEMORI_RE0011B_T3_C1_D1" 


# Create labels for phylogenetic tree
label_A1 <- paste0(sapply(sample_A1, function(x) substr(x,14,14)), "_", sapply(sample_A1, function(x) substr(x,17,17)))
label_A2 <- paste0(sapply(sample_A2, function(x) substr(x,14,14)), "_", sapply(sample_A2, function(x) substr(x,17,17)))

label_B1 <- paste0(sapply(sample_B1, function(x) substr(x,14,14)), "_", sapply(sample_B1, function(x) substr(x,17,17)))
label_B2 <- paste0(sapply(sample_B2, function(x) substr(x,14,14)), "_", sapply(sample_B2, function(x) substr(x,17,17)))
label_B3 <- paste0(sapply(sample_B3, function(x) substr(x,14,14)), "_", sapply(sample_B3, function(x) substr(x,17,17)))

patient_label <- paste0(sapply(patient_name, function(x) substr(x,8,9)), sapply(patient_name, function(x) substr(x,12,13)))

# Read in SNV file for patient from Mobster run for multi-sample analyses 
Multi_mobster_output_nontail <- read.table(paste0("~/analysis/WES/mutations/clonality_analyses/multi_sample_MOBSTER/output_tables/", patient_name, ".driver_cluster_annotated_nontail_Mobster_output.txt"), header = T, sep = "\t", stringsAsFactors = F)

# Merge NeoSNV files to add information on NeoSNV counts to phylogenetic trees
neo_snv_1 <- read.table(paste0("~/analysis/WES/neoantigens/filtering_steps/tables/neoSNV_tables/", sample_A1, ".weakly_filtered_neoSNVs"), header = T, stringsAsFactors = F, sep="\t")
neo_snv_2 <- read.table(paste0("~/analysis/WES/neoantigens/filtering_steps/tables/neoSNV_tables/", sample_A2, ".weakly_filtered_neoSNVs"), header = T, stringsAsFactors = F, sep="\t")
neo_snv_3 <- read.table(paste0("~/analysis/WES/neoantigens/filtering_steps/tables/neoSNV_tables/", sample_B1, ".weakly_filtered_neoSNVs"), header = T, stringsAsFactors = F, sep="\t")
neo_snv_4 <- read.table(paste0("~/analysis/WES/neoantigens/filtering_steps/tables/neoSNV_tables/", sample_B2, ".weakly_filtered_neoSNVs"), header = T, stringsAsFactors = F, sep="\t")
neo_snv_5 <- read.table(paste0("~/analysis/WES/neoantigens/filtering_steps/tables/neoSNV_tables/", sample_B3, ".weakly_filtered_neoSNVs"), header = T, stringsAsFactors = F, sep="\t")

final_Neo_SNV <- data.frame(rbind(neo_snv_1, neo_snv_2, neo_snv_3, neo_snv_4, neo_snv_5))
final_Neo_SNV <- final_Neo_SNV[,c("SNV", "BindLevel")]
colnames(final_Neo_SNV)[colnames(final_Neo_SNV) == "BindLevel"]  <- "NeoSNV_Bind_level"
final_Neo_SNV <- final_Neo_SNV[!duplicated(final_Neo_SNV$SNV), ]
Multi_mobster_output_nontail2 <- merge(Multi_mobster_output_nontail, final_Neo_SNV, by.x = "non_tail_SNV", by.y = "SNV", all.x = T)

# Reorder table
selected_df = Multi_mobster_output_nontail2 %>%
  select(starts_with('MEMORI'))
gene_annot = Multi_mobster_output_nontail2 %>%
  select(!starts_with('MEMORI'))
Multi_mobster_output_nontail2 <- data.frame(cbind(gene_annot, selected_df))

# Make SNV table binary (VAF> or =0.005 => SNV present; VAF<0.005 => SNV not present)
binary_df <-  ifelse(selected_df<0.005, 0, 1); binary_df <- as.data.frame(binary_df)
Multi_mobster_output_nontail_binary <- as.data.frame(cbind(Multi_mobster_output_nontail2[,c(1:6)], binary_df))

DNAbin <- binary_df  %>%
  select(starts_with('MEMORI'))
colnames(DNAbin) <- paste0(sapply(colnames(DNAbin), function(x) substr(x,14,14)), "_", sapply(colnames(DNAbin), function(x) substr(x,17,17)))

DNAbin <- cbind(DNAbin, Root = 0)

# Run phylogenetic (Maximum Parsimony) analysis
DNAphy <- as.phyDat(DNAbin,type="USER",levels=c(0,1))
DNAtreeMP <- pratchet(DNAphy);DNAtreeMP <- root(DNAtreeMP,'Root',resolve.root=T)
DNAtreeMP <- acctran(DNAtreeMP, DNAphy);set.seed(123)
DNABStrees <- bootstrap.phyDat(DNAphy, pratchet, bs = 100) # Bootstrapping

# Calculate Consistency Index and Homoplasy Index
ci <- CI(DNAtreeMP,DNAphy) # ci=1 means no homoplasy
hi <- signif(1-ci,digits=2) # so hi=0 means no homoplasy

# Table with information on driver SNVs
DNA_bin_driver_annotated <- cbind(Multi_mobster_output_nontail_binary$non_tail_SNV, Multi_mobster_output_nontail_binary$driver_label_Intogen, Multi_mobster_output_nontail_binary$driver_label_COSMIC_INTOGEN, Multi_mobster_output_nontail_binary$Cluster, Multi_mobster_output_nontail_binary$Info_SNV_synonymous, Multi_mobster_output_nontail_binary$NeoSNV_Bind_level, DNAbin)
colnames(DNA_bin_driver_annotated)[colnames(DNA_bin_driver_annotated) == "Multi_mobster_output_nontail_binary$non_tail_SNV"] <- "non_tail_SNV"
colnames(DNA_bin_driver_annotated)[colnames(DNA_bin_driver_annotated) == "Multi_mobster_output_nontail_binary$Cluster"] <- "Cluster"
colnames(DNA_bin_driver_annotated)[colnames(DNA_bin_driver_annotated) == "Multi_mobster_output_nontail_binary$driver_label_Intogen"] <- "driver_label_Intogen"
colnames(DNA_bin_driver_annotated)[colnames(DNA_bin_driver_annotated) == "Multi_mobster_output_nontail_binary$driver_label_COSMIC_INTOGEN"] <- "driver_label_COSMIC_INTOGEN"
colnames(DNA_bin_driver_annotated)[colnames(DNA_bin_driver_annotated) == "Multi_mobster_output_nontail_binary$Info_SNV_synonymous"] <- "Info_SNV_synonymous"


DNA_bin_driver_annotated$non_synon_driver_Intogen <- if_else(!(DNA_bin_driver_annotated$Info_SNV_synonymous == "synonymous SNV" | DNA_bin_driver_annotated$Info_SNV_synonymous == "unknown" | DNA_bin_driver_annotated$Info_SNV_synonymous == "."), DNA_bin_driver_annotated$driver_label_Intogen, "")
DNA_bin_driver_annotated$non_synon_driver_COSMIC_INTOGEN <- if_else(!(DNA_bin_driver_annotated$Info_SNV_synonymous == "synonymous SNV" | DNA_bin_driver_annotated$Info_SNV_synonymous == "unknown" | DNA_bin_driver_annotated$Info_SNV_synonymous == "."), DNA_bin_driver_annotated$driver_label_COSMIC_INTOGEN, "")

# List of shared and private SNVs between timepoints. Non-synonymous driver SNVs and strong or weak binding NeoSNVs are annotated in the table and manually added to the phylogenetic tree.
shared_all <- DNA_bin_driver_annotated[which(DNA_bin_driver_annotated[,label_A1] == 1 & DNA_bin_driver_annotated[,label_A2] == 1 & DNA_bin_driver_annotated[,label_B1] == 1 & DNA_bin_driver_annotated[,label_B2] == 1 & DNA_bin_driver_annotated[,label_B3] == 1), ]
present_A1 <- DNA_bin_driver_annotated[which(DNA_bin_driver_annotated[,label_A1] == 1),]
present_A2 <- DNA_bin_driver_annotated[which(DNA_bin_driver_annotated[,label_A2] == 1),]
not_present_A1_A2 <- DNA_bin_driver_annotated[which(DNA_bin_driver_annotated[,label_A1] == 0 & DNA_bin_driver_annotated[,label_A2] == 0),]

shared_all_not_B3 <- DNA_bin_driver_annotated[which(DNA_bin_driver_annotated[,label_A1] == 1 & DNA_bin_driver_annotated[,label_A2] == 1 & DNA_bin_driver_annotated[,label_B1] == 1 & DNA_bin_driver_annotated[,label_B2] == 1 & DNA_bin_driver_annotated[,label_B3] == 0), ] 

shared_A1_B1_B2 <- DNA_bin_driver_annotated[which(DNA_bin_driver_annotated[,label_A1] == 1 & DNA_bin_driver_annotated[,label_A2] == 0 & DNA_bin_driver_annotated[,label_B1] == 1 & DNA_bin_driver_annotated[,label_B2] == 1 & DNA_bin_driver_annotated[,label_B3] == 0), ] 
shared_A1_B1 <- DNA_bin_driver_annotated[which(DNA_bin_driver_annotated[,label_A1] == 1 & DNA_bin_driver_annotated[,label_A2] == 0 & DNA_bin_driver_annotated[,label_B1] == 1 & DNA_bin_driver_annotated[,label_B2] == 0 & DNA_bin_driver_annotated[,label_B3] == 0), ] 
shared_A1_A2 <- DNA_bin_driver_annotated[which(DNA_bin_driver_annotated[,label_A1] == 1 & DNA_bin_driver_annotated[,label_A2] == 1 & DNA_bin_driver_annotated[,label_B1] == 0 & DNA_bin_driver_annotated[,label_B2] == 0 & DNA_bin_driver_annotated[,label_B3] == 0), ] 

private_A1 <- DNA_bin_driver_annotated[which(DNA_bin_driver_annotated[,label_A1] == 1 & DNA_bin_driver_annotated[,label_A2] == 0 & DNA_bin_driver_annotated[,label_B1] == 0 & DNA_bin_driver_annotated[,label_B2] == 0 & DNA_bin_driver_annotated[,label_B3] == 0), ] 
private_A2 <- DNA_bin_driver_annotated[which(DNA_bin_driver_annotated[,label_A1] == 0 & DNA_bin_driver_annotated[,label_A2] == 1 & DNA_bin_driver_annotated[,label_B1] == 0 & DNA_bin_driver_annotated[,label_B2] == 0 & DNA_bin_driver_annotated[,label_B3] == 0), ] 

private_B1 <- DNA_bin_driver_annotated[which(DNA_bin_driver_annotated[,label_A1] == 0 & DNA_bin_driver_annotated[,label_A2] == 0 & DNA_bin_driver_annotated[,label_B1] == 1 & DNA_bin_driver_annotated[,label_B2] == 0 & DNA_bin_driver_annotated[,label_B3] == 0), ] 
private_B2 <- DNA_bin_driver_annotated[which(DNA_bin_driver_annotated[,label_A1] == 0 & DNA_bin_driver_annotated[,label_A2] == 0 & DNA_bin_driver_annotated[,label_B1] == 0 & DNA_bin_driver_annotated[,label_B2] == 1 & DNA_bin_driver_annotated[,label_B3] == 0), ] 
private_B3 <- DNA_bin_driver_annotated[which(DNA_bin_driver_annotated[,label_A1] == 0 & DNA_bin_driver_annotated[,label_A2] == 0 & DNA_bin_driver_annotated[,label_B1] == 0 & DNA_bin_driver_annotated[,label_B2] == 0 & DNA_bin_driver_annotated[,label_B3] == 1), ] 

View(shared_all); View(shared_all_not_B3); View(shared_A1_B1_B2); View(shared_A1_B1);
View(private_A1); View(private_A2); View(private_B1); View(private_B2); View(private_B3)


# Export SNVs from tree branches according to received therapy 
pre_treatment <- rbind(present_A1, present_A2)
unique_pre_treat_SNV<- unique(pre_treatment$non_tail_SNV)
Other_info_list <- str_split(unique_pre_treat_SNV,"_")
extractinfo <- function(lst,n) { sapply(lst,'[', n)}
pre_treat_df <- data.frame(sample = patient_name, chr = extractinfo(Other_info_list, 1), from = extractinfo(Other_info_list, 2), ref = extractinfo(Other_info_list, 3), to = extractinfo(Other_info_list, 4))
write.table(pre_treat_df, file = paste0(data_dir,"/pre_treatment_SNVs/", patient_name, ".pre_treatment.txt"), sep = "\t", quote = F, row.names = F)

chemo_induced <- not_present_A1_A2
Other_info_list <- str_split(chemo_induced$non_tail_SNV,"_")
extractinfo <- function(lst,n) { sapply(lst,'[', n)}
chemo_induc_df <- data.frame(sample = patient_name, chr = extractinfo(Other_info_list, 1), from = extractinfo(Other_info_list, 2), ref = extractinfo(Other_info_list, 3), to = extractinfo(Other_info_list, 4))
write.table(chemo_induc_df, file = paste0(data_dir,"/chemo_induced_SNVs/", patient_name, ".chemo_induced_SNVs.txt"), sep = "\t", quote = F, row.names = F)


# Plot phylogeny - comes up with error saying lots of the options aren't graphical parameters (they are...)
# Non-synonymous driver SNVs and NeoSNV counts are manually added to the plot
pdf(paste0(fig_dir, "/", patient_name, ".phyl.tree.pdf"), height = 4.8, width = 8)
plotBS(DNAtreeMP,DNABStrees,type.tree='phylogram',bs.adj=c(1.3,1.3),label.offset=15,p=0,type="phylogram", edge.color= c('#0A9086', '#E64B35FF', 'grey', '#0A9086', 'grey', '#E64B35FF', 'grey', '#0A9086',"#7E6148FF", "#7E6148FF"), edge.width=3,font=2,cex=1.5,tip.color=c('#E64B35FF', '#E64B35FF', '#0A9086', '#0A9086', '#0A9086', "#7E6148FF"))
add.scale.bar(0,01, lwd=3,font=2);mtext(paste0('HI: ',hi),side=3,font=2,cex=1.25,ask=T);show.node.label = T;
title(main = patient_label, cex.main=1.5, col.main= "black")
dev.off()

