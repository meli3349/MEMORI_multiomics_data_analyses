
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
patient_name <- "MEMORI_NE0009"
         
sample_A1 <- "MEMORI_NE0009A_T1_B1_D1"    
sample_B1 <- "MEMORI_NE0009B_T1_B1_D1"
sample_C1 <- "MEMORI_NE0009C_T1_C1_D1" 


# Create labels for phylogenetic tree
label_A1 <- paste0(sapply(sample_A1, function(x) substr(x,14,14)), "_", sapply(sample_A1, function(x) substr(x,17,17)))
label_B1 <- paste0(sapply(sample_B1, function(x) substr(x,14,14)), "_", sapply(sample_B1, function(x) substr(x,17,17)))
label_C1 <- paste0(sapply(sample_C1, function(x) substr(x,14,14)), "_", sapply(sample_C1, function(x) substr(x,17,17)))

patient_label <- paste0(sapply(patient_name, function(x) substr(x,8,8)), "R ", sapply(patient_name, function(x) substr(x,12,13)))

# Read in SNV file for patient from Mobster run for multi-sample analyses 

Multi_mobster_output_nontail <- read.table(paste0("~/analysis/WES/mutations/clonality_analyses/multi_sample_MOBSTER/output_tables/", patient_name, ".driver_cluster_annotated_nontail_Mobster_output.txt"), header = T, sep = "\t", stringsAsFactors = F)

# Merge NeoSNV files to add information on NeoSNV counts to phylogenetic trees
neo_snv_1 <- read.table(paste0("~/analysis/WES/neoantigens/filtering_steps/tables/neoSNV_tables/", sample_A1, ".weakly_filtered_neoSNVs"), header = T, stringsAsFactors = F, sep="\t")
neo_snv_2 <- read.table(paste0("~/analysis/WES/neoantigens/filtering_steps/tables/neoSNV_tables/", sample_B1, ".weakly_filtered_neoSNVs"), header = T, stringsAsFactors = F, sep="\t")
neo_snv_3 <- read.table(paste0("~/analysis/WES/neoantigens/filtering_steps/tables/neoSNV_tables/", sample_C1, ".weakly_filtered_neoSNVs"), header = T, stringsAsFactors = F, sep="\t")
final_Neo_SNV <- data.frame(rbind(neo_snv_1, neo_snv_2, neo_snv_3))
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
DNA_bin_driver_annotated <- cbind(Multi_mobster_output_nontail_binary$non_tail_SNV, Multi_mobster_output_nontail_binary$driver_label_Intogen, Multi_mobster_output_nontail_binary$driver_label_COSMIC_INTOGEN, Multi_mobster_output_nontail_binary$Cluster, Multi_mobster_output_nontail_binary$Info_SNV_synonymous,  Multi_mobster_output_nontail_binary$NeoSNV_Bind_level, DNAbin)
colnames(DNA_bin_driver_annotated)[colnames(DNA_bin_driver_annotated) == "Multi_mobster_output_nontail_binary$non_tail_SNV"] <- "non_tail_SNV"
colnames(DNA_bin_driver_annotated)[colnames(DNA_bin_driver_annotated) == "Multi_mobster_output_nontail_binary$Cluster"] <- "Cluster"
colnames(DNA_bin_driver_annotated)[colnames(DNA_bin_driver_annotated) == "Multi_mobster_output_nontail_binary$driver_label_Intogen"] <- "driver_label_Intogen"
colnames(DNA_bin_driver_annotated)[colnames(DNA_bin_driver_annotated) == "Multi_mobster_output_nontail_binary$driver_label_COSMIC_INTOGEN"] <- "driver_label_COSMIC_INTOGEN"
colnames(DNA_bin_driver_annotated)[colnames(DNA_bin_driver_annotated) == "Multi_mobster_output_nontail_binary$Info_SNV_synonymous"] <- "Info_SNV_synonymous"


DNA_bin_driver_annotated$non_synon_driver_Intogen <- if_else(!(DNA_bin_driver_annotated$Info_SNV_synonymous == "synonymous SNV" | DNA_bin_driver_annotated$Info_SNV_synonymous == "unknown" | DNA_bin_driver_annotated$Info_SNV_synonymous == "."), DNA_bin_driver_annotated$driver_label_Intogen, "")
DNA_bin_driver_annotated$non_synon_driver_COSMIC_INTOGEN <- if_else(!(DNA_bin_driver_annotated$Info_SNV_synonymous == "synonymous SNV" | DNA_bin_driver_annotated$Info_SNV_synonymous == "unknown" | DNA_bin_driver_annotated$Info_SNV_synonymous == "."), DNA_bin_driver_annotated$driver_label_COSMIC_INTOGEN, "")

# List of shared and private SNVs between timepoints. Non-synonymous driver SNVs and strong or weak binding NeoSNVs are annotated in the table and manually added to the phylogenetic tree.
shared_ABC <- DNA_bin_driver_annotated[which(DNA_bin_driver_annotated[,label_A1] == 1 & DNA_bin_driver_annotated[,label_B1] == 1 & DNA_bin_driver_annotated[,label_C1] == 1), ]

shared_AB <- DNA_bin_driver_annotated[which(DNA_bin_driver_annotated[,label_A1] == 1 & DNA_bin_driver_annotated[,label_B1] == 1 & DNA_bin_driver_annotated[,label_C1] == 0), ]

shared_AC <- DNA_bin_driver_annotated[which(DNA_bin_driver_annotated[,label_A1] == 1 & DNA_bin_driver_annotated[,label_B1] == 0 & DNA_bin_driver_annotated[,label_C1] == 1), ]

shared_BC <- DNA_bin_driver_annotated[which(DNA_bin_driver_annotated[,label_A1] == 0 & DNA_bin_driver_annotated[,label_B1] == 1 & DNA_bin_driver_annotated[,label_C1] == 1), ]

private_A <- DNA_bin_driver_annotated[which(DNA_bin_driver_annotated[,label_A1] == 1 & DNA_bin_driver_annotated[,label_B1] == 0 & DNA_bin_driver_annotated[,label_C1] == 0), ]

private_B <- DNA_bin_driver_annotated[which(DNA_bin_driver_annotated[,label_A1] == 0 & DNA_bin_driver_annotated[,label_B1] == 1 & DNA_bin_driver_annotated[,label_C1] == 0), ]

private_C <- DNA_bin_driver_annotated[which(DNA_bin_driver_annotated[,label_A1] == 0 & DNA_bin_driver_annotated[,label_B1] == 0 & DNA_bin_driver_annotated[,label_C1] == 1), ]

View(shared_ABC); View(shared_AB); View(shared_AC); View(shared_BC); View(private_A); View(private_B); View(private_C)

# Export SNVs from tree branches according to received therapy
pre_treatment <- rbind(shared_ABC, private_A)
Other_info_list <- str_split(pre_treatment$non_tail_SNV,"_")
extractinfo <- function(lst,n) { sapply(lst,'[', n)}
pre_treat_df <- data.frame(sample = patient_name, chr = extractinfo(Other_info_list, 1), from = extractinfo(Other_info_list, 2), ref = extractinfo(Other_info_list, 3), to = extractinfo(Other_info_list, 4))
write.table(pre_treat_df, file = paste0(data_dir,"/pre_treatment_SNVs/", patient_name, ".pre_treatment.txt"), sep = "\t", quote = F, row.names = F)

chemo_induced <- rbind(shared_BC, private_B)
Other_info_list <- str_split(chemo_induced$non_tail_SNV,"_")
extractinfo <- function(lst,n) { sapply(lst,'[', n)}
chemo_induc_df <- data.frame(sample = patient_name, chr = extractinfo(Other_info_list, 1), from = extractinfo(Other_info_list, 2), ref = extractinfo(Other_info_list, 3), to = extractinfo(Other_info_list, 4))
write.table(chemo_induc_df, file = paste0(data_dir,"/chemo_induced_SNVs/", patient_name, ".chemo_induced_SNVs.txt"), sep = "\t", quote = F, row.names = F)

radiotherapy_induced <- private_C
Other_info_list <- str_split(radiotherapy_induced$non_tail_SNV,"_")
extractinfo <- function(lst,n) { sapply(lst,'[', n)}
radio_induc_df <- data.frame(sample = patient_name, chr = extractinfo(Other_info_list, 1), from = extractinfo(Other_info_list, 2), ref = extractinfo(Other_info_list, 3), to = extractinfo(Other_info_list, 4))
write.table(radio_induc_df, file = paste0(data_dir,"/radiochemo_induced_SNVs/", patient_name, ".radio_induced_SNVs.txt"), sep = "\t", quote = F, row.names = F)


# Plot phylogeny - comes up with error saying lots of the options aren't graphical parameters (they are...)
# Non-synonymous driver SNVs and NeoSNV counts are manually added to the plot
pdf(paste0(fig_dir, "/", patient_name, ".phyl.tree.pdf"), height = 4.8, width = 6)
plotBS(DNAtreeMP,DNABStrees,type.tree='phylogram',bs.adj=c(1.3,1.3),label.offset=10,p=0,type="phylogram", edge.color= c('#568ebd', '#0A9086', "grey", '#E64B35FF', "#7E6148FF", "#7E6148FF"), edge.width=3,font=2,cex=1.5,tip.color=c('#E64B35FF', '#0A9086', '#568ebd', "#7E6148FF"))
add.scale.bar(0,01, lwd=3,font=2);mtext(paste0('HI: ',hi),side=3,font=2,cex=1.25,ask=T);show.node.label = T;
title(main = patient_label, cex.main=1.5, col.main= "black")
dev.off()

