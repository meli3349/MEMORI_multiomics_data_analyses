# dN/dS analysis 

library("dndscv")
library("stringr")

# Import SNV file with column "sampleID", "chr", "pos", "ref", "mut"
patient_name <- "MEMORI_RE0005"

Multi_mobster_output_nontail <- read.table(paste0("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/multi_sample_Mobster_on_CCF/Mobster_output_tables/", patient_name, "/", "bam_read_VAF/weakly_filtered/", patient_name, ".driver_cluster_annotated_nontail_Mobster_output.txt"), header = T, sep = "\t", stringsAsFactors = F)
trunc_df <- Multi_mobster_output_nontail[which(Multi_mobster_output_nontail$MEMORI_RE0005A_T1_C1_D1 != 1 &  Multi_mobster_output_nontail$MEMORI_RE0005B_T1_B1_D1 != 0 & Multi_mobster_output_nontail$MEMORI_RE0005C_T1_C1_D1 != 0), ]
Multi_mobster_output_nontail <- trunc_df

# Seperate SNVs into chr, pos, ref, mut
Other_info_list <- str_split(Multi_mobster_output_nontail$non_tail_SNV,"_")
extractinfo <- function(lst,n) { sapply(lst,'[', n)}

SNV_hg38_df <- data.frame("chr" = extractinfo(Other_info_list ,1), "pos" = extractinfo(Other_info_list ,2), 
                          "ref" = extractinfo(Other_info_list ,3), "mut" = extractinfo(Other_info_list ,4))

SNV_hg38_df$sampleID <- patient_name
#SNV_hg38_df$chr<- substr(SNV_hg38_df$chr,4,nchar(SNV_hg38_df$chr))


mutations <- SNV_hg38_df

mutations <- mutations[,c(5,1:4)]

#obj <- readRDS("/Users/melissaschmidt/Downloads/RefCDS_human_GRCh38.p12_dNdScv.0.1.0.rda")
#obj <- SCTransform(obj, vst.flavor="v2")

load("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/dNdS_analysis/covariates_hg19_hg38_epigenome_pcawg.rda")
#load("/Users/melissaschmidt/Downloads/RefCDS_human_hg19_GencodeV18_newcovariates.rda")
#covs <- as.data.frame(covs)

dndsout = dndscv(mutations, refdb = "/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/dNdS_analysis/RefCDS_human_GRCh38.p12_dNdScv.0.1.0.rda", cv=NULL)
sel_cv = dndsout$sel_cv
print(head(sel_cv), digits = 3)

# find 
signif_genes = sel_cv[sel_cv$qtrunc_cv<0.1, c("gene_name","qtrunc_cv")]
rownames(signif_genes) = NULL
print(signif_genes)

# Pirnt globale dNdS values
print(dndsout$globaldnds)


#############################################################################################
#### Run dNdS on treatment-naive SNVs, CTx induced SNVs, RCTX induced SNVs


input_per_patient_treatment_naive_RE <- list.files("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/mutation_results/phylogenetic_trees/weakly_filtered_from_bamreadcount/phyl_tree_branch_tables/pre_treatment_SNVs/",  pattern = "MEMORI_RE*", full.names = T)
input_per_patient_treatment_naive_NE <- list.files("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/mutation_results/phylogenetic_trees/weakly_filtered_from_bamreadcount/phyl_tree_branch_tables/pre_treatment_SNVs/",  pattern = "MEMORI_NE*", full.names = T)
input_per_patient_chemo_induced <- list.files("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/mutation_results/phylogenetic_trees/weakly_filtered_from_bamreadcount/phyl_tree_branch_tables/chemo_induced_SNVs/",  pattern = "MEMORI_*", full.names = T)
input_per_patient_radio_induced <- list.files("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/mutation_results/phylogenetic_trees/weakly_filtered_from_bamreadcount/phyl_tree_branch_tables/radiochemo_induced_SNVs/",  pattern = "*.radio_induced_SNVs.txt", full.names = T)

for (x in input_per_patient_treatment_naive_RE){
  input_deconsigs_treat_naive_RE <- do.call(rbind,lapply(input_per_patient_treatment_naive_RE, function(x)read.table(x, header = T ,stringsAsFactors = F)))}
colnames(input_deconsigs_treat_naive_RE) <- c("sampleID", "chr", "pos", "ref", "to")

for (x in input_per_patient_treatment_naive_NE){
  input_deconsigs_treat_naive_NE <- do.call(rbind,lapply(input_per_patient_treatment_naive_NE, function(x)read.table(x, header = T ,stringsAsFactors = F)))}
colnames(input_deconsigs_treat_naive_NE) <- c("sampleID", "chr", "pos", "ref", "to")

for (x in input_per_patient_chemo_induced){
  input_deconsigs_chemo_ind <- do.call(rbind,lapply(input_per_patient_chemo_induced, function(x)read.table(x, header = T ,stringsAsFactors = F)))}
colnames(input_deconsigs_chemo_ind) <- c("sampleID", "chr", "pos", "ref", "to")

for (x in input_per_patient_radio_induced){
  input_deconsigs_radio_ind <- do.call(rbind,lapply(input_per_patient_radio_induced, function(x)read.table(x, header = T ,stringsAsFactors = F)))}
colnames(input_deconsigs_radio_ind) <- c("sampleID", "chr", "pos", "ref", "to")

# select which df

df_to_select <- c("input_deconsigs_treat_naive_RE", "input_deconsigs_treat_naive_NE", "input_deconsigs_chemo_ind", "input_deconsigs_radio_ind")

for (i in df_to_select){
mutations <- get(i)

load("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/dNdS_analysis/covariates_hg19_hg38_epigenome_pcawg.rda")

dndsout = dndscv(mutations, refdb = "/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/dNdS_analysis/RefCDS_human_GRCh38.p12_dNdScv.0.1.0.rda", cv=NULL)
sel_cv = dndsout$sel_cv
print(head(sel_cv), digits = 3)

# find 
signif_genes = sel_cv[sel_cv$qtrunc_cv<0.1, c("gene_name","qtrunc_cv")]
rownames(signif_genes) = NULL
print(signif_genes)

# Pirnt globale dNdS values
print(dndsout$globaldnds)
write.table(dndsout$globaldnds, file = paste0("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/dNdS_analysis/output/", i, "globaldNdS.txt"), quote = F, row.names = F, sep = "\t")
}
