
########################################################
#######  Running MOBSTER for single samples       ######
########################################################
# 1. This script computes the CCF for each SNV using MOBSTER 
# 2. Based on the CCF SNVs are classified as clonal (=C1 cluster) or subclonal (= all other C-clusters or tail mutations)


library(CNAqc)
library(grid)
library("gridExtra")
CNAqc:::get_reference("GRCh38")
require(mobster)
# R Packages CNAqc and MOBSTER are available on https://github.com/caravagnalab

# prepare output dirs:
fig_dir = "~/analysis/WES/mutations/clonality_analyses/MOBSTER/plot/mobster_model"
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)
data_dir = "~/analysis/WES/mutations/clonality_analyses/MOBSTER/output_tables"
dir.create(data_dir, showWarnings=FALSE, recursive=TRUE)


# Select samples included for genetic analyses
selected_samples <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/Study_sample_lists/WES_samples_include_for_analysis.txt", sep = "\t", header = T)
selected_samples_Seq <- selected_samples[which(selected_samples$PASS_Sequenza_best_fitting_ploidy == "y"), ]
samples <- selected_samples_Seq$SampleID

# Read in SNV lists, copy number files and purity files
for (sam in samples){  
  df <- read.table(paste0("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/mutation_results/mutations_filtered/weak_filter/", sam, ".weak.filtered.SNVs.txt"), header = T, stringsAsFactors = F)
  mut_df <- df[, c("Chr", "Start", "End", "Ref", "Alt", "ExonicFunc.refGene", "totalsum_tumour", "sum_alter_tumour", "VAF_manually")]
  colnames(mut_df)[colnames(mut_df) == "Chr"] <- "chr"
  colnames(mut_df)[colnames(mut_df) == "Start"] <- "from"
  colnames(mut_df)[colnames(mut_df) == "End"] <- "to"
  colnames(mut_df)[colnames(mut_df) == "Ref"] <- "ref"
  colnames(mut_df)[colnames(mut_df) == "Alt"] <- "alt"
  colnames(mut_df)[colnames(mut_df) == "totalsum_tumour"] <- "DP"
  colnames(mut_df)[colnames(mut_df) == "sum_alter_tumour"] <- "NV"
  colnames(mut_df)[colnames(mut_df) == "VAF_manually"] <- "VAF"
  
  
  sequenza_df <- read.table(paste0("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/", sam, "_segments.txt"), header = T, stringsAsFactors = F)
  cnv_df <- sequenza_df[, c("chromosome", "start.pos", "end.pos", "A", "B")]
  colnames(cnv_df)[colnames(cnv_df) == "chromosome"] <- "chr"
  colnames(cnv_df)[colnames(cnv_df) == "start.pos"] <- "from"
  colnames(cnv_df)[colnames(cnv_df) == "end.pos"] <- "to"
  colnames(cnv_df)[colnames(cnv_df) == "A"] <- "Major"
  colnames(cnv_df)[colnames(cnv_df) == "B"] <- "minor"
  
  pur_df <- read.table(paste0("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/", sam, "_confints_CP.txt"), header = T, stringsAsFactors = F, fill = T) 
  pur_df <- pur_df[2,1]
  
  x <- CNAqc::init(mutations = mut_df, cna = cnv_df, purity= pur_df, ref = "GRCh38")


# Run MOBSTER and calculate CCF 
x = compute_CCF(x, method = 'ROUGH') 
y = CCF(x)
y = y[!is.na(y$CCF), , drop = FALSE]
y$original_VAF = y$VAF

# Calculate adjusted VAF from CCF  to assign SNVs to clonal and subclonal clusters
# As OAC shows mainly non-diploid CNS, the adjusted VAF needs to be calculated as CCF/2
y$VAF = y$CCF
y <- y %>%
  mutate(VAF = CCF / 2) %>%
  filter(VAF < 1.0, VAF > (0.025/pur_df))

snv_cna <- compute_CCF(x) 


CCF_mobster_fit = mobster_fit(y, parallel = FALSE, auto_setup = "FAST")

# Calculate and plot Mobster fit using the CCF
g <- plot(CCF_mobster_fit$best)
p <- plot_data_histogram(snv_cna, which = "CCF")

pdf(paste0(fig_dir, "/", sam, ".mobster.on.CCF.pdf"))
pushViewport(viewport(layout = grid.layout(2, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(g, vp = vplayout(1, 1:2))  # key is to define vplayout
print(p, vp = vplayout(2, 1:2))
dev.off()
 

# Export Mobster table with adjusted VAF calculated from CCF
colnames(CCF_mobster_fit$best$data)[colnames(CCF_mobster_fit$best$data) == "VAF"] <- "Mobster_adj_CCF_VAF"
write.table(CCF_mobster_fit$best$data, paste0(data_dir, "/", sam, ".Mobster.CCF.txt"), sep = "\t", row.names = F, quote = F)

}







