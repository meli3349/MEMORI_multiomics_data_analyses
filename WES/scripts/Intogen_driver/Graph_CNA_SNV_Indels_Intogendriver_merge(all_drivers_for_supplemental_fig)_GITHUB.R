
##########################################################################
## Plot CNA and mutation/indel status in Intogen drivers  ################
##########################################################################

library(data.table)
library(dplyr)
library(rlist)
library("RColorBrewer")

# prepare output dirs:
fig_dir = "~/analysis/WES/Intogen_driver/plots"
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)

# Read in Intogen list with driver genes form Esophageal Cancer and gastric cancer 
# https://www.intogen.org/ 
Intogen_EC_GC_list  <- read.table("~/analysis/General_files/external_files/Intogen_drivers.txt", stringsAsFactors = F, header = T)

# ENSEMBLE/ HUGO Gene ID and gene location list
ENSEMBLE_list<- read.table("~/analysis/General_files/external_files/compiledGeneInfo.txt", stringsAsFactors = F, header = T)
chr_filt <- paste0("chr", c(1:23, "X", "Y"))
ENSEMBLE_list <- ENSEMBLE_list[which(ENSEMBLE_list$Chr %in% chr_filt),]

# Get ENSEMBL Gene Names and gene location for Oncogenes/TSG from Intogen list
Intogen_loc <- merge(Intogen_EC_GC_list, ENSEMBLE_list, by.x = "Symbol", by.y = "Name")
Intogen_loc  <- Intogen_loc[ ,c("Symbol", "GeneID", "Chr", "Strand", "Start", "End")]

# Select samples included for genetic analyses
selected_samples <- read.table("~/analysis/General_files/MEMORI_specific/WES_samples_include_for_analysis.txt", sep = "\t", header = T)
selected_samples_Seq <- selected_samples[which(selected_samples$PASS_Sequenza_best_fitting_ploidy == "y"), ]
samples <- selected_samples_Seq$SampleID

# Table per sample with Copy number state for each driver 
List_driver_CNS <- list()
for (sam in samples) { 
  df <- read.table(paste0("~/analysis/WES/input_files/Sequenza/best_fitting_ploidies_per_sample/", sam, "_segments.txt"), header = T, stringsAsFactors = F)
  colnames(df)[colnames(df) == "chromosome"] <- "Chr"
  colnames(df)[colnames(df) == "start.pos"] <- "Start"
  colnames(df)[colnames(df) == "end.pos"] <- "End"
  df$copy_state <- paste0(df$A, "_", df$B)
  df <- df[,c("Chr", "Start", "End", "copy_state", "CNt")]
  
  Intogen_loc  <- as.data.table(Intogen_loc)
  df <- as.data.table(df)
  
  setkey(df, Chr, Start, End)
  merged_df <- foverlaps(Intogen_loc, df, type="within")
  
  export_df_1 <- merged_df[,c("copy_state")]
  List_driver_CNS[[sam]] <- export_df_1 
  genes <- merged_df[,c("Symbol")]
}

# Combine Copy number state for driver genes from each sample to a mastertable
df <- list.cbind(List_driver_CNS)
MasterTable_detailed <- as.data.frame(cbind(genes,df),stringsAsFactors=F)
colnames(MasterTable_detailed)  <- c("GeneID", samples)

# Add clincal information (patientID, timepoint, cellularity, Becker remission score)
patient <-c("Patient", sapply(colnames(MasterTable_detailed[2:ncol(MasterTable_detailed)]), function(x) substr(x,8, 13)))
Timepoint <-c("Timepoint", sapply(colnames(MasterTable_detailed[2:ncol(MasterTable_detailed)]), function(x) substr(x,14, 14)))

  for (sam in samples) { 
    cellularity_df    <- do.call(cbind,lapply(samples,function(sam)read.table(paste0("~/analysis/WES/input_files/Sequenza/best_fitting_ploidies_per_sample/confints_files/", sam, "_confints_CP.txt"), header = T, stringsAsFactors = F, fill = T)[2,1]))
    colnames(cellularity_df) <- samples
    for(sam in samples) { cellularity_df[,sam] <- as.numeric(cellularity_df[,sam])}}
x <- as.numeric(cellularity_df[1,])
cell_lev = cut(x, c(0.17,0.4,0.6,1))
levels(cell_lev) = c("lowCell","mediumCell","highCell")
cellularity <- c("Cellularity", as.character(cell_lev))

remission_df<- read.table("~/analysis/General_files/MEMORI_specific/WES_samples_include_for_analysis.txt", sep = "\t", header = T)
remission_df <- remission_df[,c("SampleID", "Becker_remission_grade")]
remission_df <- remission_df[which(remission_df$SampleID %in% samples), ]
y <- as.character(remission_df$Becker_remission_grade)
remission_score <- c("Remission",as.character(y))
  
# Add clincial information to Mastertable
MasterTable_detailed_clinical <- rbind(patient, Timepoint, cellularity, remission_score, MasterTable_detailed)

# For plotting REP and NRP tables 
melt_MasterTable_detailed_clinical <- melt(MasterTable_detailed_clinical, id="GeneID")
melt_MasterTable_detailed_clinical$colID <- melt_MasterTable_detailed_clinical$value

# Combining copy number states into copy number categories
loss <- "1_0" 
high_LOH <- c("20_0", "19_0", "18_0", "17_0", "16_0", "15_0", "14_0", "13_0", "12_0", "11_0", "10_0", "9_0", "8_0", "7_0", "6_0", "5_0", "4_0", "3_0")
normal_LOH <- "2_0"
triploid_ampl <- "2_1"
normal_diploid <- "1_1"
# all copy number states with x_1 state and x>2 are categorized as Amplification > 3

melt_MasterTable_detailed_clinical$colID[melt_MasterTable_detailed_clinical$colID %in% loss] <- "Loss 1:0"
melt_MasterTable_detailed_clinical$colID[melt_MasterTable_detailed_clinical$colID %in% normal_LOH] <- "LOH 2:0"
melt_MasterTable_detailed_clinical$colID[melt_MasterTable_detailed_clinical$colID %in% high_LOH] <- "LOH >2:0"
melt_MasterTable_detailed_clinical$colID[melt_MasterTable_detailed_clinical$colID %in% triploid_ampl] <- "Amplification = 3"
melt_MasterTable_detailed_clinical$colID[melt_MasterTable_detailed_clinical$colID %in% normal_diploid] <- "diploid"
melt_MasterTable_detailed_clinical$colID[which(grepl("_", melt_MasterTable_detailed_clinical$colID, fixed = TRUE))] <- "Amplification > 3"
melt_MasterTable_detailed_clinical$colID[is.na(melt_MasterTable_detailed_clinical$colID)] <- "unknown"

all_genes <- MasterTable_detailed$GeneID
melt_MasterTable_detailed_clinical$rows_ordered <- factor(melt_MasterTable_detailed_clinical$GeneID, levels = rev(c("Patient", "Timepoint", "Cellularity", "Remission", all_genes)))


## Add indels in driver genes to Mastertable 
# MasterTable_indels shows cancer cell fraction (CCF) harbouring indels in driver genes. NA indicates no indel in cancer driver.
# Only indels with CCF>0.1 and <1.4 are included, as these are defined as biologically plausible and relevant
MasterTable_indels <- read.table("~/analysis/WES/input_files/mutation_results/Intogen_driver/Mastertable_CCF_indels_Intogen_drivers.txt", stringsAsFactors = F, header = T)
melt_MasterTable_indels <- melt(MasterTable_indels, id="GeneID")
melt_MasterTable_indels$colID_indels <- ifelse(melt_MasterTable_indels$value> 0.1 & melt_MasterTable_indels$value <1.4, "Indel", NA)
melt_MasterTable_detailed_clinical <- merge(melt_MasterTable_detailed_clinical, melt_MasterTable_indels, by = c("GeneID", "variable"), all.x = T)


## Add SNVs in driver genes to Mastertable 
# MasterTable_SNVs shows cancer cell fraction (CCF) harbouring SNVs in driver genes. NA indicates no indel in cancer driver.
# Only SNVs with CCF>0.1 and <1.4 are included, as these are defined as biologically plausible and relevant
MasterTable_SNVs <- read.table("~/analysis/WES/input_files/mutation_results/Intogen_driver/Mastertable_CCF_SNVs_Intogen_drivers.txt", stringsAsFactors = F, header = T)
melt_MasterTable_SNVs <- melt(MasterTable_SNVs, id="GeneID")
melt_MasterTable_SNVs$colID_SNV <- ifelse((melt_MasterTable_SNVs$value> 0.1 & melt_MasterTable_SNVs$value <1.4), "SNV", NA)
melt_MasterTable_detailed_clinical <- merge(melt_MasterTable_detailed_clinical, melt_MasterTable_SNVs, by = c("GeneID", "variable"), all.x = T)


### Plot data 
# for correct colours in the plot, information of CNA in drivers with no SNV and no Indel needs to be transferred to columns colID_SNV and colID_indels
melt_MasterTable_detailed_clinical$colID_SNV <- ifelse(is.na(melt_MasterTable_detailed_clinical$colID_SNV), melt_MasterTable_detailed_clinical$colID, melt_MasterTable_detailed_clinical$colID_SNV)
melt_MasterTable_detailed_clinical$colID_indels <- ifelse(is.na(melt_MasterTable_detailed_clinical$colID_indels), melt_MasterTable_detailed_clinical$colID, melt_MasterTable_detailed_clinical$colID_indels)

# Plot NRP and REP separately
melt_MasterTable_detailed_clinical_NE <- melt_MasterTable_detailed_clinical[which(grepl("MEMORI_NE", melt_MasterTable_detailed_clinical$variable, fixed = TRUE)), ]
melt_MasterTable_detailed_clinical_RE <- melt_MasterTable_detailed_clinical[which(grepl("MEMORI_RE", melt_MasterTable_detailed_clinical$variable, fixed = TRUE)), ]

#Adjust colours
melt_MasterTable_detailed_clinical_NE$colours_ordered <- factor(melt_MasterTable_detailed_clinical_NE$colID, levels = c("A", "Amplification = 3", "Amplification > 3", "B","C", "diploid", "LOH >2:0", "LOH 2:0", "Loss 1:0", "1", "1a", "1b", "2", "3", "lowCell", "mediumCell", "highCell", "NE0005", "NE0006", "NE0007", "NE0008", "NE0009", "NE0010", "NE0012", "NE0020", "NE0021", "NE0022", "unknown"))
colours_NE <- c("gray68", "gray50", "gray32", '#e63e00', '#D6604D', '#B2182B', 'darkseagreen3', 'steelblue', "grey", "#948ab1", "yellow", "#2166AC", "#4393C3", '#92C5DE', "#cec5e7", "#b9adde", "lightblue", "orange", "khaki2", "forestgreen", "darkmagenta", "firebrick3", "chocolate4","coral2", "plum2", "slateblue4", "yellow", "lightgrey", "black", "black", "green")

melt_MasterTable_detailed_clinical_RE$colours_ordered <- factor(melt_MasterTable_detailed_clinical_RE$colID, levels = c("A", "Amplification = 3", "Amplification > 3", "B","C", "diploid", "LOH >2:0", "LOH 2:0", "Loss 1:0", "unknown", "1", "1b", "2", "3", "lowCell", "mediumCell", "highCell", "RE0002", "RE0003", "RE0004", "RE0005", "RE0006", "RE0007", "RE0008", "RE0009", "RE0010", "RE0011", "RE0012", "RE0014", "RE0021", "RE0022", "RE0023", "RE0024", "RE0025"))
colours_RE <- c("gray68", "gray50", "gray32",  '#e63e00', '#D6604D', '#B2182B',  'darkseagreen3', 'steelblue', "grey", "#948ab1", "yellow", "#2166AC", "#4393C3", '#92C5DE', "#cec5e7", "#b9adde", "lightblue", "orange", "khaki2", "forestgreen", "darkmagenta", "firebrick3", "chocolate4","coral2", "plum2", "slateblue4", "aquamarine3", "goldenrod3", "chocolate","chocolate4", "royalblue", "violetred3", "palegreen2", "yellow", "lightgrey")


pdf(paste0(fig_dir, "/NE_Intogendriver_CNA_all_drivers.pdf"), height = 25, width = 14)
p <- ggplot(melt_MasterTable_detailed_clinical_NE, aes(x = variable, y = rows_ordered)) +
  geom_tile(aes(width = 0.85,height = 0.85 , fill = colours_ordered)) +
  geom_tile(aes(width = 0.85,height = 0.2 , fill = colID_indels)) +
  geom_tile(aes(width = 0.2,height = 0.85 , fill = colID_SNV)) +
  coord_fixed() + 
  scale_fill_manual(values = colours_NE) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(color="black", size=14, face="bold"),
    legend.title = element_blank(),
    legend.text = element_text(color = "black", size = 14, face="bold"),
    legend.position = "bottom")
plot(p)
dev.off()


pdf(paste0(fig_dir, "/RE_Intogendriver_CNA_all_drivers.pdf"), height = 25, width = 14)
p <- ggplot(melt_MasterTable_detailed_clinical_RE, aes(x = variable, y = rows_ordered)) +
  geom_tile(aes(width = 0.85,height = 0.85 , fill = colours_ordered)) +
  geom_tile(aes(width = 0.85,height = 0.2 , fill = colID_indels)) +
  geom_tile(aes(width = 0.2,height = 0.85 , fill = colID_SNV)) +
  coord_fixed() + 
  scale_fill_manual(values = colours_RE) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(color="black", size=14, face="bold"),
    legend.title = element_blank(),
    legend.text = element_text(color = "black", size = 14, face="bold"),
    legend.position = "bottom")
plot(p)
dev.off()

