
#########################################################################################################
###  1. Barplot for copy number categories in Responder (REP) and Non-Responder (NRP) at each timepoint ###
###  2. Calculation of percetage of altered copy number state                                         ###
#########################################################################################################

library(reshape2)
library(ggplot2)
library(data.table)


# prepare output dirs:
fig_dir = "~/analysis/WES/CNA/plots"
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)

# Select samples to include for analysis
selected_samples <- read.table("~/analysis/General_files/MEMORI_specific/WES_samples_include_for_analysis.txt", sep = "\t", header = T)
selected_samples_Anov <- selected_samples[which(selected_samples$PASS_Sequenza_best_fitting_ploidy == "y"), ]
samples <- selected_samples_Anov$SampleID


# Create a Mastertable with 100kb binned CNS
for (sam in samples) { 
  chr_loc <- read.csv(paste0("~/analysis/WES/input_files/Sequenza/best_fitting_ploidies_100kb_binned/", samples[1], ".100kb_chrom_bin.txt"), header=TRUE,stringsAsFactors = F, sep = "\t")[,c(1:3)]
  samples <- gsub("\\S+(MEMORI_\\S+D\\d)\\S+","\\1", samples)
  df    <- do.call(cbind,lapply(samples,function(sam)paste0(read.csv(paste0("/analysis/WES/input_files/Sequenza/binned_100kb_Sequenza_files/best_fitting_Sequenza_ploidies_per_sample/", sam, ".100kb_chrom_bin.txt"), header = TRUE,stringsAsFactors = F, sep = "\t")[,5], "_", read.csv(paste0("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/multi_sample_Mobster_on_CCF/binned_100kb_Sequenza_files/best_fitting_Sequenza_ploidies_per_sample/", sam, ".100kb_chrom_bin.txt"), header = TRUE,stringsAsFactors = F, sep = "\t")[,6])))
  CNA_MasterTable   <- as.data.frame(cbind(chr_loc,df),stringsAsFactors=F)
  colnames(CNA_MasterTable)  <- c("chr", "start", "stop",samples)
}

# Eliminate NA values
CNA_MasterTable <- data.frame(lapply(CNA_MasterTable, function(x) {
  gsub("NA_NA", NA, x)
}))
CNA_MasterTable <- na.omit(CNA_MasterTable)

CNA_MasterTable$start <- NULL
CNA_MasterTable$stop <- NULL
CNA_MasterTable$chr <- NULL

Sample_ID = colnames(CNA_MasterTable)
names_OACs = as.data.frame(Sample_ID)

# Definitions of copy number categories:
# Loss = 1:0
# Balanced LOH= 2:0
# Normal diploid= 1:1
# Unbalanced amplified LOH = 3:0, 4:0, 5:0 etc
# Amplification with triploid state = 2:1
# High amplification = more than 4 alleles and no LOH-state

# Coding of CNS categories
loss <- "1_0"
normal_LOH <- "2_0"
triploid_ampl <- "2_1"
normal_diploid <- "1_1"
# high_LOH <- pattern "_0", but not "1_0" and not "2_0
# high_ampl <- total length - all other categories


# Transform exact copy number states to copy number category for each sample in each 100kb bin
Pg_CNA_OACs <- names_OACs
for (i in 1:ncol(CNA_MasterTable)) {
  samp = CNA_MasterTable[,i]
  Pg_CNA_OACs$pg_loss[i] = length(samp[samp==loss])/length(samp)
  Pg_CNA_OACs$pg_high_LOH[i] = length(samp[which(grepl("_0", samp, fixed = TRUE) & samp != "2_0" & samp != "1_0")])/length(samp)
  Pg_CNA_OACs$pg_normal_LOH[i] = length(samp[samp==normal_LOH])/length(samp)
  Pg_CNA_OACs$pg_normal_diploid[i] = length(samp[samp==normal_diploid])/length(samp)
  Pg_CNA_OACs$pg_triploid_ampl[i] = length(samp[samp==triploid_ampl])/length(samp)
  Pg_CNA_OACs$pg_high_ampl[i] = (length(samp)-length(samp[which(grepl("_0", samp, fixed = TRUE))]) - length(samp[samp==normal_diploid]) - length(samp[samp==triploid_ampl])) /length(samp)
  
}

# Calculate fraction of altered copy number state
Pg_CNA_OACs_3 <- Pg_CNA_OACs
Pg_CNA_OACs_3$pga <- 1 - Pg_CNA_OACs_3$pg_normal_diploid
mean(Pg_CNA_OACs_3$pga)
max(Pg_CNA_OACs_3$pga)
min(Pg_CNA_OACs_3$pga)

############################################################
# Create Bargraph for individual samples 
############################################################

Pg_CNA_OACs_Plot <- melt(Pg_CNA_OACs)
Pg_CNA_OACs_Plot$variable_ordered <- factor(Pg_CNA_OACs_Plot$variable , levels = c("pg_high_ampl", "pg_triploid_ampl", "pg_high_LOH", "pg_normal_LOH", "pg_loss", "pg_normal_diploid"))
colours_CNA <- c("#B2182B", "#D6604D", "#2166AC", "#4393C3", "#92C5DE","grey") 
Pg_CNA_OACs_Plot$Sample_ID <- gsub('RE', 'REP', Pg_CNA_OACs_Plot$Sample_ID)
Pg_CNA_OACs_Plot$Sample_ID <- gsub('NE', 'NRP', Pg_CNA_OACs_Plot$Sample_ID)

pdf(paste0(fig_dir, "/CNA_percentage_per_sample.pdf"), width = 18, height = 10)
p <- ggplot(Pg_CNA_OACs_Plot, aes(Sample_ID, value, fill= variable_ordered)) + geom_bar(stat= "Identity") + theme(axis.text.x = element_text(angle = 90, size = 14, face = "bold", color = "black")) + 
  theme(axis.title.y = element_text(size = 20, face = "bold"))+ theme(axis.text.y = element_text(size = 20, face = "bold", color ="black")) + theme(legend.text=element_text(size= 16, face = "bold"))+ theme(legend.title = element_text(size= 18, face = "bold")) + 
  scale_fill_manual(values = colours_CNA, labels=c("Amplification (CNt>3)", "Gain (CNt=3)", "LOH (>2:0)", "cnLOH (2:0)", "Loss (1:0)", "Diploid (1:1)")) +
  ylab("Percentage of exome") + xlab("") + labs(fill = "CN state")
plot(p)
dev.off()


############################################################
# Create Bargraph for clincical subgroups (REP and NRP for each timepoint)
############################################################

Pg_CNA_OACs_Plot_2 <-  melt(Pg_CNA_OACs)
Pg_CNA_OACs_Plot_2$Responsiveness <- sapply(Pg_CNA_OACs_Plot_2$Sample_ID, function(x) substr(x,8,9))
Pg_CNA_OACs_Plot_2$Timepoint <- sapply(Pg_CNA_OACs_Plot_2$Sample_ID, function(x) substr(x,14, 14))


Pg_CNA_OACs_Plot_2b <- Pg_CNA_OACs_Plot_2 %>% 
  select(Timepoint, Responsiveness, variable_ordered, value) %>% 
  group_by(Timepoint, Responsiveness, variable_ordered) %>% 
  summarise(adjusted_value = value/dplyr::n())

Pg_CNA_OACs_Plot_2b$adjusted_value <- as.numeric(Pg_CNA_OACs_Plot_2b$adjusted_value)
Pg_CNA_OACs_Plot_2b$variable_ordered <- factor(Pg_CNA_OACs_Plot_2b$variable_ordered , levels = c("pg_high_ampl", "pg_triploid_ampl", "pg_high_LOH", "pg_normal_LOH", "pg_loss", "pg_normal_diploid"))
colours_CNA <- c("#B2182B", "#D6604D", "#2166AC", "#4393C3", "#92C5DE","grey") 

pdf(paste0(fig_dir, "/CNA_percentage_per_group.pdf"), width = 5, height = 4) 
q <- ggplot(Pg_CNA_OACs_Plot_2b, aes(Timepoint, adjusted_value, fill= variable_ordered)) + geom_bar(stat= "Identity") + theme(axis.text.x = element_text(size = 14, face = "bold", color = "black")) + 
  theme(axis.title.y = element_text(size = 14, face = "bold"))+ theme(axis.text.y = element_text(size = 14, face = "bold", color ="black")) + theme(legend.text=element_text(size= 14, face = "bold"))+ theme(legend.title = element_text(size= 15.5, face = "bold")) + 
  scale_fill_manual(values = colours_CNA, labels = c("CNS>3", "CNS=3", ">2:0", "2:0", "1:0", "1:1")) + 
  facet_grid(~Responsiveness, scales = "free_x", space = "free_x",  labeller = as_labeller(c('NR' = "NRP", 'RE' = "REP"))) +
  theme(panel.spacing = unit(0.5, "cm"),
        strip.text = element_text(size = 14, face = "bold")) +
  ylab("percentage of exome") + xlab("") + labs(fill = "") + ylim(0,1)
plot(q)
dev.off()

# T-tests between CN categories of different timepoints for NRP and REP
tests_list_NE <- lapply(colnames(test_NE)[2:8], function(x) t_test(as.formula(paste0(x, "~ Timepoint")), data = test_NE))
result_NRP <- rbindlist(tests_list_NE)

tests_list_RE <- lapply(colnames(test_RE)[2:8], function(x) t_test(as.formula(paste0(x, "~ Timepoint")), data = test_RE))
result_REP <- rbindlist(tests_list_RE)


