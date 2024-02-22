
################################################################################################
# 1. Barplot of imumune escaped samples in REP and NRP
#### a. Overview graph
#### b. Barplot with immune escape mechanisms
#### c. Barplot with HLA LOH

# 2. Barplot of imumune escaped samples per Timepoint
#### a. Overview graph (REP and NRP)
#### a. for REP 
#### a. for NRP

# 3. Barplot of imumune escaped samples in patients with pathological remission grade 1-3 
#### a. Barplot with samples classified as escaped or non-escaped

# 4. Barplot with occurence timepoint of genetic and transcriptomic immune escape 


library(reshape2); 
library(ggplot2);
library(dplyr)


# prepare output dirs:
fig_dir = "~/analysis/multi_omic/immune_escape/plots/Bargraphs"
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)

# my themes
theme_mypub_grid <- function(base_size = 14,
                             base_family = ""){
  theme_bw(base_size = base_size, 
           base_family = base_family) %+replace%
    theme(
      panel.grid.major = element_blank(),   
      panel.grid.minor = element_blank(), 
      panel.background = element_blank(),
      complete = TRUE
    )
}

my_theme <- theme(
  plot.title = element_text(color="black", size=16, face="bold"),
  axis.title.x = element_text(color="black", size=16, face="bold"),
  axis.text.x = element_text(color="black", size=14, face="bold"),
  axis.title.y = element_text(color="black", size=16, face="bold"),
  axis.text.y = element_text(color="black", size=14, face="bold"),
  legend.title = element_text(color = "black", size = 14, face="bold"),
  legend.text = element_text(color = "black", size = 14, face="bold"),
  
)

# Read in Mastertable with immune scape mechanisms per sample
escape_df <- read.table("~/analysis/multi_omic/immune_escape/tables/MEMORI_escape_master_file_incl_PDL1.txt", stringsAsFactors = F, header = T, sep = "\t")

# Add clinical data 
escape_df$Responsiveness <- sapply(escape_df$Sample, function(x) substr(x,8,9)) 
escape_df$Timepoint <-  sapply(escape_df$Sample, function(x) substr(x,14,14)) 

# Analysis of immune escape mechanism in samples with available WES and RNA data
escape_df_complete <- escape_df[complete.cases(escape_df[ ,c("RNA_file","WES_file")]),]
escape_df_complete$final_escape_mechanism <- ifelse((escape_df_complete$PDL1_immune_escape == "PDL1_high" & escape_df_complete$HLAmut ==  "FALSE" & escape_df_complete$final_HLA_LOH  == "FALSE" & escape_df_complete$B2Mmut == "FALSE"), "PDL1_high", 
                                           ifelse((escape_df_complete$PDL1_immune_escape == "FALSE" & escape_df_complete$HLAmut ==  "mutated" & escape_df_complete$final_HLA_LOH  == "FALSE" & escape_df_complete$B2Mmut == "FALSE"), "HLAmut",
                                                  ifelse((escape_df_complete$PDL1_immune_escape == "FALSE" & escape_df_complete$HLAmut ==  "FALSE" & escape_df_complete$final_HLA_LOH  == "LOH" & escape_df_complete$B2Mmut == "FALSE"), "HLA_LOH",
                                                         ifelse((escape_df_complete$PDL1_immune_escape == "FALSE" & escape_df_complete$HLAmut ==  "FALSE" & escape_df_complete$final_HLA_LOH  == "FALSE" & escape_df_complete$B2Mmut == "mutated"), "B2Mmut", 
                                                                ifelse((escape_df_complete$PDL1_immune_escape == "FALSE" & escape_df_complete$HLAmut ==  "FALSE" & escape_df_complete$final_HLA_LOH  == "FALSE" & escape_df_complete$B2Mmut == "FALSE"), "no escape", "mixed")))))




####################################
### 1. Immune escape per responsiveness per sample

# 1a.) Overview
escape_df_responsiveness <- escape_df_complete[,c("final_escape_result", "Responsiveness")]
Sum_df <- table(escape_df_responsiveness)
Sum_df_perc <- scale(Sum_df, FALSE, colSums(Sum_df)) * 100
melt_Sum_df <- melt(Sum_df_perc)

pdf(paste0(fig_dir, "/Barplot_escape_RE_NR_per_sample.pdf"), width = 6, height = 5)                                       
q <- ggplot(melt_Sum_df, aes(Responsiveness, value, fill= final_escape_result)) + geom_bar(stat= "Identity") + theme(axis.text.x = element_text(size = 14, face = "bold", color = "black"))  +
  theme(axis.title.y = element_text(size = 14, face = "bold"))+ theme(axis.text.y = element_text(size = 14, face = "bold", color ="black")) + theme(legend.text=element_text(size= 14, face = "bold"))+ theme(legend.title = element_text(size= 15.5, face = "bold")) + 
  scale_fill_manual(values = c("gray67", "khaki2"), labels = c("no", "yes")) + 
  ylab("percentage of samples") + xlab("") + labs(fill = "Immune escaped") + ylim(0,105) + scale_x_discrete(labels=c("NE" = "NRP", "RE" = "REP")) + theme_mypub_grid() + my_theme
plot(q) 
dev.off()

# stats
chisq_NE_RE <- chisq.test(Sum_df)
dfmean_RE_NR <- escape_df_responsiveness %>% 
  group_by(final_escape_result, Responsiveness) %>% 
  dplyr::summarise(count = dplyr::n())
names(dfmean_RE_NR)[1] <- "Category_sample"
dfmean_RE_NR$Test <- "chisq"
dfmean_RE_NR$pVal <- chisq_NE_RE$p.value
write.table(dfmean_RE_NR, file = paste0(fig_dir, "/chisq_Barplot_escape_RE_NR_samples.txt"), quote = F, row.names = F)                                       

# 1b.) Detailed immune escape mechanism

escape_df_mechanism_responsiveness <- escape_df_complete[,c("final_escape_mechanism", "Responsiveness")]
Sum_df <- table(escape_df_mechanism_responsiveness)
Sum_df_perc <- scale(Sum_df, FALSE, colSums(Sum_df)) * 100
melt_Sum_df <- melt(Sum_df_perc)
melt_Sum_df$final_escape_mechanism <- factor(melt_Sum_df$final_escape_mechanism, levels = c("no escape", "mixed", "PDL1_high", "HLA_LOH"))


pdf(paste0(fig_dir, "/Barplot_escape_RE_NR_per_sample_including_escape_mechanism.pdf"), width = 7, height = 4.5)                                       
q <- ggplot(melt_Sum_df, aes(Responsiveness, value, fill= final_escape_mechanism)) + geom_bar(stat= "Identity") + theme(axis.text.x = element_text(size = 14, face = "bold", color = "black")) +
  theme(axis.title.y = element_text(size = 14, face = "bold"))+ theme(axis.text.y = element_text(size = 14, face = "bold", color ="black")) + theme(legend.text=element_text(size= 14, face = "bold"))+ theme(legend.title = element_text(size= 15.5, face = "bold")) + 
  scale_fill_manual("Immune escape mechanisms", values = c("grey80", 'darkseagreen3', "#B2182B",'steelblue'), labels = c("no escape", "mixed", "PDL1 overexpression", "HLA LOH"))  +
  ylab("percentage of samples") + xlab("")  + ylim(0,105) + scale_x_discrete(labels=c("NE" = "NRP", "RE" = "REP"))  + theme_mypub_grid() + my_theme
plot(q)
dev.off()

# stats
chisq_mechanisms<- chisq.test(Sum_df)
dfmean_mechanism <- escape_df_mechanism_responsiveness %>% 
  group_by(final_escape_mechanism, Responsiveness) %>% 
  dplyr::summarise(count = dplyr::n())
names(dfmean_mechanism)[1] <- "Category_sample"
dfmean_mechanism$Test <- "chisq"
dfmean_mechanism$pVal <- chisq_mechanisms$p.value
write.table(dfmean_mechanism, file = paste0(fig_dir, "/chisq_Barplot_escape_mechanism_RE_NR_samples.txt"), quote = F, row.names = F)                                       


# 1c.) HLA LOH in REP and NRP

escape_df_HLALOH <- escape_df[complete.cases(escape_df[c("WES_file")]), c("final_HLA_LOH", "Responsiveness") ]
Sum_df <- table(escape_df_HLALOH)
Sum_df_perc <- scale(Sum_df, FALSE, colSums(Sum_df)) * 100
melt_Sum_df <- melt(Sum_df_perc)


pdf(paste0(fig_dir, "/Barplot_escape_RE_NR_per_sample_HLA_LOH.pdf"), width = 6, height = 5)                                       
q <- ggplot(melt_Sum_df, aes(Responsiveness, value, fill= final_HLA_LOH)) + geom_bar(stat= "Identity") + theme(axis.text.x = element_text(size = 14, face = "bold", color = "black"))  +
  theme(axis.title.y = element_text(size = 14, face = "bold"))+ theme(axis.text.y = element_text(size = 14, face = "bold", color ="black")) + theme(legend.text=element_text(size= 14, face = "bold"))+ theme(legend.title = element_text(size= 15.5, face = "bold")) + 
  scale_fill_manual("HLA LOH", values = c("gray67", "khaki2"), labels = c("no", "yes")) + 
  ylab("percentage of samples") + xlab("") + labs(fill = "final_HLA_LOH") + ylim(0,105) + scale_x_discrete(labels=c("NE" = "NRP", "RE" = "REP")) + theme_mypub_grid() + my_theme
plot(q) 
dev.off()

# stats
chisq_HLALOH <- chisq.test(Sum_df)
dfmean_HLALOH <- escape_df_HLALOH %>% 
  group_by(final_HLA_LOH, Responsiveness) %>% 
  dplyr::summarise(count = dplyr::n())
names(dfmean_HLALOH)[1] <- "Category_sample"
dfmean_HLALOH$Test <- "chisq"
dfmean_HLALOH$pVal <- chisq_HLALOH$p.value
write.table(dfmean_HLALOH, file = paste0(fig_dir, "/chisq_Barplot_HLA_LOH_RE_NR_samples.txt"), quote = F, row.names = F)                                       



####################################
### 2. Immune escape per timepoint
escape_df_timepoint <- escape_df_complete[,c("final_escape_result", "Timepoint")]
Sum_df <- table(escape_df_timepoint)
Sum_df_perc <- scale(Sum_df, FALSE, colSums(Sum_df)) * 100
melt_Sum_df <- melt(Sum_df_perc)

# 2a. Overview
pdf(paste0(fig_dir, "/Barplot_escape_per_Timepoint.pdf"), width = 6, height = 5)                                       
q <- ggplot(melt_Sum_df, aes(Timepoint, value, fill= final_escape_result)) + geom_bar(stat= "Identity") + theme(axis.text.x = element_text(size = 14, face = "bold", color = "black")) + 
  theme(axis.title.y = element_text(size = 14, face = "bold"))+ theme(axis.text.y = element_text(size = 14, face = "bold", color ="black")) + theme(legend.text=element_text(size= 14, face = "bold"))+ theme(legend.title = element_text(size= 15.5, face = "bold")) + 
  scale_fill_manual(values = c("gray67", "khaki2"), labels = c("no", "yes")) + #, "lightgreen", "khaki2", "forestgreen", "plum2", "firebrick3", "firebrick1", "orange", "floralwhite", "pink", "coral2", "darkmagenta", "slateblue4", "goldenrod3", "chocolate","chocolate4", "black", "grey33", "grey")) 
  ylab("percentage of samples") + xlab("") + labs(fill = "Immune escaped") + ylim(0,105)  + theme_mypub_grid() + my_theme
plot(q)
dev.off()

# stats
chisq_TP<- chisq.test(Sum_df)
dfmean_TP <- escape_df_timepoint %>% 
  group_by(interaction(final_escape_result, Timepoint)) %>% 
  dplyr::summarise(count = dplyr::n())
names(dfmean_TP)[1] <- "Category_sample"
dfmean_TP$Test <- "chisq"
dfmean_TP$pVal <- chisq_TP$p.value
write.table(dfmean_TP, file = paste0(fig_dir, "/chisq_Barplot_escape_per_Timepoint.txt"), quote = F, row.names = F)                                       


# 2b. Immune escape per Timemoint for REP
escape_df_RE <- escape_df_complete[which(escape_df$Responsiveness == "RE"),c("final_escape_result", "Timepoint")]
Sum_df <- table(escape_df_RE)
Sum_df_perc <- scale(Sum_df, FALSE, colSums(Sum_df)) * 100
melt_Sum_df <- melt(Sum_df_perc)

pdf(paste0(fig_dir, "/Barplot_escape_RE_per_Timepoint.pdf"), width =6, height = 5)                                       
q <- ggplot(melt_Sum_df, aes(Timepoint, value, fill= final_escape_result)) + geom_bar(stat= "Identity") + theme(axis.text.x = element_text(size = 14, face = "bold", color = "black")) + 
  theme(axis.title.y = element_text(size = 14, face = "bold"))+ theme(axis.text.y = element_text(size = 14, face = "bold", color ="black")) + theme(legend.text=element_text(size= 14, face = "bold"))+ theme(legend.title = element_text(size= 16, face = "bold")) +  
  scale_fill_manual(values = c("gray67", "khaki2"), labels = c("no", "yes")) + 
  ylab("percentage of samples") + xlab("") + labs(fill = "Immune escaped") + ylim(0,105)  + theme_mypub_grid() + my_theme
plot(q)
dev.off()

# stats
chisq_RE_TP <- chisq.test(Sum_df)
dfmean_chisq_RE_TP <- escape_df_RE %>% 
  group_by(interaction(final_escape_result, Timepoint)) %>% 
  dplyr::summarise(count = dplyr::n())
names(dfmean_chisq_RE_TP)[1] <- "Category_sample"
dfmean_chisq_RE_TP$Test <- "chisq"
dfmean_chisq_RE_TP$pVal <- chisq_RE_TP$p.value
write.table(dfmean_chisq_RE_TP, file = paste0(fig_dir, "/chisq_Barplot_escape_RE_per_Timepoint.txt"), quote = F, row.names = F)                                       


# 2b. Immune escape per Timemoint for NRP
escape_df_NR <- escape_df_complete[which(escape_df$Responsiveness == "NE"),c("final_escape_result", "Timepoint")]
Sum_df <- table(escape_df_NR)
Sum_df_perc <- scale(Sum_df, FALSE, colSums(Sum_df)) * 100
melt_Sum_df <- melt(Sum_df_perc)

pdf(paste0(fig_dir, "/Barplot_escape_NR_per_Timepoint.pdf"), width = 6, height = 5)                                       
q <- ggplot(melt_Sum_df, aes(Timepoint, value, fill= final_escape_result)) + geom_bar(stat= "Identity") + theme(axis.text.x = element_text(size = 14, face = "bold", color = "black")) + 
  theme(axis.title.y = element_text(size = 14, face = "bold"))+ theme(axis.text.y = element_text(size = 14, face = "bold", color ="black")) + theme(legend.text=element_text(size= 14, face = "bold"))+ theme(legend.title = element_text(size= 16, face = "bold")) +
  scale_fill_manual(values = c("gray67", "khaki2"), labels = c("no", "yes")) + #, "lightgreen", "khaki2", "forestgreen", "plum2", "firebrick3", "firebrick1", "orange", "floralwhite", "pink", "coral2", "darkmagenta", "slateblue4", "goldenrod3", "chocolate","chocolate4", "black", "grey33", "grey"))+
  ylab("percentage of samples") + xlab("") + labs(fill = "Immune escaped") + ylim(0,105)  + theme_mypub_grid() + my_theme
plot(q)
dev.off()

# stats
chisq_NR_TP <- chisq.test(Sum_df)
dfmean_chisq_NR_TP <- escape_df_NR %>% 
  group_by(interaction(final_escape_result, Timepoint)) %>% 
  dplyr::summarise(count = dplyr::n())
names(dfmean_chisq_NR_TP)[1] <- "Category_sample"
dfmean_chisq_NR_TP$Test <- "chisq"
dfmean_chisq_NR_TP$pVal <- chisq_NR_TP$p.value
write.table(dfmean_chisq_NR_TP, file = paste0(fig_dir, "/chisq_Barplot_escape_NR_per_Timepoint.txt"), quote = F, row.names = F)                                       



####################################
#  3. Immune escape per regression grade per sample
escape_df_remission <- escape_df_complete[,c("Sample", "final_escape_result")]

regression_df<- read.table("~/analysis/General_files/MEMORI_specific/WES_samples_include_for_analysis.txt", sep = "\t", header = T)
regression_df <- regression_df[,c("SampleID", "Becker_remission_grade")]
z_remission <- merge(escape_df_remission, regression_df, by.x = "Sample", by.y ="SampleID", all.x = T)
z_remission <- z_remission[,c("final_escape_result", "Becker_remission_grade")]
Sum_df <- table(z_remission)
Sum_df_perc <- scale(Sum_df, FALSE, colSums(Sum_df)) * 100
melt_Sum_df <- melt(Sum_df_perc)

pdf(paste0(fig_dir, "/Barplot_escape_remission_grade_per_sample.pdf"), width = 6, height = 5)                                       
q <- ggplot(melt_Sum_df, aes(Becker_remission_grade, value, fill= final_escape_result)) + geom_bar(stat= "Identity") + theme(axis.text.x = element_text(size = 14, face = "bold", color = "black")) + 
  theme(axis.title.y = element_text(size = 14, face = "bold"))+ theme(axis.text.y = element_text(size = 14, face = "bold", color ="black")) + theme(legend.text=element_text(size= 14, face = "bold"))+ theme(legend.title = element_text(size= 16, face = "bold")) +
  scale_fill_manual(values = c("gray67", "khaki2")) + #, "lightgreen", "khaki2", "forestgreen", "plum2", "firebrick3", "firebrick1", "orange", "floralwhite", "pink", "coral2", "darkmagenta", "slateblue4", "goldenrod3", "chocolate","chocolate4", "black", "grey33", "grey"))+
  ylab("percentage of samples") + xlab("") + labs(fill = "Immune escaped") + theme_mypub_grid() + my_theme
plot(q)
dev.off()

# stats
chisq_remission_samples <- chisq.test(Sum_df)
dfmean_remission_samples <- z_remission %>% 
  group_by(interaction(final_escape_result, Becker_remission_grade)) %>% 
  dplyr::summarise(count = dplyr::n())
names(dfmean_remission_samples)[1] <- "Category_sample"
dfmean_remission_samples$Test <- "chisq"
dfmean_remission_samples$pVal <- chisq_remission_samples$p.value
write.table(dfmean_remission_samples, file = paste0(fig_dir, "/chisq_Barplot_escape_remission_per_sample.txt"), quote = F, row.names = F)                                        

#####################################################################################################################################################
# 4. Genetic and Non-genetic immune escaper per patient and timepoint of occurence 


# Read table with immun escape mechanism and their occurence timepoint 
# Genetic immune escape is defined as HLA-LOH or mutations in HLA or B2M; Transcriptomic immune escape is defined as PDL1 overexpression
# "early occurence" is defined as "prior to treatment" (= Timepoint A), "late occurence" is defined as during treatment (= Timepoint B or C)
immune_escape_genetics_df <- read.table("~/analysis/multi_omic/input_files/MEMORI_escape_genetic_non_genetic_per_patient.txt", header = T, stringsAsFactors = F)
melt_immune_escape_genetics_df <- melt(immune_escape_genetics_df, id = "Patient")

melt_immune_escape_genetics_df <- melt_immune_escape_genetics_df[,c("variable", "value")]
Sum_df <- t(table(melt_immune_escape_genetics_df))
Sum_df_perc <- scale(Sum_df, FALSE, colSums(Sum_df)) * 100
melt_Sum_df <- melt(Sum_df_perc)
colnames(melt_Sum_df ) <- c("escape_TP", "escape_omic", "value")
melt_Sum_df$value <- as.numeric(melt_Sum_df$value)
melt_Sum_df$escape_TP <- factor(melt_Sum_df$escape_TP, levels = c("n", "late", "early"))


pdf(paste0(fig_dir, "Barplot_escape_genetic_transcript_timepoint_per_patient.pdf"), width = 6, height = 5)                                       
q <- ggplot(melt_Sum_df, aes(escape_omic, value, fill= escape_TP)) + geom_bar(stat= "Identity") + theme(axis.text.x = element_text(size = 16, face = "bold", color = "black")) + 
  theme(axis.title.y = element_text(size = 18, face = "bold"))+ theme(axis.text.y = element_text(size = 16, face = "bold", color ="black")) + theme(axis.text.x = element_blank()) + theme(legend.text=element_text(size= 16, face = "bold"))+ theme(legend.title = element_text(size= 16, face = "bold")) + 
  scale_fill_manual("I.E. timepoint", values = c("gray50", "#B2182B", "steelblue"), labels = c("n", "late", "early")) + #, "lightgreen", "khaki2", "forestgreen", "plum2", "firebrick3", "firebrick1", "orange", "floralwhite", "pink", "coral2", "darkmagenta", "slateblue4", "goldenrod3", "chocolate","chocolate4", "black", "grey33", "grey"))+
  #facet_grid(~escape_omic, scales = "free_x", space = "free_x", labeller = as_labeller(c("genetic_immune_escape" = "Genetic I.E.", "non_genetic_immune_escape" = "Transcriptomic I.E"))) +
  #theme(panel.spacing = unit(0.5, "cm"), strip.text = element_text(size = 16, face = "bold")) +
  ylab("percentage of patients") + xlab("") + labs(fill = "escape_omic") + ylim(0,105) + theme_mypub_grid() + my_theme
plot(q)
dev.off()


chisq_immune_escape_genetics <- chisq.test(Sum_df)
dfmean_immune_escape_genetics <- melt_immune_escape_genetics_df %>% 
  group_by(interaction(value, variable)) %>% 
  dplyr::summarise(count = dplyr::n())
names(dfmean_immune_escape_genetics)[1] <- "Category_sample"
dfmean_immune_escape_genetics$Test <- "chisq"
dfmean_immune_escape_genetics$pVal <- chisq_immune_escape_genetics$p.value
write.table(dfmean_immune_escape_genetics, file = paste0(fig_dir, "/chisq_Barplot_escape_genetic_transcript_timepoint_per_patient.txt"), quote = F, row.names = F)                                       


