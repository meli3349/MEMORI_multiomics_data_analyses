
#### Create Barcharts for mutational signatures 
## 1.) COSMIC signatures per sample
## 2.) COSMIC signatures for RE and NR grouped for different TP
## 3a.) COSMIC signatures for RE and NR grouped for different TP (but only signatures that are >5% in Responsiveness_TP groups)
## 3b) Line graph for selected signatures S4 and S5 over time
## 4.) Base change signatures per sample summarized to "C>A" "C>G" "C>T" "T>A" "T>C" "T>G"
## 5.) Base change signatures per sample

library(reshape2)
library(ggplot2)
library(rstatix)

# prepare output dirs:
fig_dir = "~/analysis/WES/mutations/plots/COSMIC_signatures"
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)

# My theme 
my_theme <- theme(
  plot.title = element_text(color="black", size=16, face="bold"),
  axis.title.x = element_text(color="black", size=16, face="bold"),
  axis.text.x = element_text(color="black", size=16, face="bold"),
  axis.title.y = element_text(color="black", size=16, face="bold"),
  axis.text.y = element_text(color="black", size=16, face="bold"),
  legend.title = element_text(color = "black", size = 14, face="bold"),
  legend.text = element_text(color = "black", size = 14, face="bold"),
  
)

theme_mypub_grid <- function(base_size = 14,
                             base_family = ""){
  theme_bw(base_size = base_size, 
           base_family = base_family) %+replace%
    theme(
      #panel.grid.major = element_blank(),   
      panel.grid.minor = element_blank(),   
      
      complete = TRUE
    )
}

################################################################
#### 1.)  Barchart for mutation COSMIC signatures per sample 

sig_files <- list.files(path = "~/analysis/WES/mutations/COSMIC_signatures/output_tables/Per_sample/signature_weight/all_signatures", pattern = "_signature_weights.txt", full.names = T)
samples <- gsub("\\S+(MEMORI_\\S+D\\d)\\S+","\\1",sig_files)
sign_names <- read.table(sig_files[1], header=TRUE,stringsAsFactors = F)[,1]
df <- do.call(cbind,lapply(sig_files,function(fn)read.table(fn, header = TRUE,stringsAsFactors = F)[,2]))
Mastertable_sig_files <- as.data.frame(cbind(sign_names, df), strings.as.factor = F)
colnames(Mastertable_sig_files) <- c("Signature", samples); rownames(Mastertable_sig_files) <- Mastertable_sig_files$Signature
for(sam in samples) { Mastertable_sig_files[,sam] <- as.numeric(Mastertable_sig_files[,sam])}

# Delete non-present signatures from Mastertable
Mastertable_sig_files$Signature <- paste0("S", c(1:30))
unknown <- as.vector(1-colSums(Mastertable_sig_files[,c(2:ncol(Mastertable_sig_files))])); unknown2 <- c("unknown", unknown); rownames(Mastertable_sig_files) <- Mastertable_sig_files$Signature
Mastertable_sig_files <- as.data.frame(rbind(Mastertable_sig_files, unknown2)); for(sam in samples) { Mastertable_sig_files[,sam] <- as.numeric(Mastertable_sig_files[,sam])}
Mastertable_sig_compressed <- Mastertable_sig_files[which(rowSums(Mastertable_sig_files[, c(2:ncol(Mastertable_sig_files))]) != 0), ]


Mastertable_sig_Plot <- melt(Mastertable_sig_compressed)
Mastertable_sig_Plot$Signature_ordered <- factor(Mastertable_sig_Plot$Signature, levels = c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S19", "S20", "S21","S22", "S23", "S24", "S25", "S26","S27", "S28", "S29", "S30", "unknown"))
Mastertable_sig_Plot$PatientNumber = sapply(Mastertable_sig_Plot$variable, function(x) substr(x,8, 13)) 
Mastertable_sig_Plot$Responsiveness <- sapply(Mastertable_sig_Plot$variable, function(x) substr(x,8, 9))
Mastertable_sig_Plot$Timepoint <- sapply(Mastertable_sig_Plot$variable, function(x) substr(x,14, 14))

Mastertable_sig_Plot_1 <- Mastertable_sig_Plot
Mastertable_sig_Plot_1$variable_to_plot <- Mastertable_sig_Plot_1$variable
Mastertable_sig_Plot_1$variable_to_plot <- str_replace(Mastertable_sig_Plot_1$variable_to_plot, "NE00", "NRP")
Mastertable_sig_Plot_1$variable_to_plot <- str_replace(Mastertable_sig_Plot_1$variable_to_plot, "RE00", "REP")

# 1.) Create Bargraph for each sample 
pdf(paste0(fig_dir, "/Bargraph_per_sample_SampleID_labeled.pdf"), width = 14, height = 5)
q <- ggplot(Mastertable_sig_Plot1, aes(variable_to_plot, value, fill= Signature_ordered)) + geom_bar(stat= "Identity") + theme(axis.text.x = element_text(size = 12, face = "bold", color = "black")) + 
  theme(axis.title.y = element_text(size = 14, face = "bold"))+ theme(axis.text.y = element_text(size = 14, face = "bold", color ="black")) + theme(legend.text=element_text(size= 14, face = "bold"))+ theme(legend.title = element_text(size= 16, face = "bold")) + 
  scale_fill_manual(values = c("darkblue", "steelblue", "lightblue", "aquamarine3", "lightgreen", "khaki2", "forestgreen", "plum2", "firebrick3", "firebrick1", "orange", "floralwhite", "pink", "coral2", "darkmagenta", "slateblue4", "goldenrod3", "chocolate","chocolate4", "black", "grey33", "grey"))+
  scale_x_discrete(label = function(x) substr(x, 14, 14)) + 
  facet_grid(~Responsiveness, scales = "free_x", space = "free_x", labeller = as_labeller(c("NE" = "NR", 'RE' = "RE"))) +
  theme(panel.spacing = unit(1, "cm"),
        strip.text = element_text(size = 14, face = "bold")) +
  ylab("COSMIC signature weight") + xlab("") + labs(fill = "COSMIC Signature")
plot(q)
dev.off()

  

################################################################
#### 2.) Barchart for mutation COSMIC signatures for NRP and REP during therapy

Mastertable_sig_Plot_2 <- Mastertable_sig_Plot %>% 
  select(Timepoint, Responsiveness, Signature_ordered, value) %>% 
  group_by(Timepoint, Responsiveness, Signature_ordered) %>% 
  summarise(percent = value/dplyr::n())

pdf(paste0(fig_dir, "/Bargraph_NE_RE_per_TP.pdf"), width = 6, height = 5)                                       
q <- ggplot(Mastertable_sig_Plot_2, aes(Timepoint, percent, fill= Signature_ordered)) + geom_bar(stat= "Identity") + theme(axis.text.x = element_text(size = 14, face = "bold", color = "black")) + 
  theme(axis.title.y = element_text(size = 14, face = "bold"))+ theme(axis.text.y = element_text(size = 14, face = "bold", color ="black")) + theme(legend.text=element_text(size= 14, face = "bold"))+ theme(legend.title = element_text(size= 16, face = "bold")) + 
  scale_fill_manual(values = c("darkblue", "steelblue", "lightblue", "aquamarine3", "lightgreen", "khaki2", "forestgreen", "plum2", "firebrick3", "firebrick1", "orange", "floralwhite", "pink", "coral2", "darkmagenta", "slateblue4", "goldenrod3", "chocolate","chocolate4", "black", "grey33", "grey"))+
  facet_grid(~Responsiveness, scales = "free_x", space = "free_x", labeller = as_labeller(c("NE" = "NRP", 'RE' = "REP"))) +
  theme(panel.spacing = unit(1, "cm"),
        strip.text = element_text(size = 14, face = "bold")) +
  ylab("COSMIC signature weight") + xlab("") + labs(fill = "COSMIC Signature")
plot(q)
dev.off()


################################################################
# 3a) Calculate all signatures >5% in Responsivenss_Timepoint group

Mastertable_sig_Plot_2_sign_select <- Mastertable_sig_Plot_2
Mastertable_sig_Plot_2_sign_select$Resp_Time_group <- paste0(Mastertable_sig_Plot_2_sign_select$Responsiveness, "_", Mastertable_sig_Plot_2_sign_select$Timepoint)

Resp_Time_group <- c("NE_A", "NE_B", "NE_C", "RE_A", "RE_B", "RE_C")
Signatures <- unique(Mastertable_sig_Plot_2_sign_select$Signature_ordered)

list_sign_percentage_Re_TP <- list() 
for (r in Resp_Time_group){

  df <- Mastertable_sig_Plot_2_sign_select[which(Mastertable_sig_Plot_2_sign_select$Resp_Time_group == r), ]
  list_sign_percentage <- list()
  
  for(S in Signatures){
    df_sign <- df[which(df$Signature_ordered == S), ]
    sum_sign <- sum(df_sign$percent)
    list_sign_percentage[[S]] <- sum_sign   
  }
  RE_TP_summary <- as.data.frame(cbind(list_sign_percentage))
  colnames(RE_TP_summary) <- r
  list_sign_percentage_Re_TP[[r]] <- RE_TP_summary  
}
  
All_summary <- do.call(cbind.data.frame, list_sign_percentage_Re_TP)

sign_greater_5 <- as.data.frame(subset(All_summary, NE_A > 0.05 | NE_B > 0.05 | NE_C > 0.05 | RE_A > 0.05 | RE_B > 0.05 | RE_C > 0.05))

sign_greater_5_string <- rownames(sign_greater_5)
write.table(sign_greater_5_string, file = "~/analysis/WES/mutations/COSMIC_signatures/output_tables/Per_sample/signature_weight/More_5_percent_in_Responsiveness_Timepoint_groups.txt", sep = "\t", row.names = T)


################################################################
#### 3a) Barchart for mutation COSMIC signatures for signatures >5% in NRP and REP during therapy

sig_files_select <- list.files(path = "~/analysis/WES/mutations/COSMIC_signatures/output_tables/Per_sample/signature_weight/more_5percent_signatures", pattern = "_signature_weights_more5percent_resp_TP_group.txt", full.names = T)
samples <- gsub("\\S+(MEMORI_\\S+D\\d)\\S+","\\1",sig_files_select)
sign_names_select <- read.table(sig_files_select[1], header=TRUE,stringsAsFactors = F)[,1]
df1 <- do.call(cbind,lapply(sig_files_select,function(fn)read.table(fn, header = TRUE,stringsAsFactors = F)[,2]))
Mastertable_sig_files_select <- as.data.frame(cbind(sign_names_select, df1), strings.as.factor = F)
colnames(Mastertable_sig_files_select) <- c("Signature", samples); rownames(Mastertable_sig_files_select) <- Mastertable_sig_files_select$Signature
for(sam in samples) { Mastertable_sig_files_select[,sam] <- as.numeric(Mastertable_sig_files_select[,sam])}

#Compress table by Deleting non-present signatures  
Mastertable_sig_files_select_compressed <- Mastertable_sig_files_select[which(rowSums(Mastertable_sig_files_select[, c(2:ncol(Mastertable_sig_files_select))]) != 0), ]
seletced_sig <- as.integer(stringr::str_extract(rownames(Mastertable_sig_files_select_compressed ), "\\d+"))
seletced_sig <- paste0("S", seletced_sig)
Mastertable_sig_files_select_compressed$Signature <- seletced_sig 
unknown <- as.vector(1-colSums(Mastertable_sig_files_select_compressed[,c(2:ncol(Mastertable_sig_files_select_compressed))])); unknown2 <- c("unknown", unknown); rownames(Mastertable_sig_files) <- Mastertable_sig_files$Signature
Mastertable_sig_files_select_compressed <- rbind(Mastertable_sig_files_select_compressed, unknown2); rownames(Mastertable_sig_files_select_compressed) <- Mastertable_sig_files_select_compressed$Signature

Mastertable_sig_Plot_select <- melt(Mastertable_sig_files_select_compressed, id = "Signature")
Mastertable_sig_Plot_select$Signature_ordered <- factor(Mastertable_sig_Plot_select$Signature, levels = c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S19", "S20", "S21","S22", "S23", "S24", "S25", "S26","S27", "S28", "S29", "S30", "unknown"))
Mastertable_sig_Plot_select$PatientNumber = sapply(Mastertable_sig_Plot_select$variable, function(x) substr(x,8, 13)) 
Mastertable_sig_Plot_select$Responsiveness <- sapply(Mastertable_sig_Plot_select$variable, function(x) substr(x,8, 9))
Mastertable_sig_Plot_select$Timepoint <- sapply(Mastertable_sig_Plot_select$variable, function(x) substr(x,14, 14))
Mastertable_sig_Plot_select$value <- as.numeric(Mastertable_sig_Plot_select$value)

Mastertable_sig_Plot_2_select <- Mastertable_sig_Plot_select %>% 
  select(Timepoint, Responsiveness, Signature_ordered, value) %>% 
  group_by(Timepoint, Responsiveness, Signature_ordered) %>% 
  summarise(percent = value/dplyr::n())


pdf(paste0(fig_dir, "/Bargraph_NE_RE_per_TP_only_more_5percent_signature_N.pdf"), width = 7, height = 5)                                       
q <- ggplot(Mastertable_sig_Plot_2_select, aes(Timepoint, percent, fill= Signature_ordered)) + geom_bar(stat= "Identity") + theme(axis.text.x = element_text(size = 14, face = "bold", color = "black")) + 
  theme(axis.title.y = element_text(size = 14, face = "bold"))+ theme(axis.text.y = element_text(size = 14, face = "bold", color ="black")) + theme(legend.text=element_text(size= 14, face = "bold"))+ theme(legend.title = element_text(size= 16, face = "bold")) + 
  scale_fill_manual(values = c("darkblue", "steelblue", "orange", "darkmagenta", "khaki2", "lightblue", "aquamarine3", "pink", "goldenrod3", "grey"))+
  facet_grid(~Responsiveness, scales = "free_x", space = "free_x", labeller = as_labeller(c("NE" = "NRP", 'RE' = "REP"))) +
  theme(panel.spacing = unit(1, "cm"),
        strip.text = element_text(size = 14, face = "bold")) +
  ylab("COSMIC signature weight") + xlab("") + labs(fill = "COSMIC Signature")
plot(q)
dev.off()


# Calculate stats for signatures >5% in NRP and REP during therapy
test_df <- as.data.frame(t(Mastertable_sig_files_select_compressed))
colnames(test_df) <- test_df[1,]; test_df <- test_df[-1,]
test_df$Responsiveness <- sapply(rownames(test_df), function(x) substr(x,8, 9))
test_df$Timepoint <- sapply(rownames(test_df), function(x) substr(x,14, 14))
cols.num <- c("S1", "S3", "S4", "S5", "S6", "S12", "S17", "S18", "S29", "unknown")
test_df[cols.num] <- sapply(test_df[cols.num],as.numeric)

test_RE <- test_df[which(test_df$Responsiveness == "RE"), ]
test_NE <- test_df[which(test_df$Responsiveness == "NE"), ]


# Wilcox test for individual signatures for NR and RE 
stat.test_RE <- Mastertable_sig_Plot_select[which(Mastertable_sig_Plot_select$Responsiveness == "RE"), ] %>%
  group_by(Signature_ordered) %>%
  wilcox_test(value ~ Timepoint) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

stat.test_NR <- Mastertable_sig_Plot_select[which(Mastertable_sig_Plot_select$Responsiveness == "NE"), ] %>%
  group_by(Signature_ordered) %>%
  wilcox_test(value ~ Timepoint) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()


####################################################################################################################################
## 3b) Plot line graphs for 2 most interesting signatures S4 and S5 

Mastertable_sig_Plot_select_Sig5_NE <- Mastertable_sig_Plot_select[which(Mastertable_sig_Plot_select$Responsiveness == "NE" & Mastertable_sig_Plot_select$Signature == "S5"), ]
Mastertable_sig_Plot_select_Sig4_RE <- Mastertable_sig_Plot_select[which(Mastertable_sig_Plot_select$Responsiveness == "RE" & Mastertable_sig_Plot_select$Signature == "S4"), ]

Mastertable_sig_Plot_select_Sig5_NE$Patient_TP <- sapply(Mastertable_sig_Plot_select_Sig5_NE$variable, function(x) substr(x,8,14)) 
Mastertable_sig_Plot_select_Sig4_RE$Patient_TP <- sapply(Mastertable_sig_Plot_select_Sig4_RE$variable, function(x) substr(x,8,14)) 

Mastertable_sig_Plot_select_Sig5_NE$Signature <- as.character(Mastertable_sig_Plot_select_Sig5_NE$Signature)


# Calculate mean for multi-region samples
mean_Mastertable_sig_Plot_select_Sig5_NE = aggregate(Mastertable_sig_Plot_select_Sig5_NE, by = list(Mastertable_sig_Plot_select_Sig5_NE$Patient_TP), FUN = mean)
mean_Mastertable_sig_Plot_select_Sig4_RE = aggregate(Mastertable_sig_Plot_select_Sig4_RE, by = list(Mastertable_sig_Plot_select_Sig4_RE$Patient_TP), FUN = mean)
mean_Mastertable_sig_Plot_select_Sig5_NE$Signature <- "S5"; mean_Mastertable_sig_Plot_select_Sig5_NE$Timepoint <- sapply(mean_Mastertable_sig_Plot_select_Sig5_NE$Group.1, function(x) substr(x,7,7)); mean_Mastertable_sig_Plot_select_Sig5_NE$PatientNumber <- sapply(mean_Mastertable_sig_Plot_select_Sig5_NE$Group.1, function(x) substr(x,1,6))
mean_Mastertable_sig_Plot_select_Sig4_RE$Signature <- "S4"; mean_Mastertable_sig_Plot_select_Sig4_RE$Timepoint <- sapply(mean_Mastertable_sig_Plot_select_Sig4_RE$Group.1, function(x) substr(x,7,7)); mean_Mastertable_sig_Plot_select_Sig4_RE$PatientNumber <- sapply(mean_Mastertable_sig_Plot_select_Sig4_RE$Group.1, function(x) substr(x,1,6))


pdf(paste0(fig_dir, "/COSMIC_S4_in_RE_during_treatment.pdf"), width = 5, height = 5)
p <- ggplot(data=mean_Mastertable_sig_Plot_select_Sig4_RE, aes(x=Timepoint, y=value, group = PatientNumber)) +
  geom_line(colour = 'orange', size = 1.0)+
  geom_point(colour = 'orange', size = 1.5)  +
  theme_mypub_grid() + 
  labs(x='', y="S4 weight") + ggtitle("COSMIC S4 in REPs") + my_theme #+ stat_compare_means(comparisons=list(c('A','B'),c('A','C'), c('B','C')), method = "wilcox.test", size = 4.5, face="bold")
plot(p)
dev.off()

pdf(paste0(fig_dir, "/COSMIC_S5_in_NR_during_treatment.pdf"), width = 5, height = 5)
p <- ggplot(data=mean_Mastertable_sig_Plot_select_Sig5_NE, aes(x=Timepoint, y=value, group = PatientNumber)) +
  geom_line(colour = "darkmagenta", size = 1.0)+
  geom_point(colour = "darkmagenta", size = 1.5)  +
  theme_mypub_grid() + 
  labs(x='', y="S5 weight") + ggtitle("COSMIC S5 in NRPs") + my_theme 
plot(p)
dev.off()



################################################################
#### 6.) SNVs from phylogenetic tree were classified as treatment-naive SNVs, CTx-induced and RCTx-induced SNVs 
#### COSMIC SNV signatures were run for each SNV category


sig_base_files <- list.files(path = "Bargraph_per_SNV_inducer_of_patients_from_phylotrees.pdf", pattern = "_signature_tumours_summarized.txt", full.names = T)
sign_base_names <- read.table(sig_base_files[1], header=TRUE,stringsAsFactors = F)[,1]
df_base_sig <- do.call(cbind,lapply(sig_base_files,function(fn)read.table(fn, header = TRUE,stringsAsFactors = F)[,2]))
Mastertable_sig_base_files <- as.data.frame(cbind(sign_base_names, df_base_sig), strings.as.factor = F)
colnames(Mastertable_sig_base_files) <- c("Signature_bases", "chemo_induced", "radio_induced", "treatment_naive"); rownames(Mastertable_sig_base_files) <- Mastertable_sig_base_files$Signature_bases
#Mastertable_sig_base_files[,c("chemo_induced", "radio_induced", "treatment_naive")] <- as.numeric(Mastertable_sig_base_files[,c("chemo_induced", "radio_induced", "treatment_naive")])
SNV_inducer <- c("chemo_induced", "radio_induced", "treatment_naive")
for(sam in SNV_inducer) { Mastertable_sig_base_files[,sam] <- as.numeric(Mastertable_sig_base_files[,sam])}

Mastertable_sig_base_Plot <- melt(Mastertable_sig_base_files)
Mastertable_sig_base_Plot$Signature_ordered <- factor(Mastertable_sig_base_Plot$Signature_bases, levels = c("base_change_T>G", "base_change_T>C", "base_change_T>A", "base_change_C>T","base_change_C>G","base_change_C>A"))
Mastertable_sig_base_Plot$variable <- factor(Mastertable_sig_base_Plot$variable,levels = c("treatment_naive", "chemo_induced", "radio_induced"))

pdf(paste0(fig_dir, "/Bargraph_per_SNV_inducer_of_patients_from_phylotrees.pdf"), width = 5, height = 7)
q <- ggplot(Mastertable_sig_base_Plot, aes(variable, value, fill= Signature_ordered)) + geom_bar(stat= "Identity") + theme(axis.text.x = element_text(angle = 90, size = 14, face = "bold", color = "black")) + 
  theme(axis.title.y = element_text(size = 14, face = "bold"))+ theme(axis.text.y = element_text(size = 14, face = "bold", color ="black")) + theme(legend.text=element_text(size= 16, face = "bold"))+ theme(legend.title = element_text(size= 16, face = "bold")) + 
  scale_fill_manual(values = c("darkblue", "steelblue", "lightblue", "aquamarine3", "khaki2", "mediumvioletred"), labels = c("T>G", "T>C", "T>A", "C>T", "C>G", "C>A"))+
  theme(panel.spacing = unit(1, "cm"), 
        strip.text = element_text(size = 14, face = "bold")) +
  ylab("Signature weight") + xlab("") + 
  scale_x_discrete(labels=c("treatment_naive" = "pre-treatment", "chemo_induced" = "CTx induced", "radio_induced" = "RCTx induced")) + 
  labs(fill = "")
plot(q)
dev.off()
