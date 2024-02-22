###########################################################################################
### Heatmap indicating immune escape mechanisms for each sample

# prepare output dirs:
fig_dir = "~/analysis/multi_omic/immune_escape/plots/Heatmap"
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)

# Read in Mastertable of immune escape
hlaAlt.df <- read.delim('~analysis/multi_omic/immune_escape/tables/MEMORI_escape_master_file_incl_PDL1.txt')

# Select columns to plot 
hlaAlt.df.plot <- hlaAlt.df[order(hlaAlt.df$Sample),c('Sample','Patient','HLAmut','B2Mmut','HLAloh','HLAregionCNLOH', 'PDL1_immune_escape')]
hlaAlt.df.plot$HLA_LOH_detailed <- ifelse(hlaAlt.df.plot$HLAloh == "LOH" | hlaAlt.df.plot$HLAloh == "LOH (?)" | hlaAlt.df.plot$HLAregionCNLOH == "LOH (?)", "LOH", ifelse(hlaAlt.df.plot$HLAloh == "Homozygous", "Homozygous", "FALSE"))
hlaAlt.df.plot <- hlaAlt.df.plot[,c("Sample", "Patient", "HLAmut", "B2Mmut", "PDL1_immune_escape", "HLA_LOH_detailed")]
colnames(hlaAlt.df.plot) <- c("Sample", "Patient", "HLA mutation", "B2M", "PDL1", "HLA CN state")
hlaAlt.df.plot$PDL1 <- ifelse(is.na(hlaAlt.df.plot$PDL1), "FALSE", hlaAlt.df.plot$PDL1)
hlaAlt.df.plot <- hlaAlt.df.plot[complete.cases(hlaAlt.df.plot[ ,c(3,4,5)]),]

# Prepare colours and labels 
plot.df <- reshape2::melt(hlaAlt.df.plot,id=c('Sample','Patient'))
escapeColours <- setNames(alpha(c('yellow','steelblue', 'darkseagreen3', "#B2182B", 'grey80','grey80', 'gray50'), c(1,.8,.8,.8,.8,.8,.8)),
                          c('mutated','LOH','HLA_mutation','PDL1_high', 'FALSE','NA', 'Homozygous'))

plot.df$value <- factor(plot.df$value, levels=names(escapeColours))
plot.df$Sample <- gsub('MEMORI_','',plot.df$Sample)
plot.df$Timepoint <- sapply(plot.df$Sample, function(x) substr(x,7,7))
hlaAlt.df.plot$typeCol <- ifelse(grepl('NE',hlaAlt.df.plot$Patient),"#F39B7FFF","#00A087FF")


# Plot heatmap
pdf(paste0(fig_dir, "/Heatmap_immune_escape.pdf"), width = 10, height = 4)
xPat <- cumsum(table(hlaAlt.df.plot$Patient)) 
p <- ggplot(plot.df, aes(x=Sample, y=variable, fill=value)) +
  geom_tile(height = 0.85, width = 0.85) +
  coord_fixed(ratio = 5) +
  theme_mypub() +
  labs(x = NULL, y = NULL)+
  theme_mypub() +
  scale_fill_manual(values=escapeColours, labels = c("mutation", "LOH", "overexpression", "normal", "homozygous")) +
  labs(fill='') +
  geom_vline(xintercept = xPat+0.5) +
  theme(axis.text.x=element_text(size=10, face = "bold", colour="black"),
        axis.title.x=element_blank(),axis.title.y=element_blank(), axis.text.y=element_text(size=10, face = "bold", colour = "black"))  +
  scale_x_discrete(label=function(x) substr(x, 7, 7)) + theme(legend.position="bottom") + theme(legend.text = element_text(size=10, face = "bold"))
plot(p)
dev.off()

