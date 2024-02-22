
#################################################################
## Plotting KEGG pathway analysis for comparisions of interest ##
#################################################################
# All plots will include KEGG pathways with p< 0.084 and intersection line at siginificance level (p = 0.05) will delineate significant pathways

library(ggplot2)
library(ggrepel)

# prepare output dirs:
fig_dirKEGG = "~/analysis/RNA_Seq/pathway_analyses/KEGG_analysis/Plots"
dir.create(data_dirKEGG, showWarnings=FALSE, recursive=TRUE)

# Read in table with pathways of interest for OAC
pathways_of_interest <- read.table("~/analysis/General_files/external_files/KEGG_pathways_of_interest.txt", header = T, sep = "\t", stringsAsFactors = F)


# Create theme for all plots 
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

### for Responder TP A versus TP B
y <- read.delim("~/analysis/RNA_Seq/pathway_analyses/KEGG_analysis/Tables/KEGG_table_res_ReA_ReB.p.adj.txt", header = T, stringsAsFactors = F)
y <- y[which(y$p.adjust < 0.084), ]
y$negat_log10_p_adj <- -log10(y$p.adjust)
y$enrichment_score_color <- ifelse(y$enrichmentScore < 0, "downregulated", "upregulated")
y <- merge(y, pathways_of_interest, by.x = "Description", by.y = "KEGG.pathways.of.interest", all.x = T)
y$colour <- ifelse(y$enrichmentScore < 0, '#afa0de', '#02818a')
# No significant pathways to plot.


### for REP TP B versus REP TP C
  y <- read.delim("~/analysis/RNA_Seq/pathway_analyses/KEGG_analysis/Tables/KEGG_table_res_ReB_ReC.p.adj.txt", header = T, stringsAsFactors = F)
  y <- y[which(y$p.adjust < 0.084), ]
  y$negat_log10_p_adj <- -log10(y$p.adjust)
  y$enrichment_score_color <- ifelse(y$enrichmentScore < 0, "downregulated", "upregulated")
  y <- merge(y, pathways_of_interest, by.x = "Description", by.y = "KEGG.pathways.of.interest", all.x = T)
  y$colour <- ifelse(y$enrichmentScore < 0, '#386890',  '#80a980')
  
  
  pdf(paste0(fig_dirKEGG, "/res_ReB_ReC.pdf"), width = 6, height = 5)
  p <- ggplot(y, aes(x=enrichmentScore, y=negat_log10_p_adj, color=enrichment_score_color)) + geom_hline(yintercept=1.301, linetype='dotted', alpha = 0.7) +
    geom_point(size=3, colour = y$colour, alpha = 0.4) + scale_color_manual(values=c('#386890', 'darkseagreen3'), labels = c("upregulated in TP C", "upregulated in TP B")) +
  ggtitle("REP timepoint B versus REP timepoint C") + xlab("KEGG enrichment score") + ylab("-log10 adj p-value") + 
    geom_text_repel(aes(label = Abbreviated.pathway.of.interest, colour="white", segment.colour="black"),
                    max.overlaps = 20,
                    box.padding = 0.35, 
                    fontface = "bold",
                    size = 5,
                    colour = 'black',
                    segment.color = 'grey50') +
    theme(legend.title = element_blank()) + theme(legend.position="bottom") + theme_mypub_grid() + my_theme + 
    ylim(1.1, 1.6) + scale_x_continuous(breaks = c(-0.6, -0.3, 0, 0.3, 0.6), limits=c(-.7,0.7))
  plot(p) 
  dev.off()
  
   
  ### for Non-Responder TP B versus Non-Responder TP C
  
  y <- read.delim("~/analysis/RNA_Seq/pathway_analyses/KEGG_analysis/Tables/KEGG_table_res_NeB_NeC.p.adj.txt", header = T, stringsAsFactors = F)
  y <- y[which(y$p.adjust <= 0.084), ]
   y$negat_log10_p_adj <- -log10(y$p.adjust)
   y$enrichment_score_color <- ifelse(y$enrichmentScore < 0, "downregulated", "upregulated")
   y <- merge(y, pathways_of_interest, by.x = "Description", by.y = "KEGG.pathways.of.interest", all.x = T)
   y$colour <- ifelse(y$enrichmentScore < 0, '#386890',  '#80a980')
  
  pdf(paste0(fig_dirKEGG, "/res_NeB_NeC.pdf"), width = 6, height = 5)
  p <- ggplot(y, aes(x=enrichmentScore, y=negat_log10_p_adj, color=enrichment_score_color)) + geom_hline(yintercept=1.301, linetype='dotted', alpha = 0.7) +
    geom_point(size=3, colour = y$colour, alpha = 0.4) + scale_color_manual(values=c('#b9adde', '#02818a'), labels = c("upregulated in TP C", "upregulated in TP B")) +
    ggtitle("NEP timepoint B vs NEP timepoint C") + xlab("KEGG enrichment score") + ylab("-log10 adj p-value") + 
    geom_text_repel(aes(label = Abbreviated.pathway.of.interest, colour="white", segment.colour="black"),
                    max.overlaps = 20,
                     box.padding = 0.35, 
                     fontface = "bold",
                     size = 5,
                     colour = 'black',
                     segment.color = 'grey50') +
   theme(legend.title = element_blank()) + theme(legend.position="bottom") + theme_mypub_grid() + my_theme + 
    ylim(1.1, 1.8) + scale_x_continuous(breaks = c(-0.6, -0.3, 0, 0.3, 0.6), limits=c(-.7,0.7))
  plot(p) 
  dev.off()
  
  ### for Non-Responder TP A versus TP C
  y <- read.delim("~/analysis/RNA_Seq/pathway_analyses/KEGG_analysis/Tables/KEGG_table_res_NeA_NeC.p.adj.txt", header = T, stringsAsFactors = F)
  y <- y[which(y$p.adjust <= 0.084), ]
  y$negat_log10_p_adj <- -log10(y$p.adjust)
  y$enrichment_score_color <- ifelse(y$enrichmentScore < 0, "downregulated", "upregulated")
  y <- merge(y, pathways_of_interest, by.x = "Description", by.y = "KEGG.pathways.of.interest", all.x = T)
  y$colour <- ifelse(y$enrichmentScore < 0, '#386890', '#e63e00')
  
  pdf(paste0(fig_dirKEGG, "/res_NeA_NeC.pdf"), width = 6, height = 5)
  p <- ggplot(y, aes(x=enrichmentScore, y=negat_log10_p_adj, color=enrichment_score_color)) + geom_hline(yintercept=1.301, linetype='dotted', alpha = 0.7) +
    geom_point(size=3, colour = y$colour, alpha = 0.4) + scale_color_manual(values=c('#b9adde', '#02818a'), labels = c("upregulated in TP C", "upregulated in TP A")) +
    ggtitle("NEP timepoint A vs NEP timepoint C") + xlab("KEGG enrichment score") + ylab("-log10 adj p-value") + 
    geom_text_repel(aes(label = Abbreviated.pathway.of.interest, colour="white", segment.colour="black"),
                    max.overlaps = 20,
                     box.padding = 0.35, 
                     fontface = "bold",
                     size = 5,
                     colour = 'black',
                     segment.color = 'grey50') +
     theme(legend.title = element_blank()) + theme(legend.position="bottom") + theme_mypub_grid() + my_theme + 
    ylim(1.1, 1.65) + scale_x_continuous(breaks = c(-0.6, -0.3, 0, 0.3, 0.6), limits=c(-.7,0.7))
  plot(p) 
  dev.off()
  
  
  ### for Responder TP A versus TP C
  y <- read.delim("~/analysis/RNA_Seq/pathway_analyses/KEGG_analysis/Tables/KEGG_table_res_ReA_ReC.p.adj.txt", header = T, stringsAsFactors = F)
  y <- y[which(y$p.adjust <= 0.084), ]
  y$negat_log10_p_adj <- -log10(y$p.adjust)
  y$enrichment_score_color <- ifelse(y$enrichmentScore < 0, "downregulated", "upregulated")
  y <- merge(y, pathways_of_interest, by.x = "Description", by.y = "KEGG.pathways.of.interest", all.x = T)
  y$colour <- ifelse(y$enrichmentScore < 0, '#386890', '#e63e00')
  
  pdf(paste0(fig_dirKEGG, "/res_ReA_ReC.pdf"), width = 6, height = 5)
  p <- ggplot(y, aes(x=enrichmentScore, y=negat_log10_p_adj, color=enrichment_score_color)) + geom_hline(yintercept=1.301, linetype='dotted', alpha = 0.7) +
    geom_point(size=3, colour = y$colour, alpha = 0.4) + scale_color_manual(values=c('#b9adde', '#02818a'), labels = c("upregulated in TP C", "upregulated in TP A")) +
    ggtitle("REP timepoint A vs REP timepoint C") + xlab("KEGG enrichment score") + ylab("-log10 adj p-value") + 
    geom_text_repel(aes(label = Abbreviated.pathway.of.interest, colour="white", segment.colour="black"),
                    max.overlaps = 20,
                    box.padding = 0.35, 
                     fontface = "bold",
                     size = 5,
                     colour = 'black',
                     segment.color = 'grey50') +
    theme(legend.title = element_blank()) + theme(legend.position="bottom") + theme_mypub_grid() + my_theme + 
    ylim(1.1, 2) + scale_x_continuous(breaks = c(-0.6, -0.3, 0, 0.3, 0.6), limits=c(-.7,0.7))
  plot(p) 
  dev.off()


  ### for TP A (REP $ NEP) versus Responder TP C
  y <- read.delim("~/analysis/RNA_Seq/pathway_analyses/KEGG_analysis/Tables/KEGG_table_res_A_ReC.p.adj.txt", header = T, stringsAsFactors = F)
  y <- y[which(y$p.adjust <= 0.084), ]
  y$negat_log10_p_adj <- -log10(y$p.adjust)
  y$enrichment_score_color <- ifelse(y$enrichmentScore < 0, "downregulated", "upregulated")
  y <- merge(y, pathways_of_interest, by.x = "Description", by.y = "KEGG.pathways.of.interest", all.x = T)
  y$colour <- ifelse(y$enrichmentScore < 0, '#386890', '#e63e00')
  
  pdf(paste0(fig_dirKEGG, "/res_A_ReC.pdf"), width = 7, height = 6)
  p <- ggplot(y, aes(x=enrichmentScore, y=negat_log10_p_adj, color=enrichment_score_color)) + geom_hline(yintercept=1.301, linetype='dotted', alpha = 0.7) +
    geom_point(size=3, colour = y$colour, alpha = 0.4) + scale_color_manual(values=c('#b9adde', '#02818a'), labels = c("upregulated in TP C", "upregulated in TP A")) +
    ggtitle("Timepoint A vs REP Timepoint C") + xlab("KEGG enrichment score") + ylab("-log10 adj p-value") + 
    geom_text_repel(aes(label = Abbreviated.pathway.of.interest, colour="white", segment.colour="black"),
                    max.overlaps = Inf,
                     box.padding = 0.1, 
                     fontface = "bold",
                     size = 3,
                     colour = 'black',
                     segment.color = 'grey50') +
    theme(legend.title = element_blank()) + theme(legend.position="bottom") + theme_mypub_grid() + my_theme + 
    ylim(1.1, 2.0) + scale_x_continuous(breaks = c(-0.6, -0.3, 0, 0.3, 0.6), limits=c(-.7,0.7))
  plot(p) 
  dev.off()
  
  ### for TP A (REP $ NEP) versus Non-Responder TP C
  y <- read.delim("~/analysis/RNA_Seq/pathway_analyses/KEGG_analysis/Tables/KEGG_table_res_A_NeC.p.adj.txt", header = T, stringsAsFactors = F)
  y <- y[which(y$p.adjust <= 0.084), ]
  y$negat_log10_p_adj <- -log10(y$p.adjust)
  y$enrichment_score_color <- ifelse(y$enrichmentScore < 0, "downregulated", "upregulated")
  y <- merge(y, pathways_of_interest, by.x = "Description", by.y = "KEGG.pathways.of.interest", all.x = T)
  y$colour <- ifelse(y$enrichmentScore < 0, '#386890', '#e63e00')
  
  pdf(paste0(fig_dirKEGG, "/res_A_NeC.pdf"), width = 6, height = 5)
  p <- ggplot(y, aes(x=enrichmentScore, y=negat_log10_p_adj, color=enrichment_score_color)) + geom_hline(yintercept=1.301, linetype='dotted', alpha = 0.7) +
    geom_point(size=3, colour = y$colour, alpha = 0.4) + scale_color_manual(values=c('#b9adde', '#02818a'), labels = c("upregulated in TP C", "upregulated in TP A")) +
    ggtitle("Timepoint A vs NRP Timepoint C") + xlab("KEGG enrichment score") + ylab("-log10 adj p-value") + 
    geom_text_repel(aes(label = Abbreviated.pathway.of.interest, colour="white", segment.colour="black"),
                    max.overlaps = Inf,
                     box.padding = 0.1, 
                     fontface = "bold",
                     size = 3,
                     colour = 'black',
                     segment.color = 'grey50') +
    theme(legend.title = element_blank()) + theme(legend.position="bottom") + theme_mypub_grid() + my_theme + 
    ylim(1.1, 2.0) + scale_x_continuous(breaks = c(-0.6, -0.3, 0, 0.3, 0.6), limits=c(-.7,0.7))
  plot(p) 
  dev.off()
  
  ### for TP B (REP $ NEP) versus Responder TP C
  y <- read.delim("~/analysis/RNA_Seq/pathway_analyses/KEGG_analysis/Tables/KEGG_table_res_B_ReC.p.adj.txt", header = T, stringsAsFactors = F)
  y <- y[which(y$p.adjust <= 0.084), ]
  y$negat_log10_p_adj <- -log10(y$p.adjust)
  y$enrichment_score_color <- ifelse(y$enrichmentScore < 0, "downregulated", "upregulated")
  y <- merge(y, pathways_of_interest, by.x = "Description", by.y = "KEGG.pathways.of.interest", all.x = T)
  y$colour <- ifelse(y$enrichmentScore < 0, '#386890', '#80a980')
  
  pdf(paste0(fig_dirKEGG, "/res_B_ReC.pdf"), width = 7, height = 6)
  p <- ggplot(y, aes(x=enrichmentScore, y=negat_log10_p_adj, color=enrichment_score_color)) + geom_hline(yintercept=1.301, linetype='dotted', alpha = 0.7) +
    geom_point(size=3, colour = y$colour, alpha = 0.4) + scale_color_manual(values=c('#b9adde', '#02818a'), labels = c("upregulated in TP C", "upregulated in TP B")) +
    ggtitle("Timepoint B vs REP Timepoint C") + xlab("KEGG enrichment score") + ylab("-log10 adj p-value") + 
    geom_text_repel(aes(label = Abbreviated.pathway.of.interest, colour="white", segment.colour="black"),
                    max.overlaps = Inf,
                    box.padding = 0.1, 
                    fontface = "bold",
                    size = 3,
                    colour = 'black',
                    segment.color = 'grey50') +
    theme(legend.title = element_blank()) + theme(legend.position="bottom") + theme_mypub_grid() + my_theme + 
    ylim(1.1, 1.8) + scale_x_continuous(breaks = c(-0.6, -0.3, 0, 0.3, 0.6), limits=c(-.7,0.7))
  plot(p) 
  dev.off()
  
  ### for TP B (REP $ NEP) versus Non-Responder TP C
  y <- read.delim("~/analysis/RNA_Seq/pathway_analyses/KEGG_analysis/Tables/KEGG_table_res_B_NeC.p.adj.txt", header = T, stringsAsFactors = F)
  y <- y[which(y$p.adjust <= 0.084), ]
  y$negat_log10_p_adj <- -log10(y$p.adjust)
  y$enrichment_score_color <- ifelse(y$enrichmentScore < 0, "downregulated", "upregulated")
  y <- merge(y, pathways_of_interest, by.x = "Description", by.y = "KEGG.pathways.of.interest", all.x = T)
  y$colour <- ifelse(y$enrichmentScore < 0, '#386890', '#80a980')
  
  pdf(paste0(fig_dirKEGG, "/res_B_NeC.pdf"), width = 6, height = 5)
  p <- ggplot(y, aes(x=enrichmentScore, y=negat_log10_p_adj, color=enrichment_score_color)) + geom_hline(yintercept=1.301, linetype='dotted', alpha = 0.7) +
    geom_point(size=3, colour = y$colour, alpha = 0.4) + scale_color_manual(values=c('#b9adde', '#02818a'), labels = c("upregulated in TP C", "upregulated in TP A")) +
    ggtitle("Timepoint B vs NRP Timepoint C") + xlab("KEGG enrichment score") + ylab("-log10 adj p-value") + 
    geom_text_repel(aes(label = Abbreviated.pathway.of.interest, colour="white", segment.colour="black"),
                    max.overlaps = Inf,
                    box.padding = 0.1, 
                    fontface = "bold",
                    size = 4,
                    colour = 'black',
                    segment.color = 'grey50') +
    theme(legend.title = element_blank()) + theme(legend.position="bottom") + theme_mypub_grid() + my_theme + 
    ylim(1.1, 1.8) + scale_x_continuous(breaks = c(-0.6, -0.3, 0, 0.3, 0.6), limits=c(-.7,0.7))
  plot(p) 
  dev.off()
  
  
