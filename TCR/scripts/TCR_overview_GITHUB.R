
################################################
# TCR overview graphs
################################################
require(ggplot2); require(ggpubr); require(scales); require(reshape2); library(plyr); require(scales)


# prepare output dirs:
fig_dir = "~/analysis/TCR/plots/"
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)

# Read in TCR summary file 
TCR_summary_file <- read.table("TCR_summary.txt", header = T, stringsAsFactors = F)

# Themes
theme_mypub <- function(base_size = 14,
                        base_family = ""){
  theme_bw(base_size = base_size, 
           base_family = base_family) %+replace%
    theme(
      panel.grid.major = element_blank(),   
      panel.grid.minor = element_blank(),   
      
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



# Import TCR data with UniqueDCRsPostCollapsing	and TotalDCRsPostCollapsing

df <- final_summary[,c("SampleID", "UniqueDCRsPostCollapsing", "TotalDCRsPostCollapsing")]
df$SampleID2 <- sapply(df$SampleID, function(x) substr(x,1,22))
df$Timepoint <- sapply(df$SampleID, function(x) substr(x,13,13))
df$Responsiveness <- sapply(df$SampleID, function(x) substr(x,7,8))
df$Patient <- sapply(df$SampleID, function(x) substr(x,7,12))
df$TCR_chain <-  sapply(df$SampleID, function(x) substr(x, 24, nchar(df$SampleID)))


df_alpha <- df[which(df$TCR_chain == "alpha"),]
df_beta <- df[which(df$TCR_chain == "beta"),]
list_df <- list(df_alpha, df_beta)


## Plot Total TCR counts during treatment 

for (l in list_df){
  TCR_chain <- unique(l$TCR_chain)
  
  # For each timepoint
  pdf(paste0(fig_dir, "TotalDCRs_", TCR_chain, "_per_timepoint.pdf"), height = 4.5, width = 4.5)
  p <- ggplot2.violinplot(data=l, xName="Timepoint",yName="TotalDCRsPostCollapsing",fill = "Responsiveness",
                          groupName="Timepoint",
                          groupColors= c('#e63e00','darkseagreen3','steelblue'), showLegend=FALSE,
                          backgroundColor="white", xtitle= "Timepoint", ytitle=paste0("Total TCRs ", TCR_chain), 
                          mainTitle=paste0("Total TCRs ", TCR_chain),
                          addDot=T, dotSize=0.5,
                          legendPosition="bottom") + 
    scale_y_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x)), limits = c(1, 100500))  + 
    stat_compare_means(comparisons=list(c('A','B'),c('A','C'), c('B','C')), method = "wilcox.test", size = 3.7)
  q <- p + my_theme + theme(plot.title = element_text(color="black", size=14, face="bold")) + theme(axis.title.x = element_text(color="black", size=16, face="bold"))+ theme(axis.title.y = element_text(color="black", size=14, face="bold")) 
  plot(q)
  dev.off()
  
  
  # For NRPs and REPs at each timepoint
  pdf(paste0(fig_dir, "TotalDCRs_", TCR_chain, "_per_responsiveness_timepoint.pdf"), height = 4.5, width = 4.5)
  p <- ggplot2.violinplot(data=l, xName="Timepoint",yName="TotalDCRsPostCollapsing",fill = "Responsiveness",
                          groupName="Timepoint",
                          groupColors= c('#e63e00','darkseagreen3','steelblue'), showLegend=FALSE,
                          backgroundColor="white", xtitle= "Timepoint", ytitle=paste0("Unique TCRs ", TCR_chain), 
                          mainTitle=paste0("Unique TCRs ", TCR_chain),
                          addDot=T, dotSize=0.5,
                          legendPosition="bottom") + 
    facet_grid(~Responsiveness , scales = "free_x", space = "free_x", labeller = as_labeller(c('NE' = "NRP", 'RE' = "REP"))) +
    theme(panel.spacing = unit(0.5, "cm"),
          strip.text = element_text(size = 12, face = "bold"))  +
    scale_y_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x)), limits = c(1, 100000))  + 
    stat_compare_means(comparisons=list(c('A','B'),c('A','C'), c('B','C')), method = "wilcox.test", size = 3.7)
  q <- p + my_theme + theme(plot.title = element_text(color="black", size=14, face="bold")) + theme(axis.title.x = element_text(color="black", size=16, face="bold"))+ theme(axis.title.y = element_text(color="black", size=14, face="bold")) 
  plot(q)
  dev.off()
  
}