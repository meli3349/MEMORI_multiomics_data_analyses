# Script to plot CNA segments of multi-timepoint/multi-region tumour samples as a whole genome heatmap
library(stringr)

# Options
fig_dir = "analysis/WES/CNA/plots"
dataset_dir = "analysis/WES/CNA/datasets"

# prepare output dirs:
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(dataset_dir, showWarnings=FALSE, recursive=TRUE)

# Functions #
# Classify CNAs in terms of gain/loss
classifyCNAs <- function(segs) {
  classes <- c()
  segs <- segs[which(segs$chromosome %in% paste0('chr',c(1:22))),]
  segs <- segs[which(!is.na(segs$A) & !is.na(segs$B)),]
  for(s in c(1:nrow(segs))) {
    curseg <- segs[s,]
    if(curseg$B==0) { #If monosomy
      if(curseg$A==1) { classes <- c(classes,'Mono') # Classic monosomy
      } else if(curseg$A==2) { classes <- c(classes,'cnLOH') # cnLOH
      } else if(curseg$A>2 | curseg$A==0) { classes <- c(classes,'Del') } # Other Deletion
    } else {
      if(curseg$CNt==3) { classes <- c(classes,'Gain') # Trisomy
      } else if(curseg$CNt>=4) { classes <- c(classes,'HighGain') # Tetrasomy and above
      } else { classes <- c(classes,'None')} # No CNA
    }
  }
  segs$CNA <- classes
  return(segs)
}


# Setup ####
chroms <- paste0('chr',c(1:22))
chrInfo <- read.csv("/analysis/General_files/external_files/centromerePositions_hg38.csv", stringsAsFactors = F);chrInfo <- chrInfo[c(1:22),]
matTab <- data.frame(chrInfo[c("chromosome","plotStart")])


# Classify CNAs in terms of gain/loss for included MEMORI samples
selected_samples <- read.table("/analysis/General_files/MEMORI_specific/WES_samples_include_for_analysis.txt", sep = "\t", header = T)
selected_samples_Seq <- selected_samples[which(selected_samples$PASS_Sequenza_best_fitting_ploidy == "y"), ]
samnames <- selected_samples_Seq$SampleID

seglist <- list()
for(sam in samnames) {
  # read in CN segments files (output from Sequenza) 
  tmpseg <- read.table(paste0("/analysis/WES/input_files/Sequenza/best_fitting_ploidies_per_sample/",sam,"_segments.txt"), header=T)
  segs <- classifyCNAs(tmpseg)
  seglist[[sam]] <- segs
}

saveRDS(seglist,file=paste0(dataset_dir, "/segement_files_for_segment_plot.R"))


TPcols <- c(A='#e63e00',B='darkseagreen3',C='steelblue')


# load CNA data
MEMORI_NE0005A_T2_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0005A_T2_B1_D1_segments.txt", header=T)
MEMORI_NE0005B_T2_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0005B_T2_B1_D1_segments.txt", header=T)
MEMORI_NE0005C_T2_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0005C_T2_C1_D1_segments.txt", header=T)

MEMORI_NE0006A_T2_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0006A_T2_B1_D1_segments.txt", header=T)
MEMORI_NE0006B_T1_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0006B_T1_B1_D1_segments.txt", header=T)
MEMORI_NE0006B_T3_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0006B_T3_C1_D1_segments.txt", header=T)

MEMORI_NE0007A_T1_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0007A_T1_B1_D1_segments.txt", header=T)
MEMORI_NE0007B_T1_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0007B_T1_B1_D1_segments.txt", header=T)

MEMORI_NE0008A_T4_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0008A_T4_B1_D1_segments.txt", header=T)
MEMORI_NE0008B_T2_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0008B_T2_B1_D1_segments.txt", header=T)

MEMORI_NE0009A_T1_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0009A_T1_B1_D1_segments.txt", header=T)
MEMORI_NE0009B_T1_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0009B_T1_B1_D1_segments.txt", header=T)
MEMORI_NE0009C_T1_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0009C_T1_C1_D1_segments.txt", header=T)

MEMORI_NE0010A_T1_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0010A_T1_C1_D1_segments.txt", header=T)
MEMORI_NE0010B_T2_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0010B_T2_C1_D1_segments.txt", header=T)
MEMORI_NE0010B_T3_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0010B_T3_C1_D1_segments.txt", header=T)
MEMORI_NE0010C_T2_M1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0010C_T2_M1_D1_segments.txt", header=T)

MEMORI_NE0012A_T1_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0012A_T1_B1_D1_segments.txt", header=T)
MEMORI_NE0012B_T1_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0012B_T1_C1_D1_segments.txt", header=T)
MEMORI_NE0012C_T1_M1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0012C_T1_M1_D1_segments.txt", header=T)

MEMORI_NE0020A_T1_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0020A_T1_B1_D1_segments.txt", header=T)

MEMORI_NE0021A_T2_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0021A_T2_C1_D1_segments.txt", header=T)
MEMORI_NE0021B_T1_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0021B_T1_B1_D1_segments.txt", header=T)
MEMORI_NE0021C_T1_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0021C_T1_C1_D1_segments.txt", header=T)

MEMORI_NE0022A_T2_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0022A_T2_B1_D1_segments.txt", header=T)
MEMORI_NE0022A_T2_B1_D2 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0022A_T2_B1_D2_segments.txt", header=T)
MEMORI_NE0022B_T1_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0022B_T1_C1_D1_segments.txt", header=T)
MEMORI_NE0022C_T1_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0022C_T1_C1_D1_segments.txt", header=T)

MEMORI_RE0002A_T4_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0002A_T4_C1_D1_segments.txt", header=T)
MEMORI_RE0002B_T1_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0002B_T1_B1_D1_segments.txt", header=T)

MEMORI_RE0003A_T1_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0003A_T1_C1_D1_segments.txt", header=T)
MEMORI_RE0003B_T2_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0003B_T2_C1_D1_segments.txt", header=T)
MEMORI_RE0003C_T1_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0003C_T1_C1_D1_segments.txt", header=T)

MEMORI_RE0004A_T2_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0004A_T2_B1_D1_segments.txt", header=T)
MEMORI_RE0004B_T3_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0004B_T3_B1_D1_segments.txt", header=T)
MEMORI_RE0004C_T1_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0004C_T1_C1_D1_segments.txt", header=T)

MEMORI_RE0005A_T1_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0005A_T1_C1_D1_segments.txt", header=T)
MEMORI_RE0005B_T1_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0005B_T1_B1_D1_segments.txt", header=T)
MEMORI_RE0005C_T1_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0005C_T1_C1_D1_segments.txt", header=T)

MEMORI_RE0006A_T3_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0006A_T3_B1_D1_segments.txt", header=T)
MEMORI_RE0006B_T1_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0006B_T1_B1_D1_segments.txt", header=T)
MEMORI_RE0006C_T3_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0006C_T3_C1_D1_segments.txt", header=T)
MEMORI_RE0006C_T2_M1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0006C_T2_M1_D1_segments.txt", header=T)

MEMORI_RE0007A_T3_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0007A_T3_C1_D1_segments.txt", header=T)
MEMORI_RE0007C_T1_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0007C_T1_C1_D1_segments.txt", header=T)

MEMORI_RE0008A_T1_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0008A_T1_B1_D1_segments.txt", header=T)
MEMORI_RE0008B_T3_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0008B_T3_B1_D1_segments.txt", header=T)
MEMORI_RE0008C_T1_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0008C_T1_C1_D1_segments.txt", header=T)
MEMORI_RE0008C_T2_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0008C_T2_B1_D1_segments.txt", header=T)
MEMORI_RE0008C_T3_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0008C_T3_C1_D1_segments.txt", header=T)

MEMORI_RE0009A_T1_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0009A_T1_B1_D1_segments.txt", header=T)

MEMORI_RE0010A_T1_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0010A_T1_C1_D1_segments.txt", header=T)
MEMORI_RE0010A_T2_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0010A_T2_C1_D1_segments.txt", header=T)
MEMORI_RE0010A_T4_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0010A_T4_C1_D1_segments.txt", header=T)
MEMORI_RE0010B_T1_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0010B_T1_B1_D1_segments.txt", header=T)
MEMORI_RE0010B_T2_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0010B_T2_C1_D1_segments.txt", header=T)
MEMORI_RE0010B_T3_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0010B_T3_C1_D1_segments.txt", header=T)
MEMORI_RE0010C_T1_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0010C_T1_C1_D1_segments.txt", header=T)


MEMORI_RE0011A_T1_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0011A_T1_B1_D1_segments.txt", header=T)
MEMORI_RE0011A_T3_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0011A_T3_B1_D1_segments.txt", header=T)
MEMORI_RE0011B_T1_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0011B_T1_C1_D1_segments.txt", header=T)
MEMORI_RE0011B_T2_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0011B_T2_C1_D1_segments.txt", header=T)
MEMORI_RE0011B_T3_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0011B_T3_C1_D1_segments.txt", header=T)
MEMORI_RE0011C_T1_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0011C_T1_C1_D1_segments.txt", header=T)

MEMORI_RE0012A_T5_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0012A_T5_B1_D1_segments.txt", header=T)
MEMORI_RE0012B_T1_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0012B_T1_B1_D1_segments.txt", header=T)
MEMORI_RE0012C_T1_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0012C_T1_B1_D1_segments.txt", header=T)

MEMORI_RE0014A_T1_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0014A_T1_B1_D1_segments.txt", header=T)

MEMORI_RE0021A_T2_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0021A_T2_C1_D1_segments.txt", header=T)
MEMORI_RE0021B_T1_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0021B_T1_C1_D1_segments.txt", header=T)
MEMORI_RE0021C_T1_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0021C_T1_C1_D1_segments.txt", header=T)

MEMORI_RE0022A_T1_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0022A_T1_C1_D1_segments.txt", header=T)
MEMORI_RE0022B_T1_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0022B_T1_B1_D1_segments.txt", header=T)
MEMORI_RE0022C_T1_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0022C_T1_C1_D1_segments.txt", header=T)

MEMORI_RE0023A_T1_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0023A_T1_B1_D1_segments.txt", header=T)
MEMORI_RE0023B_T2_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0023B_T2_B1_D1_segments.txt", header=T)
MEMORI_RE0023C_T1_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0023C_T1_C1_D1_segments.txt", header=T)

MEMORI_RE0024A_T2_B1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0024A_T2_B1_D1_segments.txt", header=T)

MEMORI_RE0025A_T1_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0025A_T1_C1_D1_segments.txt", header=T)
MEMORI_RE0025C_T1_C1_D1 <- read.table("/Users/melissaschmidt/Documents/PostDoc/Lab_Graham/COLL_cluster_files/BCI-EvoCa-OAC/Melissa/WES_Batch1_2/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0025C_T1_C1_D1_segments.txt", header=T)


samnames <- names(seglist)
sortsam <- samnames

# Get unique patients
uniquepat <- unique(gsub('([RN]E\\d+)\\w_\\S+','\\1',samnames))
sortsam <- samnames
patients <- as.factor(gsub('([RN]E\\d+)\\w_\\S+','\\1',sortsam)); patCol <- patients;levels(patCol) <- 1:length(levels(patCol));patCol <- as.numeric(patCol)

regs <- sapply(samnames, function(x) substr(x,14,14))
Responsiveness <- sapply(sortsam, function(x) substr(x, 8,8))
Responsiveness <- ifelse(Responsiveness == "R", 'Responder','Non-Responder')

# Plot segments
pdf(paste0(fig_dir, "/segments_best_fit_ploidy_selected_samples.pdf"), height = 12, width = 12)

# topy defines the highest point on plot - adjust according to number of samples/patients
topy <- 123.4

# Setup blank plot with appropriate dimensions and labels
plot(x=c(0,3.2e9),y=c(0,topy+0.5), col='white',xaxt='n',yaxt='n',xlab='Chromosomes',ylab='',bty='n',yaxs='i',xaxs='i',main='',cex.lab=1,font.lab=2)

# Vectors to track positions and gaps between patients
ypos <- topy+1;patpos <- c();gappos <- c()

# Vectors to track positions and gaps between patients
patpos <- c();gappos <- c()

# For each sample
for(i in c(1:length(sortsam))) {
  sample <- sortsam[i]
  segmuts <- seglist[[sample]]
  
  # If not the first sample
  if(i!=1) {
    # If patient number is different to previous patient number
    if(patCol[i]!=patCol[i-1]) {
      # Input a gap into the plot to separate out patients
      ypos <- ypos-2;gappos <- c(gappos,ypos)
      patpos <- c(patpos,ypos-.5-length(which(gsub('([RN]E\\d+)\\w_\\S+','\\1',sortsam)==patients[i]))/2)
    }
  } else {
    # Input gap for the first sample
    patpos <- c(patpos,ypos-1-length(which(gsub('([RN]E\\d+)\\w_\\S+','\\1',sortsam)==patients[i]))/2)
  }
  
  # Minus 1 from ypos to move down the plot area in order to plot the current sample
  ypos <- ypos-1
  
  # Convert start and end positions to genomic coordinates
  chrIDs <- segmuts$chromosome
  gPositions <- matTab[['plotStart']][match(chrIDs,matTab$chromosome)]
  startPos <- segmuts$start.pos+gPositions
  endPos <- segmuts$end.pos+gPositions
  

  # Select colouts to plot gains and losses
  colour = rep(NA,length=nrow(segmuts))
  colour[which(segmuts$CNA=='None')] = "grey"
  colour[which(segmuts$CNA=='HighGain')] = "#B2182B"
  colour[which(segmuts$CNA=='Gain')] = "#D6604D"
  colour[which(segmuts$CNA=='Mono')] = "#92C5DE"
  colour[which(segmuts$CNA=='cnLOH')] = "#4393C3"
  colour[which(segmuts$CNA=='Del')] = "#2166AC"
  
  
  # Plot segments
  rect(xleft=startPos, xright=endPos, ybottom=(ypos-.4), ytop=(ypos+.4),col=colour,border=NA)
  
  # Plot annotations of region and type to right of plot area
  rect(xleft=2.9e9,xright=2.93e9,ybottom=(ypos-.5), ytop=(ypos+.5),col=TPcols[regs[i]],border=NA)
  rect(xleft=2.93e9,xright=2.96e9,ybottom=(ypos-.5), ytop=(ypos+.5),col=ifelse(Responsiveness[i]=='Responder',"#26A69A", "#b9adde"),border=NA)
}

# Add chromosome and patient labels
patient_label <- paste0((sapply(unique(patients), function(x) substr(x,8,9))), (sapply(unique(patients), function(x) substr(x,12,13))))
patient_label <- str_replace(patient_label, "NE", "NRP")
patient_label <- str_replace(patient_label, "RE", "REP")

axis(side=1,at=c(chrInfo$plotCentroStart),labels=gsub('chr(\\S+)','\\1',chrInfo$chromosome),cex.axis=1,las=1,font=2,line=0)
axis(side=2,at=patpos,labels=patient_label,cex.axis=.8,las=1,font=2,line=0)

# Add extra lines for chromosome starts and ends and centromere starts and ends
abline(v=chrInfo$plotCentroStart,lty=2,lwd=0.3);abline(v=chrInfo$plotCentroEnd,lty=2,lwd=0.3)
abline(v=chrInfo$plotStart,lty=1,lwd=0.6);abline(v=chrInfo$plotEnd,lty=1,lwd=0.6)


# Add more lines to frame the plot area
rect(xleft=0,xright=max(chrInfo$plotEnd),ytop=topy+0.52,ybottom=topy+0.48,lwd=0.6)
rect(xleft=0,xright=max(chrInfo$plotEnd),ytop=0.02,ybottom=-0.02,lwd=0.6)
rect(xleft=1e7,xright=max(chrInfo$plotEnd)-1e7,ytop=gappos+1.5,ybottom=gappos-.5,col='white',border=NA)
rect(xleft=0,xright=max(chrInfo$plotEnd),ytop=gappos+2.6,ybottom=gappos-.6,lwd=0.3)


# Add legends for CNAs and annotations
par(xpd=T,font=2)
legend(x=2.98e9,y=topy,legend=c('CNt>3','CNt=3','1:0','2:0','>2:0'),fill=c("#B2182B","#D6604D","#92C5DE","#4393C3","#2166AC"),border=NA,cex=0.8,box.lwd=0.5,title='CNS')
legend(x=2.98e9,y=topy-20,legend=c('A','B','C'),fill=TPcols,border=NA,cex=0.8,box.lwd=0.5,title='Timepoint')
legend(x=2.98e9,y=topy-35,legend=c('NRP','REP'),fill=c("#b9adde", "#26A69A"),border=NA,cex=0.8,box.lwd=0.5,title='Responsiveness')
par(xpd=F,font=1)
dev.off()








