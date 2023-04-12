# Script to plot CNA segments of multiregion tumour as a whole genome heatmap

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
chrInfo <- read.csv("/path/to/file/centromerePositions_hg38.csv", stringsAsFactors = F);chrInfo <- chrInfo[c(1:22),]
matTab <- data.frame(chrInfo[c("chromosome","plotStart")])

# Include either samnames that have been apporoved for CNA downstream analysys

selected_samples <- read.table("/WES/data_sets/WES_samples_include_for_analysis.txt", sep = "\t", header = T)
selected_samples_Seq <- selected_samples[which(selected_samples$PASS_Sequenza_best_fitting_ploidy == "y"), ]
samnames <- selected_samples_Seq$SampleID



# load example CNA data, you should read your data instead
MEMORI_NE0005A_T2_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0005A_T2_B1_D1_segments.txt", header=T)
MEMORI_NE0005B_T2_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0005B_T2_B1_D1_segments.txt", header=T)
MEMORI_NE0005C_T2_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0005C_T2_C1_D1_segments.txt", header=T)

MEMORI_NE0006A_T2_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0006A_T2_B1_D1_segments.txt", header=T)
MEMORI_NE0006B_T1_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0006B_T1_B1_D1_segments.txt", header=T)
MEMORI_NE0006B_T3_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0006B_T3_C1_D1_segments.txt", header=T)

MEMORI_NE0007A_T1_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0007A_T1_B1_D1_segments.txt", header=T)
MEMORI_NE0007B_T1_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0007B_T1_B1_D1_segments.txt", header=T)

MEMORI_NE0008A_T4_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0008A_T4_B1_D1_segments.txt", header=T)
MEMORI_NE0008B_T2_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0008B_T2_B1_D1_segments.txt", header=T)

MEMORI_NE0009A_T1_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0009A_T1_B1_D1_segments.txt", header=T)
MEMORI_NE0009B_T1_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0009B_T1_B1_D1_segments.txt", header=T)
MEMORI_NE0009C_T1_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0009C_T1_C1_D1_segments.txt", header=T)

MEMORI_NE0010A_T1_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0010A_T1_C1_D1_segments.txt", header=T)
MEMORI_NE0010B_T2_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0010B_T2_C1_D1_segments.txt", header=T)
MEMORI_NE0010B_T3_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0010B_T3_C1_D1_segments.txt", header=T)
MEMORI_NE0010C_T2_M1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0010C_T2_M1_D1_segments.txt", header=T)

MEMORI_NE0012A_T1_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0012A_T1_B1_D1_segments.txt", header=T)
MEMORI_NE0012B_T1_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0012B_T1_C1_D1_segments.txt", header=T)
MEMORI_NE0012C_T1_M1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0012C_T1_M1_D1_segments.txt", header=T)

MEMORI_NE0020A_T1_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0020A_T1_B1_D1_segments.txt", header=T)

MEMORI_NE0021A_T2_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0021A_T2_C1_D1_segments.txt", header=T)
MEMORI_NE0021B_T1_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0021B_T1_B1_D1_segments.txt", header=T)
MEMORI_NE0021C_T1_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0021C_T1_C1_D1_segments.txt", header=T)

MEMORI_NE0022A_T2_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0022A_T2_B1_D1_segments.txt", header=T)
MEMORI_NE0022A_T2_B1_D2 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0022A_T2_B1_D2_segments.txt", header=T)
MEMORI_NE0022B_T1_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0022B_T1_C1_D1_segments.txt", header=T)
MEMORI_NE0022C_T1_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_NE0022C_T1_C1_D1_segments.txt", header=T)


MEMORI_RE0002A_T4_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0002A_T4_C1_D1_segments.txt", header=T)
MEMORI_RE0002B_T1_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0002B_T1_B1_D1_segments.txt", header=T)

MEMORI_RE0003A_T1_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0003A_T1_C1_D1_segments.txt", header=T)
MEMORI_RE0003B_T2_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0003B_T2_C1_D1_segments.txt", header=T)
MEMORI_RE0003C_T1_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0003C_T1_C1_D1_segments.txt", header=T)

MEMORI_RE0004A_T2_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0004A_T2_B1_D1_segments.txt", header=T)
MEMORI_RE0004B_T3_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0004B_T3_B1_D1_segments.txt", header=T)
MEMORI_RE0004C_T1_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0004C_T1_C1_D1_segments.txt", header=T)

MEMORI_RE0005A_T1_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0005A_T1_C1_D1_segments.txt", header=T)
MEMORI_RE0005B_T1_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0005B_T1_B1_D1_segments.txt", header=T)
MEMORI_RE0005C_T1_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0005C_T1_C1_D1_segments.txt", header=T)

MEMORI_RE0006A_T3_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0006A_T3_B1_D1_segments.txt", header=T)
MEMORI_RE0006B_T1_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0006B_T1_B1_D1_segments.txt", header=T)
MEMORI_RE0006C_T3_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0006C_T3_C1_D1_segments.txt", header=T)
MEMORI_RE0006C_T2_M1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0006C_T2_M1_D1_segments.txt", header=T)

MEMORI_RE0007A_T3_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0007A_T3_C1_D1_segments.txt", header=T)
MEMORI_RE0007C_T1_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0007C_T1_C1_D1_segments.txt", header=T)

MEMORI_RE0008A_T1_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0008A_T1_B1_D1_segments.txt", header=T)
MEMORI_RE0008B_T3_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0008B_T3_B1_D1_segments.txt", header=T)
MEMORI_RE0008C_T1_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0008C_T1_C1_D1_segments.txt", header=T)
MEMORI_RE0008C_T2_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0008C_T2_B1_D1_segments.txt", header=T)
MEMORI_RE0008C_T3_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0008C_T3_C1_D1_segments.txt", header=T)

MEMORI_RE0009A_T1_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0009A_T1_B1_D1_segments.txt", header=T)

MEMORI_RE0010A_T1_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0010A_T1_C1_D1_segments.txt", header=T)
MEMORI_RE0010A_T2_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0010A_T2_C1_D1_segments.txt", header=T)
MEMORI_RE0010A_T4_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0010A_T4_C1_D1_segments.txt", header=T)
MEMORI_RE0010B_T1_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0010B_T1_B1_D1_segments.txt", header=T)
MEMORI_RE0010B_T2_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0010B_T2_C1_D1_segments.txt", header=T)
MEMORI_RE0010B_T3_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0010B_T3_C1_D1_segments.txt", header=T)
MEMORI_RE0010C_T1_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0010C_T1_C1_D1_segments.txt", header=T)


MEMORI_RE0011A_T1_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0011A_T1_B1_D1_segments.txt", header=T)
MEMORI_RE0011A_T3_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0011A_T3_B1_D1_segments.txt", header=T)
MEMORI_RE0011B_T1_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0011B_T1_C1_D1_segments.txt", header=T)
MEMORI_RE0011B_T2_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0011B_T2_C1_D1_segments.txt", header=T)
MEMORI_RE0011B_T3_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0011B_T3_C1_D1_segments.txt", header=T)
MEMORI_RE0011C_T1_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0011C_T1_C1_D1_segments.txt", header=T)

MEMORI_RE0012A_T5_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0012A_T5_B1_D1_segments.txt", header=T)
MEMORI_RE0012B_T1_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0012B_T1_B1_D1_segments.txt", header=T)
MEMORI_RE0012C_T1_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0012C_T1_B1_D1_segments.txt", header=T)

MEMORI_RE0014A_T1_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0014A_T1_B1_D1_segments.txt", header=T)

MEMORI_RE0021A_T2_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0021A_T2_C1_D1_segments.txt", header=T)
MEMORI_RE0021B_T1_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0021B_T1_C1_D1_segments.txt", header=T)
MEMORI_RE0021C_T1_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0021C_T1_C1_D1_segments.txt", header=T)

MEMORI_RE0022A_T1_C1_D1<- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0022A_T1_C1_D1_segments.txt", header=T)
MEMORI_RE0022B_T1_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0022B_T1_B1_D1_segments.txt", header=T)
MEMORI_RE0022C_T1_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0022C_T1_C1_D1_segments.txt", header=T)

MEMORI_RE0023A_T1_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0023A_T1_B1_D1_segments.txt", header=T)
MEMORI_RE0023B_T2_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0023B_T2_B1_D1_segments.txt", header=T)
MEMORI_RE0023C_T1_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0023C_T1_C1_D1_segments.txt", header=T)

MEMORI_RE0024A_T2_B1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0024A_T2_B1_D1_segments.txt", header=T)

MEMORI_RE0025A_T1_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0025A_T1_C1_D1_segments.txt", header=T)
MEMORI_RE0025C_T1_C1_D1 <- read.table("/WES/data_sets/Sequenza/final_results/best_fitting_ploidies_per_sample/MEMORI_RE0025C_T1_C1_D1_segments.txt", header=T)





# Classify as gains and losses etc.
for(sam in samnames) {
  cursegs <- get(sam)
  segs <- classifyCNAs(cursegs)
  assign(paste0(sam,"Segs"),segs)
}

# Plot segments

pdf("output_pdf", height = 12, width = 12)
par(mar=c(5.1,12.1,3.1,3.1),font=1)
ypos <- length(samnames);sampos <- c()
plot(x=c(0,3.1e9),y=c(0,ypos+3),col='white',xaxt='n',yaxt='n',xlab='Chromosomes',ylab="",bty='n',yaxs='i',xaxs='i',main='',cex.lab=1,font.lab=2)
for(sam in samnames[10]) {
  segments <- get(paste0(sam,"Segs"))
  sampos <- c(sampos,ypos)
  # Convert start and end positions to genomic coordinates
  chrIDs <- segments$chromosome
  gPositions <- matTab[['plotStart']][match(chrIDs,matTab$chromosome)]
  startPos <- segments$start.pos+gPositions
  endPos <- segments$end.pos+gPositions
  
  # Select colouts to plot gains and losses
  colour = rep(NA,length=nrow(segments))
  colour[which(segments$CNA=='None')] = "grey"
  colour[which(segments$CNA=='HighGain')] = "#B2182B"
  colour[which(segments$CNA=='Gain')] = "#D6604D"
  colour[which(segments$CNA=='Mono')] = "#92C5DE"
  colour[which(segments$CNA=='cnLOH')] = "#4393C3"
  colour[which(segments$CNA=='Del')] = "#2166AC"
  

  # Plot segments
  rect(xleft=startPos, xright=endPos, ybottom=(ypos-.4), ytop=(ypos+.4),col=colour,border=NA)
  
  ypos <- ypos-1
}
ypos <- ypos+1
# Add axis annotations
#axis(side=1,at=c(chrInfo$gCentStart),labels=gsub('chr(\\S+)','\\1',chrInfo$chromosome),cex.axis=1,las=1,font=2,line=0)
label_samnames <- sapply(samnames, function(x) substr(x, 8,17))
label_samnames <- str_replace_all(label_samnames, "NE", "NR")
axis(side=2,at=sampos,labels=label_samnames,cex.axis=.8,las=1,font=2,line=0)

# Add chromosome and centromere lines
abline(v=chrInfo$gCentStart,lty=2,lwd=0.3);abline(v=chrInfo$gCentEnd,lty=2,lwd=0.3)
abline(v=chrInfo$plotStart,lty=1,lwd=0.6);abline(v=chrInfo$plotEnd,lty=1,lwd=0.6)
par(font=2);legend('topright',xpd=T,legend=c('CNt>3','CNt=3','1:0','2:0','>2:0'),fill=c("#B2182B","#D6604D","#92C5DE","#4393C3","#2166AC"),border=NA,cex=0.79,box.lwd=0.79);par(xpd=F,font=1)

dev.off()







