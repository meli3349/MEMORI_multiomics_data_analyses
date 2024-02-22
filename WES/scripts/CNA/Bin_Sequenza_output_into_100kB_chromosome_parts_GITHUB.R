

###################################################
# Script to convert CNA segments into 100kb bins ##
###################################################
library(data.table);library(dendextend)


# Setup #
chroms <- paste0('chr',c(1:22))
chrInfo <- read.csv("/analysis/General_files/external_files/centromerePositions_hg38.csv", stringsAsFactors = F);chrInfo <- chrInfo[c(1:22),]

# Make 100kb bins across genome
bins <- as.data.frame(matrix(nrow=0,ncol=3));colnames(bins) <- c('Chromosome','gStart','gEnd')
for(i in c(1:nrow(chrInfo))) {
  curchr <- chrInfo[i,]
  start <- curchr$plotStart
  while((start+1.e5)<=curchr$plotCentroStart) {
    end <- start+1e5
    bins[(nrow(bins)+1),] <- c(curchr$chromosome,start,end)
    start <- end
  }
  start <- curchr$plotCentroEnd
  while((start+1e5)<=curchr$plotEnd) {
    end <- start+1e5
    bins[(nrow(bins)+1),] <- c(curchr$chromosome,start,end)
    start <- end
  }
}

# Create a data table and make start and end numeric
bintable <- as.data.table(bins);bintable$gStart <- as.numeric(bintable$gStart)
bintable$gEnd <- as.numeric(bintable$gEnd);setkey(bintable,gStart,gEnd)
binnedCN <- as.data.frame(bintable)


# Get Sequenza files from best fitting ploidy run 
seq_files <- list.files("/analysis/WES/input_files/Sequenza/best_fitting_ploidies_per_sample/", pattern ="*segments.txt", full.names = T)


# Assign total copy number states to every bin
for(sam in seq_files[1]) {
  segments <- read.table(sam, header = T, stringsAsFactors = F)
  segments <- segments[which(segments$chromosome %in% paste0('chr',c(1:22))),]
  sample_name <- gsub("\\S+(MEMORI_\\S+D\\d)\\S+","\\1", sam)
  gPos <- chrInfo[match(segments$chromosome,chrInfo$chromosome),"plotStart"]
  segments$gStart <- segments$start.pos+gPos;segments$gEnd <- segments$end.pos+gPos
  
  segtable <- as.data.table(segments);setkey(segtable, gStart,gEnd)
  
  copynumbins <- as.data.frame(foverlaps(bintable,segtable,type='any',mult='first'))
  row.names(copynumbins) <- c(1:nrow(copynumbins))
  exportbin <- copynumbins[,c('Chromosome','i.gStart','i.gEnd','CNt','A','B')]
  write.table(exportbin,file=paste0("/analysis/WES/data_sets/Sequenza/best_fitting_ploidy_binned_100kb/", sample_name, ".100kb_chrom_bin.txt"), quote= F, row.names=F, sep="\t")
}
