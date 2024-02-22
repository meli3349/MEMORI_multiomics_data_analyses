# Feb 2022
# Tahel Ronel

# Correlation of expanded TCRs with remission grade
# 1. Individual fold changes between timepoints A and B
# 2. All fold changes - one plot for each timepoint alpha/beta
# 3. All timepoints - one plot for each fold change alpha/beta

library(plyr)
library(dplyr)
library(ggplot2)
library(readr)
library(psych)
library(reshape2)
library(vegan)
library(stringr)
library(R.utils)
library(tidyr)
library(igraph)
library(data.table)
library(stringdist)
library(Biobase)
library(gridExtra)
library(grid)
library(ggpubr)

######################################################

# Change to the path where your My_analysis folder is:
input<-"path/to/"

clin_data<-clin_data<-read.table(paste0(input,"analysis/TCR/input_files/Clinical_data_test_and_first_TCR_batch.txt"),sep="\t",header=TRUE) 

rem_grade<-c()

load(paste0(input,"analysis/TCR/data/merged_alphas.RData"))

# First select chain, then plot for alpha and beta separately (i.e. run through rest of the script for each chain)

###########
# Alpha
chain<-"alpha"
load(paste0(input,"analysis/TCR/data/data/all_len_alpha.RData"))
all_len<-all_len_alpha

# Beta
chain<-"beta"
load(paste0(input,"analysis/TCR/data/data/all_len_beta.RData"))
all_len<-all_len_beta
###########

i<-1
for (i in 1:length(merged_alphas)){
  #which(clin_data$patients==names(merged_alphas)[i])
  rem_grade[i]<-clin_data$simplified_remission_grade[which(clin_data$patients==names(merged_alphas)[i])][1]
}
rem_grade

cor.test(rem_grade,all_len_alpha$lenAB2)
cor.test(rem_grade,all_len_alpha$lenAB4)
cor.test(rem_grade,all_len_alpha$lenAB8)
cor.test(rem_grade,all_len_alpha$lenAC2)
cor.test(rem_grade,all_len_alpha$lenAC4) # sig
cor.test(rem_grade,all_len_alpha$lenAC8)
cor.test(rem_grade,all_len_alpha$lenBC2)
cor.test(rem_grade,all_len_alpha$lenBC4)
cor.test(rem_grade,all_len_alpha$lenBC8) # almost sig

cor.test(rem_grade,all_len_beta$lenAB2)
cor.test(rem_grade,all_len_beta$lenAB4) # almost sig
cor.test(rem_grade,all_len_beta$lenAB8)
cor.test(rem_grade,all_len_beta$lenAC2)
cor.test(rem_grade,all_len_beta$lenAC4) 
cor.test(rem_grade,all_len_beta$lenAC8)
cor.test(rem_grade,all_len_beta$lenBC2)
cor.test(rem_grade,all_len_beta$lenBC4)
cor.test(rem_grade,all_len_beta$lenBC8) # sig

# Expansion from A to B for each fold change against remission grade

dat1<-data.frame(rem_grade,all_len$lenAB2,all_len$lenAB4,all_len$lenAB8,c(rep("NE",4),rep("RE",5)))
colnames(dat1)[5]<-"group"
fold<-c("x2","x4","x8")

for (i in 1:3){
  p1<-ggplot(dat1, aes(x=rem_grade, y=dat1[,i+1], group=group)) + geom_point(size=4,alpha=0.7,aes(color=group,fill=group),pch=21) +    theme_bw() + theme(text=element_text(size=26), legend.position='right')+  labs(colour='black')  + xlab("Remission grade") + ylab("Number of expanded TCRs") +  
    guides(color="none") +theme(legend.text = element_text(size=16))  + scale_x_continuous(breaks=c(1,2,3), labels=c(1,2,3))
  p1<-p1 + ggtitle(paste0(chain,", A to B, ",fold[i])) + theme(plot.title = element_text(size=22)) 
  p1
  ggsave(filename=paste0(input,"analysis/TCR/plots/",chain,"_expAB_rem","_",fold[i],"_", Sys.Date(),'.pdf', sep=''), p1, dpi=300, device='pdf',height=6.15, width=7)
}


# Remission grade against the number of expanded TCRs - all fold changes combined

# Expansion from A to B
dat1<-data.frame(rem_grade,all_len$lenAB2,c(rep("NE",4),rep("RE",5)),"x2")
dat2<-data.frame(rem_grade,all_len$lenAB4,c(rep("NE",4),rep("RE",5)),"x4")
dat3<-data.frame(rem_grade,all_len$lenAB8,c(rep("NE",4),rep("RE",5)),"x8")
colnames(dat1)<-c("rem_grade","number_expanded","group","fold_exp")
colnames(dat2)<-c("rem_grade","number_expanded","group","fold_exp")
colnames(dat3)<-c("rem_grade","number_expanded","group","fold_exp")


rem_all_fold<-data.frame(rbind(dat1,dat2,dat3))

cor.test(rem_all_fold$rem_grade,rem_all_fold$number_expanded)

dat1<-rem_all_fold
p1<-ggplot(dat1, aes(x=rem_grade, y=dat1[,2], group=fold_exp)) + geom_point(position=position_jitter(width=0.25), size=4,alpha=0.7,aes(shape=group,color=fold_exp,fill=fold_exp)) +    theme_bw() + theme(text=element_text(size=26), legend.position='right')+  labs(colour='black')  + xlab("Remission grade") + ylab("Number of expanded TCRs") +  
  guides(color=guide_legend("fold_exp"), shape = "legend",fill=guide_legend("fold_exp")) +theme(legend.text = element_text(size=16))  + scale_x_continuous(breaks=c(1,2,3), labels=c(1,2,3))
p1<-p1 + ggtitle(paste0(chain,", A to B")) + theme(plot.title = element_text(size=22)) 
p1
#ggsave(filename=paste0(input,"analysis/TCR/plots/",chain,"_expAB_allfold_rem","_", Sys.Date(),'.pdf', sep=''), p1, dpi=300, device='pdf',height=6.15, width=7)


# Expansion from A to C
dat1<-data.frame(rem_grade,all_len$lenAC2,c(rep("NE",4),rep("RE",5)),"x2")
dat2<-data.frame(rem_grade,all_len$lenAC4,c(rep("NE",4),rep("RE",5)),"x4")
dat3<-data.frame(rem_grade,all_len$lenAC8,c(rep("NE",4),rep("RE",5)),"x8")
colnames(dat1)<-c("rem_grade","number_expanded","group","fold_exp")
colnames(dat2)<-c("rem_grade","number_expanded","group","fold_exp")
colnames(dat3)<-c("rem_grade","number_expanded","group","fold_exp")

rem_all_fold<-data.frame(rbind(dat1,dat2,dat3))

cor.test(rem_all_fold$rem_grade,rem_all_fold$number_expanded)

dat1<-rem_all_fold
p1<-ggplot(dat1, aes(x=rem_grade, y=dat1[,2], group=fold_exp)) + geom_point(position=position_jitter(width=0.25), size=4,alpha=0.7,aes(shape=group,color=fold_exp,fill=fold_exp)) +    theme_bw() + theme(text=element_text(size=26), legend.position='right')+  labs(colour='black')  + xlab("Remission grade") + ylab("Number of expanded TCRs") +  
  guides(color=guide_legend("fold_exp"), shape = "legend",fill=guide_legend("fold_exp")) +theme(legend.text = element_text(size=16))  + scale_x_continuous(breaks=c(1,2,3), labels=c(1,2,3))
p1<-p1 + ggtitle(paste0(chain,", A to C")) + theme(plot.title = element_text(size=22)) #+ coord_cartesian(xlim=c(-1,8.5), ylim=c(-1,8.5)) #+ expand_limits(x=c(-1,8),y=c(-1,8))
p1
#ggsave(filename=paste0(input,"analysis/TCR/plots/",chain,"_expAC_allfold_rem","_", Sys.Date(),'.pdf', sep=''), p1, dpi=300, device='pdf',height=6.15, width=7)


# Expansion from B to C
dat1<-data.frame(rem_grade,all_len$lenBC2,c(rep("NE",4),rep("RE",5)),"x2")
dat2<-data.frame(rem_grade,all_len$lenBC4,c(rep("NE",4),rep("RE",5)),"x4")
dat3<-data.frame(rem_grade,all_len$lenBC8,c(rep("NE",4),rep("RE",5)),"x8")
colnames(dat1)<-c("rem_grade","number_expanded","group","fold_exp")
colnames(dat2)<-c("rem_grade","number_expanded","group","fold_exp")
colnames(dat3)<-c("rem_grade","number_expanded","group","fold_exp")

rem_all_fold<-data.frame(rbind(dat1,dat2,dat3))

cor.test(rem_all_fold$rem_grade,rem_all_fold$number_expanded)

dat1<-rem_all_fold
p1<-ggplot(dat1, aes(x=rem_grade, y=dat1[,2], group=fold_exp)) + geom_point(position=position_jitter(width=0.25), size=4,alpha=0.7,aes(shape=group,color=fold_exp,fill=fold_exp)) +    theme_bw() + theme(text=element_text(size=26), legend.position='right')+  labs(colour='black')  + xlab("Remission grade") + ylab("Number of expanded TCRs") +  
  guides(color=guide_legend("fold_exp"), shape = "legend",fill=guide_legend("fold_exp")) +theme(legend.text = element_text(size=16))  + scale_x_continuous(breaks=c(1,2,3), labels=c(1,2,3))
p1<-p1 + ggtitle(paste0(chain,", B to C")) + theme(plot.title = element_text(size=22)) #+ coord_cartesian(xlim=c(-1,8.5), ylim=c(-1,8.5)) #+ expand_limits(x=c(-1,8),y=c(-1,8))
p1
#ggsave(filename=paste0(input,"analysis/TCR/plots/",chain,"_expBC_allfold_rem","_", Sys.Date(),'.pdf', sep=''), p1, dpi=300, device='pdf',height=6.15, width=7)



# Remission grade against the number of expanded TCRs - all timepoints combined (separated by fold change)

# Fold change x2
dat1<-data.frame(rem_grade,all_len$lenAB2,c(rep("NE",4),rep("RE",5)),"A-B")
dat2<-data.frame(rem_grade,all_len$lenAC2,c(rep("NE",4),rep("RE",5)),"A-C")
dat3<-data.frame(rem_grade,all_len$lenBC2,c(rep("NE",4),rep("RE",5)),"B-C")
colnames(dat1)<-c("rem_grade","number_expanded","group","timepoint")
colnames(dat2)<-c("rem_grade","number_expanded","group","timepoint")
colnames(dat3)<-c("rem_grade","number_expanded","group","timepoint")

rem_all_tp<-data.frame(rbind(dat1,dat2,dat3))

cor.test(rem_all_tp$rem_grade,rem_all_tp$number_expanded) # low

dat1<-rem_all_tp
# significance (Mann-Whitney)
p1<-ggplot(dat1, aes(x=rem_grade, y=dat1[,2], group=timepoint)) + geom_point(position=position_jitter(width=0.25), size=4,alpha=0.7,aes(shape=group,color=timepoint,fill=timepoint)) +    theme_bw() + theme(text=element_text(size=26), legend.position='right')+  labs(colour='black')  + xlab("Remission grade") + ylab("Number of expanded TCRs") +  
  guides(color=guide_legend("timepoint"), shape = "legend",fill=guide_legend("timepoint")) +theme(legend.text = element_text(size=16))  + scale_x_continuous(breaks=c(1,2,3), labels=c(1,2,3))
p1<-p1 + geom_signif(comparisons = list(c(1,3)),map_signif_level=function(p)sprintf("p = %.2f", p), textsize = 6) 
p1<-p1 + ggtitle(paste0(chain,", x2")) + theme(plot.title = element_text(size=22)) 
p1
#ggsave(filename=paste0(input,"analysis/TCR/plots/",chain,"_exp2_all_tp_rem","_", Sys.Date(),'.pdf', sep=''), p1, dpi=300, device='pdf',height=6.15, width=7)

# correlation
p2<-ggplot(dat1, aes(x=dat1[,1], y=dat1[,2])) + theme_bw() + theme(text=element_text(size=22), legend.position='right')+  labs(colour='black')  + xlab("Remission grade") + ylab("Number of expanded TCRs") +  
  geom_smooth(method=lm, se=TRUE, color='#0072B2',fill='#999999') + 
  stat_cor(method="spearman",size=8) +
  scale_x_continuous(breaks=c(1,2,3), labels=c(1,2,3)) + 
  geom_point(position=position_jitter(width=0.25), size=4,alpha=0.7,aes(shape=group,color=timepoint,fill=timepoint)) + 
  theme_bw() + theme(text=element_text(size=22), axis.text =element_text(size=22))+ xlab('Remission grade') + ylab('Number of expanded TCRs') #+ 
p2<-p2 + ggtitle(paste0(chain,", x2")) + theme(plot.title = element_text(size=22)) 
p2<-p2 + guides(color=guide_legend("timepoint"), shape = "legend",fill=guide_legend("timepoint")) +theme(legend.text = element_text(size=16))  
p2
#ggsave(filename=paste0(input,"analysis/TCR/plots/",chain,"_exp2_all_tp_cor_rem","_", Sys.Date(),'.pdf', sep=''), p2, dpi=300, device='pdf',height=6.15, width=7)


# Fold change x4
dat1<-data.frame(rem_grade,all_len$lenAB4,c(rep("NE",4),rep("RE",5)),"A-B")
dat2<-data.frame(rem_grade,all_len$lenAC4,c(rep("NE",4),rep("RE",5)),"A-C")
dat3<-data.frame(rem_grade,all_len$lenBC4,c(rep("NE",4),rep("RE",5)),"B-C")
colnames(dat1)<-c("rem_grade","number_expanded","group","timepoint")
colnames(dat2)<-c("rem_grade","number_expanded","group","timepoint")
colnames(dat3)<-c("rem_grade","number_expanded","group","timepoint")

rem_all_tp<-data.frame(rbind(dat1,dat2,dat3))

cor.test(rem_all_tp$rem_grade,rem_all_tp$number_expanded) # sig

dat1<-rem_all_tp
# significance (Mann-Whitney)
p1<-ggplot(dat1, aes(x=rem_grade, y=dat1[,2], group=timepoint)) + geom_point(position=position_jitter(width=0.25), size=4,alpha=0.7,aes(shape=group,color=timepoint,fill=timepoint)) +    theme_bw() + theme(text=element_text(size=26), legend.position='right')+  labs(colour='black')  + xlab("Remission grade") + ylab("Number of expanded TCRs") +  
  guides(color=guide_legend("timepoint"), shape = "legend",fill=guide_legend("timepoint")) +theme(legend.text = element_text(size=16))  + scale_x_continuous(breaks=c(1,2,3), labels=c(1,2,3))
p1<-p1 + geom_signif(comparisons = list(c(1,3)),map_signif_level=function(p)sprintf("p = %.2f", p), textsize = 6) 
p1<-p1 + ggtitle(paste0(chain,", x4")) + theme(plot.title = element_text(size=22)) 
p1
#ggsave(filename=paste0(input,"analysis/TCR/plots/",chain,"_exp4_all_tp_rem","_", Sys.Date(),'.pdf', sep=''), p1, dpi=300, device='pdf',height=6.15, width=7)

# correlation
p2<-ggplot(dat1, aes(x=dat1[,1], y=dat1[,2])) + theme_bw() + theme(text=element_text(size=22), legend.position='right')+  labs(colour='black')  + xlab("Remission grade") + ylab("Number of expanded TCRs") +  
  geom_smooth(method=lm, se=TRUE, color='#0072B2',fill='#999999') + 
  stat_cor(method="spearman",size=8) +
  scale_x_continuous(breaks=c(1,2,3), labels=c(1,2,3)) + 
  geom_point(position=position_jitter(width=0.25), size=4,alpha=0.7,aes(shape=group,color=timepoint,fill=timepoint)) + 
  theme_bw() + theme(text=element_text(size=22), axis.text =element_text(size=22))+ xlab('Remission grade') + ylab('Number of expanded TCRs') #+ 
p2<-p2 + ggtitle(paste0(chain,", x4")) + theme(plot.title = element_text(size=22)) 
p2<-p2 + guides(color=guide_legend("timepoint"), shape = "legend",fill=guide_legend("timepoint")) +theme(legend.text = element_text(size=16))  
p2
#ggsave(filename=paste0(input,"analysis/TCR/plots/",chain,"_exp4_all_tp_cor_rem","_", Sys.Date(),'.pdf', sep=''), p2, dpi=300, device='pdf',height=6.15, width=7)


# Fold change x8
dat1<-data.frame(rem_grade,all_len$lenAB8,c(rep("NE",4),rep("RE",5)),"A-B")
dat2<-data.frame(rem_grade,all_len$lenAC8,c(rep("NE",4),rep("RE",5)),"A-C")
dat3<-data.frame(rem_grade,all_len$lenBC8,c(rep("NE",4),rep("RE",5)),"B-C")
colnames(dat1)<-c("rem_grade","number_expanded","group","timepoint")
colnames(dat2)<-c("rem_grade","number_expanded","group","timepoint")
colnames(dat3)<-c("rem_grade","number_expanded","group","timepoint")

rem_all_tp<-data.frame(rbind(dat1,dat2,dat3))

cor.test(rem_all_tp$rem_grade,rem_all_tp$number_expanded) # sig

dat1<-rem_all_tp
# significance (Mann-Whitney)
p1<-ggplot(dat1, aes(x=rem_grade, y=dat1[,2], group=timepoint)) + geom_point(position=position_jitter(width=0.25), size=4,alpha=0.7,aes(shape=group,color=timepoint,fill=timepoint)) +    theme_bw() + theme(text=element_text(size=26), legend.position='right')+  labs(colour='black')  + xlab("Remission grade") + ylab("Number of expanded TCRs") +  
  guides(color=guide_legend("timepoint"), shape = "legend",fill=guide_legend("timepoint")) +theme(legend.text = element_text(size=16))  + scale_x_continuous(breaks=c(1,2,3), labels=c(1,2,3))
p1<-p1 + geom_signif(comparisons = list(c(1,3)),map_signif_level=function(p)sprintf("p = %.2f", p), textsize = 6) 
p1<-p1 + ggtitle(paste0(chain,", x8")) + theme(plot.title = element_text(size=22)) 
p1
#ggsave(filename=paste0(input,"analysis/TCR/plots/",chain,"_exp8_all_tp_rem","_", Sys.Date(),'.pdf', sep=''), p1, dpi=300, device='pdf',height=6.15, width=7)

# correlation
p2<-ggplot(dat1, aes(x=dat1[,1], y=dat1[,2])) + theme_bw() + theme(text=element_text(size=22), legend.position='right')+  labs(colour='black')  + xlab("Remission grade") + ylab("Number of expanded TCRs") +  
  geom_smooth(method=lm, se=TRUE, color='#0072B2',fill='#999999') + 
  stat_cor(method="spearman",size=8) +
  scale_x_continuous(breaks=c(1,2,3), labels=c(1,2,3)) + 
  geom_point(position=position_jitter(width=0.25), size=4,alpha=0.7,aes(shape=group,color=timepoint,fill=timepoint)) + 
  theme_bw() + theme(text=element_text(size=22), axis.text =element_text(size=22))+ xlab('Remission grade') + ylab('Number of expanded TCRs') #+ 
p2<-p2 + ggtitle(paste0(chain,", x8")) + theme(plot.title = element_text(size=22)) 
p2<-p2 + guides(color=guide_legend("timepoint"), shape = "legend",fill=guide_legend("timepoint")) +theme(legend.text = element_text(size=16))  
p2
#ggsave(filename=paste0(input,"analysis/TCR/plots/",chain,"_exp8_all_tp_cor_rem","_", Sys.Date(),'.pdf', sep=''), p2, dpi=300, device='pdf',height=6.15, width=7)

