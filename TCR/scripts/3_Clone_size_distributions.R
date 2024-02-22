# Feb 2022
# Tahel Ronel

# Clone size distribution of samples/subsamples

library(ggplot2)
library(ggpubr)
library(plyr)
library(dplyr)
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

# Change to path to My_analysis
input<-"/path/to"

load(paste0(input,"analysis/TCR/data/merged_alphas.RData"))
load(paste0(input,"analysis/TCR/data/merged_betas.RData"))

load(paste0(input,"analysis/TCR/data/alldcrs_f_alpha.Rdata"))
load(paste0(input,"analysis/TCR/data/alldcrs_f_beta.Rdata"))

#directory to which plots are saved. 
output<-paste0(input,"analysis/TCR/plots/")

# First select chain, then plot for alpha and beta separately (i.e. run through rest of the script for each chain)

###########
# Alpha
chain<-"alpha"
samples<-merged_alphas
samples_ss<-alldcrs_f_alpha

# Beta
chain<-"beta"
samples<-merged_betas
samples_ss<-alldcrs_f_beta

###########


#plotting abundance profiles

# 1. All TCRs

# All TCRs at timepoint A

pats<-names(merged_alphas)
#group<-c(rep("NE",4),rep("RE",5))
Alist<-list()
i<-1
for (i in 1:length(pats)){
  Alist[[i]]<-data.frame(cbind(pats[i],samples[[i]][,2])) #,group[i]))
}
Alist1<-data.frame(rbindlist(Alist))
Alist1[,2]<-as.numeric(as.character(Alist1[,2]))
zeros<-which(Alist1[,2]==0)
Alist1<-Alist1[-zeros,]
colnames(Alist1)<-c("patient","count") #,"group")

# All TCRs at timepoint B

Blist<-list()
i<-1
for (i in 1:length(pats)){
  Blist[[i]]<-data.frame(cbind(pats[i],samples[[i]][,3])) #,group[i]))
}
Blist1<-data.frame(rbindlist(Blist))
Blist1[,2]<-as.numeric(as.character(Blist1[,2]))
zeros<-which(Blist1[,2]==0)
Blist1<-Blist1[-zeros,]
colnames(Blist1)<-c("patient","count") #,"group")

ncols<-c()
for (i in 1:length(pats)){
  ncols[i]<-length(samples[[i]])  
}
pats2<-which(ncols==7)
#names2<-names(merged[pats2])

# All TCRs at timepoint C

Clist<-list()
i<-1
for (i in pats2){
  Clist[[i]]<-data.frame(cbind(pats[i],samples[[i]][,4]))#,group[i]))
}
Clist1<-data.frame(rbindlist(Clist))
Clist1[,2]<-as.numeric(as.character(Clist1[,2]))
zeros<-which(Clist1[,2]==0)
Clist1<-Clist1[-zeros,]
colnames(Clist1)<-c("patient","count") #,"group")



# Plotting clone size histogram

hist1<-hist(log2(Alist1[,2]), breaks=c(0:12),plot=FALSE)$counts
#table(samples[[i]]$duplicate_count)
cum_hist1<-log10(hist1/sum(hist1))

hist2<-hist(log2(Blist1[,2]), breaks=c(0:12),plot=FALSE)$counts
cum_hist2<-log10(hist2/sum(hist2))

hist3<-hist(log2(Clist1[,2]), breaks=c(0:12),plot=FALSE)$counts
cum_hist3<-log10(hist3/sum(hist3))

Adf<-data.frame(1:12,cum_hist1,"A")
colnames(Adf)<-c("Abundance","Proportion","Timepoint")
Bdf<-data.frame(1:12,cum_hist2,"B")  
colnames(Bdf)<-c("Abundance","Proportion","Timepoint")
Cdf<-data.frame(1:12,cum_hist3,"C")  
colnames(Cdf)<-c("Abundance","Proportion","Timepoint")

dat<-data.frame(rbind(Adf,Bdf,Cdf))

ylabel<-"Proportion"
xlabel<-"Abundance (log2)"

xdat <- dat[,1]
gr <- dat[,3]

shapes<-c(1,3,5,7,9)

inf<-which(dat$Proportion==-Inf)
dat<-dat[-inf,]

p1 <- ggplot(dat, aes(x=Abundance, y=Proportion)) + geom_point(aes(color=factor(dat[,3]), shape=factor(dat[,3]), fill=factor(dat[,3])),size=5,alpha=0.6)  +  theme_bw() + theme(text=element_text(size=26), legend.position = c(.95, .95), legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6), legend.title=element_blank(), plot.title=element_text(size=26))+ xlab(xlabel) + ylab(ylabel) + labs(colour='ID',shape='ID',fill='ID') + scale_y_continuous(breaks=c(0,-1,-2,-3,-4,-5),labels = c(1,0.1,0.01,0.001,"0.0001","0.00001")) +  scale_x_continuous(breaks = c(1:13), labels=c(0,1,2,3,4,5,6,7,8,9,10,11,12)) +
  scale_color_manual(values=c("#D55E00","#009E73","grey50"))
# p1<-p1 + ggtitle(paste0(chain,", clone size distribution, combined patients"))
p1

ggsave(filename=paste(output,chain,"_cumu_prop",'_', Sys.Date(),'.pdf', sep=''), p1, dpi=300, device='pdf' ,height=6.15, width=7.5)


# Plotting all clone sizes

tab1<-table(Alist1[,2])
prop1<-tab1/(sum(Alist1[,2]))
df1<-data.frame(cbind(as.numeric(names(tab1)),log2(as.numeric(names(tab1))),prop1,"A"))
colnames(df1)<-c("abundance","abundance (log2)","proportion","timepoint")

tab2<-table(Blist1[,2])
prop2<-tab2/(sum(Blist1[,2]))
df2<-data.frame(cbind(as.numeric(names(tab2)),log2(as.numeric(names(tab2))),prop2,"B"))
colnames(df2)<-c("abundance","abundance (log2)","proportion","timepoint")

tab3<-table(Clist1[,2])
prop3<-tab3/(sum(Clist1[,2]))
df3<-data.frame(cbind(as.numeric(names(tab3)),log2(as.numeric(names(tab3))),prop3,"C"))
colnames(df3)<-c("abundance","abundance (log2)","proportion","timepoint")

dat<-data.frame(rbind(df1,df2,df3))
dat[,1]<-as.numeric(dat[,1])
dat[,2]<-as.numeric(dat[,2])
dat[,3]<-as.numeric(dat[,3])

ylabel<-"Proportion"
xlabel<-"Abundance (log2)"

xdat <- dat[,2]
gr <- dat[,4]

shapes<-c(1,3,5,7,9)

p2 <- ggplot(dat, aes(x=dat[,2], y=dat[,3])) + geom_point(aes(color=factor(dat[,4]), shape=factor(dat[,4]), fill=factor(dat[,4])),size=5,alpha=0.6)  +  theme_bw() + theme(text=element_text(size=26), legend.position = c(.95, .95), legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6), legend.title=element_blank(), plot.title=element_text(size=26))+ xlab(xlabel) + ylab(ylabel) + labs(colour='name',shape='name',fill='name') + scale_y_log10(breaks = c(1,0.1,0.01,0.001,0.0001,0),labels=c(1,0.1,0.01,0.001,0.0001,0)) +    #scale_x_continuous(breaks = c(1:13), labels=c(0,1,2,3,4,5,6,7,8,9,10,11,12)) +
  scale_color_manual(values=c("#D55E00","#009E73","grey50"))
# p2<- p2 + ggtitle(paste0(chain,", clone size distribution - all clone sizes, combined patients"))
p2

ggsave(filename=paste(output,chain,"_prop_all_", Sys.Date(),'.pdf', sep=''), p2, dpi=300, device='pdf' ,height=6.15, width=7.5)



# 2. Timepoints subsampled to the same number of TCRs

# Plotting clone size histogram
Asamples<-samples_ss[which(samples_ss$timepoint=="A"),]
Bsamples<-samples_ss[which(samples_ss$timepoint=="B"),]
Csamples<-samples_ss[which(samples_ss$timepoint=="C"),]

hist1<-hist(log2(Asamples[,1]), breaks=c(0:12),plot=FALSE)$counts
#table(samples[[i]]$duplicate_count)
cum_hist1<-log10(hist1/sum(hist1))

hist2<-hist(log2(Bsamples[,1]), breaks=c(0:12),plot=FALSE)$counts
cum_hist2<-log10(hist2/sum(hist2))

hist3<-hist(log2(Csamples[,1]), breaks=c(0:12),plot=FALSE)$counts
cum_hist3<-log10(hist3/sum(hist3))

Adf<-data.frame(1:12,cum_hist1,"A")
colnames(Adf)<-c("Abundance","Proportion","Timepoint")
Bdf<-data.frame(1:12,cum_hist2,"B")  
colnames(Bdf)<-c("Abundance","Proportion","Timepoint")
Cdf<-data.frame(1:12,cum_hist3,"C")  
colnames(Cdf)<-c("Abundance","Proportion","Timepoint")

dat<-data.frame(rbind(Adf,Bdf,Cdf))

ylabel<-"Proportion"
xlabel<-"Abundance (log2)"

xdat <- dat[,1]
gr <- dat[,3]

shapes<-c(1,3,5,7,9)

inf<-which(dat$Proportion==-Inf)
dat<-dat[-inf,]

p1 <- ggplot(dat, aes(x=Abundance, y=Proportion)) + geom_point(aes(color=factor(dat[,3]), shape=factor(dat[,3]), fill=factor(dat[,3])),size=5,alpha=0.6)  +  theme_bw() + theme(text=element_text(size=26), legend.position = c(.95, .95), legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6), legend.title=element_blank(), plot.title=element_text(size=26))+ xlab(xlabel) + ylab(ylabel) + labs(colour='ID',shape='ID',fill='ID') + scale_y_continuous(breaks=c(0,-1,-2,-3,-4,-5),labels = c(1,0.1,0.01,0.001,"0.0001","0.00001")) +  scale_x_continuous(breaks = c(1:13), labels=c(0,1,2,3,4,5,6,7,8,9,10,11,12)) +
  scale_color_manual(values=c("#D55E00","#009E73","grey50"))
# p1<-p1 + ggtitle(paste0(chain,", clone size distribution, combined patients"))
p1

ggsave(filename=paste(output,chain,"_cumu_prop_ss",'_', Sys.Date(),'.pdf', sep=''), p1, dpi=300, device='pdf' ,height=6.15, width=7.5)


# Plotting all clone sizes

tab1<-table(Asamples[,1])
prop1<-tab1/(sum(Asamples[,1]))
df1<-data.frame(cbind(as.numeric(names(tab1)),log2(as.numeric(names(tab1))),prop1,"A"))
colnames(df1)<-c("abundance","abundance (log2)","proportion","timepoint")

tab2<-table(Bsamples[,1])
prop2<-tab2/(sum(Bsamples[,1]))
df2<-data.frame(cbind(as.numeric(names(tab2)),log2(as.numeric(names(tab2))),prop2,"B"))
colnames(df2)<-c("abundance","abundance (log2)","proportion","timepoint")

tab3<-table(Csamples[,1])
prop3<-tab3/(sum(Csamples[,1]))
df3<-data.frame(cbind(as.numeric(names(tab3)),log2(as.numeric(names(tab3))),prop3,"C"))
colnames(df3)<-c("abundance","abundance (log2)","proportion","timepoint")

dat<-data.frame(rbind(df1,df2,df3))
dat[,1]<-as.numeric(dat[,1])
dat[,2]<-as.numeric(dat[,2])
dat[,3]<-as.numeric(dat[,3])

ylabel<-"Proportion"
xlabel<-"Abundance (log2)"

xdat <- dat[,2]
gr <- dat[,4]

shapes<-c(1,3,5,7,9)

p2 <- ggplot(dat, aes(x=dat[,2], y=dat[,3])) + geom_point(aes(color=factor(dat[,4]), shape=factor(dat[,4]), fill=factor(dat[,4])),size=5,alpha=0.6)  +  theme_bw() + theme(text=element_text(size=26), legend.position = c(.95, .95), legend.justification = c("right", "top"),legend.box.just = "right",legend.margin = margin(6, 6, 6, 6), legend.title=element_blank(), plot.title=element_text(size=26))+ xlab(xlabel) + ylab(ylabel) + labs(colour='name',shape='name',fill='name') + scale_y_log10(breaks = c(1,0.1,0.01,0.001,0.0001,0),labels=c(1,0.1,0.01,0.001,0.0001,0)) +    #scale_x_continuous(breaks = c(1:13), labels=c(0,1,2,3,4,5,6,7,8,9,10,11,12)) +
  scale_color_manual(values=c("#D55E00","#009E73","grey50"))
# p2<- p2 + ggtitle(paste0(chain,", clone size distribution - all clone sizes, combined patients"))
p2

ggsave(filename=paste(output,chain,"_prop_all_ss_", Sys.Date(),'.pdf', sep=''), p2, dpi=300, device='pdf' ,height=6.15, width=7.5)
