# May 2022
# Tahel Ronel

# Tracking the proportion of the expanded clones across the timepoints

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

input<-"/path/to/"

# First select chain, then plot for alpha and beta separately

# Alpha
chain<-"alpha"
load(paste0(input,"analysis/TCR/data/AB_exp_alpha_df_220904.RData"))
load(paste0(input,"analysis/TCR/data/AC_exp_alpha_df_220904.RData"))
load(paste0(input,"analysis/TCR/data/BC_exp_alpha_df_220904.RData"))
AB_exp_df<-AB_exp_alpha_df
AC_exp_df<-AC_exp_alpha_df
BC_exp_df<-BC_exp_alpha_df

# Beta
chain<-"beta"
load(paste0(input,"analysis/TCR/data/AB_exp_beta_df_220904.RData"))
load(paste0(input,"analysis/TCR/data/AC_exp_beta_df_220904.RData"))
load(paste0(input,"analysis/TCR/data/BC_exp_beta_df_220904.RData"))
AB_exp_df<-AB_exp_beta_df
AC_exp_df<-AC_exp_beta_df
BC_exp_df<-BC_exp_beta_df

##########################################################################


# plotting: choose which fold change you want to plot:

## A to B expansion ##
# x2
fold<-"x2"
dat1<-AB_exp_df[[1]]
dat1<-data.frame(dat1,stringsAsFactors = FALSE)
dat<-data.frame(rbind(cbind(dat1$TCR,dat1$prop1,dat1$name,"A"),cbind(dat1$TCR,dat1$prop2,dat1$name,"B"),cbind(dat1$TCR,dat1$prop3,dat1$name,"C")),stringsAsFactors = FALSE)
dat[,2]<-as.numeric(as.character(dat[,2]))
unique(sort(dat[,2]))[2]
dat[,2][which(dat[,2]==0)]<-0.0001


# x4
fold<-"x4"
dat1<-AB_exp_df[[2]]
dat1<-data.frame(dat1,stringsAsFactors = FALSE)
dat<-data.frame(rbind(cbind(dat1$TCR,dat1$prop1,dat1$name,"A"),cbind(dat1$TCR,dat1$prop2,dat1$name,"B"),cbind(dat1$TCR,dat1$prop3,dat1$name,"C")),stringsAsFactors = FALSE)
dat[,2]<-as.numeric(as.character(dat[,2]))
unique(sort(dat[,2]))[2]
dat[,2][which(dat[,2]==0)]<-0.0001


# x8
fold<-"x8"
dat1<-AB_exp_df[[3]]
dat1<-data.frame(dat1,stringsAsFactors = FALSE)
dat<-data.frame(rbind(cbind(dat1$TCR,dat1$prop1,dat1$name,"A"),cbind(dat1$TCR,dat1$prop2,dat1$name,"B"),cbind(dat1$TCR,dat1$prop3,dat1$name,"C")),stringsAsFactors = FALSE)
dat[,2]<-as.numeric(as.character(dat[,2]))
unique(sort(dat[,2]))[2]
dat[,2][which(dat[,2]==0)]<-0.0001
##

xdat <- factor(dat$X4,levels=c('A','B','C'))
maxy<-max(dat[,2])

p1 <- ggplot(dat, aes(x=xdat, y=as.numeric(as.character(dat[,2])), group=dat[,1])) + geom_point(size=1.5) + geom_line(aes(group=dat[,1]),size=1,col=as.numeric(as.factor(dat[,3]))) +  theme_bw() + theme(text=element_text(size=26))+ xlab('timepoint') + ylab("proportion") + labs(colour='') + scale_y_log10(breaks = c(0.0001,0.001,0.01,0.1),labels = c("0","0.001","0.01","0.1"),limits=c(0.0001,maxy)) 
p1 <- ggplot(dat, aes(x=xdat, y=as.numeric(as.character(dat[,2])), group=dat[,1], col=dat[,3])) + geom_point(size=1.5) + geom_line(aes(group=dat[,1]),size=1) +  theme_bw() + theme(text=element_text(size=26))+ xlab('timepoint') + ylab("proportion") + labs(colour='') + scale_y_log10(breaks = c(0.0001,0.001,0.01,0.1),labels = c("0","0.001","0.01","0.1"),limits=c(0.0001,maxy)) 
p1 <- p1 + ggtitle(paste0(fold," expansion, A to B, ",chain)) 
p1

#ggsave(filename=paste0(input,"/analysis/TCR/plots/",chain,"_",fold,"_track","_AB_", Sys.Date(),'.pdf', sep=''), p1, dpi=300, device='pdf',height=6.15, width=7)



## A to C expansion ##
# x2
fold<-"x2"
dat1<-AC_exp_df[[1]]
dat1<-data.frame(dat1,stringsAsFactors = FALSE)
dat<-data.frame(rbind(cbind(dat1$TCR,dat1$prop1,dat1$name,"A"),cbind(dat1$TCR,dat1$prop2,dat1$name,"B"),cbind(dat1$TCR,dat1$prop3,dat1$name,"C")),stringsAsFactors = FALSE)
dat[,2]<-as.numeric(as.character(dat[,2]))
unique(sort(dat[,2]))[2]
dat[,2][which(dat[,2]==0)]<-0.0001


# x4
fold<-"x4"
dat1<-AC_exp_df[[2]]
dat1<-data.frame(dat1,stringsAsFactors = FALSE)
dat<-data.frame(rbind(cbind(dat1$TCR,dat1$prop1,dat1$name,"A"),cbind(dat1$TCR,dat1$prop2,dat1$name,"B"),cbind(dat1$TCR,dat1$prop3,dat1$name,"C")),stringsAsFactors = FALSE)
dat[,2]<-as.numeric(as.character(dat[,2]))
unique(sort(dat[,2]))[2]
dat[,2][which(dat[,2]==0)]<-0.0001


# x8
fold<-"x8"
dat1<-AC_exp_df[[3]]
dat1<-data.frame(dat1,stringsAsFactors = FALSE)
dat<-data.frame(rbind(cbind(dat1$TCR,dat1$prop1,dat1$name,"A"),cbind(dat1$TCR,dat1$prop2,dat1$name,"B"),cbind(dat1$TCR,dat1$prop3,dat1$name,"C")),stringsAsFactors = FALSE)
dat[,2]<-as.numeric(as.character(dat[,2]))
unique(sort(dat[,2]))[2]
dat[,2][which(dat[,2]==0)]<-0.0001
##

xdat <- factor(dat$X4,levels=c('A','B','C'))
maxy<-max(dat[,2])

p1 <- ggplot(dat, aes(x=xdat, y=as.numeric(as.character(dat[,2])), group=dat[,1])) + geom_point(size=1.5) + geom_line(aes(group=dat[,1]),size=1,col=as.numeric(as.factor(dat[,3]))) +  theme_bw() + theme(text=element_text(size=26))+ xlab('timepoint') + ylab("proportion") + labs(colour='') + scale_y_log10(breaks = c(0.0001,0.001,0.01,0.1),labels = c("0","0.001","0.01","0.1"),limits=c(0.0001,maxy)) 
p1 <- ggplot(dat, aes(x=xdat, y=as.numeric(as.character(dat[,2])), group=dat[,1], col=dat[,3])) + geom_point(size=1.5) + geom_line(aes(group=dat[,1]),size=1) +  theme_bw() + theme(text=element_text(size=26))+ xlab('timepoint') + ylab("proportion") + labs(colour='') + scale_y_log10(breaks = c(0.0001,0.001,0.01,0.1),labels = c("0","0.001","0.01","0.1"),limits=c(0.0001,maxy)) 
p1 <- p1 + ggtitle(paste0(fold," expansion, A to C, ",chain)) 
p1

#ggsave(filename=paste0(input,"/analysis/TCR/",chain,"_",fold,"_track","_AC_", Sys.Date(),'.pdf', sep=''), p1, dpi=300, device='pdf',height=6.15, width=7)



## B to C expansion ##
# x2
fold<-"x2"
dat1<-BC_exp_df[[1]]
dat1<-data.frame(dat1,stringsAsFactors = FALSE)
dat<-data.frame(rbind(cbind(dat1$TCR,dat1$prop1,dat1$name,"A"),cbind(dat1$TCR,dat1$prop2,dat1$name,"B"),cbind(dat1$TCR,dat1$prop3,dat1$name,"C")),stringsAsFactors = FALSE)
dat[,2]<-as.numeric(as.character(dat[,2]))
unique(sort(dat[,2]))[2]
dat[,2][which(dat[,2]==0)]<-0.0001


# x4
fold<-"x4"
dat1<-BC_exp_df[[2]]
dat1<-data.frame(dat1,stringsAsFactors = FALSE)
dat<-data.frame(rbind(cbind(dat1$TCR,dat1$prop1,dat1$name,"A"),cbind(dat1$TCR,dat1$prop2,dat1$name,"B"),cbind(dat1$TCR,dat1$prop3,dat1$name,"C")),stringsAsFactors = FALSE)
dat[,2]<-as.numeric(as.character(dat[,2]))
unique(sort(dat[,2]))[2]
dat[,2][which(dat[,2]==0)]<-0.0001


# x8
fold<-"x8"
dat1<-BC_exp_df[[3]]
dat1<-data.frame(dat1,stringsAsFactors = FALSE)
dat<-data.frame(rbind(cbind(dat1$TCR,dat1$prop1,dat1$name,"A"),cbind(dat1$TCR,dat1$prop2,dat1$name,"B"),cbind(dat1$TCR,dat1$prop3,dat1$name,"C")),stringsAsFactors = FALSE)
dat[,2]<-as.numeric(as.character(dat[,2]))
unique(sort(dat[,2]))[2]
dat[,2][which(dat[,2]==0)]<-0.0001
##

xdat <- factor(dat$X4,levels=c('A','B','C'))
maxy<-max(dat[,2])

p1 <- ggplot(dat, aes(x=xdat, y=as.numeric(as.character(dat[,2])), group=dat[,1])) + geom_point(size=1.5) + geom_line(aes(group=dat[,1]),size=1,col=as.numeric(as.factor(dat[,3]))) +  theme_bw() + theme(text=element_text(size=26))+ xlab('timepoint') + ylab("proportion") + labs(colour='') + scale_y_log10(breaks = c(0.0001,0.001,0.01,0.1),labels = c("0","0.001","0.01","0.1"),limits=c(0.0001,maxy)) 
p1 <- ggplot(dat, aes(x=xdat, y=as.numeric(as.character(dat[,2])), group=dat[,1], col=dat[,3])) + geom_point(size=1.5) + geom_line(aes(group=dat[,1]),size=1) +  theme_bw() + theme(text=element_text(size=26))+ xlab('timepoint') + ylab("proportion") + labs(colour='') + scale_y_log10(breaks = c(0.0001,0.001,0.01,0.1),labels = c("0","0.001","0.01","0.1"),limits=c(0.0001,maxy)) 
p1 <- p1 + ggtitle(paste0(fold," expansion, B to C, ",chain)) 
p1

#ggsave(filename=paste0(input,"/analysis/TCR/plots/",chain,"_",fold,"_track","_BC_", Sys.Date(),'.pdf', sep=''), p1, dpi=300, device='pdf',height=6.15, width=7)




