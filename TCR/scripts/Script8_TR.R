# December 2022
# Tahel Ronel

# Creating fishplots of the observed expanded TCRs for each patient for which we have three timepoints

options(stringsAsFactors = FALSE)
require(ggpubr)
require(reshape2)
require(gdata)
require(Vennerable)
require(grid)
require(fishplot)

input<-"/path/to/"

output1<-paste0(input, "analysis/TCR/plots/")

load(paste0(input, "analysis/TCR/data/AB_exp_alpha_df_220904.RData"))
load(paste0(input, "analysis/TCR/data/AB_exp_beta_df_220904.RData"))
load(paste0(input, "analysis/TCR/data/AC_exp_alpha_df_220904.RData"))
load(paste0(input, "analysis/TCR/data/AC_exp_beta_df_220904.RData"))
load(paste0(input, "analysis/TCR/data/BC_exp_alpha_df_220904.RData"))
load(paste0(input, "analysis/TCR/data/BC_exp_beta_df_220904.RData"))


# Alpha set
AB_alpha<-AB_exp_alpha_df[[2]]
colnames(AB_alpha)[9]<-"timepoint1"
AC_alpha<-AC_exp_alpha_df[[2]]
colnames(AC_alpha)[9]<-"timepoint2"
BC_alpha<-BC_exp_alpha_df[[2]]
colnames(BC_alpha)[9]<-"timepoint3"

AB_AC_alpha<-merge(AB_alpha,AC_alpha,all=TRUE)
AB_AC_BC_alpha<-merge(AB_AC_alpha,BC_alpha,all=TRUE)
AB_AC_BC_alpha[is.na(AB_AC_BC_alpha)] <- 0

# Beta set
AB_beta<-AB_exp_beta_df[[2]]
colnames(AB_beta)[9]<-"timepoint1"
AC_beta<-AC_exp_beta_df[[2]]
colnames(AC_beta)[9]<-"timepoint2"
BC_beta<-BC_exp_beta_df[[2]]
colnames(BC_beta)[9]<-"timepoint3"

AB_AC_beta<-merge(AB_beta,AC_beta,all=TRUE)
AB_AC_BC_beta<-merge(AB_AC_beta,BC_beta,all=TRUE)
AB_AC_BC_beta[is.na(AB_AC_BC_beta)] <- 0



for (my_chain in c("alpha", "beta")){
  
  if (my_chain=="alpha") exp_set <- AB_AC_BC_alpha #exp_R_all_alpha_CSC 
  if (my_chain=="beta") exp_set <- AB_AC_BC_beta #exp_R_all_beta_CSC

  
  for (i in 1:4) { 
    
    c <- FALSE
    
    pat_id <- unique(exp_set$name)[i]
    print(pat_id)
    
    tcr_tmp <- subset(exp_set, name==pat_id) 
    #tcr_tmp
    
    timepoints <- c(0,150, 300) 

    ##################################################################################################################

    seg_111 <- subset(tcr_tmp, ((A>0)&(B>0)&(C>0)))
    seg_100 <- subset(tcr_tmp, ((A>0)&(B==0)&(C==0)))
    seg_010 <- subset(tcr_tmp, ((A==0)&(B>0)&(C==0)))
    seg_001 <- subset(tcr_tmp, ((A==0)&(B==0)&(C>0)))
    seg_110 <- subset(tcr_tmp, ((A>0)&(B>0)&(C==0)))
    seg_101 <- subset(tcr_tmp, ((A>0)&(B==0)&(C>0)))
    seg_011 <- subset(tcr_tmp, ((A==0)&(B>0)&(C>0)))
    
    dim_A <- nrow(subset(tcr_tmp, A>0))
    dim_B <- nrow(subset(tcr_tmp, B>0))
    dim_C <- nrow(subset(tcr_tmp, C>0))
    

    seg <- c(nrow(seg_111), nrow(seg_100), nrow(seg_010), nrow(seg_001), nrow(seg_110), nrow(seg_101), nrow(seg_011))

    my_filter <- list(
      c(1, 1, 1),
      c(1, 0, 0),
      c(0, 1, 0),
      c(0, 0, 1),
      c(1, 1, 0),
      c(1, 0, 1),
      c(0, 1, 1)  
    )

    frac.table <- data.frame()
    dim_max <- max(dim_A, dim_B, dim_C)
    
    dim_all <- c(dim_A, dim_B, dim_C)

    prop_A <- 100*dim_A/dim_max
    prop_B <- 100*dim_B/dim_max
    prop_C <- 100*dim_C/dim_max

    props <- c(prop_A, prop_B, prop_C)

    frac.table <- rbind(frac.table, unlist(lapply(props, FUN=function(x){min(100, (round(x)+1))})))
    
    colnames(frac.table) <- c("A", "B", "C")
    
    for (j in 1:7) {

      my_vec <- seg[j]/dim_all * my_filter[[j]]
      frac.table <- rbind(frac.table, my_vec * props)

    }

    frac.table <- subset(frac.table, (A>0)|(B>0)|(C>0))
    frac.table[is.na(frac.table)]<-0
    
    parents <- c(0, rep(1, nrow(frac.table)-1))
    frac_mat <- as.matrix(floor(frac.table))
    
    ##################################################################################################################

    fish_b = createFishObject(frac_mat,parents,timepoints=timepoints, fix.missing.clones=TRUE)
    fish_b = layoutClones(fish_b)

    pdf(paste0(output1, 'Fish_', pat_id, "_", my_chain, "_", Sys.Date(),'.pdf'), width=8, height=5)
    fishPlot(fish_b,shape="spline",title =paste0(pat_id,"- Occurrence of Expanded TCRs at each timepoint, ", my_chain),
                 cex.title=1, cex.vlab =1, vlines=c(0,150, 300), 
                 vlab=c("A","B", "C"), bg.type = "solid", bg.col = c("white", "white", "white"))
    dev.off()

  }
}