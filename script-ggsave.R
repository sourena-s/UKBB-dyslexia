library(data.table)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(bigsnpr)
library(ggpubr)
library(tidyr)
library(ggnewscale)
core_count=32
args = commandArgs(trailingOnly=TRUE)
cond=as.character(args[1])

a<-fread(paste0('/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/',cond,'/lasso/lasso-hapmap.pgs'),sep=",")
b<-fread(paste0('/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/',cond,'/lasso/grid-hapmap.params'),sep=",")

master_df <- read.table(paste(sep="","/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/",cond,"/component-results/master_df.txt"),header=T)
master_df[,3:20] <- sapply(master_df[,3:20],as.numeric)
#find the max tbm rsq row
master_df[which(master_df$lassosum_rsq==max(master_df[master_df$contrast=="tbm",]$lassosum_rsq)),]
lambda_max <- master_df[which(master_df$lassosum_rsq==max(master_df[master_df$contrast=="tbm",]$lassosum_rsq)),]$lambda
delta_max <- master_df[which(master_df$lassosum_rsq==max(master_df[master_df$contrast=="tbm",]$lassosum_rsq)),]$delta
max_param_tbm <- which(b$lambda==lambda_max & b$delta==delta_max)
max_ic_tbm <-  which(master_df$lassosum_rsq==max(master_df[master_df$contrast=="tbm",]$lassosum_rsq))

#find the max dmri rsq row
master_df[which(master_df$lassosum_rsq==max(master_df[master_df$contrast=="dmri",]$lassosum_rsq)),]
lambda_max <- master_df[which(master_df$lassosum_rsq==max(master_df[master_df$contrast=="dmri",]$lassosum_rsq)),]$lambda
delta_max <- master_df[which(master_df$lassosum_rsq==max(master_df[master_df$contrast=="dmri",]$lassosum_rsq)),]$delta
max_param_dmri <- which(b$lambda==lambda_max & b$delta==delta_max)
max_ic_dmri <-  which(master_df$lassosum_rsq==max(master_df[master_df$contrast=="dmri",]$lassosum_rsq))
write_val <- as.data.frame(t(c(max_param_tbm,max_param_dmri, max_ic_tbm, max_ic_dmri, cond)))
names(write_val) <- c("lassosum_param_col_tbm","lassosum_param_col_dmri","max_ic_tbm","max_ic_dmri","condition")
write.table(write_val, paste(sep="",'/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/',cond,'/lasso/optimisation_points.txt'),sep=" ",col.names=T,quote=F,row.names=F)

n_var_tbm <- fread(paste0('/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/',cond,'/lasso/beta_lassosum2-hapmap.txt'),sep=',')
n_var_tbm <- as.data.frame(n_var_tbm)[,max_param_tbm]
n_var_tbm <- length(n_var_tbm[n_var_tbm!=0])
n_var_tbm <- format(n_var_tbm, big.mark = ",")

n_var_dmri <- fread(paste0('/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/',cond,'/lasso/beta_lassosum2-hapmap.txt'),sep=',')
n_var_dmri <- as.data.frame(n_var_dmri)[,max_param_dmri]
n_var_dmri <- length(n_var_dmri[n_var_dmri!=0])
n_var_dmri <- format(n_var_dmri, big.mark = ",")


#contrast="dmri"; xca="ica";xca_dim=200;comp_id=30
#comp_ids <- which(master_df$contrast==contrast & master_df$XCA==xca & master_df$dim==xca_dim & master_df$ic_num==comp_id)

cond_capital <- paste(toupper(substr(cond, 1, 1)), substr(cond, 2, nchar(cond)), sep="")
comp_ids_tbm <- which(master_df$contrast=="tbm" & master_df$XCA==master_df[max_ic_tbm,]$XCA & master_df$dim==master_df[max_ic_tbm,]$dim & master_df$ic_num==master_df[max_ic_tbm,]$ic_num)
comp_ids_dmri <- which(master_df$contrast=="dmri" & master_df$XCA==master_df[max_ic_dmri,]$XCA & master_df$dim==master_df[max_ic_dmri,]$dim & master_df$ic_num==master_df[max_ic_dmri,]$ic_num)

cs_lines_df_tbm <- gather(master_df[comp_ids_tbm,][1,][,c("cs1_rsq","cs2_rsq","cs3_rsq","cs4_rsq","cs5_rsq","sbayesR_rsq")],key="phi",value="rsq", cs1_rsq, cs2_rsq, cs3_rsq, cs4_rsq, cs5_rsq, sbayesR_rsq)
cs_lines_df_dmri <- gather(master_df[comp_ids_dmri,][1,][,c("cs1_rsq","cs2_rsq","cs3_rsq","cs4_rsq","cs5_rsq","sbayesR_rsq")],key="phi",value="rsq", cs1_rsq, cs2_rsq, cs3_rsq, cs4_rsq, cs5_rsq, sbayesR_rsq)

max_y_bayesian <- max(cs_lines_df_tbm$rsq, cs_lines_df_dmri$rsq ,na.rm=T)


max_y <- max(master_df[comp_ids_tbm,]$lassosum_rsq, master_df[comp_ids_dmri,]$lassosum_rsq, max_y_bayesian ,na.rm=T )

comp_ids <- comp_ids_tbm

if (master_df[max_ic_tbm,]$XCA == "ica") {dim_name=paste(sep="","-",master_df[max_ic_tbm,]$dim)}else{dim_name=""}


cs_lines_df <- gather(master_df[comp_ids,][1,][,c("cs1_rsq","cs2_rsq","cs3_rsq","cs4_rsq","cs5_rsq","sbayesR_rsq")],key="phi",value="rsq", cs1_rsq, cs2_rsq, cs3_rsq, cs4_rsq, cs5_rsq, sbayesR_rsq)

gg_tbm <-  ggplot(master_df[comp_ids,]) + geom_point(size=1, aes(x=log(lambda),color=factor(round(delta,2)),y=lassosum_rsq)) + 
	geom_line(size= 1/ .pt, aes(x=log(lambda),color=factor(round(delta,2)),y=lassosum_rsq,linetype="Frequentist"   )) +
	guides(colour=guide_legend(title="Lassosum2 L₂ regularisation") ) + xlab("Lassosum2 L₁ regularisation") + ylab("IDP explained variance (r²)") + 
	ylim(c(0,1.1*max_y)) +
	ggtitle(paste(sep="", "sMRI TBM\n", toupper(master_df[max_ic_tbm,]$XCA), dim_name, " component ",master_df[max_ic_tbm,]$ic_num )) +
	   guides(linetype=guide_legend(order=1,keywidth = 3 , override.aes=list(alpha=1 ) ) ) + 
   	   scale_color_discrete(name="Lassosum2\nregularisation (L₂)") +
	   guides(color=guide_legend(order=3,keywidth = 3 , override.aes=list(alpha=1, linewidth=3 / .pt ) ) ) + 
	   new_scale_color()+
	geom_hline(data=cs_lines_df, aes(yintercept=rsq ,linetype="Bayesian" ,color=factor(phi)),linewidth=1/ .pt)+ 
	scale_color_manual("Shrinkage type",values=c("cs1_rsq"="red", "cs2_rsq"=grey(0.7), "cs3_rsq"=grey(0.5),"cs4_rsq"=grey(0.3), "cs5_rsq"=grey(0.1), "sbayesR_rsq"="Blue" ), labels=c("PRS-CS (auto)","PRS-CS (Φ=10⁻⁶)","PRS-CS (Φ=10⁻⁴)","PRS-CS (Φ=0.01)","PRS-CS (Φ=1)", "SBayesR (auto)")) +
	scale_linetype_manual("PGS type",values=c("Frequentist"="solid", "Bayesian"="11"))  + 
	guides(linetype=guide_legend(order=1,keywidth = 3, keyheight = 1  , override.aes=list(alpha=1) ) ) +
	guides(color=guide_legend(order=2,keywidth = 3 , override.aes=list(alpha=1,linetype="11" ) ) ) + 
	geom_label_repel(size=7 / .pt , data=master_df[comp_ids,][master_df[comp_ids,]$lassosum_rsq==max(master_df[comp_ids,]$lassosum_rsq),],aes(x=log(lambda),y=lassosum_rsq, label= paste("r² =", round(max(master_df[comp_ids_tbm,]$lassosum_rsq),4),"\nn =",n_var_tbm,"variants" ) ), nudge_x=0.5, nudge_y=25e-5  ) +
	theme( axis.text = element_text(size=8), axis.title = element_text(size=8), text = element_text(size=7), panel.background = element_blank(), axis.line = element_line(colour = "black", size=0.5/ .pt),legend.position="none")

comp_ids <- comp_ids_dmri

cs_lines_df <- gather(master_df[comp_ids,][1,][,c("cs1_rsq","cs2_rsq","cs3_rsq","cs4_rsq","cs5_rsq","sbayesR_rsq")],key="phi",value="rsq", cs1_rsq, cs2_rsq, cs3_rsq, cs4_rsq, cs5_rsq, sbayesR_rsq)

if (master_df[max_ic_dmri,]$XCA == "ica") {dim_name=paste(sep="","-",master_df[max_ic_dmri,]$dim)}else{dim_name=""}
gg_dmri <- ggplot(master_df[comp_ids,]) + geom_point(size=1, aes(x=log(lambda),color=factor(round(delta,2)),y=lassosum_rsq)) + 
	geom_line(size=1/ .pt , aes(x=log(lambda),color=factor(round(delta,2)),y=lassosum_rsq,linetype="Frequentist"   )) +
	guides(colour=guide_legend(title="Lassosum2 L₂ regularisation") ) + xlab("Lassosum2 L₁ regularisation") + ylab("") + #ylab("IDP explained variance (R²)") + 
	ylim(c(0,1.1*max_y)) +
	ggtitle(paste(sep="", "dMRI AFD\n", toupper(master_df[max_ic_dmri,]$XCA), dim_name, " component ",master_df[max_ic_dmri,]$ic_num )) +
	   guides(linetype=guide_legend(order=1,keywidth = 3 , override.aes=list(alpha=1 ) ) ) + 
   	   scale_color_discrete(name="Lassosum2\nregularisation (L₂)") +
	   guides(color=guide_legend(order=3,keywidth = 3 , override.aes=list(alpha=1 ) ) ) + 
	   new_scale_color()+
	geom_hline(data=cs_lines_df, aes(yintercept=rsq ,linetype="Bayesian" ,color=factor(phi)),linewidth=1/ .pt)+ 
	scale_color_manual("Shrinkage type",values=c("cs1_rsq"="red", "cs2_rsq"=grey(0.7), "cs3_rsq"=grey(0.5),"cs4_rsq"=grey(0.3), "cs5_rsq"=grey(0.1), "sbayesR_rsq"="Blue" ), labels=c("PRS-CS (auto)","PRS-CS (Φ=10⁻⁶)","PRS-CS (Φ=10⁻⁴)","PRS-CS (Φ=0.01)","PRS-CS (Φ=1)", "SBayesR (auto)")) +
	scale_linetype_manual("PGS type",values=c("Frequentist"="solid", "Bayesian"="11"))  + 
	guides(linetype=guide_legend(order=1,keywidth = 3, keyheight = 1  , override.aes=list(alpha=1, linewidth=3 / .pt ) ) ) +
	guides(color=guide_legend(order=2,keywidth = 3 , override.aes=list(alpha=1,linetype="11" ) ) ) + 
	geom_label_repel(size=7/ .pt , data=master_df[comp_ids,][master_df[comp_ids,]$lassosum_rsq==max(master_df[comp_ids,]$lassosum_rsq),],aes(x=log(lambda),y=lassosum_rsq, label= paste("R² =", round(max(master_df[comp_ids_dmri,]$lassosum_rsq),4) ,"\nn =",n_var_dmri,"variants"  ) ), nudge_x=0.5, nudge_y=25e-5  ) +
	theme( axis.text = element_text(size=8), axis.title = element_text(size=8), text = element_text(size=7),   panel.background = element_blank(), axis.line = element_line(colour = "black", size=0.5/ .pt),legend.position="none")



gg_dmri_legend <- ggplot(master_df[comp_ids,]) + geom_point(size=1, aes(x=log(lambda),color=factor(round(delta,2)),y=lassosum_rsq)) + 
	geom_line(aes(x=log(lambda),color=factor(round(delta,2)),y=lassosum_rsq,linetype="Frequentist"   )) +
	guides(colour=guide_legend(title="Lassosum2 L₂ regularisation") ) + xlab("Lassosum2 L₁ regularisation") + ylab("") + #ylab("IDP explained variance (R²)") + 
	ylim(c(0,1.1*max_y)) +
	ggtitle(paste(sep="", "dMRI AFD\n", toupper(master_df[max_ic_dmri,]$XCA), dim_name, " component ",master_df[max_ic_dmri,]$ic_num )) +
	   guides(linetype=guide_legend(order=1,keywidth = 3 , override.aes=list(alpha=1,linewidth=c(1,0.75) ) ) ) + 
   	   scale_color_discrete(name="Lassosum2\nL₂ regularisation") +
	   guides(color=guide_legend(order=3,keywidth = 1 , override.aes=list(alpha=1, linewidth=1,size=2 ) ) ) + 
	   new_scale_color()+
	geom_hline(data=cs_lines_df, aes(yintercept=rsq ,linetype="Bayesian" ,color=factor(phi)),linewidth=3/ .pt)+ 
	scale_color_manual("Shrinkage type",values=c("cs1_rsq"="red", "cs2_rsq"=grey(0.7), "cs3_rsq"=grey(0.5),"cs4_rsq"=grey(0.3), "cs5_rsq"=grey(0.1), "sbayesR_rsq"="Blue" ), labels=c("PRS-CS (auto)","PRS-CS (Φ=10⁻⁶)","PRS-CS (Φ=10⁻⁴)","PRS-CS (Φ=0.01)","PRS-CS (Φ=1)", "SBayesR (auto)")) +
	scale_linetype_manual("PGS type",values=c("Frequentist"="solid", "Bayesian"="11"))  + 
	guides(color=guide_legend(order=2,keywidth = 3 , override.aes=list(alpha=1,linetype="11" ) ) ) + 
	theme( axis.text = element_text(size=8), axis.title = element_text(size=8), text = element_text(size=7), panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.key=element_rect(fill="white"), legend.background = element_rect(fill = "white",color="grey",linewidth=0.5))


gg_legend <- get_legend(gg_dmri_legend, position = NULL)

final_plot <- plot_grid(gg_tbm,gg_dmri, gg_legend, ncol=3, align="hv",  rel_widths = c(2, 2,1))
ggsave(paste(sep="","/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/",cond,"/r2_plot.png" ) ,plot=final_plot, bg="white", width=184,height=199,unit="mm",dpi=1200)
ggsave(paste(sep="","/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/",cond,"/r2_plot.svg" ) ,plot=final_plot, bg="white", width=320,height=150,unit="mm" )

