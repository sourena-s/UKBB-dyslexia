library(data.table)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(bigsnpr)
library(ggpubr)

core_count=16
args = commandArgs(trailingOnly=TRUE)
cond=as.character(args[1])
sbayesr_data <-  read.table(paste(sep="",'/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/',cond, '/auto-sbayesr.pgs.snpRes'),header=T)
names(sbayesr_data) <- c("id", "rsid","chr", "pos","a1","a0","freq","beta","beta_se","PIP","LastSampleEff")
obj.bigSNP.all <- snp_attach("/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/processed/big40_all_subj_subset_hapmap_2k.rds")
map <- as.data.frame(cbind(obj.bigSNP.all$map$chromosome,obj.bigSNP.all$map$rsid,obj.bigSNP.all$map$physical.pos,obj.bigSNP.all$map$allele1, obj.bigSNP.all$map$allele2))
#check ref allele
names(map) <- c("chr","rsid","pos","a0","a1")
map$chr<-as.integer(map$chr)
sbayesr_match <- snp_match(sbayesr_data, join_by_pos = FALSE, map, return_flip_and_rev=TRUE, strand_flip =FALSE)
sbayesr_score <- big_prodMat(obj.bigSNP.all$genotype, as.matrix(sbayesr_match$beta), ind.col = sbayesr_match$'_NUM_ID_',  ncores = core_count, block.size=5000 )
write.table(sbayesr_score,        file=paste(sep="",'/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/',cond, '/auto-sbayesr.pgs'),quote=F,col.names=F,row.names=F)


#sbayesr_score <- read.table(paste(sep="",'/data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/prs-thresholds/',cond,'/100.txt.sorted'),header=F)[,2]
pgs_subid<-read.table('/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/list-all-subjects-genetics-subid-and-indices.txt',header=F)[,1]
#sbayesr_score <- as.numeric(read.table(paste(sep="",'/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/',cond, '/auto-sbayesr.pgs'),header=F)[,1])

#creating dMRI covariate dataframe
base="/data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/"
modality="dmri"
cov_base <- paste(sep="",base,modality,"/")
cov1 <- read.table( paste(sep="",cov_base,"c0-post-QC") ,header=F)
cov2 <- read.table( paste(sep="",cov_base,"c1-sex-age-site-sitedate-rescanDateAgecorrected") ,header=F)
cov3 <- read.table( paste(sep="",cov_base,"c2-genomic-PCs") ,header=F)
cov4 <- read.table( paste(sep="",cov_base,"c3-AgeSq-SexByAge-SexByAgeSq") ,header=F)
cov5 <- read.table( paste(sep="",cov_base,"c4-genoarray") ,header=F)
cov6 <- read.table( paste(sep="",cov_base,"c5-inHouseEddy") ,header=F)
covariates <- cbind (cov1,cov2,cov3,cov4,cov5,cov6)
covariates_dMRI <- covariates[covariates[,2]=="OK",c(1,seq(3,26))]
colnames(covariates_dMRI) <- c("subid", paste0("COV",1:24))
subids_dmri <- covariates_dMRI[,1]
dmri_pgs_indices <- match(as.integer(subids_dmri) ,as.integer(pgs_subid) )

#we extract subids from the pruned subject column (pass==OK) below
base="/data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/"
modality="tbm"
cov_base <- paste(sep="",base,modality,"/")
cov1 <- read.table( paste(sep="",cov_base,"c0-post-QC") ,header=F)
cov2 <- read.table( paste(sep="",cov_base,"c0-5-coded") ,header=F)
cov3 <- read.table( paste(sep="",cov_base,"c7-pca_ten") ,header=F)
cov4 <- read.table( paste(sep="",cov_base,"design-agesq-ageXsex-agesqXsex") ,header=F)
cov5 <- read.table( paste(sep="",cov_base,"base-genoBatch") ,header=F)
covariates <- cbind (cov1,cov2,cov3,cov4,cov5)
covariates_tbm <- covariates[covariates[,2]=="OK",c(1,seq(4,27))]
colnames(covariates_tbm) <- c("subid", paste0("COV",1:24))
subids_tbm <- covariates_tbm[,1]
tbm_pgs_indices <- match(as.integer(subids_tbm) ,as.integer(pgs_subid) )

print(paste("Loading data for condition:",cond))
a<-read.table(paste(sep="",'/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/',cond,'/lasso/lasso-hapmap.pgs'),sep=",",header=F)
b<-read.table(paste(sep="",'/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/',cond,'/lasso/grid-hapmap.params'),sep=",",header=T)
subids_tbm_multi <- read.table('/data/clusterfs/lag/users/sousoh/ukbb/lists/list-t1-all-available.txt',header=F)[,1]

master_df<-data.frame(contrast=character(0),XCA=character(0),dim=integer(0),ic_num=integer(0),lambda=double(0),delta=double(0),lassosum_rsq=double(0),lassosum_tval=double(0),sbayesR_rsq=double(0),sbayesR_tval=double(0))

for (contrast in c("tbm","dmri") ){
	print (paste("Running contrast:",contrast))
	if (contrast=="dmri"){subids<-subids_dmri}else{subids<-subids_tbm}
	pgs_img <- scale(a[ match(as.integer(subids) ,as.integer(pgs_subid) ),])
	sbayesr_pgs_img <- scale(sbayesr_score[ match(as.integer(subids) ,as.integer(pgs_subid) ) ])

	for (xca in c("pca","ica") ){
		print(paste("Running decomposition:",xca))
               for (xca_dim in c(324,200,124,76,47,29,18,11) ){
			print(paste("Running dimension",xca_dim))
			if (contrast=="dmri"){
				pheno <- read.table(paste(sep="",'/data/clusterfs/lag/users/sousoh/ukbb/genetic/design-matrix-dMRI/ica-dim-',
						  as.character(xca_dim),'-migp-8000/melodic_mix_',xca,'_dim_',as.character(xca_dim),'.txt')    ,header=F)
			}else{ 	pheno <- read.table(paste(sep="",'/data/clusterfs/lag/users/sousoh/ukbb/jac/merged/new-ica-dim-',
						  as.character(xca_dim),'-migp-8000/melodic_mix_',xca,'_dim_',as.character(xca_dim),'.txt')    ,header=F) }

if (contrast=="tbm"){ pheno<-pheno[ match(subids_tbm,subids_tbm_multi) ,] }
			for (comp_id in seq(1:xca_dim) ){
if (xca=="pca"){pid = xca_dim - comp_id +1} else {pid=comp_id}
ic_name=paste(sep="",contrast,"_",xca,"_",comp_id,"_dim_",xca_dim)
pheno_pc <- pheno[,pid]
if (contrast=="tbm"){
l1<-paste("COV", 1:24, sep = "", collapse = "+") %>%  paste0("pheno_pc","~", .) %>% as.formula  %>% lm(.,data=cbind(pheno_pc,covariates_tbm)) }
else
{ l1<-paste("COV", 1:24, sep = "", collapse = "+") %>%  paste0("pheno_pc","~", .) %>% as.formula  %>% lm(.,data=cbind(pheno_pc,covariates_dMRI)) }
ic_corr <- scale(l1$residuals)

sbayesr_lm <- summary(lm(ic_corr ~ sbayesr_pgs_img)) 
#span all grid params
for (i in 1:dim(b)[1]){
if (!is.na(a[1,i]) ) { l<-summary(lm(ic_corr ~ pgs_img[,i])) }
	master_df[nrow(master_df)+1,]<- c(contrast,xca,xca_dim,comp_id,b[i,1],b[i,2],l$r.squared,l$coefficients[2,3],sbayesr_lm$r.squared, sbayesr_lm$coefficients[2,3]) 
}
#print(paste("IC pid",pid))
comp_ids <- which(master_df$contrast==contrast & master_df$XCA==xca & master_df$dim==xca_dim & master_df$ic_num==comp_id)

print(paste("Component",comp_id,"IC pid", pid, "sbayesR rsq",format(round(sbayesr_lm$r.squared,5),nsmall=5), "sbayesR t-value",format(round(sbayesr_lm$coefficients[2,3],2),nsmall=2),
	    " MAX rsq ", format(round(max(as.double(master_df[comp_ids,]$lassosum_rsq)),5),nsmall=5), " MAX t-value ",format(round(max(abs(as.double(master_df[comp_ids,]$lassosum_tval))),2),nsmall=2) ))
} # IC number loop
comp_ids <- which(master_df$contrast==contrast & master_df$XCA==xca & master_df$dim==xca_dim )
print (paste("Finished exploring all components and grid parameters, max rsq is",      format(round(max(as.double(master_df[comp_ids,]$lassosum_rsq)),5),nsmall=5)           ))
if (xca=="pca") {break} 
} # dimension loop, we only have one dimension for PCA hence breaking out		 
} # PCA vs ICA loop
} # contrast loop (dMRI vs TBM)

master_df[,3:10] <- sapply(master_df[,3:10],as.numeric)
dir.create(file.path("/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs", cond,"component-results"), showWarnings = FALSE)
write.table(master_df,file=paste(sep="","/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/",cond,"/component-results/master_df.txt"),quote=F,row.names=F)
#master_df <- read.table(paste(sep="","/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/",cond,"/component-results/master_df.txt"),header=T)
master_df[,3:10] <- sapply(master_df[,3:10],as.numeric)
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

#contrast="dmri"; xca="ica";xca_dim=200;comp_id=30
#comp_ids <- which(master_df$contrast==contrast & master_df$XCA==xca & master_df$dim==xca_dim & master_df$ic_num==comp_id)
comp_ids_tbm <- which(master_df$contrast=="tbm" & master_df$XCA==master_df[max_ic_tbm,]$XCA & master_df$dim==master_df[max_ic_tbm,]$dim & master_df$ic_num==master_df[max_ic_tbm,]$ic_num)
comp_ids_dmri <- which(master_df$contrast=="dmri" & master_df$XCA==master_df[max_ic_dmri,]$XCA & master_df$dim==master_df[max_ic_dmri,]$dim & master_df$ic_num==master_df[max_ic_dmri,]$ic_num)
max_y <- max(master_df[comp_ids_tbm,]$lassosum_rsq, master_df[comp_ids_dmri,]$lassosum_rsq)

comp_ids <- comp_ids_tbm

if (master_df[max_ic_tbm,]$XCA == "ica") {dim_name=paste(sep="","-",master_df[max_ic_tbm,]$dim)}else{dim_name=""}

gg_tbm <- ggplot(master_df[comp_ids,]) + geom_point(aes(x=log(lambda),color=factor(round(delta,2)),y=lassosum_rsq)) + 
	geom_line(aes(x=log(lambda),color=factor(round(delta,2)),y=lassosum_rsq,linetype="Lassosum2")) +
	xlab("Lambda (L₁) regularisation") + ylab("IDP explained variance (R²)") + 
	ylim(c(0,1.1*max_y))   + theme(legend.position = "none") + 
	ggtitle(paste(sep="", "sMRI TBM\n", toupper(master_df[max_ic_tbm,]$XCA),dim_name, " component #",master_df[max_ic_tbm,]$ic_num)) + 
	geom_hline(aes(yintercept=master_df[comp_ids_tbm,]$sbayesR_rsq[1],linetype="SBayesR" ),color="black",linewidth=1,alpha=0.5)+ 
	 scale_linetype_manual("PGS type",values=c("Lassosum2"="solid","SBayesR"="dashed")) +  guides(linetype=guide_legend(keywidth = 3, keyheight = 1  )) + theme(legend.key=element_blank()) +
	geom_label_repel(data=master_df[comp_ids,][master_df[comp_ids,]$lassosum_rsq==max(master_df[comp_ids,]$lassosum_rsq),],aes(x=log(lambda),y=lassosum_rsq, label= paste("R² =", round(max(master_df[comp_ids_tbm,]$lassosum_rsq),4)) ), nudge_x=0.5, nudge_y=25e-5  ) # + guides(colour=guide_legend(title="Delta regularisation (L2)"))  

comp_ids <- comp_ids_dmri

if (master_df[max_ic_dmri,]$XCA == "ica") {dim_name=paste(sep="","-",master_df[max_ic_dmri,]$dim)}else{dim_name=""}
gg_dmri <- ggplot(master_df[comp_ids,]) + geom_point(aes(x=log(lambda),color=factor(round(delta,2)),y=lassosum_rsq)) + 
	geom_line(aes(x=log(lambda),color=factor(round(delta,2)),y=lassosum_rsq,linetype="Lassosum2"   )) +
	guides(colour=guide_legend(title="Delta (L₂) regularisation")) + xlab("Lambda (L₁) regularisation") + ylab("") + #ylab("IDP explained variance (R²)") + 
	ylim(c(0,1.1*max_y)) + theme(legend.position = "none") +
	ggtitle(paste(sep="", "dMRI AFD\n", toupper(master_df[max_ic_dmri,]$XCA), dim_name, " component #",master_df[max_ic_dmri,]$ic_num )) +
	geom_hline(aes(yintercept=master_df[comp_ids_dmri,]$sbayesR_rsq[1], linetype="SBayesR"),color="black",linewidth=1,alpha=0.5) +   
	scale_linetype_manual("PGS type",values=c("Lassosum2"="solid","SBayesR"="dashed")) +  guides(linetype=guide_legend(keywidth = 3, keyheight = 1  )) + theme(legend.key=element_blank()) +
	geom_label_repel(data=master_df[comp_ids,][master_df[comp_ids,]$lassosum_rsq==max(master_df[comp_ids,]$lassosum_rsq),],aes(x=log(lambda),y=lassosum_rsq, label= paste("R² =", round(max(master_df[comp_ids_dmri,]$lassosum_rsq),4)) ), nudge_x=0.5, nudge_y=25e-5  ) 



cond_capital <- paste(toupper(substr(cond, 1, 1)), substr(cond, 2, nchar(cond)), sep="")

gg_dmri_legend <- ggplot(master_df[comp_ids,]) + geom_point(aes(x=log(lambda),color=factor(round(delta,2)),y=lassosum_rsq)) + 
	geom_line(aes(x=log(lambda),color=factor(round(delta,2)),y=lassosum_rsq ,linetype="Lassosum2"  )) +
	guides(colour=guide_legend(title="Delta\nregularisation (L₂)"))+ 
	xlab("Lambda regularisation (L₁)") + ylab("IDP explained variance (R²)") + ylim(c(0,1.1*max_y)) + 
	geom_hline(aes(yintercept=master_df[comp_ids,]$sbayesR_rsq[1] ,linetype="SBayesR" ),color="black",linewidth=1,alpha=0.5) + 
	   scale_linetype_manual( paste(sep="", cond_capital  , "\nPGS type") ,values=c("Lassosum2"="solid","SBayesR"="dashed")) +   
	   guides(linetype=guide_legend(keywidth = 3 , override.aes=list(alpha=1, linewidth=1 ) ) ) + 
	   theme(legend.key=element_blank()) +
	geom_label_repel(data=master_df[comp_ids,][master_df[comp_ids,]$lassosum_rsq==max(master_df[comp_ids,]$lassosum_rsq),],aes(x=log(lambda),y=lassosum_rsq, label= paste("R²:", round(max(master_df[comp_ids_dmri,]$lassosum_rsq),4)) ), nudge_x=0.5, nudge_y=1e-4  ) 


gg_legend <- get_legend(gg_dmri_legend, position = NULL)

final_plot <- plot_grid(gg_tbm,gg_dmri, gg_legend, ncol=3, align="hv",  rel_widths = c(2, 2,1))
ggsave(paste(sep="","/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/",cond,"/r2_plot.png" ) ,plot=final_plot, bg="white")
ggsave(paste(sep="","/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/",cond,"/r2_plot.svg" ) ,plot=final_plot, bg="white")

