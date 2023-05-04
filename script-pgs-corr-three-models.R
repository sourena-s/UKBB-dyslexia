library(data.table)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(bigsnpr)
library(ggpubr)
library(tidyr)
library(ggnewscale)
core_count=16
args = commandArgs(trailingOnly=TRUE)
cond=as.character(args[1])
obj.bigSNP.all <- snp_attach("/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/processed/big40_all_subj_subset_hapmap_2k.rds")
map <- as.data.frame(cbind(obj.bigSNP.all$map$chromosome,obj.bigSNP.all$map$rsid,obj.bigSNP.all$map$physical.pos,obj.bigSNP.all$map$allele1, obj.bigSNP.all$map$allele2))
#check ref allele
names(map) <- c("chr","rsid","pos","a0","a1")
map$chr<-as.integer(map$chr)

#todo: check if file exists
sbayesr_data <-  read.table(paste(sep="",'/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/',cond, '/auto-sbayesr.pgs.snpRes'),header=T)
names(sbayesr_data) <- c("id", "rsid","chr", "pos","a1","a0","freq","beta","beta_se","PIP","LastSampleEff")

sbayesr_match <- snp_match(sbayesr_data, join_by_pos = FALSE, map, return_flip_and_rev=TRUE, strand_flip =FALSE)
sbayesr_score <- big_prodMat(obj.bigSNP.all$genotype, as.matrix(sbayesr_match$beta), ind.col = sbayesr_match$'_NUM_ID_',  ncores = core_count, block.size=5000 )

write.table(sbayesr_score,        file=paste(sep="",'/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/',cond, '/auto-sbayesr.pgs'),quote=F,col.names=F,row.names=F)
####end todo



pgs_subid<-read.table('/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/list-all-subjects-genetics-subid-and-indices.txt',header=F)[,1]

##PRS-CS
cs_score <- list()
cs_lm <- list()
cs_pgs_img <- list()

obj.bigSNP.all.cs <- snp_attach("/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/big40.rds")
map.cs <- as.data.frame(cbind(obj.bigSNP.all.cs$map$chromosome,obj.bigSNP.all.cs$map$rsid,obj.bigSNP.all.cs$map$physical.pos,obj.bigSNP.all.cs$map$allele1, obj.bigSNP.all.cs$map$allele2))
#check ref allele
names(map.cs) <- c("chr","rsid","pos","a0","a1")
map.cs$chr<-as.integer(map.cs$chr)
map.cs$pos<-as.integer(map.cs$pos)

for (cs_id in 1:5) {
phi <- c("auto","1e-06", "1e-04", "1e-02", "1e+00")[cs_id]
cs_data <-  read.table(paste(sep="",'/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/',cond, '/prs-cs/pgs-autosome-phi-',phi,'.txt'),header=F)
names(cs_data) <- c("chr", "rsid","pos","a1","a0","beta")
cs_match <- snp_match(cs_data, join_by_pos = FALSE, map.cs, return_flip_and_rev=T, strand_flip =F)
cs_score[[cs_id]] <- big_prodMat(obj.bigSNP.all.cs$genotype, as.matrix(cs_match$beta), ind.col = cs_match$'_NUM_ID_',  ncores = core_count, block.size=5000 )
cs_score[[cs_id]] <- as.vector(scale(cs_score[[cs_id]]) )
write.table(cs_score[[cs_id]], file=paste(sep="",'/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/',cond, '/prs-cs/subj_pgs_phi_',phi,'.pgs'),quote=F,col.names=F,row.names=F)
#cs_score[[cs_id]] <- read.table(paste(sep="",'/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/',cond, '/prs-cs/subj_pgs_phi_',phi,'.pgs'),header=F)[,1]
}

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

master_df<-data.frame(contrast=character(0),XCA=character(0),dim=integer(0),ic_num=integer(0),lambda=double(0),delta=double(0),lassosum_rsq=double(0),lassosum_tval=double(0),sbayesR_rsq=double(0),sbayesR_tval=double(0),cs1_rsq=double(0),cs1_tval=double(0) ,cs2_rsq=double(0),cs2_tval=double(0),cs3_rsq=double(0),cs3_tval=double(0),cs4_rsq=double(0),cs4_tval=double(0),cs5_rsq=double(0),cs5_tval=double(0) )

for (contrast in c("dmri","tbm") ){
	print (paste("Running contrast:",contrast))
	if (contrast=="dmri"){subids<-subids_dmri}else{subids<-subids_tbm}
	pgs_img <- scale(a[ match(as.integer(subids) ,as.integer(pgs_subid) ),])
	sbayesr_pgs_img <- scale(sbayesr_score[ match(as.integer(subids) ,as.integer(pgs_subid) ) ])
	for (cs_id in 1:5) {
	cs_pgs_img[[cs_id]] <- scale(cs_score[[cs_id]][ match(as.integer(subids) ,as.integer(pgs_subid) ) ]) }

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
for (cs_id in 1:5) {
cs_lm[[cs_id]] <-  summary(lm(ic_corr ~ cs_pgs_img[[cs_id]]))  }

#span all grid params
for (i in 1:dim(b)[1]){
if (!is.na(a[1,i]) ) { l<-summary(lm(ic_corr ~ pgs_img[,i])) }
	master_df[nrow(master_df)+1,]<- c(contrast,xca,xca_dim,comp_id,b[i,1],b[i,2],l$r.squared,l$coefficients[2,3],sbayesr_lm$r.squared, sbayesr_lm$coefficients[2,3] , cs_lm[[1]]$r.squared, cs_lm[[1]]$coefficients[2,3],  cs_lm[[2]]$r.squared, cs_lm[[2]]$coefficients[2,3],  cs_lm[[3]]$r.squared, cs_lm[[3]]$coefficients[2,3], cs_lm[[4]]$r.squared, cs_lm[[4]]$coefficients[2,3], cs_lm[[5]]$r.squared, cs_lm[[5]]$coefficients[2,3])
}
#print(paste("IC pid",pid))
comp_ids <- which(master_df$contrast==contrast & master_df$XCA==xca & master_df$dim==xca_dim & master_df$ic_num==comp_id)

print(paste("Component",comp_id,"IC pid", pid, "sbayesR rsq",format(round(sbayesr_lm$r.squared,5),nsmall=5), "sbayesR t-value",format(round(sbayesr_lm$coefficients[2,3],2),nsmall=2),
	    " MAX rsq ", format(round(max(as.double(master_df[comp_ids,]$lassosum_rsq)),5),nsmall=5), " MAX t-value ",format(round(max(abs(as.double(master_df[comp_ids,]$lassosum_tval))),2),nsmall=2)        , " PRS-CS t-values [auto & phi increments]",        format(round(max(as.double(master_df[comp_ids,]$cs1_tval)),2),nsmall=2) ,   format(round(max(as.double(master_df[comp_ids,]$cs2_tval)),2),nsmall=2) ,  format(round(max(as.double(master_df[comp_ids,]$cs3_tval)),2),nsmall=2) ,  format(round(max(as.double(master_df[comp_ids,]$cs4_tval)),2),nsmall=2) ,  format(round(max(as.double(master_df[comp_ids,]$cs5_tval)),2),nsmall=2)      ))
} # IC number loop
comp_ids <- which(master_df$contrast==contrast & master_df$XCA==xca & master_df$dim==xca_dim )
print (paste("Finished exploring all components and grid parameters, max rsq is",      format(round(max(as.double(master_df[comp_ids,]$lassosum_rsq)),5),nsmall=5)           ))
if (xca=="pca") {break} 
} # dimension loop, we only have one dimension for PCA hence breaking out		 
} # PCA vs ICA loop
} # contrast loop (dMRI vs TBM)

master_df[,3:20] <- sapply(master_df[,3:20],as.numeric)
dir.create(file.path("/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs", cond,"component-results"), showWarnings = FALSE)
write.table(master_df,file=paste(sep="","/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/",cond,"/component-results/master_df.txt"),quote=F,row.names=F)

