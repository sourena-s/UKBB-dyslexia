library(bigsnpr)
library(data.table)
library(magrittr)

options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

sumstats <- read.table("/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/dyslexia/dyslexia.filtered.ldpred2",header=T)
names(sumstats) <-    c("chr","pos","rsid","a1","a0","n_eff","beta_se","p", "OR","INFO","MAF","beta")

#becareful with AD, bec sumstats are in GR38 coordinates
#sumstats <- read.table("/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/AD/sumstats.ldpred2",header=T)
#names(sumstats) <- c("rsid","chr","pos", "a1", "a0", "MAF", "beta", "beta_se", "p", "n_eff", "INFO")
sumstats$chr <- as.character(sumstats$chr)
sumstats <- sumstats[sumstats$chr!="X",]
sumstats$chr <- as.integer(sumstats$chr)

#only keep those variants with MAX sample sizeaverage
#sumstats <- sumstats[sumstats$n_eff==max(sumstats$n_eff),]
sumstats$chr <- as.integer(sumstats$chr)

ld<-readRDS('/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/matrix_big40_ld_hapmap.rds')
corr<-readRDS('/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/matrix_big40_corr_hapmap.rds')
obj.bigSNP <- snp_attach("/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/big40_subset_2k_hapmap_plus.rds")
obj.bigSNP.all <- snp_attach("/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/big40.rds")

fam.order <- NULL
fam.order<-read.table('/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/list-all-subjects-genetics-subid-and-indices.txt',header=F)[,1]
fam.order <- as.data.frame(cbind(fam.order,fam.order,fam.order,fam.order))
names(fam.order) <- c("family.ID", "sample.ID","FID", "IID")
obj.bigSNP$fam <- fam.order
map <- obj.bigSNP$map[,c(1,3,4,5,6)]
names(map) <- c("chr", "rsid", "pos", "a1", "a0")
map$pos<- as.integer(map$pos)
map$chr<- as.integer(map$chr)
info_snp <- snp_match(sumstats, join_by_pos = FALSE, map, return_flip_and_rev=TRUE, strand_flip =FALSE)
df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]

#now find variant indices on the large big40 genotype matrix, this is to be used in Pred matrix multiplicationa AFTER PGS estimation
#note that variant count across traits (i.e. var_ids) may change depending on the sumstat file used and its overlap with hapmap3/UKBB vars
big.map<-as.data.frame(obj.bigSNP.all$map)
colnames(big.map) <- c("chr","marker","rsid","pos","a1","a0","freq","info")
big.map$chr <- as.integer(big.map$chr)
var_ids <- snp_match(info_snp,big.map, join_by_pos = FALSE, return_flip_and_rev=TRUE, strand_flip =FALSE)
var_ids <- var_ids$'_NUM_ID_'

ldsc <- snp_ldsc(ld[info_snp$'_NUM_ID_'], length(ld[info_snp$'_NUM_ID_']), chi2 = (df_beta$beta / df_beta$beta_se)^2, sample_size = df_beta$n_eff, blocks = NULL)
h2_est <- ldsc[["h2"]]
print(paste("h2_est is",h2_est))
##lassosum2
#tolerance is 1e-5 by default, increasing to see if we can have more convergence
#delta:L2  lambda:L1
beta_lassosum2 <- snp_lassosum2(corr, df_beta, ncores = 64, ind.corr=info_snp$'_NUM_ID_',
delta=seq_log(0.01,10000,10),nlambda=30,lambda.min.ratio=0.01, maxiter=1000,tol=1e-2)
pred_lasso <- big_prodMat(obj.bigSNP.all$genotype, beta_lassosum2, ind.col = var_ids,  ncores = 64, block.size=5000 )

cond="dyslexia" # or AD

saveRDS(beta_lassosum2,file=paste(sep="","/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/",cond,"/lasso/beta_lassosum2-2.rds"))
write.table(pred_lasso,paste(sep="","/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/",cond, "/lasso/lasso-2.pgs"), quote=F,row.names=F, col.names=F,sep=",")
write.table(attr(beta_lassosum2, "grid_param") ,paste(sep="","/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/",cond,"/lasso/grid-2.params"), quote=F,row.names=F,sep=",")

##auto model
set.seed(1)
n_iter=1000
n_burn=250
multi_auto_no_alpha <- snp_ldpred2_auto(corr, burn_in=n_burn, num_iter = n_iter, df_beta,verbose=TRUE, sparse=TRUE, h2_init = h2_est  , #use_MLE=FALSE, allow_jump_sign=FALSE, 
		shrink_corr = 0.95, vec_p_init = seq_log(0.01, 1, length.out =25), ncores = 64, ind.cor = info_snp$`_NUM_ID_`)

saveRDS(multi_auto_no_alpha,file=paste(sep="","/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/",cond,"/ldpred2/multi_auto_no_alpha_mle.rds"))
#multi_auto_alpha <- snp_ldpred2_auto(corr, burn_in=n_burn, num_iter = n_iter, df_beta,verbose=TRUE, sparse=TRUE, h2_init = h2_est, use_MLE=TRUE, allow_jump_sign=FALSE, 
#		shrink_corr = 0.95, vec_p_init = seq_log(0.01, 1, length.out =25), ncores = 64, ind.cor = info_snp$`_NUM_ID_`)

range <- sapply(multi_auto_no_alpha, function(auto) diff(range(auto$corr_est)))
keep <- (range > (0.90 * quantile(range, 0.90, na.rm = TRUE )))
beta_auto_dense <-  rowMeans(sapply(multi_auto_no_alpha[keep], function(auto) auto$beta_est       ))
beta_auto_sparse <- rowMeans(sapply(multi_auto_no_alpha[keep], function(auto) auto$beta_est_sparse))

pred_auto <- big_prodMat(obj.bigSNP.all$genotype, beta_auto_dense, ind.col=var_ids) 
pred_scaled <- apply(pred_auto, 2, sd)
(abs(pred_scaled - median(pred_scaled)) < 3 * mad(pred_scaled))
final_beta_auto_dense <- rowMeans(beta_auto_dense[,abs(pred_scaled - median(pred_scaled)) < 3 * mad(pred_scaled)])
final_pred <- big_prodMat(obj.bigSNP.all$genotype, final_beta_auto_dense, ind.col=var_ids)

pred_auto <- big_prodMat(obj.bigSNP.all$genotype, beta_auto_sparse, ind.col=var_ids) 
pred_scaled <- apply(pred_auto, 2, sd)
final_beta_auto <- rowMeans(beta_auto_sparse[,abs(pred_scaled - median(pred_scaled)) < 3 * mad(pred_scaled)])

ind.valpostp <- rowMeans(sapply(multi_auto_no_alpha[keep], function(auto) auto$postp_est))

#Grid
h2_seq <- round(h2_est *  c(0.1,0.2,0.3,0.4,0.7,1,1.4) , 4)
p_seq <- signif(seq_log(1e-5, 0.9, length.out = 10), 2)
#apparently each h2 and sparseness parameter (true/false) gets its own thread
params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))
beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = 64, burn_in=n_burn, num_iter = n_iter, ind.corr=info_snp$'_NUM_ID_')
pred_grid <- big_prodMat(obj.bigSNP.all$genotype, beta_grid, ind.col = var_ids,  ncores = 64, block.size=5000 )


saveRDS(multi_auto_no_alpha,file=paste(sep="","/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/",cond,"/ldpred2/multi_auto_no_alpha.rds"))
write.table(params,paste(sep="","/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/",cond,"/ldpred2/grid.params"), quote=F,row.names=F,sep=",")
write.table(beta_grid,paste(sep="","/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/",cond,"/ldpred2/grid.beta"), quote=F,row.names=F, col.names=F,sep=",")
write.table(pred_grid,paste(sep="","/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/",cond, "/ldpred2/grid.pgs"), quote=F,row.names=F, col.names=F,sep=",")

write.table(fam.order$family.ID,paste(sep="","/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/",cond,"/subject_id"), quote=F,row.names=F, col.names=F,sep=",")


params <- read.table("/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/AD/ldpred2/grid.params",header=T,sep=",")

chain_h2 <- matrix(nrow=length(multi_auto_no_alpha),ncol=n_burn+n_iter)
chain_p <- matrix(nrow=length(multi_auto_no_alpha),ncol=n_burn+n_iter)

chain_p_avg <- vector()
chain_h2_avg <- vector()

for (i in 1:length(multi_auto_no_alpha) ) { 
chain_p_avg[i] <- mean(multi_auto_no_alpha[[i]]$path_p_est[(n_burn+1):(n_burn+n_iter)])
chain_h2_avg[i] <- mean(multi_auto_no_alpha[[i]]$path_h2_est[(n_burn+1):(n_burn+n_iter)])
chain_p[i,] <-  multi_auto_no_alpha[[i]]$path_p_est
chain_h2[i,] <-  multi_auto_no_alpha[[i]]$path_h2_est}

write.table(chain_p,"/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/AD/ldpred2/chain_p.csv",sep=",",quote=F,row.names=F)
write.table(chain_h2,"/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/AD/ldpred2/chain_h2.csv",sep=",",quote=F,row.names=F)

chain_h2 <- read.csv('/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/AD/ldpred2/chain_h2.csv',header=T)
chain_p <- read.csv('/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/AD/ldpred2/chain_p.csv',header=T)
chain_df <- data.frame(matrix(nrow=0,ncol=4))
for (i in 1:dim(chain_h2)[1]) { for (j in 1:dim(chain_h2)[2]) {
chain_df <- rbind(chain_df,c(i,j,chain_h2[i,j], chain_p[i,j]) ) }}
colnames(chain_df) <- c("chain_ID","Iteration","h2","p")



(params2 <- attr(beta_lassosum2, "grid_param"))
pred_grid2 <- big_prodMat(obj.bigSNP.all$genotype, beta_lassosum2, ind.col = var_ids,  ncores = 64, block.size=5000 )
params2$score <- apply(pred_grid2[ind.val, ], 2, function(x) {
  if (all(is.na(x))) return(NA)
  summary(lm(y[ind.val] ~ x))$coef["x", 3]
  # summary(glm(y[ind.val] ~ x, family = "binomial"))$coef["x", 3]
})

library(ggplot2)
ggplot(params2, aes(x = lambda, y = score, color = as.factor(delta))) +
  theme_bigstatsr() +
  geom_point() +
  geom_line() +
  scale_x_log10(breaks = 10^(-5:0)) +
  labs(y = "GLM Z-Score", color = "delta") +
  theme(legend.position = "top") +
  guides(colour = guide_legend(nrow = 1))

library(dplyr)
best_grid_lassosum2 <- params2 %>%
  mutate(id = row_number()) %>%
  arrange(desc(score)) %>%
  print() %>% 
  slice(1) %>%
  pull(id) %>% 
  beta_lassosum2[, .]

# Choose the best among all LDpred2-grid and lassosum2 models
best_grid_overall <- 
  `if`(max(params2$score, na.rm = TRUE) > max(params$score, na.rm = TRUE),
       best_grid_lassosum2, best_beta_grid)




#In the LDpred2 paper, we proposed an automatic way of filtering bad chains by comparing the scale of the resulting predictions. We have tested a somewhat equivalent and simpler alternative since, which we recommend here:

# `range` should be between 0 and 2

pcor(pred_auto, y[ind.test], NULL)

library(ggplot2)
auto <- multi_auto[[1]]  # first chain
plot_grid(
  qplot(y = auto$path_p_est) +
    theme_bigstatsr() +
    geom_hline(yintercept = auto$p_est, col = "blue") +
    scale_y_log10() +
    labs(y = "p"),
  qplot(y = auto$path_h2_est) +
    theme_bigstatsr() +
    geom_hline(yintercept = auto$h2_est, col = "blue") +
    labs(y = "h2"),
  ncol = 1, align = "hv"
)
