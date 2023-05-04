
##auto model
set.seed(1)
n_iter=1000
n_burn=250
multi_auto_no_alpha <- snp_ldpred2_auto(corr, burn_in=n_burn, num_iter = n_iter, df_beta,verbose=TRUE, sparse=TRUE, h2_init = h2_est  , #use_MLE=FALSE, allow_jump_sign=FALSE,
                shrink_corr = 0.95, vec_p_init = seq_log(0.01, 1, length.out =25), ncores = 64, ind.cor = info_snp$`_NUM_ID_`)

saveRDS(multi_auto_no_alpha,file=paste(sep="","/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/",cond,"/ldpred2/multi_auto_no_alpha_mle.rds"))
#multi_auto_alpha <- snp_ldpred2_auto(corr, burn_in=n_burn, num_iter = n_iter, df_beta,verbose=TRUE, sparse=TRUE, h2_init = h2_est, use_MLE=TRUE, allow_jump_sign=FALSE,
#               shrink_corr = 0.95, vec_p_init = seq_log(0.01, 1, length.out =25), ncores = 64, ind.cor = info_snp$`_NUM_ID_`)

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

