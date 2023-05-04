library(bigsnpr)
library(data.table)
library(magrittr)

options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

#sumstats <- read.table("/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/dyslexia/dyslexia.filtered.ldpred2",header=T)
#names(sumstats) <-    c("chr","pos","rsid","a1","a0","n_eff","beta_se","p", "OR","INFO","MAF","beta")

#becareful with AD, bec sumstats are in GR38 coordinates
sumstats <- read.table("/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/AD/sumstats.ldpred2",header=T)
names(sumstats) <- c("rsid","chr","pos", "a1", "a0", "MAF", "beta", "beta_se", "p", "n_eff", "INFO")
# LDpred 2 require the header to follow the exact naming

sumstats$chr <- as.character(sumstats$chr)
sumstats <- sumstats[sumstats$chr!="X",]
sumstats$chr <- as.integer(sumstats$chr)

#only keep those variants with MAX sample size (for debugging purpose)
sumstats <- sumstats[sumstats$n_eff==max(sumstats$n_eff),]

# Filter out hapmap SNPs hapmap3: 
#info <- readRDS("/data/clusterfs/lag/users/sousoh/ukbb/genetic/hapmap/hapmap3plus/map_hm3_plus.rds")
#top_variants<-read.table('/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/dyslexia/dyslexia.ldpred2.top-variants.pval',header=F,colClasses = "character")[,1]
#sumstats <- sumstats[sumstats$rsid %in% c(info$rsid),]

sumstats$chr <- as.integer(sumstats$chr)

#more filtering sumstats
#sumstats <- sumstats[  (sumstats$beta > mean(sumstats$beta) - 3* sd(sumstats$beta)) & (sumstats$beta < mean(sumstats$beta) + 3* sd(sumstats$beta)  ) ,]

NCORES <- 12
# Open a temporary file
tmp <- tempfile(tmpdir = "tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
# Initialize variables for storing the LD score and LD matrix
corr <- NULL
ld <- NULL
# We want to know the ordering of samples in the bed file
fam.order <- NULL

#obj.bigSNP <- snp_attach("/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/big40.rds")
#random_2k <-read.table('/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/list-all-subjects-genetics-subid-and-indices-white-British-suff-2k.txt',header=F)[,2]
#subset_file <- snp_subset(obj.bigSNP,ind.row=random_2k, backingfile="/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/big40_subset_2k")
obj.bigSNP.2k <- snp_attach("/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/big40_subset_2k.rds")
fam.order<-read.table('/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/list-all-subjects-genetics-subid-and-indices.txt',header=F)[,1]
fam.order.2k <- read.table('/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/list-all-subjects-genetics-subid-and-indices-white-British-suff-2k.txt',header=F)[,1]
fam.order <- as.data.frame(cbind(fam.order,fam.order,fam.order,fam.order))
fam.order.2k <- as.data.frame(cbind(fam.order.2k,fam.order.2k,fam.order.2k,fam.order.2k))
names(fam.order) <- c("family.ID", "sample.ID","FID", "IID")
names(fam.order.2k) <- c("family.ID", "sample.ID","FID", "IID")
#obj.bigSNP$fam <- fam.order
obj.bigSNP.2k$fam <- fam.order.2k

# extract the SNP information from the genotype
map <- obj.bigSNP.2k$map[,c(1,3,4,5,6)]
names(map) <- c("chr", "rsid", "pos", "a1", "a0")
map$pos<- as.integer(map$pos)
map$chr<- as.integer(map$chr)
# perform SNP matching
#sumstats<-sumstats[sample(1:dim(sumstats)[1],dim(sumstats)[1] %/% 50,replace=F),]
info_snp <- snp_match(sumstats, join_by_pos = FALSE, map, return_flip_and_rev=TRUE, strand_flip =FALSE)

# better to use af of GWAS and INFO scores as well (then can use 0.7 instead of 0.5)
#sd_ldref <- with(info_snp, sqrt(2 * af_UKBB * (1 - af_UKBB)))
#sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2))
#is_bad <-
#  sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.05 | sd_ldref < 0.05

# Assign the genotype to a variable for easier downstream analysis
#subset_file <- snp_subset(obj.bigSNP.2k,ind.col=info_snp$'_NUM_ID_')
#genotype <- snp_attach(subset_file)
genotype <- snp_attach('/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/big40_subset_2k_sub1.rds')
#Brace for vertigo
genotype <- genotype$genotypes
# Rename the data structures
CHR <- map$chr
POS <- map$pos
# get the CM information from 1000 Genome
POS2 <- snp_asGeneticPos(CHR, POS, dir = "/data/clusterfs/lag/users/sousoh/ukbb/genetic/hapmap" ,type="hapmap" )
# calculate LD
for (chr in 1:22) {
    # Extract SNPs that are included in the chromosome
    ind.chr <- which(info_snp$chr == chr)
#we already filtered genotype based on NUM_ID, so no need for a second ind2 query
ind.chr2 <- ind.chr
#    ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
    corr0 <- snp_cor(genotype, ind.col = ind.chr2, ncores = 1, infos.pos = POS2[ind.chr2], size= 3/1000)
    if (chr == 1) {
    	    ld <- Matrix::colSums(corr0^2)
        corr <- as_SFBM(corr0, tmp)
    } else {
        ld <- c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
    }
    print (paste("Chromosome", chr))
}

saveRDS(ld, file = "/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/ld_corr_dysl.rds'") 
quit(save="no")

df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
ldsc <- snp_ldsc(ld, length(ld), chi2 = (df_beta$beta / df_beta$beta_se)^2, sample_size = df_beta$n_eff, blocks = NULL)
h2_est <- ldsc[["h2"]]
print("h2_est is",h2_est)
##auto model
n_iter=1000
n_burn=200
multi_auto <- snp_ldpred2_auto(corr, 
	report_step=1000,burn_in=n_burn, num_iter = n_iter,
    df_beta,verbose=TRUE, sparse=TRUE,
    h2_init = h2_est, use_MLE=FALSE, allow_jump_sign=FALSE, 
    vec_p_init = seq_log(0.02, 0.9, length.out =5),
    ncores = 4)

beta_auto <- sapply(multi_auto, function(auto)
    auto$beta_est)


pred_auto <-
    big_prodMat(genotype,
                beta_auto,
                ind.row = ind.test, ind.col = info_snp$`_NUM_ID_`)
# scale the PRS generated from AUTO
pred_scaled <- apply(pred_auto, 2, sd)
final_beta_auto <-
    rowMeans(beta_auto[,
                abs(pred_scaled -
                    median(pred_scaled)) <
                    3 * mad(pred_scaled)])
#Grid
h2_seq <- round(h2_est *  seq_log(0.2,5,5) , 4)
p_seq <- signif(seq_log(1e-5, 1, length.out = 5), 2)
params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))
#params <- params[sample(1:dim(params)[1],30,replace=F),]
beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = 1, burn_in=20, num_iter = 100)




# Reformat the phenotype file such that y is of the same order as the
# sample ordering in the genotype file
#fam.order <- as.data.table(obj.bigSNP$fam)
subids <- read.table('/data/clusterfs/lag/users/sousoh/ukbb/genetic/design-matrix-dMRI/bases/c0-subjects-path-QC-passed-only-indices',header=F)
pcdim=200
icid=192
ica_mix <- read.table(paste(sep="",'/data/clusterfs/lag/users/sousoh/ukbb/genetic/design-matrix-dMRI/ica-dim-',as.character(pcdim),'-migp-8000/melodic_mix_pca_dim_',as.character(pcdim),'.txt')	,header=F)


phenotype <- as.data.frame(cbind(subids,subids,ica_mix[,icid]))
colnames(phenotype) <- c("FID","IID", paste(sep="","PC",icid,"DIM",pcdim) )




#y <- pheno[fam.order, on = c("FID", "IID")]


pheno <- pheno[match(fam.order[,3],pheno[,1]),]
# Prepare data for grid model
p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2)
h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4)
grid.param <-
    expand.grid(p = p_seq,
            h2 = h2_seq,
            sparse = c(FALSE, TRUE))
# Get adjusted beta from grid model
beta_grid <- snp_ldpred2_grid(corr, df_beta, grid.param, ncores = NCORES)


genotype <- obj.bigSNP$genotypes
# calculate PRS for all samples
ind.test <- 1:nrow(genotype)

pred_grid <- big_prodMat(   genotype, 
                            beta_grid, 
                            ind.col = info_snp$`_NUM_ID_`)



#model performance

# Calculate the null R2
# use glm for binary trait
# (will also need the fmsb package to calculate the pseudo R2)
null.model <- paste("COV", 1:24, sep = "", collapse = "+") %>%
    paste0("IC1DIM200~", .) %>%
    as.formula %>%
    lm(., data = pheno) %>%
    summary
null.r2 <- null.model$r.squared




reg.formula <- paste("COV", 1:24, sep = "", collapse = "+") %>%
    paste0("IC1DIM200~PRS+", .) %>%
    as.formula
reg.dat <- y
max.r2 <- 0
for(i in 1:ncol(pred_grid)){
    reg.dat$PRS <- pred_grid[,i]
    grid.model <- lm(reg.formula, dat=reg.dat) %>%
        summary  
    if(max.r2 < grid.model$r.squared){
        max.r2 <- grid.model$r.squared
    }
}
(result <- data.table(
    grid = max.r2 - null.r2,
    null = null.r2
))




#In the LDpred2 paper, we proposed an automatic way of filtering bad chains by comparing the scale of the resulting predictions. We have tested a somewhat equivalent and simpler alternative since, which we recommend here:

# `range` should be between 0 and 2
(range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
##  [1] 0.1210443 0.1216797 0.1209822 0.1197119 0.1187931 0.1199233 0.1202114
##  [8] 0.1195135 0.1214119 0.1206811 0.1189464 0.1204088 0.1196642 0.1195328
## [15] 0.1198751 0.1225441 0.1210127 0.1210234 0.1196245 0.1194824 0.1213715
## [22] 0.1188433 0.1203081 0.1196867 0.1210735 0.1201303 0.1209825 0.1195834
## [29] 0.1200190 0.1197259
(keep <- (range > (0.95 * quantile(range, 0.95))))
##  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
## [15] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
## [29] TRUE TRUE
To get the final effects / predictions, you should only use chains that pass this filtering:

beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
pred_auto <- big_prodVec(G, beta_auto, ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]])
pcor(pred_auto, y[ind.test], NULL)
## [1] 0.4967147 0.3669594 0.6075102

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
