library(bigsnpr)
library(data.table)
library(magrittr)

options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

args = commandArgs(trailingOnly=TRUE)
i=as.integer(args[1])
node_cores=as.integer(args[2])
node_cores=16
#node_cores=nb_cores() - 4
print(paste("Running Lassosum2 using", node_cores,"cores.."))

all_sumstats <- vector()
all_sumstats[1] <-"/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/dyslexia/dyslexia.filtered.ldpred2"
#becareful with AD, bec sumstats are in GR38 coordinates
all_sumstats[2] <-"/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/AD-jansen-2019/AD_sumstats_Jansenetal_2019sept_no_apoe_1MB.txt"
all_sumstats[3] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/total-surf/0648.txt"
all_sumstats[4] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/putamen/0015.txt"
all_sumstats[5] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/SCZ/SCZ_EURO_2022.tsv"
all_sumstats[6] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/ADHD/ADHD.proc"
all_sumstats[7] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/ASD/ASD.proc"
all_sumstats[8] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/grade-E1/schoolgrades.E1.sumstats"
all_sumstats[9] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/grade-E2/schoolgrades.E2.sumstats"
all_sumstats[10] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/grade-E3/schoolgrades.E3.sumstats"
all_sumstats[11] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/grade-E4/schoolgrades.E4.sumstats"
all_sumstats[12] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/migraine/migraine_ihgc2021_gws_gwama_0.corrected"
all_sumstats[13] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/PD/sumstats-with-rsid.txt"
all_sumstats[14] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/BIP/bip.txt"
all_sumstats[15] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/UKBB-GCSE/categorical-6138-both_sexes-3.tsv"
all_sumstats[16] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/UKBB-NVQ/categorical-6138-both_sexes-5.tsv"
all_sumstats[17] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/UKBB-PAIN/categorical-6159-both_sexes-8.tsv"
all_sumstats[18] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/UKBB-REACTION-TIME/continuous-20023-both_sexes-irnt.tsv"
all_sumstats[19] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/UKBB-VNR/continuous-20016-both_sexes-irnt.tsv"
all_sumstats[20] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/UKBB-VNR-comp1/sumstat.tsv"
all_sumstats[21] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/UKBB-VNR-comp2/sumstat.tsv"
all_sumstats[22] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/UKBB-VNR-comp3/sumstat.tsv"
all_sumstats[23] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/UKBB-VNR-comp4/sumstat.tsv"
all_sumstats[24] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/UKBB-VNR-comp5/sumstat.tsv"
all_sumstats[25] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/UKBB-VNR-comp6/sumstat.tsv"
all_sumstats[26] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/UKBB-VNR-comp7/sumstat.tsv"
all_sumstats[27] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/UKBB-VNR-comp8/sumstat.tsv"
all_sumstats[28] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/UKBB-VNR-comp9/sumstat.tsv"
all_sumstats[29] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/UKBB-VNR-comp10/sumstat.tsv"
all_sumstats[30] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/UKBB-VNR-comp11/sumstat.tsv"
all_sumstats[31] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/UKBB-VNR-comp12/sumstat.tsv"
all_sumstats[32] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/UKBB-VNR-comp13/sumstat.tsv"
all_sumstats[33] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/UKBB-hypertension/phecode-401-both_sexes.tsv"
all_sumstats[34] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/UKBB-risk-taking/categorical-2040-both_sexes-2040.tsv"
all_sumstats[35] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/AD-kunkle/Kunkle_etal_Stage1_results_no_apoE_1MB_flank.txt"
all_sumstats[36] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/handedness/Ambidextrous_UKBB.txt"
all_sumstats[37] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/handedness-L/LeftHandedness_MetaAnalysis_UKBB_IHC.txt"
all_sumstats[38] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/else-WR/proc-sum"
all_sumstats[39] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/else-SP/proc-sum"
all_sumstats[40] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/else-PA/proc-sum"
all_sumstats[41] <- "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/else-NREAD/proc-sum"


cond <-vector()
cond[1] <- "dyslexia"
cond[2] <- "AD-jansen-2019"
cond[3] <- "total-surf"
cond[4] <- "putamen"
cond[5] <- "SCZ"
cond[6] <- "ADHD"
cond[7] <- "ASD"
cond[8] <- "grade-E1"
cond[9] <- "grade-E2"
cond[10] <- "grade-E3"
cond[11] <- "grade-E4"
cond[12] <- "migraine"
cond[13] <- "PD"
cond[14] <- "BIP"
cond[15] <- "UKBB-GCSE"
cond[16] <- "UKBB-NVQ"
cond[17] <- "UKBB-PAIN"
cond[18] <- "UKBB-REACTION-TIME"
cond[19] <- "UKBB-VNR"
cond[20] <- "UKBB-VNR-comp1"
cond[21] <- "UKBB-VNR-comp2"
cond[22] <- "UKBB-VNR-comp3"
cond[23] <- "UKBB-VNR-comp4"
cond[24] <- "UKBB-VNR-comp5"
cond[25] <- "UKBB-VNR-comp6"
cond[26] <- "UKBB-VNR-comp7"
cond[27] <- "UKBB-VNR-comp8"
cond[28] <- "UKBB-VNR-comp9"
cond[29] <- "UKBB-VNR-comp10"
cond[30] <- "UKBB-VNR-comp11"
cond[31] <- "UKBB-VNR-comp12"
cond[32] <- "UKBB-VNR-comp13"
cond[33] <- "UKBB-hypertension"
cond[34] <- "UKBB-risk-taking"
cond[35] <- "AD-kunkle"
cond[36] <- "handedness"
cond[37] <- "handedness-L"
cond[38] <- "else-WR"
cond[39] <- "else-SP"
cond[40] <- "else-PA"
cond[41] <- "else-NREAD"

cond_name <- cond[i]

ld<-readRDS('/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/processed/matrix_big40_ld_hapmap.rds')
corr<-readRDS("/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/processed/matrix_orig_corr_hapmap.rds")
hapmap <- readRDS("/data/clusterfs/lag/users/sousoh/ukbb/genetic/hapmap/hapmap3plus/map_hm3_plus.rds")
ld <- hapmap$ld

obj.bigSNP.all <- snp_attach("/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/processed/big40_all_subj_subset_hapmap_2k.rds")
map <- obj.bigSNP.all$map[,c(1,3,4,5,6)]
names(map) <- c("chr", "rsid", "pos", "a1", "a0")
map$pos<- as.integer(map$pos)
map$chr<- as.integer(map$chr)
map$beta <- 1
fam.order<-read.table('/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/list-all-subjects-genetics-subid-and-indices.txt',header=F)[,1]
fam.order <- as.data.frame(cbind(fam.order,fam.order,fam.order,fam.order))
names(fam.order) <- c("family.ID", "sample.ID","FID", "IID")
obj.bigSNP.all$fam <- fam.order

sumstats <- readRDS(paste(sep=".",all_sumstats[i],"rds") )
     #leftover from a small bug, the below two lines can now be removed
     sumstats$a0 <- toupper(sumstats$a0)
     sumstats$a1 <- toupper(sumstats$a1)

cond_name <- cond[i]
#single overlap: sumstat in the space of hapmap
info_snp <- snp_match(sumstats, join_by_pos = FALSE, hapmap, return_flip_and_rev=TRUE, strand_flip =FALSE,  match.min.prop=0.05)
#double overlap: sumstat and big40 in the space of hapmap
info_snp_2 <- snp_match(map, join_by_pos = FALSE, hapmap[info_snp$"_NUM_ID_",], return_flip_and_rev=TRUE, strand_flip =FALSE)
#extract the overlapping variant indices in the hapmap space
info_snp_triple_overlap <- snp_match(info_snp_2[,c("chr","rsid","a0","a1","beta")], hapmap,join_by_pos = FALSE,  return_flip_and_rev=TRUE, strand_flip =FALSE)
overlap_vars <- info_snp_triple_overlap$"_NUM_ID_"
#get allele-corrected GWAS beta estimate dataframe
info_snp_3 <- snp_match(sumstats, join_by_pos = FALSE, hapmap[overlap_vars,], return_flip_and_rev=TRUE, strand_flip =FALSE)
###find var indices on big40.all for the prediction step: note: all should match
info_snp_big40 <- snp_match(info_snp_3[,c("chr","rsid","a0","a1","beta")] , join_by_pos = FALSE, map, return_flip_and_rev=TRUE, strand_flip =FALSE)

df_beta <- info_snp_3[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]

beta_lassosum2 <- snp_lassosum2(corr, df_beta, ncores = node_cores , ind.corr = overlap_vars ,
delta=seq_log(0.01,10000,10),nlambda=30,lambda.min.ratio=0.01, maxiter=1000,tol=1e-2)

pred_lasso <- big_prodMat(obj.bigSNP.all$genotype, beta_lassosum2, ind.col = info_snp_big40$"_NUM_ID_",  ncores = node_cores, block.size=5000 )

dir.create(file.path("/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs", cond_name), showWarnings = FALSE)
dir.create(file.path("/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs", cond_name, "lasso" ) , showWarnings = FALSE)

write.table(beta_lassosum2,paste(sep="","/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/",cond_name,"/lasso/beta_lassosum2-hapmap.txt"), quote=F,row.names=F, col.names=F,sep=",")
write.table(pred_lasso,paste(sep="","/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/",cond_name, "/lasso/lasso-hapmap.pgs"), quote=F,row.names=F, col.names=F,sep=",")
write.table(attr(beta_lassosum2, "grid_param") ,paste(sep="","/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/",cond_name,"/lasso/grid-hapmap.params"), quote=F,row.names=F,sep=",")
write.table(info_snp_triple_overlap[,1:4],paste(sep="","/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/",cond_name, "/lasso/overlap-variants-in-pgs-beta.txt"), quote=F,row.names=F, col.names=F, sep=" ")

