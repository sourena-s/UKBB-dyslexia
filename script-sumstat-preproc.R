library(bigsnpr)
library(data.table)
library(magrittr)

options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

args = commandArgs(trailingOnly=TRUE)
i=as.integer(args[1])


all_sumstats <- vector()
all_sumstats[1] <-"/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/dyslexia/dyslexia.filtered.ldpred2"
#becareful with AD, bec sumstats are in GR38 coordinates
#all_sumstats[2] <-"/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/AD/sumstats.ldpred2"
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
all_sumstat_colnames <- list()


colnames_ukbb_categorical <- c("chr","pos","a0","a1","af_cases_meta_hq","af_controls_meta_hq","beta_meta_hq","se_meta_hq","neglog10_pval_meta_hq","neglog10_pval_heterogeneity_hq","af_cases_meta","af_controls_meta","beta_meta","se_meta","neglog10_pval_meta","neglog10_pval_heterogeneity","af_cases_AFR","af_cases_AMR","af_cases_CSA","af_cases_EAS","af_cases","af_cases_MID","af_controls_AFR","af_controls_AMR","af_controls_CSA","af_controls_EAS","af_controls","af_controls_MID","beta_AFR","beta_AMR","beta_CSA","beta_EAS","beta","beta_MID","se_AFR","se_AMR","se_CSA","se_EAS","beta_se","se_MID","neglog10_pval_AFR","neglog10_pval_AMR","neglog10_pval_CSA","neglog10_pval_EAS","neglog10_pval","neglog10_pval_MID","low_confidence_AFR","low_confidence_AMR","low_confidence_CSA","low_confidence_EAS","low_confidence","low_confidence_MID")

colnames_ukbb_continuous <- c("chr","pos","a0","a1","af_meta_hq","beta_meta_hq","se_meta_hq","neglog10_pval_meta_hq","neglog10_pval_heterogeneity_hq","af_meta","beta_meta","se_meta","neglog10_pval_meta","neglog10_pval_heterogeneity","af_AFR","af_AMR","af_CSA","af_EAS","MAF","beta_AFR","beta_AMR","beta_CSA","beta_EAS","beta","se_AFR","se_AMR","se_CSA","se_EAS","beta_se","neglog10_pval_AFR","neglog10_pval_AMR","neglog10_pval_CSA","neglog10_pval_EAS","neglog10_pval","low_confidence_AFR","low_confidence_AMR","low_confidence_CSA","low_confidence_EAS","low_confidence")

all_sumstat_colnames[[1]] <- c("chr","pos","rsid","a1","a0","n_eff","beta_se","p", "OR","INFO","MAF","beta")
#all_sumstat_colnames[[2]] <- c("rsid","chr","pos", "a1", "a0", "MAF", "beta", "beta_se", "p", "n_eff", "INFO")
all_sumstat_colnames[[2]] <- c("uniqID.a1a2","chr","pos","a1","a0","rsid","z","p","Nsum","n_eff","dir","EAF","beta","beta_se")
#in BIG40 imaging genetic GWAS, the second allele is effect (a2)
all_sumstat_colnames[[3]] <- c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "minuslogp")
all_sumstat_colnames[[4]] <- all_sumstat_colnames[[3]]
all_sumstat_colnames[[5]] <- c("chr","rsid","pos","a1","a0","FCAS","FCON","IMPINFO","beta","beta_se","p","NCAS","NCON","n_eff")
all_sumstat_colnames[[6]] <- c("n_eff" ,"freq" ,"beta", "chr" ,"rsid" ,"pos" ,"a1" ,"a0" ,"FRQ_A_38691" ,"FRQ_U_186843" ,"info" ,"OR" ,"beta_se" ,"p" ,"Direction" ,"Nca" ,"Nco")
all_sumstat_colnames[[7]] <- c("beta", "chr", "rsid", "pos", "a1", "a0", "info", "OR", "beta_se", "p")
all_sumstat_colnames[[8]] <- c("rsid", "chr", "pos", "a1", "a0", "beta", "beta_se", "p", "info", "n_eff")
all_sumstat_colnames[[9]] <- all_sumstat_colnames[[8]] 
all_sumstat_colnames[[10]] <- all_sumstat_colnames[[8]] 
all_sumstat_colnames[[11]] <- all_sumstat_colnames[[8]] 
all_sumstat_colnames[[12]] <- c("rsid","chr","pos","a1","a0","freq","beta","beta_se","beta_95L","beta_95U","z","p","minus_log10_p","q_statistic","q_p_value","i2","n_studies","n_samples","effects","n_eff")
all_sumstat_colnames[[13]] <- c("chr","pos","a0","a1","freq","beta","beta_se","p","N_cases","N_controls","n_eff","rsid")
all_sumstat_colnames[[14]] <- c("chr","pos","rsid","a1","a0","beta","beta_se","p","NGT","FCAS","FCON","info","n_eff_div_2","NCAS","NCON","DIRE")

#<old> won't be used (see the if statement below)
all_sumstat_colnames[[15]] <- colnames_ukbb_categorical
all_sumstat_colnames[[16]] <- colnames_ukbb_categorical
all_sumstat_colnames[[17]] <- c("chr","pos","a0","a1","af_cases_meta_hq","af_controls_meta_hq","beta_meta_hq","se_meta_hq","neglog10_pval_meta_hq","neglog10_pval_heterogeneity_hq","af_cases_meta","af_controls_meta","beta_meta","se_meta","neglog10_pval_meta","neglog10_pval_heterogeneity","af_cases_AFR","af_cases_CSA","af_cases_EUR","af_cases_MID","af_controls_AFR","af_controls_CSA","af_controls_EUR","af_controls_MID","beta_AFR","beta_CSA","beta","beta_MID","se_AFR","se_CSA","beta_se","se_MID","neglog10_pval_AFR","neglog10_pval_CSA","neglog10_pval","neglog10_pval_MID","low_confidence_AFR","low_confidence_CSA","low_confidence","low_confidence_MID")

all_sumstat_colnames[[18]] <- c("chr","pos","a0","a1","af_meta_hq","beta_meta_hq","se_meta_hq","neglog10_pval_meta_hq","neglog10_pval_heterogeneity_hq","af_meta","beta_meta","se_meta","neglog10_pval_meta","neglog10_pval_heterogeneity","af_AFR","af_AMR","af_CSA","af_EAS","af_EUR","af_MID","beta_AFR","beta_AMR","beta_CSA","beta_EAS","beta","beta_MID","se_AFR","se_AMR","se_CSA","se_EAS","beta_se","se_MID","neglog10_pval_AFR","neglog10_pval_AMR","neglog10_pval_CSA","neglog10_pval_EAS","neglog10_pval","neglog10_pval_MID","low_confidence_AFR","low_confidence_AMR","low_confidence_CSA","low_confidence_EAS","low_confidence","low_confidence_MID")

all_sumstat_colnames[[19]] <- c("chr","pos","a0","a1","af_meta_hq","beta_meta_hq","se_meta_hq","neglog10_pval_meta_hq","neglog10_pval_heterogeneity_hq","af_meta","beta_meta","se_meta","neglog10_pval_meta","neglog10_pval_heterogeneity","af_AFR","af_AMR","af_CSA","af_EAS","af_EUR","beta_AFR","beta_AMR","beta_CSA","beta_EAS","beta","se_AFR","se_AMR","se_CSA","se_EAS","beta_se","neglog10_pval_AFR","neglog10_pval_AMR","neglog10_pval_CSA","neglog10_pval_EAS","neglog10_pval","low_confidence_AFR","low_confidence_AMR","low_confidence_CSA","low_confidence_EAS","low_confidence")

all_sumstat_colnames[[20]] <- c("chr","pos","a0","a1","af_cases_meta","af_controls_meta","beta_meta","se_meta","neglog10_pval_meta","neglog10_pval_heterogeneity","af_cases_AFR","af_cases_AMR","af_cases_CSA","af_cases_EAS","af_cases_EUR","af_cases_MID","af_controls_AFR","af_controls_AMR","af_controls_CSA","af_controls_EAS","af_controls_EUR","af_controls_MID","beta_AFR","beta_AMR","beta_CSA","beta_EAS","beta","beta_MID","se_AFR","se_AMR","se_CSA","se_EAS","beta_se","se_MID","neglog10_pval_AFR","neglog10_pval_AMR","neglog10_pval_CSA","neglog10_pval_EAS","neglog10_pval","neglog10_pval_MID","low_confidence_AFR","low_confidence_AMR","low_confidence_CSA","low_confidence_EAS","low_confidence","low_confidence_MID")
all_sumstat_colnames[[21]] <- all_sumstat_colnames[[20]] 
all_sumstat_colnames[[22]] <- c("chr","pos","a1","a0","af_cases_meta_hq","af_controls_meta_hq","beta_meta_hq","se_meta_hq","neglog10_pval_meta_hq","neglog10_pval_heterogeneity_hq","af_cases_meta","af_controls_meta","beta_meta","se_meta","neglog10_pval_meta","neglog10_pval_heterogeneity","af_cases_AFR","af_cases_AMR","af_cases_CSA","af_cases_EAS","af_cases_EUR","af_cases_MID","af_controls_AFR","af_controls_AMR","af_controls_CSA","af_controls_EAS","af_controls_EUR","af_controls_MID","beta_AFR","beta_AMR","beta_CSA","beta_EAS","beta","beta_MID","se_AFR","se_AMR","se_CSA","se_EAS","beta_se","se_MID","neglog10_pval_AFR","neglog10_pval_AMR","neglog10_pval_CSA","neglog10_pval_EAS","neglog10_pval","neglog10_pval_MID","low_confidence_AFR","low_confidence_AMR","low_confidence_CSA","low_confidence_EAS","low_confidence","low_confidence_MID") 
all_sumstat_colnames[[23]] <- all_sumstat_colnames[[20]] 
all_sumstat_colnames[[24]] <- all_sumstat_colnames[[20]] 
all_sumstat_colnames[[25]] <- all_sumstat_colnames[[20]] 
all_sumstat_colnames[[26]] <- colnames_ukbb_categorical
all_sumstat_colnames[[27]] <- colnames_ukbb_categorical
all_sumstat_colnames[[28]] <- c("chr","pos","a0","a1","af_cases_meta_hq","af_controls_meta_hq","beta_meta_hq","se_meta_hq","neglog10_pval_meta_hq","neglog10_pval_heterogeneity_hq","af_cases_meta","af_controls_meta","beta_meta","se_meta","neglog10_pval_meta","neglog10_pval_heterogeneity","af_cases_AFR","af_cases_CSA","af_cases_EUR","af_controls_AFR","af_controls_CSA","af_controls_EUR","beta_AFR","beta_CSA","beta","se_AFR","se_CSA","beta_se","neglog10_pval_AFR","neglog10_pval_CSA","neglog10_pval","low_confidence_AFR","low_confidence_CSA","low_confidence")
all_sumstat_colnames[[29]] <- c("chr","pos","a0","a1","af_cases_meta","af_controls_meta","beta_meta","se_meta","neglog10_pval_meta","neglog10_pval_heterogeneity","af_cases_CSA","af_cases_EUR","af_controls_CSA","af_controls_EUR","beta_CSA","beta","se_CSA","beta_se","neglog10_pval_CSA","neglog10_pval","low_confidence_CSA","low_confidence")
all_sumstat_colnames[[30]] <- all_sumstat_colnames[[29]] 
all_sumstat_colnames[[31]] <- c("chr","pos","a0","a1","af_cases_EUR","af_controls_EUR","beta","beta_se","neglog10_pval","low_confidence")
all_sumstat_colnames[[32]] <- all_sumstat_colnames[[31]] 
all_sumstat_colnames[[35]] <- c("chr","pos","rsid","a1","a0","beta","beta_se","p")
all_sumstat_colnames[[36]] <- c("rsid","chr","pos","GENPOS","a1","a0","freq","info","CHISQ_LINREG","P_LINREG","beta","beta_se","CHISQ_BOLT_LMM_INF","p","n","n_eff")
all_sumstat_colnames[[37]] <- c("rsid","a1","a0","freq","FreqSE","MinFreq","MaxFreq","n_eff","z","p","Direction","HetISq","HetChiSq","HetDf","HetPVal","chr","pos")

#The standardized beta (i.e assuming both Y and X are transformed to have unit variance and mean zero) = Zscore*sqrt(Var(Y|X)/N)

#Var(Y|X) = 1/(1 + (Zscore*Zscore)/N)
#Var(beta) = Var(Y|X)/N

#Beta = z / sqrt(2p(1− p)(n + z^2)) and
#SE =1 / sqrt(2p(1− p)(n + z^2))
#Where p is the frequency of the imputed SNP, 

##<end old/>
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
#dummy beta
map$beta <- 1
fam.order<-read.table('/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/list-all-subjects-genetics-subid-and-indices.txt',header=F)[,1]
fam.order <- as.data.frame(cbind(fam.order,fam.order,fam.order,fam.order))
names(fam.order) <- c("family.ID", "sample.ID","FID", "IID")
obj.bigSNP.all$fam <- fam.order

#if (file.exists (paste(sep=".",all_sumstats[i],"rds")) ) {
#	sumstats <- readRDS( paste(sep=".",all_sumstats[i],"rds")  )

#}else {
print(paste("Reading summary stat file:",all_sumstats[i]))
	sumstats <- read.table(all_sumstats[i],header=T)

	#check if UK Biobabnk header
	if(i>=15 & i<=34) {
		header <- colnames(sumstats)
		header[which(header=="alt")] <- "a1"
		header[which(header=="ref")] <- "a0"
		header[which(header=="beta_EUR")] <- "beta"
		header[which(header=="se_EUR")] <- "beta_se"
		header[which(header=="neglog10_pval_EUR")] <- "neglog10_pval"
		header[which(header=="low_confidence_EUR")] <- "low_confidence"
		names(sumstats) <- header

	}else{

		names(sumstats) <- all_sumstat_colnames[[i]]
	}



	if (i==3 | i==4) {sumstats$n_eff <- 33000
		sumstats$p <- 10 ^ (-1 * as.numeric(sumstats$minuslogp))
	}
	if (i==7) {sumstats$n_eff <-  46350 }
	if (i==14) {sumstats$n_eff <-  sumstats$n_eff_div_2 * 2  }
	
	if (i==35) {sumstats$n_eff <-  57693}

	if (i>=15 & i<=34){ sumstats$n_eff <- c(
				       195589 + 220727, #GCSE
				       78286  + 338030, #NVQ
				       6603   + 412985, #pain
				       417660, #reaction time
				       rep(135088,14), #VNR (fluid intelligence) and 13 subcomponents				       
				       94311  + 325488, #hypertension
				       104809 + 301022 #risk taking
				       ) [i-14]
		sumstats$p <- 10 ^ (-1 * sumstats$neglog10_pval)
		variants<-read.table('/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/UKBB-VNR/full_variant_qc_metrics.txt',header=T,sep="\t")
		sumstats$rsid<- variants$rsid
		sumstats$info <- variants$info
		sumstats$freq <- variants$af_EUR
		sumstats <- sumstats[ ! sumstats$low_confidence == "true" & !is.na(sumstats$beta) ,]
		}	
	if (i==37) {
	sumstats$z <- as.numeric(sumstats$z)
	sumstats$n_eff <- as.numeric (sumstats$n_eff)
	sumstats$var_y_x = 1/(1 + (sumstats$z ^ 2)/ sumstats$n_eff)
	sumstats$beta <- sumstats$z * sqrt(sumstats$var_y_x / sumstats$n_eff)
	sumstats$beta_se = sumstats$var_y_x  / sumstats$n_eff
	sumstats$a0 <- toupper(sumstats$a0)
	sumstats$a1 <- toupper(sumstats$a1)
	}

	sumstats$chr <- as.character(sumstats$chr)
	sumstats <- sumstats[sumstats$chr!="X",]
	sumstats$chr <- as.integer(sumstats$chr)
	#sumstats <- sumstats[sumstats$n_eff==max(sumstats$n_eff),]
	#sumstats <- sumstats %>% drop_na()
	saveRDS(sumstats,file=paste(sep=".",all_sumstats[i],"rds") )
#}

#single overlap: sumstat in the space of hapmap
info_snp <- snp_match(sumstats, join_by_pos = F, hapmap, return_flip_and_rev=T, strand_flip =F,  match.min.prop=0.05)
#double overlap: sumstat and big40 in the space of hapmap
info_snp_2 <- snp_match(map, join_by_pos = F, hapmap[info_snp$"_NUM_ID_",], return_flip_and_rev=T, strand_flip =F)
#extract the overlapping variant indices in the hapmap space
info_snp_triple_overlap <- snp_match(info_snp_2[,c("chr","rsid","a0","a1","beta")], hapmap,join_by_pos = F,  return_flip_and_rev=T, strand_flip =F)
overlap_vars <- info_snp_triple_overlap$"_NUM_ID_"
#get allele-corrected GWAS beta estimate dataframe
info_snp_3 <- snp_match(sumstats, join_by_pos = F, hapmap[overlap_vars,], return_flip_and_rev=T, strand_flip =F)
###find var indices on big40.all for the prediction step: note: all should match
info_snp_big40 <- snp_match(info_snp_3[,c("chr","rsid","a0","a1","beta")] , join_by_pos = F, map, return_flip_and_rev=T, strand_flip =F)

sumstat_ma <- paste(sep="","/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/",cond_name,"/auto-sumstats.ma")
#if (! file.exists(sumstat_ma) ) {
sum_export <- info_snp_3[,c("rsid","a1","a0","af_UKBB","beta","beta_se","p","n_eff")]
names(sum_export)<-c("SNP", "A1", "A2", "freq" , "b", "se", "p", "N")
write.table(sum_export,file=sumstat_ma,quote=F,row.names=F,col.names=T,sep=" ") 
print(paste("Writing .ma corrected sumstat file for",cond_name)) 
#}

df_beta <- info_snp_3[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]

ldsc <- snp_ldsc(ld[overlap_vars], length(ld[overlap_vars]), chi2 = (df_beta$beta / df_beta$beta_se)^2, sample_size = df_beta$n_eff, blocks = NULL)
h2_est <- ldsc[["h2"]]
print(paste("h2_est is",h2_est))

write.table(h2_est,file=paste(sep=".",sumstat_ma,"h2"),quote=F,row.names=F,col.names=T,sep=" ")


