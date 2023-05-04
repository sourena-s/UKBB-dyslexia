library(bigsnpr)
library(data.table)
library(magrittr)

options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)
args = commandArgs(trailingOnly=TRUE)
m_chr=args[1]
print(paste("loading chromosome", m_chr))

#sumstats <- read.table("/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/dyslexia/dyslexia.filtered.ldpred2",header=T)
#names(sumstats) <-    c("chr","pos","rsid","a1","a0","n_eff","beta_se","p", "OR","INFO","MAF","beta")

#becareful with AD, bec sumstats are in GR38 coordinates
#umstats <- read.table("/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/AD/sumstats.ldpred2",header=T)
#ames(sumstats) <- c("rsid","chr","pos", "a1", "a0", "MAF", "beta", "beta_se", "p", "n_eff", "INFO")
# LDpred 2 require the header to follow the exact naming

#umstats$chr <- as.character(sumstats$chr)
#umstats <- sumstats[sumstats$chr!="X",]
#umstats$chr <- as.integer(sumstats$chr)

#only keep those variants with MAX sample size (for debugging purpose)
#umstats <- sumstats[sumstats$n_eff==max(sumstats$n_eff),]

# Filter out hapmap SNPs hapmap3: 
info <- readRDS("/data/clusterfs/lag/users/sousoh/ukbb/genetic/hapmap/hapmap3plus/map_hm3_plus.rds")
#top_variants<-read.table('/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/dyslexia/dyslexia.ldpred2.top-variants.pval',header=F,colClasses = "character")[,1]
#sumstats <- sumstats[sumstats$rsid %in% c(info$rsid),]

#umstats$chr <- as.integer(sumstats$chr)

#more filtering sumstats
#sumstats <- sumstats[  (sumstats$beta > mean(sumstats$beta) - 3* sd(sumstats$beta)) & (sumstats$beta < mean(sumstats$beta) + 3* sd(sumstats$beta)  ) ,]

#CORES <- 12
# Open a temporary file
tmp <- tempfile(tmpdir = paste(sep="","tmp-data-ld-processed")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
# Initialize variables for storing the LD score and LD matrix
corr <- NULL
ld <- NULL
# We want to know the ordering of samples in the bed file

obj.bigSNP <- snp_attach("/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/big40.rds")
fam.order<-read.table('/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/list-all-subjects-genetics-subid-and-indices.txt',header=F)[,1]
fam.order <- as.data.frame(cbind(fam.order,fam.order,fam.order,fam.order))
names(fam.order) <- c("family.ID", "sample.ID","FID", "IID")
obj.bigSNP$fam <- fam.order

info <- readRDS("/data/clusterfs/lag/users/sousoh/ukbb/genetic/hapmap/hapmap3plus/map_hm3_plus.rds")
hapmap_plus_ids <- which(obj.bigSNP$map$rsid %in% c(info$rsid))
random_2k <-read.table('/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/list-all-subjects-genetics-subid-and-indices-white-British-suff-2k.txt',header=F)[,1]
subj_2k_ids <- which(obj.bigSNP$fam$sample.ID %in% random_2k)
subset_file <- snp_subset(obj.bigSNP,ind.row=subj_2k_ids,ind.col=hapmap_plus_ids, backingfile="/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/processed/big40_subset_hapmap_2k")
obj.bigSNP.2k <- snp_attach("/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/processed/big40_subset_hapmap_2k.rds")


subset_file_all_40k <- snp_subset(obj.bigSNP,ind.col=hapmap_plus_ids, backingfile="/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/processed/big40_all_subj_subset_hapmap_2k")
#subset_file <- snp_subset(obj.bigSNP,ind.row=random_2k, backingfile="/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/big40_subset_2k")

#switch the below lines for (not) filtering based on hapmap plus variants
#obj.bigSNP.2k <- snp_attach("/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/big40_subset_2k.rds")

#this is old...................
obj.bigSNP.2k <- snp_attach("/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/big40_subset_2k_hapmap_plus.rds")

#how hapmap filtering was performed
#hapmap_plus_ids <- which(obj.bigSNP.2k$map$rsid %in% c(info$rsid))
#subset_file <- snp_subset(obj.bigSNP.2k,ind.col=hapmap_plus_ids, backingfile="/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/big40_subset_2k_hapmap_plus")


#fam.order<-read.table('/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/list-all-subjects-genetics-subid-and-indices.txt',header=F)[,1]
fam.order.2k <- read.table('/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/list-all-subjects-genetics-subid-and-indices-white-British-suff-2k.txt',header=F)[,1]
#fam.order <- as.data.frame(cbind(fam.order,fam.order,fam.order,fam.order))
fam.order.2k <- as.data.frame(cbind(fam.order.2k,fam.order.2k,fam.order.2k,fam.order.2k))
#names(fam.order) <- c("family.ID", "sample.ID","FID", "IID")
names(fam.order.2k) <- c("family.ID", "sample.ID","FID", "IID")
#obj.bigSNP$fam <- fam.order
obj.bigSNP.2k$fam <- fam.order.2k
#end old.................................



# extract the SNP information from the genotype
#map <- obj.bigSNP.2k$map[,c(1,3,4,5,6)]
#names(map) <- c("chr", "rsid", "pos", "a1", "a0")
#map$pos<- as.integer(map$pos)
#map$chr<- as.integer(map$chr)
# perform SNP matching
#sumstats<-sumstats[sample(1:dim(sumstats)[1],dim(sumstats)[1] %/% 50,replace=F),]
#info_snp <- snp_match(sumstats, join_by_pos = FALSE, map, return_flip_and_rev=TRUE, strand_flip =FALSE)

# better to use af of GWAS and INFO scores as well (then can use 0.7 instead of 0.5)
#sd_ldref <- with(info_snp, sqrt(2 * af_UKBB * (1 - af_UKBB)))
#sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2))
#is_bad <-
#  sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.05 | sd_ldref < 0.05

#subset_file <- snp_subset(obj.bigSNP.2k,ind.col=info_snp$'_NUM_ID_')
# genotype <- snp_attach('/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/big40_subset_verysmall.rds')#genotype <- snp_attach(subset_file)

#genotype <- snp_attach('/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/big40_subset_2k_sub1.rds')
#genotype <- snp_attach("/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/big40_subset_small.rds")
genotype <- obj.bigSNP.2k
POS2 <- snp_asGeneticPos(as.integer(genotype$map$chromosome), as.integer(genotype$map$physical.pos), dir = "/data/clusterfs/lag/users/sousoh/ukbb/genetic/hapmap" ,type="hapmap" )

#uncomment for batch submnission of chromosomes
#for (chr in m_chr) {
#    ind.chr <- which(as.integer(genotype$map$chromosome) == chr)
#print (paste("calculating LD for ", length(ind.chr), "variants in chromosome", chr))
#    corr0 <- snp_cor(genotype$genotype, ind.col = ind.chr, ncores = 1, infos.pos = POS2[ind.chr], size= 3/1000)
#    	ld <- Matrix::colSums(corr0^2)
#        corr <- as_SFBM(corr0, tmp)
#    print (paste("Chromosome", chr,"processed."))
#}

#saveRDS(ld, file = paste(sep="","/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/matrix_big40_ld_chr_",m_chr,".rds")) 
#saveRDS(corr, file = paste(sep="","/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/matrix_big40_corr_chr_",m_chr,".rds")) 

for (chr in 1:22) {
    ind.chr <- which(as.integer(genotype$map$chromosome) == chr)
print (paste("calculating LD for ", length(ind.chr), "variants in chromosome", chr))
    corr0 <- snp_cor(genotype$genotype, ind.col = ind.chr, ncores = 64, infos.pos = POS2[ind.chr], size= 3/1000)

if (chr == 1) {
    	ld <- Matrix::colSums(corr0^2)
        corr <- as_SFBM(corr0, tmp)
    } else {
        ld <- c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
    }
    print (paste("Chromosome", chr,"processed."))
}

print ("done calculating, now saving..")

saveRDS(ld, file = "/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/processed/matrix_big40_ld_hapmap.rds") 
saveRDS(corr, file = "/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/processed/matrix_big40_corr_hapmap.rds") 
print ("saved LD matrices.")
