

library(bigsnpr)
ukbb_genotypes <- snp_attach('/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/big40.rds')


a<-read.table('/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/british-ancestry-subjects',header=F)
a<-as.data.frame(cbind(a,0,0,0,-9))
colnames(a) <- c("family.ID","sample.ID","paternal.ID","maternal.ID","sex","affection")

ukbb_genotypes$fam <- a

bed_file <- tempfile(fileext = ".bed")
bed <- snp_writeBed(ukbb_genotypes, bed_file)
