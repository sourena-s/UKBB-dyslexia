


library(bigsnpr)
ukbb_genotypes <- snp_attach('/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/big40.rds')

for (thresh in c(0.001,0.01,0.1,100,10,15,1,20,25,30,35,40,45,50,55,5,60,65,70,75,80,85,90,95)) {
print(paste("analyzing thresh ", thresh))
#sum_stat_file='/data/clusterfs/lag/users/sousoh/ukbb/genetic/sbayesR/prs-AD-with-APOE.snpRes'
sum_stat_file=paste(sep="","/data/clusterfs/lag/users/sousoh/ukbb/genetic/sbayesR/prs-dyslexia.snpRes.stratified/top-",thresh,"-percentile")
my_sumstat<-read.table(sum_stat_file,header=T)
proc_sumstats <- my_sumstat[,c(3,4,6,5,8)]
colnames(proc_sumstats) <- c("chr", "pos", "a0", "a1" , "beta")
info_snp <- as.data.frame(ukbb_genotypes$map[,c(1,4,5,6)])
colnames(info_snp) <- c("chr", "pos", "a0" , "a1" )
info_snp$chr <- as.numeric(info_snp$chr)
results <- snp_match(proc_sumstats,info_snp,return_flip_and_rev=TRUE)

subject_PRS <- snp_PRS(ukbb_genotypes$genotypes, betas.keep=results$beta, ind.keep=results[,9])

m_path="/data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/prs-thresholds/dyslexia/"
m_filename=paste(sep="",m_path,thresh,".txt")

subs<-read.table('/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/list-all-subjects-genetics-subid-and-indices.txt',header=F)[,1]
write.table(cbind(subs,subject_PRS), m_filename,  row.names = FALSE, col.names=FALSE, quote=FALSE)

system(paste(sep="", "sort -k1,1 ", m_filename, paste(sep=""," >",m_filename,".sorted")))

}
