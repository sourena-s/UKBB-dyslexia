library(bigsnpr)
convert_flag=F
if (convert_flag) {

#bgens=vector(mode="character",length=22)
bgens=vector(mode="character",length=1)
bgen_snps=list()
#my_backing="/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/big40"
my_backing="/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/big40_chrX"

subjects=read.table('/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/list-all-subjects-genetics-indices.txt',header=F,colClasses="numeric")[,1]
subjects_x= read.table('list-all-subjects-genetics-subid-and-indices-chrx-and-autosome',header=F,colClasses="numeric")[,2]

#for (i in 1:22) { bgens[i]=paste(sep="","/data/workspaces/lag/workspaces/lg-ukbiobank/primary_data/genetic_data/snp/snp_release_v3/data/imp/ukb_imp_chr", i,"_v3.bgen")}
#for (i in 1:22) { bgen_snps[[i]]=read.table(paste(sep="","/data/clusterfs/lag/users/sousoh/ukbb/genetic/snp-qc/snp-post-qc-chr-", i), header=F,colClasses="character")[,1]}

i="X"
bgens[1]=paste(sep="","/data/workspaces/lag/workspaces/lg-ukbiobank/primary_data/genetic_data/snp/snp_release_v3/data/imp/ukb_imp_chr", i,"_v3.bgen")
bgen_snps[[1]]=read.table(paste(sep="","/data/clusterfs/lag/users/sousoh/ukbb/genetic/snp-qc/snp-post-qc-chr-", i), header=F,colClasses="character")[,1]


snp_readBGEN(bgens, ind_row=subjects_x, backingfile=my_backing,  read_as = "dosage", list_snp_id=bgen_snps, ncores=16 )
snp_readBGEN(bgens, ind_row=subjects, backingfile=my_backing,  read_as = "dosage", list_snp_id=bgen_snps, ncores=16 )


}
ukbb_genotypes <- snp_attach('/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/big40.rds')
#ukbb_genotypes <- snp_attach('/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/big40_chrX.rds')


#or sum_stat_file='/data/clusterfs/lag/users/sousoh/ukbb/genetic/sbayesR/prs-dyslexia.snpRes.stratified/top-100-percentile'
#or sum_stat_file='/data/clusterfs/lag/users/sousoh/ukbb/genetic/sbayesR/prs-AD-with-APOE.snpRes'
#raw dyslexia sumstats for C+P method:
#sum_stat_file='/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/dyslexia-top10k/2nd-clump.snpRes'
#sum_stat_file='/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/dyslexia/dyslexia.filtered.2.dat.overlap.big40.with.indel.alleles.snpRes'
sum_stat_file='/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/dyslexia/dyslexia.filtered.ldpred2'
#sum_stat_file='/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/big40.map.TUBB4B.snpRes'
#sum_stat_file='/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/dyslexia/dyslexia.filtered.2.dat.overlap.big40.with.indel.alleles.chrX.snpRes'
#sum_stat_file='/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/AD/top-ad-migraine.snpRes'
my_sumstat<-read.table(sum_stat_file,header=T)

#In the SNP weight sumstat file (SBayesR .snpRes output file), allele codings are A1, A2 and A1-eff,
#i.e., the effect (beta) is wrt to A1 which is the first of the two allele columns
#in contrast, bigSNPr expects a0 and a1 coding, and the effect allele corresponds to the second column (a1)
#hence the column orders have been reversed below: c(6,5) to match SbayesR and bigSNPr effect allele conventions
proc_sumstats <- my_sumstat[,c(1,2,5,4,12)]
colnames(proc_sumstats) <- c("chr", "pos", "a0", "a1" , "beta")
info_snp <- as.data.frame(ukbb_genotypes$map[,c(1,4,5,6)])
colnames(info_snp) <- c("chr", "pos", "a0" , "a1" )
info_snp$chr <- as.numeric(info_snp$chr)
info_snp$pos <- as.numeric(info_snp$pos)

proc_sumstats$chr <- as.numeric(proc_sumstats$chr)
proc_sumstats$pos <- as.numeric(proc_sumstats$pos)

#disabling strand flips so we won't lose too many variants being marked as ambiguous
#only 1 variant was flipped in the previous run
results <- snp_match(proc_sumstats,info_snp,return_flip_and_rev=TRUE, strand_flip =FALSE)

#subject_PRS <- snp_PRS(ukbb_genotypes$genotypes, betas.keep=results$beta, ind.keep=results[,9])

subs<-read.table('/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/list-all-subjects-genetics-subid-and-indices.txt',header=F)[,1]
#subs<-read.table('/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/list-all-subjects-genetics-subid-and-indices-chrx-and-autosome',header=F)[,1]
#write.table(cbind(subs,subject_PRS), "/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/dyslexia-top10k/2nd-clump-PRS.txt",  row.names = FALSE, col.names=FALSE, quote=FALSE)

b<- as.data.frame(ukbb_genotypes$map)
colnames(b)<- c("chr","marker_id","rsid","pos","a1","a2","freq","info")
b$chr <- as.numeric(b$chr)
#b$chr <- 23
b$pos <- as.numeric(b$pos)
merged_data <- merge(b, results,by=c("chr","pos"))
#snp_list=merged_data$rsid

#snp_list=read.table('/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/dyslexia-top10k/clumped-variants.txt', colClasses="character", header=F)[,1]
#snp_list=read.table('/data/clusterfs/lag/users/sousoh/ukbb/genetic/1000GP/test', header=F)[,1]
#snp_list=read.table('/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/dyslexia/var-list-clumped-p-1percent',header=F)[,1]
#snp_list=read.table('/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs/dyslexia/var-list-clumped-p-1percent.chrX',header=T)[,1]
snp_list=read.table('/data/clusterfs/lag/users/sousoh/ukbb/genetic/UKBB-clump/var-clump-1e-3.txt',header=F)[,1]

for (snp in snp_list){

snp_id=which(merged_data$rsid==snp)

if ((length(snp_id)==1)  && !(file.exists(paste(sep="","/data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/variant-dosage/",snp,".txt.sorted")))) {

snp_risk_allele_dosage <- ukbb_genotypes$genotypes[,merged_data[snp_id,15]] * ifelse(merged_data[snp_id,]$beta < 0, 1,-1) +  ifelse(merged_data[snp_id,]$beta < 0, 0, 2)
snp_data <- as.data.frame(cbind(ukbb_genotypes$map[merged_data[snp_id,15],],merged_data[snp_id,]$beta,merged_data[snp_id,14]))
colnames(snp_data)[9] <- "corrected_gwas_beta_allele2"
colnames(snp_data)[10] <- "reverse"

m_path="/data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/variant-dosage/"
m_filename=paste(sep="",m_path,snp,".txt")
write.table(cbind(subs,snp_risk_allele_dosage), m_filename,  row.names = FALSE, col.names=FALSE, quote=FALSE)

system(paste(sep="", "sort -k1,1 ", m_filename, paste(sep=""," >",m_filename,".sorted")))

m_filename=paste(sep="",m_path,snp,".variant")
write.table(snp_data, m_filename,  row.names = FALSE, col.names=T, quote=FALSE)
print(paste("Writing design mat for variant:", snp))
#system(paste('bash /data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/scripts/script-variant-design-matrix-generator.sh', m_filename))
system(paste('bash /data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/design-mat/scripts/script-variant-design-matrix-generator.sh', snp))
}
}


