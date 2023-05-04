tmp <- tempfile(tmpdir = paste(sep="","tmp-data-ld-hapmap3-processed"))
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
corr <- NULL

        ld <- readRDS("/data/clusterfs/lag/users/sousoh/ukbb/genetic/hapmap/hapmap3plus/map_hm3_plus.rds")
	ld<-ld$ld

for (chr in 1:22) {

if (chr == 1) {
print ("loading chr 1")
temp_mat <- readRDS(paste(sep="","/data/clusterfs/lag/users/sousoh/ukbb/genetic/hapmap/hapmap3plus/LD/LD_with_blocks_chr",chr,".rds"))
corr <- as_SFBM(temp_mat, tmp)
    } else {
	print(paste("loading chr",chr))
	temp_mat <- readRDS(paste(sep="","/data/clusterfs/lag/users/sousoh/ukbb/genetic/hapmap/hapmap3plus/LD/LD_with_blocks_chr",chr,".rds"))
    	corr$add_columns(temp_mat, nrow(corr))
    }
    print (paste("Chromosome", chr,"processed."))
}

saveRDS(corr, file = "/data/clusterfs/lag/users/sousoh/ukbb/genetic/bigsnp/processed/matrix_orig_corr_hapmap.rds")
print ("saved LD matrices.")

