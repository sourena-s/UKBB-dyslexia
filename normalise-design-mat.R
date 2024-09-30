args = commandArgs(trailingOnly=TRUE)

m_filename=args[1]
a<-read.table(m_filename,header=F)

#normalise all columns except the last which is the intercept
a<- cbind (scale(a[,1:dim(a)[2]-1])  , a[, dim(a)[2] ]  )
write.table(a,file=paste(sep="", m_filename, ".norm"), quote=F,col.names=F,row.names=F)
