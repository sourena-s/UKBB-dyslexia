library(ggplot2)
library(scales)

a<-read.table('/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-genetic-corr/all.matrix',header=F)
h2<-read.table('/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-genetic-corr/all.h2snp',header=F)

b<-a

pheno_code <- levels(factor(c(b$V1,b$V2)))
b$V1 <- factor(a$V1, levels=pheno_code)
b$V2<-factor(a$V2,levels=pheno_code)
h2$condition <-factor(h2$V1,  levels=pheno_code)

 levels(b$V1) <- 1:34
 levels(b$V2) <- 1:34
 levels(h2$condition) <- 1:34

 b$col1 <- pmax(as.numeric(as.character(b$V1)),as.numeric(as.character(b$V2)))
 b$col2 <- pmin(as.numeric(as.character(b$V1)),as.numeric(as.character(b$V2)))
h2$condition <- as.numeric(as.character(h2$condition))



 ggplot()+ geom_point(data=b,shape=22,aes(x=35-col1,y=col2,size=pmin(1,abs(V3)),fill=V3 ) )+ 
 geom_point(data=b,shape=22,aes(y=col1,x=35-col2,size=pmin(1,abs(V3)),fill=V3 ) )+
 geom_point(data=h2,shape=21,aes(y=condition,x=35-condition,size=abs(as.numeric(V2)),fill=abs(as.numeric(V2)) ) )+
 theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), axis.title.x=element_blank(),axis.title.y=element_blank() ,
       panel.background = element_rect(fill = "gray95",
                                colour = "darkgrey",
                                size = 1, linetype = "solid"),
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                colour = "darkgrey") )+
scale_y_continuous(labels=    as.character(pheno_code) , breaks=1:34, minor_breaks=1.5:33.5, guide = guide_axis(angle = 45), limits=c(0.5,34.5),expand=c(0,0)) +  
	 scale_x_continuous(labels=rev(as.character(pheno_code)), breaks=1:34, minor_breaks=1.5:33.5, guide = guide_axis(angle = 45), limits=c(0.5,34.5),expand=c(0,0)) +
	 scale_fill_gradientn( name="Genetic correlation (r)",limits = c(-1,1), colours=c("blue",rgb(0,191/255,255/255), "white","yellow","red"), space = "Lab", breaks=c(-1,-0.5,0,0.5,1)  ,oob=squish) + scale_size_continuous( name="|r|", breaks=c(0,0.25,0.5,1))

