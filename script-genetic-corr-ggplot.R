library(ggplot2)
library(scales)

#a<-read.table('/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-genetic-corr/all.matrix',header=F)
#h2<-read.table('/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-genetic-corr/all.h2snp',header=F)
a<-read.table('/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-genetic-corr/new/all-corr.mat',header=F)
h2<-read.table('/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-genetic-corr/new/all.h2snp',header=F)[,c(1,6)]

#selected_phenos <- h2$V1[c(2,16:18,20:22,23:29,31,32,33)]
#selected_phenos <- h2$V1[c(16,20:22,26:29,32,33)]

selected_phenos <- h2$V1

a <- a[a$V1 %in% selected_phenos & a$V2 %in% selected_phenos,]

h2 <- h2 [h2$V1 %in% selected_phenos,]

b<-a

pheno_code <- levels(factor(c(b$V1,b$V2)))
reorder_ids <-  match(c("dyslexia","UKBB-REACTION-TIME","total-surf","grade-E1","grade-E2","grade-E3","grade-E4","AD","PD","ASD","ADHD","SCZ","BIP","UKBB-GCSE","UKBB-risk-taking","UKBB-PAIN","UKBB-VNR"), pheno_code)
reorder_ids <-  match(c("dyslexia","grade-E1","AD","PD","ASD","ADHD","SCZ","BIP","UKBB-GCSE","UKBB-VNR"), pheno_code)
reorder_ids <-  match(c("dyslexia","grade-E1","ADHD","UKBB-GCSE","UKBB-VNR","else-NREAD","else-PA","else-WR","else-SP"), pheno_code)

pheno_code<-pheno_code[reorder_ids]

b$V1 <- factor(a$V1, levels=pheno_code)
b$V2<-factor(a$V2,levels=pheno_code)
h2$condition <-factor(h2$V1,  levels=pheno_code)

n_pheno <-  dim(h2)[1]

 levels(b$V1) <- 1:( n_pheno +1)
 levels(b$V2) <- 1:( n_pheno +1)

 levels(h2$condition) <- 1:( n_pheno +1)


 b$col1 <- pmax(as.numeric(as.character(b$V1)),as.numeric(as.character(b$V2)))
 b$col2 <- pmin(as.numeric(as.character(b$V1)),as.numeric(as.character(b$V2)))
h2$condition <- as.numeric(as.character(h2$condition))

pheno_code[which(pheno_code=="UKBB-VNR")] <- "Fluid intelligence"
pheno_code[which(pheno_code=="UKBB-risk-taking")] <- "Risk taking"
pheno_code[which(pheno_code=="UKBB-REACTION-TIME")] <- "Reaction time"
pheno_code[which(pheno_code=="UKBB-PAIN")] <- "Pain"
pheno_code[which(pheno_code=="UKBB-GCSE")] <- "Education (GCSE)"
pheno_code[which(pheno_code=="total-surf")] <- "Surface area (total)"
pheno_code[which(pheno_code=="grade-E1")] <- "School grade"
pheno_code[which(pheno_code=="grade-E2")] <- "School grade PC2"
pheno_code[which(pheno_code=="grade-E3")] <- "School grade PC3"
pheno_code[which(pheno_code=="grade-E4")] <- "School grade PC4"
pheno_code[which(pheno_code=="else-WR")] <- "Word reading"
pheno_code[which(pheno_code=="else-PA")] <- "Phonemic awareness"
pheno_code[which(pheno_code=="else-NREAD")] <- "Non-word reading"
pheno_code[which(pheno_code=="else-SP")] <- "Spelling"


# geom_point(data=b,shape=22,aes(x=col1,y=n_pheno-col2+1,size=pmin(1,abs(V3)),fill=V3 ) )+
my_plot <-  ggplot()+ geom_point(data=b,shape=22,size=12, aes(y=n_pheno-col1+1,x=col2,fill=V3 ) )+ 
 geom_point(data=b,shape=22,size=12, aes(x=col1,y=n_pheno-col2+1,fill=V3 ) )+
 geom_point(data=h2,shape=21,aes(x=condition,y=n_pheno-condition+1,size=abs(as.numeric(V6)),fill=abs(as.numeric(V6)) ) )+
 theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), axis.title.x=element_blank(),axis.title.y=element_blank() ,
       panel.background = element_rect(fill = "gray95",
                                colour = "darkgrey",
                                size = 1, linetype = "solid"),
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                colour = "darkgrey") )+
scale_y_continuous(labels=    as.character(rev(pheno_code)) , breaks=1:(n_pheno), minor_breaks=1.5:(n_pheno+0.5), guide = guide_axis(angle = 45), limits=c(0.5,n_pheno+0.5),expand=c(0,0)) +  
	 scale_x_continuous(labels=as.character(pheno_code), breaks=1:n_pheno, minor_breaks=1.5:(n_pheno+0.5), guide = guide_axis(angle = 45), limits=c(0.5,n_pheno+0.5),expand=c(0,0)) +
	 scale_fill_gradientn( name="Genetic\ncorrelation (r)",limits = c(-1,1), colours=c("blue",rgb(0,191/255,255/255), "white","yellow","red"), space = "Lab", breaks=c(-1,-0.5,0,0.5,1)  ,oob=squish)  +
guides(size=FALSE)+
annotate(geom="rect",xmin = 1.5+0.02, xmax = 0.5+n_pheno-0.02, ymin = 0.5+n_pheno  + 0.02 - 1, ymax = 0.5+n_pheno - 0.02,alpha=0.2,fill="white",color="red",size=1)

ggsave("/data/clusterfs/lag/users/sousoh/ukbb/genetic/sbayesR/genetic_corr.png",my_plot,width=150,height=119,units="mm",dpi=1200)

#+	 scale_size_continuous( name="|r|", breaks=c(0,0.25,0.5,1), legend.position="na")





