library(ggplot2)
library(ggrepel)
 m_df<- read.table('/data/clusterfs/lag/users/sousoh/ukbb/genetic/vaiant-maps/impact-mode-merged-var-maps/dyslexia-top-1percent/tbm/ica-dim-10/melodic_mix_norm',header=F)


#for structural MRI
#label_df<- t(c("MAPT\nrs12150530",max(m_df$V8),0))
#label_df<- rbind(label_df, c("RSPO3\nrs13220350",min(m_df$V8),0) )

#for diffusion MRI
label_df<- c("NEUROD2\nrs12453682",max(m_df$V10),0)
label_df<- rbind(label_df, c("SLC39A8\nrs35518360",min(m_df$V10),0) )
label_df<-as.data.frame(label_df)
colnames(label_df) <- c("txt", "x" ,"y")

label_df$x<-as.numeric(label_df$x)
label_df$y<-as.numeric(label_df$y)

m_plot <- ggplot(m_df) + geom_histogram(size=1,aes(x=V10),bins=100,color='black',fill='white') + geom_label_repel(size=20,data=label_df, aes(label=txt,x=x,y=y),nudge_y=300) + xlab("Impact mode weight (z-score)") + ylab("Frequency") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=1))  
ggsave("plot_weigth.png",m_plot, width=60,height=25,units="cm",dpi=300)
