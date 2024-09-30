 library(ggplot2)
 library(data.table)
 library(cowplot)
 library(ggnewscale)


m_k=c(-1,-1,1,1,1,1,1,1,1)
id=0

for (trait in c("dyslexia","ADHD","UKBB-VNR", "UKBB-GCSE", "grade-E1","else-PA","else-WR","else-NREAD","else-SP")){
print(paste('analysing',trait,'please standby'))
	id=id+1
k=m_k[id] # minus one for 'negative' phenotypes, e.g. dyslexia, ADHD, positive for 'positive' phenotype, e.g. reading/reasoning performance
contrast="tbm"
 a<-read.table(paste0("/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/",trait,"/component-results/master_df.txt"),header=T)

a<-a[a$contrast==contrast,]

max_stat <- max(c(a$lassosum_tval,a$sbayesR_tval,a$cs1_tval,a$cs2_tval,a$cs3_tval,a$cs4_tval,a$cs5_tval))
min_stat<- min(c(a$lassosum_tval,a$sbayesR_tval,a$cs1_tval,a$cs2_tval,a$cs3_tval,a$cs4_tval,a$cs5_tval))
max_stat<-max(max_stat,abs(min_stat))
my_plots <- list()
 b=1
m_margin=c(1,1,1,1,0,0,0,0)
for (i in c(11, 18, 29, 47, 76, 124, 200 , 324)) {
my_plots[[b]] <- ggplot(a[a$contrast==contrast & a$XCA=="ica" & a$dim==i ,]) + 
	geom_point(shape=18, aes(x=ic_num,y=k*lassosum_tval,color=log(lambda) ,size=abs(lassosum_tval) ) ) +  
	geom_hline(yintercept=0,color="black",size=0.5) +
	ylim(-max_stat,max_stat) + xlab(paste ("ICA",i) ) + 
	scale_size(range=c(0.3,1.5)) + new_scale('size') +
	geom_point(shape=23, aes(ic_num,y=k*sbayesR_tval,size=abs(sbayesR_tval)) ,fill="green") + 
	geom_point(shape=23, aes(x=ic_num,y=k*cs1_tval,  size=abs(cs1_tval)    ) ,fill="red") + 
	scale_size(range=c(0.3,3)) + 
	theme(panel.border = element_blank(), 
			  axis.line.x = element_line(size = 1, linetype = "solid", colour = "grey"), 
			  panel.background = element_blank(), 
			  axis.title.x=element_text(angle = 90,hjust=0,vjust=0.5),  
			  legend.position = "none",
			  axis.title.y=element_blank(),
			  axis.text.y=element_blank(),
			  axis.ticks=element_blank(),
			  axis.text.x=element_blank() ,
			  plot.margin = unit(c(0,0,0,0), "mm")  ) +
		expand_limits(x=c(  m_margin[b]*( 1- (i*0.2)) , m_margin[b]* i*1.2))



b=b+1
 }

my_plots[[b]] <-  ggplot(a[a$contrast==contrast & a$XCA=="pca" & a$dim==324 ,]) +
	geom_point(shape=18, aes(x=ic_num,y=-k*lassosum_tval,color=log(lambda) ,size=abs(lassosum_tval) ) ) +  
	xlab("PCA 324") + ylab("t-value")  + 
	geom_hline(yintercept=0,color="black",size=0.5) +
	ylim(-max_stat,max_stat) + 
	scale_size(range=c(0.3,1.5)) + new_scale('size') +
	geom_point(shape=23, aes(ic_num,y=-k*sbayesR_tval,size=abs(sbayesR_tval)) ,fill="green") + 
	geom_point(shape=23, aes(x=ic_num,y=-k*cs1_tval,  size=abs(cs1_tval)    ) ,fill="red") + 
	scale_size(range=c(0.3,3)) + 
	theme(panel.border = element_blank(), 
			  axis.line.x = element_line(size = 1, linetype = "solid", colour = "grey"), 
			  panel.background = element_blank(), 
			  axis.title.x=element_text(angle = 90,hjust=0,vjust=0.5),  
			  legend.position = "none",
			  axis.ticks.x=element_blank(),
			  axis.text.x=element_blank() ,
			  plot.margin = unit(c(0,0,0,0), "cm")  ) #+
#		expand_limits(x=c(1- (324*0.2) , 324*1.2)


my_final_plot <- plot_grid(plotlist=my_plots[9:1],align="h",nrow=1, rel_widths=rev(c(47, 47, 47, 47, 76, 124, 200, 324, 324)) )
ggsave(paste0("/data/clusterfs/lag/users/sousoh/ukbb/misc-GWASs-pgs/",trait,"/idp_phewas_plot_",contrast,".png" ) ,plot=my_final_plot, bg="white", width=400,height=150,unit="mm",dpi=600)
 }
