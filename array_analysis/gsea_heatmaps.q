FWER.cutoff<-0.01

sigGeneSets<-unique(data$NAME[data$FWER.p.val<=FWER.cutoff])

data$sigGeneSet<-F

data$sigGeneSet[data$NAME %in% sigGeneSets]<-T

data.sig<-subset(data,sigGeneSet==T)
data.sig$NAME<-factor(data.sig$NAME,levels=sigGeneSets)

#data.melt<-melt(data.sig,id=c(1,3,11))

p<-ggplot(data.sig) + scale_y_discrete(expand = c(0, 0)) + scale_x_discrete(expand=c(0,0)) +
	opts(
	       plot.title = theme_text(size = 15), 
	       axis.title.x = theme_text(size = 12),
           axis.text.x = theme_text(size=10,angle = -90, hjust=0),
	       axis.title.y = theme_text(size = 12, angle = 90), 
	       strip.text.x = theme_text(size = 12), 
	       strip.text.y = theme_text(size = 12, angle = -90),
	       axis.text.y = theme_text(size=5,hjust=1),
	       	#axis.text.y = theme_blank(),
			panel.background = theme_rect(colour=NA),
	       panel.grid.major = theme_line(colour=NA),
	       panel.grid.minor = theme_line(colour=NA),
	       panel.border = theme_rect(colour=NA)
	)

pdf("NES_heatmap_vs_GFP.pdf")
 p + geom_tile(aes(x=SAMPLE,y=NAME,fill=NES,alpha=-log10(FWER.p.val+0.0001)))  + opts(title="Significant GeneSet Enrichment\nby clone overexpression") +
	scale_fill_gradient2("Enrichment",low="blue",high="red")

dev.off()	



###############
#Reordering
###############
NES.melt<-melt(data.sig[,c(1,2,5)])
NES.cast<-cast(NES.melt,...~SAMPLE)
NES.rows<-NES.cast$NAME
NES.cast<-NES.cast[,-c(1:2)]

rowOrder<-order.dendrogram(as.dendrogram(hclust(dist(NES.cast))))
colOrder<-order.dendrogram(as.dendrogram(hclust(dist(t(NES.cast)))))

data.sig$NAME<-factor(data.sig$NAME,levels=NES.rows[rowOrder])
data.sig$SAMPLE<-factor(data.sig$SAMPLE,levels=colnames(NES.cast)[colOrder])

q<-ggplot(data.sig) + scale_y_discrete(expand = c(0, 0)) + scale_x_discrete(expand=c(0,0)) +
	opts(
	       plot.title = theme_text(size = 15), 
	       axis.title.x = theme_text(size = 12),
           axis.text.x = theme_text(size=10,angle = -90, hjust=0),
	       axis.title.y = theme_text(size = 12, angle = 90), 
	       strip.text.x = theme_text(size = 12), 
	       strip.text.y = theme_text(size = 12, angle = -90),
	       axis.text.y = theme_text(size=4,hjust=1),
	       	#axis.text.y = theme_blank(),
			panel.background = theme_rect(colour=NA),
	       panel.grid.major = theme_line(colour=NA),
	       panel.grid.minor = theme_line(colour=NA),
	       panel.border = theme_rect(colour=NA)
	)

pdf("NES_heatmap_vs_GFP_ordered.pdf")
q + geom_tile(aes(x=SAMPLE,y=NAME,fill=NES,alpha=-log10(FWER.p.val+0.0001)))  + opts(title="Significant GeneSet Enrichment\nby clone overexpression") +
	scale_fill_gradient("Enrichment",low="blue",high="red") + scale_alpha("Significance")


dev.off()