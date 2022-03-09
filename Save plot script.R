
DFhtsp2<-rbind(DFhtsp[[1]], DFhtsp[[2]], DFhtsp[[3]], DFhtsp[[4]], DFhtsp[[5]], DFhtsp[[6]], DFhtsp[[7]], DFhtsp[[8]])

DFhtsp2$htsp2<-plyr::revalue(DFhtsp2$htsp, c("ATL"="Atlantic Forest", "AND"="Tropical Andes", "TUM"="Tumbes-Choco-Magdalena", "MES"="Mesoamerica", "EAS"="Eastern Afromontane", "GHA"="Western Ghats and Sri Lanka", "IND"="Indo-Burma", "SUN"="Sundaland"))
DFhtsp2$htsp2<-factor(DFhtsp2$htsp2, c("Atlantic Forest", "Tropical Andes", "Tumbes-Choco-Magdalena", "Mesoamerica", "Eastern Afromontane", "Western Ghats and Sri Lanka", "Indo-Burma", "Sundaland"))

# Hotspot plots
GHtsp<-list()
HisHt<-list()

for(HT in 1:8){
  hotname<-c("ATL", "AND", "TUM", "MES", "EAS", "GHA", "IND", "SUN")[HT]
  
  HisHt[[HT]]<-ggplot(chlist[chlist$htsp==hotname,])+geom_histogram(aes(x=hfp), breaks=seq(0,50,2))+ theme_void()+scale_x_continuous(limits=c(0,50), expand=c(0,0))+scale_y_continuous(expand=c(0,0))
  
  GHtsp[[HT]]<-ggplot(DFhtsp2[DFhtsp2$htsp==hotname,])+geom_line(aes(x=hfp, y=Rich4), col=as.character(colour[HT]), size=1.7)+  scale_y_continuous(limits=c(0,1.2*max(DFhtsp2$Rich4)))+  xlab("") + ylab("")+theme_minimal()+
    annotation_custom(ggplotGrob(HisHt[[HT]]), xmin=0, xmax=50, ymin=0, ymax=0.2*max(DFhtsp2$Rich4[DFhtsp2$htsp==hotname]))+ggtitle(hotname)
}

Ghot<-gridExtra::grid.arrange(GHtsp[[1]], GHtsp[[2]], GHtsp[[3]], GHtsp[[4]], GHtsp[[5]], GHtsp[[6]], GHtsp[[7]], GHtsp[[8]], ncol=4)



# Main plot
Gtot<-gridExtra::grid.arrange(G1, 
                              Ghot,
                              GPred, GBBS,
                              layout_matrix=matrix(c(1,1,1,1,1, 5,5,5, 1,1,1,1,1, 5,5,5, 1,1,1,1,1, 5,5,5,rep(2,24), 3,3,3,3,4,4,4,4,3,3,3,3,4,4,4,4), ncol=8, byrow=T))


