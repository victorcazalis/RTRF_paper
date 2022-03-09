


Routes<-subset(Routes, is.na(Alt_mean)==F & is.na(npp)==F & is.na(hfp)==F)



### Model
Routes$Alt_mean<-as.numeric(scale(Routes$Alt_mean))
Routes$Alt_mean2<-Routes$Alt_mean^2
Routes$npp<-as.numeric(scale(Routes$npp))
Routes$npp2<-Routes$npp^2
Routes$hfp2<-Routes$hfp^2
Routes$Latitude<-as.numeric(scale(Routes$Latitude))
Routes$Latitude2<-Routes$Latitude^2

M0bbs<-gam(sp_richness.NONNAT ~ Alt_mean + Alt_mean2 + npp + npp2 + Latitude + Latitude2, data=Routes, family="nb")
Routes$Residus<-residuals(M0bbs)

library(chngpt)
M4bbs<-chngptm (formula.1=Residus ~ 1, formula.2= ~hfp, Routes, type="segmented", family="gaussian")

summary(M4bbs)


# Save coefficients
Mbbs_lin<-lm(Residus ~ hfp, data=Routes)
coef[coef$Index=="Overall" & coef$Zone=="BBS", ColNum]<-save_coef(Routes, M0bbs, M4bbs, Mbbs_lin, "Residus", "Threshold")


## Predict
DFbbs<-data.frame(hfp=seq(0,50,0.5), Alt_mean=median(Routes$Alt_mean, na.rm=T), Alt_mean2=median(Routes$Alt_mean, na.rm=T)^2, npp=median(Routes$npp, na.rm=T), npp2=median(Routes$npp, na.rm=T)^2, Latitude=median(Routes$Latitude, na.rm=T), Latitude2=median(Routes$Latitude, na.rm=T)^2)


DFbbs$Rich4<-exp(predict(M0bbs, newdata=DFbbs, type="link") + predict(M4bbs, newdata=DFbbs, type="response"))


Hisbbs<-ggplot(Routes)+
  geom_histogram(aes(x=hfp), breaks=seq(0,50,2))+
  theme_void()+
  scale_x_continuous(limits=c(0,50), expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))


GBBS<-ggplot()+
  geom_line(aes(x=hfp, y=Rich4), data=DFbbs, col="#0606f0ff", size=1.3)+
  scale_y_continuous(limits=c(0, 1.2*max(DFbbs$Rich4, na.rm=T)))+
  xlab("Human footprint") + ylab("Species richness")+ggtitle("BBS database (N=3016)")+
  theme_minimal()+
  annotation_custom(ggplotGrob(Hisbbs), xmin=0, xmax=50, ymin=0, ymax=0.2*max(DFbbs$Rich4))



