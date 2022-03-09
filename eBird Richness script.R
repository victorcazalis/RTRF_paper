chlist$alt<-as.numeric(scale(chlist$alt))
chlist$alt2<-chlist$alt^2
chlist$npp<-as.numeric(scale(chlist$npp))
chlist$npp2<-chlist$npp^2
chlist<-subset(chlist, is.na(npp)==F & is.na(hfp)==F & is.na(alt)==F & year>=2010)
chlist$Exotic[chlist$htsp %in% c("EAS", "IND")]<-NA

DF<-data.frame(hfp=seq(0,50,0.5), Rich1=NA, Rich2=NA, Rich3=NA, ForDep=NA, Endem=NA, Threat=NA, duration=median(chlist$duration, na.rm=T), alt=median(chlist$alt, na.rm=T), alt2=median(chlist$alt, na.rm=T)^2, npp=median(chlist$npp, na.rm=T), npp2=median(chlist$npp, na.rm=T)^2, day=median(chlist$day, na.rm=T), obsKelling=median(chlist$obsKelling, na.rm=T), year=median(chlist$year, na.rm=T), N_obs=median(chlist$N_obs, na.rm=T), htsp="ATL")


### Supplementary descriptive plots
# Bird diversity
Ga<-ggplot(chlist)+geom_histogram(aes(x=rich))+xlab("Species richness")+theme_minimal()
Gb<-ggplot(chlist)+geom_histogram(aes(x=Ass.specialisation1))+xlab("High-dependency")+theme_minimal()
Gc<-ggplot(chlist)+geom_histogram(aes(x=(Ass.specialisation2-Ass.specialisation1)))+xlab("Medium-dependency")+theme_minimal()
Gd<-ggplot(chlist)+geom_histogram(aes(x=LowForest))+xlab("Low-dependency")+theme_minimal()
Ge<-ggplot(chlist)+geom_histogram(aes(x=NonForest))+xlab("Non-forest")+theme_minimal()
Gf<-ggplot(chlist)+geom_histogram(aes(x=endem.index))+xlab("Endemic species")+theme_minimal()
Gg<-ggplot(chlist)+geom_histogram(aes(x=LargeRange))+xlab("Large-range species")+theme_minimal()
Gh<-ggplot(chlist)+geom_histogram(aes(x=iucn.index))+xlab("Threatened + NT species")+theme_minimal()
Gi<-ggplot(chlist)+geom_histogram(aes(x=Exotic))+xlab("Non-native species")+theme_minimal()
Gj<-ggplot(chlist)+geom_histogram(aes(x=SensiHigh))+xlab("High-sensitivity species")+theme_minimal()
Gk<-ggplot(chlist)+geom_histogram(aes(x=SensiTol))+xlab("Tolerant species")+theme_minimal()
Gl<-ggplot(chlist)+geom_histogram(aes(x=SensiAnthropo))+xlab("Anthropophilic species")+theme_minimal()

Gdistri<-cowplot::plot_grid(Ga,Gb,Gc,Gd,Ge,Gf,Gg,Gh,Gi,Gj,Gk,Gl, ncol=4, labels="AUTO")
cowplot::save_plot("Figures/Supplementary/Distribution of diversity indices.pdf", Gdistri, base_width=10, base_height=5)




### Ecological model
chlist$htsp2<-plyr::revalue(chlist$htsp, c("ATL"="Atlantic Forest", "AND"="Tropical Andes", "TUM"="Tumbes-Choco-Magdalena", "MES"="Mesoamerica", "EAS"="Eastern Afromontane", "GHA"="Western Ghats and Sri Lanka", "IND"="Indo-Burma", "SUN"="Sundaland"))
chlist$htsp2<-factor(chlist$htsp2, c("Atlantic Forest", "Tropical Andes", "Tumbes-Choco-Magdalena", "Mesoamerica", "Eastern Afromontane", "Western Ghats and Sri Lanka", "Indo-Burma", "Sundaland"))

M0<-gam(rich ~ alt + alt2 + npp + npp2 + s(day)+ s(duration, k=4) + s(obsKelling, k=4) + year + s(N_obs, k=4) + htsp, data=chlist, family="nb")


### Footprint model on residuals
# Extract residuals
chlist$Residus<-residuals(M0)

# GAM (for SI)
M3<-gam(Residus ~ s(hfp,k=6), data=chlist, family="gaussian")
DF$Rich3<-exp(predict(M0, newdata=DF, type="link") + predict(M3, newdata=DF, type="response"))

# Fit threshold regression
library(chngpt)
M4<-chngptm(formula.1=Residus ~ 1, formula.2= ~hfp, chlist, type="segmented", family="gaussian")
DF$Rich4<-exp(predict(M0, newdata=DF, type="link") + predict(M4, newdata=DF, type="response"))

# Save coefficients
M4_lin<-lm(Residus ~ hfp, data=chlist)
coef[coef$Index=="Overall" & coef$Zone=="Global", ColNum]<-save_coef(chlist, M0, M4, M4_lin, "Residus", "Threshold")


### Fig.1A
His1<-ggplot(chlist)+
  geom_histogram(aes(x=hfp), breaks=seq(0,50,2))+
  theme_void()+
  scale_x_continuous(limits=c(0,50), expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))

G1<-ggplot(DF, aes(x=hfp))+
  geom_line(aes(y=Rich4), col="#0606f0ff", size=1.5)+
  scale_y_continuous(limits=c(0,1.2*max(DF$Rich4)))+
  xlab("Human footprint") + ylab("Species richness")+
  theme_bw()+
  ggtitle("eBird (N=65,465)")+
  annotation_custom(ggplotGrob(His1), xmin=0, xmax=50, ymin=0, ymax=0.2*max(DF$Rich4))



### Per hotspot ###
DFhtsp<-list()
colour<-NA

for(i in 1:8){
  HT=c("ATL", "AND", "TUM", "MES", "EAS", "GHA", "IND", "SUN")[i]
  
  chlisthtsp<-subset(chlist, chlist$htsp==HT)
  
  DFhtsp[[i]]<-data.frame(hfp=seq(0,50,0.5), Rich1=NA, Rich2=NA, Rich3=NA, ForDep=NA, Endem=NA, Threat=NA, duration=median(chlisthtsp$duration, na.rm=T), alt=median(chlisthtsp$alt, na.rm=T), alt2=median(chlisthtsp$alt, na.rm=T)^2, npp=median(chlisthtsp$npp, na.rm=T), npp2=median(chlisthtsp$npp, na.rm=T)^2, day=median(chlisthtsp$day, na.rm=T), obsKelling=median(chlisthtsp$obsKelling, na.rm=T), year=median(chlisthtsp$year, na.rm=T), N_obs=median(chlisthtsp$N_obs, na.rm=T))
  
  DFhtsp[[i]]$htsp<-HT
  
  # Ecological model
  M0htsp<-gam(rich ~ alt + alt2 + npp + npp2 + s(day)+ s(duration, k=4) + s(obsKelling, k=4) + year + s(N_obs, k=4), data=chlisthtsp, family="nb")
  chlisthtsp$Residus<-residuals(M0htsp)
  
  # Best model (threshold if significant, linear otherwise)
  M4htsp<-chngptm (formula.1=Residus ~ 1, formula.2= ~hfp, chlisthtsp, type="segmented", family="gaussian")
  
  # Save coefficients
  M1htsp<-lm(Residus ~ hfp, data=chlisthtsp)
  Best_mod_value=revalue(as.character(as.data.frame(summary(M4htsp)[1])[,5][3]<0.05), c("TRUE"="Threshold", "FALSE"="Linear"))
  coef[coef$Index=="Overall" & coef$Zone==HT, ColNum]<-save_coef(chlisthtsp, M0htsp, M4htsp, M1htsp, "Residus", Best_mod_value)
  
  
  if(as.numeric(coef$P2[coef$Index=="Overall" & coef$Zone==HT])<0.05){ # Stock curve for hotspots with significant threshold
    DFhtsp[[i]]$Rich4<-exp(predict(M0htsp, newdata=DFhtsp[[i]], type="link") + predict(M4htsp, newdata=DFhtsp[[i]], type="response"))
  } else { # Stock curve from linear trend for non-significant threshold
    M1htsp<-gam(Residus ~ hfp, data=chlisthtsp)
    coef$Coef1[coef$Index=="Overall" & coef$Zone==paste0(HT, "lin")]<-round(M1htsp$coefficients["hfp"],3)
    coef$P1[coef$Index=="Overall" & coef$Zone==paste0(HT, "lin")]<-summary(M1htsp)$p.pv["hfp"]
    DFhtsp[[i]]$Rich4<-exp(predict(M0htsp, newdata=DFhtsp[[i]], type="link") + predict(M1htsp, newdata=DFhtsp[[i]], type="response"))
  }
  
  # Smoothed model (Supplementary)
  M3htsp<-gam(Residus ~ s(hfp, k=6), data=chlisthtsp)
  DFhtsp[[i]]$Rich3<-exp(predict(M0htsp, newdata=DFhtsp[[i]], type="link") + predict(M3htsp, newdata=DFhtsp[[i]], type="response"))
  
  # Save colour to plot
  if(as.numeric(coef$P2[coef$Index=="Overall" & coef$Zone==HT])<0.05){colour[i]<-"#0606f0ff"} else {colour[i]<-"#7570b3ff"}
  
}