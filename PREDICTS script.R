

# Calculate square variables
dataPRED$npp_sc<-scale(dataPRED$npp)
dataPRED$npp_sc2<-dataPRED$npp_sc^2
dataPRED$Latitud_sc<-scale(dataPRED$Latitud)
dataPRED$Latitud_sc2<-dataPRED$Latitud_sc^2
dataPRED$elevation_sc<-scale(dataPRED$elevation)
dataPRED$elevation_sc2<-dataPRED$elevation_sc^2


### Model 
#stock<-dataPRED
#dataPRED<-stock[stock$taxon_of_interest=="Invertebrates",] 


library(lme4)

M0PRED<-glmer(Species_richness ~ elevation_sc + elevation_sc2 + npp_sc + npp_sc2 + Latitud_sc + Latitud_sc2 + (1|SS)+ (1|SSB)  + (1|SSBS), data=dataPRED, family="poisson")

dataPRED$ResidusPRED<-residuals(M0PRED)


library(chngpt)
M4PRED<-chngptm (formula.1=ResidusPRED ~ 1, formula.2= ~hfp, dataPRED, type="segmented", family="gaussian")
summary(M4PRED)

M1PRED<-lm(ResidusPRED ~ hfp, data=dataPRED)
summary(M1PRED)


DFPRED<-data.frame(hfp=seq(0,50,0.5), elevation_sc=median(dataPRED$elevation_sc, na.rm=T), elevation_sc2=median(dataPRED$elevation_sc, na.rm=T)^2, SS=dataPRED$SS[1], SSB=dataPRED$SSB[1], SSBS=dataPRED$SSBS[1], Latitud_sc=median(dataPRED$Latitud_sc, na.rm=T), Latitud_sc2=median(dataPRED$Latitud_sc, na.rm=T)^2, npp_sc=median(dataPRED$npp_sc, na.rm=T), npp_sc2=median(dataPRED$npp_sc, na.rm=T)^2)

DFPRED$Rich0<-predict(M0PRED, DFPRED, type="link")
DFPRED$Rich4<-exp(DFPRED$Rich0 + predict(M4PRED, DFPRED, type="response"))
DFPRED$Rich1<-exp(DFPRED$Rich0 + predict(M1PRED, DFPRED, type="response"))

# Save coefficients
coef[coef$Index=="Overall" & coef$Zone=="PREDICTS", ColNum]<-save_coef(dataPRED, M0PRED, M4PRED, M1PRED, "ResidusPRED", "Linear")



# Plot
HisPRED<-ggplot(dataPRED)+
  geom_histogram(aes(x=hfp), breaks=seq(0,50,2))+
  theme_void()+
  scale_x_continuous(limits=c(0,50), expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))


GPred<-ggplot()+
  geom_line(aes(x=hfp, y=Rich1), data=DFPRED, col="#7570b3ff", size=1.3)+
  scale_y_continuous(limits=c(0, 1.2*max(DFPRED$Rich1, na.rm=T)))+
  xlab("Human footprint") + ylab("Species richness")+ggtitle("PREDICTS database (N=6382)")+
  theme_minimal()+
  annotation_custom(ggplotGrob(HisPRED), xmin=0, xmax=50, ymin=0, ymax=0.2*max(DFPRED$Rich4))




### SUPPLEMENTARY SCRIPT BIT to isolate PREDICTS data on birds, or tropical birds
# pred=readRDS("database.rds") # Table downloaded here to complete the Gray et al dataset with more precise taxonomic data: https://data.nhm.ac.uk/dataset/the-2016-release-of-the-predicts-database
# pred=subset(pred, pred$Source_ID %in% dataPRED$Source_ID) # Only keep studies present in Gray database
# studies_classes=ddply(pred, .(Source_ID, Study_number, Site_number), function(x){data.frame(Birds=("Aves" %in% x$Class), Classes=paste(unique(x$Class), collapse=" "))}) # Isolate sites working only on birds
# aves_studies=studies_classes[studies_classes$Classes=="Aves",]
# dataPRED_birds=dataPRED[paste(dataPRED$Source_ID, dataPRED$Study_number, dataPRED$Site_number, sep=" ") %in% paste(aves_studies$Source_ID, aves_studies$Study_number, aves_studies$Site_number, sep=" "),] # Keep Gray et al data only for those sites
# dataPRED_TropBirds=subset(dataPRED_birds, dataPRED_birds$Zone=="Tropical")
