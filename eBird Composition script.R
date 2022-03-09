### FOREST-DEPENDENCY
## High
M_spe_nat<-gam(Ass.specialisation1 ~ alt + alt2 + npp + npp2 + s(day)+ s(duration, k=4) + s(obsKelling, k=4) + year + s(N_obs, k=4) + htsp, data=chlist, family="nb") # NB fits better
chlist$Res_spe<-residuals(M_spe_nat)
M_spe4<-chngptm(formula.1=Res_spe ~ 1, formula.2= ~hfp, chlist, type="segmented", family="gaussian") ; summary(M_spe4)
M_spe4b<-lm(Res_spe~hfp, data=chlist)
if((as.data.frame(summary(M_spe4)[1])[,5][3])<0.05){M_spePRED<-M_spe4} else {M_spePRED<-M_spe4b}
DF$ForHigh4<-exp(predict(M_spe_nat, newdata=DF, type="link") + predict(M_spePRED, newdata=DF, type="response"))
M_spe4free<-gam(Res_spe ~ s(hfp,k=6), data=chlist)
DF$ForHigh4free<-exp(predict(M_spe_nat, newdata=DF, type="link") + predict(M_spe4free, newdata=DF, type="response"))

# Save coefficients
coef[coef$Index=="HighFor" & coef$Zone=="Global", ColNum]<-save_coef(chlist, M_spe_nat, M_spe4, M_spe4b, "Res_spe", "Linear")



## Medium
chlist$MediumForest<-chlist$Ass.specialisation2-chlist$Ass.specialisation1
M_formed_nat<-gam(MediumForest ~ alt + alt2 + npp + npp2 + s(day)+ s(duration, k=4) + s(obsKelling, k=4) + year + s(N_obs, k=4) + htsp, data=chlist, family="nb") # NB fits better
chlist$Res_formed<-residuals(M_formed_nat)
M_formed4<-chngptm(formula.1=Res_formed ~ 1, formula.2= ~hfp, chlist, type="segmented", family="gaussian") ; summary(M_formed4)
M_formed4b<-lm(Res_formed~hfp, data=chlist)
if((as.data.frame(summary(M_formed4)[1])[,5][3])<0.05){M_formedPRED<-M_formed4} else {M_formedPRED<-M_formed4b}
DF$ForMed4<-exp(predict(M_formed_nat, newdata=DF, type="link") + predict(M_formedPRED, newdata=DF, type="response"))
M_formed4free<-gam(Res_formed ~ s(hfp,k=6), data=chlist)
DF$ForMed4free<-exp(predict(M_formed_nat, newdata=DF, type="link") + predict(M_formed4free, newdata=DF, type="response"))

# Save coefficients
coef[coef$Index=="MedFor" & coef$Zone=="Global", ColNum]<-save_coef(chlist, M_formed_nat, M_formed4, M_formed4b, "Res_formed", "Threshold")



## Low
M_forlow_nat<-gam(LowForest ~ alt + alt2 + npp + npp2 + s(day)+ s(duration, k=4) + s(obsKelling, k=4) + year + s(N_obs, k=4) + htsp, data=chlist, family="nb") # NB fits better
chlist$Res_forlow<-residuals(M_forlow_nat)
M_forlow4<-chngptm(formula.1=Res_forlow ~ 1, formula.2= ~hfp, chlist, type="segmented", family="gaussian") ; summary(M_forlow4)
M_forlow4b<-lm(Res_forlow~hfp, data=chlist)
if((as.data.frame(summary(M_forlow4)[1])[,5][3])<0.05){M_forlowPRED<-M_forlow4} else {M_forlowPRED<-M_forlow4b}
DF$ForLow4<-exp(predict(M_forlow_nat, newdata=DF, type="link") + predict(M_forlowPRED, newdata=DF, type="response"))
M_forlow4free<-gam(Res_forlow ~ s(hfp,k=6), data=chlist)
DF$ForLow4free<-exp(predict(M_forlow_nat, newdata=DF, type="link") + predict(M_forlow4free, newdata=DF, type="response"))

# Save coefficients
coef[coef$Index=="LowFor" & coef$Zone=="Global", ColNum]<-save_coef(chlist, M_forlow_nat, M_forlow4, M_forlow4b, "Res_forlow", "Threshold")



## Non-forest
M_NF_nat<-gam(NonForest ~ alt + alt2 + npp + npp2 + s(day)+ s(duration, k=4) + s(obsKelling, k=4) + year + s(N_obs, k=4) + htsp, data=chlist, family="nb") # NB fits better
chlist$Res_NF<-residuals(M_NF_nat)
M_NF4<-chngptm(formula.1=Res_NF ~ 1, formula.2= ~hfp, chlist, type="segmented", family="gaussian") ; summary(M_NF4)
M_NF4b<-lm(Res_NF~hfp, data=chlist)
if((as.data.frame(summary(M_NF4)[1])[,5][3])<0.05){M_NFPRED<-M_NF4} else {M_NFPRED<-M_NF4b}
DF$NF4<-exp(predict(M_NF_nat, newdata=DF, type="link") + predict(M_NFPRED, newdata=DF, type="response"))
M_NF4free<-gam(Res_NF ~ s(hfp,k=6), data=chlist)
DF$NF4free<-exp(predict(M_NF_nat, newdata=DF, type="link") + predict(M_NF4free, newdata=DF, type="response"))

# Save coefficients
coef[coef$Index=="NoFor" & coef$Zone=="Global", ColNum]<-save_coef(chlist, M_NF_nat, M_NF4, M_NF4b, "Res_NF", "Threshold")



## Plot
G2a<-ggplot(DF, aes(x=hfp))+
  geom_line(aes(y=ForHigh4), col="#1b9e77", size=1.5)+
  geom_line(aes(y=ForMed4), col="#b8e186", size=1.5)+
  geom_line(aes(y=ForLow4), col="lightsalmon", size=1.5)+
  geom_line(aes(y=NF4), col="#d95f02", size=1.5)+
  scale_y_continuous(trans="log10")+
  xlab("Human footprint") + ylab("Species richness")+
  ggtitle("Forest specialisation")

G3a<-ggplot(DF, aes(x=hfp))+
  geom_line(aes(y=ForHigh4free), col="#1b9e77", size=1.5)+
  geom_line(aes(y=ForMed4free), col="#b8e186", size=1.5)+
  geom_line(aes(y=ForLow4free), col="lightsalmon", size=1.5)+
  geom_line(aes(y=NF4free), col="#d95f02", size=1.5)+
  scale_y_continuous(trans="log10")+
  xlab("Human footprint") + ylab("Species richness")+
  ggtitle("Forest specialisation")




### ENDEMISM ###

### Endemic
M_endem_nat<-gam(endem.index ~ alt + alt2 + npp + npp2 + s(day)+ s(duration, k=4) + s(obsKelling, k=4) + year + s(N_obs, k=4) + htsp, data=chlist, family="nb") # Fits better in NB
chlist$Res_endem<-residuals(M_endem_nat)
M_endem4<-chngptm(formula.1=Res_endem ~ 1, formula.2= ~hfp, chlist, type="segmented", family="gaussian") ; summary(M_endem4)
M_endem4b<-lm(Res_endem~hfp, data=chlist)
if((as.data.frame(summary(M_endem4)[1])[,5][3])<0.05){M_endemPRED<-M_endem4} else {M_endemPRED<-M_endem4b}
DF$Endem4<-exp(predict(M_endem_nat, newdata=DF, type="link") + predict(M_endemPRED, newdata=DF, type="response"))
M_endem4free<-gam(Res_endem ~ s(hfp,k=6), data=chlist)
DF$Endem4free<-exp(predict(M_endem_nat, newdata=DF, type="link") + predict(M_endem4free, newdata=DF, type="response"))

# Save coefficients
coef[coef$Index=="Endemic" & coef$Zone=="Global", ColNum]<-save_coef(chlist, M_endem_nat, M_endem4, M_endem4b, "Res_endem", "Threshold")



## Large range
M_large_nat<-gam(LargeRange ~ alt + alt2 + npp + npp2 + s(day)+ s(duration, k=4) + s(obsKelling, k=4) + year + s(N_obs, k=4) + htsp, data=chlist, family="nb") # NB fits better
chlist$Res_large<-residuals(M_large_nat)
M_large4<-chngptm(formula.1=Res_large ~ 1, formula.2= ~hfp, chlist, type="segmented", family="gaussian") ; summary(M_large4)
M_large4b<-lm(Res_large~hfp, data=chlist)
if((as.data.frame(summary(M_large4)[1])[,5][3])<0.05){M_largePRED<-M_large4} else {M_largePRED<-M_large4b}
DF$large4<-exp(predict(M_large_nat, newdata=DF, type="link") + predict(M_largePRED, newdata=DF, type="response"))
M_large4free<-gam(Res_large ~ s(hfp,k=6), data=chlist)
DF$large4free<-exp(predict(M_large_nat, newdata=DF, type="link") + predict(M_large4free, newdata=DF, type="response"))

# Save coefficients
coef[coef$Index=="LargeRange" & coef$Zone=="Global", ColNum]<-save_coef(chlist, M_large_nat, M_large4, M_large4b, "Res_large", "Threshold")


## Plot
G2b<-ggplot(DF, aes(x=hfp))+
  geom_line(aes(y=Endem4), col="#1b9e77", size=1.5)+
  geom_line(aes(y=large4), col="#d95f02", size=1.5)+
  scale_y_continuous(trans="log10")+
  xlab("Human footprint") + ylab("Species richness")+
  ggtitle("Endemism")

G3b<-ggplot(DF, aes(x=hfp))+
  geom_line(aes(y=Endem4free), col="#1b9e77", size=1.5)+
  geom_line(aes(y=large4free), col="#d95f02", size=1.5)+
  scale_y_continuous(trans="log10")+
  xlab("Human footprint") + ylab("Species richness")+
  ggtitle("Endemism")



### Threat / Exotic ###

## Threatened
M_threat_nat<-gam(iucn.index ~ alt + alt2 + npp + npp2 + s(day)+ s(duration, k=4) + s(obsKelling, k=4) + year + s(N_obs, k=4) + htsp, data=chlist, family="nb") # Fits better in NB
chlist$Res_threat<-residuals(M_threat_nat)
M_threat4<-chngptm(formula.1=Res_threat ~ 1, formula.2= ~hfp, chlist, type="segmented", family="gaussian") ; summary(M_threat4)
M_threat4b<-lm(Res_threat~hfp, data=chlist)
if((as.data.frame(summary(M_threat4)[1])[,5][3])<0.05){M_threatPRED<-M_threat4} else {M_threatPRED<-M_threat4b}
DF$threat4<-exp(predict(M_threat_nat, newdata=DF, type="link") + predict(M_threatPRED, newdata=DF, type="response"))
M_threat4free<-gam(Res_threat ~ s(hfp,k=6), data=chlist)
DF$threat4free<-exp(predict(M_threat_nat, newdata=DF, type="link") + predict(M_threat4free, newdata=DF, type="response"))

# Save coefficients
coef[coef$Index=="thrNT" & coef$Zone=="Global", ColNum]<-save_coef(chlist, M_threat_nat, M_threat4, M_threat4b, "Res_threat", "Threshold")



### Exotic
M_exotic_nat<-gam(Exotic ~ alt + alt2 + npp + npp2 + s(day)+ s(duration, k=4) + s(obsKelling, k=4) + year + s(N_obs, k=4) + htsp, data=chlist, family="nb") # NB fits better
chlist$Res_exotic[is.na(chlist$Exotic)==F]<-residuals(M_exotic_nat)
M_exotic4<-chngptm(formula.1=Res_exotic ~ 1, formula.2= ~hfp, chlist, type="segmented", family="gaussian") ; summary(M_exotic4)
M_exotic4b<-lm(Res_exotic~hfp, data=chlist)
if((as.data.frame(summary(M_exotic4)[1])[,5][3])<0.05){M_exoticPRED<-M_exotic4} else {M_exoticPRED<-M_exotic4b}
DF$exotic4<-exp(predict(M_exotic_nat, newdata=DF, type="link") + predict(M_exoticPRED, newdata=DF, type="response"))
M_exotic4free<-gam(Res_exotic ~ s(hfp,k=6), data=chlist)
DF$exotic4free<-exp(predict(M_exotic_nat, newdata=DF, type="link") + predict(M_exotic4free, newdata=DF, type="response"))

# Save coefficients
coef[coef$Index=="Exotic" & coef$Zone=="Global", ColNum]<-save_coef(chlist, M_exotic_nat, M_exotic4, M_exotic4b, "Res_exotic", "Threshold")



## Plot
G2c<-ggplot(DF, aes(x=hfp))+
  geom_line(aes(y=threat4), col="#1b9e77", size=1.5)+
  geom_line(aes(y=exotic4), col="#d95f02", size=1.5)+
  scale_y_continuous(trans="log10")+
  xlab("Human footprint") + ylab("Species richness")+
  ggtitle("Threatened / Exotic")

G3c<-ggplot(DF, aes(x=hfp))+
  geom_line(aes(y=threat4free), col="#1b9e77", size=1.5)+
  geom_line(aes(y=exotic4free), col="#d95f02", size=1.5)+
  scale_y_continuous(trans="log10")+
  xlab("Human footprint") + ylab("Species richness")+
  ggtitle("Threatened / Exotic")


### Sensitivity ###

## High-sensitivity
M_sensi_nat<-gam(SensiHigh ~ alt + alt2 + npp + npp2 + s(day)+ s(duration, k=4) + s(obsKelling, k=4) + year + s(N_obs, k=4) + htsp, data=chlist, family="nb") # NB fits better
chlist$Res_sensi[is.na(chlist$SensiHigh)==F]<-residuals(M_sensi_nat)
M_sensi4<-chngptm(formula.1=Res_sensi ~ 1, formula.2= ~hfp, chlist, type="segmented", family="gaussian") ; summary(M_sensi4)
M_sensi4b<-lm(Res_sensi~hfp, data=chlist)
if((as.data.frame(summary(M_sensi4)[1])[,5][3])<0.05){M_sensiPRED<-M_sensi4} else {M_sensiPRED<-M_sensi4b}
DF$sensi4<-exp(predict(M_sensi_nat, newdata=DF, type="link") + predict(M_sensiPRED, newdata=DF, type="response"))
M_sensi4free<-gam(Res_sensi ~ s(hfp,k=6), data=chlist)
DF$sensi4free<-exp(predict(M_sensi_nat, newdata=DF, type="link") + predict(M_sensi4free, newdata=DF, type="response"))

# Save coefficients
coef[coef$Index=="SensiHigh" & coef$Zone=="Global", ColNum]<-save_coef(chlist, M_sensi_nat, M_sensi4, M_sensi4b, "Res_sensi", "Linear")


## Tolerant
M_sensitol_nat<-gam(SensiTol ~ alt + alt2 + npp + npp2 + s(day)+ s(duration, k=4) + s(obsKelling, k=4) + year + s(N_obs, k=4) + htsp, data=chlist, family="nb") # NB fits better
chlist$Res_sensitol[is.na(chlist$SensiTol)==F]<-residuals(M_sensitol_nat)
M_sensitol4<-chngptm(formula.1=Res_sensitol ~ 1, formula.2= ~hfp, chlist, type="segmented", family="gaussian") ; summary(M_sensitol4)
M_sensitol4b<-lm(Res_sensitol~hfp, data=chlist)
if((as.data.frame(summary(M_sensitol4)[1])[,5][3])<0.05){M_sensitolPRED<-M_sensitol4} else {M_sensitolPRED<-M_sensitol4b}
DF$sensitol4<-exp(predict(M_sensitol_nat, newdata=DF, type="link") + predict(M_sensitolPRED, newdata=DF, type="response"))
M_sensitol4free<-gam(Res_sensitol ~ s(hfp,k=6), data=chlist)
DF$sensitol4free<-exp(predict(M_sensitol_nat, newdata=DF, type="link") + predict(M_sensitol4free, newdata=DF, type="response"))

# Save coefficients
coef[coef$Index=="SensiTol" & coef$Zone=="Global", ColNum]<-save_coef(chlist, M_sensitol_nat, M_sensitol4, M_sensitol4b, "Res_sensitol", "Threshold")


## Anthropophilic
M_sensiant_nat<-gam(SensiAnthropo ~ alt + alt2 + npp + npp2 + s(day)+ s(duration, k=4) + s(obsKelling, k=4) + year + s(N_obs, k=4) + htsp, data=chlist, family="nb") # NB fits better
chlist$Res_sensiant[is.na(chlist$SensiAnthropo)==F]<-residuals(M_sensiant_nat)
M_sensiant4<-chngptm(formula.1=Res_sensiant ~ 1, formula.2= ~hfp, chlist, type="segmented", family="gaussian") ; summary(M_sensiant4)
M_sensiant4b<-lm(Res_sensiant ~ hfp, data=chlist) 
if((as.data.frame(summary(M_sensiant4)[1])[,5][3])<0.05){M_sensiantPRED<-M_sensiant4} else {M_sensiantPRED<-M_sensiant4b}
DF$sensiant4<-exp(predict(M_sensiant_nat, newdata=DF, type="link") + predict(M_sensiantPRED, newdata=DF, type="response"))
M_sensiant4free<-gam(Res_sensiant ~ s(hfp,k=6), data=chlist)
DF$sensiant4free<-exp(predict(M_sensiant_nat, newdata=DF, type="link") + predict(M_sensiant4free, newdata=DF, type="response"))

# Save coefficients
coef[coef$Index=="SensiAnt" & coef$Zone=="Global", ColNum]<-save_coef(chlist, M_sensiant_nat, M_sensiant4, M_sensiant4b, "Res_sensiant", "Threshold")


## Plot
G2d<-ggplot(DF, aes(x=hfp))+
  geom_line(aes(y=sensi4), col="#1b9e77", size=1.5)+
  geom_line(aes(y=sensitol4), col="lightsalmon", size=1.5)+
  geom_line(aes(y=sensiant4), col="#d95f02", size=1.5)+
  scale_y_continuous(trans="log10")+
  xlab("Human footprint") + ylab("Species richness")+
  ggtitle("Sensitivity")

G3d<-ggplot(DF, aes(x=hfp))+
  geom_line(aes(y=sensi4free), col="#1b9e77", size=1.5)+
  geom_line(aes(y=sensitol4free), col="lightsalmon", size=1.5)+
  geom_line(aes(y=sensiant4free), col="#d95f02", size=1.5)+
  scale_y_continuous(trans="log10")+
  xlab("Human footprint") + ylab("Species richness")+
  ggtitle("Sensitivity")

GT<-gridExtra::grid.arrange(G2a, G2b, G2c, G2d, ncol=2)
GTsupp<-gridExtra::grid.arrange(G3a, G3b, G3c, G3d, ncol=2)