library(mgcv) ; library(ggplot2) ; library(dplyr) ; library(gtools) ; library(plyr)
setwd("H:/RTRF_project")



### Create a table to store model results
coef<-data.frame(Index = c(rep("Overall", 11), "HighFor", "MedFor", "LowFor", "NoFor", "Endemic", "LargeRange", "thrNT", "Exotic", "SensiHigh", "SensiTol", "SensiAnt", rep("Overall", 8)), 
               Zone = c("Global", "ATL", "AND", "TUM", "MES", "EAS", "GHA", "IND", "SUN", "PREDICTS", "BBS", rep("Global", 11), "ATLlin", "ANDlin", "TUMlin", "MESlin", "EASlin", "GHAlin", "INDlin", "SUNlin"))
coef$R2_footprint<-coef$R2_ecol<-coef$AIC_lin<-coef$AIC_thresh<-coef$SSE_linear<-coef$SSE_thresh<-coef$SE_thresh<-coef$SE_inter<-coef$SE_hft<-coef$P2<-coef$Coef2<-coef$Threshold<-coef$P1<-coef$Coef1<-NA
               


### Create a function to calculate SSE of linear and threshold models
ColNum<-c(which(names(coef)=="Coef1"), which(names(coef)=="Coef2"), which(names(coef)=="Threshold"), which(names(coef)=="P1"), which(names(coef)=="P2"), which(names(coef)=="SE_hft"), which(names(coef)=="SE_inter"), which(names(coef)=="SE_thresh"), which(names(coef)=="SSE_thresh"), which(names(coef)=="SSE_linear"), which(names(coef)=="AIC_thresh"), which(names(coef)=="AIC_lin"), which(names(coef)=="R2_ecol"), which(names(coef)=="R2_footprint")) # Order of columns of coef in which I'll store the different values


save_coef<-function(data, model_ecol, model_thre, model_lin, Name_residuals, Best_model){
  
  data$Res_SSE<-data[,which(names(data)==Name_residuals)]
  
  data$Pred_thre<-predict(model_thre, data)
  data$Pred_lin<-predict(model_lin, data)
  
  data$SSE_thre<-(data$Res_SSE-data$Pred_thre)^2
  data$SSE_lin<-(data$Res_SSE-data$Pred_lin)^2
  
  # Calculate Rsquared of the ecological model (different for PREDICTS as I use glmer)
  if(Name_residuals != "ResidusPRED"){
    R2ecol=summary(model_ecol)$r.sq
  } else {
    R2ecol=as.numeric(piecewiseSEM::rsquared(model_ecol)$Conditional)
  }
  
  # Calculate Rsquared of the footprint model
  if(Best_model=="Linear"){ # Linear models
    R2footprint=1-(sum(residuals(model_lin)^2, na.rm=T) / sum((data$Res_SSE-mean(data$Res_SSE, na.rm=T))^2, na.rm=T))
  } 
  if(Best_model=="Threshold"){ # Threshold models
    R2footprint=1-(sum(residuals(model_thre)^2, na.rm=T) / sum((data$Res_SSE-mean(data$Res_SSE, na.rm=T))^2, na.rm=T))
  }
  
  return(as.numeric(c(
    round(model_thre$coefficients["hfp"],3),
    round(model_thre$coefficients["(hfp-chngpt)+"],3),
    round(model_thre$coefficients["chngpt"],2),
    as.data.frame(summary(model_thre)[1])[,5][2],
    as.data.frame(summary(model_thre)[1])[,5][3],
    as.data.frame(summary(model_thre)[1])[,2][2],
    as.data.frame(summary(model_thre)[1])[,2][3],
    as.data.frame(summary(model_thre)[2])[2,1],
    sum(data$SSE_thre, na.rm=T), 
    sum(data$SSE_lin, na.rm=T), 
    AIC(model_thre), 
    AIC(model_lin),
    R2ecol,
    R2footprint
  )))
}






##########################
### PREPARE EBIRD DATA ###
##########################

### Charge eBird data
chlist<-readRDS("Intermediate files/Data_RTRF_eBird.rds")

# Choose the HF to use (here for the main text, see end of the script for the SI)
chlist$hfp<-(chlist$hf13+chlist$hf09)/2



#####################################
### ANALYSIS 1 : SPECIES RICHNESS ###
#####################################

### eBird
source("eBird Richness script.R")




### PREDICTS
# Charge the data
dataPRED<-read.csv("PREDICTS/Data_RTRF_PREDICTS.csv")
dataPRED$year<-as.Date(dataPRED$Sample_start_earliest, "%d/%m/%Y") %>% format(., "%Y")
dataPRED<-subset(dataPRED, dataPRED$year >= 2000 & is.na(dataPRED$npp)==F) # Remove studies that started before 2000

# Choose HF value
dataPRED$hfp<-(dataPRED$hf05+dataPRED$hf09)/2

# Run the analysis
source("PREDICTS script.R")





### BBS
# Charge data
Routes<-read.csv("BBS/Data_RTRF_BBS.csv")
Routes<-Routes[order(Routes$RtCode),]

# Choose HF value
Routes$hfp<-(Routes$hf13+Routes$hf09)/2 

# Run the analysis
source("BBS script.R")





### Save plots 
source("Save plot script.R")
cowplot::save_plot("Figures/Main/Species.richness2.raw.svg", Gtot, base_width=7.5, base_height=8)



### Supplementary plot (smoothed terms)
Gtotsupp<-cowplot::plot_grid(
  ggplot(DF, aes(x=hfp))+
    geom_line(aes(y=Rich3), col="#7570b3ff", size=1.3)+
    scale_y_continuous(limits=c(0,1.2*max(DF$Rich3)))+
    xlab("Human footprint") + ylab("Species richness")+
    theme_bw(),
  
  ggplot(DFhtsp2)+
    geom_line(aes(x=hfp, y=Rich3), col="#7570b3ff", size=1.7)+
    scale_y_continuous(limits=c(0,1.2*max(DFhtsp2$Rich4)))+
    facet_wrap(~ htsp2, ncol=4)+
    xlab("Human footprint") + ylab("Species richness")+
    theme_minimal(),
  
  ncol=1, labels="AUTO")
cowplot::save_plot("Figures/Supplementary/Species.richnessSUPP.raw.pdf", Gtotsupp, base_width=7.5, base_height=6)











##########################################
### Analysis 2: ASSEMBLAGE COMPOSITION ###
##########################################


### Run the analysis
source("eBird Composition script.R")


### Save plots
cowplot::save_plot("Figures/Main/Species.replacement.raw.svg", GT, base_width=10, base_height=6.5)
cowplot::save_plot("Figures/Supplementary/Species.replacement.SUPP.raw.svg", GTsupp, base_width=10, base_height=6.5)


### Save Table

# Round values
for(COL in c("P1", "P2")){coef[,COL]<-round(coef[,COL], 3)}

# Transform p-values in stars
for(i in which(is.na(coef$P2)==F)){
  if(as.numeric(coef$P1[i])<0.05){coef$P1[i]<-stars.pval(as.numeric(coef$P1[i]))} else{coef$P1[i]<-round(as.numeric(coef$P1[i]),2)}
  if(as.numeric(coef$P2[i])<0.05){coef$P2[i]<-stars.pval(as.numeric(coef$P2[i]))} else{coef$P2[i]<-round(as.numeric(coef$P2[i]),2)}
}

# Save
write.csv(coef, "Figures/Main/Table1.csv", row.names=F) 






#######################################################################################################
### SUPPLEMENTARY ANALYSIS: Fig.2 with human footprint from 2009 and with human footprint from 2013 ###
#######################################################################################################


######################################################
### HUMAN FOOTPRINT 2009 (from Venter et al. 2016) ###
######################################################

rm(list=setdiff(ls(), c("save_coef", "ColNum"))) ; gc()
setwd("H:/RTRF_project")

# Create coef
coef<-data.frame(Index = c(rep("Overall", 11), "HighFor", "MedFor", "LowFor", "NoFor", "Endemic", "LargeRange", "thrNT", "Exotic", "SensiHigh", "SensiTol", "SensiAnt", rep("Overall", 8)), 
                 Zone = c("Global", "ATL", "AND", "TUM", "MES", "EAS", "GHA", "IND", "SUN", "PREDICTS", "BBS", rep("Global", 11), "ATLlin", "ANDlin", "TUMlin", "MESlin", "EASlin", "GHAlin", "INDlin", "SUNlin"))
coef$R2_footprint<-coef$R2_ecol<-coef$AIC_lin<-coef$AIC_thresh<-coef$SSE_linear<-coef$SSE_thresh<-coef$SE_thresh<-coef$SE_inter<-coef$SE_hft<-coef$P2<-coef$Coef2<-coef$Threshold<-coef$P1<-coef$Coef1<-NA


# Run richness for eBird (total and per hotspot)
chlist<-readRDS("Intermediate files/Data_RTRF_eBird.rds")
chlist$hfp<-chlist$hf09 # Choose the year of HF to use

source("eBird Richness script.R")


# Run richness for PREDICTS
dataPRED<-read.csv("PREDICTS/Data_RTRF_PREDICTS.csv")
dataPRED$year<-as.Date(dataPRED$Sample_start_earliest, "%d/%m/%Y") %>% format(., "%Y")
dataPRED<-subset(dataPRED, dataPRED$year >= 2000 & is.na(dataPRED$npp)==F) # Remove studies that started before 2000
dataPRED$hfp<-dataPRED$hf09 # Choose the year of HF to use

source("PREDICTS script.R")


# Run richness for BBS
Routes<-read.csv("BBS/Data_RTRF_BBS.csv")
Routes<-Routes[order(Routes$RtCode),]
Routes$hfp<-Routes$hf09 # Choose the year of HF to use

source("BBS script.R")



# Save plot
source("Save plot script.R")
cowplot::save_plot("Figures/Supplementary/Species.richness2.raw_HF09.svg", Gtot, base_width=7.5, base_height=8)

# Composition analysis
source("eBird Composition script.R")
cowplot::save_plot("Figures/Supplementary/Species.replacement.raw_HF09.svg", GT, base_width=10, base_height=6.5)

# Save table
for(COL in c("P1", "P2")){coef[,COL]<-round(coef[,COL], 3)}

for(i in which(is.na(coef$P2)==F)){
  if(as.numeric(coef$P1[i])<0.05){coef$P1[i]<-stars.pval(as.numeric(coef$P1[i]))}
  if(as.numeric(coef$P2[i])<0.05){coef$P2[i]<-stars.pval(as.numeric(coef$P2[i]))}
}
write.csv(coef, "Figures/Supplementary/Table1_HF09.csv", row.names=F) 






########################################################
### HUMAN FOOTPRINT 2013 (from Williams et al. 2020) ###
########################################################

rm(list=setdiff(ls(), c("save_coef", "ColNum"))) ; gc()
setwd("H:/RTRF_project")

# Create coef
coef<-data.frame(Index = c(rep("Overall", 11), "HighFor", "MedFor", "LowFor", "NoFor", "Endemic", "LargeRange", "thrNT", "Exotic", "SensiHigh", "SensiTol", "SensiAnt", rep("Overall", 8)), 
                 Zone = c("Global", "ATL", "AND", "TUM", "MES", "EAS", "GHA", "IND", "SUN", "PREDICTS", "BBS", rep("Global", 11), "ATLlin", "ANDlin", "TUMlin", "MESlin", "EASlin", "GHAlin", "INDlin", "SUNlin"))
coef$R2_footprint<-coef$R2_ecol<-coef$AIC_lin<-coef$AIC_thresh<-coef$SSE_linear<-coef$SSE_thresh<-coef$SE_thresh<-coef$SE_inter<-coef$SE_hft<-coef$P2<-coef$Coef2<-coef$Threshold<-coef$P1<-coef$Coef1<-NA






# Run richness for eBird (total and per hotspot)
chlist<-readRDS("Intermediate files/Data_RTRF_eBird.rds")
chlist$hfp<-chlist$hf13 # Choose the year of HF to use

source("eBird Richness script.R")


# Run richness for PREDICTS
dataPRED<-read.csv("PREDICTS/Data_RTRF_PREDICTS.csv")
dataPRED$year<-as.Date(dataPRED$Sample_start_earliest, "%d/%m/%Y") %>% format(., "%Y")
dataPRED<-subset(dataPRED, dataPRED$year >= 2000 & is.na(dataPRED$npp)==F) # Remove studies that started before 2000
dataPRED$hfp<-dataPRED$hf13 # Choose the year of HF to use

source("PREDICTS script.R")


# Run richness for BBS
Routes<-read.csv("BBS/Data_RTRF_BBS.csv")
Routes<-Routes[order(Routes$RtCode),]
Routes$hfp<-Routes$hf13 # Choose the year of HF to use

source("BBS script.R")


# Save plot
source("Save plot script.R")
cowplot::save_plot("Figures/Supplementary/Species.richness2.raw_HF13.svg", Gtot, base_width=7.5, base_height=8)

# Composition analysis
source("eBird Composition script.R")
cowplot::save_plot("Figures/Supplementary/Species.replacement.raw_HF13.svg", GT, base_width=10, base_height=6.5)

# Save table
for(COL in c("P1", "P2")){coef[,COL]<-round(coef[,COL], 3)}

for(i in which(is.na(coef$P2)==F)){
  if(as.numeric(coef$P1[i])<0.05){coef$P1[i]<-stars.pval(as.numeric(coef$P1[i]))}
  if(as.numeric(coef$P2[i])<0.05){coef$P2[i]<-stars.pval(as.numeric(coef$P2[i]))}
}
write.csv(coef, "Figures/Supplementary/Table1_HF13.csv", row.names=F) 


