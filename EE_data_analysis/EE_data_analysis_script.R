###############################################################################################################################################
###############################################################################################################################################
######################################Inhibition of soil extracelular enzyme activity - monophosphoesterases###################################
###############################################################################################################################################
###############################################################################################################################################

#############################################################Loading libraries#################################################################
library(readODS)
library(openxlsx)
library(ggplot2)
library(minpack.lm)
library(dplyr)
library(lamW)
library(bbmle)
library(nls2)
#############################################################GGPLOT THEME######################################################################
theme_min<-theme(axis.text.x=element_text(vjust=0.2, size=18, colour="black"),
                 axis.text.y=element_text(hjust=0.2, size=18, colour="black"),
                 axis.title=element_text(size=18, colour="black"),
                 axis.line=element_line(size=0.5, colour="black"),
                 strip.text=element_text(size=18, face="bold"),
                 axis.ticks=element_line(size=1, colour="black"),
                 axis.ticks.length=unit(-0.05, "cm"),
                 panel.background=element_rect(colour="black", fill="white"),
                 panel.grid=element_line(linetype=0),
                 legend.text=element_text(size=14, colour="black"),
                 legend.title=element_text(size=14, colour="black"),
                 legend.position=c("right"),
                 legend.key.size=unit(1, "cm"),
                 strip.background=element_rect(fill="grey98", colour="black"),
                 legend.key=element_rect(fill="white", size=1.2),
                 legend.spacing=unit(0.5, "cm"),
                 plot.title=element_text(size=18, face="bold", hjust=-0.05))
###############################################################################################################################################
#Not run
################################################################Raw data#######################################################################
# #loading the script to recalculate raw data
# source("../EE_data_analysis/E_calc.R")
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Plesne catchment~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# #Litter horizon 
# ##Measured on 14.4.2020
# ###uploading the file with reaction time
# po_time<-read_ods(path="../Raw_data//14.4.2020/casy14.4.2020.ods", sheet = 2, col_names=TRUE)
# po_all<-E_calc(dataset = c("/mnt/580CBE2464C5F83D/pracovni/data_statistika/Junior/Enzyme_inhibition/Raw_data/14.4.2020/MUB_14.4.2020.xlsx"),
#                  MUBconc = c(0, 5, 25, 50, 125, 250),
#                  APconc = c(50, 100, 250, 500, 1000),
#                  Iconc = c(0, 20, 80, 120, 240),
#                  Nmeasure = 25, empty = 1,
#                  Times = as.numeric(po_time[po_time$Ctr=="NE", "Time"]))
# ###Extract the data
# po<-po_all$data
# ###Do the correction against the control
# po$Pcorr<-po$P-pmax(0, po$Pcontr)
# ###Remove negative values
# po$Pcorr2<-pmax(0, as.numeric(po$Pcorr))
# 
# ###All units are umol/L
# ###recalculate to umol/g(DW) - fresh soil:buffer - 1:100 (w/w)
# po$Product<-po$Pcorr2/10/0.3195
# po$Substrate<-po$C_AP/10/0.3195
# po$Inhibitor<-po$C_I/10/0.3195
# 
# ###Add the sample description
# po$Catchment<-"Plesne"
# po$Horizon<-"Litter"
# 
# #Organic topsoil horizon
# ###Measured 24.3.2020
# ###uploading the file with reaction time
# pa_time<-read_ods(path="../Raw_data/24.3.2020/casy24.3.2020.ods", sheet = 2, col_names=TRUE)
# pa_all<-E_calc(dataset = c("/mnt/580CBE2464C5F83D/pracovni/data_statistika/Junior/Enzyme_inhibition/Raw_data/24.3.2020/MUB_24.3.2020.xlsx"),
#                  MUBconc = c(0, 5, 25, 50, 125, 250),
#                  APconc = c(50, 100, 250, 500, 1000),
#                  Iconc = c(0, 20, 80, 120, 240),
#                  Nmeasure = 25, empty = 1,
#                  Times = as.numeric(pa_time[pa_time$Ctr=="NE", "Time"]))
# ###Extract the data
# pa<-pa_all$data
# ###Do the correction against the control
# pa$Pcorr<-pa$P-pmax(0, pa$Pcontr)
# ###Remove negative values
# pa$Pcorr2<-pmax(0, as.numeric(pa$Pcorr))
# 
# ###All units are umol/L
# ###recalculate to umol/g(DW) - fresh soil:buffer - 1:100 (w/w)
# pa$Product<-pa$Pcorr2/10/0.2992
# pa$Substrate<-pa$C_AP/10/0.2992
# pa$Inhibitor<-pa$C_I/10/0.2992
# 
# pa$Catchment<-"Plesne"
# pa$Horizon<-"Organic topsoil"
# 
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Certovo catchment~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# #Litter horizon 
# ##Measured on 6.4.2020
# ###uploading the file with reaction time
# co_time<-read_ods(path="../Raw_data/6.4.2020/casy6.4.2020.ods", sheet = 2, col_names=TRUE)
# co_all<-E_calc(dataset = c("/mnt/580CBE2464C5F83D/pracovni/data_statistika/Junior/Enzyme_inhibition/Raw_data/6.4.2020/MUB_6.4.2020.xlsx"),
#                  MUBconc = c(0, 5, 25, 50, 125, 250),
#                  APconc = c(50, 100, 250, 500, 1000),
#                  Iconc = c(0, 20, 80, 120, 240),
#                  Nmeasure = 25, empty = 1,
#                  Times = as.numeric(co_time[co_time$Ctr=="NE", "Time"]))
# 
# ###Extract the data
# co<-co_all$data
# ###Do the correction against the control
# co$Pcorr<-co$P-pmax(0, co$Pcontr)
# ###Remove negative values
# co$Pcorr2<-pmax(0, as.numeric(co$Pcorr))
# 
# ###All units are umol/L
# ###recalculate to umol/g(DW) - fresh soil:buffer - 1:100 (w/w)
# co$Product<-co$Pcorr2/10/0.2911
# co$Substrate<-co$C_AP/10/0.2911
# co$Inhibitor<-co$C_I/10/0.2911
# 
# ###Add the sample description
# co$Catchment<-"Certovo"
# co$Horizon<-"Litter"
# 
# #Organic topsoil horizon 
# ##Measured on 30.3.2020
# ###uploading the file with reaction time
# ca_time<-read_ods(path="../Raw_data/30.3.2020/casy30.3.2020.ods", sheet = 2, col_names=TRUE)
# ca_all<-E_calc(dataset = c("/mnt/580CBE2464C5F83D/pracovni/data_statistika/Junior/Enzyme_inhibition/Raw_data/30.3.2020/MUB_30.3.2020.xlsx"),
#                  MUBconc = c(0, 5, 25, 50, 125, 250),
#                  APconc = c(50, 100, 250, 500, 1000),
#                  Iconc = c(0, 20, 80, 120, 240),
#                  Nmeasure = 25, empty = 1,
#                  Times = as.numeric(ca_time[ca_time$Ctr=="NE", "Time"]))
# ###Extract the data
# ca<-ca_all$data
# ###Do the correction against the control
# ca$Pcorr<-ca$P-pmax(0, ca$Pcontr)
# ###Remove negative values
# ca$Pcorr2<-pmax(0, as.numeric(ca$Pcorr))
# 
# ###All units are umol/L
# ###recalculate to umol/g(DW) - fresh soil:buffer - 1:100 (w/w)
# ca$Product<-ca$Pcorr2/10/0.3122
# ca$Substrate<-ca$C_AP/10/0.3122
# ca$Inhibitor<-ca$C_I/10/0.3122
# 
# ###Add the sample description
# ca$Catchment<-"Certovo"
# ca$Horizon<-"Organic topsoil"
# 
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Merge and export the data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# edata<-rbind(po,pa,co,ca)
# write.csv(edata, file = c("/mnt/580CBE2464C5F83D/pracovni/data_statistika/Junior/Enzyme_inhibition/enzyme_data.csv"))
# #Product, Substrate, Inhibitor - units umol/g(DW), Time - min
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Script starts here
##########################################################Reading all data############################################################
#Enzyme data
edata<-read.csv("../enzyme_data.csv")
#add lumped factor for graphics
edata$Legend<-NA
for(i in 1:nrow(edata)){
  if(edata$Catchment[i]=="Plesne" & edata$Horizon[i]=="Litter"){
    edata$Legend[i]<-"Plesne - Litter"
  }else{
    if(edata$Catchment[i]=="Plesne" & edata$Horizon[i]=="Organic topsoil"){
      edata$Legend[i]<-"Plesne - Organic topsoil"
    }else{
      if(edata$Catchment[i]=="Certovo" & edata$Horizon[i]=="Litter"){
        edata$Legend[i]<-"Certovo - Litter"
      }else{
        edata$Legend[i]<-"Certovo - Organic topsoil"
      }
    }
  }
}
#P sorption data
pdata<-read.csv("../sorption_data.csv")
pdata$Catchment<-factor(pdata$Catchment, levels=c("Plesne", "Certovo"))
#add lumped factor for graphics
pdata$Legend<-NA
for(i in 1:nrow(pdata)){
  if(pdata$Catchment[i]=="Plesne" & pdata$Horizon[i]=="Litter"){
    pdata$Legend[i]<-"Plesne - Litter"
  }else{
    if(pdata$Catchment[i]=="Plesne" & pdata$Horizon[i]=="Organic topsoil"){
      pdata$Legend[i]<-"Plesne - Organic topsoil"
    }else{
      if(pdata$Catchment[i]=="Certovo" & pdata$Horizon[i]=="Litter"){
        pdata$Legend[i]<-"Certovo - Litter"
      }else{
        pdata$Legend[i]<-"Certovo - Organic topsoil"
      }
    }
  }
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Initial product (P-PO4 measured as SRP) concentration~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Measured vs added SRP
ggplot(pdata[pdata$Time==0, ], aes(SRP_a, SRP_o)) + geom_point(cex=6, aes(fill=Legend, shape=Legend)) +
  theme_min + scale_fill_manual(values = c("black", "black", "grey", "grey")) + 
  scale_shape_manual(values = c(21, 22, 21, 22)) + 
  theme(legend.title = element_blank(), legend.position = c(0.1, 0.8)) +
  ylim(0, 20) + xlim(0, 20) + geom_abline(intercept = 0, slope=1, color="black") +
  stat_smooth(method="lm", lty=2, lwd=1.5, se=F, colour="black") +
  ylab(expression(paste("Measured SRP (", mu, "mol ", g^{-1}, ")"))) + 
  xlab(expression(paste("Added SRP (", mu, "mol ", g^{-1}, ")")))

#Define 4 different relationships between the added and measured initial SRP concentrations
lm_po<-lm(SRP_o~SRP_a, pdata, subset = Catchment=="Plesne" & Horizon=="Litter")
lm_pa<-lm(SRP_o~SRP_a, pdata, subset = Catchment=="Plesne" & Horizon=="Organic topsoil")
lm_co<-lm(SRP_o~SRP_a, pdata, subset = Catchment=="Certovo" & Horizon=="Litter")
lm_ca<-lm(SRP_o~SRP_a, pdata, subset = Catchment=="Certovo" & Horizon=="Organic topsoil")

summary(lm_po)
summary(lm_pa)
summary(lm_co)
summary(lm_ca)

#Calculating the initial concentration of SRP/Inhibitor
edata$SRP_a<-edata$Inhibitor
edata$InhibitorSRP<-NA
edata[(edata$Catchment=="Plesne" & edata$Horizon=="Litter"), "InhibitorSRP"]<-
  predict(lm_po, newdata = as.data.frame(edata[(edata$Catchment=="Plesne" & edata$Horizon=="Litter"), ]))
edata[(edata$Catchment=="Plesne" & edata$Horizon=="Organic topsoil"), "InhibitorSRP"]<-
  predict(lm_pa, newdata = as.data.frame(edata[(edata$Catchment=="Plesne" & edata$Horizon=="Organic topsoil"), ]))
edata[(edata$Catchment=="Certovo" & edata$Horizon=="Litter"), "InhibitorSRP"]<-
  predict(lm_co, newdata = as.data.frame(edata[(edata$Catchment=="Plesne" & edata$Horizon=="Litter"), ]))
edata[(edata$Catchment=="Certovo" & edata$Horizon=="Organic topsoil"), "InhibitorSRP"]<-
  predict(lm_ca, newdata = as.data.frame(edata[(edata$Catchment=="Plesne" & edata$Horizon=="Organic topsoil"), ]))

#############################################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Testing the product inhibition~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#############################################################################################################################################
#Defining the non-linear functions whose fit will be compared
##Without product inhibition (wi)
wi<-"time ~ -1/Vmax*(Km*log((Substrate-Product)/(Substrate))-Product)"
##Competitive inhibition (ci)
ci<-"time ~ -1/Vmax*(Km*((Substrate+InhibitorSRP)/Kic+1)*log((Substrate-Product)/(Substrate))+
                      (1-Km/Kic)*-Product)"
##Non competitive inhibition (nci)
nci<-"time ~ -1/Vmax*(Km*((Substrate+InhibitorSRP)/Kiu+1)*log((Substrate-Product)/(Substrate))+
                      (1-(Km+Substrate+InhibitorSRP)/Kiu)*-Product-((Substrate-Product)^2-(Substrate)^2)/2/Kiu)"
##Uncompetitive inhibition (uci)
uci<-"time ~ -1/Vmax*(Km*log((Substrate-Product)/(Substrate))+
                      (1+(Km+Substrate+InhibitorSRP)/Kiu)*-Product-((Substrate-Product)^2-(Substrate)^2)/2/Kiu)"

#Plesne catchment
##Litter horizon
plo_wi<-nls(wi, data = edata, subset = Catchment=="Plesne" & Horizon=="Litter", start=list(Vmax=1, Km=1))
plo_ci<-nls(ci, data = edata, subset = Catchment=="Plesne" & Horizon=="Litter", start=list(Vmax=1, Km=1, Kic=1))
plo_nci<-nls(nci, data = edata, subset = Catchment=="Plesne" & Horizon=="Litter", start=list(Vmax=0.1, Km=10, Kiu=10))
plo_uci<-nls(uci, data = edata, subset = Catchment=="Plesne" & Horizon=="Litter", start=list(Vmax=0.1, Km=10, Kiu=20))
###AIC test
AICtab(plo_wi, plo_ci, plo_nci, plo_uci, weights=T, sort=T, base=T, logLik=T)
###Ftest
anova(plo_wi, plo_ci)
anova(plo_wi, plo_nci)
anova(plo_wi, plo_uci)

###Parameters
coef(plo_wi)
coef(plo_ci)
coef(plo_nci)
coef(plo_uci)

##Organic topsoil horizon
pla_wi<-nls(wi, data = edata, subset = Catchment=="Plesne" & Horizon=="Organic topsoil", start=list(Vmax=1, Km=1))
pla_ci<-nls(ci, data = edata, subset = Catchment=="Plesne" & Horizon=="Organic topsoil", start=list(Vmax=1, Km=1, Kic=1))
pla_nci<-nls(nci, data = edata, subset = Catchment=="Plesne" & Horizon=="Organic topsoil", start=list(Vmax=0.1, Km=10, Kiu=10))
pla_uci<-nlsLM(nci, data = edata, subset = Catchment=="Plesne" & Horizon=="Organic topsoil", start=list(Vmax=0.1, Km=10, Kiu=10))
###AIC test
AICtab(pla_wi, pla_ci, pla_nci, pla_uci, weights=T, sort=T, base=T, logLik=T)
###Ftest
anova(pla_wi, pla_ci)
anova(pla_wi, pla_nci)
anova(pla_wi, pla_uci)

###Parameters
coef(pla_wi)
coef(pla_ci)
coef(pla_nci)
coef(pla_uci)

#Certovo catchment
##Litter horizon
co_wi<-nls(wi, data = edata, subset = Catchment=="Certovo" & Horizon=="Litter", start=list(Vmax=1, Km=1))
co_ci<-nls(ci, data = edata, subset = Catchment=="Certovo" & Horizon=="Litter", start=list(Vmax=1, Km=1, Kic=1))
co_nci<-nls(nci, data = edata, subset = Catchment=="Certovo" & Horizon=="Litter", start=list(Vmax=0.1, Km=10, Kiu=10))
co_uci<-nls(uci, data = edata, subset = Catchment=="Certovo" & Horizon=="Litter", start=list(Vmax=0.1, Km=10, Kiu=20))
###AIC test
AICtab(co_wi, co_ci, co_nci, co_uci, weights=T, sort=T, base=T, logLik=T)
###Ftest
anova(co_wi, co_ci)
anova(co_wi, co_nci)
anova(co_wi, co_uci)

###Parameters
coef(co_wi)
coef(co_ci)
coef(co_nci)
coef(co_uci)

##Organic topsoil horizon
ca_wi<-nls(wi, data = edata, subset = Catchment=="Certovo" & Horizon=="Organic topsoil", start=list(Vmax=1, Km=1))
ca_ci<-nls(ci, data = edata, subset = Catchment=="Certovo" & Horizon=="Organic topsoil", start=list(Vmax=1, Km=1, Kic=1))
ca_nci<-nls(nci, data = edata, subset = Catchment=="Certovo" & Horizon=="Organic topsoil", start=list(Vmax=0.1, Km=10, Kiu=10))
ca_uci<-nlsLM(nci, data = edata, subset = Catchment=="Certovo" & Horizon=="Organic topsoil", start=list(Vmax=0.1, Km=10, Kiu=10))
###AIC test
AICtab(ca_wi, ca_ci, ca_nci, ca_uci, weights=T, sort=T, base=T, logLik=T)
###Ftest
anova(ca_wi, ca_ci)
anova(ca_wi, ca_nci)
anova(ca_wi, ca_uci)

###Parameters
coef(ca_wi)
coef(ca_ci)
coef(ca_nci)
coef(ca_uci)

#Calculate the rates as a slope of linear regression
##For first 5 minutes
###sem budu ukladat data
erates<-data.frame(Substrate=c(rep(unique(edata[(edata$Catchment=="Plesne" & edata$Horizon=="Litter"), "Substrate"]), each=5),
                               rep(unique(edata[(edata$Catchment=="Plesne" & edata$Horizon=="Organic topsoil"), "Substrate"]), each=5),
                               rep(unique(edata[(edata$Catchment=="Certovo" & edata$Horizon=="Litter"), "Substrate"]), each=5),
                               rep(unique(edata[(edata$Catchment=="Certovo" & edata$Horizon=="Organic topsoil"), "Substrate"]), each=5)),
                   InhibitorSRP=c(rep(unique(edata[(edata$Catchment=="Plesne" & edata$Horizon=="Litter"), "InhibitorSRP"]), times=5),
                                  rep(unique(edata[(edata$Catchment=="Plesne" & edata$Horizon=="Organic topsoil"), "InhibitorSRP"]), times=5),
                                  rep(unique(edata[(edata$Catchment=="Certovo" & edata$Horizon=="Litter"), "InhibitorSRP"]), times=5),
                                  rep(unique(edata[(edata$Catchment=="Certovo" & edata$Horizon=="Organic topsoil"), "InhibitorSRP"]), times=5)),
                   Catchment=c(rep("Plesne", times=50), rep("Certovo", times=50)),
                   Horizon=c(rep("Litter", times=25), rep("Organic topsoil", times=25),
                             rep("Litter", times=25), rep("Organic topsoil", times=25)))
erates$v<-NA
erates$v.se<-NA

# for(i in unique(ed3$C_AP)){
#   for(n in unique(ed3$C_I)){
#     lmv<-lm(Pcorr2~time-1, ed3[(ed3$C_AP==i & ed3$C_I==n & ed3$time<5 & ed3$outlier=="NO"), ])
#     ed3_v[(ed3_v$C_AP==i & ed3_v$C_I==n), "v"]<-summary(lmv)$coefficients[1]
#     ed3_v[(ed3_v$C_AP==i & ed3_v$C_I==n), "v.se"]<-summary(lmv)$coefficients[2]
#   }
# }