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
library(FME)
library(DEoptim)
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
##Target data frame
erates<-data.frame(Substrate=c(rep(unique(edata[(edata$Catchment=="Plesne" & edata$Horizon=="Litter"), "Substrate"]), each=5),
                               rep(unique(edata[(edata$Catchment=="Plesne" & edata$Horizon=="Organic topsoil"), "Substrate"]), each=5),
                               rep(unique(edata[(edata$Catchment=="Certovo" & edata$Horizon=="Litter"), "Substrate"]), each=5),
                               rep(unique(edata[(edata$Catchment=="Certovo" & edata$Horizon=="Organic topsoil"), "Substrate"]), each=5)),
                   InhibitorSRP=c(rep(unique(edata[(edata$Catchment=="Plesne" & edata$Horizon=="Litter"), "InhibitorSRP"]), times=5),
                                  rep(unique(edata[(edata$Catchment=="Plesne" & edata$Horizon=="Organic topsoil"), "InhibitorSRP"]), times=5),
                                  rep(unique(edata[(edata$Catchment=="Certovo" & edata$Horizon=="Litter"), "InhibitorSRP"]), times=5),
                                  rep(unique(edata[(edata$Catchment=="Certovo" & edata$Horizon=="Organic topsoil"), "InhibitorSRP"]), times=5)),
                   Inhibitor=c(rep(unique(edata[(edata$Catchment=="Plesne" & edata$Horizon=="Litter"), "Inhibitor"]), times=5),
                                  rep(unique(edata[(edata$Catchment=="Plesne" & edata$Horizon=="Organic topsoil"), "Inhibitor"]), times=5),
                                  rep(unique(edata[(edata$Catchment=="Certovo" & edata$Horizon=="Litter"), "Inhibitor"]), times=5),
                                  rep(unique(edata[(edata$Catchment=="Certovo" & edata$Horizon=="Organic topsoil"), "Inhibitor"]), times=5)),
                   Catchment=c(rep("Plesne", times=50), rep("Certovo", times=50)),
                   Horizon=c(rep("Litter", times=25), rep("Organic topsoil", times=25),
                             rep("Litter", times=25), rep("Organic topsoil", times=25)))
##For first 5 minutes
erates$v5<-NA
erates$v5.se<-NA

for(l in unique(edata$Horizon)){
  for(k in unique(edata$Catchment)){
    for(i in unique(edata[(edata$Horizon==l & edata$Catchment==k), "Substrate"])){
      for(n in unique(edata[(edata$Horizon==l & edata$Catchment==k), "InhibitorSRP"])){
        lmv<-lm(Product~time-1, edata[(edata$Substrate==i & edata$InhibitorSRP==n & 
                                         edata$time<5 & edata$Horizon==l & edata$Catchment==k), ])
        erates[(erates$Substrate==i & erates$InhibitorSRP==n & 
                  erates$Horizon==l & erates$Catchment==k), "v5"]<-summary(lmv)$coefficients[1]
        erates[(erates$Substrate==i & erates$InhibitorSRP==n & 
                  erates$Horizon==l & erates$Catchment==k), "v5.se"]<-summary(lmv)$coefficients[2]
      }
    }
  }
}

# ##For first 10 minutes
# erates$v10<-NA
# erates$v10.se<-NA
# 
# for(l in unique(edata$Horizon)){
#   for(k in unique(edata$Catchment)){
#     for(i in unique(edata[(edata$Horizon==l & edata$Catchment==k), "Substrate"])){
#       for(n in unique(edata[(edata$Horizon==l & edata$Catchment==k), "InhibitorSRP"])){
#         lmv<-lm(Product~time-1, edata[(edata$Substrate==i & edata$InhibitorSRP==n & 
#                                          edata$time<10 & edata$Horizon==l & edata$Catchment==k), ])
#         erates[(erates$Substrate==i & erates$InhibitorSRP==n & 
#                   erates$Horizon==l & erates$Catchment==k), "v10"]<-summary(lmv)$coefficients[1]
#         erates[(erates$Substrate==i & erates$InhibitorSRP==n & 
#                   erates$Horizon==l & erates$Catchment==k), "v10.se"]<-summary(lmv)$coefficients[2]
#       }
#     }
#   }
# }
# 
# ##For first 20 minutes
# erates$v20<-NA
# erates$v20.se<-NA
# 
# for(l in unique(edata$Horizon)){
#   for(k in unique(edata$Catchment)){
#     for(i in unique(edata[(edata$Horizon==l & edata$Catchment==k), "Substrate"])){
#       for(n in unique(edata[(edata$Horizon==l & edata$Catchment==k), "InhibitorSRP"])){
#         lmv<-lm(Product~time-1, edata[(edata$Substrate==i & edata$InhibitorSRP==n & 
#                                          edata$time<20 & edata$Horizon==l & edata$Catchment==k), ])
#         erates[(erates$Substrate==i & erates$InhibitorSRP==n & 
#                   erates$Horizon==l & erates$Catchment==k), "v20"]<-summary(lmv)$coefficients[1]
#         erates[(erates$Substrate==i & erates$InhibitorSRP==n & 
#                   erates$Horizon==l & erates$Catchment==k), "v20.se"]<-summary(lmv)$coefficients[2]
#       }
#     }
#   }
# }

##For 2 hours
erates$v120<-NA
erates$v120.se<-NA

for(l in unique(edata$Horizon)){
  for(k in unique(edata$Catchment)){
    for(i in unique(edata[(edata$Horizon==l & edata$Catchment==k), "Substrate"])){
      for(n in unique(edata[(edata$Horizon==l & edata$Catchment==k), "InhibitorSRP"])){
        lmv<-lm(Product~time-1, edata[(edata$Substrate==i & edata$InhibitorSRP==n & 
                                         edata$time>30 & edata$time<120 & edata$Horizon==l & edata$Catchment==k), ])
        erates[(erates$Substrate==i & erates$InhibitorSRP==n & 
                  erates$Horizon==l & erates$Catchment==k), "v120"]<-summary(lmv)$coefficients[1]
        erates[(erates$Substrate==i & erates$InhibitorSRP==n & 
                  erates$Horizon==l & erates$Catchment==k), "v120.se"]<-summary(lmv)$coefficients[2]
      }
    }
  }
}

##Calculate the relative activity (in respect to substrate without the inhibitor)
erates$v5rel<-NA
erates$v5rel.se<-NA

for(l in unique(edata$Horizon)){
 for(k in unique(edata$Catchment)){
   for(i in unique(edata[(edata$Horizon==l & edata$Catchment==k), "Substrate"])){
     for(n in unique(edata[(edata$Horizon==l & edata$Catchment==k), "Inhibitor"])){
       erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==n), "v5rel"]<-
         erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==n), "v5"]/
         erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==0), "v5"]*100
       erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==n), "v5rel.se"]<-
         sqrt((erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==n), "v5.se"]/
         erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==n), "v5"])^2 + 
           (erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==0), "v5.se"]/
              erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==0), "v5"])^2)* 
         erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==n), "v5rel"]
     }
   }
 }
}

# erates$v10rel<-NA
# erates$v10rel.se<-NA
# 
# for(l in unique(edata$Horizon)){
#   for(k in unique(edata$Catchment)){
#     for(i in unique(edata[(edata$Horizon==l & edata$Catchment==k), "Substrate"])){
#       for(n in unique(edata[(edata$Horizon==l & edata$Catchment==k), "Inhibitor"])){
#         erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==n), "v10rel"]<-
#           erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==n), "v10"]/
#           erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==0), "v10"]*100
#         erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==n), "v10rel.se"]<-
#           sqrt((erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==n), "v10.se"]/
#                   erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==n), "v10"])^2 + 
#                  (erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==0), "v10.se"]/
#                     erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==0), "v10"])^2)* 
#           erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==n), "v10rel"]
#       }
#     }
#   }
# }
# 
# erates$v20rel<-NA
# erates$v20rel.se<-NA
# 
# for(l in unique(edata$Horizon)){
#   for(k in unique(edata$Catchment)){
#     for(i in unique(edata[(edata$Horizon==l & edata$Catchment==k), "Substrate"])){
#       for(n in unique(edata[(edata$Horizon==l & edata$Catchment==k), "Inhibitor"])){
#         erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==n), "v20rel"]<-
#           erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==n), "v20"]/
#           erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==0), "v20"]*100
#         erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==n), "v20rel.se"]<-
#           sqrt((erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==n), "v20.se"]/
#                   erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==n), "v20"])^2 + 
#                  (erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==0), "v20.se"]/
#                     erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==0), "v20"])^2)* 
#           erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==n), "v20rel"]
#       }
#     }
#   }
# }

erates$v120rel<-NA
erates$v120rel.se<-NA

for(l in unique(erates$Horizon)){
  for(k in unique(erates$Catchment)){
    for(i in unique(erates[(erates$Horizon==l & erates$Catchment==k), "Substrate"])){
      for(n in unique(erates[(erates$Horizon==l & erates$Catchment==k), "Inhibitor"])){
        erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==n), "v120rel"]<-
          erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==n), "v120"]/
          erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==0), "v120"]*100
        erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==n), "v120rel.se"]<-
          sqrt((erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==n), "v120.se"]/
                  erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==n), "v120"])^2 + 
                 (erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==0), "v120.se"]/
                    erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==0), "v120"])^2)* 
          erates[(erates$Horizon==l & erates$Catchment==k & erates$Substrate==i & erates$Inhibitor==n), "v120rel"]
      }
    }
  }
}

##Plot
###add lumped factor for graphics
erates$Legend<-NA
for(i in 1:nrow(erates)){
  if(erates$Catchment[i]=="Plesne" & erates$Horizon[i]=="Litter"){
    erates$Legend[i]<-"Plesne - Litter"
  }else{
    if(erates$Catchment[i]=="Plesne" & erates$Horizon[i]=="Organic topsoil"){
      erates$Legend[i]<-"Plesne - Organic topsoil"
    }else{
      if(erates$Catchment[i]=="Certovo" & erates$Horizon[i]=="Litter"){
        erates$Legend[i]<-"Certovo - Litter"
      }else{
        erates$Legend[i]<-"Certovo - Organic topsoil"
      }
    }
  }
}

erates$Substrate2<-round(erates$Substrate, 0)
erates$Substrate3<-rep(erates$Substrate2[1:25], times=4)

###Predictions of the product inhibition function
epredts<-data.frame(Substrate=c(rep(unique(edata[(edata$Catchment=="Plesne" & edata$Horizon=="Litter"), "Substrate"]), each=50),
                               rep(unique(edata[(edata$Catchment=="Plesne" & edata$Horizon=="Organic topsoil"), "Substrate"]), each=50),
                               rep(unique(edata[(edata$Catchment=="Certovo" & edata$Horizon=="Litter"), "Substrate"]), each=50),
                               rep(unique(edata[(edata$Catchment=="Certovo" & edata$Horizon=="Organic topsoil"), "Substrate"]), each=50)),
                   InhibitorSRP=rep(seq(0, 16, length.out = 50), times=4),
                   Catchment=c(rep("Plesne", times=50*5*2), rep("Certovo", times=50*5*2)),
                   Horizon=c(rep("Litter", times=50*5), rep("Organic topsoil", times=50*5),
                             rep("Litter", times=50*5), rep("Organic topsoil", times=50*5)))


epredts$v<-NA
for(i in 1:nrow(epredts)){
  if(epredts$Catchment[i]=="Plesne" & epredts$Horizon[i]=="Litter"){
    epredts$v[i]<-coef(plo_ci)[1]*epredts$Substrate[i]/(coef(plo_ci)[2]*(1+epredts$InhibitorSRP[i]/coef(plo_ci)[3]) + epredts$Substrate[i])
  }else{
    if(epredts$Catchment[i]=="Plesne" & epredts$Horizon[i]=="Organic topsoil"){
      epredts$v[i]<-coef(pla_ci)[1]*epredts$Substrate[i]/(coef(pla_ci)[2]*(1+epredts$InhibitorSRP[i]/coef(pla_ci)[3]) + epredts$Substrate[i])
    }else{
      if(epredts$Catchment[i]=="Certovo" & epredts$Horizon[i]=="Litter"){
        epredts$v[i]<-coef(co_ci)[1]*epredts$Substrate[i]/(coef(co_ci)[2]*(1+epredts$InhibitorSRP[i]/coef(co_ci)[3]) + epredts$Substrate[i])
      }else{
        epredts$v[i]<-coef(ca_ci)[1]*epredts$Substrate[i]/(coef(ca_ci)[2]*(1+epredts$InhibitorSRP[i]/coef(ca_ci)[3]) + epredts$Substrate[i])
      }
    }
  }
}

###relative rates
for(l in unique(epredts$Horizon)){
  for(k in unique(epredts$Catchment)){
    for(i in unique(epredts[(epredts$Horizon==l & epredts$Catchment==k), "Substrate"])){
      for(n in unique(epredts[(epredts$Horizon==l & epredts$Catchment==k), "InhibitorSRP"])){
        epredts[(epredts$Horizon==l & epredts$Catchment==k & epredts$Substrate==i & epredts$InhibitorSRP==n), "vrel"]<-
          epredts[(epredts$Horizon==l & epredts$Catchment==k & epredts$Substrate==i & epredts$InhibitorSRP==n), "v"]/
          epredts[(epredts$Horizon==l & epredts$Catchment==k & epredts$Substrate==i & epredts$InhibitorSRP==0), "v"]*100
      }
    }
  }
}

###add lumped factor for graphics
epredts$Legend<-NA
for(i in 1:nrow(epredts)){
  if(epredts$Catchment[i]=="Plesne" & epredts$Horizon[i]=="Litter"){
    epredts$Legend[i]<-"Plesne - Litter"
  }else{
    if(epredts$Catchment[i]=="Plesne" & epredts$Horizon[i]=="Organic topsoil"){
      epredts$Legend[i]<-"Plesne - Organic topsoil"
    }else{
      if(epredts$Catchment[i]=="Certovo" & epredts$Horizon[i]=="Litter"){
        epredts$Legend[i]<-"Certovo - Litter"
      }else{
        epredts$Legend[i]<-"Certovo - Organic topsoil"
      }
    }
  }
}

epredts$Substrate2<-round(epredts$Substrate, 0)
epredts$Substrate3<-rep(epredts$Substrate2[1:250], times=4)


ggplot(erates, aes(InhibitorSRP, v5rel)) + geom_point(cex=6, pch=21, aes(fill = as.factor(Substrate3))) +
  facet_grid(.~Legend) + theme_min + scale_fill_manual(values=c("black", "grey30", "grey60", "grey90", "white")) +
  scale_y_continuous(limits = c(0, 110), breaks = c(0, 20, 40, 60, 80, 100)) +
  xlim(0, 16) + geom_errorbar(aes(ymin=v5rel-v5rel.se, ymax=v5rel+v5rel.se), width=0.01) +
  geom_line(data=epredts, aes(InhibitorSRP, vrel, color=as.factor(Substrate3), linetype=as.factor(Substrate3)), lwd=1.2) +
  scale_color_manual(values=c("black", "grey30", "grey60", "grey90", "grey90")) +
  scale_linetype_manual(values=c(rep("solid", 4), "longdash"))

# ggplot(erates, aes(InhibitorSRP, v120rel)) + geom_point(cex=6, pch=21, aes(fill = as.factor(Substrate3))) +
# facet_grid(.~Legend) + theme_min + scale_fill_manual(values=c("black", "grey30", "grey60", "grey90", "white")) +
# scale_y_continuous(limits = c(0, 110), breaks = c(0, 20, 40, 60, 80, 100)) +
# xlim(0, 16) + geom_errorbar(aes(ymin=v120rel-v120rel.se, ymax=v120rel+v120rel.se), width=0.01) +
# geom_line(data=epredts, aes(InhibitorSRP, vrel, color=as.factor(Substrate3), linetype=as.factor(Substrate3)), lwd=1.2) +
# scale_color_manual(values=c("black", "grey30", "grey60", "grey90", "grey90")) +
# scale_linetype_manual(values=c(rep("solid", 4), "longdash"))

#############################################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Acknkowledging the SRP sorption/release kinetic~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#############################################################################################################################################
pdata$SRP_a2<-round(pdata$SRP_a, 0)
#Calculate relative change of SRP concentration in respects to time 0
##Experiment1
pdata_e1<-subset(pdata, Experiment==1)
pdata_e1$SRP_orel<-NA
for(i in unique(pdata_e1$Catchment)){
  for(n in unique(pdata_e1$Horizon)){
    for(k in unique(pdata_e1[(pdata_e1$Catchment==i & pdata_e1$Horizon==n), "SRP_a"])){
      for(l in unique(pdata_e1[(pdata_e1$Catchment==i & pdata_e1$Horizon==n & pdata_e1$SRP_a==k), "Time"])){
        pdata_e1[(pdata_e1$Catchment==i & pdata_e1$Horizon==n & pdata_e1$SRP_a==k & pdata_e1$Time==l), "SRP_orel"]<-
          pdata_e1[(pdata_e1$Catchment==i & pdata_e1$Horizon==n & pdata_e1$SRP_a==k & pdata_e1$Time==l), "SRP_o"]/
          pdata_e1[(pdata_e1$Catchment==i & pdata_e1$Horizon==n & pdata_e1$SRP_a==k & pdata_e1$Time==0), "SRP_o"]*100
      }
    }
  }
}

##Experiment2
pdata_e2<-subset(pdata, Experiment==2)
pdata_e2$SRP_orel<-NA
for(i in unique(pdata_e2$Catchment)){
  for(n in unique(pdata_e2$Horizon)){
    for(k in unique(pdata_e2[(pdata_e2$Catchment==i & pdata_e2$Horizon==n), "SRP_a"])){
      for(l in unique(pdata_e2[(pdata_e2$Catchment==i & pdata_e2$Horizon==n & pdata_e2$SRP_a==k), "Time"])){
        pdata_e2[(pdata_e2$Catchment==i & pdata_e2$Horizon==n & pdata_e2$SRP_a==k & pdata_e2$Time==l), "SRP_orel"]<-
          pdata_e2[(pdata_e2$Catchment==i & pdata_e2$Horizon==n & pdata_e2$SRP_a==k & pdata_e2$Time==l), "SRP_o"]/
          pdata_e2[(pdata_e2$Catchment==i & pdata_e2$Horizon==n & pdata_e2$SRP_a==k & pdata_e2$Time==0), "SRP_o"]*100
      }
    }
  }
}

#Visualize the results
ggplot(rbind(pdata_e1, pdata_e2), aes(Time, SRP_orel)) + geom_point(cex=6, pch=21, aes(fill = as.factor(SRP_a2))) +
  theme_min + #scale_fill_manual(values=c("black", "grey30", "grey60", "grey90", "white")) +
  stat_smooth(method = lm, se=F, aes(colour=as.factor(SRP_a2)))

ggplot(pdata, aes(Time, log(SRP_o))) + geom_point(cex=6, pch=21, aes(fill = as.factor(SRP_a2))) +
  facet_grid(.~Legend) + theme_min + #scale_fill_manual(values=c("black", "grey30", "grey60", "grey90", "white")) +
  stat_smooth(method = lm, se=F, aes(colour=as.factor(SRP_a2)))

#Mark the outliers
pdata$outlier<-c("NO")
pdata[(pdata$Catchment=="Certovo" & pdata$Horizon=="Litter" & pdata$SRP_a2==0 & pdata$Time==0 & log(pdata$SRP_o)>1), "outlier"]<-c("YES")
pdata[(pdata$Catchment=="Certovo" & pdata$Horizon=="Organic topsoil" & pdata$SRP_a2==0 & log(pdata$SRP_o)<0.4), "outlier"]<-c("YES")
pdata[(pdata$Catchment=="Certovo" & pdata$Horizon=="Organic topsoil" & pdata$SRP_a2==1 & log(pdata$SRP_o)<0.4), "outlier"]<-c("YES")
pdata[(pdata$Catchment=="Plesne" & pdata$Horizon=="Organic topsoil" & pdata$SRP_a2==0 & pdata$Experiment==2), "outlier"]<-c("YES")


#Calculate slopes of change of SRP concentration over time for each initial concentration of inhibitor, catchment and horizon
##Initialize the dataframe storing the results
SRPslopes<-data_frame(Catchment=character(), Horizon=character(), SRP_a=numeric(), slope=numeric(), se=numeric())

for(i in unique(pdata$Catchment)){
  for(n in unique(pdata$Horizon)){
    for(k in unique(pdata[(pdata$Catchment==i & pdata$Horizon==n), "SRP_a"])){
      tryCatch({
        SRPslopes<-rbind(SRPslopes,
                         data.frame(Catchment = pdata[(pdata$Catchment==i & pdata$Horizon==n & pdata$SRP_a==k), "Catchment"][1],
                                    Horizon = pdata[(pdata$Catchment==i & pdata$Horizon==n & pdata$SRP_a==k), "Horizon"][1],
                                    SRP_a=k,
                                    slope=coef(lm(SRP_o~Time, data=pdata[(pdata$Catchment==i & pdata$Horizon==n & pdata$SRP_a==k & pdata$outlier=="NO"),]))[2],
                                    se=summary(lm(SRP_o~Time, data=pdata[(pdata$Catchment==i & pdata$Horizon==n & pdata$SRP_a==k & pdata$outlier=="NO"),]))[[4]][[4]]
                         ))
      }, error = function(e){print("Nejde kurva")})
      
    }
  }
}
##Visualize
###add lumped factor for graphics
SRPslopes$Legend<-NA
for(i in 1:nrow(SRPslopes)){
  if(SRPslopes$Catchment[i]=="Plesne" & SRPslopes$Horizon[i]=="Litter"){
    SRPslopes$Legend[i]<-"Plesne - Litter"
  }else{
    if(SRPslopes$Catchment[i]=="Plesne" & SRPslopes$Horizon[i]=="Organic topsoil"){
      SRPslopes$Legend[i]<-"Plesne - Organic topsoil"
    }else{
      if(SRPslopes$Catchment[i]=="Certovo" & SRPslopes$Horizon[i]=="Litter"){
        SRPslopes$Legend[i]<-"Certovo - Litter"
      }else{
        SRPslopes$Legend[i]<-"Certovo - Organic topsoil"
      }
    }
  }
}
ggplot(SRPslopes, aes(SRP_a, slope)) + geom_point(cex=6, aes(fill=Legend, shape=Legend)) + theme_min +
  scale_fill_manual(values = c("black", "black", "grey", "grey")) + 
  scale_shape_manual(values = c(21, 22, 21, 22)) + geom_errorbar(aes(ymin=slope-se, ymax=slope+se))+
  theme(legend.title = element_blank(), legend.position = c(0.8, 0.8)) +
  stat_smooth(method=lm, lwd=1.2, color="black", lty=2) +
  ylab(expression(paste("SRP relative rate of increase (", min^{-1}, ")"))) + 
  xlab(expression(paste("Added SRP (", mu, "mol ", g^{-1}, ")"))) 


#Add to original data frame
##Generate the general relationship
SRPslopes$Inhibitor<-SRPslopes$SRP_a
slope_calc<-lm(slope~Inhibitor, SRPslopes)
summary(slope_calc)
##Apply to original data
edata$slope<-predict(slope_calc, newdata=edata)

#Run the competitive inhibition model again with SRP release over time accounted for
##Load the function
source("CI_SRP.R")
##Plesne catchment
###Litter horizon
plo_cisrp<-CI_SRP(data=subset(edata, Catchment=="Plesne" & Horizon=="Litter"))
plo_cisrp$Parameters
plo_cisrp$Goodness$Gfit

###Organic topsoil
pla_cisrp<-CI_SRP(data=subset(edata, Catchment=="Plesne" & Horizon=="Organic topsoil"))
pla_cisrp$Parameters
pla_cisrp$Goodness$Gfit

##Certovo catchment
###Litter horizon
co_cisrp<-CI_SRP(data=subset(edata, Catchment=="Certovo" & Horizon=="Litter"))
co_cisrp$Parameters
co_cisrp$Goodness$Gfit

###Organic topsoil
ca_cisrp<-CI_SRP(data=subset(edata, Catchment=="Certovo" & Horizon=="Organic topsoil"))
ca_cisrp$Parameters
ca_cisrp$Goodness$Gfit

#Add results to a epredts data frame
epredts$v_srp<-NA
for(i in 1:nrow(epredts)){
  if(epredts$Catchment[i]=="Plesne" & epredts$Horizon[i]=="Litter"){
    epredts$v_srp[i]<-plo_cisrp$Parameters[1]*epredts$Substrate[i]/(plo_cisrp$Parameters[2]*(1+epredts$InhibitorSRP[i]/plo_cisrp$Parameters[3]) + epredts$Substrate[i])
  }else{
    if(epredts$Catchment[i]=="Plesne" & epredts$Horizon[i]=="Organic topsoil"){
      epredts$v_srp[i]<-pla_cisrp$Parameters[1]*epredts$Substrate[i]/(pla_cisrp$Parameters[2]*(1+epredts$InhibitorSRP[i]/pla_cisrp$Parameters[3]) + epredts$Substrate[i])
    }else{
      if(epredts$Catchment[i]=="Certovo" & epredts$Horizon[i]=="Litter"){
        epredts$v_srp[i]<-co_cisrp$Parameters[1]*epredts$Substrate[i]/(co_cisrp$Parameters[2]*(1+epredts$InhibitorSRP[i]/co_cisrp$Parameters[3]) + epredts$Substrate[i])
      }else{
        epredts$v_srp[i]<-ca_cisrp$Parameters[1]*epredts$Substrate[i]/(ca_cisrp$Parameters[2]*(1+epredts$InhibitorSRP[i]/ca_cisrp$Parameters[3]) + epredts$Substrate[i])
      }
    }
  }
}

###relative rates
epredts$v_srprel<-NA
for(l in unique(epredts$Horizon)){
  for(k in unique(epredts$Catchment)){
    for(i in unique(epredts[(epredts$Horizon==l & epredts$Catchment==k), "Substrate"])){
      for(n in unique(epredts[(epredts$Horizon==l & epredts$Catchment==k), "InhibitorSRP"])){
        epredts[(epredts$Horizon==l & epredts$Catchment==k & epredts$Substrate==i & epredts$InhibitorSRP==n), "v_srprel"]<-
          epredts[(epredts$Horizon==l & epredts$Catchment==k & epredts$Substrate==i & epredts$InhibitorSRP==n), "v_srp"]/
          epredts[(epredts$Horizon==l & epredts$Catchment==k & epredts$Substrate==i & epredts$InhibitorSRP==0), "v_srp"]*100
      }
    }
  }
}

ggplot(erates, aes(InhibitorSRP, v5rel)) + geom_point(cex=6, pch=21, aes(fill = as.factor(Substrate3))) +
  facet_grid(.~Legend) + theme_min + scale_fill_manual(values=c("black", "grey30", "grey60", "grey90", "white")) +
  scale_y_continuous(limits = c(0, 110), breaks = c(0, 20, 40, 60, 80, 100)) +
  xlim(0, 16) + geom_errorbar(aes(ymin=v5rel-v5rel.se, ymax=v5rel+v5rel.se), width=0.01) +
  geom_line(data=epredts, aes(InhibitorSRP, v_srprel, color=as.factor(Substrate3), linetype=as.factor(Substrate3)), lwd=1.2) +
  scale_color_manual(values=c("black", "grey30", "grey60", "grey90", "grey90")) +
  scale_linetype_manual(values=c(rep("solid", 4), "longdash"))

#############################################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Acknkowledging different initial Porg concentration~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#############################################################################################################################################
#Vizualize
ggplot(subset(pdata, Time!=0), aes(SRP_a, log(Porg))) + geom_point(cex=6, aes(colour = as.factor(Horizon), shape=Catchment)) +
  theme_min + stat_smooth(method=lm, se=F, aes(color=Horizon))

#Generate the relationship for each horizon separately
pdata$Inhibitor<-pdata$SRP_a
porg_o<-lm(log(Porg)~Inhibitor, pdata, subset = Horizon=="Litter")
summary(porg_o)

porg_a<-lm(log(Porg)~Inhibitor, pdata, subset = Horizon=="Organic topsoil")
summary(porg_a)

#Add to data frame
edata$Porg<-NA
edata[edata$Horizon=="Litter", "Porg"]<-exp(predict(porg_o, newdata=edata[edata$Horizon=="Litter", ]))
edata[edata$Horizon!="Litter", "Porg"]<-exp(predict(porg_a, newdata=edata[edata$Horizon!="Litter", ]))

#Run the competitive inhibition model again with competition between fluorescent substrate and Porg
##Load the function
source("CI_Porg2.R") #2. scenario
##Plesne catchment
###Litter horizon
plo_cip2<-CI_Porg2(data=subset(edata, Catchment=="Plesne" & Horizon=="Litter"))
plo_cip2$Parameters
plo_cip2$Goodness$Gfit

###Organic topsoil
pla_cip2<-CI_Porg2(data=subset(edata, Catchment=="Plesne" & Horizon=="Organic topsoil"))
pla_cip2$Parameters
pla_cip2$Goodness$Gfit

##Certovo catchment
###Litter horizon
co_cip2<-CI_Porg2(data=subset(edata, Catchment=="Certovo" & Horizon=="Litter"))
co_cip2$Parameters
co_cip2$Goodness$Gfit

###Organic topsoil
ca_cip2<-CI_Porg2(data=subset(edata, Catchment=="Certovo" & Horizon=="Organic topsoil"))
ca_cip2$Parameters
ca_cip2$Goodness$Gfit

#Add results to a epredts data frame
epredts$v_p2<-NA
epredts$Porg<-NA
epredts$Inhibitor<-epredts$InhibitorSRP
epredts[epredts$Horizon=="Litter", "Porg"]<-exp(predict(porg_o, newdata=epredts[epredts$Horizon=="Litter", ]))
epredts[epredts$Horizon!="Litter", "Porg"]<-exp(predict(porg_a, newdata=epredts[epredts$Horizon!="Litter", ]))

for(i in 1:nrow(epredts)){
  if(epredts$Catchment[i]=="Plesne" & epredts$Horizon[i]=="Litter"){
    epredts$v_p2[i]<-plo_cip2$Parameters[1]*epredts$Substrate[i]/
      (plo_cip2$Parameters[2]*(1+epredts$InhibitorSRP[i]/plo_cip2$Parameters[3])*(1+epredts$Porg[i]/plo_cip2$Parameters[4]) + epredts$Substrate[i])
  }else{
    if(epredts$Catchment[i]=="Plesne" & epredts$Horizon[i]=="Organic topsoil"){
      epredts$v_p2[i]<-pla_cip2$Parameters[1]*epredts$Substrate[i]/
        (pla_cip2$Parameters[2]*(1+epredts$InhibitorSRP[i]/pla_cip2$Parameters[3])*(1+epredts$Porg[i]/pla_cip2$Parameters[4]) + epredts$Substrate[i])
    }else{
      if(epredts$Catchment[i]=="Certovo" & epredts$Horizon[i]=="Litter"){
        epredts$v_p2[i]<-co_cip2$Parameters[1]*epredts$Substrate[i]/
          (co_cip2$Parameters[2]*(1+epredts$InhibitorSRP[i]/co_cip2$Parameters[3])*(1+epredts$Porg[i]/co_cip2$Parameters[4]) + epredts$Substrate[i])
      }else{
        epredts$v_p2[i]<-ca_cip2$Parameters[1]*epredts$Substrate[i]/
          (ca_cip2$Parameters[2]*(1+epredts$InhibitorSRP[i]/ca_cip2$Parameters[3])*(1+epredts$Porg[i]/ca_cip2$Parameters[4]) + epredts$Substrate[i])
      }
    }
  }
}

###relative rates
epredts$v_prel2<-NA
for(l in unique(epredts$Horizon)){
  for(k in unique(epredts$Catchment)){
    for(i in unique(epredts[(epredts$Horizon==l & epredts$Catchment==k), "Substrate"])){
      for(n in unique(epredts[(epredts$Horizon==l & epredts$Catchment==k), "InhibitorSRP"])){
        epredts[(epredts$Horizon==l & epredts$Catchment==k & epredts$Substrate==i & epredts$InhibitorSRP==n), "v_prel2"]<-
          epredts[(epredts$Horizon==l & epredts$Catchment==k & epredts$Substrate==i & epredts$InhibitorSRP==n), "v_p2"]/
          epredts[(epredts$Horizon==l & epredts$Catchment==k & epredts$Substrate==i & epredts$InhibitorSRP==0), "v_p2"]*100
      }
    }
  }
}

ggplot(erates, aes(InhibitorSRP, v5rel)) + geom_point(cex=6, pch=21, aes(fill = as.factor(Substrate3))) +
  facet_grid(.~Legend) + theme_min + scale_fill_manual(values=c("black", "grey30", "grey60", "grey90", "white")) +
  scale_y_continuous(limits = c(0, 110), breaks = c(0, 20, 40, 60, 80, 100)) +
  xlim(0, 16) + geom_errorbar(aes(ymin=v5rel-v5rel.se, ymax=v5rel+v5rel.se), width=0.01) +
  geom_line(data=epredts, aes(InhibitorSRP, v_prel2, color=as.factor(Substrate3), linetype=as.factor(Substrate3)), lwd=1.2) +
  scale_color_manual(values=c("black", "grey30", "grey60", "grey90", "grey90")) +
  scale_linetype_manual(values=c(rep("solid", 4), "longdash"))

#############################################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Testing substrate preference~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#############################################################################################################################################
##Load the function
source("Substrate_preference.R")
##Plesne catchment
###Litter horizon
plo_sp<-Substrate_preference(data=subset(edata, Catchment=="Plesne" & Horizon=="Litter"))
plo_sp$Parameters
plo_sp$Goodness$Gfit

###Organic topsoil
pla_sp<-Substrate_preference(data=subset(edata, Catchment=="Plesne" & Horizon=="Organic topsoil"))
pla_sp$Parameters
pla_sp$Goodness$Gfit

##Certovo catchment
###Litter horizon
co_sp<-Substrate_preference(data=subset(edata, Catchment=="Certovo" & Horizon=="Litter"))
co_sp$Parameters
co_sp$Goodness$Gfit

###Organic topsoil
ca_sp<-Substrate_preference(data=subset(edata, Catchment=="Certovo" & Horizon=="Organic topsoil"))
ca_sp$Parameters
ca_sp$Goodness$Gfit


#Add results to a epredts data frame
epredts$v_sp<-NA

for(i in 1:nrow(epredts)){
  if(epredts$Catchment[i]=="Plesne" & epredts$Horizon[i]=="Litter"){
    epredts$v_sp[i]<-plo_sp$Parameters[1]*epredts$Substrate[i]/
      (plo_sp$Parameters[2]*(1+epredts$Porg[i]/plo_sp$Parameters[3]) + epredts$Substrate[i])
  }else{
    if(epredts$Catchment[i]=="Plesne" & epredts$Horizon[i]=="Organic topsoil"){
      epredts$v_sp[i]<-pla_sp$Parameters[1]*epredts$Substrate[i]/
        (pla_sp$Parameters[2]*(1+epredts$Porg[i]/pla_sp$Parameters[3]) + epredts$Substrate[i])
    }else{
      if(epredts$Catchment[i]=="Certovo" & epredts$Horizon[i]=="Litter"){
        epredts$v_sp[i]<-co_sp$Parameters[1]*epredts$Substrate[i]/
          (co_sp$Parameters[2]*(1+epredts$Porg[i]/co_sp$Parameters[3]) + epredts$Substrate[i])
      }else{
        epredts$v_sp[i]<-ca_sp$Parameters[1]*epredts$Substrate[i]/
          (ca_sp$Parameters[2]*(1+epredts$Porg[i]/ca_sp$Parameters[3]) + epredts$Substrate[i])
      }
    }
  }
}

###relative rates
epredts$v_sprel<-NA
for(l in unique(epredts$Horizon)){
  for(k in unique(epredts$Catchment)){
    for(i in unique(epredts[(epredts$Horizon==l & epredts$Catchment==k), "Substrate"])){
      for(n in unique(epredts[(epredts$Horizon==l & epredts$Catchment==k), "InhibitorSRP"])){
        epredts[(epredts$Horizon==l & epredts$Catchment==k & epredts$Substrate==i & epredts$InhibitorSRP==n), "v_sprel"]<-
          epredts[(epredts$Horizon==l & epredts$Catchment==k & epredts$Substrate==i & epredts$InhibitorSRP==n), "v_sp"]/
          epredts[(epredts$Horizon==l & epredts$Catchment==k & epredts$Substrate==i & epredts$InhibitorSRP==0), "v_sp"]*100
      }
    }
  }
}

ggplot(erates, aes(InhibitorSRP, v5rel)) + geom_point(cex=6, pch=21, aes(fill = as.factor(Substrate3))) +
  facet_grid(.~Legend) + theme_min + scale_fill_manual(values=c("black", "grey30", "grey60", "grey90", "white")) +
  scale_y_continuous(limits = c(0, 110), breaks = c(0, 20, 40, 60, 80, 100)) +
  xlim(0, 16) + geom_errorbar(aes(ymin=v5rel-v5rel.se, ymax=v5rel+v5rel.se), width=0.01) +
  geom_line(data=epredts, aes(InhibitorSRP, v_sprel, color=as.factor(Substrate3), linetype=as.factor(Substrate3)), lwd=1.2) +
  scale_color_manual(values=c("black", "grey30", "grey60", "grey90", "grey90")) +
  scale_linetype_manual(values=c(rep("solid", 4), "longdash"))

#############################################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Testing the two pools of enzyme~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#############################################################################################################################################
#Load the function
source("two_epools.R")
##Plesne catchment
###Litter horizon
plo_2ci<-two_epools(data=subset(edata, Catchment=="Plesne" & Horizon=="Litter"))
plo_2ci$Parameters
plo_2ci$Goodness$Gfit

###Organic topsoil
pla_2ci<-two_epools(data=subset(edata, Catchment=="Plesne" & Horizon=="Organic topsoil"))
pla_2ci$Parameters
pla_2ci$Goodness$Gfit

##Certovo catchment
###Litter horizon
co_2ci<-two_epools(data=subset(edata, Catchment=="Certovo" & Horizon=="Litter"))
co_2ci$Parameters
co_2ci$Goodness$Gfit

###Organic topsoil
ca_2ci<-two_epools(data=subset(edata, Catchment=="Certovo" & Horizon=="Organic topsoil"))
ca_2ci$Parameters
ca_2ci$Goodness$Gfit

###Add to a previous data frame
epredts$v_two<-NA
for(i in 1:nrow(epredts)){
  if(epredts$Catchment[i]=="Plesne" & epredts$Horizon[i]=="Litter"){
    epredts$v_two[i]<-plo_2ci$Parameters[1]*epredts$Substrate[i]/(plo_2ci$Parameters[2] + epredts$Substrate[i])+
      plo_2ci$Parameters[3]*epredts$Substrate[i]/(plo_2ci$Parameters[4]*(1+epredts$InhibitorSRP[i]/plo_2ci$Parameters[5]) + epredts$Substrate[i])
  }else{
    if(epredts$Catchment[i]=="Plesne" & epredts$Horizon[i]=="Organic topsoil"){
      epredts$v_two[i]<-pla_2ci$Parameters[1]*epredts$Substrate[i]/(pla_2ci$Parameters[2] + epredts$Substrate[i])+
        pla_2ci$Parameters[3]*epredts$Substrate[i]/(pla_2ci$Parameters[4]*(1+epredts$InhibitorSRP[i]/pla_2ci$Parameters[5]) + epredts$Substrate[i])
    }else{
      if(epredts$Catchment[i]=="Certovo" & epredts$Horizon[i]=="Litter"){
        epredts$v_two[i]<-co_2ci$Parameters[1]*epredts$Substrate[i]/(co_2ci$Parameters[2] + epredts$Substrate[i])+
          co_2ci$Parameters[3]*epredts$Substrate[i]/(co_2ci$Parameters[4]*(1+epredts$InhibitorSRP[i]/co_2ci$Parameters[5]) + epredts$Substrate[i])
      }else{
        epredts$v_two[i]<-ca_2ci$Parameters[1]*epredts$Substrate[i]/(ca_2ci$Parameters[2] + epredts$Substrate[i])+
          ca_2ci$Parameters[3]*epredts$Substrate[i]/(ca_2ci$Parameters[4]*(1+epredts$InhibitorSRP[i]/ca_2ci$Parameters[5]) + epredts$Substrate[i])
      }
    }
  }
}

###relative rates
epredts$v_tworel<-NA
for(l in unique(epredts$Horizon)){
  for(k in unique(epredts$Catchment)){
    for(i in unique(epredts[(epredts$Horizon==l & epredts$Catchment==k), "Substrate"])){
      for(n in unique(epredts[(epredts$Horizon==l & epredts$Catchment==k), "InhibitorSRP"])){
        epredts[(epredts$Horizon==l & epredts$Catchment==k & epredts$Substrate==i & epredts$InhibitorSRP==n), "v_tworel"]<-
          epredts[(epredts$Horizon==l & epredts$Catchment==k & epredts$Substrate==i & epredts$InhibitorSRP==n), "v_two"]/
          epredts[(epredts$Horizon==l & epredts$Catchment==k & epredts$Substrate==i & epredts$InhibitorSRP==0), "v_two"]*100
      }
    }
  }
}

ggplot(erates, aes(InhibitorSRP, v5rel)) + geom_point(cex=6, pch=21, aes(fill = as.factor(Substrate3))) +
  facet_grid(.~Legend) + theme_min + scale_fill_manual(values=c("black", "grey30", "grey60", "grey90", "white")) +
  scale_y_continuous(limits = c(0, 110), breaks = c(0, 20, 40, 60, 80, 100)) +
  xlim(0, 16) + geom_errorbar(aes(ymin=v5rel-v5rel.se, ymax=v5rel+v5rel.se), width=0.01) +
  geom_line(data=epredts, aes(InhibitorSRP, v_tworel, color=as.factor(Substrate3), linetype=as.factor(Substrate3)), lwd=1.2) +
  scale_color_manual(values=c("black", "grey30", "grey60", "grey90", "grey90")) +
  scale_linetype_manual(values=c(rep("solid", 4), "longdash"))

