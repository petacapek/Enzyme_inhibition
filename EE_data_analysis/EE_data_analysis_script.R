###############################################################################################################################################
###############################################################################################################################################
################################################Biochemical inhibition of soil phosphatase activity############################################
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
library(rgenoud)
library(ABCoptim)
library(gridExtra)
library(reshape)
library(segmented)
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
    pdata$Legend[i]<-"Plešné - Litter"
  }else{
    if(pdata$Catchment[i]=="Plesne" & pdata$Horizon[i]=="Organic topsoil"){
      pdata$Legend[i]<-"Plešné - Organic topsoil"
    }else{
      if(pdata$Catchment[i]=="Certovo" & pdata$Horizon[i]=="Litter"){
        pdata$Legend[i]<-"Čertovo - Litter"
      }else{
        pdata$Legend[i]<-"Čertovo - Organic topsoil"
      }
    }
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Initial product (P-PO4 measured as SRP) concentration~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
pdata$Legend<-factor(pdata$Legend, levels=c("Plešné - Litter", "Plešné - Organic topsoil", "Čertovo - Litter", "Čertovo - Organic topsoil"))
#Measured vs added SRP
ggplot(pdata[pdata$Time==0, ], aes(SRP_a, SRP_o)) + geom_point(cex=6, aes(fill=Legend, shape=Legend)) +
  theme_min + scale_fill_manual(values = c("black", "black", "grey", "grey")) + 
  scale_shape_manual(values = c(21, 22, 21, 22)) + 
  theme(legend.title = element_blank(), legend.position = c(0.2, 0.8)) +
  ylim(0, 20) + xlim(0, 20) + geom_abline(intercept = 0, slope=1, color="black") +
  stat_smooth(method="lm", lty=2, lwd=1.5, se=F, colour="black") +
  ylab(expression(paste("Measured SRP (", mu, "mol ", g^{-1}, ")"))) + 
  xlab(expression(paste("Added P-P", O[4]," (", mu, "mol ", g^{-1}, ")")))

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
plo_wi<-nls(wi, data = edata, subset = Catchment=="Plesne" & Horizon=="Litter" & time<120, start=list(Vmax=1, Km=1))
plo_ci<-nls(ci, data = edata, subset = Catchment=="Plesne" & Horizon=="Litter" & time<120, start=list(Vmax=1, Km=1, Kic=1))
plo_nci<-nls(nci, data = edata, subset = Catchment=="Plesne" & Horizon=="Litter" & time<120, start=list(Vmax=0.1, Km=10, Kiu=10))
plo_uci<-nls(uci, data = edata, subset = Catchment=="Plesne" & Horizon=="Litter" & time<120, start=list(Vmax=0.1, Km=10, Kiu=20))
###AIC test
AICtab(plo_wi, plo_ci, plo_nci, plo_uci, weights=T, sort=T, base=T, logLik=T)
###Likelihood ratio test
-2*(as.numeric(logLik(plo_wi) - as.numeric(logLik(plo_ci))))
pchisq(-2*(as.numeric(logLik(plo_wi) - as.numeric(logLik(plo_ci)))), df=1, lower.tail=F)
-2*(as.numeric(logLik(plo_wi) - as.numeric(logLik(plo_nci))))
pchisq(-2*(as.numeric(logLik(plo_wi) - as.numeric(logLik(plo_nci)))), df=1, lower.tail=F)
-2*(as.numeric(logLik(plo_wi) - as.numeric(logLik(plo_uci))))
pchisq(-2*(as.numeric(logLik(plo_wi) - as.numeric(logLik(plo_uci)))), df=1, lower.tail=F)

###Parameters
coef(plo_wi)
coef(plo_ci)
coef(plo_nci)
coef(plo_uci)

##Organic topsoil horizon
pla_wi<-nls(wi, data = edata, subset = Catchment=="Plesne" & Horizon=="Organic topsoil" & time<120, start=list(Vmax=1, Km=1))
pla_ci<-nls(ci, data = edata, subset = Catchment=="Plesne" & Horizon=="Organic topsoil" & time<120, start=list(Vmax=1, Km=1, Kic=1))
pla_nci<-nls(nci, data = edata, subset = Catchment=="Plesne" & Horizon=="Organic topsoil" & time<120, start=list(Vmax=0.1, Km=10, Kiu=10))
pla_uci<-nls(nci, data = edata, subset = Catchment=="Plesne" & Horizon=="Organic topsoil" & time<120, start=list(Vmax=0.1, Km=10, Kiu=10))
###AIC test
AICtab(pla_wi, pla_ci, pla_nci, pla_uci, weights=T, sort=T, base=T, logLik=T)
###Likelihood ratio test
-2*(as.numeric(logLik(pla_wi) - as.numeric(logLik(pla_ci))))
pchisq(-2*(as.numeric(logLik(pla_wi) - as.numeric(logLik(pla_ci)))), df=1, lower.tail=F)
-2*(as.numeric(logLik(pla_wi) - as.numeric(logLik(pla_nci))))
pchisq(-2*(as.numeric(logLik(pla_wi) - as.numeric(logLik(pla_nci)))), df=1, lower.tail=F)
-2*(as.numeric(logLik(pla_wi) - as.numeric(logLik(pla_uci))))
pchisq(-2*(as.numeric(logLik(pla_wi) - as.numeric(logLik(pla_uci)))), df=1, lower.tail=F)

###Parameters
coef(pla_wi)
coef(pla_ci)
coef(pla_nci)
coef(pla_uci)

#Certovo catchment
##Litter horizon
co_wi<-nls(wi, data = edata, subset = Catchment=="Certovo" & Horizon=="Litter" & time<120, start=list(Vmax=1, Km=1))
co_ci<-nls(ci, data = edata, subset = Catchment=="Certovo" & Horizon=="Litter" & time<120, start=list(Vmax=1, Km=1, Kic=1))
co_nci<-nls(nci, data = edata, subset = Catchment=="Certovo" & Horizon=="Litter" & time<120, start=list(Vmax=0.1, Km=10, Kiu=10))
co_uci<-nlsLM(uci, data = edata, subset = Catchment=="Certovo" & Horizon=="Litter" & time<120, start=list(Vmax=1.859265, Km=210.954127, Kiu=30.958534))
###AIC test
AICtab(co_wi, co_ci, co_nci, co_uci, weights=T, sort=T, base=T, logLik=T)
###Likelihood ratio test
-2*(as.numeric(logLik(co_wi) - as.numeric(logLik(co_ci))))
pchisq(-2*(as.numeric(logLik(co_wi) - as.numeric(logLik(co_ci)))), df=1, lower.tail=F)
-2*(as.numeric(logLik(co_wi) - as.numeric(logLik(co_nci))))
pchisq(-2*(as.numeric(logLik(co_wi) - as.numeric(logLik(co_nci)))), df=1, lower.tail=F)
-2*(as.numeric(logLik(co_wi) - as.numeric(logLik(co_uci))))
pchisq(-2*(as.numeric(logLik(co_wi) - as.numeric(logLik(co_uci)))), df=1, lower.tail=F)

###Parameters
coef(co_wi)
coef(co_ci)
coef(co_nci)
coef(co_uci)

##Organic topsoil horizon
ca_wi<-nls(wi, data = edata, subset = Catchment=="Certovo" & Horizon=="Organic topsoil" & time<120, start=list(Vmax=1, Km=1))
ca_ci<-nls(ci, data = edata, subset = Catchment=="Certovo" & Horizon=="Organic topsoil" & time<120, start=list(Vmax=1, Km=1, Kic=1))
ca_nci<-nls(nci, data = edata, subset = Catchment=="Certovo" & Horizon=="Organic topsoil" & time<120, start=list(Vmax=0.1, Km=10, Kiu=10))
ca_uci<-nlsLM(uci, data = edata, subset = Catchment=="Certovo" & Horizon=="Organic topsoil" & time<120, start=list(Vmax=0.1, Km=10, Kiu=10))
###AIC test
AICtab(ca_wi, ca_ci, ca_nci, ca_uci, weights=T, sort=T, base=T, logLik=T)
###Likelihood ratio test
-2*(as.numeric(logLik(ca_wi) - as.numeric(logLik(ca_ci))))
pchisq(-2*(as.numeric(logLik(ca_wi) - as.numeric(logLik(ca_ci)))), df=1, lower.tail=F)
-2*(as.numeric(logLik(ca_wi) - as.numeric(logLik(ca_nci))))
pchisq(-2*(as.numeric(logLik(ca_wi) - as.numeric(logLik(ca_nci)))), df=1, lower.tail=F)
-2*(as.numeric(logLik(ca_wi) - as.numeric(logLik(ca_uci))))
pchisq(-2*(as.numeric(logLik(ca_wi) - as.numeric(logLik(ca_uci)))), df=1, lower.tail=F)
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

##Plot
###add lumped factor for graphics
erates$Legend<-NA
for(i in 1:nrow(erates)){
  if(erates$Catchment[i]=="Plesne" & erates$Horizon[i]=="Litter"){
    erates$Legend[i]<-"Plešné - Litter"
  }else{
    if(erates$Catchment[i]=="Plesne" & erates$Horizon[i]=="Organic topsoil"){
      erates$Legend[i]<-"Plešné - Organic topsoil"
    }else{
      if(erates$Catchment[i]=="Certovo" & erates$Horizon[i]=="Litter"){
        erates$Legend[i]<-"Čertovo - Litter"
      }else{
        erates$Legend[i]<-"Čertovo - Organic topsoil"
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
    epredts$Legend[i]<-"Plešné - Litter"
  }else{
    if(epredts$Catchment[i]=="Plesne" & epredts$Horizon[i]=="Organic topsoil"){
      epredts$Legend[i]<-"Plešné - Organic topsoil"
    }else{
      if(epredts$Catchment[i]=="Certovo" & epredts$Horizon[i]=="Litter"){
        epredts$Legend[i]<-"Čertovo - Litter"
      }else{
        epredts$Legend[i]<-"Čertovo - Organic topsoil"
      }
    }
  }
}

epredts$Substrate2<-round(epredts$Substrate, 0)
epredts$Substrate3<-rep(epredts$Substrate2[1:250], times=4)

erates$Legend<-factor(erates$Legend, levels=c("Plešné - Litter", "Plešné - Organic topsoil", "Čertovo - Litter", "Čertovo - Organic topsoil"))
epredts$Legend<-factor(epredts$Legend, levels=c("Plešné - Litter", "Plešné - Organic topsoil", "Čertovo - Litter", "Čertovo - Organic topsoil"))

ggplot(erates, aes(InhibitorSRP, v5rel)) + geom_point(cex=6, pch=21, aes(fill = as.factor(Substrate3))) +
  facet_grid(.~Legend) + theme_min + scale_fill_manual(values=c("black", "grey30", "grey60", "grey90", "white")) +
  scale_y_continuous(limits = c(0, 110), breaks = c(0, 20, 40, 60, 80, 100)) +
  xlim(0, 16) + geom_errorbar(aes(ymin=v5rel-v5rel.se, ymax=v5rel+v5rel.se), width=0.01) +
  geom_line(data=epredts, aes(InhibitorSRP, vrel, color=as.factor(Substrate3), linetype=as.factor(Substrate3)), lwd=1.2, show.legend = F) +
  scale_color_manual(values=c("black", "grey30", "grey60", "grey90", "grey90")) +
  scale_linetype_manual(values=c(rep("solid", 4), "longdash")) +
  ylab("Relative potential phosphatase activity  (%)") +
  xlab(expression(paste("Initial SRP (", mu, "mol ", g^{-1}, ")"))) +
  labs(fill=expression(paste("MUB-P (", mu, "mol ", g^{-1}, ")"))) +
  theme(legend.position = c(0.35,0.3))

#Piecewise linear regression
##Run for saturating concentrations of substrate and no inhibitor
edata_sat<-subset(edata, Substrate>60 & Inhibitor==0)

##Define the breakpoints
edata_sat$bp<-NA
###Plesne
####Litter horizon
edata_sat[(edata_sat$Catchment=="Plesne" & edata_sat$Horizon=="Litter"), "bp"]<-
  summary(segmented(lm(Product~time-1, subset(edata_sat, Catchment=="Plesne" & Horizon=="Litter"))))$psi[2]

slope(segmented(lm(Product~time-1, subset(edata_sat, Catchment=="Certovo" & Horizon=="Organic topsoil"))))

####Organic topsoil
edata_sat[(edata_sat$Catchment=="Plesne" & edata_sat$Horizon=="Organic topsoil"), "bp"]<-
  summary(segmented(lm(Product~time-1, subset(edata_sat, Catchment=="Plesne" & Horizon=="Organic topsoil"))))$psi[2]

###Certovo
####Litter horizon
edata_sat[(edata_sat$Catchment=="Certovo" & edata_sat$Horizon=="Litter"), "bp"]<-
  summary(segmented(lm(Product~time-1, subset(edata_sat, Catchment=="Certovo" & Horizon=="Litter"))))$psi[2]

####Organic topsoil
edata_sat[(edata_sat$Catchment=="Certovo" & edata_sat$Horizon=="Organic topsoil"), "bp"]<-
  summary(segmented(lm(Product~time-1, subset(edata_sat, Catchment=="Certovo" & Horizon=="Organic topsoil"))))$psi[2]

##Do the linear regression for defined time intervals 
erates_sat1<-data.frame(Catchment=c(rep("Plesne", times=2), rep("Certovo", times=2)),
                   Horizon=c(rep("Litter", times=1), rep("Organic topsoil", times=1),
                             rep("Litter", times=1), rep("Organic topsoil", times=1)))
erates_sat1$v<-NA
erates_sat1$v.se<-NA
erates_sat1$Interval<-c("To the breakpoint")

for(i in unique(edata_sat$Catchment)){
  for(n in unique(edata_sat$Horizon)){
    lmv<-lm(Product~time-1, edata_sat[(edata_sat$Catchment==i & edata_sat$Horizon==n & 
                                     edata_sat$time>=0 & 
                                       edata_sat$time<unique(edata_sat[(edata_sat$Catchment==i & edata_sat$Horizon==n), "bp"])), ])
    erates_sat1[(erates_sat1$Catchment==i & erates_sat1$Horizon==n), "v"]<-summary(lmv)$coefficients[1]
    erates_sat1[(erates_sat1$Catchment==i & erates_sat1$Horizon==n), "v.se"]<-summary(lmv)$coefficients[2]
  }
}

erates_sat2<-data.frame(Catchment=c(rep("Plesne", times=2), rep("Certovo", times=2)),
                        Horizon=c(rep("Litter", times=1), rep("Organic topsoil", times=1),
                                  rep("Litter", times=1), rep("Organic topsoil", times=1)))
erates_sat2$v<-NA
erates_sat2$v.se<-NA
erates_sat2$Interval<-c("From the breakpoint")

for(i in unique(edata_sat$Catchment)){
  for(n in unique(edata_sat$Horizon)){
    lmv<-lm(Product~time-1, edata_sat[(edata_sat$Catchment==i & edata_sat$Horizon==n & 
                                         edata_sat$time>unique(edata_sat[(edata_sat$Catchment==i & edata_sat$Horizon==n), "bp"])), ])
    erates_sat2[(erates_sat2$Catchment==i & erates_sat2$Horizon==n), "v"]<-summary(lmv)$coefficients[1]
    erates_sat2[(erates_sat2$Catchment==i & erates_sat2$Horizon==n), "v.se"]<-summary(lmv)$coefficients[2]
  }
}

Erts<-rbind(erates_sat1, erates_sat2)
Erts$Catchment2<-ifelse(Erts$Catchment=="Plesne", "Plešné", "Čertovo")
Erts$Catchment2<-factor(Erts$Catchment2, levels = c("Plešné", "Čertovo"))
Erts$Interval<-factor(Erts$Interval, levels = c("To the breakpoint", "From the breakpoint"))
Erts$bp<-c("11 min", "67 min", "23 min", "62 min", NA, NA, NA, NA)
Erts$bpx<-c(0.1, 0.04, 0.1, 0.04, 0.1, 0.04, 0.1, 0.04)

ggplot(Erts, aes(Catchment2, v)) + geom_bar(aes(fill=Interval), stat = "identity", position=position_dodge(), color="black") +
  facet_wrap(~Horizon, scales="free_y") + theme_min + 
  scale_fill_manual(values = c("white", "grey")) +
  geom_errorbar(aes(ymin=v, ymax=v+v.se, color=Interval), position=position_dodge(), show.legend = F) +
  scale_color_manual(values = c("black", "black")) +
  theme(axis.title.x = element_blank(), legend.title = element_blank(), legend.position = c(0.12, 0.9)) +
  ylab(expression(paste("Potential phosphatase activity (", mu, "mol ", g^{-1}, min^{-1}, ")"))) +
  geom_text(aes(Catchment2, bpx, label=bp, hjust=1.1, vjust=8, fontface="italic"), cex=8)
  

#Calculate the explained variability in concentration of product over time
source("CI_model.R")
##Plesne catchment
###Litter horizon
plo_ci2<-CI_model(data=subset(edata, Catchment=="Plesne" & Horizon=="Litter" & time<120), parameters = coef(plo_ci))
plo_ci2$Gfit

###Organic topsoil
pla_ci2<-CI_model(data=subset(edata, Catchment=="Plesne" & Horizon=="Organic topsoil" & time<120), parameters = coef(pla_ci))
pla_ci2$Gfit

##Certovo catchment
###Litter horizon
co_ci2<-CI_model(data=subset(edata, Catchment=="Certovo" & Horizon=="Litter" & time<120), parameters = coef(co_ci))
co_ci2$Gfit

###Organic topsoil
ca_ci2<-CI_model(data=subset(edata, Catchment=="Certovo" & Horizon=="Organic topsoil" & time<120), parameters = coef(ca_ci))
ca_ci2$Gfit

#############################################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Acknkowledging different initial Porg concentration~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#############################################################################################################################################
#Vizualize
ggplot(subset(pdata, Time!=0), aes(SRP_a, log(Porg))) + geom_point(cex=6, aes(colour = as.factor(Horizon), shape=Catchment)) +
  theme_min + stat_smooth(method=lm, se=F, aes(color=Horizon))

pdata_org<-subset(pdata, !is.na(Porg))
pdata_org$SRP_a2
pdata_org[pdata_org$SRP_a2==9, "SRP_a2"]<-8
pdata_org[pdata_org$SRP_a2==17, "SRP_a2"]<-16

grid.arrange(
  ggplot(pdata[pdata$Time==0, ], aes(SRP_a, SRP_o)) + geom_point(cex=6, aes(fill=Legend, shape=Legend)) +
               theme_min + scale_fill_manual(values = c("black", "black", "grey", "grey")) + 
               scale_shape_manual(values = c(21, 22, 21, 22)) + 
               theme(legend.title = element_blank(), legend.position = c(0.28, 0.8)) +
               ylim(0, 20) + xlim(0, 20) + geom_abline(intercept = 0, slope=1, color="black") +
               stat_smooth(method="lm", lty=2, lwd=1.5, se=F, colour="black") +
               ylab(expression(paste("Measured SRP (", mu, "mol ", g^{-1}, ")"))) + 
               xlab(expression(paste("Added P-P", O[4]," (", mu, "mol ", g^{-1}, ")"))) +
               ggtitle("A)"),
  pdata_org %>% group_by(Horizon, SRP_a2) %>% summarize(y=mean(Porg), y.sd=sd(Porg)) %>%
    ggplot(aes(as.factor(SRP_a2), y)) + geom_bar(aes(fill=Horizon), stat = "identity", position=position_dodge(), color="black") +
    theme_min + scale_fill_manual(values = c("black", "grey")) + 
    geom_errorbar(aes(ymin=y, ymax=y+y.sd, color=Horizon), position=position_dodge(), show.legend = F) +
    scale_color_manual(values = c("black", "black")) +
    theme(legend.title = element_blank(), legend.position = c(0.2, 0.85)) +
    ylab(expression(paste("Measured DOP (", mu, "mol ", g^{-1}, ")"))) + 
    xlab(expression(paste("Added P-P", O[4]," (", mu, "mol ", g^{-1}, ")"))) +
    ggtitle("B)"), nrow=1
)

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
source("CI_Porg.R")
##Plesne catchment
###Litter horizon
plo_cip<-CI_Porg(data=subset(edata, Catchment=="Plesne" & Horizon=="Litter" & time<120))
plo_cip$Parameters
plo_cip$Goodness$Gfit

###Organic topsoil
pla_cip<-CI_Porg(data=subset(edata, Catchment=="Plesne" & Horizon=="Organic topsoil" & time<120))
pla_cip$Parameters
pla_cip$Goodness$Gfit

##Certovo catchment
###Litter horizon
co_cip<-CI_Porg(data=subset(edata, Catchment=="Certovo" & Horizon=="Litter" & time<120))
co_cip$Parameters
co_cip$Goodness$Gfit

###Organic topsoil
ca_cip<-CI_Porg(data=subset(edata, Catchment=="Certovo" & Horizon=="Organic topsoil" & time<120))
ca_cip$Parameters
ca_cip$Goodness$Gfit

#Add results to a epredts data frame
epredts$v_p<-NA
epredts$Porg<-NA
epredts$Inhibitor<-epredts$InhibitorSRP
epredts[epredts$Horizon=="Litter", "Porg"]<-exp(predict(porg_o, newdata=epredts[epredts$Horizon=="Litter", ]))
epredts[epredts$Horizon!="Litter", "Porg"]<-exp(predict(porg_a, newdata=epredts[epredts$Horizon!="Litter", ]))

for(i in 1:nrow(epredts)){
  if(epredts$Catchment[i]=="Plesne" & epredts$Horizon[i]=="Litter"){
    epredts$v_p[i]<-plo_cip$Parameters[1]*epredts$Substrate[i]/
      (plo_cip$Parameters[2]*(1+epredts$InhibitorSRP[i]/plo_cip$Parameters[3])*(1+epredts$Porg[i]/plo_cip$Parameters[4]) + epredts$Substrate[i])
  }else{
    if(epredts$Catchment[i]=="Plesne" & epredts$Horizon[i]=="Organic topsoil"){
      epredts$v_p[i]<-pla_cip$Parameters[1]*epredts$Substrate[i]/
        (pla_cip$Parameters[2]*(1+epredts$InhibitorSRP[i]/pla_cip$Parameters[3])*(1+epredts$Porg[i]/pla_cip$Parameters[4]) + epredts$Substrate[i])
    }else{
      if(epredts$Catchment[i]=="Certovo" & epredts$Horizon[i]=="Litter"){
        epredts$v_p[i]<-co_cip$Parameters[1]*epredts$Substrate[i]/
          (co_cip$Parameters[2]*(1+epredts$InhibitorSRP[i]/co_cip$Parameters[3])*(1+epredts$Porg[i]/co_cip$Parameters[4]) + epredts$Substrate[i])
      }else{
        epredts$v_p[i]<-ca_cip$Parameters[1]*epredts$Substrate[i]/
          (ca_cip$Parameters[2]*(1+epredts$InhibitorSRP[i]/ca_cip$Parameters[3])*(1+epredts$Porg[i]/ca_cip$Parameters[4]) + epredts$Substrate[i])
      }
    }
  }
}

###relative rates
epredts$v_prel<-NA
for(l in unique(epredts$Horizon)){
  for(k in unique(epredts$Catchment)){
    for(i in unique(epredts[(epredts$Horizon==l & epredts$Catchment==k), "Substrate"])){
      for(n in unique(epredts[(epredts$Horizon==l & epredts$Catchment==k), "InhibitorSRP"])){
        epredts[(epredts$Horizon==l & epredts$Catchment==k & epredts$Substrate==i & epredts$InhibitorSRP==n), "v_prel"]<-
          epredts[(epredts$Horizon==l & epredts$Catchment==k & epredts$Substrate==i & epredts$InhibitorSRP==n), "v_p"]/
          epredts[(epredts$Horizon==l & epredts$Catchment==k & epredts$Substrate==i & epredts$InhibitorSRP==0), "v_p"]*100
      }
    }
  }
}

#############################################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Testing sole DOP inhibition~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#############################################################################################################################################
##Load the function
source("Substrate_preference.R")
##Plesne catchment
###Litter horizon
plo_sp<-Substrate_preference(data=subset(edata, Catchment=="Plesne" & Horizon=="Litter" & time<120))
plo_sp$Parameters
plo_sp$Goodness$Gfit

###Organic topsoil
pla_sp<-Substrate_preference(data=subset(edata, Catchment=="Plesne" & Horizon=="Organic topsoil" & time<120))
pla_sp$Parameters
pla_sp$Goodness$Gfit

##Certovo catchment
###Litter horizon
co_sp<-Substrate_preference(data=subset(edata, Catchment=="Certovo" & Horizon=="Litter" & time<120))
co_sp$Parameters
co_sp$Goodness$Gfit

###Organic topsoil
ca_sp<-Substrate_preference(data=subset(edata, Catchment=="Certovo" & Horizon=="Organic topsoil" & time<120))
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

#############################################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Models comparison~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#############################################################################################################################################
##Plesne catchment
###Litter horizon
plo_ci2$Gfit
plo_sp$Goodness$Gfit
plo_cip$Goodness$Gfit

####logLike test
#####CI vs SP
-2*(plo_ci2$Gfit[["ll"]]-plo_sp$Goodness$Gfit[["ll"]])
pchisq(-2*(plo_ci2$Gfit[["ll"]]-plo_sp$Goodness$Gfit[["ll"]]), df=1, lower.tail=F)
#####CI vs CIP
-2*(plo_ci2$Gfit[["ll"]]-plo_cip$Goodness$Gfit[["ll"]])
pchisq(-2*(plo_ci2$Gfit[["ll"]]-plo_cip$Goodness$Gfit[["ll"]]), df=1, lower.tail=F)

###Organic topsoil
pla_ci2$Gfit
pla_sp$Goodness$Gfit
pla_cip$Goodness$Gfit

####logLike test
#####CI vs SP
-2*(pla_ci2$Gfit[["ll"]]-pla_sp$Goodness$Gfit[["ll"]])
pchisq(-2*(pla_ci2$Gfit[["ll"]]-pla_sp$Goodness$Gfit[["ll"]]), df=1, lower.tail=F)
#####CI vs CIP
-2*(pla_ci2$Gfit[["ll"]]-pla_cip$Goodness$Gfit[["ll"]])
pchisq(-2*(pla_ci2$Gfit[["ll"]]-pla_cip$Goodness$Gfit[["ll"]]), df=1, lower.tail=F)

##Certovo catchment
###Litter horizon
co_ci2$Gfit
co_sp$Goodness$Gfit
co_cip$Goodness$Gfit

####logLike test
#####CI vs SP
-2*(co_ci2$Gfit[["ll"]]-co_sp$Goodness$Gfit[["ll"]])
pchisq(-2*(co_ci2$Gfit[["ll"]]-co_sp$Goodness$Gfit[["ll"]]), df=1, lower.tail=F)
#####CI vs CIP
-2*(co_ci2$Gfit[["ll"]]-co_cip$Goodness$Gfit[["ll"]])
pchisq(-2*(co_ci2$Gfit[["ll"]]-co_cip$Goodness$Gfit[["ll"]]), df=1, lower.tail=F)

###Organic topsoil
ca_ci2$Gfit
ca_sp$Goodness$Gfit
ca_cip$Goodness$Gfit

####logLike test
#####CI vs SP
-2*(ca_ci2$Gfit[["ll"]]-ca_sp$Goodness$Gfit[["ll"]])
pchisq(-2*(ca_ci2$Gfit[["ll"]]-ca_sp$Goodness$Gfit[["ll"]]), df=1, lower.tail=F)
#####CI vs CIP
-2*(ca_ci2$Gfit[["ll"]]-ca_cip$Goodness$Gfit[["ll"]])
pchisq(-2*(ca_ci2$Gfit[["ll"]]-ca_cip$Goodness$Gfit[["ll"]]), df=1, lower.tail=F)

#Vizualization
viz_inh<-data.frame(SRP=rep(rep(seq(0.1, 10, by=0.1), 100), 4),
                    DOP=rep(rep(seq(0.1, 10, by=0.1), each=100), 4),
                    Legend=c(rep("Plešné - Litter", 1e4),
                             rep("Plešné - Organic topsoil", 1e4),
                             rep("Čertovo - Litter", 1e4),
                             rep("Čertovo - Organic topsoil", 1e4)),
                    Catchment=c(rep("Plešné", 2e4), rep("Čertovo", 2e4)),
                    Horizon=c(rep("Litter", 1e4), rep("Organic topsoil", 1e4),
                              rep("Litter", 1e4), rep("Organic topsoil", 1e4)))
viz_inh$activity<-NA

for(i in 1:nrow(viz_inh)){
  if(viz_inh$Legend[i]=="Plešné - Litter"){
    viz_inh$activity[i]<-(100*plo_cip$Parameters[1]/(plo_cip$Parameters[2]*
                                     (1+viz_inh$SRP[i]/plo_cip$Parameters[3])*
                                     (1+viz_inh$DOP[i]/plo_cip$Parameters[4]) + 100))/(100*plo_cip$Parameters[1]/(plo_cip$Parameters[2]+ 100))*100
  }else{
    if(viz_inh$Legend[i]=="Plešné - Organic topsoil"){
      viz_inh$activity[i]<-(100*pla_cip$Parameters[1]/(pla_cip$Parameters[2]*
                                       (1+viz_inh$SRP[i]/pla_cip$Parameters[3])*
                                       (1+viz_inh$DOP[i]/pla_cip$Parameters[4]) + 100))/(100*pla_cip$Parameters[1]/(pla_cip$Parameters[2]+ 100))*100
    }else{
      if(viz_inh$Legend[i]=="Čertovo - Litter"){
        viz_inh$activity[i]<-(100*co_cip$Parameters[1]/(co_cip$Parameters[2]*
                                         (1+viz_inh$SRP[i]/co_cip$Parameters[3])*
                                         (1+viz_inh$DOP[i]/co_cip$Parameters[4]) + 100))/(100*co_cip$Parameters[1]/(co_cip$Parameters[2]+ 100))*100
      }else{
        viz_inh$activity[i]<-(100*ca_cip$Parameters[1]/(ca_cip$Parameters[2]*
                                         (1+viz_inh$SRP[i]/ca_cip$Parameters[3])*
                                         (1+viz_inh$DOP[i]/ca_cip$Parameters[4]) + 100))/(100*ca_cip$Parameters[1]/(ca_cip$Parameters[2]+ 100))*100
      }
    }
  }
}

viz_inh$Catchment<-factor(viz_inh$Catchment, levels = c("Plešné", "Čertovo"))

grid.arrange(ggplot(subset(viz_inh, DOP==0.1 | DOP==10), aes(SRP, activity)) + 
               geom_line(lwd=1.2, aes(linetype=Horizon, colour=as.factor(DOP))) +
               theme_min + ylim(50, 100) + facet_grid(.~Catchment) + 
               scale_linetype_manual(values = c("solid", "dotdash")) +
               scale_color_manual(values = c("black", "grey60")) +
               labs(color=expression(paste("DOP (", mu, "mol ", g^{-1}, ")")), linetype=c(" ")) +
               ggtitle("A)") + ylab(expression(atop("Potential phosphatase activity", paste("(% of the actual activity)")))) +
               xlab(expression(paste("SRP (", mu, "mol ", g^{-1}, ")"))) +
               theme(legend.direction = c("horizontal"), legend.box = c("vertical"), legend.position = c(0.25, 0.3),
                     panel.spacing = unit(2, "lines")),
             ggplot(subset(viz_inh, SRP==0.1 | SRP==10), aes(DOP, activity)) + 
               geom_line(lwd=1.2, aes(linetype=Horizon, colour=as.factor(SRP))) +
               theme_min + ylim(50, 100) + facet_grid(.~Catchment) + 
               scale_linetype_manual(values = c("solid", "dotdash")) +
               scale_color_manual(values = c("black", "grey60")) +
               labs(color=expression(paste("SRP (", mu, "mol ", g^{-1}, ")")), linetype=c(" ")) +
               ggtitle("B)") + ylab(expression(atop("Potential phosphatase activity", paste("(% of the actual activity)")))) +
               xlab(expression(paste("DOP (", mu, "mol ", g^{-1}, ")"))) +
               theme(legend.direction = c("horizontal"), legend.box = c("vertical"), legend.position = c(0.25, 0.3),
                     panel.spacing = unit(2, "lines")),
             nrow=2)


#############################################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Parameter values predictors~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#############################################################################################################################################
#Values
pars_all<-as.data.frame(rbind(plo_cip$Parameters, pla_cip$Parameters,
                              co_cip$Parameters, ca_cip$Parameters))
pars_all$Catchment<-c("Plešné", "Plešné", "Čertovo", "Čertovo")
pars_all$Horizon<-c("Litter", "Organic topsoil", "Litter", "Organic topsoil")
pars_all$Legend<-c("Plešné - Litter", "Plešné - Organic topsoil", 
                   "Čertovo - Litter", "Čertovo - Organic topsoil")
pars_all$Legend<-factor(pars_all$Legend, levels = c("Plešné - Litter", "Plešné - Organic topsoil", 
                                                    "Čertovo - Litter", "Čertovo - Organic topsoil"))
#Predictors
pars_all$pH<-c(4.51, 3.58, 4.06, 3.61)
pars_all$MBP<-c(12.1, 9.9, 14.4, 5.8)
pars_all$Pnahco3<-c(2.7, 1, 2.2, 0.6)
pars_all$MBC<-c(863.4, 765.8, 851, 627.4)
pars_all$MBN<-c(42.2, 28.3, 40.8, 28.7)
pars_all$DOC<-c(71.1, 64.9, 90.4, 41.7)
pars_all$TN<-c(111.8, 49.6, 127.7, 39.3)
pars_all$Porg<-c(3.97, 2.42, 3.15,  2.23)
pars_all$SRP<-c(2.83, 2.04, 2.19, 1.41)
pars_all$DON<-c(24.2, 13.3, 26.8, 10.4)
pars_all$NH4<-c(83.9, 33.2, 93.4, 23.6)
pars_all$NO3<-c(3.69, 3.15, 7.48, 5.37)
pars_all$Ctot<-c(47, 43, 49, 38)#in %
pars_all$Ntot<-c(1.86, 1.58, 2.03, 1.53)#in %
#Vmax
ggplot(pars_all, aes(pH, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(Pnahco3, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(Porg, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(SRP, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(MBP, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(MBC, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(MBN, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(Ctot, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(Ntot, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(DOC, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(TN, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(DON, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(NH4, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(NO3, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(NH4+NO3, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(MBC/MBP, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(MBN/MBP, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(Porg/MBP, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(Porg/SRP, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(MBC/MBN, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(DOC/Porg, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(DOC/SRP, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(TN/Porg, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(TN/SRP, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(DON/Porg, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(DON/SRP, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(NH4/Porg, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(NH4/SRP, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(NO3/Porg, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(NO3/SRP, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes((NH4+NO3)/Porg, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes((NH4+NO3)/SRP, Vmax)) + geom_point(cex=6, pch=21) + theme_min
#Kmorg
ggplot(pars_all, aes(pH, Kmorg)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(Pnahco3, Kmorg)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(Porg, Kmorg)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(SRP, Kmorg)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(MBP, Kmorg)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(MBC, Kmorg)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(MBN, Kmorg)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(Ctot, Kmorg)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(Ntot, Kmorg)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(DOC, Kmorg)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(TN, Kmorg)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(DON, Kmorg)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(NH4, Kmorg)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(NO3, log(Kmorg))) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(NH4+NO3, Kmorg)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(MBC/MBP, Kmorg)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(MBN/MBP, Kmorg)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(Porg/MBP, Kmorg)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes((Porg/SRP), (Kmorg))) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(MBC/MBN, Kmorg)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(DOC/Porg, Kmorg)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(DOC/SRP, Kmorg)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(TN/Porg, Kmorg)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(TN/SRP, Kmorg)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(DON/Porg, Kmorg)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(DON/SRP, Kmorg)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(NH4/Porg, Kmorg)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(NH4/SRP, Kmorg)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(NO3/Porg, Kmorg)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(NO3/SRP, Kmorg)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes((NH4+NO3)/Porg, Kmorg)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(((NH4+NO3)/SRP), log(Kmorg))) + geom_point(cex=6, pch=21) + theme_min
#Kmf
ggplot(pars_all, aes(pH, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(Pnahco3, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(Porg, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(SRP, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(MBP, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(MBC, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(MBN, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(Ctot, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(Ntot, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(DOC, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(TN, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(DON, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(NH4, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(NO3, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(NH4+NO3, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(MBC/MBP, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(MBN/MBP, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(Porg/MBP, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(Porg/SRP, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(MBC/MBN, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(DOC/Porg, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(DOC/SRP, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(TN/Porg, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(TN/SRP, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(DON/Porg, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(DON/SRP, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(NH4/Porg, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(NH4/SRP, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(NO3/Porg, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(NO3/SRP, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes((NH4+NO3)/Porg, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes((NH4+NO3)/SRP, Kmf)) + geom_point(cex=6, pch=21) + theme_min
#Kic
ggplot(pars_all, aes(pH, Kic)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(Pnahco3, Kic)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(Porg, Kic)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(SRP, Kic)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(MBP, Kic)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(MBC, Kic)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(MBN, Kic)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(Ctot, Kic)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(Ntot, Kic)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(DOC, Kic)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(TN, Kic)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(DON, Kic)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(NH4, Kic)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(NO3, Kic)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(NH4+NO3, Kic)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(MBC/MBP, Kic)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(MBN/MBP, (Kic))) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(Porg/MBP, Kic)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(Porg/SRP, Kic)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(MBC/MBN, log(Kic))) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(DOC/Porg, Kic)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(DOC/SRP, Kic)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(TN/Porg, Kic)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(TN/SRP, Kic)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(DON/Porg, Kic)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(DON/SRP, Kic)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(NH4/Porg, Kic)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(NH4/SRP, Kic)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(NO3/Porg, Kic)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(NO3/SRP, Kic)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes((NH4+NO3)/Porg, Kic)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes((NH4+NO3)/SRP, Kic)) + geom_point(cex=6, pch=21) + theme_min

ggplot(pars_all, aes(TN/SRP, Vmax)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(MBP, Kmf)) + geom_point(cex=6, pch=21) + theme_min
ggplot(pars_all, aes(pH, Kic)) + geom_point(cex=6, pch=21, aes(fill=Legend)) + theme_min

grid.arrange(ggplot(pars_all, aes((NH4+NO3)/SRP, Vmax)) + geom_point(cex=6, aes(fill=Legend, shape=Legend), show.legend = F) + theme_min +
               scale_fill_manual(values = c("black", "black", "grey", "grey")) + 
               scale_shape_manual(values = c(21, 22, 21, 22)) +
               stat_smooth(method=lm, se=F, color="grey30") +
               ylim(0, 0.20) +
               theme(legend.title = element_blank(), legend.position = c(0.7, 0.3)) +
               ylab(expression(paste(V[MAX], " (",mu, "mol ", g^{-1},min^{-1}, ")" ))) +
               xlab(expression(paste("(N",H[4]^{"+"}," + N", O[3]^{"-"}, ")/SRP (mol/mol)"))) + ggtitle(("A)")),
             ggplot(pars_all, aes(MBP, Kmf)) + geom_point(cex=6, aes(fill=Legend, shape=Legend), show.legend = F) + theme_min +
               scale_fill_manual(values = c("black", "black", "grey", "grey")) + 
               scale_shape_manual(values = c(21, 22, 21, 22)) +
               stat_smooth(method=lm, se=F, color="grey30") +
               ylim(0, 16) + xlim(5, 15) +
               theme(legend.title = element_blank(), legend.position = c(0.7, 0.3)) +
               ylab(expression(paste(K[M-MUB-P], " (",mu, "mol ", g^{-1}, ")" ))) +
               xlab(expression(paste("MBP (", mu, "mol ", g^{-1}, ")"))) + ggtitle(("B)")),
             ggplot(pars_all, aes(MBC/MBN, Kic)) + geom_point(cex=6, aes(fill=Legend, shape=Legend), show.legend = T) + theme_min +
               scale_fill_manual(values = c("black", "black", "grey", "grey")) + 
               scale_shape_manual(values = c(21, 22, 21, 22)) +
               stat_smooth(method=lm, se=F, color="grey30") +
               scale_y_log10(limits=c(10, 200), breaks=c(10, 30, 70, 200))+
               theme(legend.title = element_blank(), legend.position = c(0.7, 0.75)) +
               ylab(expression(paste("Kic (",mu, "mol ", g^{-1}, ")" ))) +
               xlab("MBC/MBN (mol/mol)") + ggtitle(("C)")), nrow=3)
  


# ggplot(pars_all, aes(DOC/Porg, Kic)) + geom_point(cex=6, aes(fill=Legend, shape=Legend), show.legend = F) + theme_min +
#  scale_fill_manual(values = c("black", "black", "grey", "grey")) + 
#  scale_shape_manual(values = c(21, 22, 21, 22)) + scale_y_log10(limits=c(3, 20))+
#  theme(legend.title = element_blank(), legend.position = c(0.8, 0.2)) + ggtitle("B)") 


#############################################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SRP release kinetic~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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

ggplot(subset(pdata, outlier=="NO"), aes(Time, log(SRP_o))) + geom_point(cex=6, pch=21, aes(fill = as.factor(SRP_a2))) +
  facet_grid(.~Legend) + theme_min + #scale_fill_manual(values=c("black", "grey30", "grey60", "grey90", "white")) +
  stat_smooth(method = lm, se=F, aes(colour=as.factor(SRP_a2)), show.legend = F) +
  ylab(expression(paste("Measured SRP (", mu, "mol ", g^{-1}, ")"))) + xlab("Time (hours)") + 
  labs(fill=expression(paste("Added P-P",O[4], "(", mu, "mol ", g^{-1}, ")")))


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
                                    slope=coef(lm(log(SRP_o)~Time, data=pdata[(pdata$Catchment==i & pdata$Horizon==n & pdata$SRP_a==k & pdata$outlier=="NO"),]))[2],
                                    se=summary(lm(log(SRP_o)~Time, data=pdata[(pdata$Catchment==i & pdata$Horizon==n & pdata$SRP_a==k & pdata$outlier=="NO"),]))[[4]][[4]]
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

ggplot(SRPslopes, aes(SRP_a, slope*60)) + geom_point(cex=6, aes(fill=Legend, shape=Legend)) + theme_min +
  scale_fill_manual(values = c("black", "black", "grey", "grey")) + 
  scale_shape_manual(values = c(21, 22, 21, 22)) + geom_errorbar(aes(ymin=(slope-se)*60, ymax=(slope+se)*60))+
  theme(legend.title = element_blank(), legend.position = c(0.7, 0.7)) +
  ylab(expression(paste("SRP increase rate (", h^{-1}, ")"))) + 
  xlab(expression(paste("Added P-P", O[4]," (", mu, "mol ", g^{-1}, ")"))) +
  ggtitle("A)")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~End of the script~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#############################################################################################################################################
#############################################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Additional tests~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#############################################################################################################################################
#############################################################################################################################################
#Not run
#############################################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Testing the two pools of enzyme~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#############################################################################################################################################
# #Load the function
# source("two_epools.R")
# source("two_epools2.R")
# ##Plesne catchment
# ###Litter horizon
# plo_2ci<-two_epools(data=subset(edata, Catchment=="Plesne" & Horizon=="Litter" & time<120))
# plo_2ci$Parameters
# plo_2ci$Goodness$Gfit
# plo_ci2$Gfit
# plo_sp$Goodness$Gfit
# plo_cip$Goodness$Gfit
# plo_cie$Goodness$Gfit
# 
# ####logLike test
# #####CI vs SP
# -2*(plo_ci2$Gfit[["ll"]]-plo_sp$Goodness$Gfit[["ll"]])
# pchisq(-2*(plo_ci2$Gfit[["ll"]]-plo_sp$Goodness$Gfit[["ll"]]), df=1, lower.tail=F)
# #####CI vs CIP
# -2*(plo_ci2$Gfit[["ll"]]-plo_cip$Goodness$Gfit[["ll"]])
# pchisq(-2*(plo_ci2$Gfit[["ll"]]-plo_cip$Goodness$Gfit[["ll"]]), df=1, lower.tail=F)
# 
# ####additional calculations
# plo_2ci$Goodness<-two_epools2(data=subset(edata, Catchment=="Plesne" & Horizon=="Litter" & time<120), parameters = plo_2ci$Parameters)
# 
# ###Organic topsoil
# pla_2ci<-two_epools(data=subset(edata, Catchment=="Plesne" & Horizon=="Organic topsoil" & time<120))
# pla_2ci$Parameters
# pla_2ci$Goodness$Gfit
# pla_ci2$Gfit
# pla_sp$Goodness$Gfit
# pla_cip$Goodness$Gfit
# pla_cie$Goodness$Gfit
# 
# ####logLike test
# #####CI vs SP
# -2*(pla_ci2$Gfit[["ll"]]-pla_sp$Goodness$Gfit[["ll"]])
# pchisq(-2*(pla_ci2$Gfit[["ll"]]-pla_sp$Goodness$Gfit[["ll"]]), df=1, lower.tail=F)
# #####CI vs CIP
# -2*(pla_ci2$Gfit[["ll"]]-pla_cip$Goodness$Gfit[["ll"]])
# pchisq(-2*(pla_ci2$Gfit[["ll"]]-pla_cip$Goodness$Gfit[["ll"]]), df=1, lower.tail=F)
# 
# ####additional calculations
# pla_2ci$Goodness<-two_epools2(data=subset(edata, Catchment=="Plesne" & Horizon=="Organic topsoil" & time<120), parameters = pla_2ci$Parameters)
# 
# ##Certovo catchment
# ###Litter horizon
# co_2ci<-two_epools(data=subset(edata, Catchment=="Certovo" & Horizon=="Litter" & time<120))
# co_2ci$Parameters
# co_2ci$Goodness$Gfit
# co_ci2$Gfit
# co_sp$Goodness$Gfit
# co_cip$Goodness$Gfit
# co_cie$Goodness$Gfit
# 
# ####logLike test
# #####CI vs SP
# -2*(co_ci2$Gfit[["ll"]]-co_sp$Goodness$Gfit[["ll"]])
# pchisq(-2*(co_ci2$Gfit[["ll"]]-co_sp$Goodness$Gfit[["ll"]]), df=1, lower.tail=F)
# #####CI vs CIP
# -2*(co_ci2$Gfit[["ll"]]-co_cip$Goodness$Gfit[["ll"]])
# pchisq(-2*(co_ci2$Gfit[["ll"]]-co_cip$Goodness$Gfit[["ll"]]), df=1, lower.tail=F)
# 
# #####SP vs CIP
# -2*(co_sp$Goodness$Gfit[["ll"]]-co_cip$Goodness$Gfit[["ll"]])
# pchisq(-2*(co_sp$Goodness$Gfit[["ll"]]-co_cip$Goodness$Gfit[["ll"]]), df=1, lower.tail=F)
# 
# ####additional calculations
# co_2ci$Goodness<-two_epools2(data=subset(edata, Catchment=="Certovo" & Horizon=="Litter" & time<120), parameters = co_2ci$Parameters)
# 
# ###Organic topsoil
# ca_2ci<-two_epools(data=subset(edata, Catchment=="Certovo" & Horizon=="Organic topsoil" & time<120))
# ca_2ci$Parameters
# ca_2ci$Goodness$Gfit
# ca_ci2$Gfit
# ca_sp$Goodness$Gfit
# ca_cip$Goodness$Gfit
# ca_cie$Goodness$Gfit
# 
# ####logLike test
# #####CI vs SP
# -2*(ca_ci2$Gfit[["ll"]]-ca_sp$Goodness$Gfit[["ll"]])
# pchisq(-2*(ca_ci2$Gfit[["ll"]]-ca_sp$Goodness$Gfit[["ll"]]), df=1, lower.tail=F)
# #####CI vs CIP
# -2*(ca_ci2$Gfit[["ll"]]-ca_cip$Goodness$Gfit[["ll"]])
# pchisq(-2*(ca_ci2$Gfit[["ll"]]-ca_cip$Goodness$Gfit[["ll"]]), df=1, lower.tail=F)
# 
# ####additional calculations
# ca_2ci$Goodness<-two_epools2(data=subset(edata, Catchment=="Certovo" & Horizon=="Organic topsoil" & time<120), parameters = ca_2ci$Parameters)
# 
# 
# ###Add to a previous data frame
# epredts$v_two<-NA
# for(i in 1:nrow(epredts)){
#   if(epredts$Catchment[i]=="Plesne" & epredts$Horizon[i]=="Litter"){
#     epredts$v_two[i]<-plo_2ci$Parameters[1]*epredts$Substrate[i]/(plo_2ci$Parameters[2] + epredts$Substrate[i])+
#       plo_2ci$Parameters[3]*epredts$Substrate[i]/(plo_2ci$Parameters[4]*(1+epredts$InhibitorSRP[i]/plo_2ci$Parameters[5]) + epredts$Substrate[i])
#   }else{
#     if(epredts$Catchment[i]=="Plesne" & epredts$Horizon[i]=="Organic topsoil"){
#       epredts$v_two[i]<-pla_2ci$Parameters[1]*epredts$Substrate[i]/(pla_2ci$Parameters[2] + epredts$Substrate[i])+
#         pla_2ci$Parameters[3]*epredts$Substrate[i]/(pla_2ci$Parameters[4]*(1+epredts$InhibitorSRP[i]/pla_2ci$Parameters[5]) + epredts$Substrate[i])
#     }else{
#       if(epredts$Catchment[i]=="Certovo" & epredts$Horizon[i]=="Litter"){
#         epredts$v_two[i]<-co_2ci$Parameters[1]*epredts$Substrate[i]/(co_2ci$Parameters[2] + epredts$Substrate[i])+
#           co_2ci$Parameters[3]*epredts$Substrate[i]/(co_2ci$Parameters[4]*(1+epredts$InhibitorSRP[i]/co_2ci$Parameters[5]) + epredts$Substrate[i])
#       }else{
#         epredts$v_two[i]<-ca_2ci$Parameters[1]*epredts$Substrate[i]/(ca_2ci$Parameters[2] + epredts$Substrate[i])+
#           ca_2ci$Parameters[3]*epredts$Substrate[i]/(ca_2ci$Parameters[4]*(1+epredts$InhibitorSRP[i]/ca_2ci$Parameters[5]) + epredts$Substrate[i])
#       }
#     }
#   }
# }
# 
# ###relative rates
# epredts$v_tworel<-NA
# for(l in unique(epredts$Horizon)){
#   for(k in unique(epredts$Catchment)){
#     for(i in unique(epredts[(epredts$Horizon==l & epredts$Catchment==k), "Substrate"])){
#       for(n in unique(epredts[(epredts$Horizon==l & epredts$Catchment==k), "InhibitorSRP"])){
#         epredts[(epredts$Horizon==l & epredts$Catchment==k & epredts$Substrate==i & epredts$InhibitorSRP==n), "v_tworel"]<-
#           epredts[(epredts$Horizon==l & epredts$Catchment==k & epredts$Substrate==i & epredts$InhibitorSRP==n), "v_two"]/
#           epredts[(epredts$Horizon==l & epredts$Catchment==k & epredts$Substrate==i & epredts$InhibitorSRP==0), "v_two"]*100
#       }
#     }
#   }
# }
# 
# ggplot(erates, aes(InhibitorSRP, v5rel)) + geom_point(cex=6, pch=21, aes(fill = as.factor(Substrate3))) +
#   facet_grid(.~Legend) + theme_min + scale_fill_manual(values=c("black", "grey30", "grey60", "grey90", "white")) +
#   scale_y_continuous(limits = c(0, 110), breaks = c(0, 20, 40, 60, 80, 100)) +
#   xlim(0, 16) + geom_errorbar(aes(ymin=v5rel-v5rel.se, ymax=v5rel+v5rel.se), width=0.01) +
#   geom_line(data=epredts, aes(InhibitorSRP, v_tworel, color=as.factor(Substrate3), linetype=as.factor(Substrate3)), lwd=1.2) +
#   scale_color_manual(values=c("black", "grey30", "grey60", "grey90", "grey90")) +
#   scale_linetype_manual(values=c(rep("solid", 4), "longdash"))
# 
# #Combine all predictions
# Yhat_all<-rbind(plo_2ci$Goodness$Yhat, pla_2ci$Goodness$Yhat,
#                 co_2ci$Goodness$Yhat, ca_2ci$Goodness$Yhat)
# ##Vizualize
# ###add lumped factor for graphics
# for(i in 1:nrow(Yhat_all)){
#   if(Yhat_all$Catchment[i]=="Plesne" & Yhat_all$Horizon[i]=="Litter"){
#     Yhat_all$Legend[i]<-"Plesne - Litter"
#   }else{
#     if(Yhat_all$Catchment[i]=="Plesne" & Yhat_all$Horizon[i]=="Organic topsoil"){
#       Yhat_all$Legend[i]<-"Plesne - Organic topsoil"
#     }else{
#       if(Yhat_all$Catchment[i]=="Certovo" & Yhat_all$Horizon[i]=="Litter"){
#         Yhat_all$Legend[i]<-"Certovo - Litter"
#       }else{
#         Yhat_all$Legend[i]<-"Certovo - Organic topsoil"
#       }
#     }
#   }
# }
# 
# Yhat_all$InhibitorSRP2<-round(Yhat_all$InhibitorSRP, 0)
# Yhat_all$Substrate2<-round(Yhat_all$Substrate, 0)
# Yhat_all$InhibitorSRP3<-rep(Yhat_all[c(1:nrow(plo_2ci$Goodness$Yhat)), "InhibitorSRP2"], 4)
# Yhat_all$Substrate3<-rep(Yhat_all[c(1:nrow(plo_2ci$Goodness$Yhat)), "Substrate2"], 4)
# 
# ggplot(subset(Yhat_all, Substrate3==63 & InhibitorSRP3==4 & time<90), aes(time, Product)) +
#   geom_point(aes(shape=Legend, color=Legend), cex=6) + theme_min + 
#   geom_line(aes(time, Pred, color=Legend)) +
#   stat_smooth(mehtod=lm, aes(color=Legend), lty=2)
# 
#############################################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Testing the variable concentration of enzyme pool~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#############################################################################################################################################
# ##Load the function
# source("CIE.R")
# ##Plesne catchment
# ###Litter horizon
# plo_cie<-CIE(data=subset(edata, Catchment=="Plesne" & Horizon=="Litter" & time<120))
# plo_cie$Parameters
# plo_cie$Goodness$Gfit
# 
# ###Organic topsoil
# pla_cie<-CIE(data=subset(edata, Catchment=="Plesne" & Horizon=="Organic topsoil" & time<120))
# pla_cie$Parameters
# pla_cie$Goodness$Gfit
# 
# ##Certovo catchment
# ###Litter horizon
# co_cie<-CIE(data=subset(edata, Catchment=="Certovo" & Horizon=="Litter" & time<120))
# co_cie$Parameters
# co_cie$Goodness$Gfit
# 
# ###Organic topsoil
# ca_cie<-CIE(data=subset(edata, Catchment=="Certovo" & Horizon=="Organic topsoil" & time<120))
# ca_cie$Parameters
# ca_cie$Goodness$Gfit
# 
# #Add results to a epredts data frame
# epredts$v_e<-NA
# epredts$Porg<-NA
# 
# for(i in 1:nrow(epredts)){
#   if(epredts$Catchment[i]=="Plesne" & epredts$Horizon[i]=="Litter"){
#     epredts$v_e[i]<-plo_cie$Parameters[1]*plo_cie$Parameters[5]*epredts$Substrate[i]/
#       (plo_cie$Parameters[2]*(1+epredts$InhibitorSRP[i]/plo_cie$Parameters[3]) + epredts$Substrate[i])
#   }else{
#     if(epredts$Catchment[i]=="Plesne" & epredts$Horizon[i]=="Organic topsoil"){
#       epredts$v_e[i]<-pla_cie$Parameters[1]*pla_cie$Parameters[5]*epredts$Substrate[i]/
#         (pla_cie$Parameters[2]*(1+epredts$InhibitorSRP[i]/pla_cie$Parameters[3]) + epredts$Substrate[i])
#     }else{
#       if(epredts$Catchment[i]=="Certovo" & epredts$Horizon[i]=="Litter"){
#         epredts$v_e[i]<-co_cie$Parameters[1]*co_cie$Parameters[5]*epredts$Substrate[i]/
#           (co_cie$Parameters[2]*(1+epredts$InhibitorSRP[i]/co_cie$Parameters[3]) + epredts$Substrate[i])
#       }else{
#         epredts$v_e[i]<-ca_cie$Parameters[1]*ca_cie$Parameters[5]*epredts$Substrate[i]/
#           (ca_cie$Parameters[2]*(1+epredts$InhibitorSRP[i]/ca_cie$Parameters[3]) + epredts$Substrate[i])
#       }
#     }
#   }
# }
# 
# ###relative rates
# epredts$v_erel<-NA
# for(l in unique(epredts$Horizon)){
#   for(k in unique(epredts$Catchment)){
#     for(i in unique(epredts[(epredts$Horizon==l & epredts$Catchment==k), "Substrate"])){
#       for(n in unique(epredts[(epredts$Horizon==l & epredts$Catchment==k), "InhibitorSRP"])){
#         epredts[(epredts$Horizon==l & epredts$Catchment==k & epredts$Substrate==i & epredts$InhibitorSRP==n), "v_erel"]<-
#           epredts[(epredts$Horizon==l & epredts$Catchment==k & epredts$Substrate==i & epredts$InhibitorSRP==n), "v_e"]/
#           epredts[(epredts$Horizon==l & epredts$Catchment==k & epredts$Substrate==i & epredts$InhibitorSRP==0), "v_e"]*100
#       }
#     }
#   }
# }
# 
# ggplot(erates, aes(InhibitorSRP, v5rel)) + geom_point(cex=6, pch=21, aes(fill = as.factor(Substrate3))) +
#   facet_grid(.~Legend) + theme_min + scale_fill_manual(values=c("black", "grey30", "grey60", "grey90", "white")) +
#   scale_y_continuous(limits = c(0, 110), breaks = c(0, 20, 40, 60, 80, 100)) +
#   xlim(0, 16) + geom_errorbar(aes(ymin=v5rel-v5rel.se, ymax=v5rel+v5rel.se), width=0.01) +
#   geom_line(data=epredts, aes(InhibitorSRP, v_erel, color=as.factor(Substrate3), linetype=as.factor(Substrate3)), lwd=1.2) +
#   scale_color_manual(values=c("black", "grey30", "grey60", "grey90", "grey90")) +
#   scale_linetype_manual(values=c(rep("solid", 4), "longdash"))

# ggplot(erates, aes(InhibitorSRP, v90rel)) + geom_point(cex=6, pch=21, aes(fill = as.factor(Substrate3))) +
#   facet_grid(.~Legend) + theme_min + scale_fill_manual(values=c("black", "grey30", "grey60", "grey90", "white")) +
#   scale_y_continuous(limits = c(0, 110), breaks = c(0, 20, 40, 60, 80, 100)) +
#   xlim(0, 16) + geom_errorbar(aes(ymin=v90rel-v90rel.se, ymax=v90rel+v90rel.se), width=0.01) +
#   geom_line(data=epredts, aes(InhibitorSRP, v_erel, color=as.factor(Substrate3), linetype=as.factor(Substrate3)), lwd=1.2) +
#   scale_color_manual(values=c("black", "grey30", "grey60", "grey90", "grey90")) +
#   scale_linetype_manual(values=c(rep("solid", 4), "longdash"))

 
# #############################################################################################################################################
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Testing the two pools of enzyme with substrate preference~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# #############################################################################################################################################
# #Assumption: one enzyme pool is highly specific and inhibited by the product, fluoresecent substrate and Porg competes for enzymes
# #Load the function  
# source("two_eppools.R")
# ##Plesne catchment
# ###Litter horizon
# plo_2epci<-two_eppools(data=subset(edata, Catchment=="Plesne" & Horizon=="Litter" & time<90))
# plo_2epci$Parameters
# plo_2epci$Goodness$Gfit
# 
# ###Organic topsoil
# pla_2epci<-two_eppools(data=subset(edata, Catchment=="Plesne" & Horizon=="Organic topsoil" & time<90))
# pla_2epci$Parameters
# pla_2epci$Goodness$Gfit
# 
# ##Certovo catchment
# ###Litter horizon
# co_2epci<-two_eppools(data=subset(edata, Catchment=="Certovo" & Horizon=="Litter" & time<90))
# co_2epci$Parameters
# co_2epci$Goodness$Gfit
#  
# ###Organic topsoil
# ca_2epci<-two_eppools(data=subset(edata, Catchment=="Certovo" & Horizon=="Organic topsoil" & time<90))
# ca_2epci$Parameters
# ca_2epci$Goodness$Gfit
# 
# ###Add to a previous data frame
# epredts$v_eptwo<-NA
# for(i in 1:nrow(epredts)){
#  if(epredts$Catchment[i]=="Plesne" & epredts$Horizon[i]=="Litter"){
#    epredts$v_eptwo[i]<-plo_2epci$Parameters[1]*epredts$Substrate[i]/
#      (plo_2epci$Parameters[2]*(1+epredts$Porg[i]/plo_2epci$Parameters[3]) + epredts$Substrate[i])+
#      plo_2epci$Parameters[4]*epredts$Substrate[i]/
#      (plo_2epci$Parameters[5]*(1+epredts$Porg[i]/plo_2epci$Parameters[6])*
#      (1+epredts$InhibitorSRP[i]/plo_2epci$Parameters[7]) + epredts$Substrate[i])
#  }else{
#    if(epredts$Catchment[i]=="Plesne" & epredts$Horizon[i]=="Organic topsoil"){
#      epredts$v_eptwo[i]<-pla_2epci$Parameters[1]*epredts$Substrate[i]/
#        (pla_2epci$Parameters[2]*(1+epredts$Porg[i]/pla_2epci$Parameters[3]) + epredts$Substrate[i])+
#        pla_2epci$Parameters[4]*epredts$Substrate[i]/
#       (pla_2epci$Parameters[5]*(1+epredts$Porg[i]/pla_2epci$Parameters[6])*
#        (1+epredts$InhibitorSRP[i]/pla_2epci$Parameters[7]) + epredts$Substrate[i])
#    }else{
#      if(epredts$Catchment[i]=="Certovo" & epredts$Horizon[i]=="Litter"){
#        epredts$v_eptwo[i]<-co_2epci$Parameters[1]*epredts$Substrate[i]/
#          (co_2epci$Parameters[2]*(1+epredts$Porg[i]/co_2epci$Parameters[3]) + epredts$Substrate[i])+
#          co_2epci$Parameters[4]*epredts$Substrate[i]/
#          (co_2epci$Parameters[5]*(1+epredts$Porg[i]/co_2epci$Parameters[6]) *
#          (1+epredts$InhibitorSRP[i]/co_2epci$Parameters[7]) + epredts$Substrate[i])
#      }else{
#        epredts$v_eptwo[i]<-ca_2epci$Parameters[1]*epredts$Substrate[i]/
#          (ca_2epci$Parameters[2]*(1+epredts$Porg[i]/ca_2epci$Parameters[3]) + epredts$Substrate[i])+
#          ca_2epci$Parameters[4]*epredts$Substrate[i]/
#          (ca_2epci$Parameters[5]*(1+epredts$Porg[i]/ca_2epci$Parameters[6]) *
#          (1+epredts$InhibitorSRP[i]/ca_2epci$Parameters[7]) + epredts$Substrate[i])
#      }
#    }
#  }
# }
# 
# ###relative rates
# epredts$v_eptworel<-NA
# for(l in unique(epredts$Horizon)){
#  for(k in unique(epredts$Catchment)){
#    for(i in unique(epredts[(epredts$Horizon==l & epredts$Catchment==k), "Substrate"])){
#      for(n in unique(epredts[(epredts$Horizon==l & epredts$Catchment==k), "InhibitorSRP"])){
#        epredts[(epredts$Horizon==l & epredts$Catchment==k & epredts$Substrate==i & epredts$InhibitorSRP==n), "v_eptworel"]<-
#          epredts[(epredts$Horizon==l & epredts$Catchment==k & epredts$Substrate==i & epredts$InhibitorSRP==n), "v_eptwo"]/
#          epredts[(epredts$Horizon==l & epredts$Catchment==k & epredts$Substrate==i & epredts$InhibitorSRP==0), "v_eptwo"]*100
#      }
#    }
#  }
# }
# 
# ggplot(erates, aes(InhibitorSRP, v5rel)) + geom_point(cex=6, pch=21, aes(fill = as.factor(Substrate3))) +
#  facet_grid(.~Legend) + theme_min + scale_fill_manual(values=c("black", "grey30", "grey60", "grey90", "white")) +
#  scale_y_continuous(limits = c(0, 110), breaks = c(0, 20, 40, 60, 80, 100)) +
#  xlim(0, 16) + geom_errorbar(aes(ymin=v5rel-v5rel.se, ymax=v5rel+v5rel.se), width=0.01) +
#  geom_line(data=epredts, aes(InhibitorSRP, v_eptworel, color=as.factor(Substrate3), linetype=as.factor(Substrate3)), lwd=1.2) +
#  scale_color_manual(values=c("black", "grey30", "grey60", "grey90", "grey90")) +
#  scale_linetype_manual(values=c(rep("solid", 4), "longdash"))
# 
# source("CI_hill.R")
# ##Plesne catchment
# ###Litter horizon
# plo_hill<-CI_hill(data=subset(edata, Catchment=="Plesne" & Horizon=="Litter" & time<120))
# plo_hill$Parameters
# plo_hill$Goodness$Gfit
# 
# ###Organic topsoil
# pla_hill<-CI_hill(data=subset(edata, Catchment=="Plesne" & Horizon=="Organic topsoil" & time<120))
# pla_hill$Parameters
# pla_hill$Goodness$Gfit
# 
# ##Certovo catchment
# ###Litter horizon
# co_hill<-CI_hill(data=subset(edata, Catchment=="Certovo" & Horizon=="Litter" & time<120))
# co_hill$Parameters
# co_hill$Goodness$Gfit
# 
# ###Organic topsoil
# ca_hill<-CI_hill(data=subset(edata, Catchment=="Certovo" & Horizon=="Organic topsoil" & time<120))
# ca_hill$Parameters
# ca_hill$Goodness$Gfit
# 
# source("CI_SI.R")
# ##Plesne catchment
# ###Litter horizon
# plo_si<-CI_SI(data=subset(edata, Catchment=="Plesne" & Horizon=="Litter" & time<120))
# plo_si$Parameters
# plo_si$Goodness$Gfit
# 
# ###Organic topsoil
# pla_si<-CI_SI(data=subset(edata, Catchment=="Plesne" & Horizon=="Organic topsoil" & time<120))
# pla_si$Parameters
# pla_si$Goodness$Gfit
# 
# ##Certovo catchment
# ###Litter horizon
# co_si<-CI_SI(data=subset(edata, Catchment=="Certovo" & Horizon=="Litter" & time<120))
# co_si$Parameters
# co_si$Goodness$Gfit
# 
# ###Organic topsoil
# ca_si<-CI_SI(data=subset(edata, Catchment=="Certovo" & Horizon=="Organic topsoil" & time<120))
# ca_si$Parameters
# ca_si$Goodness$Gfit
