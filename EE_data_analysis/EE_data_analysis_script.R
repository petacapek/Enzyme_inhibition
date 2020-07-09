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

################################################################Raw data#######################################################################
#loading the script to recalculate raw data
source("../EE_data_analysis/E_calc.R")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Plesne catchment~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Litter horizon 
##Measured on 14.4.2020
###uploading the file with reaction time
po_time<-read_ods(path="../Raw_data//14.4.2020/casy14.4.2020.ods", sheet = 2, col_names=TRUE)
po_all<-E_calc(dataset = c("/mnt/580CBE2464C5F83D/pracovni/data_statistika/Junior/Enzyme_inhibition/Raw_data/14.4.2020/MUB_14.4.2020.xlsx"),
                 MUBconc = c(0, 5, 25, 50, 125, 250),
                 APconc = c(50, 100, 250, 500, 1000),
                 Iconc = c(0, 20, 80, 120, 240),
                 Nmeasure = 25, empty = 1,
                 Times = as.numeric(po_time[po_time$Ctr=="NE", "Time"]))
###Extract the data
po<-po_all$data
###Do the correction against the control
po$Pcorr<-po$P-pmax(0, po$Pcontr)
###Remove negative values
po$Pcorr2<-pmax(0, as.numeric(po$Pcorr))

###All units are umol/L
###recalculate to umol/g(DW) - fresh soil:buffer - 1:100 (w/w)
po$Product<-po$Pcorr2/10/0.3195
po$Substrate<-po$C_AP/10/0.3195
po$Inhibitor<-po$C_I/10/0.3195

###Add the sample description
po$Catchment<-"Plesne"
po$Horizon<-"Litter"

#Organic topsoil horizon
###Measured 24.3.2020
###uploading the file with reaction time
pa_time<-read_ods(path="../Raw_data/24.3.2020/casy24.3.2020.ods", sheet = 2, col_names=TRUE)
pa_all<-E_calc(dataset = c("/mnt/580CBE2464C5F83D/pracovni/data_statistika/Junior/Enzyme_inhibition/Raw_data/24.3.2020/MUB_24.3.2020.xlsx"),
                 MUBconc = c(0, 5, 25, 50, 125, 250),
                 APconc = c(50, 100, 250, 500, 1000),
                 Iconc = c(0, 20, 80, 120, 240),
                 Nmeasure = 25, empty = 1,
                 Times = as.numeric(pa_time[pa_time$Ctr=="NE", "Time"]))
###Extract the data
pa<-pa_all$data
###Do the correction against the control
pa$Pcorr<-pa$P-pmax(0, pa$Pcontr)
###Remove negative values
pa$Pcorr2<-pmax(0, as.numeric(pa$Pcorr))

###All units are umol/L
###recalculate to umol/g(DW) - fresh soil:buffer - 1:100 (w/w)
pa$Product<-pa$Pcorr2/10/0.2992
pa$Substrate<-pa$C_AP/10/0.2992
pa$Inhibitor<-pa$C_I/10/0.2992

pa$Catchment<-"Plesne"
pa$Horizon<-"Organic topsoil"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Certovo catchment~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Litter horizon 
##Measured on 6.4.2020
###uploading the file with reaction time
co_time<-read_ods(path="../Raw_data/6.4.2020/casy6.4.2020.ods", sheet = 2, col_names=TRUE)
co_all<-E_calc(dataset = c("/mnt/580CBE2464C5F83D/pracovni/data_statistika/Junior/Enzyme_inhibition/Raw_data/6.4.2020/MUB_6.4.2020.xlsx"),
                 MUBconc = c(0, 5, 25, 50, 125, 250),
                 APconc = c(50, 100, 250, 500, 1000),
                 Iconc = c(0, 20, 80, 120, 240),
                 Nmeasure = 25, empty = 1,
                 Times = as.numeric(co_time[co_time$Ctr=="NE", "Time"]))

###Extract the data
co<-co_all$data
###Do the correction against the control
co$Pcorr<-co$P-pmax(0, co$Pcontr)
###Remove negative values
co$Pcorr2<-pmax(0, as.numeric(co$Pcorr))

###All units are umol/L
###recalculate to umol/g(DW) - fresh soil:buffer - 1:100 (w/w)
co$Product<-co$Pcorr2/10/0.2911
co$Substrate<-co$C_AP/10/0.2911
co$Inhibitor<-co$C_I/10/0.2911

###Add the sample description
co$Catchment<-"Certovo"
co$Horizon<-"Litter"

#Organic topsoil horizon 
##Measured on 30.3.2020
###uploading the file with reaction time
ca_time<-read_ods(path="../Raw_data/30.3.2020/casy30.3.2020.ods", sheet = 2, col_names=TRUE)
ca_all<-E_calc(dataset = c("/mnt/580CBE2464C5F83D/pracovni/data_statistika/Junior/Enzyme_inhibition/Raw_data/30.3.2020/MUB_30.3.2020.xlsx"),
                 MUBconc = c(0, 5, 25, 50, 125, 250),
                 APconc = c(50, 100, 250, 500, 1000),
                 Iconc = c(0, 20, 80, 120, 240),
                 Nmeasure = 25, empty = 1,
                 Times = as.numeric(ca_time[ca_time$Ctr=="NE", "Time"]))
###Extract the data
ca<-ca_all$data
###Do the correction against the control
ca$Pcorr<-ca$P-pmax(0, ca$Pcontr)
###Remove negative values
ca$Pcorr2<-pmax(0, as.numeric(ca$Pcorr))

###All units are umol/L
###recalculate to umol/g(DW) - fresh soil:buffer - 1:100 (w/w)
ca$Product<-ca$Pcorr2/10/0.3122
ca$Substrate<-ca$C_AP/10/0.3122
ca$Inhibitor<-ca$C_I/10/0.3122

###Add the sample description
ca$Catchment<-"Certovo"
ca$Horizon<-"Organic topsoil"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Merge and export the data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
edata<-rbind(po,pa,co,ca)
write.csv(edata, file = c("/mnt/580CBE2464C5F83D/pracovni/data_statistika/Junior/Enzyme_inhibition/enzyme_data.csv"))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#