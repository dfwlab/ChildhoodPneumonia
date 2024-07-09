setwd("E:/Rworkspace/ChildhoodPneumonia")#change your url
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(tidyverse)
library(plyr)

patients<-read.csv("./PICdata/Pneu_patients.csv")
#####################
#  lab data filter  #
#####################
Lab<-read.csv("./PICdata/LABEVENTS.csv")
Lab_item<-read.csv("./PICdata/D_LABITEMS.csv")
#----24h data--#
#patient Lab data
Lab_pneu_patient<-Lab%>% filter(Lab$SUBJECT_ID %in% patients$SUBJECT_ID)
# order
Lab_pneu_patient<-Lab_pneu_patient[order(Lab_pneu_patient$SUBJECT_ID,
                                         Lab_pneu_patient$ITEMID,
                                         Lab_pneu_patient$CHARTTIME),]
#24h selected
setDT(Lab_pneu_patient)
index<-duplicated(Lab_pneu_patient,by=c("SUBJECT_ID","ITEMID"))
first_lab<-Lab_pneu_patient[which(index==F),]
first_lab<-na.omit(first_lab)

#######################
# data cleaning       #
#######################
lab_index<-unique(first_lab$ITEMID)
lab_px<-table(first_lab$ITEMID)
lab_px<-lab_px[order(lab_px,decreasing = T)]

#----AGE distribution Figure s1---#
sdf<-data.frame(Age=fData$ADMI_AGE,time=fData$Time/30,survival=as.factor(fData$EXPIRE_FLAG))
sdf<-sdf[which(sdf$Age/365 <5),]#under 5 years old
scolor<-c("#00b9f1","#f9320c")
p=ggscatter(sdf, x = "Age", y = "time",color = "survival",shape = 21,
            palette = scolor)+
  xlab("Age")+ylab("Time (month)")+
  scale_x_continuous(breaks=seq(0, 365*5, 365))+
  theme_bw() + theme(panel.grid=element_blank())
pdf("./Figs/s. age-time.pdf",width=10, height=6)
p
dev.off()


#----select age<=3 year-------#
patient<-patients[which(patients$AGE_TYPE=="Infancy"|patients$AGE_TYPE=="Toddler"),] #749 patients
patient<-patient[order(patient$SUBJECT_ID),]
length(patient$SUBJECT_ID)
write.csv(patient[,2:17],"./Temps/1.PatientsData.csv",row.names = F)

#---lab data filter--#
Lab_pneu_patient<-Lab%>% filter(Lab$SUBJECT_ID %in% patient$SUBJECT_ID)
Lab_pneu_patient<-Lab_pneu_patient[order(Lab_pneu_patient$SUBJECT_ID,
                                         Lab_pneu_patient$ITEMID,
                                         Lab_pneu_patient$CHARTTIME),]
setDT(Lab_pneu_patient)
index<-duplicated(Lab_pneu_patient,by=c("SUBJECT_ID","ITEMID"))
first_lab<-Lab_pneu_patient[which(index==F),]
length(unique(first_lab$SUBJECT_ID)) #patients num 633  have lab marker

lab_index<-unique(first_lab$ITEMID)
lab_px<-table(first_lab$ITEMID)
lab_px<-lab_px[order(lab_px,decreasing = T)]
# delete lab marker which missing>20%
lab_useful<-as.data.frame(lab_px[which(lab_px/633>=0.8)])
Lab_select<-lab_useful$Var1 
#delete wrong data (for a given feature, unite not same)
funits<-as.data.frame(table(first_lab$ITEMID,first_lab$VALUEUOM)) 
funits<-funits[-which(funits$Freq==0),]
funits<-funits[-which(funits$Var2==""),]
a<-as.data.frame(table(funits$Var1))
wunits<-a$Var1[which(a$Freq>1)] 
Lab_select<-Lab_select[-which(Lab_select%in%wunits)]
fName<-Lab_item[which(Lab_item$ITEMID%in%Lab_select),]
length(fName$ITEMID)# 85 features
write.csv(fName[,2:7],"./Temps/1.featureName.csv",row.names = F)


#-----final lab data------#
finalData<-patient%>%select(SUBJECT_ID,HADM_ID,ADMI_AGE,GENDER,Time,EXPIRE_FLAG)
finalData<-finalData[order(finalData$SUBJECT_ID),]
finalData<-finalData %>% filter(finalData$SUBJECT_ID %in% unique(first_lab$SUBJECT_ID))
pipei<-first_lab[which(first_lab$SUBJECT_ID%in%patient$SUBJECT_ID),c(2,4,6)]
pipei<-pipei[which(pipei$ITEMID%in%Lab_select),]
data<-dcast(pipei,SUBJECT_ID~ITEMID)
finalData<-join(finalData,data,match="first")
dim(finalData)
colnames(finalData)[7:ncol(finalData)]<-paste(fName$ITEMID,fName$LABEL,sep = ".")
head(finalData)
write.csv(finalData,
        file = "./Temps/1.labDataSelected.csv",
        row.names = F,fileEncoding = "UTF-8")

