setwd("E:/Rworkspace/ChildhoodPneumonia")#change your url
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(tidyverse)
library(plyr)
#################
#  Finial data  #
#################
labdata<-read.csv("./Temps/1.labDataSelected.csv") 
#factor cols
factor_item<-Filter(function(u) any(c('Negative','Positive','Clear','Pasty','Soft','Yellow','type O') %in% u), labdata) #选择包含negative的列
factor_item$GENDER <-labdata$GENDER
#ransform chrs
which.nonnum <- function(x) {
  which(is.na(suppressWarnings(as.numeric(as.character(x)))))
}
colIndex_fac<-which(colnames(labdata) %in% colnames(factor_item))
colIndex_num<- which(!colnames(labdata) %in% colnames(factor_item))
for (i in colIndex_num) {  #num col but have chr value
  labdata[which.nonnum(labdata[,i]),i]<-NA
}

labdata[,colIndex_num]<-as.data.frame(lapply(labdata[,colIndex_num],as.numeric))
labdata[,colIndex_fac]<-as.data.frame(lapply(labdata[,colIndex_fac],as.factor))
fac_data<-as.data.frame(lapply(labdata[,colIndex_fac],as.factor))
test1<-fac_data %>% select_if(is.numeric)
test2<-fac_data %>% select_if(is.character)
labdata[,which(colnames(labdata) %in% colnames(test1))]<-test1
labdata[,which(colnames(labdata) %in% colnames(test2))]<-test2
write.csv(labdata,
          file = "./Temps/2.finalData.csv",
          row.names = F,fileEncoding = "UTF-8")#then check by hand

data<-read.csv("./Temps/2.finalData.csv") 
test<-data %>% select_if(is.character)
colIndex_fac<-which(colnames(data) %in% colnames(test))
fac_data<-as.data.frame(lapply(data[,colIndex_fac],as.factor))
fac_data<-fac_data[,-1]
num_data<-data[,-colIndex_fac]
num_data<-num_data[,-(1:5)]
length(fac_data)#15 factor indicators
length(num_data)#70 numeric indicators
Cox_sx<-cbind(SUBJECT_ID=data$SUBJECT_ID,ADMI_AGE=data$ADMI_AGE,
              GENDER=data$GENDER,Time=data$Time,EXPIRE_FLAG=data$EXPIRE_FLAG,
              fac_data,num_data)
write.csv(Cox_sx,"./Results/1.dataBeforeCox.csv",row.names = F)

##############################
# data distrubution Table s1 #
##############################
labdata[,colIndex_num]<-as.data.frame(lapply(labdata[,colIndex_num],as.numeric))
labdata[,colIndex_fac]<-as.data.frame(lapply(labdata[,colIndex_fac],as.factor))
fac_data<-as.data.frame(lapply(labdata[,colIndex_fac],as.factor))
test1<-labdata %>% select_if(is.numeric)
test2<-labdata %>% select_if(is.factor)
test2$EXPIRE_FLAG<-labdata$EXPIRE_FLAG
##alive data
a<-summary(test1[which(test1$EXPIRE_FLAG==0),])
a<-as.matrix.data.frame(a,nrow=7,ncol=75)
colnames(a)<-colnames(test1)
rownames(a)<-c("Min", "1st Qu", "Median", "Mean", "3rd Qu", "Max","NA")
a<-t(a)
summary(test2[which(test1$EXPIRE_FLAG==0),])
write.csv(a,"./Results/2.1NumFeatureDistribution0.csv")

##death data
b<-summary(test1[which(test1$EXPIRE_FLAG==1),])
b<-as.matrix.data.frame(b,nrow=7)
colnames(b)<-colnames(test1)
rownames(b)<-c("Min", "1st Qu", "Median", "Mean", "3rd Qu", "Max","NA")
b<-t(b)
summary(test2[which(test1$EXPIRE_FLAG==1),])
write.csv(b,"./Results/2.1NumFeatureDistribution1.csv")