setwd("E:/Rworkspace/ChildhoodPneumonia")#change your url
library(dplyr)
library(survival)
library(plyr)
library(MASS)
library(forestplot)
library(pheatmap)
##########
# uniCOX #
##########
data<-read.csv("./Temps/2.finalData.csv",stringsAsFactors = F)
Cox_sx<-read.csv("./Results/1.dataBeforeCox.csv",stringsAsFactors = F)
UniCox_Data<-Cox_sx[,-c(1,4,5)]
UniCox_Data$GENDER<-as.factor(UniCox_Data$GENDER)
surv<-Surv(time = Cox_sx$Time,event = Cox_sx$EXPIRE_FLAG)
UniCox_Data$surv=surv
#unicox function
UniCox<-function(x){
  FML<-as.formula(paste0("surv~",x))
  Cox<-coxph(FML,data = UniCox_Data) 
  Sum<-summary(Cox)
  CI<-paste0(round(Sum$conf.int[,3:4],3),collapse = "-") 
  Pvalue<-round(Sum$coefficients[,5],3)
  HR<-round(Sum$coefficients[,2],3)
  Unicox<-data.frame("Characteristics"=x,
                     "Hazard Ratio"=HR,
                     "CI95"=CI,
                     "P value"=Pvalue)
  return(Unicox)
}

varNames<-colnames(UniCox_Data)[1:87] 
length(varNames)
UniVar<-list()
for (nm in varNames){
  temp <- tryCatch(
    { UniCox(nm) },
    error = function(e) { message('Error @ ',nm) ; return(NA) },
    finally = { message('next...') 
    }
  )
  UniVar[[nm]] <- temp
}

UniVar<-ldply(UniVar,data.frame)
UniVar<-UniVar[,2:5]
UniVar<-na.omit(UniVar)
write.csv(UniVar,"./Results/3.unicoxResults.csv",row.names = F)

Cox_Item<-UniVar$Characteristics[UniVar$P.value<0.05]
Cox_data<-data[,c(which(colnames(data) %in% Cox_Item))]
Cox_data<-Cox_data %>%filter(!if_all(.fns = is.na))
p_id<-data$SUBJECT_ID[as.numeric(rownames(Cox_data))] 

Cox_all<-cbind(SUBJECT_ID=data$SUBJECT_ID, Time=data$Time, EXPIRE_FLAG=data$EXPIRE_FLAG)
Cox_all<-as.data.frame(Cox_all[which(Cox_all[,1]%in% p_id),]) 
heat_fac <- Cox_data %>% select_if(is.character)
heat_num <- Cox_data %>% select_if(is.numeric)
Cox_all<-cbind(Cox_all,heat_fac,heat_num)
table(Cox_all$EXPIRE_FLAG)

##############
#  mutiCOX   #
##############
UniCox_var<-cbind( Time=Cox_all$Time, EXPIRE_FLAG=Cox_all$EXPIRE_FLAG,
                   X5409.Hardness=Cox_all[,4],
                   Cox_all[,c(5:27)])
UniCox_var<-na.omit(UniCox_var)
items<-colnames(UniCox_var)[-(1:2)]
s<-paste(items,collapse="+")
FML<-as.formula(paste0("Surv(Time, EXPIRE_FLAG)~",s))
res.cox <- coxph(FML, data = UniCox_var)
mucox<-summary(res.cox)

formatFit<-function(fit){
  p_val<-summary(fit)$coefficients[,5]
  wald<-summary(fit)$coefficients[,4]^2
  valueB<-summary(fit)$coefficients[,1]
  valueOR<-summary(fit)$coefficients[,2]
  confitOR<-exp(confint(fit))
  data.frame(
    OR_with_CI=paste(round(valueOR,3),"(",
                     round(confitOR[,1],3),"-",round(confitOR[,2],3),")",sep=""),
    P=round(p_val,3),
    coef=round(valueB,3),
    Wald=round(wald,3),
    exp.coef=round(valueOR,3),
    lower=round(confitOR[,1],3),
    upper=round(confitOR[,2],3)
  )
}
summary(res.cox)$coefficients
data_f=formatFit(res.cox)
write.csv(data_f,"./Temps/3. forestplot_muti.csv")#your need to add a title

index=rownames(data_f)[which(data_f$`Pr(>|z|)`<0.05)] 
df=data[,c(1:6,which(colnames(data)%in%index[-1]))]

x_RowsAllNa_removed =  df[apply(df[7:11], 1, function(y) any(!is.na(y))),]
write.csv(x_RowsAllNa_removed,"./Results/3.modelData.csv",row.names = F)


##########################
# forestplot figure 2B  #
#########################
data=read.csv("./Temps/2. forestplot_muti.csv")
f<-forestplot(as.matrix(data[,c(1:3)]),data$exp.coef,data$lower,
              data$upper,zero = 1,xlog = F,
              clip = c(0,5),
              colgap = unit(5,"mm"),graphwidth=unit(60,"mm"),
              lineheight = unit(0.8,"cm"),graph.pos = 3,
              col = fpColors(box="blue", lines="blue", zero = "grey"),boxsize = 0.3,
              ci.vertices = T,ci.vertices.height = 0.2,
              lty.ci = 7,lwd.zero=0.4,lwd.ci = 3,
              txt_gp=fpTxtGp(label = gpar(cex=0.8),
                             ticks = gpar(cex=0.8)) ,
              is.summary=c(TRUE,rep(FALSE,100))
)
f
pdf("./Figs/2. forestPlot.pdf",width=10, height=12)
f
dev.off()

#####################
# heatmap figure 2A #
#####################
Cox_all<-na.omit(Cox_all)
Cox_all<-Cox_all[order(Cox_all$EXPIRE_FLAG),]
heat_cli=Cox_all[,3:4]#factor
heat_exp=Cox_all[,5:27]#num
annotation_col = data.frame(Survival=as.factor(heat_cli[,1]),
                            Hardness = as.factor(heat_cli[,2]))
summary(annotation_col)
rownames(annotation_col) <- Cox_all$SUBJECT_ID
ann_colors = list(
  Survival=c("0"="#2E94B9","1"="#de4307"),
  Hardness = c( "Pasty"="#A6CEE3" ,"Soft"="#1F78B4","Thin"="#B2DF8A"))
red <- rgb(255,0,0,maxColorValue = 255)
blue <- rgb(0,0,255,maxColorValue = 255)
white <- rgb(255,255,255,maxColorValue = 255)
linshi <- scale(heat_exp[,1:23])
rownames(linshi) <- Cox_all$SUBJECT_ID
linshi<-t(linshi)
hist(linshi)
linshi[linshi>2] <- 2
linshi[linshi<(-2)] <- c(-2)
p<-pheatmap(linshi,fontsize=6,cutree_col = 4,#cellheight = 2,cellwidth = 1 ,
            color  = colorRampPalette(c(blue,white,red))(100),
            annotation_col = annotation_col,
            annotation_colors = ann_colors,
            clustering_method = "ward.D2",
            border_color = "grey60",
            cluster_cols = F, cluster_rows = F,
            show_rownames = T, show_colnames = F
)
p
pdf("./Figs/2. heatmap.pdf")
p
dev.off()

