setwd("E:/Rworkspace/ChildhoodPneumonia")#change your url
library(readxl)
library(Matrix)
library(verification)
library(pROC)
library(ggplot2)
library(ggpubr)
library(caret)
library(DiagrammeR)
library(tidyverse)
library(pacman)
library(cowplot)
library(grid)
library(survival)
library(survminer)
library(regplot)
library(pheatmap)
p_load(xgboost, ParBayesianOptimization, mlbench, dplyr, skimr, recipes, resample)

################################
#        markers  validation   #
################################
data <- read.csv("./Results/3.modelData.csv",stringsAsFactors = F)
###########-----cox--------##########
summary(data)
coxData<-data
res.cox <- coxph(Surv(Time, survival)~., data = coxData)
mucox<-summary(res.cox)
A<-as.matrix(res.cox$coefficients,ncol=1) 
B=as.matrix(data[,3:7])
score<-B%*%A
coxData$score=score[,1]
score.cox<-coxph(Surv(Time, survival)~scale(score),data = coxData)
summary(score.cox)


res.cut <- surv_cutpoint(coxData, #data
                         time = "Time", 
                         event = "survival", 
                         variables = "score") 

summary(res.cut) 
coxData$fenzu<-"ALow"
coxData$fenzu[which(coxData$score>(-1.521))]<-"High"
table(coxData$fenzu,coxData$survival)
coxData$score<-coxData$score+1.521

tcolor=c("#f1ac9d","#1c96a6")
z1<-ggplot(coxData, aes(x = reorder(rownames(coxData),score), y = score)) +
  geom_col(aes(fill = as.factor(fenzu)))+
  scale_fill_manual(values =tcolor)+
  xlab("")+ylab("Score")+
  theme_bw() + theme(panel.grid=element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x=element_blank())
#dot
sdf<-data.frame(id=reorder(rownames(coxData),score),time=coxData$Time,survival=as.factor(coxData$survival),fenzu=coxData$fenzu)
scolor<-c("#ffa200","#128fce")
z2<-ggscatter(sdf, x = "id", y = "time", color = "survival",shape = 21,
              palette = scolor)+
  ylim(0,100)+
  xlab("")+ylab("Time (day)")+
  theme_bw() + theme(panel.grid=element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x=element_blank())

z3 <- plot_grid(z1,z2, ncol = 1, align = "v")
z3

#survival plot 
fit <- survfit(Surv(Time,survival)~fenzu,data = coxData)
summary(coxph(Surv(Time,survival)~fenzu,data = coxData))
# coef exp(coef) se(coef)     z Pr(>|z|)    
# fenzuHigh 1.5143    4.5461   0.2206 6.866 6.61e-12 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# fenzuHigh     4.546       0.22     2.951     7.004

tcolor=c("#1c96a6","#f1ac9d")
ggsurvplot(fit,
           pval = TRUE,#P-value
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           palette = tcolor,
           conf.int = T,
           conf.int.alpha=0.1, #显示置信区间
           surv.median.line = "hv", # 增加中位生存时间
           xlim=c(0,100),
           break.x.by = 10,
           legend.labs=c("Low","High"))

#nomogram
FML<-as.formula(Surv(Time,survival==1)~X5083.Uric.Acid..Urine+X5033.Urea+
                  X5252.Oxygen.Saturation+X5099.Hemoglobin+X6473.Lipase)
Coxfit<-coxph(FML,data=coxData)
regplot(Coxfit, plots=c("density","boxes"), observation=T, failtime=c(7,30),
        title="Risk Nomogram", prfail=T, clickable=TRUE, points=TRUE)

#####################################
#  validation by 2019 PICU dataset  #
#####################################
library(yardstick)
library(ggplot2)
library(tinyarray)
library(vegan)
######## Batch effect #########
#If needed, please email the corresponding author to request the data set.
patients=read_xlsx("./data/pnue19.xlsx")
summary(patients)
vali=as.data.frame(lapply(patients[,c(2,5:11)], function(x) as.numeric(as.character(x))))
vali$gender=patients$Gender
y.vali=na.omit(vali)

colnames(coxData)
colnames(y.vali)
alldf=rbind(coxData[,3:7],y.vali[,4:8])
alldf$group=c(rep("PICData",631),rep("validata",129))
p=draw_pca(exp = t(alldf[,1:5]), group_list = factor(alldf$group))
pca_score= p$data$coord
panosim <- anosim(pca_score, factor(alldf$group), permutations=999)
summary(panosim)
# ANOSIM statistic R: -0.01933 
# Significance: 0.787 

######## model validation ##########
median(y.vali$Age_month)#8[4-21]
IQR = quantile(y.vali$Age_month, 0.75) - quantile(y.vali$Age_month, 0.25)
table(y.vali$survival)
table(y.vali$gender)
validata=as.matrix(y.vali[,4:8]) 
vali_pred <- predict(mxgb5m,validata,type="response")

xgboost_roc=roc(y.vali$survival,vali_pred)
aa<-roc(y.vali$survival,vali_pred, plot=TRUE, print.thres=TRUE, ci=TRUE,
        print.auc=TRUE,legacy.axes = TRUE,col = "#1687A7")

sp.obj <- ci.sp(aa, sensitivities=seq(0, 1, .01), boot.n=100) 
plot(sp.obj, type="shape", col="gray60") 

pred=numeric(length(vali_pred)) 
pred[which(vali_pred>0.206)]=1

##########################
#   results Figure plot  #
##########################
library(survival)
library(survminer)
library(pheatmap)
library(grid)
library(cowplot)
#######--split group--#########
y.vali$score=vali_pred
y.vali$fenzu="Low"
y.vali$fenzu[which(y.vali$score>0.206)]="High"
#bar plot
tcolor=c("#f1ac9d","#1c96a6")
p1<-ggplot(y.vali, aes(x = reorder(rownames(y.vali),score-0.206), y = score-0.206)) +
  geom_col(aes(fill = as.factor(fenzu)))+
  scale_fill_manual(values =tcolor)+
  xlab("")+ylab("Score")+
  theme_bw() + theme(panel.grid=element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x=element_blank())
#scatter plot
sdf<-data.frame(id=reorder(rownames(y.vali),y.vali$score),
                time=y.vali$Time,survival=as.factor(y.vali$survival),
                fenzu=y.vali$fenzu,
                score=y.vali$score)
scolor<-c("#ffa200","#128fce")
p2<-ggscatter(sdf, x = "id", y = "time", color = "survival",shape = 21,
              palette = scolor)+
  xlab("id")+ylab("Time (day)")+
  theme_bw() + theme(panel.grid=element_blank())

#heatmap
y.vali<-y.vali[order(y.vali$score),]
red <- rgb(255,0,0,maxColorValue = 255)
blue <- rgb(0,0,255,maxColorValue = 255)
white <- rgb(255,255,255,maxColorValue = 255)
summary(y.vali)
heat_cli=y.vali[,c(10,2)]#risk and sur
heat_cli<-heat_cli[order(y.vali$score),]
summary(heat_cli)
heat_exp=y.vali[,c(4:8)]
heat_exp<-heat_exp[order(y.vali$score),]
summary(heat_exp)
annotation_col = data.frame(Risk=as.factor(heat_cli[,1]),
                            survival=as.factor(heat_cli[,2]))
summary(annotation_col)
rownames(annotation_col) <- rownames(y.vali)
ann_colors = list(Risk=c("Low"="#003399","High"="#BC3C28"),
                  survival=c("0"="#2E94B9","1"="#de4307"))

linshi <- scale(heat_exp[,1:5])
rownames(linshi) <- rownames(y.vali)
linshi<-t(linshi)
hist(linshi)
linshi[linshi>2] <- 2
linshi[linshi<(-2)] <- c(-2)
p3<-pheatmap(linshi,fontsize=6,cutree_col = 4,cellheight = 5,cellwidth = 0.5 ,
             color  = colorRampPalette(c(blue,white,red))(100),
             annotation_col = annotation_col,
             annotation_colors = ann_colors,
             clustering_method = "ward.D2",
             border_color = "grey60",
             cluster_cols = F, cluster_rows = F,
             show_rownames = T, show_colnames = F
)

score=y.vali$score[order(y.vali$score)]
score<-t(score)
p4<-pheatmap(score,fontsize=6,cutree_col = 4,cellheight = 5,cellwidth = 0.5 ,
             color  = colorRampPalette(c("#dadada","#0B0B0B"))(100),
             clustering_method = "ward.D2",
             border_color = "grey60",
             cluster_cols = F, cluster_rows = F,
             show_rownames = T, show_colnames = T
)

p_vali <- plot_grid(p4$gtable,p3$gtable,p1,p2, ncol = 1, align = "v")


########## boxplot ############
data2 <- data.frame(score = y.vali$score, 
                    cluster = as.factor(y.vali$fenzu)) 
a2<-ggplot(data = data2, aes(x = cluster , y = score,fill=cluster)) + 
  stat_compare_means(label.y = 1, label.x = 1)  +  
  geom_violin(trim=F,color="white",width = 0.9) +
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.1,fill = "#FFFFFF")+ 
  scale_fill_manual(values = tcolor)+ 
  theme_bw()+ 
  theme(axis.text.x=element_text(hjust = 0,colour="black",size=10), 
        axis.text.y=element_text(size=10,colour="black"), 
        axis.title.y=element_text(size = 10,colour="black"), 
        legend.text=element_text(colour="black",  
                                 size=10),
        legend.title=element_text(colour="black", 
                                  size=10),
        panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())+ 
  ylab("Predict score")


######  KM  ########
fit <- survfit(Surv(Time,survival)~fenzu,data = y.vali)
testSur<-ggsurvplot(fit,
                    pval = TRUE,#P-value
                    risk.table = TRUE, # Add risk table
                    risk.table.col = "strata", # Change risk table color by groups
                    palette = c("#f1ac9d","#1c96a6"),
                    conf.int = T,
                    conf.int.alpha=0.1, #显示置信区间
                    surv.median.line = "hv", # 增加中位生存时间
                    xlim=c(0,80),
                    break.x.by = 20,
                    legend.labs=c("High","Low")
)
# HR
summary(coxph(Surv(y.vali$Time,y.vali$survival) ~ scale(y.vali$score)))
#P=3.06e-06 *** 
# exp(coef) exp(-coef) lower .95 upper .95
# scale(y.vali$score)     2.107     0.4745     1.541     2.882
