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
score.cox<-coxph(Surv(Time, survival)~score,data = coxData)
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
# > summary(coxph(Surv(Time,survival)~fenzu,data = coxData))
# Call:
#   coxph(formula = Surv(Time, survival) ~ fenzu, data = coxData)
# 
# n= 631, number of events= 95 
# 
# coef exp(coef) se(coef)     z Pr(>|z|)    
# fenzuHigh 1.5143    4.5461   0.2206 6.866 6.61e-12 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# fenzuHigh     4.546       0.22     2.951     7.004
# 
# Concordance= 0.647  (se = 0.028 )
# Likelihood ratio test= 37.73  on 1 df,   p=8e-10
# Wald test            = 47.14  on 1 df,   p=7e-12
# Score (logrank) test = 56.72  on 1 df,   p=5e-14
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
regplot(Coxfit, plots=c("density","boxes"), observation=T, failtime=c(7,30,100),
        title="Risk Nomogram", prfail=T, clickable=TRUE, points=TRUE)

#####################################
#  validate data  3marker remodeled #
#####################################
vali.data=read_excel("./Data/validata/zjyu.xlsx",sheet = 1)#ALL DATA
y.vali=data.frame(X5033.Urea=vali.data$`Blood urea nitrogen (mmol/L)`,
                  X5099.Hemoglobin=vali.data$`Hemoglobin (g/dL)`,
                  X5252.Oxygen.Saturation=vali.data$`Oxygen saturation`*100,
                  survival=vali.data$`Outcome(1=survival,0=death)`,#这里标错了，1是死亡
                  Time=vali.data$`LOS-Hospitalization time(days)`)
vdata1 <- data.matrix(y.vali[,-c(4:5)]) 
vali_label<-y.vali$survival
vdata2 <- Matrix(vdata1,sparse=T)  ## 
vdata3 <- factor(vali_label,levels = c(0,1))   ### 
vdata4 <- list(data=vdata2,label=vali_label)  ### candidate training data

dvali <- xgb.DMatrix(data = vdata4$data, label = vdata4$label)
#use original model parameters
best_vparams<-as.list(c(eta=0.2,max_depth=5,min_child_weight=1,
                        subsample=0.8,colsample_bytree=1,
                        objective="binary:logistic",eval_metric="logloss"))
# Fit final model
set.seed(1)
xgbcv <- xgb.cv(params = best_vparams,
                data = dvali,
                nround = 100,
                nfold = 5,
                prediction = TRUE,
                early_stopping_rounds = 5,
                verbose = 0,
                maximize = F)

vnrounds = xgbcv$best_iteration


#final model
set.seed(1)
mxgb_vali <- xgboost(data = dvali, 
                     params = best_vparams, 
                     maximize = F, 
                     nrounds = vnrounds,
                     early_stopping_rounds = 5, 
                     verbose = 0)

# Evaluate on validate set
vali_pred <- predict(mxgb_vali,vdata1,type = "response")

roc(vali_label,vali_pred)
## [1] 0.909
aa<-roc(vali_label,vali_pred, plot=TRUE, print.thres=TRUE, ci=TRUE,
        print.auc=TRUE,legacy.axes = TRUE,col = "#1687A7")#cutoff=0.498

sp.obj <- ci.sp(aa, sensitivities=seq(0, 1, .01), boot.n=100) 
plot(sp.obj, type="shape", col="gray60") 

roc(vali_label,y.vali$X5033.Urea)
roc(vali_label,y.vali$X5099.Hemoglobin)
roc(vali_label,y.vali$X5252.Oxygen.Saturation)


new.vali=as.data.frame(vdata1)
new.vali$survival=y.vali$survival
new.vali$Time=y.vali$Time
new.vali$score=vali_pred
new.vali$fenzu="Low"
new.vali$fenzu[which(new.vali$score>0.498)]="High"

tcolor=c("#f1ac9d","#1c96a6")
p1<-ggplot(new.vali, aes(x = reorder(rownames(new.vali),score), y = score-0.502)) +
  geom_col(aes(fill = as.factor(fenzu)))+
  scale_fill_manual(values =tcolor)+
  xlab("")+ylab("Score")+
  theme_bw() + theme(panel.grid=element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x=element_blank())

sdf<-data.frame(id=reorder(rownames(new.vali),new.vali$score),
                time=new.vali$Time,survival=as.factor(new.vali$survival),
                fenzu=new.vali$fenzu,
                score=new.vali$score)
scolor<-c("#ffa200","#128fce")
p2<-ggscatter(sdf, x = "id", y = "time", color = "survival",shape = 21,
              palette = scolor)+
  xlab("id")+ylab("Time (day)")+
  theme_bw() + theme(panel.grid=element_blank())


heat_exp=new.vali[,1:3]
heat_exp<-heat_exp[order(new.vali$score),]
red <- rgb(255,0,0,maxColorValue = 255)
blue <- rgb(0,0,255,maxColorValue = 255)
white <- rgb(255,255,255,maxColorValue = 255)
linshi <- scale(heat_exp[,1:3])#归一化
linshi<-t(linshi)
hist(linshi)
linshi[linshi>2] <- 2
linshi[linshi<(-2)] <- -2

p3<-pheatmap(linshi,fontsize=6,cutree_col = 4,cellheight = 6,cellwidth = 2,
             color  = colorRampPalette(c(blue,white,red))(100),
             clustering_method = "ward.D2",
             border_color = "grey60",
             cluster_cols = F, cluster_rows = F,
             show_rownames = T, show_colnames = F
)
p4 <- plot_grid(p1,p2,p3$gtable, ncol = 1, align = "v")
p4
######KM########
red="#BC3C28"
blue="#003399"
fit <- survfit(Surv(Time,survival)~fenzu,data = new.vali)
testSur<-ggsurvplot(fit,
                    pval = TRUE,#P-value
                    risk.table = TRUE, # Add risk table
                    risk.table.col = "strata", # Change risk table color by groups
                    palette = c(red,blue),
                    conf.int = T,
                    conf.int.alpha=0.1, 
                    surv.median.line = "hv", 
                    xlim=c(0,100),
                    break.x.by = 20,
                    legend.labs=c("High","Low")
)
# HR
summary(coxph(Surv(new.vali$Time,new.vali$survival) ~ vali_pred))
#1403[ 149.8     13141]

pdf("./Figs/4. valiSur-KM.pdf",height=6,width=5)
testSur
dev.off()
#####diff boxplot#######
pdata=melt(new.vali,
           id.vars='fenzu', 
           measure.vars=colnames(new.vali)[1:3],
           variable.name = "variety", 
           value.name = "value")

#pdata=data.frame(variety,fenzu,value)
pdata$fenzu<-factor(pdata$fenzu,levels=c("High","Low"))
testDiff=ggplot(pdata, aes(x=variety, y=value, fill=fenzu)) + 
                stat_compare_means(label = "p.format",method = "wilcox.test")+
                geom_boxplot(outlier.shape = NA) +
                scale_fill_manual(values = c("#D94E48","#5175A4"))+
                xlab("")+ylab("")+
                facet_wrap(~variety, scale="free",nrow=1)
pdf("./Figs/4. boxplot.pdf",height=3,width=5)
testDiff
dev.off()
