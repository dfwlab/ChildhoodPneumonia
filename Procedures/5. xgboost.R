setwd("E:/Rworkspace/ChildhoodPneumonia")#change your url
library(xgboost)
library(Matrix)
library(verification)
library(pROC)
library(ggplot2)
library(ggpubr)
library(caret)
library(DiagrammeR)
library(tidyverse)
library(pacman)
library(SHAPforxgboost)
library(dcurves)

data <- read.csv("./Results/3.modelData.csv",stringsAsFactors = F)
dead <- data[which(data$survival  == 1),]
alive <- data[which(data$survival  == 0),]
######################
# split train & test #
######################
set.seed(1)
dead_random <-sample(1:95,76)
alive_random <-sample(1:536,429)

train <- rbind(dead[dead_random,],alive[alive_random,])
test <- rbind(dead[-dead_random,],alive[-alive_random,])
train_data <- as.matrix(train[,-(1:2)]) 
test_data <- as.matrix(test[,-(1:2)])## Xgboost training
train_label <- train$survival  ### sample clinical information
traindata1 <- data.matrix(train[,-(1:2)]) 
traindata2 <- Matrix(traindata1,sparse=T)  ## 
traindata3 <- factor(train_label,levels = c(0,1))   ### 
traindata4 <- list(data=traindata2,label=train_label)  ### candidate training data

dtrain <- xgb.DMatrix(data = traindata4$data, label = traindata4$label)

##############
# best param #
##############
grd <- expand.grid(
  eta = seq(0.001, 0.2, length.out = 5),
  max_depth = seq(2L, 5L, by = 1),
  min_child_weight = seq(1, 25, length.out = 3),
  subsample = c(0.2,0.6,0.8),
  colsample_bytree=c(0.25, 0.5, 1))

dim(grd)#540 5

grd_out <- apply(grd, 1, function(par){
  
  par <- append(par, list(objective = "binary:logistic",eval_metric = "logloss"))
  mdl <- xgboost(data = dtrain, params = par, nrounds = 100, early_stopping_rounds = 5, maximize = F, verbose = 0)
  lst <- data.frame(par, score = mdl$best_score)
  
  return(lst)
})
grd_out <- do.call(rbind, grd_out)

best_par <- grd_out %>%
  data.frame() %>%
  arrange(score) %>%
  .[1,]
best_par
# eta	max_depth	min_child_weight	subsample	colsample_bytree	objective	eval_metric	score
# 0.2	   5	     1	               0.8	      1	       binary:logistic	logloss	0.06350403

best_params<-as.list(c(eta=0.2,max_depth=5,min_child_weight=1,
                       subsample=0.8,colsample_bytree=1,
                       objective="binary:logistic",eval_metric="logloss"))

#################
# xgboost model #
#################
set.seed(978)
xgbcv <- xgb.cv(params = best_params,
                data = dtrain,
                nround = 100,
                nfold = 5,
                prediction = TRUE,
                early_stopping_rounds = 5,
                verbose = 0,
                maximize = F)

nround = xgbcv$best_iteration
nround#14

#final model
set.seed(978)
mxgb5m <- xgboost(data = dtrain, 
                  params = best_params, 
                  maximize = F, 
                  nrounds = nround,
                  early_stopping_rounds = 5, 
                  verbose = 0)


#####################
# model performance #
#####################
#importance
importance_matrix <- xgb.importance(colnames(dtrain), model = mxgb5m)
xgb.plot.importance(importance_matrix, rel_to_first = TRUE, xlab = "Importance")

#shap
shap_values <- shap.values(xgb_model = mxgb5m, X_train = train_data)
shap_contrib <- shap_values$shap_score
shap_long <- shap.prep(xgb_model = mxgb5m, X_train = train_data, shap_contrib = shap_contrib)
shap.plot.summary(shap_long)

######---AUC-----#######
#train
train_pred <- predict(mxgb5m,train_data,type = "response")
roc(train_label,train_pred)# 0.9425
aa<-roc(train_label,train_pred, plot=TRUE, print.thres=TRUE, ci=TRUE,
        print.auc=TRUE,legacy.axes = TRUE,col = "#1687A7")#cutoff=0.206
sp.obj <- ci.sp(bb, sensitivities=seq(0, 1, .01), boot.n=100) 
plot(sp.obj, type="shape", col="gray60") #add ci

#test
test_pred <- predict(mxgb5m,test_data,type = "response")
roc(test$survival,test_pred)#0.8763
bb<-roc(test$survival,test_pred, plot=TRUE, print.thres=TRUE, ci=TRUE,
        print.auc=TRUE,legacy.axes = TRUE,col = "#1687A7")
sp.obj <- ci.sp(bb, sensitivities=seq(0, 1, .01), boot.n=100) 
plot(sp.obj, type="shape", col="gray60") 

#######---DCA&CIC---######
dca_data <- data.frame(pred = c(train_pred,test_pred), 
                       outcome = c(train_label,test$survival))
dca_result <- decision_curve(outcome ~ pred, data = dca_data, 
                             thresholds = seq(0, 1, by = 0.01))
#DCA
plot_decision_curve(dca_result,
                    curve.names = "model",
                    confidence.intervals = F)
#CIC
plot_clinical_impact(dca_result,col=c('#1687A7','#E64B35FF'))


##########################
#   results Figure plot  #
##########################
library(survival)
library(survminer)
library(pheatmap)
library(grid)
library(cowplot)
#######--train--#########
train$fenzu<-NA
train$fenzu[which(train_pred>0.206)]<-"High"
train$fenzu[which(train_pred<=0.206)]<-"Low"
blue <- "#003399"
red <- "#BC3C28"
train_surv <- survfit(Surv(train$Time,train$survival)~train$fenzu,data = train)
survdiff(Surv(train$Time,train$survival)~train$fenzu,data = train)
ggsurvplot(train_surv,
           pval = TRUE,#P-value
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           palette = c(red,blue),
           conf.int = T,
           conf.int.alpha=0.1, #CI
           surv.median.line = "hv", # median survival time
           xlim=c(0,100),
           break.x.by = 10,
           legend.labs=c("High","Low")
)
# HR
summary(coxph(Surv(train$Time,train$survival) ~ train_pred))  ## HR 4202   [1431-12337]

#heatmap
red <- rgb(255,0,0,maxColorValue = 255)
blue <- rgb(0,0,255,maxColorValue = 255)
white <- rgb(255,255,255,maxColorValue = 255)
summary(train)
train$score=train_pred
heat_cli=train[,c(8,1)]#risk and sur
heat_cli<-heat_cli[order(train$score),]
summary(heat_cli)
heat_exp=train[,c(3:7)]
heat_exp<-heat_exp[order(train$score),]
summary(heat_exp)
annotation_col = data.frame(Risk=as.factor(heat_cli[,1]),
                            Survival=as.factor(heat_cli[,2]))
summary(annotation_col)
rownames(annotation_col) <- rep(1:nrow(train))
ann_colors = list(Risk=c("Low"="#003399","High"="#BC3C28"),
                  Survival=c("0"="#2E94B9","1"="#de4307"))

linshi <- scale(heat_exp[,1:5])#归一化
rownames(linshi) <- rep(1:nrow(train))
linshi<-t(linshi)
hist(linshi)
linshi[linshi>2] <- 2
linshi[linshi<(-2)] <- c(-2)
p1<-pheatmap(linshi,fontsize=6,cutree_col = 4,cellheight = 5,cellwidth = 0.5 ,
             color  = colorRampPalette(c(blue,white,red))(100),
             annotation_col = annotation_col,
             annotation_colors = ann_colors,
             clustering_method = "ward.D2",
             border_color = "grey60",
             cluster_cols = F, cluster_rows = F,
             show_rownames = T, show_colnames = F
)
score=train$score[order(train$score)]
score<-t(score)
p2<-pheatmap(score,fontsize=6,cutree_col = 4,cellheight = 5,cellwidth = 0.5 ,
             color  = colorRampPalette(c("#dadada","#0B0B0B"))(100),
             clustering_method = "ward.D2",
             border_color = "grey60",
             cluster_cols = F, cluster_rows = F,
             show_rownames = T, show_colnames = T
)

p_train <- plot_grid(p2$gtable,p1$gtable, nrow = 2, align = "v")

pdf("./Figs/3. heatmapTrain.pdf")
p_train
dev.off()
#----boxplot compare High and Low risk
data2 <- data.frame(score = train_pred, 
                    cluster = as.factor(train$fenzu)) 
a2<-ggplot(data = data2, aes(x = cluster , y = score,fill=cluster)) + 
  stat_compare_means(label.y = 1, label.x = 1)  +  
  geom_violin(trim=F,color="white",width = 0.9) +
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.1,fill = "#FFFFFF")+ 
  scale_fill_manual(values = c("#BC3C28","#003399"))+ 
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

#----boxplot compare High and Low risk 5 markers--#
# 5 markers
marker_traindata<-data.frame(marker=rep(c("Uric Acid","Urea","Oxygen Saturation","Hemoglobin","Lipase"),each=nrow(train)),
                             value=c(train$X5083.Uric.Acid..Urine,train$X5033.Urea,train$X5252.Oxygen.Saturation,
                                     train$X5099.Hemoglobin,train$X6473.Lipase),
                             cluster=as.factor(rep(train$fenzu,5))) 
summary(marker_traindata)
a3<-ggplot(marker_traindata, aes(x=marker, y=value, fill=cluster)) + 
  stat_compare_means(label = "p.format",method = "wilcox.test")+
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("#BC3C28","#003399"))+
  xlab("")+ylab("")+
  facet_wrap(~marker, scale="free",nrow=1)+
  theme_bw() + theme(panel.grid=element_blank())

grid.newpage()  ##new plotpage
pushViewport(viewport(layout = grid.layout(2,2))) 
vplayout <- function(x,y){ viewport(layout.pos.row = x, layout.pos.col = y)}
print(NULL, vp = vplayout(1,1))   
print(a2, vp = vplayout(1,2))   
print(a3, vp = vplayout(2,1:2))  

###########################
#       test KM           #
###########################
test$fenzu<-NA
test$fenzu[which(test_pred>0.206)]<-"High"
test$fenzu[which(test_pred<=0.206)]<-"Low"
test_surv <- survfit(Surv(test$Time,test$survival)~test$fenzu,data = test)
survdiff(Surv(test$Time,test$survival)~test$fenzu,data = test)
blue <- "#003399"
red <- "#BC3C28"
ggsurvplot(test_surv,
           pval = TRUE,#P-value
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           palette = c(red,blue),
           conf.int = T,
           conf.int.alpha=0.1, 
           surv.median.line = "hv", 
           xlim=c(0,100),
           break.x.by = 10,
           legend.labs=c("High","Low")
)
# HR
summary(coxph(Surv(test$Time,test$survival) ~ test_pred))  

#-------heatmap---------#
summary(test)
test$score=test_pred
heat_cli=test[,c(8,1)]#risk and sur
heat_cli<-heat_cli[order(test$score),]
summary(heat_cli)
heat_exp=test[,c(3:7)]
heat_exp<-heat_exp[order(test$score),]
summary(heat_exp)
annotation_col = data.frame(Risk=as.factor(heat_cli[,1]),
                            Survival=as.factor(heat_cli[,2]))
summary(annotation_col)
rownames(annotation_col) <- rep(1:nrow(test))

linshi <- scale(heat_exp[,1:5])
rownames(linshi) <- rep(1:nrow(test))
linshi<-t(linshi)
hist(linshi)
linshi[linshi>2] <- 2
linshi[linshi<(-2)] <- c(-2)
p3<-pheatmap(linshi,fontsize=6,cutree_col = 4,cellheight = 5,cellwidth = 0.5,
             color  = colorRampPalette(c(blue,white,red))(100),
             annotation_col = annotation_col,
             annotation_colors = ann_colors,
             clustering_method = "ward.D2",
             border_color = "grey60",
             cluster_cols = F, cluster_rows = F,
             show_rownames = T, show_colnames = F
)
score=test$score[order(test$score)]
score<-t(score)
p4<-pheatmap(score,fontsize=6,cutree_col = 4,cellheight = 5,cellwidth =0.5,
             color  = colorRampPalette(c("#dadada","#0B0B0B"))(100),
             clustering_method = "ward.D2",
             border_color = "grey60",
             cluster_cols = F, cluster_rows = F,
             show_rownames = T, show_colnames = T
)

p_test<- plot_grid(p4$gtable,p3$gtable, nrow = 2, align = "v")

pdf("./Figs/3. heatmapTest.pdf")
p_test
dev.off()

#----boxplot compare High and Low risk
data3 <- data.frame(score = test_pred, 
                    cluster = as.factor(test$fenzu)) 
a4<-ggplot(data = data3, aes(x = cluster , y = score,fill=cluster)) + 
  stat_compare_means(label.y = 1, label.x = 1)  +  
  geom_violin(trim=F,color="white",width = 0.9) +
  geom_boxplot(outlier.size=NA, size=0.01, outlier.shape = NA,width = 0.1,fill = "#FFFFFF")+ 
  scale_fill_manual(values = c("#BC3C28","#003399"))+ 
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

#----boxplot compare High and Low risk 5 markers--#
marker_testdata<-data.frame(marker=rep(c("Uric Acid","Urea","Oxygen Saturation","Hemoglobin","Lipase"),each=nrow(test)),
                            value=c(test$X5083.Uric.Acid..Urine,test$X5033.Urea,test$X5252.Oxygen.Saturation,
                                    test$X5099.Hemoglobin,test$X6473.Lipase),
                            cluster=as.factor(rep(test$fenzu,5))) 
a5<-ggplot(marker_testdata, aes(x=marker, y=value, fill=cluster)) + 
  stat_compare_means(label = "p.format",method = "wilcox.test")+
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("#BC3C28","#003399"))+
  xlab("")+ylab("")+
  facet_wrap(~marker, scale="free",nrow=1)+
  theme_bw() + theme(panel.grid=element_blank())

grid.newpage()  
pushViewport(viewport(layout = grid.layout(2,2))) 
vplayout <- function(x,y){ viewport(layout.pos.row = x, layout.pos.col = y)}
print(NULL, vp = vplayout(1,1))   
print(a4, vp = vplayout(1,2))   
print(a5, vp = vplayout(2,1:2)) 

