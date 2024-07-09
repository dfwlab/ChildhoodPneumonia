setwd("E:/Rworkspace/ChildhoodPneumonia")#change your url
library(rms)
library(pROC)
library(e1071)
library(randomForest)
library(nnet)
library(xgboost)
library(Matrix)
library(ggsci)
library(scales)
##########logistics#######
com.df=na.omit(mdata)
colnames(mdata)
log_mod <- glm(survival~.,data=train[,-2])
log_mod
log_pred <- predict(log_mod,test[,3:7],type = "response")
lg_roc=roc(test$survival,log_pred)#0.8043

#### -----SVM-------######
svm_model <- svm(as.factor(survival) ~ ., data = train[,-2],probability=T,
                 kernel = "radial")
tune_result <- tune(svm, as.factor(survival) ~ ., 
                    data = train[,-2], kernel = "radial",probability=T,
                    ranges = list(cost = c(0.1, 1, 10), gamma = c(0.1, 1, 10)))
svm_mod <- tune_result$best.model
summary(svm_mod)
svm_pred=predict(svm_mod,test[,3:7],probability=T)
pred_prob=attr(svm_pred,"probabilities")
na_rows <- !complete.cases(test)
svm_roc=roc(test$survival[-which(na_rows)],pred_prob[,1])#0.6197
#######--Random forest--########
rf_mod <- randomForest(as.factor(survival) ~ ., data =na.omit(train[,-2]))
rf_pred<- predict(rf_mod, test[,3:7],"prob")[,2]
rf_roc=roc(test$survival, rf_pred) # 0.8224

########--KNN--#######
nn_model <- nnet(survival ~ ., 
                 data = train[,-2], size = 5, decay = 0.01, maxit = 200, trace = FALSE)

nn_probs <- predict(nn_model, test[,3:7], type = "raw")

knn_roc=roc(test$survival, nn_probs) #0.6057

########--xgboost--#####
params <- list(
  booster = "gbtree",
  objective = "binary:logistic",
  eval_metric = "auc"
)

xgb_mod <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 10,
  verbose = 1
)
pred <- predict(xgb_mod, test_data)
xgb_roc=roc(test$survival, pred)#0.8357

######ROC plot########
mypal<-pal_lancet("lanonc",alpha=0.6)(9)
show_col(mypal)

plot(lg_roc, col = mypal[1], main = "ROC Curves for Multiple Models")
plot(knn_roc, col = mypal[2], add = TRUE)
plot(rf_roc, col = mypal[3], add = TRUE)
plot(svm_roc, col = mypal[4], add = TRUE)
plot(xgb_roc, col = mypal[5], add = TRUE)

legend("bottomright", legend = c("logistic: AUC=0.804",
                                 "Neural Network: AUC=0.606",
                                 "random forest: AUC=0.822",
                                 "SVM: AUC=0.620", 
                                 "XGBoost: AUC=0.836"),
       col = mypal[1:5], lty = 1)
