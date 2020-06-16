wdbc <- read.csv("wdbc.csv", header = T)
features <- c("radius", "texture", "perimeter", "area", "smoothness", "compactness",
              "concavity", "concave_points", "symmetry", "fractal_dimension")
names(wdbc) <- c("id", "diagnosis", paste0(features,"_mean"), paste0(features,"_se"),
                 paste0(features,"_worst"))
wdbc.data <- wdbc[,c(3:32)]
row.names(wdbc.data) <- wdbc$id
wdbc_raw <- cbind(wdbc.data, as.factor(as.numeric(wdbc$diagnosis)-1))
colnames(wdbc_raw)[31] <- "diagnosis"

####################################################################
library(C50)
library(caret)
library(MASS)
library(ggplot2)
library(RColorBrewer)
library(ROCR)
library(tidyverse)
library(kernlab)
library(class)
library(boot)

class <- read.csv("Classification.csv", header = TRUE)

# Summaries
summary(class)
par(mfrow=c(1,2))
hist(class$X1[class$Group == 1],main="Group 1",xlab="")
hist(class$X2[class$Group == 0],main="Group 0",xlab="")

# Splitting data - split into groups based on even response group membership 
ii <- createDataPartition(class$Group, p=.8, list=F)  
xTrain <- class[ii,1:2]; yTrain <- class[ii,3]
xTest <- class[-ii,1:2]; yTest <- class[-ii,3]

################## Linear discriminant analysis #####################

# Similar approach to PCA as it finds linear combinations of variables in the data that best explain the data, 
# but instead of transforming the axes to minimise the variance, LDA also 
# finds the axes which maximise the separation between classes - generally favoured where there are more
# than two response classes
# It produces a 'discriminant score' for each class, and assigns to the group which results in the largest score.
# When there are more than 1 predictors, the observations are assumed to be from a multivariate gaussian distribution,
# so the class specific variance is changed to a covariance function that is common to all classes, with class-specific mean vectors

lda <- lda(yTrain~X1+X2,data=xTrain)
lda
# Preditions
lda.pred <- predict(lda, newdata = cbind(xTest,yTest))
# Evaluating
confusionMatrix(as.factor(lda.pred$class),as.factor(yTest))
# Quite poor predictions - accuracy is ok but sensitivity (true positives) very low despite 
# specificity (true negatives) being high. Poor pos predictive value also suggests the model not
# good at predicting group 1. 

twoClassColor <- brewer.pal(3,'Set1')[1:2]
names(twoClassColor) <- c('0','1')

# Decision boundary
nbp <- 250;
PredA <- seq(min(class$X1), max(class$X1), length = nbp)
PredB <- seq(min(class$X2), max(class$X2), length = nbp)
Grid <- expand.grid(X1 = PredA, X2 = PredB)
ggplot(data = class, aes(x = X1, y = X2,
                            color = as.factor(Group))) +
  geom_contour(data = cbind(Grid, Group = predict(lda,Grid)$class),
  aes(z = as.numeric(Group)), color = "red", breaks = c(1.5),size=1) + geom_point(size=2,alpha=0.8) +
  ggtitle("LDA Decision Boundary") + theme(legend.text = element_text(size = 10)) + 
  scale_colour_manual(name = 'Group', values = twoClassColor) + theme_bw()

# It is observable that the variance of group 1 in respect of the X2 variable, meaning the quadratic element
# of the QDA should be able to take into account this variance

# ROC
prediction(lda.pred$posterior[,2], yTest) %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot(main="LDA: AUC 0.732")
# AUC
prediction(lda.pred$posterior[,2], yTest) %>%
  performance(measure = "auc") %>%
  .@y.values


################## Quadratic discriminant analysis #####################

# LDA assumed that each class has the same covariance - so each class has its own covariance matrix unlike LDA, 
# so \Sigma^{-1} in the expression becomes \Sigma_k^{-1} means that the boundaries become quadratic
var(class$X2[class$Group==0])
var(class$X2[class$Group==1]) # Large disparity especially in the X2 explanatory variable
# 

qda <- qda(yTrain~X1+X2,data=xTrain)
qda
# Preditions
qda.pred <- predict(qda, newdata = cbind(xTest,yTest))
# Evaluating
confusionMatrix(as.factor(qda.pred$class),as.factor(yTest))

# Decision boundary
ggplot(data = class, aes(x = X1, y = X2,
                         color = as.factor(Group))) +
  geom_contour(data = cbind(Grid, Group = predict(qda,Grid)$class),
               aes(z = as.numeric(Group)), color = "red", breaks = c(1.5),size=1) + geom_point(size=2,alpha=0.8) +
  ggtitle("QDA Decision Boundary") + theme(legend.text = element_text(size = 10)) + 
  scale_colour_manual(name = 'Group', values = twoClassColor) + theme_bw()

# ROC
prediction(qda.pred$posterior[,2], yTest) %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot(main="QDA: AUC 0.830")
# AUC
prediction(qda.pred$posterior[,2], yTest) %>%
  performance(measure = "auc") %>%
  .@y.values

# On the basis of each of these measures, it is clear that quadratic discriminant analysis is more 
# appropriate for the classifications 

################# Logistic regression ###################################

# Generally can outperform the LDA and QDA methods when the covariances are not Gaussian, which is an 
# assumption of the previous methods. Rather than calculating a discriminant score which calculates straight
# away which group the class should be included in, logistic regression is a form of linear regression where 
# the link function is on the logit scale and classifies based on a cut-off probability like 0.5 when there
# are two explanatory variables. Predictions on the reponse variable level as converts log-odds to probabilities 

log.reg <- glm(yTrain~X1+X2,data=xTrain,family = "binomial")
a <- predict(log.reg,newdata=cbind(xTest,yTest),type = "response")
b <- ifelse(a>0.5,1,0)
confusionMatrix(as.factor(b),as.factor(yTest))

# Decision boundary 
ggplot(data = class, aes(x = X1, y = X2,
                         color = as.factor(Group))) +
  geom_abline(intercept=(coef(log.reg)[1]/(-coef(log.reg)[3])),slope=(coef(log.reg)[2]/(-coef(log.reg)[3])),
              color = "red",size=1) + 
  geom_point(size=2,alpha=0.8) +
  ggtitle("LR Decision Boundary") + theme(legend.text = element_text(size = 10)) + 
  scale_colour_manual(name = 'Group', values = twoClassColor) + theme_bw()

# ROC
prediction(a, yTest) %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot()
# AUC
prediction(a, yTest) %>%
  performance(measure = "auc") %>%
  .@y.values

# Very similar to LDA result, only slightly less favourable 

################### LOOCV #######################################

cv.error=rep(0,10)
# For polynomials 1,2,3... 10 fit a LOOCV model
for (i in 1:10){
  glm.fit=glm(Group~poly(X1+X2,i),data=class)
  cv.error[i]=cv.glm(class,glm.fit)$delta[1]
}
cv.error
# Create and display a plot
folds <- seq(1,10)
df <- data.frame(folds,cvError=cv.error)
ggplot(df,aes(x=folds,y=cvError)) + geom_point() +geom_line(color=twoClassColor[2],lwd=1.01) +
  xlab("Degree of Polynomial") + ylab("Cross Validation Error") +
  ggtitle("Leave one out Cross Validation") +
  theme_bw()

# May be worth including a second or third order term but no improvement with any higher terms

#################### Support vector machines ############################

# Support vector classifiers aim to maximise the distance between the two classes (optimal separating hyperplane) and works
# well with more than 2 predictors. This optimmisation is also conditional on a certain level of error which
# is allowed. Support vector machines use basis functions to to find a linear optimal hyperplane in the transformed space which 
# when transformed back will form a non-linear decision boundary. This non-linearity is controlled by different
# kernel functions. Primarily concerned with binary classification.

mdl <- train(x=xTrain,y=as.factor(yTrain), method = "svmLinear",
             trControl = trainControl("cv", number = 5),
             tuneGrid = expand.grid(C = seq(0, 2, length = 20)))
# Plot model accuracy vs different values of Cost
plot(mdl)
mdl$bestTune

mdl1 <- train(x=xTrain,y=as.factor(yTrain), method = "svmPoly")
              #trControl = trainControl("cv", number = 5))
             # tuneGrid = expand.grid(C = seq(0, 2, length = 20)))
mdl1$bestTune

svp <- ksvm(cbind(as.matrix(xTrain)[,2],as.matrix(xTrain)[,1]), yTrain, type = "C-svc", kernel = "vanilladot",C=mdl$bestTune,cross=5)
plot(data=cbind(as.matrix(xTrain)[,2],as.matrix(xTrain)[,1]),svp) # Support vectors in bold

# Decision boundary
w <- colSums(coef(svp)[[1]] * class[alphaindex(svp)[[1]],c('X1','X2')]) 
b <- b(svp)
class[, 1:2] <- sapply(class[,1:2], scale)
ggplot(data=class,aes(X1, X2, color=as.factor(Group))) +
  geom_point(size=2,alpha=0.8) + geom_abline(intercept=b/w[1], slope=-w[2]/w[1]) + 
  geom_abline(intercept=(b+1)/w[1], slope=-w[2]/w[1], linetype=2) + 
  geom_abline(intercept=(b-1)/w[1], slope=-w[2]/w[1], linetype=2) +
  theme_bw() + scale_colour_manual(name = 'Group', values = twoClassColor) +
  ggtitle("SVM Linear Kernel: Optimal Separating Hyperplane") + theme(legend.text = element_text(size = 10))

# Confusion matrix
yTestPred <- predict(svp, newdata=xTest)
confusionMatrix(as.factor(yTestPred), as.factor(yTest)) 

# Radial basis kernel plot
svp <- ksvm(cbind(as.matrix(xTrain)[,2],as.matrix(xTrain)[,1]), yTrain, type = "C-svc", kernel = "rbfdot")
plot(data=cbind(as.matrix(xTrain)[,2],as.matrix(xTrain)[,1]),svp,xlab="X1",ylab="X2")

svp <- ksvm(as.matrix(xTrain), yTrain, type = "C-svc", kernel = "rbfdot")

# ROC
prediction(predict(mdl1,newdata=xTest), yTest) %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot()
# AUC
prediction(predict(mdl1,newdata=xTest), yTest) %>%
  performance(measure = "auc") %>%
  .@y.values


################### K nearest neighbour ###################################

# This is a method of classifying points based on the group of its k nearest neighbours. A smaller value of 
# k will mean the method is more flexible and vice versa. Whilst the training error will consistently decline
# with a decrease in k, this will reach a point where there is an adverse effect on the test error rate due to
# the higher variance countering the reduction in bias.

class <- read.csv("Classification.csv", header = TRUE)
opts <- trainControl(method='repeatedcv', number=5, repeats=10, p=0.8)
# Optimal k (model)
set.seed(1040) 
knn <- train(x=xTrain, y=as.factor(yTrain), method='knn', 
             trControl=opts, tuneGrid=data.frame(k=seq(2, 20))) 
plot(knn)

yTestPred_knn <- as.factor(predict(knn, newdata = xTest))
confusionMatrix(yTestPred_knn, as.factor(yTest))

knn1 <- knn(xTrain,yTrain,cl=1,k=5)

# Decision boundary
library(ggpubr)
k19 <- ggplot(data = class, aes(x = X1, y = X2,
                         color = as.factor(Group))) +
  geom_contour(data = cbind(Grid, Group = predict(knn,Grid)),
               aes(z = as.numeric(Group)), color = "red", breaks = c(1.5),size=1) + geom_point(size=2,alpha=0.8) +
  ggtitle("K=15") + theme(legend.text = element_text(size = 10)) + 
  scale_colour_manual(name = 'Group', values = twoClassColor) + theme_bw()


get_best_result <- function(caret_fit) {
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  best_result
}

get_best_result(knn)

# K=5...
k5 <- ggplot(data = class, aes(x = X1, y = X2,
                         color = as.factor(Group))) +
  geom_contour(data = cbind(Grid, Group = knn(xTrain,Grid,cl=yTrain,k=5)),
               aes(z = as.numeric(Group)), color = "red", breaks = c(1.5),size=1) + geom_point(size=2,alpha=0.8) +
  ggtitle("K=5") + theme(legend.text = element_text(size = 10)) + 
  scale_colour_manual(name = 'Group', values = twoClassColor) + theme_bw()
figure1 <- ggarrange(k19,k5,ncol = 2, nrow = 1,common.legend = TRUE,legend = "right")
annotate_figure(figure1,top = text_grob("KNN Decision Boundary",face = "bold", size = 14))



# Another confusion matrix 
confusionMatrix(knn(xTrain,xTest,cl=yTrain,k=5),as.factor(yTest))

# Training and test errors

v <- vector("numeric", 20)  
for(i in 1:20){
  test_prediction <- knn(train=xTrain,test=xTest,cl=yTrain,k=i) 
  test_err <- mean(yTest != test_prediction)
  v[i] <- test_err
}

d <- vector("numeric", 20)  
for(i in 1:20){
  tr_prediction <- knn(train=xTrain,test=xTrain,cl=yTrain,k=i) 
  tr_err <- mean(yTrain != tr_prediction)
  d[i] <- tr_err
}

ggplot() + geom_point(aes(x=c(1:20),y=v)) + geom_line(aes(x=c(1:20),y=v,colour="Test Error"),lwd=1.01) + geom_point(aes(x=c(1:20),y=d)) + geom_line(aes(x=c(1:20),y=d,colour="Training Error"),lwd=1.01) + 
  theme_bw() + scale_colour_manual(name="",values=c(`Test Error`=twoClassColor[1],`Training Error`=twoClassColor[2])) + theme_bw() + 
  xlab("Number of Neighbours") + ylab("Error Rate") 
 

# Cross-validation
cMat <- NULL
predMat <- NULL
neighbors <- 2:20
set.seed(1040)
for(i in seq_along(neighbors)){
  fit = knn(xTrain,xTest,yTrain,k=i)
  table(fit,yTest)
  a <- confusionMatrix(as.factor(fit),as.factor(yTest))
  err <- mean(yTest != fit)
  cMat[i] <- a$overall[1]
  predMat[i] <- err
}
cMat
predMat

df <- data.frame(neighbors,Accuracy=cMat)
ggplot(df,aes(x=neighbors,y=Accuracy)) + geom_point() + geom_line(color=twoClassColor[2],lwd=1.01) +
  xlab("Number of Neighbours") + ylab("Accuracy") +
  ggtitle("KNN Regression - Accuracy vs Number of Neighbours") + theme_bw()


# Initially the optimal number of neighbours chosen using repeated cross validation was 19. When considering
# the test and training error rates however, the mininum error rate is shown at k=5 whilst the training error 
# increases relatively constantly. There appears to be little benefit to increasing the number of neighbours past
# 5. When assessing the confusion matrices, the option with fewer neighbours outperforms the higher number, again
# suggesting that k=5 appears more appropriate for the data. The decision boundaries do illustrate this difference,
# as the model with k=5 visually appears more flexible and captures the separation of the two classes best. Finally,
# the accuracy of each value of k was plotted against model accuract was shown to have the highest accuracy just above k=5. 
# These results suggest that the number of neighbours which will achieve the most accurate predictions will be 5. 

######################## True classifications ##########################################

class_true <- read.csv("ClassificationTrue.csv", header = TRUE)

# LDA
c1 <- confusionMatrix(predict(lda, newdata = class_true)$class,as.factor(class_true$Group))
yTestPred <- predict(lda, newdata = class_true)$class
e1 <- mean(class_true$Group != yTestPred)

# QDA
c2 <- confusionMatrix(predict(qda, newdata = class_true)$class,as.factor(class_true$Group))
yTestPred <- predict(qda, newdata = class_true)$class
e2 <- mean(class_true$Group != yTestPred)

# Logistic 
a <- predict(log.reg,newdata = class_true, type = "response")
b <- ifelse(a>0.5,1,0)
c3 <- confusionMatrix(as.factor(b),as.factor(class_true$Group))
e3 <- mean(class_true$Group != b)

# No improvement found when increasing the order of the polynomial 

# SVM
# Linear kernel
svp1 <- ksvm(as.matrix(xTrain), yTrain, type = "C-svc", kernel = "vanilladot",C=mdl$bestTune,cross=5)
yTestPred1 <- predict(svp1, newdata=cbind(class_true$X1,class_true$X2))
confusionMatrix(as.factor(yTestPred1), as.factor(class_true$Group)) 
# Radial
svp1 <- ksvm(as.matrix(xTrain), yTrain, type = "C-svc", kernel = "rbfdot")
yTestPred1 <- predict(svp1, newdata=cbind(class_true$X1,class_true$X2))
c4 <- confusionMatrix(as.factor(yTestPred1), as.factor(class_true$Group)) 
e4 <- mean(class_true$Group != yTestPred1)

# KNN
yTestPred2 <- knn(xTrain,cbind(class_true$X1,class_true$X2),cl=yTrain,k=5,prob=TRUE)
c5 <- confusionMatrix(knn(xTrain,cbind(class_true$X1,class_true$X2),cl=yTrain,k=5),as.factor(class_true$Group))
e5 <- mean(class_true$Group != yTestPred2)

comp <- data.frame(LDA=c1$overall[1],QDA=c2$overall[1],LR=c3$overall[1],SVM=c4$overall[1],KNN=c5$overall[1])
comp <- t(comp)
comp <- data.frame(Accuracy=comp[order(comp),])
acc <- data.frame(acc=t(data.frame(e1,e2,e3,e4,e5)))

ggplot(comp,aes(x=reorder(c("LDA","LR","KNN","SVM","QDA"), comp$Accuracy), y=comp$Accuracy)) +
  geom_segment(stat="identity",aes(xend=c("LDA","LR","KNN","SVM","QDA"), yend=0.5),lwd=1.01,color=twoClassColor[2]) +
  geom_point(size=4,color="orange") + ylim(0.5,1) + theme_bw() +
  xlab("Method") + ylab("Prediction Accuracy") + coord_flip() +
  ggtitle("Accuracy of each method in predicting the true values")

# Sensitivity

################### AUC ########################

# LDA
library(verification)
prediction(predict(lda, newdata = class_true)$posterior[,2], class_true$Group) %>%
  performance(measure = "auc") %>%
  .@y.values # 0.883
roc.area(class_true$Group,as.numeric(yTestPred))$A ## 0.679

# QDA
prediction(predict(qda, newdata = class_true)$posterior[,2], class_true$Group) %>%
  performance(measure = "auc") %>%
  .@y.values # 0.956
roc.area(class_true$Group,as.numeric(yTestPred))$A ## 0.838

# LR
prediction(a, class_true$Group) %>%
  performance(measure = "auc") %>%
  .@y.values # 0.883
roc.area(class_true$Group,as.numeric(b))$A ## 0.684


# SVM
roc.area(class_true$Group,as.numeric(yTestPred1))$A ## 0.836


# KNN
yTestPred2 <- as.numeric(yTestPred2) - 1
roc.area(class_true$Group,yTestPred2)$A ## 0.794
prediction(attr(yTestPred2,"prob"), class_true$Group) %>%
  performance(measure = "auc") %>%
  .@y.values 

prediction(attributes(yTestPred2)$prob, class_true$Group) %>%
  performance(measure = "auc") %>%
  .@y.values # 0.883


