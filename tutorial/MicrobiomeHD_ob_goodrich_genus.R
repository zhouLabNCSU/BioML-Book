library(randomForest)
library(xgboost)
library(e1071)
library(caret)
library(glmnet)
library(MASS)
library(neuralnet)
library(pROC)
library(stringr)

#Reference: http://dx.doi.org/10.1016/j.cell.2014.09.053
#load metadata and processed OTU table from MicrobiomeHD
metadata <- read.table("microbiome/MicrobiomeHD/ob_goodrich.metadata.txt",sep="\t",header=T,row.names=1)
otu_table <- read.table("microbiome/MicrobiomeHD/ob_goodrich.otu_table.100.denovo.rdp_assigned",sep="\t",header=T,row.names=1)
#vector of discrete trait coded 1=OB cases and 0=Healthy controls
metadata <- subset(metadata, n_sample==0 & DiseaseState %in% c("H","OB"))
#select one sample per individual and one individual per twin pair
metadata <- metadata[!duplicated(metadata$familyid),]
y <- as.numeric(metadata$DiseaseState=="OB")
#check that order of samples in metadata and OTU table are identical
otu_table <- otu_table[,rownames(metadata)]
#remove samples with fewer than 100 reads
otu_table <- otu_table[,which(colSums(otu_table)>=100)]
#remove OTUs with fewer than 10 reads
otu_table <- otu_table[which(rowSums(otu_table)>=10),]
#remove OTUs which were present in fewer than 5% of samples
otu_table <- otu_table[which(rowSums(otu_table>0)>=ncol(otu_table)*.05),]
#calculate relative abundance of each OTU by dividing its value by the total reads per sample
otu_table <- sweep(otu_table,2,colSums(otu_table),"/")
#collapse OTUs to genus level by summing their respective relative abundances, discarding any OTUs which were un-annotated at the genus level
tax_table <- str_split_fixed(rownames(otu_table),";",n=8)
otu_table <- data.frame(genus=tax_table[,6],otu_table)
otu_table <- subset(otu_table, !(genus=="g__"))
otu_table <- aggregate(otu_table[,-1], by=list(otu_table$genus), sum)
rownames(otu_table) <- otu_table[,1]
otu_table <- otu_table[,-1]
x <- data.matrix(otu_table)
x <- t(x)
id <- seq(1:length(y))
numsets <- 100

#actual y representing 100 random samples of sample IDs
y.sample.matrix <- matrix(nrow=length(y),ncol=numsets)
y.test.matrix <- matrix(nrow=length(y),ncol=numsets)
#predictions are probabilitities that y=1 for each sample
rf.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
xgboost.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
svm.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
lasso.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
ridge.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
enet.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
knn.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
neural.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
lda.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
hclust.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
kmeans.pred.matrix <- matrix(nrow=length(y),ncol=numsets)
set.seed(1000)
for (n in 1:numsets) {
  print(n)
  #randomly split samples into five groups for 5-fold cross-validation
  subset <- matrix(0,nrow=83,ncol=5)
  subset[1:83,1] <- sort(sample(id,83))
  subset[1:83,2] <- sort(sample(which(!(id %in% subset)),83))
  subset[1:83,3] <- sort(sample(which(!(id %in% subset)),83))
  subset[1:83,4] <- sort(sample(which(!(id %in% subset)),83))
  subset[1:82,5] <- which(!(id %in% subset))
  y.test <- NULL
  rf.pred <- NULL
  xgboost.pred <- NULL
  svm.pred <- NULL
  lasso.pred <- NULL
  ridge.pred <- NULL
  enet.pred <- NULL
  knn.pred <- NULL
  neural.pred <- NULL
  lda.pred <- NULL
  #use 4 groups as training set and 1 group as testing set for each supervised method, and repeat until each group is tested
  for (i in 1:5) {
    print(i)
    y.test <- c(y.test, y[subset[,i]])
    y.train <- y[-subset[,i]]
    x.train <- x[-subset[,i],]
    #remove OTUs that have all zeroes in training set
    x.train <- subset(x.train,select=apply(x.train,2,var)>0)
    x.test <- x[subset[,i],colnames(x.train)]
    #random forest
    rf.train <- randomForest(x.train,y=as.factor(y.train))
    rf.pred <- c(rf.pred,predict(rf.train,x.test,type="prob")[,2])
    #gradient boosting
    xgboost.train <- xgboost(x.train,label=y.train,nrounds=10,objective="binary:logistic",verbose=0)
    xgboost.pred <- c(xgboost.pred,predict(xgboost.train,x.test))
    #SVM
    svm.train <- svm(x.train,as.factor(y.train),probability=T)
    svm.model <- predict(svm.train,x.test,probability=T)
    svm.pred <- c(svm.pred,attr(svm.model,"probabilities")[,1])
    #lasso
    lasso.train <- train(x.train,as.factor(y.train),method="glmnet",tuneGrid=expand.grid(.alpha=1,.lambda=seq(0,1,by=0.1)))
    lasso.pred <- c(lasso.pred,predict(lasso.train,x.test,type="prob")[,2])
    #ridge
    ridge.train <- train(x.train,as.factor(y.train),method="glmnet",tuneGrid=expand.grid(.alpha=0,.lambda=seq(0,1,by=0.1)))
    ridge.pred <- c(ridge.pred,predict(ridge.train,x.test,type="prob")[,2])
    #elastic net
    enet.train <- train(x.train,as.factor(y.train),method="glmnet")
    enet.pred <- c(enet.pred,predict(enet.train,x.test,type="prob")[,2])
    #k-nearest neighbors
    knn.train <- knn3(x.train,as.factor(y.train))
    knn.pred <- c(knn.pred,predict(knn.train,x.test,type="prob")[,2])
    #neural networks
    train <- data.frame(y.train,x.train)
    names <- names(train)
    formula <- as.formula(paste("y.train ~", paste(names[!names %in% "y.train"], collapse = " + ")))
    neural.train <- neuralnet(formula,data=train,hidden=1,linear.output=F)
    neural.pred <- c(neural.pred,compute(neural.train,x.test)$net.result)
    #LDA
    lda.train <- lda(x.train,y.train,tol=0)
    lda.pred <- c(lda.pred,predict(lda.train,x.test)$posterior[,2])
  }
  y.sample <- as.vector(subset)
  y.sample <- y.sample[which(y.sample>0)]
  y.sample.matrix[,n] <- y.sample
  y.test.matrix[,n] <- y.test
  rf.pred.matrix[,n] <- rf.pred
  xgboost.pred.matrix[,n] <- xgboost.pred
  svm.pred.matrix[,n] <- svm.pred
  lasso.pred.matrix[,n] <- lasso.pred
  ridge.pred.matrix[,n] <- ridge.pred
  enet.pred.matrix[,n] <- enet.pred
  knn.pred.matrix[,n] <- knn.pred
  neural.pred.matrix[,n] <- neural.pred
  lda.pred.matrix[,n] <- lda.pred
  x.sample <- x[y.sample,]
  #hierarchial clustering
  hclust.pred.matrix[,n] <- cutree(hclust(dist(x.sample)),k=2)
  #k-means clustering
  kmeans.pred.matrix[,n] <- kmeans(x.sample,centers=2)$cluster
}
save(y.sample.matrix,y.test.matrix,rf.pred.matrix,xgboost.pred.matrix,svm.pred.matrix,lasso.pred.matrix,ridge.pred.matrix,enet.pred.matrix,knn.pred.matrix,neural.pred.matrix,lda.pred.matrix,hclust.pred.matrix,kmeans.pred.matrix,file="microbiome/paper/MicrobiomeHD_ob_goodrich_genus.Rdata")


#calculate average predicted y by sample
load("microbiome/paper/MicrobiomeHD_ob_goodrich_genus.Rdata")
for (n in 1:numsets) {
  y.test.temp <- cbind(y.sample.matrix[,n],y.test.matrix[,n])
  y.test.temp <- y.test.temp[order(y.test.temp[,1]),]
  y.test.matrix[,n] <- y.test.temp[,2]
  rf.pred.temp <- cbind(y.sample.matrix[,n],rf.pred.matrix[,n])
  rf.pred.temp <- rf.pred.temp[order(rf.pred.temp[,1]),]
  rf.pred.matrix[,n] <- rf.pred.temp[,2]
  xgboost.pred.temp <- cbind(y.sample.matrix[,n],xgboost.pred.matrix[,n])
  xgboost.pred.temp <- xgboost.pred.temp[order(xgboost.pred.temp[,1]),]
  xgboost.pred.matrix[,n] <- xgboost.pred.temp[,2]
  svm.pred.temp <- cbind(y.sample.matrix[,n],svm.pred.matrix[,n])
  svm.pred.temp <- svm.pred.temp[order(svm.pred.temp[,1]),]
  svm.pred.matrix[,n] <- svm.pred.temp[,2]
  lasso.pred.temp <- cbind(y.sample.matrix[,n],lasso.pred.matrix[,n])
  lasso.pred.temp <- lasso.pred.temp[order(lasso.pred.temp[,1]),]
  lasso.pred.matrix[,n] <- lasso.pred.temp[,2]
  ridge.pred.temp <- cbind(y.sample.matrix[,n],ridge.pred.matrix[,n])
  ridge.pred.temp <- ridge.pred.temp[order(ridge.pred.temp[,1]),]
  ridge.pred.matrix[,n] <- ridge.pred.temp[,2]
  enet.pred.temp <- cbind(y.sample.matrix[,n],enet.pred.matrix[,n])
  enet.pred.temp <- enet.pred.temp[order(enet.pred.temp[,1]),]
  enet.pred.matrix[,n] <- enet.pred.temp[,2]
  knn.pred.temp <- cbind(y.sample.matrix[,n],knn.pred.matrix[,n])
  knn.pred.temp <- knn.pred.temp[order(knn.pred.temp[,1]),]
  knn.pred.matrix[,n] <- knn.pred.temp[,2]
  neural.pred.temp <- cbind(y.sample.matrix[,n],neural.pred.matrix[,n])
  neural.pred.temp <- neural.pred.temp[order(neural.pred.temp[,1]),]
  neural.pred.matrix[,n] <- neural.pred.temp[,2]
  lda.pred.temp <- cbind(y.sample.matrix[,n],lda.pred.matrix[,n])
  lda.pred.temp <- lda.pred.temp[order(lda.pred.temp[,1]),]
  lda.pred.matrix[,n] <- lda.pred.temp[,2]
  hclust.pred.temp <- cbind(y.sample.matrix[,n],hclust.pred.matrix[,n])
  hclust.pred.temp <- hclust.pred.temp[order(hclust.pred.temp[,1]),]
  hclust.pred.matrix[,n] <- hclust.pred.temp[,2]
  kmeans.pred.temp <- cbind(y.sample.matrix[,n],kmeans.pred.matrix[,n])
  kmeans.pred.temp <- kmeans.pred.temp[order(kmeans.pred.temp[,1]),]
  kmeans.pred.matrix[,n] <- kmeans.pred.temp[,2]
}
identical(y,y.test.matrix[,1])
rf.pred.avg <- rowMeans(rf.pred.matrix)
rf <- roc(y,rf.pred.avg)
rf$auc
xgboost.pred.avg <- rowMeans(xgboost.pred.matrix)
xgboost <- roc(y,xgboost.pred.avg)
xgboost$auc
svm.pred.avg <- rowMeans(svm.pred.matrix)
svm <- roc(y,svm.pred.avg)
svm$auc
lasso.pred.avg <- rowMeans(lasso.pred.matrix)
lasso <- roc(y,lasso.pred.avg)
lasso$auc
ridge.pred.avg <- rowMeans(ridge.pred.matrix)
ridge <- roc(y,ridge.pred.avg)
ridge$auc
enet.pred.avg <- rowMeans(enet.pred.matrix)
enet <- roc(y,enet.pred.avg)
enet$auc
knn.pred.avg <- rowMeans(knn.pred.matrix)
knn <- roc(y,knn.pred.avg)
knn$auc
neural.pred.avg <- rowMeans(neural.pred.matrix)
neural <- roc(y,neural.pred.avg)
neural$auc
lda.pred.avg <- rowMeans(lda.pred.matrix)
lda <- roc(y,lda.pred.avg)
lda$auc
hclust.pred.avg <- rowMeans(hclust.pred.matrix)
hclust <- roc(y,hclust.pred.avg)
hclust$auc
kmeans.pred.avg <- rowMeans(kmeans.pred.matrix)
kmeans <- roc(y,kmeans.pred.avg)
kmeans$auc

png(file="microbiome/paper/MicrobiomeHD_ob_goodrich_genus.png")
plot(rf$specificities,rf$sensitivities,type="n",xlim=c(1,0),xlab="Specificity",ylab="Sensitivity",main="Goodrich (135/279 x 99 OTUs collapsed at genus level)")
lines(rf$specificities,rf$sensitivities,col="black")
par(new=T)
plot(xgboost$specificities,xgboost$sensitivities,type="n",xlim=c(1,0),xlab=NA,ylab=NA,axes=F)
lines(xgboost$specificities,xgboost$sensitivities,col="blue")
par(new=T)
plot(svm$specificities,svm$sensitivities,type="n",xlim=c(1,0),xlab=NA,ylab=NA,axes=F)
lines(svm$specificities,svm$sensitivities,col="red")
par(new=T)
plot(lasso$specificities,lasso$sensitivities,type="n",xlim=c(1,0),xlab=NA,ylab=NA,axes=F)
lines(lasso$specificities,lasso$sensitivities,col="green")
par(new=T)
plot(ridge$specificities,ridge$sensitivities,type="n",xlim=c(1,0),xlab=NA,ylab=NA,axes=F)
lines(ridge$specificities,ridge$sensitivities,col="orange")
par(new=T)
plot(enet$specificities,enet$sensitivities,type="n",xlim=c(1,0),xlab=NA,ylab=NA,axes=F)
lines(enet$specificities,enet$sensitivities,col="purple")
par(new=T)
plot(knn$specificities,knn$sensitivities,type="n",xlim=c(1,0),xlab=NA,ylab=NA,axes=F)
lines(knn$specificities,knn$sensitivities,col="brown")
par(new=T)
plot(neural$specificities,neural$sensitivities,type="n",xlim=c(1,0),xlab=NA,ylab=NA,axes=F)
lines(neural$specificities,neural$sensitivities,col="pink")
par(new=T)
plot(lda$specificities,lda$sensitivities,type="n",xlim=c(1,0),xlab=NA,ylab=NA,axes=F)
lines(lda$specificities,lda$sensitivities,col="gold")
abline(1,-1,col="gray")
legend("bottomright",c("RF (0.65)","Gboost (0.66)","SVM (0.61)","Lasso (0.62)","Ridge (0.65)","Enet (0.66)","k-NN (0.53)","Neural (0.63)","LDA (0.63)"),col=c("black","blue","red","green","orange","purple","brown","pink","gold"),lty=1)
dev.off()
