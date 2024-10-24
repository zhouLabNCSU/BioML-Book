library(randomForest)
library(xgboost)
library(e1071)
library(caret)
library(glmnet)
library(MASS)
library(neuralnet)

#Reference: http://dx.doi.org/10.1016/j.cell.2014.09.053
#load metadata and processed OTU table from MicrobiomeHD
metadata <- read.table("microbiome/MicrobiomeHD/ob_goodrich.metadata.txt",sep="\t",header=T,row.names=1)
otu_table <- read.table("microbiome/MicrobiomeHD/ob_goodrich.otu_table.100.denovo.rdp_assigned",sep="\t",header=T,row.names=1)
#vector of BMI as continuous trait
metadata <- subset(metadata, n_sample==0 & DiseaseState %in% c("H","OB"))
#select one sample per individual and one individual per twin pair
metadata <- metadata[!duplicated(metadata$familyid),]
y <- metadata$body_mass_index
#check that order of samples in metadata and OTU table are identical
otu_table <- otu_table[,rownames(metadata)]
#remove samples with fewer than 100 reads
otu_table <- otu_table[,which(colSums(otu_table)>=100)]
#remove OTUs with fewer than 10 reads
otu_table <- otu_table[which(rowSums(otu_table)>=10),]
#remove OTUs which were present in fewer than 5% of samples
otu_table <- otu_table[which(rowSums(otu_table>0)>=ncol(otu_table)*.05),]
rownames(otu_table) <- paste0("OTU",seq(1,nrow(otu_table)))
x <- data.matrix(otu_table)
#calculate relative abundance of each OTU by dividing its value by the total reads per sample
x <- sweep(x,2,colSums(x),"/")
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
    rf.train <- randomForest(x.train,y=y.train)
    rf.pred <- c(rf.pred, predict(rf.train,x.test))
    #gradient boosting
    xgboost.train <- xgboost(x.train,label=y.train,nrounds=10,objective="reg:linear",verbose=0)
    xgboost.pred <- c(xgboost.pred,predict(xgboost.train,x.test))
    #SVM
    svm.train <- svm(x.train, y.train)
    svm.pred <- c(svm.pred,predict(svm.train,x.test))
    #lasso
    lasso.train <- train(x.train,y.train,method="glmnet",tuneGrid=expand.grid(.alpha=1,.lambda=seq(0,1,by=0.1)))
    lasso.pred <- c(lasso.pred,predict(lasso.train,x.test))
    #ridge
    ridge.train <- train(x.train,y.train,method="glmnet",tuneGrid=expand.grid(.alpha=0,.lambda=seq(0,1,by=0.1)))
    ridge.pred <- c(ridge.pred,predict(ridge.train,x.test))
    #elastic net
    enet.train <- train(x.train,y.train,method="glmnet")
    enet.pred <- c(enet.pred,predict(enet.train,x.test))
    #k-nearest neighbors
    knn.train <- knnreg(x.train,y.train)
    knn.pred <- c(knn.pred,predict(knn.train,x.test))
    #neural networks
    train <- data.frame(y.train,x.train)
    names <- names(train)
    formula <- as.formula(paste("y.train ~", paste(names[!names %in% "y.train"], collapse = " + ")))
    neural.train <- neuralnet(formula,data=train,hidden=1,linear.output=T)
    neural.pred <- c(neural.pred,compute(neural.train,x.test)$net.result)
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
}
save(y.sample.matrix,y.test.matrix,rf.pred.matrix,xgboost.pred.matrix,svm.pred.matrix,lasso.pred.matrix,ridge.pred.matrix,enet.pred.matrix,knn.pred.matrix,neural.pred.matrix,file="microbiome/paper/MicrobiomeHD_bmi_goodrich_predict.Rdata")


load("microbiome/paper/MicrobiomeHD_bmi_goodrich_predict.Rdata")
#calculate average predicted y by sample, then generate ROC between predicted and actual y
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
}
identical(y,y.test.matrix[,1])
rf.pred.avg <- rowMeans(rf.pred.matrix)
rf <- cor(y,rf.pred.avg)
xgboost.pred.avg <- rowMeans(xgboost.pred.matrix)
xgboost <- cor(y,xgboost.pred.avg)
svm.pred.avg <- rowMeans(svm.pred.matrix)
svm <- cor(y,svm.pred.avg)
lasso.pred.avg <- rowMeans(lasso.pred.matrix)
lasso <- cor(y,lasso.pred.avg)
ridge.pred.avg <- rowMeans(ridge.pred.matrix)
ridge <- cor(y,ridge.pred.avg)
enet.pred.avg <- rowMeans(enet.pred.matrix)
enet <- cor(y,enet.pred.avg)
knn.pred.avg <- rowMeans(knn.pred.matrix)
knn <- cor(y,knn.pred.avg)
neural.pred.avg <- rowMeans(neural.pred.matrix)
neural <- cor(y,neural.pred.avg)
corr <- c(rf,xgboost,svm,lasso,ridge,enet,knn,neural)
corr

png(file="microbiome/paper/MicrobiomeHD_bmi_goodrich.png")
barplot(corr,names.arg=c("RF","Gboost","SVM","Lasso","Ridge","Enet","k-NN","Neural"),col=c("black","blue","red","green","orange","purple","brown","pink"),ylab="Average r",main="Goodrich (414 x 11,225)",ylim=c(0,0.4))
dev.off()
