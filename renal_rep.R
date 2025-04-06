library(readxl)
library("SpatialNP")
library(parallel)
library(MASS)
library("SPCAvRP")
library("mvtnorm")
library("mnormt")
library(ICSNP)
library(ICS)
library("MNM")
library(spatstat)
library(glasso)
library(flare)
library(glassoFast)
library(datamicroarray)
data <- read.csv("undergrad_design_project/Real Data Analysis/renal/Renal_GSE53757.csv")
data_renal <- data
contaminate_data <- function(X, r) {
  contaminated_X <- X 
  n <- nrow(X) 
  d <- ncol(X)  
  n_contaminate <- floor(n * r)  
  
 
  for (i in 1:d) {
    contaminated_entries <- sample(1:n, n_contaminate)  
    contaminated_X[contaminated_entries, i] <- sample(c(10, -10), n_contaminate, replace = TRUE)
  }
  
  return(contaminated_X)  
}


# Preprocessing
X <- data_renal[,-1:-2]
Y <- data_renal[,2]
X <- as.matrix(X)
Y <- as.numeric(Y == "ccRCC")



p_values <- numeric(ncol(X))

# Perform the t-test for each column (feature)
for (i in 1:ncol(X)) {
  # For each feature, perform a t-test comparing the two classes
  t_test_result <- t.test(X[which(Y == 0), i], X[which(Y == 1), i])
  p_values[i] <- t_test_result$p.value
}

# selected feartures
p <- 100
final_selected_features <- order(p_values)[1:p]
X <- X[, final_selected_features]
# normalization
X <- scale(X, center = TRUE, scale = TRUE)
# contamination
X <- contaminate_data(X,0.1)
# number of repeats
num_iterations <- 20  
methods <- c("CLIME", "SCLIME", "GLASSO", "SGLASSO")
results <- list()

# initialization
for (method in methods) {
  results[[method]] <- data.frame(
    Misclassification = numeric(num_iterations),
    Specificity = numeric(num_iterations),
    Sensitivity = numeric(num_iterations),
    MCC = numeric(num_iterations)
  )
}

for (iter in 1:num_iterations) {
  print(paste("Running iteration:", iter))
  
  # redivide training set and test set
  class_1_indices <- which(Y == 1)
  class_0_indices <- which(Y == 0)
  
  test_class_1_indices <- sample(class_1_indices, 50)
  test_class_0_indices <- sample(class_0_indices, 50)
  test_indices <- c(test_class_1_indices, test_class_0_indices)
  
  X_test <- X[test_indices, ]
  Y_test <- Y[test_indices]
  train_indices <- setdiff(1:nrow(X), test_indices)
  X_train <- X[train_indices, ]
  Y_train <- Y[train_indices]
  
  # compute hatSigma and hatS
  n0 <- sum(Y_train == 0)
  n1 <- sum(Y_train == 1)
  n <- n0 + n1
  hatSigma <- 1 / n * (n0 * cov(X_train[Y_train == 0, ]) + 
                         n1 * cov(X_train[Y_train == 1, ]))
  hatS <- 1 / n * (n0 * SCov(X_train[Y_train == 0, ]) + 
                     n1 * SCov(X_train[Y_train == 1, ]))
  
  LL <- 50
  BB <- 3
  
  ###### section: CLIME ######
  if (TRUE) {
    ## Cross-Validation
    cv1 <- rep(0, LL)
    cvb1 <- matrix(0, LL, BB)
    
    lambda_grid <-  exp(seq(log(0.01), log(1), length.out = LL)) 
    folds <- sample(1:BB, nrow(X_train), replace = TRUE)
    
    for (kb in 1:BB){
      X_cv_train = X_train[folds != kb,]
      Y_cv_train = Y_train[folds != kb]
      X_cv_test = X_train[folds == kb,]
      Y_cv_test = Y_train[folds == kb]
      n0_cv_train <- sum(Y_cv_train==0)
      n1_cv_train <- sum(Y_cv_train==1)
      n0_cv_test <- sum(Y_cv_test==0)
      n1_cv_test <- sum(Y_cv_test==1)
      # check
      if (n0_cv_test <= 1 | n1_cv_test <= 1 | n0_cv_train <= 1 | n1_cv_train <= 1) {
        next  # skip
      }
      hatmud=colMeans(X_cv_train[Y_cv_train==1,])-colMeans(X_cv_train[Y_cv_train==0,])
      hatmua=(colMeans(X_cv_train[Y_cv_train==1,])+colMeans(X_cv_train[Y_cv_train==0,]))/2
      cov_cv_train <- 1 / (n0_cv_train+n1_cv_train) * (n0_cv_train * cov(X_cv_train[Y_cv_train == 0, ]) + 
                                                         n1_cv_train * cov(X_cv_train[Y_cv_train == 1, ]))
      cov_cv_test <- 1 / (n0_cv_test+n1_cv_test) * (n0_cv_test * cov(X_cv_test[Y_cv_test == 0, ]) + 
                                                      n1_cv_test * cov(X_cv_test[Y_cv_test == 1, ]))
      out_cv =sugm(cov_cv_train, method = "clime", lambda = lambda_grid)
      for (ll in 1:LL)
      {
        ox1 <- out_cv$icov[[ll]]
        #result <- ((X_cv_test-matrix(rep(hatmua, nrow(X_cv_test)), ncol = length(hatmua), byrow = TRUE)) 
        #           %*% (ox1 %*% hatmud)+log(nrow(X_cv_train[Y_cv_train==1,])/nrow(X_cv_train[Y_cv_train==0,]))>0)
        #misclrt <- sum(result!= Y_cv_test)/length(Y_cv_test)
        misclrt <- sum(diag(ox1 %*% cov_cv_test)) - log(det(ox1))
        cvb1[ll,kb]<-misclrt
      }
    }
    cv1 <- rowMeans(cvb1)
    BK1 <- which.min(cv1)
    V1 <- sugm(hatSigma, method = "clime", lambda = lambda_grid[BK1])$icov[[1]]
    
    ## CLIME 
    ## Naive LDA for regular sample covariance matrix
    hatmud=colMeans(X_train[Y_train==1,])-colMeans(X_train[Y_train==0,])
    hatmua=(colMeans(X_train[Y_train==1,])+colMeans(X_train[Y_train==0,]))/2
    result1 <- ((X_test-matrix(rep(hatmua, nrow(X_test)), ncol = length(hatmua), byrow = TRUE)) 
                %*% (V1 %*% hatmud)+log(nrow(X_cv_train[Y_cv_train==1,])/nrow(X_train[Y_train==0,]))>0)
    TP1 <- sum((result1 == Y_test) & (Y_test > 0))
    TN1 <- sum((result1 == Y_test) & (Y_test == 0))
    FP1 <- sum((result1!= Y_test) & (Y_test == 0))
    FN1 <- sum((result1!= Y_test) & (Y_test > 0))
    misclrt1=sum(result1!= Y_test)/length(Y_test)
    Specificity1 = TN1/(TN1+FP1)
    Sensitivity1 = TP1/(TP1+FN1)
    MCC1 = (TP1*TN1-FP1*FN1)/(sqrt((TP1+FP1)*(TP1+FN1)*(TN1+FP1)*(TN1+FN1)))
    results[["CLIME"]][iter, ] <- c(misclrt1, Specificity1, Sensitivity1, MCC1)
  }
  
  ###### section: SCLIME ######
  if (TRUE) {
    ## Cross-Validation
    cv2 <- rep(0, LL)
    cvb2 <- matrix(0, LL, BB)
    
    lambda_grid <- exp(seq(log(0.01), log(1), length.out = LL)) 
    folds <- sample(1:BB, nrow(X_train), replace = TRUE)
    
    for (kb in 1:BB){
      X_cv_train = X_train[folds != kb,]
      Y_cv_train = Y_train[folds != kb]
      X_cv_test = X_train[folds == kb,]
      Y_cv_test = Y_train[folds == kb]
      n0_cv_train <- sum(Y_cv_train==0)
      n1_cv_train <- sum(Y_cv_train==1)
      n0_cv_test <- sum(Y_cv_test==0)
      n1_cv_test <- sum(Y_cv_test==1)
      if (n0_cv_test <= 1 | n1_cv_test <= 1 | n0_cv_train <= 1 | n1_cv_train <= 1) {
        next  # skip
      }
      hatmusd = spatial.median(X_cv_train[Y_cv_train == 1, ]) - spatial.median(X_cv_train[Y_cv_train == 0, ])
      hatmusa = (spatial.median(X_cv_train[Y_cv_train == 1, ]) + spatial.median(X_cv_train[Y_cv_train == 0, ])) / 2
      cov_cv_train <- p / (n0_cv_train+n1_cv_train) * (n0_cv_train * SCov(X_cv_train[Y_cv_train == 0, ]) + 
                                                         n1_cv_train * SCov(X_cv_train[Y_cv_train == 1, ]))
      cov_cv_test <- p / (n0_cv_test+n1_cv_test) * (n0_cv_test * SCov(X_cv_test[Y_cv_test == 0, ]) + 
                                                      n1_cv_test * SCov(X_cv_test[Y_cv_test == 1, ]))
      out_cv =sugm(cov_cv_train, method = "clime", lambda = lambda_grid)
      for (ll in 1:LL)
      {
        ox2 <- out_cv$icov[[ll]]
        #result <- ((X_cv_test-matrix(rep(hatmusa, nrow(X_cv_test)), ncol = length(hatmusa), byrow = TRUE)) 
        #           %*% (ox2 %*% hatmusd)+log(nrow(X_cv_train[Y_cv_train==1,])/nrow(X_cv_train[Y_cv_train==0,]))>0)
        #misclrt <- sum(result!= Y_cv_test)/length(Y_cv_test)
        misclrt <- sum(diag(ox2 %*% cov_cv_test)) - log(det(ox2))
        cvb2[ll,kb]<-misclrt
      }
      
    }
    
    cv2 <- rowMeans(cvb2)
    BK2 <- which.min(cv2)
    V2 <- sugm(p * hatS, method = "clime", lambda = lambda_grid[BK2])$icov[[1]]
    
    ## SCLIME
    ## Naive LDA for regular sample covariance matrix
    hatmusd=spatial.median(X_train[Y_train==1,])-spatial.median(X_train[Y_train==0,])
    hatmusa=(spatial.median(X_train[Y_train==1,])+spatial.median(X_train[Y_train==0,]))/2
    result2 <- ((X_test-matrix(rep(hatmusa, nrow(X_test)), ncol = length(hatmusa), byrow = TRUE)) 
                %*% (V2 %*% hatmusd)+log(nrow(X_cv_train[Y_cv_train==1,])/nrow(X_train[Y_train==0,]))>0)
    TP2 <- sum((result2 == Y_test) & (Y_test > 0))
    TN2 <- sum((result2 == Y_test) & (Y_test == 0))
    FP2 <- sum((result2!= Y_test) & (Y_test == 0))
    FN2 <- sum((result2!= Y_test) & (Y_test > 0))
    misclrt2=sum(result2!= Y_test)/length(Y_test)
    Specificity2 = TN2/(TN2+FP2)
    Sensitivity2 = TP2/(TP2+FN2)
    MCC2 = (TP2*TN2-FP2*FN2)/(sqrt((TP2+FP2)*(TP2+FN2)*(TN2+FP2)*(TN2+FN2)))
    results[["SCLIME"]][iter, ] <- c(misclrt2, Specificity2, Sensitivity2, MCC2)
  }
  
  ###### section: GLASSO ######
  if (TRUE) {
    ## Cross-Validation
    cv3 <- rep(0, LL)
    cvb3 <- matrix(0, LL, BB)
    
    lambda_grid <- exp(seq(log(0.01), log(2), length.out = LL)) 
    folds <- sample(1:BB, nrow(X_train), replace = TRUE)
    
    for (ll in 1:LL)
    {
      lambda <- lambda_grid[ll]
      
      for (kb in 1:BB){
        
        X_cv_train = X_train[folds != kb,]
        Y_cv_train = Y_train[folds != kb]
        X_cv_test = X_train[folds == kb,]
        Y_cv_test = Y_train[folds == kb]
        n0_cv_train <- sum(Y_cv_train==0)
        n1_cv_train <- sum(Y_cv_train==1)
        n0_cv_test <- sum(Y_cv_test==0)
        n1_cv_test <- sum(Y_cv_test==1)
        if (n0_cv_test <= 1 | n1_cv_test <= 1 | n0_cv_train <= 1 | n1_cv_train <= 1) {
          next  # skip
        }
        cov_cv_train <- 1 / (n0_cv_train+n1_cv_train) * (n0_cv_train * cov(X_cv_train[Y_cv_train == 0, ]) + 
                                                           n1_cv_train * cov(X_cv_train[Y_cv_train == 1, ]))
        cov_cv_test <- 1 / (n0_cv_test+n1_cv_test) * (n0_cv_test * cov(X_cv_test[Y_cv_test == 0, ]) + 
                                                        n1_cv_test * cov(X_cv_test[Y_cv_test == 1, ]))
        ox3 <- glassoFast(cov_cv_train,lambda)$wi
        hatmud=colMeans(X_cv_train[Y_cv_train==1,])-colMeans(X_cv_train[Y_cv_train==0,])
        hatmua=(colMeans(X_cv_train[Y_cv_train==1,])+colMeans(X_cv_train[Y_cv_train==0,]))/2
        #result <- ((X_cv_test-matrix(rep(hatmua, nrow(X_cv_test)), ncol = length(hatmua), byrow = TRUE)) 
        #           %*% (ox3 %*% hatmud)+ log(nrow(X_cv_train[Y_cv_train==1,])/nrow(X_cv_train[Y_cv_train==0,]))>0)
        #misclrt <- sum(result!= Y_cv_test)/length(Y_cv_test)
        misclrt <- sum(diag(ox3 %*% cov_cv_test)) - log(det(ox3))
        cvb3[ll,kb]<-misclrt
      }
    }
    
    cv3 <- rowMeans(cvb3)
    BK3 <- which.min(cv3)
    V3 <- glassoFast(hatSigma, lambda_grid[BK3])$wi
    
    ## GLASSO
    ## Naive LDA for regular sample covariance matrix
    hatmud=colMeans(X_train[Y_train==1,])-colMeans(X_train[Y_train==0,])
    hatmua=(colMeans(X_train[Y_train==1,])+colMeans(X_train[Y_train==0,]))/2
    result3 <- ((X_test-matrix(rep(hatmua, nrow(X_test)), ncol = length(hatmua), byrow = TRUE)) 
                %*% (V3 %*% hatmud)+log(nrow(X_cv_train[Y_cv_train==1,])/nrow(X_train[Y_train==0,]))>0)
    TP3 <- sum((result3 == Y_test) & (Y_test > 0))
    TN3 <- sum((result3 == Y_test) & (Y_test == 0))
    FP3 <- sum((result3!= Y_test) & (Y_test == 0))
    FN3 <- sum((result3!= Y_test) & (Y_test > 0))
    misclrt3=sum(result3!= Y_test)/length(Y_test)
    Specificity3 = TN3/(TN3+FP3)
    Sensitivity3 = TP3/(TP3+FN3)
    MCC3 = (TP3*TN3-FP3*FN3)/(sqrt((TP3+FP3)*(TP3+FN3)*(TN3+FP3)*(TN3+FN3)))
    results[["GLASSO"]][iter, ] <- c(misclrt3, Specificity3, Sensitivity3, MCC3)
  }
  
  ###### section: SGLASSO ######
  if (TRUE) {
    ## Cross-Validation
    cv4 <- rep(0, LL)
    cvb4 <- matrix(0, LL, BB)
    
    lambda_grid <- exp(seq(log(0.01), log(1), length.out = LL)) 
    folds <- sample(1:BB, nrow(X_train), replace = TRUE)
    
    for (ll in 1:LL)
    {
      lambda <- lambda_grid[ll]
      
      for (kb in 1:BB){
        
        X_cv_train = X_train[folds != kb,]
        Y_cv_train = Y_train[folds != kb]
        X_cv_test = X_train[folds == kb,]
        Y_cv_test = Y_train[folds == kb]
        n0_cv_train <- sum(Y_cv_train==0)
        n1_cv_train <- sum(Y_cv_train==1)
        n0_cv_test <- sum(Y_cv_test==0)
        n1_cv_test <- sum(Y_cv_test==1)
        if (n0_cv_test <= 1 | n1_cv_test <= 1 | n0_cv_train <= 1 | n1_cv_train <= 1) {
          next  # skip 
        }
        cov_cv_train <- p / (n0_cv_train+n1_cv_train) * (n0_cv_train * SCov(X_cv_train[Y_cv_train == 0, ]) + 
                                                           n1_cv_train * SCov(X_cv_train[Y_cv_train == 1, ]))
        cov_cv_test <- p / (n0_cv_test+n1_cv_test) * (n0_cv_test * SCov(X_cv_test[Y_cv_test == 0, ]) + 
                                                        n1_cv_test * SCov(X_cv_test[Y_cv_test == 1, ]))
        ox4 <- glassoFast(cov_cv_train,lambda)$wi
        hatmusd = spatial.median(X_cv_train[Y_cv_train == 1, ]) - spatial.median(X_cv_train[Y_cv_train == 0, ])
        hatmusa = (spatial.median(X_cv_train[Y_cv_train == 1, ]) + spatial.median(X_cv_train[Y_cv_train == 0, ])) / 2
        #result <- ((X_cv_test-matrix(rep(hatmusa, nrow(X_cv_test)), ncol = length(hatmusa), byrow = TRUE)) 
        #           %*% (ox4 %*% hatmusd)+ log(nrow(X_cv_train[Y_cv_train==1,])/nrow(X_cv_train[Y_cv_train==0,]))>0)
        #misclrt <- sum(result!= Y_cv_test)/length(Y_cv_test)
        misclrt <- sum(diag(ox4 %*% cov_cv_test)) - log(det(ox4))
        cvb4[ll,kb]<-misclrt
      }
    }
    
    cv4 <- rowMeans(cvb4)
    BK4 <- which.min(cv4)
    V4 <- glassoFast(p * hatS, lambda_grid[BK4])$wi
    
    ## SGLASSO
    ## Naive LDA for regular sample covariance matrix
    hatmusd=spatial.median(X_train[Y_train==1,])-spatial.median(X_train[Y_train==0,])
    hatmusa=(spatial.median(X_train[Y_train==1,])+spatial.median(X_train[Y_train==0,]))/2
    result4 <- ((X_test-matrix(rep(hatmusa, nrow(X_test)), ncol = length(hatmusa), byrow = TRUE)) 
                %*% (V4 %*% hatmusd)+log(nrow(X_cv_train[Y_cv_train==1,])/nrow(X_train[Y_train==0,]))>0)
    TP4 <- sum((result4 == Y_test) & (Y_test > 0))
    TN4 <- sum((result4 == Y_test) & (Y_test == 0))
    FP4 <- sum((result4!= Y_test) & (Y_test == 0))
    FN4 <- sum((result4!= Y_test) & (Y_test > 0))
    misclrt4=sum(result4!= Y_test)/length(Y_test)
    Specificity4 = TN4/(TN4+FP4)
    Sensitivity4 = TP4/(TP4+FN4)
    MCC4 = (TP4*TN4-FP4*FN4)/(sqrt((TP4+FP4)*(TP4+FN4)*(TN4+FP4)*(TN4+FN4)))
    results[["SGLASSO"]][iter, ] <- c(misclrt4, Specificity4, Sensitivity4, MCC4)
  }
}

# compute mean and standard error
for (method in methods) {
  print(paste("\n", method))
  print(colMeans(results[[method]]))
  print(apply(results[[method]], 2, sd))
}
