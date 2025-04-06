## Load required packages
library(flare)
library(parallel)
library(MASS)
library("SPCAvRP")
library("mvtnorm")
library("mnormt")
library(ICSNP)
library(ICS)
library("SpatialNP")
library("MNM")
library(spatstat)
library(glasso)
library(glassoFast)

generate_Omega <- function(p) {
  
  # Step 1: Generate a symmetric matrix B
  B <- matrix(0, nrow = p, ncol = p)
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      B[i, j] <- ifelse(runif(1) < 0.1, 0.5, 0)
      B[j, i] <- B[i, j]  # Ensure symmetry
    }
  }

  # Step 2: Ensure positive definiteness by choosing delta
  eig_B <- eigen(B, symmetric = TRUE)
  lambda_max <- max(eig_B$values)
  lambda_min <- min(eig_B$values)
  
  # Solve for delta: (lambda_max + delta) / (lambda_min + delta) = p
  delta <- (p * lambda_min - lambda_max) / (1 - p)
  
  # Construct Omega
  Omega <- B + delta * diag(p)
  if(max(svd(Omega)$d)<0) Omega <- -Omega
  
  # Step 3: Standardize to have unit diagonals
  D_inv <- diag(1 / sqrt(diag(Omega)))
  Omega <- D_inv %*% Omega %*% D_inv

  
  return(Omega)
}

off_diagonal <- function(rho,p){
  A= matrix(rho, nrow=p, ncol =p)
  return(A-diag(diag(A)))
}
# Generate data
generate_data <- function(n0, n1, mu0, mu1, Sigma, mi) {
  X0 <- rmvnorm(n0, mean = mu0, sigma = Sigma)
  X1 <- rmvnorm(n1, mean = mu1, sigma = Sigma)
  
  if (mi == 2) {
    X0 <- rmvt(n0, sigma = Sigma, df = 3, delta = mu0)/sqrt(3)
    X1 <- rmvt(n1, sigma = Sigma, df = 3, delta = mu1)/sqrt(3)
  }
  if (mi == 3) {
    xr0 <- rbinom(n0, 1, 0.8)
    xr1 <- rbinom(n1, 1, 0.8)
    X0<-(X0*xr0+X0*3*(1-xr0))/sqrt(2.6)
    X1<-(X1*xr1+X1*3*(1-xr1))/sqrt(2.6)
  }
  
  X <- rbind(X0, X1)
  Y <- c(rep(0, n0), rep(1, n1))
  return(list(X = X, Y = Y))
}

# Cross-validation function
cross_validate <- function(method, LL, BB, lambda_grid, train_data, n0, n1, mu0, mu1, Sigma, mi) {
  p <- ncol(train_data$X)
  cvb <- matrix(0, LL, BB)
  
  # Model fitting and covariance usage
  if (method == "CLIME") {
    out1 <- sugm(train_data$hatSigma, method = "clime", ,nlambda=LL,lambda.min.ratio=0.05)
    #out2 <- sugm(train_data$hatS * p, method = "clime", lambda = lambda_grid)
    for (ll in 1:LL) {
      lambda <- lambda_grid[ll]
      for (bb in 1:BB) {
        # Generate independent CV data
        val_data <- generate_data(n0, n1, mu0, mu1, Sigma, mi)
        X_cv <- val_data$X
        Y_cv <- val_data$Y
        
        ox <- out1$icov[[ll]]
        hatmud <- colMeans(train_data$X[train_data$Y == 1, ]) - colMeans(train_data$X[train_data$Y == 0, ])
        hatmua <- (colMeans(train_data$X[train_data$Y == 1, ]) + colMeans(train_data$X[train_data$Y == 0, ])) / 2
        predictions <- ((X_cv - matrix(rep(hatmua, nrow(X_cv)), ncol = p, byrow = TRUE)) %*% 
                          (ox %*% hatmud)+log(n1/n0)) > 0
        cvb[ll, bb] <- sum(predictions != Y_cv) / length(Y_cv)
        
      }
    }   
    cv <- rowMeans(cvb)
    BK <-  which.min(cv)
    V <- out1$icov[[BK]]
    return(V)
  } else if (method == "SCLIME") {
    #out1 <- sugm(train_data$hatSigma, method = "clime", lambda = lambda_grid)
    out2 <- sugm(train_data$hatS * p, method = "clime", ,nlambda=LL,lambda.min.ratio=0.01)
    for (ll in 1:LL) {
      lambda <- lambda_grid[ll]
      for (bb in 1:BB) {
        # Generate independent CV data
        val_data <- generate_data(n0, n1, mu0, mu1, Sigma, mi)
        X_cv <- val_data$X
        Y_cv <- val_data$Y
        
        ox <- out2$icov[[ll]]
        
        hatmusd <- spatial.median(train_data$X[train_data$Y == 1, ]) - spatial.median(train_data$X[train_data$Y == 0, ])
        hatmusa <- (spatial.median(train_data$X[train_data$Y == 1, ]) + spatial.median(train_data$X[train_data$Y == 0, ])) / 2
        predictions <- ((X_cv - matrix(rep(hatmusa, nrow(X_cv)), ncol = p, byrow = TRUE)) %*% 
                          (ox %*% hatmusd)+log(n1/n0)) > 0
        cvb[ll, bb] <- sum(predictions != Y_cv) / length(Y_cv)
        
      }
    }
    cv <- rowMeans(cvb)
    BK <-  which.min(cv)
    V <- out2$icov[[BK]]
    return(V)    
  } else if (method == "Glasso") {
    for (ll in 1:LL) {
      lambda <- lambda_grid[ll]
      for (bb in 1:BB) {
        # Generate independent CV data
        val_data <- generate_data(n0, n1, mu0, mu1, Sigma, mi)
        X_cv <- val_data$X
        Y_cv <- val_data$Y
        
        
        ox <- glassoFast(train_data$hatSigma, rho = lambda)$wi
        hatmud <- colMeans(train_data$X[train_data$Y == 1, ]) - colMeans(train_data$X[train_data$Y == 0, ])
        hatmua <- (colMeans(train_data$X[train_data$Y == 1, ]) + colMeans(train_data$X[train_data$Y == 0, ])) / 2
        predictions <- ((X_cv - matrix(rep(hatmua, nrow(X_cv)), ncol = p, byrow = TRUE)) %*% 
                          (ox %*% hatmud)+log(n1/n0)) > 0
        cvb[ll, bb] <- sum(predictions != Y_cv) / length(Y_cv)
        
      }
    }   
    cv <- rowMeans(cvb)
    BK <-  which.min(cv)
    V <- glassoFast(train_data$hatSigma,lambda_grid[BK])$wi
    return(V)
    }else if (method == "SGlasso") {
      for (ll in 1:LL) {
        lambda <- lambda_grid[ll]
        for (bb in 1:BB) {
          # Generate independent CV data
          val_data <- generate_data(n0, n1, mu0, mu1, Sigma, mi)
          X_cv <- val_data$X
          Y_cv <- val_data$Y
          
          
          ox <- glassoFast(train_data$hatS * p, rho = lambda)$wi
          hatmusd <- spatial.median(train_data$X[train_data$Y == 1, ]) - spatial.median(train_data$X[train_data$Y == 0, ])
          hatmusa <- (spatial.median(train_data$X[train_data$Y == 1, ]) + spatial.median(train_data$X[train_data$Y == 0, ])) / 2
          predictions <- ((X_cv - matrix(rep(hatmusa, nrow(X_cv)), ncol = p, byrow = TRUE)) %*% 
                            (ox %*% hatmusd)+log(n1/n0)) > 0
          cvb[ll, bb] <- sum(predictions != Y_cv) / length(Y_cv)
          
        }
      }   
      cv <- rowMeans(cvb)
      BK <-  which.min(cv)
      V <- glassoFast(train_data$hatS * p,lambda_grid[BK])$wi
      return(V)
  }
  
  # Compute misclassification rate on CV data
  
  
  
}

# Calculate metrics
calculate_metrics <- function(predictions, Y_test) {
  misclrt <- sum(predictions != Y_test) / length(Y_test)
  TP <- sum((predictions == Y_test) & (Y_test > 0))
  TN <- sum((predictions == Y_test) & (Y_test == 0))
  FP <- sum((predictions != Y_test) & (Y_test == 0))
  FN <- sum((predictions != Y_test) & (Y_test > 0))
  Specificity <- TN / (TN + FP)
  Sensitivity <- TP / (TP + FN)
  MCC <- (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  return(list(misclrt = misclrt, Specificity = Specificity, Sensitivity = Sensitivity, MCC = MCC))
}

# Main experiment function
run_experiment <- function(n = 200, p = 60, s = 10, a = 0.05, mi = 1, LL = 50, BB = 10, pi = 0.5, lambda_grid = NULL) {
  # Generate data and prepare training and test sets
  mu0 <- rep(0, p)
  mu1 <- c(rep(a, s), rep(0, p - s)) * 10
  Omega <- matrix(0, p, p)
  for (i in 1:p) {
   for (j in 1:p) {
      Omega[i, j] <- 0.6 ^ abs(i - j)
    }
  }
  Omega <- generate_Omega(p)
  Sigma <- solve(Omega)
  #Sigma <-Omega
 
  
  n0 <- rbinom(1, n, pi)
  n1 <- n - n0
  train_data <- generate_data(n0, n1, mu0, mu1, Sigma, mi)
  
  
  # Precompute hatSigma and hatS for training data
  train_data$hatSigma <- 1 / n * (n0 * cov(train_data$X[train_data$Y == 0, ]) + 
                                    n1 * cov(train_data$X[train_data$Y == 1, ]))
  train_data$hatS <- 1 / n * (n0 * SCov(train_data$X[train_data$Y == 0, ]) + 
                                n1 * SCov(train_data$X[train_data$Y == 1, ]))
  
  # Default lambda grid if not provided
  if (is.null(lambda_grid)) {
    lambda_grid <- exp(seq(log(0.005), log(1), length.out = LL))
  }
  
  # Perform Cross-Validation for each method
  clime_cv <- cross_validate("CLIME", LL, BB, lambda_grid, train_data, n0, n1, mu0, mu1, Sigma, mi)
  sclime_cv <- cross_validate("SCLIME", LL, BB, lambda_grid, train_data, n0, n1, mu0, mu1, Sigma, mi)
  glasso_cv <- cross_validate("Glasso", LL, BB, lambda_grid, train_data, n0, n1, mu0, mu1, Sigma, mi)
  sglasso_cv <- cross_validate("SGlasso", LL, BB, lambda_grid, train_data, n0, n1, mu0, mu1, Sigma, mi)
  # Check for NaN or Inf in CV results, skip if found
  if (any(is.nan(clime_cv)) || any(is.infinite(clime_cv)) ||
      any(is.nan(sclime_cv)) || any(is.infinite(sclime_cv)) ||
      any(is.nan(glasso_cv)) || any(is.infinite(glasso_cv))||
      any(is.nan(sglasso_cv)) || any(is.infinite(sglasso_cv))) {
    return(NULL)  # Return NULL if any NaN or Inf found
  }
  
  # Compute results for this replication
  n0_test <- rbinom(1, n, pi)
  n1_test <- n - n0_test
  test_data <- generate_data(n0_test, n1_test, mu0, mu1, Sigma, mi)
  X_test <- test_data$X
  Y_test <- test_data$Y
  
  results <- list()
  
  # Metrics calculation for each method
  # CLIME
  V_clime <- clime_cv
  hatmud <- colMeans(train_data$X[train_data$Y == 1, ]) - colMeans(train_data$X[train_data$Y == 0, ])
  hatmua <- (colMeans(train_data$X[train_data$Y == 1, ]) + colMeans(train_data$X[train_data$Y == 0, ])) / 2
  predictions <- ((X_test - matrix(rep(hatmua, nrow(X_test)), ncol = p, byrow = TRUE)) %*% 
                    (V_clime %*% hatmud)+log(n1/n0)) > 0
  results$CLIME <- calculate_metrics(predictions, Y_test)
  
  # SCLIME
  V_sclime <- sclime_cv
  hatmusd <- spatial.median(train_data$X[train_data$Y == 1, ]) - spatial.median(train_data$X[train_data$Y == 0, ])
  hatmusa <- (spatial.median(train_data$X[train_data$Y == 1, ]) + spatial.median(train_data$X[train_data$Y == 0, ])) / 2
  predictions <- ((X_test - matrix(rep(hatmusa, nrow(X_test)), ncol = p, byrow = TRUE)) %*% 
                    (V_sclime %*% hatmusd)+log(n1/n0)) > 0
  results$SCLIME <- calculate_metrics(predictions, Y_test)
  
  # Glasso
  V_glasso <- glasso_cv
  predictions <- ((X_test - matrix(rep(hatmua, nrow(X_test)), ncol = p, byrow = TRUE)) %*% 
                    (V_glasso %*% hatmud)+log(n1/n0)) > 0
  results$Glasso <- calculate_metrics(predictions, Y_test)
  
  # SGlasso
  V_sglasso <- sglasso_cv
  predictions <- ((X_test - matrix(rep(hatmusa, nrow(X_test)), ncol = p, byrow = TRUE)) %*% 
                    (V_sglasso %*% hatmusd)+log(n1/n0)) > 0
  results$SGlasso <- calculate_metrics(predictions, Y_test)
  
  return(results)
}

# Repeat function to run experiments until valid replicates are reached
repeat_experiment <- function(n_repeats, n = 200, p = 60, s = 10, a = 0.05, mi = 1, LL = 50, BB = 10, pi = 0.5, lambda_grid = NULL) {
  valid_repeats <- 0
  results_list <- list()
  
  while (valid_repeats < n_repeats) {
    result <- run_experiment(n, p, s, a, mi, LL, BB, pi, lambda_grid)
    
    if (!is.null(result)) {
      valid_repeats <- valid_repeats + 1
      results_list[[valid_repeats]] <- result
    }
  }
  
  return(results_list)
}

# Example usage:
num_repeats <- 100
experiment_results <- repeat_experiment(n_repeats = num_repeats, mi = 2, p = 30)


clime_misclrt <- sapply(experiment_results, function(res) res$CLIME$misclrt)
clime_specificity <- sapply(experiment_results, function(res) res$CLIME$Specificity)
clime_sensitivity <- sapply(experiment_results, function(res) res$CLIME$Sensitivity)
clime_mcc <- sapply(experiment_results, function(res) res$CLIME$MCC)

# Extract specific metrics for SCLIME
sclime_misclrt <- sapply(experiment_results, function(res) res$SCLIME$misclrt)
sclime_specificity <- sapply(experiment_results, function(res) res$SCLIME$Specificity)
sclime_sensitivity <- sapply(experiment_results, function(res) res$SCLIME$Sensitivity)
sclime_mcc <- sapply(experiment_results, function(res) res$SCLIME$MCC)

# Extract specific metrics for Glasso
glasso_misclrt <- sapply(experiment_results, function(res) res$Glasso$misclrt)
glasso_specificity <- sapply(experiment_results, function(res) res$Glasso$Specificity)
glasso_sensitivity <- sapply(experiment_results, function(res) res$Glasso$Sensitivity)
glasso_mcc <- sapply(experiment_results, function(res) res$Glasso$MCC)

# Extract specific metrics for SGlasso
sglasso_misclrt <- sapply(experiment_results, function(res) res$SGlasso$misclrt)
sglasso_specificity <- sapply(experiment_results, function(res) res$SGlasso$Specificity)
sglasso_sensitivity <- sapply(experiment_results, function(res) res$SGlasso$Sensitivity)
sglasso_mcc <- sapply(experiment_results, function(res) res$SGlasso$MCC)

# Calculate mean and standard deviation for each metric
cat("CLIME - Mean Misclassification Rate:", mean(clime_misclrt), "SD:", sd(clime_misclrt)/sqrt(num_repeats), "\n")
cat("CLIME - Mean Specificity:", mean(clime_specificity), "SD:", sd(clime_specificity)/sqrt(num_repeats), "\n")
cat("CLIME - Mean Sensitivity:", mean(clime_sensitivity), "SD:", sd(clime_sensitivity)/sqrt(num_repeats), "\n")
cat("CLIME - Mean MCC:", mean(clime_mcc), "SD:", sd(clime_mcc)/sqrt(num_repeats), "\n\n")

cat("SCLIME - Mean Misclassification Rate:", mean(sclime_misclrt), "SD:", sd(sclime_misclrt)/sqrt(num_repeats), "\n")
cat("SCLIME - Mean Specificity:", mean(sclime_specificity), "SD:", sd(sclime_specificity)/sqrt(num_repeats), "\n")
cat("SCLIME - Mean Sensitivity:", mean(sclime_sensitivity), "SD:", sd(sclime_sensitivity)/sqrt(num_repeats), "\n")
cat("SCLIME - Mean MCC:", mean(sclime_mcc), "SD:", sd(sclime_mcc)/sqrt(num_repeats), "\n\n")

cat("Glasso - Mean Misclassification Rate:", mean(glasso_misclrt), "SD:", sd(glasso_misclrt)/sqrt(num_repeats), "\n")
cat("Glasso - Mean Specificity:", mean(glasso_specificity), "SD:", sd(glasso_specificity)/sqrt(num_repeats), "\n")
cat("Glasso - Mean Sensitivity:", mean(glasso_sensitivity), "SD:", sd(glasso_sensitivity)/sqrt(num_repeats), "\n")
cat("Glasso - Mean MCC:", mean(glasso_mcc), "SD:", sd(glasso_mcc)/sqrt(num_repeats), "\n\n")

cat("SGlasso - Mean Misclassification Rate:", mean(sglasso_misclrt), "SD:", sd(sglasso_misclrt)/sqrt(num_repeats), "\n")
cat("SGlasso - Mean Specificity:", mean(sglasso_specificity), "SD:", sd(sglasso_specificity)/sqrt(num_repeats), "\n")
cat("SGlasso - Mean Sensitivity:", mean(sglasso_sensitivity), "SD:", sd(sglasso_sensitivity)/sqrt(num_repeats), "\n")
cat("SGlasso - Mean MCC:", mean(sglasso_mcc), "SD:", sd(sglasso_mcc)/sqrt(num_repeats), "\n")