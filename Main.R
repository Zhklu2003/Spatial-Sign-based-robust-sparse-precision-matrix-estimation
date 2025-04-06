library(flare)
library(parallel)
library(MASS)
library("SPCAvRP")
library("mvtnorm")
library("mnormt")
library(ICSNP)
library("ICS")
library("SpatialNP")
library("MNM")
library(glasso)


## Function to Generate Data
generate_data <- function(n, p, mi, Sigma) {
  X <- rmvnorm(n, rep(0, p), Sigma)
  if (mi == 2) X <- rmt(n, rep(0, p), Sigma, 3)/sqrt(3)
  if (mi == 3) {
    xr <- rbinom(n, 1, 0.8)
    X <- (X * xr + X * 3 * (1 - xr)) / sqrt(2.6)
  }
  return(X)
}

## Function to Run Estimation for a Single Replication
run_estimation <- function(n=100, d, mi, LL=50, BB=10, method = c("CLIME", "SCLIME", "Glasso", "SCAD")) {
  
  # Generate precision matrix Omega and covariance Sigma
  Omega <- matrix(0, d, d)
  for (i in 1:d) {
    for (j in 1:d) {
      Omega[i, j] <- 0.6 ^ abs(i - j)
    }
  }
  Omegas <- Omega / sum(diag(Omega)) * d
  Sigma <- solve(Omega)
  
  # Generate initial data
  X <- generate_data(n, d, mi, Sigma)
  hatSigma <- cov(X)
  
  # Initialize results
  results <- list()
  
  ## CLIME
  if ("CLIME" %in% method) {
    out1 <- sugm(X, method = "clime", nlambda = LL, lambda.min.ratio = 0.05)
    cv1 <- rep(0, LL)
    cvb1 <- matrix(0, LL, BB)
    for (ll in 1:LL) {
      for (bb in 1:BB) {
        X2 <- generate_data(n, d, mi, Sigma)
        ox1 <- out1$icov[[ll]]
        sx2 <- cov(X2)
        cvb1[ll, bb] <- sum(diag(ox1 %*% sx2)) - log(det(ox1))
      }
    }
    cv1 <- rowMeans(cvb1)
    BK1 <- which.min(cv1)
    eoc <- out1$icov[[BK1]] / sum(diag(out1$icov[[BK1]])) * d
    
    # Check for NaN or Inf in eoc
    if (any(is.nan(eoc)) || any(is.infinite(eoc))) {
      return(NULL)  # Drop this replication if NaN or Inf is detected
    }
    results$CLIME <- list(
      operator_norm = max(svd(eoc - Omegas)$d),
      L1_norm = max(colSums(abs(eoc - Omegas))),
      frobenius_norm = sqrt(sum((eoc - Omegas) ^ 2))
    )
  }
  
  ## SCLIME
  if ("SCLIME" %in% method) {
    SX <- SCov(X)
    out2 <- sugm(SX * p, method = "clime", nlambda = LL, lambda.min.ratio = 0.01)
    cv2 <- rep(0, LL)
    cvb2 <- matrix(0, LL, BB)
    for (ll in 1:LL) {
      for (bb in 1:BB) {
        X2 <- generate_data(n, d, mi, Sigma)
        ox1 <- out2$icov[[ll]]
        sx2 <- SCov(X2) * p
        cvb2[ll, bb] <- sum(diag(ox1 %*% sx2)) - log(det(ox1))
      }
    }
    cv2 <- rowMeans(cvb2)
    BK2 <- which.min(cv2)
    eos <- out2$icov[[BK2]] / sum(diag(out2$icov[[BK2]])) * d
    # Check for NaN or Inf in eoc
    if (any(is.nan(eos)) || any(is.infinite(eos))) {
      return(NULL)  # Drop this replication if NaN or Inf is detected
    }
    results$SCLIME <- list(
      operator_norm = max(svd(eos - Omegas)$d),
      L1_norm = max(colSums(abs(eos - Omegas))),
      frobenius_norm = sqrt(sum((eos - Omegas) ^ 2))
    )
  }
  
  ## Glasso
  if ("Glasso" %in% method) {
    lambda_grid <- exp(seq(log(0.01), log(1), length.out = LL))
    cv3 <- rep(0, LL)
    cvb3 <- matrix(0, LL, BB)
    for (ll in 1:LL) {
      for (bb in 1:BB) {
        X3 <- generate_data(n, d, mi, Sigma)
        sx3 <- cov(X3)
        ox3 <- glasso(hatSigma, lambda_grid[ll])$wi
        cvb3[ll, bb] <- sum(diag(ox3 %*% sx3)) - log(det(ox3))
      }
    }
    cv3 <- rowMeans(cvb3)
    BK3 <- which.min(cv3)
    eog <- glasso(hatSigma, lambda_grid[BK3])$wi / sum(diag(glasso(hatSigma, lambda_grid[BK3])$wi)) * d
    # Check for NaN or Inf in eoc
    if (any(is.nan(eog)) || any(is.infinite(eog))) {
      return(NULL)  # Drop this replication if NaN or Inf is detected
    }
    results$Glasso <- list(
      operator_norm = max(svd(eog - Omegas)$d),
      L1_norm = max(colSums(abs(eog - Omegas))),
      frobenius_norm = sqrt(sum((eog - Omegas) ^ 2))
    )
  }
  
  ## SCAD with Cross-Validation
  if ("SCAD" %in% method) {
    lambda_grid <- seq(0.1, 1, length.out = LL)
    cv4 <- rep(0, LL)
    cvb4 <- matrix(0, LL, BB)
    a <- 3.7
    tol <- 1e-5
    max_iter <- 100
    
    for (ll in 1:LL) {
      lambda <- lambda_grid[ll]
      for (bb in 1:BB) {
        X4 <- generate_data(n, d, mi, Sigma)
        ox4 <- matrix(0, d, d)
        for (iter in 1:max_iter) {
          weights <- matrix(0, nrow = d, ncol = d)
          for (i in 1:d) {
            for (j in 1:d) {
              omega <- abs(ox4[i, j])
              if (omega <= lambda) {
                weights[i, j] <- lambda
              } else if (omega <= a * lambda) {
                weights[i, j] <- (a * lambda - omega) / (a - 1)
              } else {
                weights[i, j] <- 0
              }
            }
          }
          glasso_result <- glasso(hatSigma, rho = weights)
          new_Omega <- glasso_result$wi
          if (max(abs(new_Omega - ox4)) < tol) {
            break
          }
          ox4 <- new_Omega
        }
        sx4 <- cov(X4)
        cvb4[ll, bb] <- sum(diag(ox4 %*% sx4)) - log(det(ox4))
      }
    }
    cv4 <- rowMeans(cvb4)
    BK4 <- which.min(cv4)
    lambda_opt <- lambda_grid[BK4]
    
    # Final SCAD estimation with best lambda
    eoad <- matrix(0, d, d)
    for (iter in 1:max_iter) {
      weights <- matrix(0, nrow = d, ncol = d)
      for (i in 1:d) {
        for (j in 1:d) {
          omega <- abs(eoad[i, j])
          if (omega <= lambda_opt) {
            weights[i, j] <- lambda_opt
          } else if (omega <= a * lambda_opt) {
            weights[i, j] <- (a * lambda_opt - omega) / (a - 1)
          } else {
            weights[i, j] <- 0
          }
        }
      }
      glasso_result <- glasso(hatSigma, rho = weights)
      new_Omega <- glasso_result$wi
      if (max(abs(new_Omega - eoad)) < tol) {
        break
      }
      eoad <- new_Omega
    }
    eoad <- eoad / sum(diag(eoad)) * d
    results$SCAD <- list(
      operator_norm = max(svd(eoad - Omegas)$d),
      L1_norm = max(colSums(abs(eoad - Omegas))),
      frobenius_norm = sqrt(sum((eoad - Omegas) ^ 2))
    )
  }
  
  return(results)
}

## Replicate the Experiment 100 Times
num_repeats <- 2
method <- c("CLIME", "SCLIME", "Glasso")  # Specify methods to run

# Initialize an empty list to store results
results <- list()

# Keep track of the number of valid replications
valid_repeats <- 0

# Run the replications
while (valid_repeats < num_repeats) {
  res <- run_estimation(d = 90, mi = 2, method = method)
  
  if (!is.null(res)) {
    results[[valid_repeats + 1]] <- res  # Store the valid result
    valid_repeats <- valid_repeats + 1    # Increment the valid replication counter
  }
}


# Extract norms for CLIME
clime_frobenius_norm <- sapply(results, function(res) res$CLIME$frobenius_norm)
clime_l1_norm <- sapply(results, function(res) res$CLIME$L1_norm)
clime_operator_norm <- sapply(results, function(res) res$CLIME$operator_norm)

# Extract norms for SCLIME
sclime_frobenius_norm <- sapply(results, function(res) res$SCLIME$frobenius_norm)
sclime_l1_norm <- sapply(results, function(res) res$SCLIME$L1_norm)
sclime_operator_norm <- sapply(results, function(res) res$SCLIME$operator_norm)

# Extract norms for Glasso
glasso_frobenius_norm <- sapply(results, function(res) res$Glasso$frobenius_norm)
glasso_l1_norm <- sapply(results, function(res) res$Glasso$L1_norm)
glasso_operator_norm <- sapply(results, function(res) res$Glasso$operator_norm)

# Helper function to calculate mean and standard deviation
calculate_mean_sd <- function(metric_values) {
  list(mean = mean(metric_values), sd = sd(metric_values))
}

# Calculate mean and SD for all norms
clime_metrics <- list(
  Frobenius = calculate_mean_sd(clime_frobenius_norm),
  L1 = calculate_mean_sd(clime_l1_norm),
  Operator = calculate_mean_sd(clime_operator_norm)
)

sclime_metrics <- list(
  Frobenius = calculate_mean_sd(sclime_frobenius_norm),
  L1 = calculate_mean_sd(sclime_l1_norm),
  Operator = calculate_mean_sd(sclime_operator_norm)
)

glasso_metrics <- list(
  Frobenius = calculate_mean_sd(glasso_frobenius_norm),
  L1 = calculate_mean_sd(glasso_l1_norm),
  Operator = calculate_mean_sd(glasso_operator_norm)
)

# Display results
cat("CLIME Metrics:\n")
cat("Frobenius Norm - Mean:", clime_metrics$Frobenius$mean, "SD:", clime_metrics$Frobenius$sd/sqrt(num_repeats), "\n")
cat("L1 Norm - Mean:", clime_metrics$L1$mean, "SD:", clime_metrics$L1$sd/sqrt(num_repeats), "\n")
cat("Operator Norm - Mean:", clime_metrics$Operator$mean, "SD:", clime_metrics$Operator$sd/sqrt(num_repeats), "\n\n")

cat("SCLIME Metrics:\n")
cat("Frobenius Norm - Mean:", sclime_metrics$Frobenius$mean, "SD:", sclime_metrics$Frobenius$sd/sqrt(num_repeats), "\n")
cat("L1 Norm - Mean:", sclime_metrics$L1$mean, "SD:", sclime_metrics$L1$sd/sqrt(num_repeats), "\n")
cat("Operator Norm - Mean:", sclime_metrics$Operator$mean, "SD:", sclime_metrics$Operator$sd/sqrt(num_repeats), "\n\n")

cat("Glasso Metrics:\n")
cat("Frobenius Norm - Mean:", glasso_metrics$Frobenius$mean, "SD:", glasso_metrics$Frobenius$sd/sqrt(num_repeats), "\n")
cat("L1 Norm - Mean:", glasso_metrics$L1$mean, "SD:", glasso_metrics$L1$sd/sqrt(num_repeats), "\n")
cat("Operator Norm - Mean:", glasso_metrics$Operator$mean, "SD:", glasso_metrics$Operator$sd/sqrt(num_repeats), "\n")


