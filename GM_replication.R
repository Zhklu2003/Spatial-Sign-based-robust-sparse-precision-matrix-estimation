# Load required libraries
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
library(igraph)
library(glasso)
library(glassoFast)

off_diagonal <- function(rho,p){
  A= matrix(rho, nrow=p, ncol =p)
  return(A-diag(diag(A)))
}
# Function to generate the graphical model
generate_graph <- function(d, max_degree) {
  Y <- matrix(runif(2 * d), nrow = 2, ncol = d)  # Step 1: Generate bivariate data points
  g <- graph.empty(d, directed = FALSE)          # Step 2: Initialize empty graph
  dist_matrix <- as.matrix(dist(t(Y)))           # Step 3: Calculate distances
  
  # Add edges to the graph
  for (i in 1:(d - 1)) {
    for (j in (i + 1):d) {
      dist_ij <- dist_matrix[i, j]
      prob <- exp(-dist_ij^2 / 0.25) / sqrt(2 * pi)
      if (runif(1) < prob && degree(g, i) < max_degree && degree(g, j) < max_degree) {
        g <- add_edges(g, c(i, j))
      }
    }
  }
  
  # Ensure degree constraint
  while (any(degree(g) > max_degree)) {
    for (i in 1:d) {
      if (degree(g, i) > max_degree) {
        edges_to_remove <- which(E(g)$from == i | E(g)$to == i)
        g <- delete_edges(g, edges_to_remove[1])  # Remove one edge
      }
    }
  }
  
  return(g)
}

# Function to generate the precision matrix Omega
generate_precision_matrix <- function(g, d) {
  Omega <- matrix(0, d, d)
  for (edge in E(g)) {
    vertices <- ends(g, edge)
    i <- vertices[1]
    j <- vertices[2]
    Omega[i, j] <- 0.145
    Omega[j, i] <- 0.145
  }
  diag(Omega) <- 1
  return(Omega)
}

# Function to generate contaminated data
generate_data <- function(n, d, Sigma, mi, r) {
  X <- rmvnorm(n, rep(0, d), Sigma)
  if (mi == 2) {
    X <- rmt(n, rep(0, d), Sigma, df = 3)/sqrt(3)
  } else if (mi == 3) {
    xr <- rbinom(n, 1, 0.8)
    X <- (X * xr + X * 3 * (1 - xr)) / sqrt(2.6)
  }
  
  # Add contamination
  n_contaminate <- floor(n * r)
  for (i in 1:d) {
    contaminated_entries <- sample(1:n, n_contaminate)
    X[contaminated_entries, i] <- sample(c(5, -5), n_contaminate, replace = TRUE)
  }
  
  return(X)
}

# Function to calculate ROC
calculate_roc <- function(results, Omega, gamma) {
  n_methods <- length(results)
  ROC <- matrix(0, nrow = 2, ncol = n_methods)
  
  for (i in 1:n_methods) {
    graph <- abs(results[[i]]) > gamma
    TP <- sum((graph > 0) & (Omega > 0))
    TN <- sum((graph == 0) & (Omega == 0))
    FP <- sum((graph > 0) & (Omega == 0))
    FN <- sum((graph == 0) & (Omega > 0))
    FPR <- FP / (FP + TN)
    TPR <- TP / (TP + FN)
    ROC[1, i] <- FPR
    ROC[2, i] <- TPR
  }
  
  return(ROC)
}

# Function to plot ROC curves
plot_roc_curves <- function(roc_list, legend_labels, colors) {
  plot(roc_list[[1]][1, ], roc_list[[1]][2, ], type = "b", col = colors[1],
       xlab = "False Positive Rate", ylab = "True Positive Rate", xlim = c(0, 1), ylim = c(0, 1),
       main = "ROC Curves for Different Methods")
  
  for (i in 2:length(roc_list)) {
    lines(roc_list[[i]][1, ], roc_list[[i]][2, ], type = "b", col = colors[i])
  }
  
  legend("bottomright", legend = legend_labels, col = colors, lty = 1, pch = 1)
}

# Main Experiment Function with Replication
run_experiment <- function(replications = 100, n = 400, d = 100, max_degree = 4, LL = 50, mi = 1, r = 0, gamma = 1e-5) {
  # Generate graph and precision matrix
  g <- generate_graph(d, max_degree)
  Omega <- generate_precision_matrix(g, d)
  Sigma <- solve(Omega)
  
  # Store the ROC results
  avg_ROC_clime <- matrix(0, 2, LL)
  avg_ROC_sclime <- matrix(0, 2, LL)
  avg_ROC_glasso <- matrix(0, 2, LL)
  avg_ROC_sglasso <- matrix(0, 2, LL)
  
  valid_replications <- 0  # Counter for valid replications
  
  for (rep in 1:replications) {
    # Generate data
    X <- generate_data(n, d, Sigma, mi, r)
    
    # Run CLIME, SCLIME, and Glasso
    lambda_grid <- exp(seq(log(0.005), log(1), length.out = LL))
    clime_results <- sugm(X, method = "clime", lambda = lambda_grid)$icov
    sclime_results <- sugm(p * SCov(X), method = "clime", lambda = lambda_grid)$icov
    glasso_results <- lapply(lambda_grid, function(lambda) glassoFast(cov(X), lambda)$wi)
    sglasso_results <- lapply(lambda_grid, function(lambda) glassoFast(p*SCov(X), lambda)$wi)
    
    # Check if any element in the graph is NaN or Inf for each method
    valid_clime <- all(sapply(clime_results, function(res) !any(is.nan(res)) && !any(is.infinite(res))))
    valid_sclime <- all(sapply(sclime_results, function(res) !any(is.nan(res)) && !any(is.infinite(res))))
    valid_glasso <- all(sapply(glasso_results, function(res) !any(is.nan(res)) && !any(is.infinite(res))))
    valid_sglasso <- all(sapply(sglasso_results, function(res) !any(is.nan(res)) && !any(is.infinite(res))))
    # If any result contains NaN or Inf, skip the replication
    if (valid_clime && valid_sclime && valid_glasso) {
      # Calculate ROC if the graph is valid
      ROC_clime <- calculate_roc(clime_results, Omega, gamma)
      ROC_sclime <- calculate_roc(sclime_results, Omega, gamma)
      ROC_glasso <- calculate_roc(glasso_results, Omega, gamma)
      ROC_sglasso <- calculate_roc(sglasso_results, Omega, gamma)
      # Accumulate ROC values
      avg_ROC_clime <- avg_ROC_clime + ROC_clime
      avg_ROC_sclime <- avg_ROC_sclime + ROC_sclime
      avg_ROC_glasso <- avg_ROC_glasso + ROC_glasso
      avg_ROC_sglasso <- avg_ROC_sglasso + ROC_sglasso
      
      # Increment valid replication count
      valid_replications <- valid_replications + 1
    }
  }
  
  # Check if the valid replications count matches the expected number
  if (valid_replications == 0) {
    stop("No valid replications found. Please check the data or model setup.")
  }
  
  # Calculate average ROC for each method
  avg_ROC_clime <- avg_ROC_clime / valid_replications
  avg_ROC_sclime <- avg_ROC_sclime / valid_replications
  avg_ROC_glasso <- avg_ROC_glasso / valid_replications
  avg_ROC_sglasso <- avg_ROC_sglasso / valid_replications
  
  return(list(avg_ROC_clime = avg_ROC_clime, avg_ROC_sclime = avg_ROC_sclime, avg_ROC_glasso = avg_ROC_glasso,avg_ROC_sglasso = avg_ROC_sglasso))
}

# Run Experiment and Analyze Results
result <- run_experiment(replications = 100, mi = 2)

# Plot ROC Curves
plot_roc_curves(list(result$avg_ROC_sclime, result$avg_ROC_sglasso, result$avg_ROC_clime, result$avg_ROC_glasso),
                legend_labels = c("SCLIME", "SGLASSO","CLIME", "GLASSO"),
                colors = c("red", "black","blue", "green"))    
