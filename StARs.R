lambda_call<- function(X, method ="glasso", nlambda = NULL, lambda.min.ratio = NULL)
{
  n <- nrow(X)
  p <- ncol(X)
  if(is.null(nlambda))
    nlambda = 10
  if(is.null(lambda.min.ratio))
    lambda.min.ratio = 0.1
  if(method =="sclime"|method =="clime"){
    lambda.max = 1
  }
  if(method =="sglasso"|method =="glasso"){
    X =scale(X)
    S = cor(X)
  lambda.max = max(max(S-diag(p)),-min(S-diag(p)))
  }
  lambda.min = lambda.min.ratio*lambda.max
  lambda_grid = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
}

cor_of_cov <- function(X){
  D <- diag(1 / sqrt(diag(X)))
  standardized_X <- D %*% X %*% D
  return(standardized_X)
}


precision_matrix_StARs <- 
  function(data, method, stars.thresh = 0.1,lambda_grid = NULL, nlambda = 10, stars.subsample.ratio = NULL, rep.num =20, verbose = TRUE){
  n <- nrow(data)
  p <- ncol(data)
  if(is.null(stars.subsample.ratio))
  {
    if(n>144) stars.subsample.ratio = 10*sqrt(n)/n
    if(n<=144) stars.subsample.ratio = 0.8
  }

  merge = list()
  if(method == "clime"|method =="sclime"){
     lambda_list = lambda_call(data, method ="clime", nlambda=nlambda)
  }
  if(method == "glasso"| method =="sglasso"){
    lambda_list = lambda_call(data, method ="glasso", nlambda=nlambda)
  }
 
  for(i in 1:nlambda) merge[[i]] = Matrix(0,p,p)
  
  for(j in 1:rep.num)
  {
    if(verbose)
    {
      mes <- paste(c("Conducting Subsampling....in progress:", floor(100*j/rep.num), "%"), collapse="")
      cat(mes, "\r")
      flush.console()
    }
    ind.sample = sample(c(1:n), floor(n*stars.subsample.ratio), replace=FALSE)
    if(method == "clime"){
      cov = cor(data[ind.sample,])
      tmp = sugm(cov, method = "clime",lambda = lambda_list)$icov
    }
    if(method == "sclime"){
      Scov = p*SCov(data[ind.sample,])
      Scov = cor_of_cov(Scov)
      tmp = sugm(Scov, method = "clime",lambda = lambda_list)$icov
    }
    if(method == "glasso"){
      cov = cor(data[ind.sample,])
      tmp = lapply(lambda_list, function(lambda) glassoFast(cov, rho = lambda)$wi)
    }
    if(method == "sglasso"){
      Scov = p*SCov(data[ind.sample,])
      Scov = cor_of_cov(Scov)
      tmp = lapply(lambda_list, function(lambda) glassoFast(Scov, rho = lambda)$wi)
    }
    for(i in 1:nlambda){
      merge[[i]] = merge[[i]] + (abs(tmp[[i]]) !=0)
    }
    rm(ind.sample,tmp)
    gc()
  }
  
  if(verbose)
  {
    mes = "Conducting Subsampling....done.                 "
    cat(mes, "\r")
    cat("\n")
    flush.console()
  }
 
  variability = rep(0,nlambda)
  for(i in 1:nlambda){
    merge[[i]] = merge[[i]]/rep.num
    variability[i] = 4*sum(merge[[i]]*(1-merge[[i]]))/(p*(p-1))
  }
  
  opt.index = max(which.max(variability >= stars.thresh)[1]-1,1)
  opt.lambda = lambda_list[opt.index]
  if (method =="clime"){
    cov = cor(data)
    best_precision_matrix = sugm(cov, method = "clime",lambda = opt.lambda)$icov[[1]]
  }
  if(method == "sclime"){
    Scov = p*SCov(data)
    Scov <- cor_of_cov(Scov)
    best_precision_matrix = sugm(Scov, method = "clime",lambda = opt.lambda)$icov[[1]]
  }
  if (method =="glasso"){
    cov1 = cor(data)
    best_precision_matrix = glassoFast(cov1, opt.lambda )$wi
  }
  if (method =="sglasso"){
    Scov = p*SCov(data)
    Scov <- cor_of_cov(Scov)
    best_precision_matrix = glassoFast(Scov, opt.lambda )$wi
  }

  return(list(lambda=opt.lambda, index = opt.index, pm = best_precision_matrix))
}

contaminate_data <- function(X, r) {
  contaminated_X <- X  
  n <- nrow(X)  # sample
  d <- ncol(X)  # variable
  n_contaminate <- floor(n * r)  # compute the number of contaminated points
  
  # column-wise contamination
  for (i in 1:d) {
    contaminated_entries <- sample(1:n, n_contaminate)  # random contamination
    contaminated_X[contaminated_entries, i] <- sample(c(0.13, -0.13), n_contaminate, replace = TRUE)
  }
  
  return(contaminated_X)  # 
}
#####################################################################
# Load required libraries
library(huge)
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
library(ggraph)
library(withr)

# read data
log_returns <- read.csv("undergrad_design_project/Real Data Analysis/GM_stock/sp500_log_returns_2015_2019.csv", row.names = 1)
gics_data <- read.csv("undergrad_design_project/Real Data Analysis/GM_stock/sp500_gics_sectors.csv")

# Ensure log_returns is a matrix instead of a DataFrame
log_returns_matrix <- as.matrix(log_returns)
log_returns_matrix <- with_seed(124, contaminate_data(log_returns_matrix, 0.005))

StARs <- precision_matrix_StARs(log_returns_matrix, method = "sglasso")
best_precision_matrix = StARs$pm  
opt.lambda = StARs$lambda
opt.index = StARs$index

# transformation
adj_matrix <- (abs(best_precision_matrix_clime_contaminate) !=0) * 1

# create labels of GICS for every stock
gics_map <- setNames(gics_data$GICS.Sector, gics_data$Symbol)
valid_tickers <- colnames(log_returns)
gics_for_valid <- gics_map[valid_tickers]
gics_for_valid <- factor(gics_for_valid)

# create graph
graph <- graph_from_adjacency_matrix(adj_matrix, mode="undirected", diag=FALSE)
# label every single sector
V(graph)$sector <- gics_for_valid

# define colors
sector_colors <- c("red", "blue", "green", "purple", "orange", "brown", "pink", "yellow", "gray", "cyan","black")
names(sector_colors) <- unique(gics_for_valid)
V(graph)$color <- sector_colors[V(graph)$sector]

# grid
ggraph(graph, layout = "fr") +
  geom_edge_link(alpha = 0.5, color = "darkgray") +
  geom_node_point(aes(color = sector), size = 1.5) +
  scale_color_manual(values = sector_colors) +
  theme_void() +
  ggtitle("CLIME")+
  theme(plot.title = element_text(hjust = 0.5))  


# Recording codes(Please ignore them)

#adj_matrix_sglasso_contaminate <- adj_matrix
#adj_matrix_glasso_contaminate <- adj_matrix
#best_precision_matrix <- best_precision_matrix_clime
#best_precision_matrix_clime_normal <- best_precision_matrix
#best_lambda_clime_normal = opt.lambda
#best_precision_matrix_clime_contaminate <- best_precision_matrix
#best_lambda_clime_contaminate = opt.lambda
#adj_matrix_clime_normal <- (abs(best_precision_matrix_clime_normal) !=0) * 1
#adj_matrix_clime_contaminate <- (abs(best_precision_matrix_clime_contaminate) !=0) * 1
#sum(adj_matrix_clime_normal != adj_matrix_clime_contaminate)




#layout_SGLASSO <- create_layout(graph_SGLASSO_contaminate, layout = "fr")
#graph_CLIME_contaminate <- graph_from_adjacency_matrix(adj_matrix_clime_contaminate, mode="undirected", diag=FALSE)
#V(graph_CLIME_contaminate)$sector <- gics_for_valid
#graph_SGLASSO_contaminate <- graph_from_adjacency_matrix(adj_matrix_sglasso_contaminate, mode="undirected", diag=FALSE)
#V(graph_SGLASSO_contaminate)$sector <- gics_for_valid
#xmin = min(layout_SGLASSO$x)
#xmax = max(layout_SGLASSO$x)
#ymin = min(layout_SGLASSO$y)
#ymax = max(layout_SGLASSO$y)
#library(Matrix)
#sparse_adj_matrix_clime_contaminate <- Matrix(adj_matrix_clime_contaminate, sparse=TRUE)
#writeMM(sparse_adj_matrix_clime_contaminate, "clime_contaminate_matrix.mtx")
