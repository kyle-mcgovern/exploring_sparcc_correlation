library(abind)


#######################
#### Function Code ####
#######################


#' Estimates absolute abundance correlations given relative abundances,
#' scale variance, and correlations between each taxon and scale.
#' 
#' @param rel_abund a D x N matrix. D taxa, N samples, relative
#'        abundances cols sum to 1
#' @param scale_var A number. The variane in the scale
#' @param scale_rel_corr A vector. A D-length vector of correlation
#'        between each taxon's abundance and scale.
#' @return An estimated DxD absolute abundance correlation
#'         matrix
absCorrEst <- function(rel_abund, scale_var, scale_rel_corr) {
  D <- nrow(rel_abund)
  # D relative abundance variances
  rel_var <- cbind(apply(log(rel_abund), 1, var))
  # Vector of 1s
  v1 <- cbind(rep(1, D))
  # Relative abundance covariance
  Lambda <- cov(log(t(rel_abund)))
  # Scale/rel sd
  scale_sigma <- sqrt(scale_var)
  rel_sigma <- sqrt(rel_var) 
  # Absolute abundance covariance estimate
  abs_cov <- Lambda + (scale_sigma * (rel_sigma * scale_rel_corr))%*%t(v1)
  abs_cov <- abs_cov + v1%*%t(scale_sigma * (rel_sigma * scale_rel_corr))
  abs_cov <- abs_cov + scale_var*matrix(1, D, D)
  # Return correlation matrix
  return(cov2cor(abs_cov))
}

bayesianAbsCorr <- function(Y, scale_var_sampler, scale_var_args=list(),
			    nsample=1000) {
  all_abs_corr <- c() 
  for(i in 1:nsample) {
    # Sample scale args
    scale_var <- do.call(scale_var_sampler, c(list(n = 1), scale_var_args))
    scale_rel_corr <- do.call(scale_rel_sampler, c(list(n=1), scale_rel_args))
    # Dirichlet Sample Relative Abundances
    P_sample <- apply(Y_genus, 2, function(col) rdirichlet(1, col+0.5))
    # Estimate absolute correlations
    abs_corr <- cov2cor(matrix_cov(P_sample, sigma_perp2, x))
    all_abs_corr <- abind(all_res, R, along=3)
  }
  return(all_abs_corr)
}

sparcc_basis <- function(P) {
  taxa_names <- row.names(P)
  T_mat <- matrix(0, nrow=nrow(P), ncol=nrow(P))
  for (i in 1:nrow(P)) {
    for (j in 1:nrow(P)) {
      T_mat[i, j] <- var(log(P[i,]/P[j,]))
    }
  }
  T_sum <- apply(T_mat, 1, sum)
  A <- c()
  for(i in 1:nrow(P)) {
    D <- nrow(P)
    A[i] <- ((T_sum[i] - (sum(T_sum) / (2*D-2))) / (D-2))
  }
  corr_mat <- matrix(0, nrow=nrow(P), ncol=nrow(P))
  for (i in 1:nrow(P)) {
    for (j in 1:nrow(P)) {
      corr_mat[i, j] <- (A[i] + A[j] - T_mat[i,j]) / (2*sqrt(A[i])*sqrt(A[j]))
    }
  }
  row.names(corr_mat) <- taxa_names
  colnames(corr_mat) <- taxa_names
  return(corr_mat)
}

sparcc <- function(Y, nsample=1000) {
  corr_mats <- c()
  for (iter in 1:nsample) {
    P_sample <- apply(Y, 2, function(col) rdirichlet(1, col+prior))
    corr_mat <- sparcc_basis(P_sample)
    corr_mats <- abind(corr_mats, corr_mat, along=3)
  }
  return(corr_mats)
}

#######################
#### ANALYSIS CODE ####
#######################


## Read Sequence Counts and Metadata
metadata <- read.csv("../data/41586_2017_BFnature24460_MOESM10_ESM.csv")
metadata <- metadata[metadata$Health.status=="Control",]
row.names(metadata) <- metadata$Sample
taxonomy <- read.csv("../data/taxa_assignments.csv", row.names=1, header=T)
Y <- apply(t(read.csv("../data/OTU_nochim.csv", row.names=1)), c(1,2), as.numeric)
Y <- Y[row.names(Y)%in%row.names(taxonomy),]
Y <- Y[,colnames(Y)%in%metadata$Sample]
genus <- taxonomy[row.names(Y), "Genus"]
genus[is.na(genus)] <- "unclassified"
row.names(Y) <- genus
# Combine counts if more than 80% sparsity
Y_genus <- rowsum(Y, group=rownames(Y))
other <- colSums(Y_genus[(rowSums(Y_genus==0)/ncol(Y_genus))>0.8,])
Y_genus <- Y_genus[(rowSums(Y_genus==0)/ncol(Y_genus))<=0.8,]
Y_genus <- rbind(Y_genus, other)
Y_genus <- Y_genus[,colnames(Y_genus)%in%row.names(metadata)]
metadata <- metadata[colnames(Y_genus),]

## SparCC Analysis
sparcc_corr_res <- sparcc(Y_genus, nsample=2000)
print(sparcc_corr_res)
