library(abind)
library(rBeta2009)
set.seed(12754)

#######################
#### Function Code ####
#######################


closure <- function(X) {
  return(apply(X, 2, function(col) col/sum(col)))
}

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

#' BpimC Code
#'
#' @param Y D x N matrix of sequence counts. D taxa, N samples.
#' @param scale_var_sampler A function. Returns a sampled estimate
#'        of scalar variance in scale
#' @param scale_rel_sampler A function. Returns a D-length vector
#'        of estimates of correlation between relative abundances
#'        and scale.
#' @param scale_var_args A list. Arguments to pass to scale_var_sampler 
#' @param scale_rel_args A list. Arguments to pass to scale_rel_sampler 
#' @param nsample A numeric. The number of estimates of corr to sample.
#' @return A D x D x S array. S estimates of D x D correlation matrices
bpimC <- function(Y, scale_var_sampler, scale_rel_sampler,
			    scale_var_args=list(), scale_rel_args=list(),
			    nsample=1000) {
  taxa_names <- row.names(Y)
  all_abs_corr <- c() 
  for(i in 1:nsample) {
    # Sample scale args
    scale_var <- do.call(scale_var_sampler, c(list(n = 1), scale_var_args))
    scale_rel_corr <- do.call(scale_rel_sampler, c(list(n=1), scale_rel_args))
    # Dirichlet Sample Relative Abundances
    P_sample <- apply(Y_genus, 2, function(col) rdirichlet(1, col+0.5))
    # Estimate absolute correlations
    abs_corr <- cov2cor(absCorrEst(P_sample, scale_var, scale_rel_corr))
    all_abs_corr <- abind(all_abs_corr, abs_corr, along=3)
  }
  row.names(all_abs_corr) <- taxa_names
  colnames(all_abs_corr) <- taxa_names
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
    P_sample <- apply(Y, 2, function(col) rdirichlet(1, col+0.5))
    corr_mat <- sparcc_basis(P_sample)
    corr_mats <- abind(corr_mats, corr_mat, along=3)
  }
  return(corr_mats)
}

write_df <- function(taxon_i, abs_corr, res) {
  # Write_dfs
  i <- which(rownames(abs_corr) == taxon_i)
  post <- t(res[i, , ])  # samples x taxa
  df <- data.frame(
    taxon = row.names(abs_corr),
    median = apply(post, 2, median),
    lower = apply(post, 2, quantile, 0.025),
    upper = apply(post, 2, quantile, 0.975),
    true = abs_corr[taxon_i, ]
  )
  return(df)
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

## Absolute correlations
rel_abund <- closure(Y_genus + 0.5)
abs_abund <- sweep(rel_abund, 2, metadata$Average.cell.count..per.gram.of.frozen.feces., FUN="*")
abs_corr <- cor(t(log(abs_abund)))

## SparCC Analysis
sparcc_corr_res <- sparcc(Y_genus, nsample=2000)

# Write_dfs
df <- write_df("Oscillibacter", abs_corr, sparcc_corr_res)
write.csv(df, "../results/sparcc_oscillibacter.csv", row.names=F)
df <- write_df("Prevotella", abs_corr, sparcc_corr_res)
write.csv(df, "../results/sparcc_prevotella.csv", row.names=F)

## Bayesian Sampling Correlation Analysis
scale_var_sampler <- function(n, s_l, s_u) {
  return(runif(n, s_l, s_u))
}
scale_rel_sampler <- function(n, D, m_l, m_u, r_l, r_u) {
  m <- runif(n, m_l, m_u)
  scale_rel_corr <- runif(D, m+r_l, m+r_u)
  return(scale_rel_corr)
}
scale_var_args <- list(s_l=0.4, 0.8)
scale_rel_args <- list(m_l=0, m_u=0.15, r_l=-0.3, r_u=0.3, D=nrow(abs_abund))
bayes_corr_res <- bpimC(Y_genus, scale_var_sampler,
				  scale_rel_sampler, scale_var_args,
				  scale_rel_args, nsample=2000)

df <- write_df("Oscillibacter", abs_corr, bayes_corr_res)
write.csv(df, "../results/bayes_oscillibacter.csv", row.names=F)
df <- write_df("Prevotella", abs_corr, bayes_corr_res)
write.csv(df, "../results/bayes_prevotella.csv", row.names=F)
