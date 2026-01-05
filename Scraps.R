############################################################
# Fri Nov 14 08:22:03 2025
# Single subject PCA with a regularization towards group maps
# All ICA replaced by PCA/SVD everywhere
############################################################

library(CICA)
library(multiway)

############################################################
##### 1. Generate true group spatial maps S_true (V x K)
############################################################

set.seed(1234)
Xe <- Sim_CICA(Nr = 15, Q = 5, R = 2, voxels = 100, timepoints = 10,
               E = 0.4, overlap = .50, externalscore = TRUE)

Sr_true <- Xe$Sr[[1]] + Xe$Sr[[2]]    # average = group maps
K <- ncol(Sr_true)
S_subj <- length(Xe$X)

############################################################
##### 2. Group template S_G (in practice from group-PCA/SCA)
############################################################

# Use Sr_true as the template
S_G <- Sr_true

############################################################
##### 3. PCA-based penalized factorization for a single subject
############################################################

pca_ls_spatial_penalty <- function(X, S_G, lambda, K = ncol(S_G), n_iter = 20, ortho = TRUE) {
  
  V  <- nrow(X)
  Tn <- ncol(X)
  
  #### INITIALIZATION BY PCA/SVD ####
  sv <- svd(X, nu = K, nv = K)
  S  <- sv$u[,1:K] %*% diag(sv$d[1:K])     # V x K spatial components
  A  <- t(sv$v[,1:K])                      # K x T temporal scores
  
  for (it in 1:n_iter) {
    
    #### 1) Update S with penalty ####
    AA_t <- A %*% t(A)              # K x K
    XA_t <- X %*% t(A)              # V x K
    
    # S = (X A^T + λ S_G)(A A^T + λ I)^(-1)
    S <- (XA_t + lambda * S_G) %*% solve(AA_t + lambda * diag(K))
    
    #### 2) Update A (plain LS) ####
    StS <- t(S) %*% S
    SX  <- t(S) %*% X
    A   <- solve(StS, SX)
  }
  
  ### Do orthogonalisation
  
  if(ortho == TRUE){
    qrS <- qr(S)
    Q <- qr.Q(qrS)
    R <- qr.R(qrS)
    
    S <- Q
    A <- R %*% A
  }
  
  return(list(S = S, A = A))
}

#### quick test for ortho #####

test <- pca_ls_spatial_penalty(X = Xe$X[[1]], S_G = S_G, lambda = 100,
                               K = K, ortho = TRUE)

t(test$S) %*% test$S
round(t(test$S) %*% test$S, 5)
round( max(abs( t(test$S) %*% test$S - diag(ncol(test$S)) )) , 10)


############################################################
##### 4. Run model for several lambda values
############################################################

set.seed(42)
alphas  <- c(0, 0.1, 10, 1000)
results <- vector("list", length(alphas))
names(results) <- paste0("lambda_", alphas)

frobenius_norm <- function(M) sqrt(sum(M^2))

for (j in seq_along(alphas)) {
  
  lambda <- alphas[j]
  
  S_hat_list    <- vector("list", S_subj)
  spatial_err   <- numeric(S_subj)
  recon_err     <- numeric(S_subj)
  
  for (i in 1:S_subj) {
    
    X_i <- Xe$X[[i]]
    #fit_i <- penalized_pca_ica(X = X_i, S_G = S_G, K = Q, lambda = lambda)
    fit_i <- pca_ls_spatial_penalty(X_i, S_G, lambda, K = K, n_iter = 50, ortho = T)
    S_hat <- fit_i$S
    A_hat <- fit_i$A
    
    S_hat_list[[i]] <- S_hat
    
    # Spatial error relative to true group
    spatial_err[i] <- frobenius_norm(S_hat - S_G)
    
    # Reconstruction error
    recon_err[i] <- mean((X_i - S_hat %*% A_hat)^2)
  }
  
  results[[j]] <- list(
    S_hat_list  = S_hat_list,
    spatial_err = spatial_err,
    recon_err   = recon_err
  )
}

############################################################
##### 5. Inspect results
############################################################

cat("=== Penalized PCA/SVD with group spatial maps ===\n")
for (j in seq_along(alphas)) {
  lambda <- alphas[j]
  cat("\nlambda =", lambda, "\n")
  cat("  Mean spatial error (Frobenius): ",
      round(mean(results[[j]]$spatial_err), 4), "\n")
  cat("  Mean reconstruction MSE: ",
      round(mean(results[[j]]$recon_err), 4), "\n")
}

############################################################
##### 6. compute modified RV matrices across lambdas
##### change names depending on chosen lambda
############################################################

Ss <- c(list(S_G), results$lambda_0$S_hat_list)
mod_0 <- CICA::computeRVmat(Ss, dist = TRUE, verbose = TRUE)

Ss <- c(list(S_G), results$lambda_0.1$S_hat_list)
mod_01 <- CICA::computeRVmat(Ss, dist = TRUE, verbose = TRUE)

Ss <- c(list(S_G), results$lambda_10$S_hat_list)
mod_10 <- CICA::computeRVmat(Ss, dist = TRUE, verbose = TRUE)

Ss <- c(list(S_G), results$lambda_1000$S_hat_list)
mod_1000 <- CICA::computeRVmat(Ss, dist = TRUE, verbose = TRUE)

round(mod_0, 5)
round(mod_01, 5)
round(mod_10, 5)
round(mod_1000, 5)

############################################################
##### 7. Congruence examples
############################################################

round(congru(S_G, results$lambda_0$S_hat_list[[1]]), 3)
round(congru(S_G, results$lambda_1000$S_hat_list[[1]]), 3)
round(congru(S_G, results$lambda_10$S_hat_list[[1]]), 3)

############################################################
##### 8. MDS / clustering plots
############################################################

par(mfrow=c(2,2))
plot(cmdscale(mod_0), col = c(1, rep(2,30)),
     xlim=c(-0.5,0.5), ylim=c(-0.5,0.5))
plot(cmdscale(mod_01), col = c(1, rep(2,30)),
     xlim=c(-0.5,0.5), ylim=c(-0.5,0.5))
plot(cmdscale(mod_10), col = c(1, rep(2,30)),
     xlim=c(-0.5,0.5), ylim=c(-0.5,0.5))
plot(cmdscale(mod_1000), col = c(1, rep(2,30)),
     xlim=c(-0.5,0.5), ylim=c(-0.5,0.5))

plot(cmdscale(mod_1000), col = c(1, rep(2,30)))

par(mfrow=c(1,1))
hcl <- hclust(mod_1000)
plot(hcl)
