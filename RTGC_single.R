# Mon Jan  5 10:47:54 2026 ------------------------------
# Single subject RTGC function

#### input arguments: ####
# X: a single matrix
# Q: number of components (uses the number of components from the template)
# S_G: template components (e.g., group component matrix)
# lambda 
# n_iter: number of ALS type of iterations for updating S and A
# ortho: orthogonalization of the components
# center: centering of X
# scale: scaling of X

RTGC_single <- function(X, Q = ncol(S_G), S_G, lambda, n_iter = 50, ortho = TRUE, center = TRUE, scale = TRUE){
  
  # Check if centering or scaling is necessary
  if (center || scale) {
    X <- scale(X, center = center, scale = scale)
  }
  
  
  # Initialization
  
  sv <- svd(X, nu = Q, nv = Q)
  S  <- sv$u[,1:Q] %*% diag(sv$d[1:Q])     # V x Q spatial components
  A  <- t(sv$v[,1:Q])                      # Q x T temporal components
  
  
  for (it in 1:n_iter){
    
    #### 1) Update S with penalty ####
    AA_t <- A %*% t(A)              # Q x Q
    XA_t <- X %*% t(A)              # V x Q
    
    # S = (X A^T + λ S_G)(A A^T + λ I)^(-1)
    S <- (XA_t + lambda * S_G) %*% solve(AA_t + lambda * diag(Q))
    
    
    #### 2) Update A (plain LS) ####
    StS <- t(S) %*% S
    SX  <- t(S) %*% X
    A   <- solve(StS, SX)
    
  }
  
  if(ortho == TRUE){
    qrS <- qr(S)
    Q <- qr.Q(qrS)
    R <- qr.R(qrS)
    
    S <- Q
    A <- R %*% A
  }
  
  
  return(list(S = S, A = A))
  
}



pca_ls_spatial_penalty <- function(X, S_G, lambda, K = ncol(S_G), n_iter = 50, ortho = TRUE) {
  
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


##### minimal test ######
library(CICA)
library(multiway)

set.seed(1234)
Xe <- Sim_CICA(Nr = 15, Q = 5, R = 2, voxels = 100, timepoints = 10,
               E = 0.4, overlap = .50, externalscore = TRUE)

Sr_true <- Xe$Sr[[1]] + Xe$Sr[[2]]    # average = group maps
S_G <- scale(Sr_true)

test <- RTGC_single(X = Xe$X[[1]], S_G = S_G, lambda = 1)

t(test$S) %*% test$S
round(t(test$S) %*% test$S, 5)
round( max(abs( t(test$S) %*% test$S - diag(ncol(test$S)) )) , 10)

sum((test$S - S_G)^2)

seq <- c(0, .001, .01, 1, 10, 100, 1000, 10000)

for(i in 1:length(seq)){
  test <- RTGC_single(X = Xe$X[[1]], S_G = S_G, lambda = seq[i], ortho = T)
  
  t(test$S) %*% test$S
  round(t(test$S) %*% test$S, 5)
  round( max(abs( t(test$S) %*% test$S - diag(ncol(test$S)) )) , 10)
  
  print( sum((test$S - S_G)^2) )
}

