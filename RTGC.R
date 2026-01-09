# Mon Jan  5 13:56:44 2026 ------------------------------
# Main function for RTGC

RTGC <- function(Xlist, Q = ncol(S_G), S_G, lambda, ortho = TRUE, 
                 center = TRUE, scale = TRUE){
  
  results <- list()
  
  for(i in 1:length(Xlist)){
    results[[i]] <- RTGC_single(Xlist[[i]], Q = Q, S_G = S_G, lambda = lambda,
                           ortho = ortho, center = center, scale = scale)
  }
  return(results)
}

Xl <- do.call(cbind, Xe$X)

sv <- svd(Xl, nu = 5)
S_G <- sv$u[,1:5] %*% diag(sv$d[1:5])
#S_G <- scale(S_G)

test <- RTGC(Xlist = Xe$X, S_G = S_G, lambda = 0.00001, ortho = T, center = T, scale = T)

Shats <- lapply(1:length(test), function(anom) test[[anom]]$S)
Shats <- c(Shats, list(S_G))
mod <- computeRVmat(Shats,dist = T)
#heatmap(mod)
plot(cmdscale(mod), asp = TRUE, col = c(rep(1, 30), 2))
