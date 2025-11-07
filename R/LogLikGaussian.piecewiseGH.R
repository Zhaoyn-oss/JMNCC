LogLikGaussian.piecewisefJM <- function(thetas){
  thetas <- relist(thetas, skeleton = list.thetas)
  gammas <- thetas$gammas
  sigma <- exp(thetas$log.sigma)
  betas <- thetas$betas 
  alphas <- thetas$alphas
  xi <- if(!CompRisk) exp(thetas$log.xi) else lapply(thetas$log.xi, exp) 
  D <- thetas$D
  D <- if (!diag.D) exp(D) else chol.transf(D)
  
  L <- chol(D) 
  bm <- sqrt(2) * t(L %*% t(tao)) 
  dimnames(bm) <- NULL
  
  # log f(Yi|bm, Xi)
  mu.y <- drop(X%*%gammas)+Z%*%t(bm)
  logNorm <- dnorm(y, mu.y, sigma, log=TRUE)
  logNorm[is.na(logNorm)] <- 0 
  log.p.yb <- rowsum(logNorm, id, reorder = FALSE); dimnames(log.p.yb) <- NULL
  
  #log f(Ti, deltai|bm, Xi)
  Ytime <- drop(Xtime %*% gammas) + Ztime%*%t(bm)
  eta.tw <- lapply(seq_len(nRisks), function(k) if (!is.null(WW)) drop(WW %*% alphas[[k]]) else 0) 
  exp.eta.tw <- lapply(eta.tw, exp)
  log.hazard <- lapply(seq_len(nRisks), function(k) (delta == k) * (log(xi[[k]][ind.D]) + eta.tw[[k]] + betas[k] * Ytime))  # log lambda_k(Ti|bm)
  
  Ys <- drop(Xs %*% gammas) + Zs%*%t(bm)
  eta.s <- lapply(seq_len(nRisks), function(k) betas[k] * Ys)
  exp.eta.s <- lapply(eta.s, exp) 
  xi.wkP.s <- lapply(seq_len(nRisks), function(k) xi[[k]][ind.K] * wkP * exp.eta.s[[k]]) 
  log.survival <- lapply(seq_len(nRisks), function(k) -exp.eta.tw[[k]]*rowsum(xi.wkP.s[[k]], id.GK, reorder = FALSE)) 
  log.p.tb <- Reduce(`+`, log.hazard)+Reduce(`+`, log.survival)
  
  log.wGH <- log(wGH)
  log.wGH.mat <- matrix(rep(log.wGH, each = n), nrow = n, ncol = M)
  log.res <- log.p.tb+log.p.yb+log.wGH.mat
  mx <- apply(log.res, 1, max) 
  p.yt <- rowSums(exp(log.res - mx)) 
  log.p.yt <- log(p.yt)+ mx
  
  res <- -sum(log.p.yt[is.finite(log.p.yt)], na.rm = TRUE) 
  return(res) 
}


LogLikGaussian.piecewisewJM <- function(thetas){
  thetas <- relist(thetas, skeleton = list.thetas)
  gammas <- thetas$gammas
  sigma <- exp(thetas$log.sigma)
  betas <- thetas$betas 
  alphas <- thetas$alphas
  xi <- if(!CompRisk) exp(thetas$log.xi) else lapply(thetas$log.xi, exp) 
  D <- thetas$D
  D <- if (!diag.D) exp(D) else chol.transf(D)
  
  L <- chol(D) 
  bm <- sqrt(2) * t(L %*% t(tao)) 
  dimnames(bm) <- NULL
  
  # log f(Yi|bm, Xi)
  mu.y <- drop(X%*%gammas)+Z%*%t(bm)
  logNorm <- dnorm(y, mu.y, sigma, log=TRUE)
  logNorm[is.na(logNorm)] <- 0 
  log.p.yb <- rowsum(logNorm, id, reorder = FALSE); dimnames(log.p.yb) <- NULL
  
  #log f(Ti, deltai|bm, Xi)
  Ytime <- drop(Xtime %*% gammas) + Ztime%*%t(bm)
  eta.tw <- lapply(seq_len(nRisks), function(k) if (!is.null(WW)) drop(WW %*% alphas[[k]]) else 0) 
  exp.eta.tw <- lapply(eta.tw, exp)
  log.hazard <- lapply(seq_len(nRisks), function(k) (delta == k) * (log(xi[[k]][ind.D]) + eta.tw[[k]] + betas[k] * Ytime))  # log lambda_k(Ti|bm)
  
  Ys <- drop(Xs %*% gammas) + Zs%*%t(bm)
  eta.s <- lapply(seq_len(nRisks), function(k) betas[k] * Ys)
  exp.eta.s <- lapply(eta.s, exp) 
  xi.wkP.s <- lapply(seq_len(nRisks), function(k) xi[[k]][ind.K] * wkP * exp.eta.s[[k]]) 
  log.survival <- lapply(seq_len(nRisks), function(k) -exp.eta.tw[[k]]*rowsum(xi.wkP.s[[k]], id.GK, reorder = FALSE)) 
  log.p.tb <- Reduce(`+`, log.hazard)+Reduce(`+`, log.survival)
  
  log.wGH <- log(wGH)
  log.wGH.mat <- matrix(rep(log.wGH, each = n), nrow = n, ncol = M)
  log.res <- log.p.tb+log.p.yb+log.wGH.mat
  mx <- apply(log.res, 1, max) 
  p.yt <- rowSums(exp(log.res - mx)) 
  log.p.yt <- (log(p.yt)+ mx)*weight
  
  res <- -sum(log.p.yt[is.finite(log.p.yt)], na.rm = TRUE) 
  return(res) 
}
