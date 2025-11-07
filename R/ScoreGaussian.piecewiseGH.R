ScoreGaussian.piecewisefJM <- function(thetas){
  thetas <- relist(thetas, skeleton = list.thetas)
  betas <- thetas$betas
  sigma <- exp(thetas$log.sigma)
  gammas <- thetas$gammas
  alphas <- thetas$alphas 
  xi <- if(!CompRisk) exp(thetas$log.xi) else lapply(thetas$log.xi, exp) 
  D <- thetas$D
  D <- if (!diag.D) exp(D) else chol.transf(D)
  
  L <- chol(D) 
  bm <- sqrt(2) * t(L %*% t(tao)) 
  dimnames(bm) <- NULL
  
  mu.y <- drop(X%*%gammas)+Z%*%t(bm)
  logNorm <- dnorm(y, mu.y, sigma, log=TRUE)
  logNorm[is.na(logNorm)] <- 0 
  log.p.yb <- rowsum(logNorm, id, reorder = FALSE); dimnames(log.p.yb) <- NULL
  
  Ytime <- drop(Xtime %*% gammas) + Ztime%*%t(bm)
  eta.tw <- lapply(seq_len(nRisks), function(k) if (!is.null(WW)) drop(WW %*% alphas[[k]]) else 0) 
  exp.eta.tw <- lapply(eta.tw, exp)
  log.hazard <- lapply(seq_len(nRisks), function(k) (delta == k) * (log(xi[[k]][ind.D]) + eta.tw[[k]] + betas[k] * Ytime))  # log lambda_k(Ti|bm)
  
  Ys <- drop(Xs %*% gammas) + Zs%*%t(bm)
  eta.s <- lapply(seq_len(nRisks), function(k) betas[k] * Ys)
  exp.eta.s <- lapply(eta.s, exp) 
  xi.wkP.s <- lapply(seq_len(nRisks), function(k) xi[[k]][ind.K] * wkP * exp.eta.s[[k]]) 
  log.survival <- lapply(seq_len(nRisks), function(k) exp.eta.tw[[k]]*rowsum(xi.wkP.s[[k]], id.GK, reorder = FALSE)) 
  log.p.tb <- Reduce(`+`, log.hazard)-Reduce(`+`, log.survival)
   
  
  log.wGH <- log(wGH)
  log.wGH.mat <- matrix(rep(log.wGH, each = n), nrow = n, ncol = M)
  log.res <- log.p.tb+log.p.yb+log.wGH.mat
  mx <- apply(log.res, 1, max) 
  p.ytb <- exp(log.res - mx)
  p.yt <- rowSums(p.ytb) 
  p.byt <- p.ytb/p.yt 
  
  # for gamma
  sc1 <- numeric(ncx)
  for (i in seq_len(ncx)) { 
    rs <- (X[, i]*(y-mu.y))/sigma^2
    rs[is.na(rs)] <- 0
    tmp <- rowsum(rs, id, reorder = FALSE) 
    sc1[i] <- sum(tmp * p.byt, na.rm = TRUE)
  } 
  ki.haz <- 0
  for(k in seq_len(nRisks)){
    ki.haz <- ki.haz+ (betas[k]*Xtime)*(delta==k) 
  } 
  sc2 <- colSums(ki.haz, na.rm = TRUE)
  
  sc3 <- numeric(ncx)
  for(i in seq_len(ncx)){
    tmp <- 0
    for(k in 1:nRisks){
      rs <- exp.eta.tw[[k]]*rowsum(xi.wkP.s[[k]]*betas[k]*Xs[, i], id.GK, reorder = FALSE)
      tmp <- tmp+rs * p.byt
    } 
    sc3[i] <- sum(tmp, na.rm = TRUE)
  }
  score.gammas <- -sc1-sc2+sc3
   
  # sigma: in longitudinal sub-model
  b2 <- if (ncz == 1) bm * bm else t(apply(bm, 1, function (x) x %o% x))
  post.b <- p.byt %*% bm
  Zb <- if (ncz == 1) post.b[id] else rowSums(Z * post.b[id, , drop = FALSE], na.rm = TRUE)  
  
  post.vb <- if (ncz == 1) {     
    drop(p.byt %*% b2) - drop(post.b^2) 
  } else {
    (p.byt %*%b2 ) - t(apply(post.b, 1, function(x) x %o% x))
  }
  eta.yx <- drop(X %*% gammas)
  tr.tZZvarb <- rowSums(ZtZ*post.vb, na.rm = TRUE) 
  resid <- (y - eta.yx-Zb)^2; resid[is.na(resid)] <- 0
  resid.subject <- rowsum(resid, id, reorder = FALSE)
  resid.subject[sub.id] <- NA
  scsigma2 <- (resid.subject+tr.tZZvarb)/sigma^2 
  score.sigma <- length(obs.id)-sum(scsigma2, na.rm=TRUE) 
  score.y <- c(score.gammas, score.sigma)
  
  # score for beta and alpha: in survival sub-model
  scbeta <- lapply(seq_len(nRisks), function(k){
    ki <- (delta==k)*Ytime- exp.eta.tw[[k]]*rowsum(xi.wkP.s[[k]]* Ys, id.GK, reorder = FALSE)
    res <- -p.byt*ki
    return(res)
  }) 
  score.beta <- unlist(lapply(scbeta, sum, na.rm=TRUE))
  
  scalpha <- lapply(seq_len(nRisks), function(k) {
    if (!is.null(WW)) -colSums(WW*((delta==k)-rowSums(p.byt * (exp.eta.tw[[k]] * rowsum(xi.wkP.s[[k]], id.GK, reorder = FALSE)), na.rm = TRUE)), na.rm = TRUE) else NULL
  }) 
  score.alpha <- unlist(scalpha)
 
  scxi <- lapply(seq_len(nRisks), function(k){ 
    scxi <- matrix(0, nrow=n, ncol=Q) 
    for(i in 1:Q){
      i1 <- ind.D == i
      i2 <- ind.K == i
      i3 <- ind.D >= i
      ki <- rowSums((p.byt[i3, ] * (exp.eta.tw[[k]][i3] * rowsum((wkP*exp.eta.s[[k]])[i2, ], id.GK[i2], reorder = FALSE))), na.rm = TRUE)
      kk <- numeric(n); kk[i3] <- ki 
      par1 <- numeric(n); par1[i1] <- ((delta == k))[i1]
      scxi[, i] <- - xi[[k]][i] *(par1/xi[[k]][i]-kk) 
    }
    return(scxi)
  })
  
  scxi <- lapply(scxi, function(x){
    res <- if(is.null(dim(x))) sum(x, na.rm = TRUE) else colSums(x, na.rm = TRUE)
    return(res)
  }) 
  score.xi <- unlist(scxi)
  score.t <- c(score.beta, score.alpha, score.xi)
  
  # for D  
  V  <- matrix(0, ncz, M)
  for(i in 1:ncz){
    rs <- (Z[, i]*(y-mu.y))/sigma^2
    rs[is.na(rs)] <- 0
    rs_long <- rowsum(rs, id, reorder = FALSE) * p.byt
    V[i, ] <- V[i, ] - colSums(rs_long , na.rm = TRUE)
   
    haz_i <- 0
    for (k in seq_len(nRisks)) haz_i <- haz_i + (betas[k] * Ztime[, i]) * (delta == k)    
    haz_i_mat <- matrix(rep(haz_i, M), n, M)* p.byt
    V[i, ] <- V[i, ] - colSums(haz_i_mat , na.rm = TRUE)
   
    surv_i_mat <- matrix(0, n, M)
    for (k in seq_len(nRisks)) {  
      tmp <- matrix(0, n, M)
      for (m in seq_len(M)) {
        w_km  <- xi.wkP.s[[k]][, m]    
        tmp[, m] <- rowsum(w_km * betas[k] * Zs[, i], id.GK, reorder = FALSE)
      }
      surv_i_mat <- surv_i_mat + (exp.eta.tw[[k]] * tmp) * p.byt
    }
    V[i, ] <- V[i, ] + colSums(surv_i_mat , na.rm = TRUE)
  }
  dL <- sqrt(2) * (V %*% tao) 
   
  if (!diag.D) {
    score.b <- 0.5 * diag(dL) * diag(L)
  } else {
    idx <- upper.tri(L, diag = TRUE)
    dL_vec <- dL[idx]
    J <- jacobian2(t(L)) 
    score.b <- drop(dL_vec %*% J)
  }
   
  c(score.y, score.t, score.b)
}
 

ScoreGaussian.piecewisewJM <- function(thetas){
  thetas <- relist(thetas, skeleton = list.thetas)
  betas <- thetas$betas
  sigma <- exp(thetas$log.sigma)
  gammas <- thetas$gammas
  alphas <- thetas$alphas 
  xi <- if(!CompRisk) exp(thetas$log.xi) else lapply(thetas$log.xi, exp) 
  D <- thetas$D
  D <- if (!diag.D) exp(D) else chol.transf(D)
  
  L <- chol(D) 
  bm <- sqrt(2) * t(L %*% t(tao)) 
  dimnames(bm) <- NULL
  
  mu.y <- drop(X%*%gammas)+Z%*%t(bm)
  logNorm <- dnorm(y, mu.y, sigma, log=TRUE)
  logNorm[is.na(logNorm)] <- 0 
  log.p.yb <- rowsum(logNorm, id, reorder = FALSE); dimnames(log.p.yb) <- NULL
  
  Ytime <- drop(Xtime %*% gammas) + Ztime%*%t(bm)
  eta.tw <- lapply(seq_len(nRisks), function(k) if (!is.null(WW)) drop(WW %*% alphas[[k]]) else 0) 
  exp.eta.tw <- lapply(eta.tw, exp)
  log.hazard <- lapply(seq_len(nRisks), function(k) (delta == k) * (log(xi[[k]][ind.D]) + eta.tw[[k]] + betas[k] * Ytime))  # log lambda_k(Ti|bm)
  
  Ys <- drop(Xs %*% gammas) + Zs%*%t(bm)
  eta.s <- lapply(seq_len(nRisks), function(k) betas[k] * Ys)
  exp.eta.s <- lapply(eta.s, exp) 
  xi.wkP.s <- lapply(seq_len(nRisks), function(k) xi[[k]][ind.K] * wkP * exp.eta.s[[k]]) 
  log.survival <- lapply(seq_len(nRisks), function(k) exp.eta.tw[[k]]*rowsum(xi.wkP.s[[k]], id.GK, reorder = FALSE)) 
  log.p.tb <- Reduce(`+`, log.hazard)-Reduce(`+`, log.survival)
  
  log.wGH <- log(wGH)
  log.wGH.mat <- matrix(rep(log.wGH, each = n), nrow = n, ncol = M)
  log.res <- log.p.tb+log.p.yb+log.wGH.mat
  mx <- apply(log.res, 1, max) 
  p.ytb <- exp(log.res - mx)
  p.yt <- rowSums(p.ytb) 
  p.byt <- p.ytb/p.yt
  
  # for gamma
  sc1 <- numeric(ncx)
  for (i in seq_len(ncx)) { 
    rs <- (X[, i]*(y-mu.y))/sigma^2
    rs[is.na(rs)] <- 0
    tmp <- rowsum(rs, id, reorder = FALSE) 
    sc1[i] <- sum(tmp * p.byt*weight, na.rm = TRUE)
  } 
  ki.haz <- 0
  for(k in seq_len(nRisks)){
    ki.haz <- ki.haz+ (betas[k]*Xtime)*(delta==k) 
  } 
  sc2 <- colSums(weight*ki.haz, na.rm = TRUE)
  
  sc3 <- numeric(ncx)
  for(i in seq_len(ncx)){
    tmp <- 0
    for(k in 1:nRisks){
      rs <- exp.eta.tw[[k]]*rowsum(xi.wkP.s[[k]]*betas[k]*Xs[, i], id.GK, reorder = FALSE)
      tmp <- tmp+rs * p.byt*weight
    } 
    sc3[i] <- sum(tmp, na.rm = TRUE)
  }
  score.gammas <- -sc1-sc2+sc3
  
  # sigma: in longitudinal sub-model
  b2 <- if (ncz == 1) bm * bm else t(apply(bm, 1, function (x) x %o% x))
  post.b <- p.byt %*% bm
  Zb <- if (ncz == 1) post.b[id] else rowSums(Z * post.b[id, , drop = FALSE], na.rm = TRUE)  
  
  post.vb <- if (ncz == 1) {     
    drop(p.byt %*% b2) - drop(post.b^2) 
  } else {
    (p.byt %*%b2 ) - t(apply(post.b, 1, function(x) x %o% x))
  }
  eta.yx <- drop(X %*% gammas)
  tr.tZZvarb <- rowSums(ZtZ*post.vb, na.rm = TRUE) 
  weights.obs <- weight[id]
  weights.obs[-obs.id] <- 0
  scsigma1 <- rowsum(weights.obs, id, reorder = FALSE)
  resid <- (y - eta.yx-Zb)^2; resid[is.na(resid)] <- 0
  resid.subject <- rowsum(resid, id, reorder = FALSE)
  resid.subject[sub.id] <- NA
  scsigma2 <- weight*(resid.subject+tr.tZZvarb)/sigma^2 
  score.sigma <- sum(scsigma1-scsigma2, na.rm=TRUE)
  score.y <- c(score.gammas, score.sigma)
  
  # score for beta and alpha: in survival sub-model
  scbeta <- lapply(seq_len(nRisks), function(k){
    ki <- (delta==k)*Ytime- exp.eta.tw[[k]]*rowsum(xi.wkP.s[[k]]* Ys, id.GK, reorder = FALSE)
    res <- -weight*(p.byt*ki)
    return(res)
  }) 
  score.beta <- unlist(lapply(scbeta, sum, na.rm=TRUE))
  
  scalpha <- lapply(seq_len(nRisks), function(k) {
    if (!is.null(WW)) -colSums(weight*WW*((delta==k)-rowSums(p.byt * (exp.eta.tw[[k]] * rowsum(xi.wkP.s[[k]], id.GK, reorder = FALSE)), na.rm = TRUE)), na.rm = TRUE) else NULL
  }) 
  score.alpha <- unlist(scalpha)
  
  scxi <- lapply(seq_len(nRisks), function(k){ 
    scxi <- matrix(0, nrow=n, ncol=Q) 
    for(i in 1:Q){
      i1 <- ind.D == i
      i2 <- ind.K == i
      i3 <- ind.D >= i
      ki <- rowSums((p.byt[i3, ] * (exp.eta.tw[[k]][i3] * rowsum((wkP*exp.eta.s[[k]])[i2, ], id.GK[i2], reorder = FALSE))), na.rm = TRUE)
      kk <- numeric(n); kk[i3] <- ki 
      par1 <- numeric(n); par1[i1] <- (weight*(delta == k))[i1]
      scxi[, i] <- - xi[[k]][i] *(par1/xi[[k]][i]-weight*kk) 
    }
    return(scxi)
  })
  
  scxi <- lapply(scxi, function(x){
    res <- if(is.null(dim(x))) sum(x, na.rm = TRUE) else colSums(x, na.rm = TRUE)
    return(res)
  }) 
  score.xi <- unlist(scxi)
  score.t <- c(score.beta, score.alpha, score.xi)
  
  # for D  
  V  <- matrix(0, ncz, M)
  for(i in 1:ncz){
    rs <- (Z[, i]*(y-mu.y))/sigma^2
    rs[is.na(rs)] <- 0
    rs_long <- rowsum(rs, id, reorder = FALSE) * p.byt*weight 
    V[i, ] <- V[i, ] - colSums(rs_long , na.rm = TRUE) 
    
    haz_i <- 0
    for (k in seq_len(nRisks)) haz_i <- haz_i + (betas[k] * Ztime[, i]) * (delta == k)   
    haz_i_mat <- matrix(rep(haz_i, M), n, M)* p.byt*weight
    V[i, ] <- V[i, ] - colSums(haz_i_mat, na.rm = TRUE)
    
    surv_i_mat <- matrix(0, n, M)
    for (k in seq_len(nRisks)) { 
      tmp <- matrix(0, n, M)
      for (m in seq_len(M)) {
        w_km  <- xi.wkP.s[[k]][, m]   
        tmp[, m] <- rowsum(w_km * betas[k] * Zs[, i], id.GK, reorder = FALSE)
      }
      surv_i_mat <- surv_i_mat + (exp.eta.tw[[k]] * tmp) * p.byt*weight
    }
    V[i, ] <- V[i, ] + colSums(surv_i_mat , na.rm = TRUE)
  }
  dL <- sqrt(2) * (V %*% tao) 
  
  if (!diag.D) {
    score.b <- 0.5 * diag(dL) * diag(L)
  } else {  
    idx <- upper.tri(L, diag = TRUE)
    dL_vec <- dL[idx]
    J <- jacobian2(t(L)) 
    score.b <- drop(dL_vec %*% J)
  }
  
  c(score.y, score.t, score.b)
}


