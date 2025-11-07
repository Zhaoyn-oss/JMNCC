ScoreGaussianMat.piecewisefJM <- function(thetas){
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
  
  #p.ytb <- exp(log.p.tb+log.p.yb)
  #p.yt <- c(p.ytb%*%wGH) 
  #wGH.mat <- matrix(rep(wGH, each = n), nrow = n, ncol = M)
  #p.byt <- (p.ytb * wGH.mat)/p.yt
  
  log.wGH <- log(wGH)
  log.wGH.mat <- matrix(rep(log.wGH, each = n), nrow = n, ncol = M)
  log.res <- log.p.tb+log.p.yb+log.wGH.mat
  mx <- apply(log.res, 1, max) 
  p.ytb <- exp(log.res - mx)
  p.yt <- rowSums(p.ytb) 
  p.byt <- p.ytb/p.yt
 
   # for gamma
  sc1 <- sc2 <- sc3 <- matrix(0, nrow = n, ncol = ncx)
  for (i in seq_len(ncx)) { 
    rs <- (X[, i]*(y-mu.y))/sigma^2
    rs[is.na(rs)] <- 0
    tmp <- rowsum(rs, id, reorder = FALSE) * p.byt
    sc1[, i] <- rowSums(tmp, na.rm = TRUE) 
  }  
  for(k in seq_len(nRisks)){
    sc2 <- sc2+ (betas[k]*Xtime)*(delta==k) 
  }  
  for(i in seq_len(ncx)){
    tmp <- 0
    for(k in 1:nRisks){
      rs <- exp.eta.tw[[k]]*rowsum(xi.wkP.s[[k]]*betas[k]*Xs[, i], id.GK, reorder = FALSE)
      tmp <- tmp+rs * p.byt
    } 
    sc3[, i] <- rowSums(tmp, na.rm = TRUE)
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
  weights <- rep(1, n)
  weights.obs <- weights[id]
  weights.obs[-obs.id] <- 0
  scsigma1 <- rowsum(weights.obs, id, reorder = FALSE)
  resid <- (y - eta.yx-Zb)^2; resid[is.na(resid)] <- 0
  resid.subject <- rowsum(resid, id, reorder = FALSE)
  resid.subject[sub.id] <- NA
  scsigma2 <- (resid.subject+tr.tZZvarb)/sigma^2
  scsigma <- scsigma1-scsigma2 
  scsigma[is.na(scsigma)] <- 0   
  score.sigma <- scsigma
  score.y <- cbind(score.gammas, score.sigma)
  
  # score for beta and alpha: in survival sub-model
  scbeta <- lapply(seq_len(nRisks), function(k){
    ki <- (delta==k)*Ytime- exp.eta.tw[[k]]*rowsum(xi.wkP.s[[k]]* Ys, id.GK, reorder = FALSE)
    res <- -rowSums(p.byt*ki, na.rm = TRUE) 
    return(res)
  }) 
  score.beta <- do.call(cbind, scbeta)
  
  scalpha <- lapply(seq_len(nRisks), function(k) {
    if (!is.null(WW)){
      -WW*((delta==k)-rowSums(p.byt * (exp.eta.tw[[k]] * rowsum(xi.wkP.s[[k]], id.GK, reorder = FALSE)), na.rm = TRUE))
    } else NULL
  }) 
  score.alpha <- do.call(cbind, scalpha)
  
  scxi <- lapply(seq_len(nRisks), function(k){ 
    scxi <- matrix(0, nrow=n, ncol=Q) 
    for(i in 1:Q){
      i1 <- ind.D == i
      i2 <- ind.K == i
      i3 <- ind.D >= i
      ki <- rowSums(p.byt[i3, ] * (exp.eta.tw[[k]][i3] * rowsum((wkP*exp.eta.s[[k]])[i2, ], id.GK[i2], reorder = FALSE)), na.rm=TRUE)
      kk <- numeric(n); kk[i3] <- ki 
      par1 <- numeric(n); par1[i1] <- ((delta == k))[i1]
      scxi[, i] <- - xi[[k]][i] *(par1/xi[[k]][i]-kk) 
    }
    return(scxi)
  })
  
  score.xi <- do.call(cbind, scxi) 
  
  score.t <- cbind(score.beta, score.alpha, score.xi)
  
  # for D   
  V_i <- array(0, dim = c(n, ncz, M))   
  for (j in seq_len(ncz)) {
    rs <- (Z[, j]*(y - mu.y)) / sigma^2
    rs[is.na(rs)] <- 0
    rs_long <- rowsum(rs, id, reorder = FALSE) * p.byt  
    V_i[, j, ] <- V_i[, j, ] - rs_long
    
    haz_j <- 0
    for (k in seq_len(nRisks)) {
      haz_j <- haz_j + (betas[k] * Ztime[, j]) * (delta == k)
    }
    haz_j_mat <- matrix(rep(haz_j, M), n, M) * p.byt  
    V_i[, j, ] <- V_i[, j, ] - haz_j_mat
    
    surv_j_mat <- matrix(0, n, M)
    for (k in seq_len(nRisks)) {
      tmp <- matrix(0, n, M)
      for (m in seq_len(M)) {
        w_km  <- xi.wkP.s[[k]][, m]
        tmp[, m] <- rowsum(w_km * betas[k] * Zs[, j], id.GK, reorder = FALSE)
      }
      surv_j_mat <- surv_j_mat + (exp.eta.tw[[k]] * tmp) * p.byt 
    }
    V_i[, j, ] <- V_i[, j, ] + surv_j_mat
  }
  
  dL_i_list <- vector("list", n)
  for (ii in seq_len(n)) {
    V_sub <- V_i[ii, , , drop = FALSE]   # 1 × ncz × M
    V_sub <- matrix(V_sub, ncz, M)
    dL_i_list[[ii]] <- sqrt(2) * (V_sub %*% tao)
  }
  
  score.b_i <- vector("list", n)
  for (ii in seq_len(n)) {
    dL_i <- dL_i_list[[ii]]
    if (!diag.D) {
      score.b_i[[ii]] <- 0.5 * diag(dL_i) * diag(L)
    } else {
      idx <- upper.tri(L, diag = TRUE)
      dL_vec_i <- dL_i[idx]
      J <- jacobian2(t(L))
      score.b_i[[ii]] <- drop(dL_vec_i %*% J)
    }
  }
  score.b <- do.call(rbind, score.b_i)
  cbind(score.y, score.t, score.b)
}


ScoreGaussianMat.piecewisewJM <- function(thetas){
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
  sc1 <- sc2 <- sc3 <- matrix(0, nrow = n, ncol = ncx)
  for (i in seq_len(ncx)) { 
    rs <- (X[, i]*(y-mu.y))/sigma^2
    rs[is.na(rs)] <- 0
    tmp <- rowsum(rs, id, reorder = FALSE) * p.byt*weight
    sc1[, i] <- rowSums(tmp, na.rm = TRUE) 
  }  
  for(k in seq_len(nRisks)){
    sc2 <- sc2+ (betas[k]*Xtime)*(delta==k)*weight 
  }  
  for(i in seq_len(ncx)){
    tmp <- 0
    for(k in 1:nRisks){
      rs <- exp.eta.tw[[k]]*rowsum(xi.wkP.s[[k]]*betas[k]*Xs[, i], id.GK, reorder = FALSE)
      tmp <- tmp+rs * p.byt*weight
    } 
    sc3[, i] <- rowSums(tmp, na.rm = TRUE)
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
  scsigma <- scsigma1-scsigma2 
  scsigma[is.na(scsigma)] <- 0   
  score.sigma <- scsigma
  score.y <- cbind(score.gammas, score.sigma)
  
  # score for beta and alpha: in survival sub-model
  scbeta <- lapply(seq_len(nRisks), function(k){
    ki <- (delta==k)*Ytime- exp.eta.tw[[k]]*rowsum(xi.wkP.s[[k]]* Ys, id.GK, reorder = FALSE)
    res <- -rowSums(weight*p.byt*ki, na.rm = TRUE) 
    return(res)
  }) 
  score.beta <- do.call(cbind, scbeta)
  
  scalpha <- lapply(seq_len(nRisks), function(k) {
    if (!is.null(WW)){
      -weight*WW*((delta==k)-rowSums(p.byt * (exp.eta.tw[[k]] * rowsum(xi.wkP.s[[k]], id.GK, reorder = FALSE)), na.rm = TRUE))
    } else NULL
  }) 
  score.alpha <- do.call(cbind, scalpha)
  
  scxi <- lapply(seq_len(nRisks), function(k){ 
    scxi <- matrix(0, nrow=n, ncol=Q) 
    for(i in 1:Q){
      i1 <- ind.D == i
      i2 <- ind.K == i
      i3 <- ind.D >= i
      ki <- rowSums(p.byt[i3, ] * (exp.eta.tw[[k]][i3] * rowsum((wkP*exp.eta.s[[k]])[i2, ], id.GK[i2], reorder = FALSE)), na.rm=TRUE)
      kk <- numeric(n); kk[i3] <- ki 
      par1 <- numeric(n); par1[i1] <- (weight*(delta == k))[i1]
      scxi[, i] <- - xi[[k]][i] *(par1/xi[[k]][i]-weight*kk) 
    }
    return(scxi)
  })
  
  score.xi <- do.call(cbind, scxi) 
  
  score.t <- cbind(score.beta, score.alpha, score.xi)
  
  # for D  
  V_i <- array(0, dim = c(n, ncz, M))   
  for (j in seq_len(ncz)) {
    rs <- (Z[, j]*(y - mu.y)) / sigma^2
    rs[is.na(rs)] <- 0
    rs_long <- rowsum(rs, id, reorder = FALSE) * p.byt *weight 
    V_i[, j, ] <- V_i[, j, ] - rs_long
    
    haz_j <- 0
    for (k in seq_len(nRisks)) {
      haz_j <- haz_j + (betas[k] * Ztime[, j]) * (delta == k)
    }
    haz_j_mat <- matrix(rep(haz_j, M), n, M) * p.byt*weight
    V_i[, j, ] <- V_i[, j, ] - haz_j_mat
    
    surv_j_mat <- matrix(0, n, M)
    for (k in seq_len(nRisks)) {
      tmp <- matrix(0, n, M)
      for (m in seq_len(M)) {
        w_km  <- xi.wkP.s[[k]][, m]
        tmp[, m] <- rowsum(w_km * betas[k] * Zs[, j], id.GK, reorder = FALSE)
      }
      surv_j_mat <- surv_j_mat + (exp.eta.tw[[k]] * tmp) * p.byt *weight
    }
    V_i[, j, ] <- V_i[, j, ] + surv_j_mat
  }
  
  dL_i_list <- vector("list", n)
  for (ii in seq_len(n)) {
    V_sub <- V_i[ii, , , drop = FALSE]   # 1 × ncz × M
    V_sub <- matrix(V_sub, ncz, M)
    dL_i_list[[ii]] <- sqrt(2) * (V_sub %*% tao)
  }
  
  score.b_i <- vector("list", n)
  for (ii in seq_len(n)) {
    dL_i <- dL_i_list[[ii]]
    if (!diag.D) {
      score.b_i[[ii]] <- 0.5 * diag(dL_i) * diag(L)
    } else {
      idx <- upper.tri(L, diag = TRUE)
      dL_vec_i <- dL_i[idx]
      J <- jacobian2(t(L))
      score.b_i[[ii]] <- drop(dL_vec_i %*% J)
    }
  }
  score.b <- do.call(rbind, score.b_i)
  
  cbind(score.y, score.t, score.b)
}
