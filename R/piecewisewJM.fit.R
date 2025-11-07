piecewisewJM.fit <- function(obsdata, surdata, id, initial.values, control){
  
  #response vectors
  y <- as.vector(obsdata$y)
  obs.id <- which(!is.na(y))
  ind.D <- surdata$ind.D 
  
  # design matrix
  X <- obsdata$X
  Z <- obsdata$Z 
  WW <- surdata$X2
  delta <- surdata$delta
  weight <- surdata$weight
  Xtime <- surdata$Xtime
  Ztime <- surdata$Ztime
  Xs <- obsdata$Xs
  Zs <- obsdata$Zs 
  
  X <- dropAttr(X); Z <- dropAttr(Z)
  Xtime <- dropAttr(Xtime); Ztime <- dropAttr(Ztime);
  Xs <- dropAttr(Xs); Zs <- dropAttr(Zs)
  
  # sample size
  ncx <- ncol(X)             
  ncz <- ncol(Z)
  ni <- as.vector(tapply(id, id, length))       # the number of records for each subject
  n <- length(ni)                               # subject sizes
  nRisks <- surdata$nRisks                      # competing event number
  CompRisk <- nRisks>1
  Q <- obsdata$Q                                
  
  ZtZ <- lapply(split(seq_along(y), id),function(idx){
    obs_idx <- idx[!is.na(y[idx])]
    if(length(obs_idx)==0){
      Zi <- matrix(0, ncz, ncz)
    }else{
      Zi <- matrix(Z[obs_idx, ], ncol=ncz)
      Zi <- crossprod(Zi)
    }
    return(Zi)  
  }) 
  names(ZtZ) <- NULL
  ZtZ <- matrix(unlist(ZtZ), n, ncz * ncz, TRUE)  # nsubject x 4  
  sub.id <- lapply(split(seq_along(y), id),function(idx){
    obs_idx <- idx[!is.na(y[idx])]
    if(length(obs_idx)==0) sub.id <- unique(id[idx])
  })
  sub.id <- unlist(sub.id)
  
  fam <- control$family
  # Gauss-Hermite quadrature
  GH <- statmod::gauss.quad(control$GHk, kind = "hermite")   
  tao <- as.matrix(expand.grid(rep(list(GH$nodes), ncz))) 
  M <- nrow(tao) 
  wGH <- as.matrix(expand.grid(rep(list(GH$weights), ncz))) 
  wGH <- pi^(-ncz/2)*apply(wGH, 1, prod)
  
  #bm <- sqrt(2) * t(control$inv.chol.VC %*% t(tao)) 
  #dimnames(bm) <- NULL 
  #M <- nrow(tao)                                # the number of notes 
  
  # Gauss-Kronrod quadrature 
  nk <- control$GKk
  id.GK <- obsdata$id.GK    
  P <- as.vector(obsdata$P)  
  st <- obsdata$st 
  ind.K <- rep(unlist(lapply(ind.D, seq_len)), each = nk)
  wk <- unlist(lapply(ind.D, function (n) rep(obsdata$wk, n)))
  wkP <- wk * rep(P, each = nk)
  
  # initial values
  gammas <- as.vector(initial.values$gammas)
  sigma <- initial.values$sigma
  D <- initial.values$D
  betas <- initial.values$betas
  alphas <- if (!is.null(WW)) initial.values$alphas else NULL    
  xi <- initial.values$xi
  
  diag.D <- is.matrix(D)
  if (diag.D) dimnames(D) <- NULL else names(D) <- NULL  
  if(!CompRisk) log.xi <- log(xi) else log.xi <- lapply(xi, log)
  old <- options(warn = (-1))
  on.exit(options(old))
  if(fam=="gaussian"){
    list.thetas <- list(gammas = gammas, log.sigma = log(sigma), 
                        betas = betas, alphas=alphas, log.xi = log.xi, 
                        D = if (!diag.D) log(D) else chol.transf(D)) 
  }
  if(fam=="poisson"){
    list.thetas <- list(gammas = gammas, betas = betas, alphas=alphas, log.xi = log.xi, 
                        D = if (!diag.D) log(D) else chol.transf(D)) 
  }
  list.thetas <- list.thetas[!sapply(list.thetas, is.null)]
  thetas <- unlist(as.relistable(list.thetas)) 
  
  # fix environments for functions
  environment(LogLikGaussian.piecewisewJM) <- environment(ScoreGaussian.piecewisewJM) <- environment(ScoreGaussianMat.piecewisewJM) <- environment()
  environment(LogLikPoisson.piecewisewJM) <- environment(ScorePoisson.piecewisewJM) <-  environment(ScorePoissonMat.piecewisewJM) <-environment()
 
  #--------------quasi-Newton iterations---------------------------
  if (is.null(control$parscale)){
    parscale <- rep(0.01, length(thetas))
  }else if(length(control$parscale)==1){
    parscale <- rep(control$parscale, length(thetas))
  }else{
    parscale <- control$parscale
  }
  if (control$verbose)  cat("\n\nquasi-Newton iterations start.\n\n")
  
  if(fam=="gaussian"){
    out <- optim(thetas, LogLikGaussian.piecewisewJM, ScoreGaussian.piecewisewJM, method = "BFGS", hessian = TRUE,
                 control = list(maxit = control$iter.qN, parscale = parscale, trace = 10 * control$verbose))
    thetas <- relist(out$par, skeleton = list.thetas)
    gammas <- thetas$gammas
    sigma <- thetas$log.sigma
    sigma <- exp(sigma)
    D <- thetas$D
    D <- if (!diag.D) exp(D) else chol.transf(D)
    betas <- thetas$betas
    alphas <- thetas$alphas
    xi <- thetas$log.xi
    conv <- out$convergence
    score <- ScoreGaussian.piecewisewJM(unlist(thetas))
    score.mat <- ScoreGaussianMat.piecewisewJM(unlist(thetas)) 
  } 
  if(fam=="poisson"){
    out <- optim(thetas, LogLikPoisson.piecewisewJM, ScorePoisson.piecewisewJM, method = "BFGS", hessian = TRUE,
                 control = list(maxit = control$iter.qN, parscale = parscale, trace = 10 * control$verbose))
    
    thetas <- relist(out$par, skeleton = list.thetas)
    gammas <- thetas$gammas 
    D <- thetas$D
    D <- if (!diag.D) exp(D) else chol.transf(D)
    betas <- thetas$betas
    alphas <- thetas$alphas
    xi <- thetas$log.xi
    conv <- out$convergence
    score <- ScorePoisson.piecewisewJM(unlist(thetas))
    score.mat <- ScorePoissonMat.piecewisewJM(unlist(thetas)) 
  }  
  Hessian <- out$hessian
  lgLik <- out$value
  names(gammas) <- names(initial.values$gammas)
  if (diag.D) dimnames(D) <- dimnames(initial.values$D) else names(D) <- names(initial.values$D)
  names(betas) <- paste0("Event", 1:nRisks) 
  alphas <- unlist(alphas)
  names(alphas) <- if(!CompRisk) colnames(surdata$X2) else{unlist(lapply(1:nRisks, function(i) paste0(colnames(surdata$X2), "(", i, ")")))}
  xi <- unlist(xi)
  names(xi) <- if (!CompRisk) paste("xi", 1:Q, sep = "") else { 
    paste("xi", sapply(Q, seq_len), "(", rep(1:nRisks, each=Q), ")", sep = "")
  }
  if(fam=="gaussian"){
    nams <- c(paste("Y.", c(names(gammas), "sigma"), sep = ""), 
              paste("T.", c(names(betas), names(alphas), names(xi)), sep = ""),
              paste("B.", paste("D", seq(1, ncz * (ncz + 1) / 2), sep = ""), sep = "")) 
  }
  if(fam=="poisson"){
    nams <- c(paste("Y.", c(names(gammas)), sep = ""), 
              paste("T.", c(names(betas), names(alphas), names(xi)), sep = ""),
              paste("B.", paste("D", seq(1, ncz * (ncz + 1) / 2), sep = ""), sep = ""))
  }
  dimnames(Hessian) <- list(nams, nams)
  colnames(score.mat) <- nams
  names(score) <- nams
  
  coefficients <- list(gammas = gammas, betas = betas, alphas=alphas, xi = xi) 
  res <- list(coefficients =coefficients, sigma = sigma, D = as.matrix(D), score=score, 
              Hessian = Hessian, score.mat=score.mat, logLik = lgLik, 
              convergence = conv)
}