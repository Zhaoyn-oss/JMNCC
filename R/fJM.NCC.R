#' fJM-NCC: Likelihood inference with full maximum likelihood approach
#'
#' The fJM-NCC method incorporates all observed data and survival outcome, treating the unobserved
#' longitudinal outcomes as missing at random. The longitudinal sub-model must be a
#' generalized linear mixed-effects model (`lmerMod` or `glmerMod`) and the survival sub-model can
#' be a Cox model (`coxph`) .
#'
#' @param obsObject An `lmerMod` or `glmerMod` object (longitudinal sub-model).
#' @param survObject An object inheriting from class \code{coxph} or class \code{survreg}. In the call to \code{coxph()}
#' or \code{survreg()}, you need to specify the argument \code{x = TRUE} such that the design matrix is
#'  contained in the object fit. See \bold{Examples}..
#' @param timeVar Character string for the time variable in the longitudinal data.
#' @param deltaVar Character string for the event indicator column in `interFact$data`,
#' 0 representing control, k=1, 2... representing competing event k.
#' @param CompRisk Logical, whether to include competing risks (default `TRUE`).
#' @param interFact List with `value` (formula) and `data` (data frame) for additional covariates or strata.
#' @param control A list of control parameters with the following elements:
#' \describe{
#'   \item{iter.qN}{The number of quasi-Newton iterations. Default is 300.}
#'   \item{tol}{tolerance value for convergence in the log-likelihood.
#'    Default is \code{sqrt(.Machine$double.eps)}.}
#'   \item{parscale}{The \code{parscale} control argument for \code{optim()}.
#'   It should be a numeric vector of length equal to the number of parameters.
#'   Default is 0.01 for all parameters.}
#'   \item{knots}{A numeric vector of the knot positions for the piecewise-constant baseline risk function.}
#'   \item{ObsTimes.knots}{Logical; if \code{TRUE} (default), the positions of the knots are specified
#'   based on the observed event times; otherwise, they are based only on the true event times.}
#'   \item{Q}{The number of internal knots, 4 default.}
#'   \item{GHk}{The number of Gauss–Hermite quadrature points used to approximate the integrals over the random effects.
#'   The default is 15 for one-, two-, or three-dimensional integration when \eqn{N < 2000}, and 9 otherwise,
#'   for the simple Gauss–Hermite rule.}
#'   \item{GKk}{The number of Gauss–Kronrod points used to approximate the integral involved in the calculation
#'   of the survival function. Two options are available: 7 (default) or 15.}
#' }
#'
#' @return A list of class `"JM.NCC"` with estimated coefficients, variance components, Hessian matrix, and other elements.
#'
#' @details
#' Constructs design matrices for longitudinal and survival sub-models, sets up
#' piecewise-constant baseline hazards, and performs joint likelihood estimation
#' using Gauss-Hermite and Gauss-Kronrod quadrature.
#'
#' @importFrom stats coef model.frame model.matrix model.response na.pass terms as.formula
#' @importFrom utils getFromNamespace
#' @importFrom lme4 fixef findbars getME ranef
#' @importFrom survival coxph survreg Surv
#'
#' @examples
#' \dontrun{
#' # load example data for linear mixed-effects model
#' data(obs_full_linear)
#' data(sur_full_linear)
#'
#' # Convert survival data to long format for competing risks
#' full.surdata.long <- crLong(data = sur_full_linear, statusVar="event",
#' censLevel = 1, nameStrata = "CR", nameStatus = "comp.event")
#'
#' # Fit Cox model
#' survObject.full <- coxph(Surv(sur.time, comp.event) ~ gender+strata(CR),
#' data = full.surdata.long, x = TRUE)
#'
#' # Fit longitudinal model
#' obsObject.full <- lmer(obs.value~sample.time+(1+sample.time|subject),
#' data=obs_full_linear)
#'
#' # Fit joint model
#' res.fjm <- fJM.NCC(obsObject=obsObject.full, survObject=survObject.full,
#' timeVar = "sample.time", deltaVar="delta", CompRisk=TRUE,
#' interFact = list(value = ~ CR, data = full.surdata.long), control = list(Q=4))
#'
#'
#' # load example data for Poisson generalized linear mixed-effects model
#' data(obs_full_poisson)
#' data(sur_full_poisson)
#'
#' obsObject.full <- glmer(obs.value~sample.time+(1+sample.time|subject),
#' data=obs_full_poisson, family=poisson(link = "log"))
#' full.surdata.long <- crLong(data = sur_full_poisson, statusVar="event",
#' censLevel = 1, nameStrata = "CR", nameStatus = "comp.event")
#' survObject.full <- coxph(Surv(sur.time, comp.event) ~ gender+strata(CR),
#' data = full.surdata.long, x = TRUE)
#' res.fjm <- fJM.NCC(obsObject=obsObject.full, survObject=survObject.full,
#' timeVar = "sample.time", deltaVar="delta", CompRisk=TRUE,
#' interFact = list(value = ~ CR, data = full.surdata.long), control = list(Q=4))
#'
#' }
#' @export
fJM.NCC <- function(obsObject, survObject, timeVar, deltaVar,
                    CompRisk=TRUE, interFact, control=list()){

  #--------------------------------------------------
  cl <- match.call()
  # Check model class and link/family
  if (inherits(obsObject, "glmerMod")) {
    fam <- family(obsObject)$family
    linkfun <- family(obsObject)$link
    if (!(fam == "poisson" && linkfun == "log")) {
      stop("\nInvalid GLMM specification.
             For glmer models, use: family = poisson(link = 'log').")
    }
  } else if (inherits(obsObject, "lmerMod")) {
    # lmer always uses gaussian(identity)
    # but we check the link in case user tried something exotic
    fam <- "gaussian"
    linkfun <- "identity"
    if (!(fam == "gaussian" && linkfun == "identity")) {
      stop("\nInvalid LMM specification.
             For lmer models, Gaussian identity link is required.")
    }
  } else {
    stop("\nobsObject must be either an lmerMod or glmerMod model.")
  }

  if (length(getME(obsObject, "flist")) > 1)
    stop("\nMultiple or nested random effects detected.\nThis function only supports a single random grouping factor.")

  if (!inherits(survObject, "coxph") && !inherits(survObject, "survreg"))
    stop("\n'survObject' must inherit from class coxph or class survreg.")

  if (!is.matrix(survObject$x))
    stop("\nuse argument 'x = TRUE' in ",
         if (inherits(survObject, "coxph")) "'coxph()'." else "'survreg()'.")

  if (length(timeVar) != 1 || !is.character(timeVar))
    stop("\n'timeVar' must be a character string.")

  if (length(deltaVar) != 1 || !is.character(deltaVar))
    stop("\n'deltaVar' must be a character string.")

  if (!deltaVar %in% names(interFact$data))
    stop("\n'deltaVar' does not correspond to one of the columns in the 'interFact:data'.")

  #----------------survival sub-model-------------------
  formT <- formula(survObject)
  if (inherits(survObject, "coxph")) {
    X2 <- survObject$x                     # baseline covariates associated with competing event in survival sub-model (alpha)
    keepX2 <- suppressWarnings(!is.na(survObject$coefficients))  #
    X2 <- X2[, keepX2, drop = FALSE]       # Select effective regression coefficients (alpha)
    if (CompRisk) {
      nRisks <- length(unique(survObject$strata))
    } else {
      nRisks <- 1
    }
    surv <- survObject$y                   # survival information (event time and status)
    if (attr(surv, "type") == "right") {
      LongFormat <- FALSE
      Time <- survObject$y[, 1]             #  event time for each subject
      d <- survObject$y[, 2]                # event indicator: 1 experience any type competing risk, 0 control
    } else if (attr(surv, "type") == "counting") {
      LongFormat <- TRUE
      if (is.null(survObject$model))
        stop("\nplease refit the Cox model including in the ",
             "call to coxph() the argument 'model = TRUE'.")
      Time <- survObject$y[, 2]
      d <- survObject$y[, 3]
    }
    idT <- if (!is.null(survObject$model$`(cluster)`)) {
      as.vector(unclass(survObject$model$`(cluster)`))
    } else {
      if (!CompRisk) seq_along(Time)
      else rep(seq_len(length(Time)/nRisks), each = nRisks)
    }
    idT <- match(idT, unique(idT))        # index for each subject 1, 1, 2, 2, 3, 3, 4, 4,....
  } else {
    X2 <- survObject$x[, -1, drop = FALSE]
    Time <- exp(survObject$y[, 1])
    d <- survObject$y[, 2]
    idT <- seq_along(Time)
    LongFormat <- FALSE
    nRisks <- 1
  }

  nT <- length(unique(idT))   # subject size in survival sub-model
  if (LongFormat && is.null(survObject$model$`(cluster)`))
    stop("\nuse argument 'model = TRUE' and cluster() in coxph().")
  if (!length(X2)) X2 <- NULL
  if (sum(d) < 5) warning("\nmore than 5 events are required.")

  WintF.vl <- as.matrix(rep(1, length(Time)))            # design matrix for competing event (strata)
  if (!is.null(interFact)) {
    if (!is.null(interFact$value)){
      WintF.vl <- if (is.null(survObject$model)) {
        model.matrix(interFact$value, data = interFact$data)
      } else {
        model.matrix(interFact$value, data = survObject$model)
      }
    }
  }

  delta <- interFact$data[[deltaVar]]
  delta <- delta[!duplicated(idT)]
  delta <- as.numeric(delta)     # competing event k=0, 1, 2, 3....

  #----------------longitudinal sub-model--------------
  data <- eval(getCall(obsObject)$data, envir = environment(obsObject))  # observation for all subjects, including missing
  randomVar <- names(getME(obsObject, "flist"))                          # random effcts variates
  group.factor <- factor(data[[randomVar]])                              # group ID
  id <- match(group.factor, unique(group.factor))                        # subject id: 1, 2, 3, 4, 4, 4, 4.....
  nY <- length(unique(id))                                               # subject size in longitudinal sub-model

  # fixed effect part
  TermsX <- terms(obsObject)
  formYx <- lme4::nobars(formula(obsObject))                              # fixed coefficient formula
  mfX <- model.frame(TermsX, data=data, na.action=na.pass)                # fixed effects design matrix
  X <- model.matrix(formYx, mfX, na.action=na.pass)
  y.long <- model.response(mfX, "numeric")

  # random effect part
  reTerm <- lme4::findbars(formula(obsObject))[[1]]
  formYz <- as.formula(paste("~", deparse(reTerm[[2]])))
  mfZ <- model.frame(terms(formYz), data = data, na.action=na.pass)
  TermsZ <- attr(mfZ, "terms")
  Z <- model.matrix(formYz, mfZ, na.action = na.pass)  # random effect design matrxi

  b_list <- ranef(obsObject)  # random effects
  b <- do.call(rbind, lapply(names(b_list), function(grp) b_list[[grp]]))  # random effects for each subject
  fitted.subject.lel <- rownames(b)                                    # subjects with observation
  subject.unique <- group.factor[!duplicated(group.factor)]            # full cohort subject
  match_index <- match(subject.unique, fitted.subject.lel)
  b <- b[match_index, , drop = FALSE]
  b <- data.matrix(b)                                                 # random effects for full cohort subjects
  dimnames(b) <- NULL
  long <- c(X %*% fixef(obsObject)) + rowSums(Z * b[id, ])            # longitudinal sub-model: X(1)*r+Z*b

  if (!timeVar %in% names(data))
    stop("\n'timeVar' does not correspond to one of the columns in the model.frame of 'obsObject'.")

  #----------------data in longitudinal and survival sub-model--------------

  obsdata <- list(y=y.long, X=X, Z=Z)
  surdata <- list(X2=X2[!duplicated(idT), ,drop=FALSE], delta=delta, nRisks=nRisks)

  con <- list(iter.qN = 300, tol = if (!CompRisk) sqrt(.Machine$double.eps) else 1e-09,
              parscale = NULL, knots = NULL, ObsTimes.knots = TRUE,
              Q = 3, GHk = if (ncol(Z) < 3 && nrow(Z) < 2000) 15 else 9,
              GKk = 7, verbose = FALSE)
  con[(namc <- names(control))] <- control
  con$family <- fam

  #----------------split of time scale qs: 0=v0<v1<v2<..<vQ-----------
  event.time <- Time[!duplicated(idT)]                                       # event time for subjects
  if (is.null(con$knots) || !is.numeric(con$knots)) {
    Q <- con$Q                                        #the number of time intervals for the piecewise baseline hazard
    qs <- if (con$ObsTimes.knots) {
      unique(quantile(unique(event.time), seq(0, 1, len = Q + 1), names = FALSE)[-c(1, Q + 1)])
    } else {
      unique(quantile(unique(Time[d == 1]), seq(0, 1, len = Q - 1), names = FALSE))
    }
    qs <- qs + 1e-06
    if (max(qs) > max(event.time)) qs[which.max(qs)] <- max(event.time) - 1e-06
    con$knots <- qs
    qs <- c(0, qs, max(event.time))  #
    Q <- length(qs) - 1
  } else {
    qs <- c(0, sort(con$knots) + 1e-06, max(event.time))
    Q <- length(qs) - 1
  }

  ind <- findInterval(event.time, qs, rightmost.closed = TRUE)
  Tiq <- outer(event.time, qs, pmin)
  Lo <- Tiq[, 1:Q]
  Up <- Tiq[, 2:(Q+1)]
  P <- (Up - Lo) / 2
  P[P < con$tol] <- as.numeric(NA)
  P1 <- (Up + Lo) / 2

  #-----------------------GK notes and weights----------------------
  wk <- gaussKronrod(con$GKk)$wk
  sk <- gaussKronrod(con$GKk)$sk
  nk <- con$GKk

  st <- matrix(0, nrow(P), nk*Q)                       # All the integral point positions of each individual
  skQ <- rep(sk, Q)
  for (i in seq_along(event.time)) {
    st[i, ] <- rep(P[i, ], each = nk) * skQ + rep(P1[i, ], each = nk)
  }

  #-----------for E(Yi|Ti)---------------
  data.Ti <- data[!duplicated(id), ]
  data.Ti[[timeVar]] <- pmax(event.time, 0)
  mfX.Ti <- model.frame(TermsX, data = data.Ti, na.action=na.pass)
  mfZ.Ti <- model.frame(TermsZ, data = data.Ti, na.action=na.pass)
  Xtime <- model.matrix(formYx, mfX.Ti)                  # for E(Yi|Ti)
  Ztime <- model.matrix(formYz, mfX.Ti)
  surdata <- c(surdata, list(Xtime=Xtime, Ztime=Ztime))

  #-----------for E(Yi|st)-----------------
  id.GK <- rep(seq_along(event.time), rowSums(!is.na(st)))
  P <- c(t(P))
  data.st <- data.Ti[rep(seq_along(event.time), each = nk*Q), ]
  data.st[[timeVar]] <- pmax(c(t(st)), 0)
  data.st <- data.st[!is.na(data.st[[timeVar]]), ]
  surdata <- c(surdata, list(ind.D = ind))

  mfX <- model.frame(TermsX, data = data.st, na.action=na.pass)
  mfZ <- model.frame(TermsZ, data = data.st, na.action=na.pass)
  Xs <- model.matrix(formYx, mfX)
  Zs <- model.matrix(formYz, mfZ)
  obsdata <- c(obsdata, list(P = P[!is.na(P)], st = st[!is.na(st)], wk = wk, id.GK = id.GK, Q = Q, Xs = Xs, Zs = Zs))

  #----------------initial values in longitudinal and survival sub-model--------------
  # covariance matrix for random effects
  VC_list <- VarCorr(obsObject)
  VC1 <- VC_list[[1]]
  sd <- attr(VC1, "stddev")
  cor <- attr(VC1, "correlation")
  VC <- diag(sd, nrow = length(sd), ncol = length(sd)) %*% cor %*% diag(sd, nrow = length(sd), ncol = length(sd))
  if (all(VC[upper.tri(VC)] == 0)) VC <- diag(VC)

  #initial value for longitudinal sub-model
  gammas <- fixef(obsObject)
  if(fam=="gaussian"){
    sigma <- sigma(obsObject)
  }
  if(fam == "poisson"){
    sigma <- NULL
  }
  D <- VC

  # initial value for survival sub-model
  init.surv <- initial.surv(Time=Time, d=d, X2=X2, WintF.vl=WintF.vl, id=id,
                            times = data[[timeVar]], long = long,
                            extra = list(idT = idT, strata = survObject$strata, control=con),
                            LongFormat = CompRisk | length(Time) > nT)
  betas <- init.surv$beta
  betas <- c(betas[1], betas[-1]+betas[1])
  alphas <- replicate(nRisks, init.surv$alpha, simplify = FALSE)
  xi <- replicate(nRisks, init.surv$xi, simplify = FALSE)
  initial.values <- list(gammas=gammas, sigma=sigma, D=D, betas=betas, alphas=alphas, xi=xi)

  #--------------------------results-------------------------------
  rmObjs <- c("X","Z","P", "st","wk", "id.GK","Q","Xs","Zs", "y.long", "mfX", "mfZ", "data.st",
              "delta","nRisks", "WintF.vl", "Xtime" , "Ztime")
  rm(list = rmObjs); gc()
  control=con

  out <- piecewisefJM.fit(obsdata=obsdata, surdata = surdata, id=id,
                           initial.values=initial.values, control=con)

  # check for problems with the Hessian at convergence
  H <- out$Hessian
  if (any(is.na(H) | !is.finite(H))) {
    warning("infinite or missing values in Hessian at convergence.\n")
  } else {
    eig <- eigen(H, symmetric = TRUE)
    if (!all(eig$values >= 0)){
      warning("Hessian matrix at convergence is not positive definite.\n")
    }
  }
  out$coefficients <- out$coefficients[!sapply(out$coefficients, is.null)]
  class(out) <- "JM.NCC"

  H.inverse.try <- try(solve(H), silent = TRUE)
  vmat <- if (!inherits(H.inverse.try, "try-error")){
    structure(H.inverse.try, dimnames = dimnames(H))
  }else{
    structure(ginv(H), dimnames = dimnames(H))
  }
  H.inverse <- (vmat + t(vmat)) / 2
  sd.hessian <- sqrt(diag(H.inverse))
  score.mat <- out$score.mat
  score.mat[is.na(score.mat)] <- 0
  S <- crossprod(score.mat)
  Cov.sand <- H.inverse %*% S %*% H.inverse
  Cov.sand <- (Cov.sand +t(Cov.sand))/2
  sd.sand <- sqrt(diag(Cov.sand))

  if(fam=="gaussian"){
    gammas <- out$coefficients$gammas; p_gamma <- length(gammas)
    log.sigma <- log(out$sigma); p_sigma <- length(log.sigma)
    betas <- out$coefficients$betas; p_beta <- length(betas)
    alphas <- unlist(out$coefficients$alphas); p_alpha <- length(alphas)
    log.xi <- unlist(out$coefficients$xi); p_xi <- length(log.xi)
    D <- out$D

    idx_gammas <- 1:p_gamma
    idx_sigma  <- (max(idx_gammas) + 1):(max(idx_gammas) + p_sigma)
    idx_betas  <- (max(idx_sigma)  + 1):(max(idx_sigma)  + p_beta)
    idx_alphas <- (max(idx_betas)  + 1):(max(idx_betas)  + p_alpha)
    idx_xi     <- (max(idx_alphas) + 1):(max(idx_alphas) + p_xi)
  }
  if(fam=="poisson"){
    gammas <- out$coefficients$gammas; p_gamma <- length(gammas)
    betas <- out$coefficients$betas; p_beta <- length(betas)
    alphas <- unlist(out$coefficients$alphas); p_alpha <- length(alphas)
    log.xi <- unlist(out$coefficients$xi); p_xi <- length(log.xi)
    D <- out$D

    idx_gammas <- 1:p_gamma
    idx_betas  <- (max(idx_gammas) + 1):(max(idx_gammas) + p_beta)
    idx_alphas <- (max(idx_betas)  + 1):(max(idx_betas)  + p_alpha)
    idx_xi     <- (max(idx_alphas) + 1):(max(idx_alphas) + p_xi)
  }

  coeff <- c(gammas, betas, alphas, log.xi)
  idx.coeff <- c(idx_gammas, idx_betas, idx_alphas, idx_xi)

  sd.hessian <- sd.sand[idx.coeff]
  sd.sand <- sd.sand[idx.coeff]
  Zval.hessian <- abs(coeff/sd.hessian)
  Zval.sand <- abs(coeff/sd.sand)
  pval.hessian <- 2*(1-pnorm(Zval.hessian))
  pval.sand <- 2*(1-pnorm(Zval.sand))

  coefficients <- cbind(Est=coeff, Se.sand=sd.sand, Z.sand=Zval.sand, P.sand=pval.sand)
  out$coefficients <- coefficients
  out
}
