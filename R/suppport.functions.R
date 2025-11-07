gaussKronrod <- function (k = 15) {
  sk <- c(-0.949107912342758524526189684047851, -0.741531185599394439863864773280788, -0.405845151377397166906606412076961, 0,
          0.405845151377397166906606412076961, 0.741531185599394439863864773280788, 0.949107912342758524526189684047851, -0.991455371120812639206854697526329,
          -0.864864423359769072789712788640926, -0.586087235467691130294144838258730, -0.207784955007898467600689403773245, 0.207784955007898467600689403773245,
          0.586087235467691130294144838258730, 0.864864423359769072789712788640926, 0.991455371120812639206854697526329)
  wk15 <- c(0.063092092629978553290700663189204, 0.140653259715525918745189590510238, 0.190350578064785409913256402421014,
            0.209482141084727828012999174891714, 0.190350578064785409913256402421014, 0.140653259715525918745189590510238, 0.063092092629978553290700663189204,
            0.022935322010529224963732008058970, 0.104790010322250183839876322541518, 0.169004726639267902826583426598550, 0.204432940075298892414161999234649,
            0.204432940075298892414161999234649, 0.169004726639267902826583426598550, 0.104790010322250183839876322541518, 0.022935322010529224963732008058970)
  wk7 <- c(0.129484966168869693270611432679082, 0.279705391489276667901467771423780, 0.381830050505118944950369775488975,
           0.417959183673469387755102040816327, 0.381830050505118944950369775488975, 0.279705391489276667901467771423780, 0.129484966168869693270611432679082)
  if (k == 7)
    list(sk = sk[1:7], wk = wk7)
  else
    list(sk = sk, wk = wk15)
}


dropAttr <-
  function (mat) {
    d <- dim(mat)
    mat <- as.vector(mat)
    dim(mat) <- d
    mat
  }


#' @keywords internal
#' @noRd
chol.transf <-
  function (x) {
    if (any(is.na(x) | !is.finite(x)))
      stop("NA or infinite values in 'x'.\n")
    if (is.matrix(x)) {
      k <- nrow(x)
      U <- chol(x)
      U[cbind(1:k, 1:k)] <- log(U[cbind(1:k, 1:k)])
      U[upper.tri(U, TRUE)]
    } else {
      nx <- length(x)
      k <- round((- 1 + sqrt(1 + 8 * nx)) / 2)
      mat <- matrix(0, k, k)
      mat[upper.tri(mat, TRUE)] <- x
      mat[cbind(1:k, 1:k)] <- exp(mat[cbind(1:k, 1:k)])
      res <- crossprod(mat)
      attr(res, "L") <- t(mat)[lower.tri(mat, TRUE)]
      res
    }
  }

#' @keywords internal
#' @noRd
initial.surv <- function(Time, d, X2, WintF.vl, id, times, long=NULL, extra=NULL,  LongFormat){

  #---------basic check-----------
  stopifnot(length(Time) == length(d))
  if (!LongFormat) stopifnot(length(Time) == length(id), length(Time) == length(times))
  if (!is.null(long)) stopifnot(length(long) == length(times))

  #-------pre-processing the last long for each id----------
  if (!is.null(long)) {
    long.id <- tapply(long, id, tail, 1)       #  last value for each subject
  } else {
    long.id <- NULL
  }

  #-----------construct WW, the second covariate for Cox-----
  idT <- extra$idT
  strata <- extra$strata

  if (LongFormat) {
    WW <- if (is.null(long.id)) {
      as.matrix(X2)
    } else {
      cbind(as.matrix(X2), long_last = as.vector(long.id[idT]))
    }
  } else {
    WW <- if (is.null(long.id)) {
      as.matrix(X2[id, , drop = FALSE])
    } else {
      cbind(as.matrix(X2[id, , drop = FALSE]), long_last = as.vector(long.id))
    }
  }
  #----------construct the data DD for the first Cox---------
  if (!LongFormat) {
    DD <- data.frame(id = id, Time = Time[id], d = d[id], times = times, check.names = FALSE)
    k <- 0L
    if (!is.null(long)) {
      L <- as.matrix(long) * WintF.vl[id, , drop = FALSE]
      k <- ncol(L)
      colnames(L) <- paste0("long_", seq_len(k))
      DD <- cbind(DD, as.data.frame(L, check.names = FALSE))
    }

    # baseline covariates
    dW <- as.data.frame(X2[id, , drop = FALSE], check.names = FALSE)
    if (ncol(dW) > 0L) {
      colnames(dW) <- paste0("W", seq_len(ncol(dW)))
      DD <- cbind(DD, dW)
    }
    # Start-Stop & event
    DD$start <- DD$times
    # for each id：stop = c(times[-1], Time[1])
    ord <- order(DD$id, DD$start)
    DD <- DD[ord, ]
    idx <- with(DD, ave(seq_len(nrow(DD)), id, FUN = function(ii) {
      c(ii[-1], ii[1])
    }))
    DD$stop <- DD$Time[idx]

    # event
    DD$event <- with(DD, ave(d, id, FUN = function(x) {
      if (length(x) == 1L) x else { x[-length(x)] <- 0; x }
    }))
  } else {
    DD <- data.frame(Time = Time, d = d, check.names = FALSE)
    k <- 0L
    if (!is.null(long)) {
      L <- as.vector(long.id[idT]) * WintF.vl
      L <- as.matrix(L)
      k <- ncol(L)
      colnames(L) <- paste0("long_", seq_len(k))
      DD <- cbind(DD, as.data.frame(L, check.names = FALSE))
    }

    dW <- as.data.frame(X2, check.names = FALSE)
    if (ncol(dW) > 0L) {
      colnames(dW) <- paste0("W", seq_len(ncol(dW)))
      DD <- cbind(DD, dW)
    }
    DD$strata <- strata
  }

  #-------formulate for the first Cox--------
  lhs <- if (!LongFormat) "Surv(start, stop, event)" else "Surv(Time, d)"
  rhs_terms <- character(0)
  if (k > 0L) rhs_terms <- c(rhs_terms, paste0("long_", seq_len(k)))
  if (ncol(dW) > 0L) rhs_terms <- c(rhs_terms, paste0("W", seq_len(ncol(dW))))
  if (!is.null(DD$strata)) rhs_terms <- c(rhs_terms, "strata(strata)")

  if (length(rhs_terms) == 0L) rhs_terms <- "1"
  form <- as.formula(paste(lhs, "~", paste(rhs_terms, collapse = " + ")))
  cph <- survival::coxph(form, data = DD, x = FALSE)
  coefs <- stats::coef(cph)

  beta_idx <- grep("^long_", names(coefs))
  alpha_idx <- setdiff(seq_along(coefs), beta_idx)

  out <- list(
    beta  = if (length(beta_idx)) unname(coefs[beta_idx]) else numeric(0),
    alpha = if (length(alpha_idx)) unname(coefs[alpha_idx]) else numeric(0)
  )

  # ---- the second Cox: for piecewise-constanr----
  dat <- data.frame(Time = Time, d = d)
  cph2 <- survival::coxph(Surv(Time, d) ~ WW, data = dat, x = TRUE)

  init.fit <- piecewiseExp.ph(cph2, knots = extra$control$knots)
  cf2 <- init.fit$coefficients
  xi_idx <- grep("xi", names(cf2))
  out$xi <- exp(cf2[xi_idx])

  out
}

#' @keywords internal
#' @noRd
piecewiseExp.ph <-
  function (coxObject, knots = NULL, length.knots = 6) {
    Time <- coxObject$y[, 1]
    d <- coxObject$y[, 2]
    n <- length(Time)
    if (is.null(knots)) {
      Q <- length.knots + 1
      knots <- unique(quantile(Time, seq(0, 1, len = Q + 1), names = FALSE)[-c(1, Q + 1)])
      knots <- knots + 1e-06
      if (max(knots) > max(Time))
        knots[which.max(knots)] <- max(Time) - 1e-06
    }
    knots <- c(0, sort(knots), max(Time) + 1)
    Q <- length(knots) - 1
    ind <- findInterval(Time, knots, rightmost.closed = TRUE)
    D <- matrix(0, n, Q)
    D[cbind(seq_along(ind), ind)] <- 1
    D <- c(D * d)
    Tiq <- outer(Time, knots, pmin)
    T <- c(Tiq[, 2:(Q+1)] - Tiq[, 1:Q])
    X <- coxObject$x[rep(seq_len(n), Q), ]
    ND <- suppressWarnings(data.frame(Time = T, D = D, X, xi = gl(Q, n),
                                      check.names = FALSE)[T > 0, ])
    glm(D ~ . + offset(log(Time)) - Time - 1,
        family = poisson, data = ND)
  }


jacobian2 <- function(L) {
  idx <- lower.tri(L, diag = TRUE)
  vecl_L <- L[idx]
  p <- length(vecl_L)
  J <- diag(1, p)
  r <- row(L)[idx]; c <- col(L)[idx]
  diag_pos <- which(r == c)
  diag(J)[diag_pos] <- vecl_L[diag_pos]
  return(J)
}

