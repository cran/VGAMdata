# Copyright (C) 1998-2020. T. W. Yee    All rights reserved.
# 20201207: family.vd2.R is for \pkg{VGAMdata}.
#   This file contains any 'old' functions from \pkg{VGAM} such
#   as ones made obsolete by new GAIT family functions.



# Last modified:
# 20201207: Trying to put [dpqr]tikuv(), etc. here.
#   Thats because family.vd1.R is getting big.






# ==========================================================


# 20060525 [dp]tikuv look fine.


dtikuv <- function(x, d, mean = 0, sigma = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  if (!is.Numeric(d, length.arg = 1) ||
      max(d) >= 2)
    stop("bad input for argument 'd'")

  L <- max(length(x), length(mean), length(sigma))
  if (length(x)     != L) x     <- rep_len(x,     L)
  if (length(mean)  != L) mean  <- rep_len(mean,  L)
  if (length(sigma) != L) sigma <- rep_len(sigma, L)


  hh <- 2 - d
  KK <- 1 / (1 + 1/hh + 0.75/hh^2)
  logden <- dnorm(x = x, mean = mean, sd = sigma, log = TRUE) +
    log(KK) + 2 * log1p(((x-mean)/sigma)^2 / (2*hh))
  logden[is.infinite(x)] <- log(0)  # 20141209 KaiH
  if (log.arg) logden else exp(logden)
}  # dtikuv



ptikuv <- function(q, d, mean = 0, sigma = 1,
                   lower.tail = TRUE, log.p = FALSE) {
  if (!is.Numeric(d, length.arg = 1) ||
      max(d) >= 2)
    stop("bad input for argument 'd'")

  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.arg <- log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")
  rm(log.p)  # 20141231 KaiH

  L <- max(length(q), length(mean), length(sigma))
  if (length(q)     != L) q     <- rep_len(q,     L)
  if (length(mean)  != L) mean  <- rep_len(mean,  L)
  if (length(sigma) != L) sigma <- rep_len(sigma, L)

  zedd1 <- 0.5 * ((q - mean) / sigma)^2
  ans <- q*0 + 0.5
  hh <- 2 - d
  KK <- 1 / (1 + 1/hh + 0.75/hh^2)
  if (any(lhs <- q < mean)) {
    ans[lhs] <- ( KK/(2*sqrt(pi))) * (
    gamma(0.5) * (1 - pgamma(zedd1[lhs], 0.5)) +
    2 * gamma(1.5) * (1 - pgamma(zedd1[lhs], 1.5)) / hh +
    gamma(2.5) * (1 - pgamma(zedd1[lhs], 2.5)) / hh^2)
  }
  if (any(rhs <- q > mean)) {
    ans[rhs] <- 1.0 - Recall(q = (2*mean[rhs] - q[rhs]), d = d,
                             mean = mean[rhs], sigma = sigma[rhs])
  }

# 20141231 KaiH
  if (lower.tail) {
    if (log.arg) log(ans) else ans
  } else {
    if (log.arg) log1p(-ans) else 1 - ans
  }
}  # ptikuv




qtikuv <- function(p, d, mean = 0, sigma = 1,
                   lower.tail = TRUE, log.p = FALSE, ...) {
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")
# rm(log.p)   # 20150102 KaiH

# if (!is.Numeric(p, positive = TRUE) || max(p) >= 1)
#   stop("bad input for argument 'p'")
  if (!is.Numeric(d, length.arg = 1) || max(d) >= 2)
    stop("bad input for argument 'd'")
# if (!is.Numeric(mean))
#   stop("bad input for argument 'mean'")
# if (!is.Numeric(sigma))
#   stop("bad input for argument 'sigma'")

### 20150102 KaiH
  orig.p <- p
  if (lower.tail) {
    if (log.p) p <- exp(p)
  } else {
    p <- if (log.p) -expm1(p) else 1 - p
  }

  L <- max(length(p), length(mean), length(sigma))
  if (length(p)     != L) p     <- rep_len(p,     L)
  if (length(mean)  != L) mean  <- rep_len(mean,  L)
  if (length(sigma) != L) sigma <- rep_len(sigma, L)
  ans <- rep_len(0.0, L)


  myfun <- function(x, d, mean = 0, sigma = 1, p)
    ptikuv(q = x, d = d, mean = mean, sigma = sigma) - p

  for (ii in 1:L) {
    Lower <- ifelse(p[ii] <= 0.5, mean[ii] - 3 * sigma[ii], mean[ii])
    while (ptikuv(q = Lower, d = d, mean = mean[ii],
                  sigma = sigma[ii]) > p[ii])
      Lower <- Lower - sigma[ii]
    Upper <- ifelse(p[ii] >= 0.5, mean[ii] + 3 * sigma[ii], mean[ii])
    while (ptikuv(q = Upper, d = d, mean = mean[ii],
                  sigma = sigma[ii]) < p[ii])
      Upper <- Upper + sigma[ii]
#print("c(Lower,Upper)")
#print( c(Lower,Upper) )
    ans[ii] <- uniroot(f = myfun, lower = Lower, upper = Upper,
                       d = d, p = p[ii],
                       mean = mean[ii], sigma = sigma[ii], ...)$root
  }


  if (log.p) {
    ans[orig.p > 0] <- NaN
  } else {
    ans[orig.p < 0] <- NaN
    ans[orig.p > 1] <- NaN
  }

  ans
}  # qtikuv



# 20060526; Uses the rejection method.
# Note: mean and sigma (and d) must be of length 1.
rtikuv <- function(n, d, mean = 0, sigma = 1, Smallno = 1.0e-6) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
            stop("bad input for argument 'n'") else n
# if (!is.Numeric(n, positive = TRUE, integer.valued = TRUE))
#   stop("bad input for argument 'n'")
  if (!is.Numeric(d, length.arg = 1) || max(d) >= 2)
    stop("bad input for argument 'd'")
  if (!is.Numeric(mean, length.arg = 1))
    stop("bad input for argument 'mean'")
  if (!is.Numeric(sigma, length.arg = 1))
    stop("bad input for argument 'sigma'")
  if (!is.Numeric(Smallno, positive = TRUE, length.arg = 1) ||
      Smallno > 0.01 ||
      Smallno < 2 * .Machine$double.eps)
      stop("bad input for argument 'Smallno'")
  ans <- rep_len(0.0, use.n)

  ptr1 <- 1; ptr2 <- 0
  hh <- 2 - d
  KK <- 1 / (1 + 1/hh + 0.75/hh^2)
# if (4 - 2*hh > 0) cat("bimodal\n") else cat("unimodal\n")
  ymax <- ifelse(hh < 2,
                 dtikuv(x = mean + sigma*sqrt(4 - 2*hh),
                        d = d, mean = mean, sigma = sigma),
                 KK / (sqrt(2 * pi) * sigma))
  while (ptr2 < use.n) {
    Lower <- mean - 5 * sigma
    while (ptikuv(q = Lower, d = d, mean = mean, sigma = sigma) > Smallno)
      Lower <- Lower - sigma
    Upper <- mean + 5 * sigma
    while (ptikuv(q = Upper, d = d,
                  mean = mean, sigma = sigma) < 1-Smallno)
      Upper <- Upper + sigma
#print("c(Lower,Upper)")
#print( c(Lower,Upper) )
    x <- runif(2*use.n, min = Lower, max = Upper)
    index <- runif(2*use.n, max = ymax) <
             dtikuv(x, d = d, mean = mean, sigma = sigma)
    sindex <- sum(index)
    if (sindex) {
      ptr2 <- min(use.n, ptr1 + sindex - 1)
      ans[ptr1:ptr2] <- (x[index])[1:(1+ptr2-ptr1)]
      ptr1 <- ptr2 + 1
    }
  }
  ans
}  # rtikuv




# 20060524
# Reference: TEST Akkaya and Tiku, in press.
# Works.
 tikuv <- function(d, lmean = "identitylink", lsigma = "loglink",
#                  emean = list(), esigma = list(),
                   isigma = NULL, zero = "sigma") {
# if (mode(lmean) != "character" && mode(lmean) != "name")
#   lmean = as.character(substitute(lmean))
# if (mode(lsigma) != "character" && mode(lsigma) != "name")
#   lsigma = as.character(substitute(lsigma))
# if (!is.list(emean)) emean = list()
# if (!is.list(esigma)) esigma = list()


  lmean <- as.list(substitute(lmean))
  emean <- link2list(lmean)
  lmean <- attr(emean, "function.name")

  lsigma <- as.list(substitute(lsigma))
  esigma <- link2list(lsigma)
  lsigma <- attr(esigma, "function.name")



# if (length(zero) &&
#    (!is.Numeric(zero, integer.valued = TRUE, positive = TRUE) ||
#    max(zero) > 2))
#   stop("bad input for argument 'zero'")
  if (!is.Numeric(d, length.arg = 1) || max(d) >= 2)
      stop("bad input for argument 'd'")



  new("vglmff",
  blurb = c("Short-tailed symmetric [Tiku and Vaughan (1999)] ",
            "distribution\n",
          "Link:     ",
          namesof("mean",  lmean,  earg = emean), ", ",
          namesof("sigma", lsigma, earg = esigma),
          "\n", "\n",
          "Mean:     mean"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
#   dotzero <- .zero
#   M1 <- 2
#   eval(negzero.expression.VGAM)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
#        dpqrfun = "zz",
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("mean", "sigma"),
         zero = .zero )
  }, list( .zero = zero ))),

  initialize = eval(substitute(expression({
#   if (NCOL(y) != 1)
#     stop("the response must be a vector or one-column matrix")

# 20120601:
#   temp5 <-
    w.y.check(w = w, y = y)
#             Is.positive.y = TRUE,
#             Is.nonnegative.y = TRUE,
#             ncol.w.max = Inf,
#             ncol.y.max = Inf,
#             Is.integer.y = TRUE,
#             out.wy = TRUE,
#             colsyperw = 1,
#             maximize = TRUE,
#   w <- temp5$w
#   y <- temp5$y
#print("head(w)")
#print( head(w) )
#print("head(y)")
#print( head(y) )


    predictors.names <-
      c(namesof("mean",  .lmean  , earg = .emean  , tag = FALSE),
        namesof("sigma", .lsigma , earg = .esigma , tag = FALSE))


    if (!length(etastart)) {
      sigma.init <- if (length( .isigma )) rep_len( .isigma , n) else {
# Eqn (2.3)
        hh <- 2 - .d
        KK <- 1 / (1 + 1/hh + 0.75/hh^2)
        K2 <- 1 + 3/hh + 15/(4*hh^2)
        rep_len(sqrt(var(y) / (KK*K2)), n)
      }
      mean.init <- rep_len(weighted.mean(y, w), n)
      etastart <-
        cbind(theta2eta(mean.init,  .lmean  , earg = .emean  ),
              theta2eta(sigma.init, .lsigma , earg = .esigma ))
    }
  }),list( .lmean = lmean, .lsigma = lsigma,
                           .isigma = isigma, .d = d,
           .emean = emean, .esigma = esigma ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta[, 1], .lmean , earg = .emean )
  }, list( .lmean = lmean,
           .emean = emean, .esigma = esigma ))),
  last = eval(substitute(expression({
      misc$link <-    c("mean" = .lmean , "sigma"= .lsigma )

      misc$earg <- list("mean" = .emean , "sigma"= .esigma )

      misc$expected <- TRUE
      misc$d <- .d
  }), list( .lmean = lmean, .lsigma = lsigma, .d = d,
            .emean = emean, .esigma = esigma ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    mymu  <- eta2theta(eta[, 1], .lmean  , earg = .emean  )
    sigma <- eta2theta(eta[, 2], .lsigma , earg = .esigma )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dtikuv(x = y, d = .d , mean = mymu,
                               sigma = sigma, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lmean = lmean, .lsigma = lsigma, .d = d,
           .emean = emean, .esigma = esigma ))),
  vfamily = c("tikuv"),
# ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
  validparams = eval(substitute(function(eta, y, extra = NULL) {
# 20160618;
    mymu  <- eta2theta(eta[, 1], .lmean  , earg = .emean  )
    sigma <- eta2theta(eta[, 2], .lsigma , earg = .esigma )
    dee   <- .d
#print("hi5a")
    okay1 <- all(is.finite(mymu )) &&
             all(is.finite(sigma)) && all(0 < sigma) &&
             all(is.finite(dee  )) && all(0 < dee & dee < 2)
#print("okay1 in @validparams in tikuv()")
#print( okay1 )
    okay1
  }, list( .lmean = lmean, .lsigma = lsigma, .d = d,
           .emean = emean, .esigma = esigma ))),
# ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,





# simslot = eval(substitute(
# function(object, nsim) {
##20140102; does not work because 'mean' argument must be scalar
#
#   pwts <- if (length(pwts <- object@prior.weights) > 0)
#             pwts else weights(object, type = "prior")
#   if (any(pwts != 1))
#     warning("ignoring prior weights")
#   eta <- predict(object)
#   mymu  <- eta2theta(eta[, 1], .lmean  , earg = .emean  )
#   sigma <- eta2theta(eta[, 2], .lsigma , earg = .esigma )
#   rtikuv(nsim * length(mymu), d = .d , mean = mymu,
#          sigma = sigma)
# }, list( .lmean = lmean, .lsigma = lsigma, .d = d,
#          .emean = emean, .esigma = esigma ))),










  deriv = eval(substitute(expression({
    mymu  <- eta2theta(eta[, 1], .lmean  , earg = .emean  )
    sigma <- eta2theta(eta[, 2], .lsigma , earg = .esigma )

    dmu.deta <- dtheta.deta(mymu, .lmean , earg = .emean )
    dsigma.deta <- dtheta.deta(sigma, .lsigma, earg = .esigma)

    zedd <- (y - mymu) / sigma
    hh <- 2 - .d
    gzedd <- zedd / (1 + 0.5*zedd^2 / hh)

    dl.dmu <- zedd / sigma - 2 * gzedd / (hh*sigma)
    dl.dsigma <- (zedd^2 - 1 - 2 * zedd * gzedd / hh) / sigma
#   dl.dsigma <- -1/sigma + zedd^2 / sigma - 2 * zedd * gzedd / (hh*sigma)

    c(w) * cbind(dl.dmu    * dmu.deta,
                 dl.dsigma * dsigma.deta)
  }), list( .lmean = lmean, .lsigma = lsigma, .d = d,
            .emean = emean, .esigma = esigma ))),
  weight = eval(substitute(expression({
# Notation: ayy=a, Dnos=D, DDDD=D^*
    ayy <- 1 / (2*hh)
    Dnos <- 1 - (2/hh) * (1 - ayy) / (1 + 2*ayy + 3*ayy^2)
    Dstar <- -1 + 3 * (1 + 2*ayy + 11*ayy^2) / (1 + 2*ayy + 3*ayy^2)

    ned2l.dmymu2 <- Dnos / sigma^2
    ned2l.dnu2   <- Dstar / sigma^2

    wz <- matrix(NA_real_, n, M)  # diagonal matrix
    wz[, iam(1, 1, M)] <- ned2l.dmymu2 * dmu.deta^2
    wz[, iam(2, 2, M)] <- ned2l.dnu2 * dsigma.deta^2
    c(w) * wz
  }), list( .lmean = lmean, .lsigma = lsigma,
            .emean = emean, .esigma = esigma ))))
}  # tikuv



# ==========================================================


# Works;
# Notes:
# 1. NR = Fisher scoring only if the scale parameter
#    enters in the form exp(- a y2).
# 2. Implemented has the scale parameter enters in
#    in the form exp(- y2 / a).
#    Then NR is not FS.

# 20050806
# Last modified: 20100107,



 bigamma.mckay <-
    function(lscale = "loglink",
             lshape1 = "loglink",
             lshape2 = "loglink",
             iscale = NULL,
             ishape1 = NULL,
             ishape2 = NULL,
             imethod = 1,
             zero = "shape") {
#            zero = 2:3
  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")

  lshape1 <- as.list(substitute(lshape1))
  eshape1 <- link2list(lshape1)
  lshape1 <- attr(eshape1, "function.name")

  lshape2 <- as.list(substitute(lshape2))
  eshape2 <- link2list(lshape2)
  lshape2 <- attr(eshape2, "function.name")


  if (!is.null(iscale))
    if (!is.Numeric(iscale, positive = TRUE))
      stop("argument 'iscale' must be positive or NULL")
  if (!is.null(ishape1))
    if (!is.Numeric(ishape1, positive = TRUE))
      stop("argument 'ishape1' must be positive or NULL")
  if (!is.null(ishape2))
    if (!is.Numeric(ishape2, positive = TRUE))
      stop("argument 'ishape2' must be positive or NULL")

  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2.5)
    stop("argument 'imethod' must be 1 or 2")



  new("vglmff",
  blurb = c("Bivariate gamma: McKay's distribution\n",
            "Links:    ",
            namesof("scale",  lscale,  earg = escale ), ", ",
            namesof("shape1", lshape1, earg = eshape1), ", ",
            namesof("shape2", lshape2, earg = eshape2)),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x,
                                .zero , M = M,
                      predictors.names = predictors.names,
                                M1 = 3)
#   dotzero <- .zero
#   M1 <- 3
#   Q1 <- 2
#   eval(negzero.expression.VGAM)
  }), list( .zero = zero ))),



  infos = eval(substitute(function(...) {
    list(M1 = 3,
         Q1 = 2,
#        dpqrfun = "zz",
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("scale", "shape1", "shape2"),
         lscale  = .lscale  ,
         lshape1 = .lshape1 ,
         lshape2 = .lshape2 ,
         zero = .zero )
  },
  list( .zero = zero,
        .lscale  = lscale ,
        .lshape1 = lshape1,
        .lshape2 = lshape2 ))),


  initialize = eval(substitute(expression({
#   if (!is.matrix(y) || ncol(y) != 2)
#     stop("the response must be a 2 column matrix")

    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
#             Is.nonnegative.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = 2,
              ncol.y.min = 2,
#             Is.integer.y = TRUE,
              out.wy = TRUE,
              colsyperw = 2,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


# For @linkinv
    extra$colnames.y  <- colnames(y)

    if (any(y[, 1] >= y[, 2]))
      stop("the second column minus the first column must ",
           "be a vector of positive values")


    predictors.names <-
      c(namesof("scale",  .lscale,  .escale,  short = TRUE),
        namesof("shape1", .lshape1, .eshape1, short = TRUE),
        namesof("shape2", .lshape2, .eshape2, short = TRUE))

    if (!length(etastart)) {
      momentsY <- if ( .imethod == 1) {
        cbind(median(y[, 1]),  # This may not be monotonic
              median(y[, 2])) + 0.01
      } else {
        cbind(weighted.mean(y[, 1], w),
              weighted.mean(y[, 2], w))
      }

      mcg2.loglik <- function(thetaval, y, x, w, extraargs) {
        ainit <- a <- thetaval
        momentsY <- extraargs$momentsY
          p <- (1/a) * abs(momentsY[1]) + 0.01
          q <- (1/a) * abs(momentsY[2] - momentsY[1]) + 0.01
          sum(c(w) * (-(p+q)*log(a) - lgamma(p) - lgamma(q) +
               (p - 1)*log(y[, 1]) +
               (q - 1)*log(y[, 2]-y[, 1]) - y[, 2] / a ))
      }

      a.grid <- if (length( .iscale )) c( .iscale ) else
                c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5,
                  10, 20, 50, 100)
      extraargs <- list(momentsY = momentsY)
      ainit <- grid.search(a.grid, objfun = mcg2.loglik,
                           y = y, x = x, w = w,
                           maximize = TRUE,
                           extraargs = extraargs)
      ainit <- rep_len(if (is.Numeric( .iscale ))
                           .iscale else ainit, n)
      pinit <- (1/ainit) * abs(momentsY[1]) + 0.01
      qinit <- (1/ainit) *
          abs(momentsY[2] - momentsY[1]) + 0.01

      pinit <- rep_len(if (is.Numeric( .ishape1 ))
                           .ishape1 else pinit, n)
      qinit <- rep_len(if (is.Numeric( .ishape2 ))
                           .ishape2 else qinit, n)

#print("c(ainit[1], pinit[1], qinit[1])")
#print( c(ainit[1], pinit[1], qinit[1]) )
      etastart <-
        cbind(theta2eta(ainit, .lscale),
              theta2eta(pinit, .lshape1),
              theta2eta(qinit, .lshape2))
    }
  }),
  list(
    .lscale = lscale, .lshape1 = lshape1, .lshape2 = lshape2,
    .escale = escale, .eshape1 = eshape1, .eshape2 = eshape2,
    .iscale = iscale, .ishape1 = ishape1, .ishape2 = ishape2,
    .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    NOS <- NCOL(eta) / c(M1 = 3)
    a <- eta2theta(eta[, 1], .lscale  ,  .escale )
    p <- eta2theta(eta[, 2], .lshape1 , .eshape1 )
    q <- eta2theta(eta[, 3], .lshape2 , .eshape2 )
    fv.mat <-  cbind("y1" = p*a,
                     "y2" = (p+q)*a)  # Overwrite the colnames:
    label.cols.y(fv.mat, colnames.y = extra$colnames.y,
                 NOS = NOS)
  },
  list(
      .lscale = lscale, .lshape1 = lshape1, .lshape2 = lshape2,
      .escale = escale, .eshape1 = eshape1, .eshape2 = eshape2
       ))),
  last = eval(substitute(expression({
    misc$link <-    c("scale"  = .lscale ,
                      "shape1" = .lshape1 ,
                      "shape2" = .lshape2 )

    misc$earg <- list("scale"  = .escale ,
                      "shape1" = .eshape1 ,
                      "shape2" = .eshape2 )

    misc$ishape1 <- .ishape1
    misc$ishape2 <- .ishape2
    misc$iscale <- .iscale
    misc$expected <- TRUE
    misc$multipleResponses <- FALSE
  }),
  list(
      .lscale = lscale, .lshape1 = lshape1, .lshape2 = lshape2,
      .escale = escale, .eshape1 = eshape1, .eshape2 = eshape2,
      .iscale = iscale, .ishape1 = ishape1, .ishape2 = ishape2,
      .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    a <- eta2theta(eta[, 1], .lscale  ,  .escale )
    p <- eta2theta(eta[, 2], .lshape1 , .eshape1 )
    q <- eta2theta(eta[, 3], .lshape2 , .eshape2 )

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        c(w) * (-(p+q)*log(a) - lgamma(p) - lgamma(q) +
                (p - 1)*log(y[, 1]) +
                (q - 1)*log(y[, 2]-y[, 1]) -
               y[, 2] / a)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
    },
  list(
      .lscale = lscale, .lshape1 = lshape1, .lshape2 = lshape2,
      .escale = escale, .eshape1 = eshape1, .eshape2 = eshape2
       ))),
  vfamily = c("bigamma.mckay"),
# ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
  validparams = eval(substitute(function(eta, y, extra = NULL) {
# 20160618;
    aparam <- eta2theta(eta[, 1], .lscale  ,  .escale )
    shape1 <- eta2theta(eta[, 2], .lshape1 , .eshape1 )
    shape2 <- eta2theta(eta[, 3], .lshape2 , .eshape2 )
#print("hi5a")
    okay1 <- all(is.finite(aparam)) && all(0 < aparam) &&
             all(is.finite(shape1)) && all(0 < shape1) &&
             all(is.finite(shape2)) && all(0 < shape2)
#print("okay1 in @validparams in bigamma.mckay()")
#print( okay1 )
    okay1
  },
  list(
    .lscale = lscale, .lshape1 = lshape1, .lshape2 = lshape2,
    .escale = escale, .eshape1 = eshape1, .eshape2 = eshape2
       ))),
# ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
  deriv = eval(substitute(expression({
    aparam <- eta2theta(eta[, 1], .lscale  ,  .escale )
    shape1 <- eta2theta(eta[, 2], .lshape1 , .eshape1 )
    shape2 <- eta2theta(eta[, 3], .lshape2 , .eshape2 )

    dl.da <- (-(shape1+shape2) + y[, 2] / aparam) / aparam
    dl.dshape1 <- -log(aparam) - digamma(shape1) + log(y[, 1])
    dl.dshape2 <- -log(aparam) - digamma(shape2) + log(y[, 2] -
                                                       y[, 1])

    c(w) * cbind(dl.da      * dtheta.deta(aparam, .lscale),
                 dl.dshape1 * dtheta.deta(shape1, .lshape1),
                 dl.dshape2 * dtheta.deta(shape2, .lshape2))
  }),
  list(
      .lscale = lscale, .lshape1 = lshape1, .lshape2 = lshape2,
      .escale = escale, .eshape1 = eshape1, .eshape2 = eshape2
       ))),
  weight = eval(substitute(expression({
    d11 <- (shape1+shape2) / aparam^2
    d22 <- trigamma(shape1)
    d33 <- trigamma(shape2)
    d12 <- 1 / aparam
    d13 <- 1 / aparam
    d23 <- 0

    wz <- matrix(0, n, dimm(M))
    wz[, iam(1, 1, M)] <- dtheta.deta(aparam, .lscale  )^2 * d11
    wz[, iam(2, 2, M)] <- dtheta.deta(shape1, .lshape1 )^2 * d22
    wz[, iam(3, 3, M)] <- dtheta.deta(shape2, .lshape2 )^2 * d33
    wz[, iam(1, 2, M)] <- dtheta.deta(aparam, .lscale  ) *
                          dtheta.deta(shape1, .lshape1 ) * d12
    wz[, iam(1, 3, M)] <- dtheta.deta(aparam, .lscale  ) *
                          dtheta.deta(shape2, .lshape2 ) * d13
    wz[, iam(2, 3, M)] <- dtheta.deta(shape1, .lshape1 ) *
                          dtheta.deta(shape2, .lshape2 ) * d23

    c(w) * wz
  }),
  list( .lscale = lscale, .lshape1 = lshape1,
                          .lshape2 = lshape2 ))))
}  # bigamma.mckay



# ==========================================================



# ==========================================================




# ==========================================================




# ==========================================================



# ==========================================================



