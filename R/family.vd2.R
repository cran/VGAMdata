# Copyright (C) 1998-2020. T. W. Yee    All rights reserved.
# 20201207: family.vd2.R is for \pkg{VGAMdata}.
#   This file contains any 'old' functions from \pkg{VGAM} such
#   as ones made obsolete by new GAIT family functions.



# Last modified:
# 20201207: Trying to put [dpqr]tikuv(), etc. here.
#   Thats because family.vd1.R is getting big.






# ==================================================================


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



# ==================================================================





# ==================================================================



# ==================================================================




# ==================================================================




# ==================================================================



# ==================================================================



