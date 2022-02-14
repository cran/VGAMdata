# These functions are
# Copyright (C) 1998-2022 T.W. Yee, University of Auckland.
# All rights reserved.















dgenpois <- function(x, lambda = 0, theta, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(lambda), length(theta))
  if (length(x)      != LLL) x      <- rep_len(x,      LLL)
  if (length(lambda) != LLL) lambda <- rep_len(lambda, LLL)
  if (length(theta)  != LLL) theta  <- rep_len(theta,  LLL)

  llans <- -x*lambda - theta + (x-1) * log(theta + x*lambda) +
           log(theta) - lgamma(x+1)
  llans[x < 0] <- log(0)
  llans[x != round(x)] <- log(0)  # x should be integer-valued
  llans[lambda > 1] <- NaN
  if (any(ind1 <- (lambda < 0))) {
    epsilon <- 1.0e-9  # Needed to handle a "<" rather than a "<=".
    mmm <- pmax(4, floor(theta/abs(lambda) - epsilon))
    llans[ind1 & mmm < pmax(-1, -theta/mmm)] <- NaN
    llans[ind1 & mmm < x] <- log(0)  # probability 0, not NaN
  }
  if (log.arg) {
    llans
  } else {
    exp(llans)
  }
}  # dgenpois








 genpoisson <-
  function(llambda = "rhobitlink",
           ltheta = "loglink",
           ilambda = NULL, itheta = NULL,  # use.approx = TRUE,
           imethod = 1,
           ishrinkage = 0.95,
           zero = "lambda") {



  llambda <- as.list(substitute(llambda))
  elambda <- link2list(llambda)
  llambda <- attr(elambda, "function.name")

  ltheta <- as.list(substitute(ltheta))
  etheta <- link2list(ltheta)
  ltheta <- attr(etheta, "function.name")

  if (!is.Numeric(ishrinkage, length.arg = 1) ||
     ishrinkage < 0 ||
     ishrinkage > 1)
    stop("bad input for argument 'ishrinkage'")

  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")

  new("vglmff",
  blurb = c("Generalized Poisson distribution\n\n",
            "Links:    ",
            namesof("lambda", llambda, earg = elambda), ", ",
            namesof("theta",  ltheta,  earg = etheta ), "\n",
            "Mean:     theta / (1-lambda)\n",
            "Variance: theta / (1-lambda)^3"),
 constraints = eval(substitute(expression({

    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = FALSE,
         multipleResponses = TRUE,
         parameters.names = c("lambda", "theta"),
         imethod = .imethod ,
         zero = .zero )
  }, list( .zero = zero,
           .imethod = imethod ))),

  initialize = eval(substitute(expression({
    temp5 <-
    w.y.check(w = w, y = y,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,  # 1,
              ncol.y.max = Inf,  # 1,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y

    extra$ncoly <- ncoly <- NOS <- ncol(y)
    extra$M1 <- M1 <- 2
    M <- M1 * ncoly
    mynames1 <- param.names("lambda", NOS, skip1 = TRUE)
    mynames2 <- param.names("theta",  NOS, skip1 = TRUE)

    predictors.names <-
       c(namesof(mynames1, .llambda , earg = .elambda , tag = FALSE),
         namesof(mynames2, .ltheta  , earg = .etheta  , tag = FALSE))
    predictors.names <- predictors.names[interleave.VGAM(M, M1 = M1)]

    init.lambda <- init.theta <- matrix(0, n, NOS)
    for (spp. in 1: NOS) {
      init.lambda[, spp.] <- if ( .imethod == 1) {
        min(max(0.05,
                1 - sqrt(weighted.mean(y[, spp.],
                                       w[, spp.]) / var(y[, spp.]))),
            0.95)
      } else if ( .imethod == 2) {
        runif(n, max = 0.1)
      } else {
        runif(n, max = 0.7)
      }

      init.theta[, spp.]  <- if ( .imethod == 2) {
        (y[, spp.] + weighted.mean(y[, spp.], w[, spp.])) / 2
      } else if ( .imethod == 3) {
        (y[, spp.] + median(y[, spp.])) / 2
      } else {
        (1 - .ishrinkage ) * y[, spp.] +
             .ishrinkage   * weighted.mean(y[, spp.], w[, spp.])
      }
    }

    if (!length(etastart)) {
      init.lambda <- if (length( .ilambda ))
                       matrix( .ilambda , n, NOS, byrow = TRUE) else
                       init.lambda
      init.theta  <- if (length( .itheta ))
                       matrix( .itheta  , n, NOS, byrow = TRUE) else
                       init.theta
      etastart <-
        cbind(theta2eta(init.lambda, .llambda , earg = .elambda ),
              theta2eta(init.theta,  .ltheta  , earg = .etheta  ))
      etastart <- etastart[, interleave.VGAM(M, M1 = M1), drop = FALSE]
    }
  }), list( .ltheta = ltheta, .llambda = llambda,
            .etheta = etheta, .elambda = elambda,
            .itheta = itheta, .ilambda = ilambda,
            .imethod = imethod, .ishrinkage = ishrinkage)) ),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    lambda <- eta2theta(eta[, c(TRUE, FALSE)], .llambda , earg = .elambda )
    theta  <- eta2theta(eta[, c(FALSE, TRUE)], .ltheta  , earg = .etheta  )
    theta / (1 - lambda)
  }, list( .ltheta = ltheta, .llambda = llambda,
           .etheta = etheta, .elambda = elambda ))),
  last = eval(substitute(expression({
    M1 <- extra$M1

    temp.names <- c(mynames1, mynames2)
    temp.names <- temp.names[interleave.VGAM(M1 * ncoly, M1 = M1)]

    misc$link <- rep_len( .llambda , M1 * ncoly)
    misc$earg <- vector("list", M1 * ncoly)
    names(misc$link) <-
    names(misc$earg) <- temp.names
    for (ii in 1:ncoly) {
      misc$link[ M1*ii-1 ] <- .llambda
      misc$link[ M1*ii   ] <- .ltheta
      misc$earg[[M1*ii-1]] <- .elambda
      misc$earg[[M1*ii  ]] <- .etheta
    }
  }), list( .ltheta = ltheta, .llambda = llambda,
            .etheta = etheta, .elambda = elambda,
            .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    lambda <- eta2theta(eta[, c(TRUE, FALSE)], .llambda , earg = .elambda )
    theta  <- eta2theta(eta[, c(FALSE, TRUE)], .ltheta  , earg = .etheta  )
    index <- (y == 0)
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- dgenpois(x = y, lambda = lambda, theta = theta,
                          log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .ltheta = ltheta, .llambda = llambda,
           .etheta = etheta, .elambda = elambda ))),
   vfamily = c("genpoisson"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    lambda <- eta2theta(eta[, c(TRUE, FALSE)], .llambda , earg = .elambda )
    theta  <- eta2theta(eta[, c(FALSE, TRUE)], .ltheta  , earg = .etheta  )
    mmm <- ifelse(lambda < 0, floor(-theta/lambda), Inf)
    if (any(mmm < 4)) {
      warning("the lower bound is less than 4; choosing 4")
      mmm <- pmax(mmm, 4)
    }
    Lbnd <- pmax(-1, -theta / mmm)
    okay1 <- all(is.finite(lambda)) && all(Lbnd < lambda & lambda < 1) &&
             all(is.finite(theta )) && all(0 < theta)
    okay1
  }, list( .ltheta = ltheta, .llambda = llambda,
           .etheta = etheta, .elambda = elambda ))),
  deriv = eval(substitute(expression({
    M1  <- 2
    NOS <- ncol(eta)/M1

    lambda <- eta2theta(eta[, c(TRUE, FALSE)], .llambda , earg = .elambda )
    theta  <- eta2theta(eta[, c(FALSE, TRUE)], .ltheta  , earg = .etheta  )
    dl.dlambda <- -y + y*(y-1) / (theta+y*lambda)
    dl.dtheta  <- -1 +   (y-1) / (theta+y*lambda) + 1/theta
    dTHETA.deta  <- dtheta.deta(theta,  .ltheta  , earg = .etheta  )
    dlambda.deta <- dtheta.deta(lambda, .llambda , earg = .elambda )
    myderiv <- c(w) * cbind(dl.dlambda * dlambda.deta,
                            dl.dtheta  * dTHETA.deta )
    myderiv[, interleave.VGAM(M, M1 = M1)]
  }), list( .ltheta = ltheta, .llambda = llambda,
            .etheta = etheta, .elambda = elambda ))),
  weight = eval(substitute(expression({
    wz <- matrix(0, n, M + M-1)  # Tridiagonal
    ned2l.dlambda2 <- theta / (1 - lambda) +
                      2 * theta / (theta + 2 * lambda)
    ned2l.dtheta2 <- 1 / theta - lambda / (theta + 2 * lambda)
    ned2l.dthetalambda <- theta / (theta + 2 * lambda)
    wz[, M1*(1:NOS) - 1    ] <- ned2l.dlambda2 * dlambda.deta^2
    wz[, M1*(1:NOS)        ] <- ned2l.dtheta2 * dTHETA.deta^2
    wz[, M1*(1:NOS) + M - 1] <- ned2l.dthetalambda * dTHETA.deta *
                                                     dlambda.deta
    wz <- w.wz.merge(w = w, wz = wz, n = n, M = M + (M - 1),
                     ndepy = NOS)
    wz
  }), list( .ltheta = ltheta, .llambda = llambda,
            .etheta = etheta, .elambda = elambda ))))
}  # genpoisson










 yip88 <- function(link = "loglink", n.arg = NULL, imethod = 1) {








  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  new("vglmff",
  blurb = c("Zero-inflated Poisson (based on Yip (1988))\n\n",
            "Link:     ",
            namesof("lambda", link, earg), "\n",
            "Variance: (1 - pstr0) * lambda"),
  first = eval(substitute(expression({
    zero <- y == 0
    if (any(zero)) {
      if (length(extra)) extra$sumw <- sum(w) else
        extra <- list(sumw=sum(w))
      if (is.numeric(.n.arg) && extra$sumw != .n.arg)
        stop("value of 'n.arg' conflicts with data ",
             "(it need not be specified anyway)")
      warning("trimming out the zero observations")


      axa.save <-  attr(x, "assign")
      x <- x[!zero,, drop = FALSE]
      attr(x, "assign") <- axa.save    # Don't lose these!!
      w <- w[!zero]
      y <- y[!zero]
    } else {
      if (!is.numeric(.n.arg))
        stop("n.arg must be supplied")
    }

  }), list( .n.arg = n.arg ))),

  initialize = eval(substitute(expression({
    narg <- if (is.numeric(.n.arg)) .n.arg else extra$sumw
    if (sum(w) > narg)
      stop("sum(w) > narg")

    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1)


    predictors.names <-
      namesof("lambda", .link, list(theta = NULL), tag = FALSE)

    if (!length(etastart)) {
      lambda.init <- Init.mu(y = y, w = w, imethod = .imethod ,  # x = x,
                             pos.only = FALSE)
      etastart <- theta2eta(lambda.init, .link , earg = .earg )
    }
    if (length(extra)) {
      extra$sumw <- sum(w)
      extra$narg <- narg   # For @linkinv
    } else {
      extra <- list(sumw = sum(w), narg = narg)
    }
  }), list( .link = link, .earg = earg,
            .n.arg = n.arg, .imethod = imethod ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    lambda <- eta2theta(eta, .link, .earg)
    temp5 <- exp(-lambda)
    pstr0 <- (1 - temp5 - extra$sumw/extra$narg) / (1 - temp5)
    if (any(pstr0 <= 0))
      stop("non-positive value(s) of pstr0")
    (1 - pstr0) * lambda
  }, list( .link = link, .earg = earg ))),

  last = eval(substitute(expression({
    misc$link <-    c(lambda = .link )

    misc$earg <- list(lambda = .earg )

    if (intercept.only) {
      suma <- extra$sumw
      pstr0 <- (1 - temp5[1] - suma / narg) / (1 - temp5[1])
      pstr0 <- if (pstr0 < 0 || pstr0 > 1) NA else pstr0
      misc$pstr0 <- pstr0
    }
  }), list( .link = link, .earg = earg ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    lambda <- eta2theta(eta, .link)
    temp5 <- exp(-lambda)
    pstr0 <- (1 - temp5 - extra$sumw / extra$narg) / (1 - temp5)
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) *
               dzipois(x = y, pstr0 = pstr0, lambda = lambda, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .link = link, .earg = earg ))),

  vfamily = c("yip88"),

  validparams = eval(substitute(function(eta, y, extra = NULL) {
    lambda <- eta2theta(eta, .link , earg = .earg )
    okay1 <- all(is.finite(lambda)) && all(0 < lambda)
    okay1
  }, list( .link = link, .earg = earg ))),

  deriv = eval(substitute(expression({
    lambda <- eta2theta(eta, .link , earg = .earg )
    temp5 <- exp(-lambda)
    dl.dlambda <- -1 + y/lambda - temp5/(1-temp5)
    dlambda.deta <- dtheta.deta(lambda, .link , earg = .earg )
    w * dl.dlambda * dlambda.deta
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({
    d2lambda.deta2 <- d2theta.deta2(lambda, .link , earg = .earg )
    d2l.dlambda2 <- -y / lambda^2 + temp5 / (1 - temp5)^2
    -w * (d2l.dlambda2*dlambda.deta^2 + dl.dlambda*d2lambda.deta2)
  }), list( .link = link, .earg = earg ))))
}  # yip88









doiposbinom <- function(x, size, prob, pstr1 = 0, log = FALSE) {



  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(size), length(prob), length(pstr1))
  if (length(x)      != LLL) x      <- rep_len(x,      LLL)
  if (length(size)   != LLL) size   <- rep_len(size,   LLL)
  if (length(prob)   != LLL) prob   <- rep_len(prob,   LLL)
  if (length(pstr1)  != LLL) pstr1  <- rep_len(pstr1,  LLL)

  ans <- x  # + prob + pstr1
  index1 <- (x == 1)
  if (log.arg) {
    ans[ index1] <- log(pstr1[ index1] + (1 - pstr1[ index1]) *
                        dposbinom(x[ index1], size[ index1],
                                  prob[ index1]))
    ans[!index1] <- log1p(-pstr1[!index1]) +
                    dposbinom(x[!index1], size[!index1],
                              prob[!index1], log = TRUE)
  } else {
    ans[ index1] <-      pstr1[ index1] + (1 - pstr1[ index1]) *
                        dposbinom(x[ index1], size[ index1], prob[ index1])
    ans[!index1] <- (1 - pstr1[!index1]) *
                        dposbinom(x[!index1], size[!index1], prob[!index1])
  }


  deflat.limit <- size * prob / (1 + (size-1) * prob -
                                 1 / (1-prob)^(size-1))
  ans[pstr1 < deflat.limit] <- NaN
  ans[1 < pstr1] <- NaN

  ans
}  # doiposbinom



poiposbinom <- function(q, size, prob, pstr1 = 0) {

  LLL <- max(length(q), length(size), length(prob), length(pstr1))
  if (length(q)      != LLL) q      <- rep_len(q,      LLL)
  if (length(size)   != LLL) size   <- rep_len(size,   LLL)
  if (length(prob)   != LLL) prob   <- rep_len(prob,   LLL)
  if (length(pstr1)  != LLL) pstr1  <- rep_len(pstr1,  LLL)
  ans <- rep_len(NA_real_, LLL)

  ans <- pposbinom(q, size, prob)  # lower.tail=lower.tail, log.p=log.p
  ans <- ifelse(q < 1, 0, pstr1 + (1 - pstr1) * ans)

  deflat.limit <- size * prob / (1 + (size-1) * prob -
                                 1 / (1-prob)^(size-1))
  ans[pstr1 < deflat.limit] <- NaN
  ans[1 < pstr1] <- NaN

  ans
}




qoiposbinom <- function(p, size, prob, pstr1 = 0) {

  LLL <- max(length(p), length(size), length(prob), length(pstr1))
  if (length(p)      != LLL) p      <- rep_len(p,      LLL)
  if (length(size)   != LLL) size   <- rep_len(size,   LLL)
  if (length(prob)   != LLL) prob   <- rep_len(prob,   LLL)
  if (length(pstr1)  != LLL) pstr1  <- rep_len(pstr1,  LLL)
  ans <- rep_len(NA_real_, LLL)
  deflat.limit <- size * prob / (1 + (size-1) * prob -
                                 1 / (1-prob)^(size-1))

  ans[p <= pstr1] <- 1
  pindex <- (deflat.limit <= pstr1) & (pstr1 < p)
  ans[pindex] <-
    qposbinom((p[pindex] - pstr1[pindex]) / (1 - pstr1[pindex]),
               size = size[pindex],
               prob = prob[pindex])

  ans[p == 0] <- 1
  ans[prob == 0] <- NaN

  ans[pstr1 < deflat.limit] <- NaN
  ans[1 < pstr1] <- NaN

  ans[p < 0] <- NaN
  ans[1 < p] <- NaN
  ans
}  # qoiposbinom



roiposbinom <- function(n, size, prob, pstr1 = 0) {

  qoiposbinom(runif(n), size, prob, pstr1 = pstr1)
}  # roiposbinom



 oiposbinomial <-
  function(lpstr1 = "logitlink", lprob = "logitlink",
           type.fitted = c("mean", "prob", "pobs1", "pstr1", "onempstr1"),
           iprob = NULL,
    gpstr1 = ppoints(9),  # (1:19)/20,
    gprob  = ppoints(9),  # (1:19)/20,  # 20160613; grid for finding prob
           multiple.responses = FALSE,
           zero = NULL) {

  gprobb <- gprob

  lpstr1 <- as.list(substitute(lpstr1))
  epstr1 <- link2list(lpstr1)
  lpstr1 <- attr(epstr1, "function.name")

  lprobb <- as.list(substitute(lprob))
  eprobb <- link2list(lprobb)
  lprobb <- attr(eprobb, "function.name")



  type.fitted <- match.arg(type.fitted,
                   c("mean", "prob", "pobs1", "pstr1", "onempstr1"))[1]

  iprobb <- iprob
  if (length(iprobb))
    if (!is.Numeric(iprobb, positive = TRUE) ||
        any(iprobb >= 1))
      stop("argument 'iprob' values must be in (0, 1)")


  new("vglmff",
  blurb = c("One-inflated positive binomial\n\n",
            "Links:    ",
            namesof("pstr1", lpstr1, earg = epstr1 ), ", ",
            namesof("prob",  lprobb, earg = eprobb ), "\n",
            "Mean:     pstr1 + (1 - pstr1) * ",
                      "size * prob / (1 - (1-prob)^size)"),

  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = .multiple.responses ,  # FALSE,  # TRUE,
         parameters.names = c("pstr1", "prob"),
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .multiple.responses = multiple.responses,
           .type.fitted = type.fitted
         ))),
  initialize = eval(substitute(expression({
    M1 <- 2
    multiple.responses <- .multiple.responses
    y <- as.matrix(y)
    w.orig <- as.matrix(w)  # zz this may be of a weird dimension

    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              Is.integer.y =  if (multiple.responses) FALSE else FALSE,
              ncol.w.max =  if (multiple.responses) ncol(y) else 1,
              ncol.y.max =  if (multiple.responses) ncol(y) else ncol(y),
              out.wy = TRUE,
              colsyperw = if (multiple.responses) 1 else ncol(y),
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    if (multiple.responses) {
      if (!all(w == round(w)))
        stop("the 'weights' argument must be integer valued")
      if (min(y) < 0 || max(y) > 1)
        stop("the response must be a proportion")
      Nvec <- w
    } else {
      if (ncol(y) > 1) {
        Nvec <- rowSums(y)
        y[, 1] <- y[, 1] / Nvec
        y <- y[, 1, drop = FALSE]
        w[, 1] <- w[, 1] * Nvec  # == w.orig * Nvec
        w <- w[, 1, drop = FALSE]
      } else {
        Nvec <- w  # rep_len(1, nrow(x))
        if (!all(Nvec == round(Nvec)))
          stop("number of trials is not integer-valued")
      }
    }
    extra$Nvec <- Nvec
    w.orig <- matrix(w.orig, n, ncol(y))



    NOS <- ncoly <- ncol(y)
    extra$ncoly <- ncoly
    M <- M1 * ncoly
    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)

    mynames1 <- param.names("pstr1", ncoly, skip1 = TRUE)
    mynames2 <- param.names("prob",  ncoly, skip1 = TRUE)
    predictors.names <-
        c(namesof(mynames1, .lpstr1 , earg = .epstr1 , tag = FALSE),
          namesof(mynames2, .lprobb , earg = .eprobb , tag = FALSE))[
          interleave.VGAM(M, M1 = M1)]


    if (!length(etastart)) {
      probb.init <-
      pstr1.init <- matrix(NA_real_, n, NOS)
      gpstr1 <- .gpstr1
      gprobb <- .gprobb
      iprobb <- .iprobb
      if (length(iprobb))
        gprobb <- iprobb

      oiposbinom.Loglikfun <- function(pstr1, prob, y, x, w, extraargs) {
      sum(c(w) * doiposbinom(x = y, pstr1 = pstr1, size = extraargs$size,
                             prob = prob, log = TRUE))
      }


      for (jay in 1:NOS) {  # For each response 'y_jay'... do:


        try.this <-
          grid.search2(gpstr1, gprobb,
                       objfun = oiposbinom.Loglikfun,
                       y = round(y[, jay] * Nvec[, jay]),
                     w = 1,  # w.orig[, jay], or 1, or w[, jay], possibly
                       extraargs = list(size = Nvec),
                       ret.objfun = TRUE)  # Last value is the loglik
        pstr1.init[, jay] <-  try.this["Value1"]
        probb.init[, jay] <-  try.this["Value2"]
      }  # for (jay ...)

      etastart <- cbind(theta2eta(pstr1.init, .lpstr1 , earg = .epstr1 ),
                        theta2eta(probb.init, .lprobb , earg = .eprobb ))[,
                        interleave.VGAM(M, M1 = M1)]
      mustart <- NULL  # Since etastart has been computed.
    }  # End of !length(etastart)
  }), list( .lpstr1 = lpstr1, .lprobb = lprobb,
            .epstr1 = epstr1, .eprobb = eprobb,
                              .iprobb = iprobb,
            .gpstr1 = gpstr1, .gprobb = gprobb,
            .multiple.responses = multiple.responses,
            .type.fitted = type.fitted ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    NOS <- ncol(eta) / c(M1 = 2)
    type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "prob", "pobs1", "pstr1", "onempstr1"))[1]

    pstr1 <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr1 , earg = .epstr1 )
    probb <- eta2theta(eta[, c(FALSE, TRUE)], .lprobb , earg = .eprobb )
    Nvec  <- extra$Nvec
    if (!is.numeric(Nvec))
      stop("something gone wrong with 'Nvec'")

    ans <-
      switch(type.fitted,
             "mean"      = pstr1 + (1 - pstr1) * Nvec *
                           probb / (1 - (1-probb)^Nvec),
             "prob"      = probb,
             "pobs1"     = doiposbinom(1, prob = probb,
                                   size = Nvec, pstr1 = pstr1), # Pr(Y=1)
             "pstr1"     =     pstr1,
             "onempstr1" = 1 - pstr1)
    label.cols.y(ans, colnames.y = extra$colnames.y, NOS = NOS)
  }, list( .lpstr1 = lpstr1, .lprobb = lprobb,
           .epstr1 = epstr1, .eprobb = eprobb
         ))),
  last = eval(substitute(expression({
    misc$link <-
      c(rep_len( .lpstr1 , NOS),
        rep_len( .lprobb , NOS))[interleave.VGAM(M, M1 = M1)]
    temp.names <- c(mynames1, mynames2)[interleave.VGAM(M, M1 = M1)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for (ii in 1:ncoly) {
      misc$earg[[M1*ii-1]] <- .epstr1
      misc$earg[[M1*ii  ]] <- .eprobb
    }
  }), list( .lpstr1 = lpstr1, .lprobb = lprobb,
            .epstr1 = epstr1, .eprobb = eprobb ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    pstr1 <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr1 , earg = .epstr1 )
    probb <- eta2theta(eta[, c(FALSE, TRUE)], .lprobb , earg = .eprobb )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * doiposbinom(x = round(extra$Nvec * y),
                                    size = extra$Nvec,  # w,
                                    pstr1 = pstr1, prob = probb,
                                    log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lpstr1 = lpstr1, .lprobb = lprobb,
           .epstr1 = epstr1, .eprobb = eprobb ))),
  vfamily = c("oiposbinomial"),



  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    pstr1 <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr1 , earg = .epstr1 )
    probb <- eta2theta(eta[, c(FALSE, TRUE)], .lprobb , earg = .eprobb )
    Nvec <- object@extra$Nvec
    roiposbinom(nsim * length(probb), size = Nvec,
                probb = probb, pstr1 = pstr1)
  }, list( .lpstr1 = lpstr1, .lprobb = lprobb,
           .epstr1 = epstr1, .eprobb = eprobb ))),

  validparams = eval(substitute(function(eta, y, extra = NULL) {
    pstr1 <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE], .lpstr1 ,
                        earg = .epstr1 )
    probb <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE], .lprobb ,
                        earg = .eprobb )
    size <- extra$Nvec
    okay1 <- all(is.finite(pstr1)) && all(pstr1 < 1) &&
             all(is.finite(probb)) && all(0 < probb & probb < 1)
    deflat.limit <- size * probb / (1 + (size-1) * probb -
                                    1 / (1-probb)^(size-1))
    okay2.deflat <- TRUE
    if (okay1 && !(okay2.deflat <- all(deflat.limit < pstr1)))
      warning("parameter 'pstr1' is too negative even allowing for ",
              "1-deflation.")
    okay1 && okay2.deflat
  }, list( .lpstr1 = lpstr1, .lprobb = lprobb,
           .epstr1 = epstr1, .eprobb = eprobb ))),






  deriv = eval(substitute(expression({
    M1 <- 2
    NOS <- M / M1
    pstr1 <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE], .lpstr1 ,
                       earg = .epstr1 )
    probb <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE], .lprobb ,
                       earg = .eprobb )

    size <- extra$Nvec

    Qn <- function(n, prob) (1 - prob)^n
    pmf1 <- size * probb * Qn(size-1, probb) / (1 - Qn(size, probb))


    onempmf1 <- 1 - pmf1  # doiposbinom(1, probb = probb, pstr1 = pstr1)
    pobs1 <- pstr1 + (1 - pstr1) * pmf1
    index1 <- as.matrix(round(w * y) == 1)

    dl.dpstr1 <- onempmf1 / pobs1
    dl.dpstr1[!index1] <- -1 / (1 - pstr1[!index1])


    d3 <- deriv3( ~ size * probb * ((1 - probb)^(size-1))
                    / (1 - (1 - probb)^size),
                  c("probb"), hessian = TRUE)
    eval.d3 <- eval(d3)
    dpmf1.dprobb   <- attr(eval.d3, "gradient")  # For checking only
    d2pmf1.dprobb2 <- attr(eval.d3, "hessian")   #
    dim(dpmf1.dprobb)   <- c(n, NOS)  # Matrix it, even for NOS==1
    dim(d2pmf1.dprobb2) <- c(n, NOS)  # Matrix it, even for NOS==1



    dl.dprobb <-  size *
      (      y  /      probb   -
       (1  - y) / (1 - probb)  -
       Qn(size-1, probb) / (1 - Qn(size, probb)))
    dl.dprobb[index1] <- (1 - pstr1[index1]) *
                         dpmf1.dprobb[index1] / pobs1[index1]

    dpstr1.deta <- dtheta.deta(pstr1, .lpstr1 , earg = .epstr1 )
    dprobb.deta <- dtheta.deta(probb, .lprobb , earg = .eprobb )

    myderiv <- cbind(dl.dpstr1 * dpstr1.deta,  # * c(w),
                     dl.dprobb * dprobb.deta)
    myderiv[, interleave.VGAM(M, M1 = M1)]
  }), list( .lpstr1 = lpstr1, .lprobb = lprobb,
            .epstr1 = epstr1, .eprobb = eprobb ))),




  weight = eval(substitute(expression({


    d4 <- deriv3( ~ size * ((1 - probb)^(size-1)) /
                    (1 - (1 - probb)^size),
                 c("probb"), hessian = FALSE)
    eval.d4 <- eval(d4)
    d2logonempmf0.dprobb2 <- attr(eval.d4, "gradient")
    dim(d2logonempmf0.dprobb2) <- c(n, NOS)  # Matrix it, even for NOS==1


    E2 <- function(size, prob) {
      size *
      prob * (1 - Qn(size-1, prob)) /
     (1 - Qn(size, prob) - size * prob * Qn(size-1, prob))
    }

    E2mat <- E2(size, probb)
    RHS <- onempmf1 * (        E2mat  /    probb^2 +
                       (size - E2mat) / (1-probb)^2 +
                       d2logonempmf0.dprobb2)


    LHS <- -d2pmf1.dprobb2 + ((1-pstr1) / pobs1) * dpmf1.dprobb^2


    ned2l.dpstr12 <- onempmf1 / ((1 - pstr1) * pobs1)
    ned2l.dpstr1probb <- dpmf1.dprobb / pobs1
    ned2l.dprobb2 <- (1 - pstr1) * (LHS + RHS)


    wz <- array(c(ned2l.dpstr12 * dpstr1.deta^2,
                  ned2l.dprobb2 * dprobb.deta^2,
                  ned2l.dpstr1probb * dpstr1.deta * dprobb.deta),
                dim = c(n, M / M1, 3))
    wz <- arwz2wz(wz, M = M, M1 = M1)
    wz
  }), list( .lprobb = lprobb, .eprobb = eprobb ))))
}  # oiposbinomial











deflat.limit.oilog  <- function(shape) {
  if (any(shape <= 0 | 1 <= shape ))
    stop("argument 'shape' must be in (0, 1)")
  ans <- 1 / (1 - 1 / dlog(1, shape))
  ans
}



doilog <- function(x, shape, pstr1 = 0, log = FALSE) {



  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(shape), length(pstr1))
  if (length(x)     != LLL) x     <- rep_len(x,     LLL)
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  if (length(pstr1) != LLL) pstr1 <- rep_len(pstr1, LLL)

  ans <- rep(NA_real_, LLL)
  index1 <- (x == 1)
  if (log.arg) {
    ans[ index1] <- log(pstr1[ index1] + (1 - pstr1[ index1]) *
                       dlog(x[ index1], shape[ index1]))
    ans[!index1] <- log1p(-pstr1[!index1]) +
                          dlog(x[!index1], shape[!index1], log = TRUE)
  } else {
    ans[ index1] <-      pstr1[ index1] + (1 - pstr1[ index1]) *
                        dlog(x[ index1], shape[ index1])
    ans[!index1] <- (1 - pstr1[!index1]) *
                        dlog(x[!index1], shape[!index1])
  }


  ans[pstr1 < deflat.limit.oilog(shape) | 1 < pstr1] <- NaN
  ans[shape <= 0 | 1 <= shape] <- NaN
  ans
}  # doilog




poilog <- function(q, shape, pstr1 = 0) {

  LLL <- max(length(q), length(shape), length(pstr1))
  if (length(q)     != LLL) q     <- rep_len(q,     LLL)
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  if (length(pstr1) != LLL) pstr1 <- rep_len(pstr1, LLL)
  ans <- rep_len(NA_real_, LLL)
  deflat.limit <- deflat.limit.oilog(shape)

  ans <- plog(q, shape)  #, lower.tail = lower.tail, log.p = log.p
  ans <- ifelse(q < 1, 0, pstr1 + (1 - pstr1) * ans)

  ans[pstr1 < deflat.limit] <- NaN
  ans[1 < pstr1] <- NaN
  ans[shape <= 0] <- NaN
  ans[1 <= shape] <- NaN

  ans
}  # poilog





qoilog <- function(p, shape, pstr1 = 0) {

  LLL <- max(length(p), length(shape), length(pstr1))
  if (length(p)     != LLL) p     <- rep_len(p,     LLL)
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  if (length(pstr1) != LLL) pstr1 <- rep_len(pstr1, LLL)
  ans <- rep_len(NA_real_, LLL)
  deflat.limit <- deflat.limit.oilog(shape)

  ans[p <= pstr1] <- 1
  pindex <- (deflat.limit <= pstr1) & (pstr1 < p)
  ans[pindex] <-
    qlog((p[pindex] - pstr1[pindex]) / (1 - pstr1[pindex]),
         shape = shape[pindex])

  ans[pstr1 < deflat.limit] <- NaN
  ans[1 < pstr1] <- NaN

  ans[p < 0] <- NaN
  ans[1 < p] <- NaN
  ans[shape <= 0] <- NaN
  ans[1 <= shape] <- NaN

  ans
}  # qoilog



roilog <- function(n, shape, pstr1 = 0) {
  qoilog(runif(n), shape, pstr1 = pstr1)
}





 oilog <-
  function(lpstr1 = "logitlink", lshape = "logitlink",
       type.fitted = c("mean", "shape", "pobs1", "pstr1", "onempstr1"),
           ishape = NULL,
           gpstr1 = ppoints(8),
           gshape = ppoints(8),
           zero = NULL) {

  lpstr1 <- as.list(substitute(lpstr1))
  epstr1 <- link2list(lpstr1)
  lpstr1 <- attr(epstr1, "function.name")

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")



  type.fitted <- match.arg(type.fitted,
                 c("mean", "shape", "pobs1", "pstr1", "onempstr1"))[1]


  if (length(ishape))
    if (!is.Numeric(ishape, positive = TRUE))
      stop("argument 'ishape' values must be positive")


  new("vglmff",
  blurb = c("One-inflated logarithmic distribution\n\n",
            "Links:    ",
            namesof("pstr1", lpstr1, earg = epstr1 ), ", ",
            namesof("shape", lshape, earg = eshape ), "\n",
            "Mean:     pstr1 + (1 - pstr1) * a * shape / (1 - shape), ",
                       "a = -1 / log(1-shape), 0 < shape < 1"),

  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = TRUE,
         parameters.names = c("pstr1", "shape"),
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted
         ))),
  initialize = eval(substitute(expression({
    M1 <- 2

    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    NOS <- ncoly <- ncol(y)
    extra$ncoly <- ncoly
    M <- M1 * ncoly
    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)

    mynames1 <- param.names("pstr1", ncoly, skip1 = TRUE)
    mynames2 <- param.names("shape", ncoly, skip1 = TRUE)
    predictors.names <-
        c(namesof(mynames1, .lpstr1 , earg = .epstr1 , tag = FALSE),
          namesof(mynames2, .lshape , earg = .eshape , tag = FALSE))[
          interleave.VGAM(M, M1 = M1)]


    if (!length(etastart)) {

      shape.init <-
      pstr1.init <- matrix(NA_real_, n, NOS)
      gpstr1 <- .gpstr1
      gshape <- .gshape

      oilog.Loglikfun <- function(pstr1, shape, y, x, w, extraargs) {
        sum(c(w) * doilog(x = y, pstr1 = pstr1,
                          shape = shape, log = TRUE))
      }


      for (jay in 1:NOS) {  # For each response 'y_jay'... do:


        try.this <-
          grid.search2(gpstr1, gshape,
                       objfun = oilog.Loglikfun,
                       y = y[, jay],  # x = x[TFvec, , drop = FALSE],
                       w = w[, jay],
                       ret.objfun = TRUE)  # Last value is the loglik
        pstr1.init[, jay] <-  try.this["Value1"]
        shape.init[, jay] <-  try.this["Value2"]
      }  # for (jay ...)

      etastart <- cbind(theta2eta(pstr1.init, .lpstr1 , earg = .epstr1 ),
                        theta2eta(shape.init, .lshape , earg = .eshape ))[,
                        interleave.VGAM(M, M1 = M1)]
      mustart <- NULL  # Since etastart has been computed.
    }  # End of !length(etastart)
  }), list( .lpstr1 = lpstr1, .lshape = lshape,
            .epstr1 = epstr1, .eshape = eshape,
                              .ishape = ishape,
            .gpstr1 = gpstr1,
            .gshape  = gshape,
            .type.fitted = type.fitted ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    NOS <- ncol(eta) / c(M1 = 2)
    type.fitted <- if (length(extra$type.fitted))
                     extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                   c("mean", "shape", "pobs1", "pstr1", "onempstr1"))[1]

    pstr1 <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr1 , earg = .epstr1 )
    shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , earg = .eshape )

    Meanfun <- function(shape) {
      aa <- -1 / log1p(-shape)
      Mean <- aa * shape / (1 - shape)
      Mean[shape <= 0 | 1 <= shape] <- NaN
      Mean
    }

    ans <-
      switch(type.fitted,
             "mean"      = pstr1 + (1 - pstr1) * Meanfun(shape),
             "shape"     = shape,
           "pobs1" = doizeta(1, shape = shape, pstr1 = pstr1),  # P(Y=1)
             "pstr1"     =     pstr1,
             "onempstr1" = 1 - pstr1)
    label.cols.y(ans, colnames.y = extra$colnames.y, NOS = NOS)
  }, list( .lpstr1 = lpstr1, .lshape = lshape,
           .epstr1 = epstr1, .eshape = eshape ))),
  last = eval(substitute(expression({
    misc$link <-
      c(rep_len( .lpstr1 , NOS),
        rep_len( .lshape , NOS))[interleave.VGAM(M, M1 = M1)]
    temp.names <- c(mynames1, mynames2)[interleave.VGAM(M, M1 = M1)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for (ii in 1:ncoly) {
      misc$earg[[M1*ii-1]] <- .epstr1
      misc$earg[[M1*ii  ]] <- .eshape
    }
  }), list( .lpstr1 = lpstr1, .lshape = lshape,
            .epstr1 = epstr1, .eshape = eshape ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    pstr1 <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr1 , earg = .epstr1 )
    shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , earg = .eshape )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * doilog(x = y, pstr1 = pstr1, shape = shape,
                               log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lpstr1 = lpstr1, .lshape = lshape,
           .epstr1 = epstr1, .eshape = eshape ))),
  vfamily = c("oilog"),



  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    pstr1 <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr1 , earg = .epstr1 )
    shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , earg = .eshape )
    roilog(nsim * length(shape), shape = shape, pstr1 = pstr1)
  }, list( .lpstr1 = lpstr1, .lshape = lshape,
           .epstr1 = epstr1, .eshape = eshape ))),

  validparams = eval(substitute(function(eta, y, extra = NULL) {
    pstr1 <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE], .lpstr1 ,
                       earg = .epstr1 )
    shape <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE], .lshape  ,
                       earg = .eshape )
    okay1 <- all(is.finite(shape)) && all(0 < shape & shape < 1) &&
             all(is.finite(pstr1)) && all(pstr1 < 1)
    deflat.limit <- deflat.limit.oizeta(shape)
    okay2.deflat <- TRUE
    if (okay1 && !(okay2.deflat <- all(deflat.limit < pstr1)))
      warning("parameter 'pstr1' is too negative even allowing for ",
              "1-deflation.")
    okay1 && okay2.deflat
  }, list( .lpstr1 = lpstr1, .lshape = lshape,
           .epstr1 = epstr1, .eshape = eshape ))),





  deriv = eval(substitute(expression({
    M1 <- 2
    NOS <- M / M1
    pstr1 <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE], .lpstr1 ,
                       earg = .epstr1 )
    shape <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE], .lshape ,
                       earg = .eshape )

    pmf1 <- dlog(1, shape)
    onempmf1 <- 1 - pmf1  # dozeta(1, shape = shape, pstr1 = pstr1)
    pobs1 <- pstr1 + (1 - pstr1) * pmf1
    index1 <- as.matrix(y == 1)

    mraa <- log1p(-shape)
    aaaa <- -1 / mraa

    dl.dpstr1 <- onempmf1 / pobs1
    dl.dpstr1[!index1] <- -1 / (1 - pstr1[!index1])




    dpmf1.dshape <- -1 / mraa - shape / ((1 - shape) * mraa^2)
    d2pmf1.dshape2 <- -2 / ((1 - shape) * mraa^2) -
                      shape * (2 + mraa) / ((1 - shape)^2 * mraa^3)

    dl.dshape <- (1 - pstr1) * dpmf1.dshape / pobs1  #
    dl.dshape[!index1] <- y[!index1] / shape[!index1] +
                         1 / ((1 - shape[!index1]) * mraa[!index1])

    dpstr1.deta <- dtheta.deta(pstr1, .lpstr1 , earg = .epstr1 )
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )

    myderiv <- c(w) * cbind(dl.dpstr1 * dpstr1.deta,
                            dl.dshape * dshape.deta)
    myderiv[, interleave.VGAM(M, M1 = M1)]
  }), list( .lpstr1 = lpstr1, .lshape = lshape,
            .epstr1 = epstr1, .eshape = eshape ))),
  weight = eval(substitute(expression({

    EY.y.gt.1 <- aaaa * shape^2 / ((1 - shape) * (1 - aaaa * shape))
    LHS <- ((1 - pstr1) / pobs1) * dpmf1.dshape^2 - d2pmf1.dshape2
    RHS <- EY.y.gt.1 / shape^2 - (1 + mraa) / ((1 - shape) * mraa)^2
    ned2l.dpstr12 <- onempmf1 / ((1 - pstr1) * pobs1)  #
    ned2l.dpstr1shape <- dpmf1.dshape / pobs1  #
    ned2l.dshape2 <- (1 - pstr1) * (LHS + (1 - pmf1) * RHS)

    wz <- array(c(c(w) * ned2l.dpstr12 * dpstr1.deta^2,
                  c(w) * ned2l.dshape2 * dshape.deta^2,
                  c(w) * ned2l.dpstr1shape * dpstr1.deta * dshape.deta),
                dim = c(n, M / M1, 3))
    wz <- arwz2wz(wz, M = M, M1 = M1)
    wz
  }), list( .lshape = lshape, .eshape = eshape ))))
}  # oilog









dotlog <- function(x, shape, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  if (log.arg) {
    ans <- dlog(x, shape, log = log.arg) - log1p(-dlog(1, shape))
    ans[x == 1] <- log(0)
  } else {
    ans <- dlog(x, shape) / (1 - dlog(1, shape))
    ans[x == 1] <- 0
  }
  ans
}  # dotlog



potlog  <- function(q, shape, log.p = FALSE) {
  ans <- if (log.p) log(plog(q, shape) - dlog(1, shape)) -
      log1p(-dlog(1, shape)) else
    (plog(q, shape) - dlog(1, shape)) / (1 - dlog(1, shape))
  ans[q < 1] <- if (log.p) log(0) else 0
  ans
}






 qotlog <- function(p, shape) {

  ans <- qlog((1 - dlog(1, shape)) * p + dlog(1, shape), shape = shape)

  ans[p == 1] <- Inf
  ans[p < 0] <- NaN
  ans[1 < p] <- NaN

  ans[shape < 0 | 1 < shape] <- NaN
  ans
}  # qotlog



rotlog <- function(n, shape) {
  qotlog(runif(n), shape)
}



 otlog <-
  function(lshape = "logitlink", gshape = ppoints(8), zero = NULL) {

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")



  new("vglmff",
  blurb = c("One-truncated logarithmic distribution ",
            "f(y) = shape^y / ((-shape - log1p(-shape)) * y), ",
             "y = 2, 3,...,\n",
             "            0 < shape < 1,\n\n",
             "Link:    ", namesof("shape", lshape, earg = eshape)),
  constraints = eval(substitute(expression({
    dotzero <- .zero
    M1 <- 1
    eval(negzero.expression.VGAM)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = TRUE,
         parameters.names = "shape",
         zero = .zero )
  }, list( .zero = zero ))),


  initialize = eval(substitute(expression({
    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
    if (any(y <= 1))
      stop("cannot have any 1s in the response")


    ncoly <- ncol(y)
    M1 <- 1
    extra$ncoly <- ncoly
    M <- M1 * ncoly


    mynames1  <- param.names("shape", ncoly, skip1 = TRUE)
    predictors.names <- namesof(mynames1, .lshape , earg = .eshape ,
                                tag = FALSE)


    if (!length(etastart)) {
      dotlog.Loglikfun <- function(shapeval, y, x, w, extraargs) {
        sum(c(w) * dotlog(x = y, shape = shapeval, log = TRUE))
      }
      Init.shape <- matrix(0, n, M)
      shape.grid <- .gshape

      for (ilocal in 1:ncoly) {
        Init.shape[, ilocal] <- grid.search(shape.grid,
                                            objfun = dotlog.Loglikfun,
                                            y = y[, ilocal],  # x = x,
                                            w = w[, ilocal])
      }  # for
      etastart <- theta2eta(Init.shape, .lshape , earg = .eshape )
    }
  }), list( .lshape = lshape, .eshape = eshape, .gshape = gshape ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    aa <- -1 / log1p(-shape)
    ((aa * shape / (1 - shape)) - dlog(1, shape)) / (1 - dlog(1, shape))
  }, list( .lshape = lshape, .eshape = eshape ))),

  last = eval(substitute(expression({
    misc$link <- c(rep_len( .lshape , ncoly))
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for (ii in 1:ncoly) {
      misc$earg[[ii]] <- .eshape
    }
  }), list( .lshape = lshape, .eshape = eshape ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dotlog(x = y, shape = shape, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lshape = lshape, .eshape = eshape ))),
  vfamily = c("otlog"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    okay1 <- all(is.finite(shape)) && all(0 < shape & shape < 1)
    okay1
  }, list( .lshape = lshape, .eshape = eshape ))),


  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    rotlog(nsim * length(shape), shape = shape)
  }, list( .lshape = lshape, .eshape = eshape ))),



  deriv = eval(substitute(expression({
    M1 <- 1
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    aa <- -1 / log1p(-shape)
    dl.dshape <- y / shape +
                shape / ((1 - shape) * (shape + log1p(-shape)))
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )
    c(w) * dl.dshape * dshape.deta
  }), list( .lshape = lshape, .eshape = eshape ))),
  weight = eval(substitute(expression({
    EY.logff <-  aa * shape / (1 - shape)

    d3 <- deriv3( ~ shape / ((1 - shape) * (shape + log(1 - shape))),
                  c("shape"), hessian = FALSE)
    eval.d3 <- eval(d3)
    d2pmf1.dshape2 <- c(attr(eval.d3, "gradient"))

    ned2l.dshape2 <-
      (EY.logff - dlog(1, shape)) / ((1 - dlog(1, shape)) * shape^2) -
      d2pmf1.dshape2
    wz <- c(w) * ned2l.dshape2 * dshape.deta^2
    wz
  }), list( .lshape = lshape, .eshape = eshape ))))
}  # otlog






dotpospois <- function(x, lambda, log = FALSE) {
  if (!is.logical(larg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  if (larg) {
    ans <- dpospois(x, lambda, log = larg) - log1p(-dpospois(1, lambda))
    ans[x == 1] <- log(0)
  } else {
    ans <- dpospois(x, lambda) / (1 - dpospois(1, lambda))
    ans[x == 1] <- 0
  }
  ans
}  # dotpospois



potpospois  <- function(q, lambda, log.p = FALSE) {
  if (log.p) log(ppospois(q, lambda) - dpospois(1, lambda)) -
      log1p(-dpospois(1, lambda)) else
    (ppospois(q, lambda) - dpospois(1, lambda)) / (1-dpospois(1, lambda))
}



 qotpospois <- function(p, lambda) {
  ans <- qpospois((1 - dpospois(1, lambda)) * p +
                  dpospois(1, lambda), lambda = lambda)

  ans[p == 1 & 0 < lambda] <- Inf
  ans[p < 0] <- NaN
  ans[1 < p] <- NaN

  ans[lambda < 0] <- NaN
  ans
}  # qotpospois



rotpospois <- function(n, lambda) {
  qotpospois(runif(n), lambda)
}




 otpospoisson <-
    function(llambda = "loglink",
             type.fitted = c("mean", "lambda", "prob0", "prob1"),
             ilambda = NULL, imethod = 1, zero = NULL) {

  llambda <- as.list(substitute(llambda))
  elambda <- link2list(llambda)
  llambda <- attr(elambda, "function.name")


  if (length( ilambda) && !is.Numeric(ilambda, positive = TRUE))
    stop("bad input for argument 'ilambda'")

  type.fitted <- match.arg(type.fitted,
                           c("mean", "lambda", "prob0", "prob1"))[1]


  new("vglmff",
  blurb = c("One-truncated Positive-Poisson distribution\n\n",
            "Links:    ",
            namesof("lambda", llambda, earg = elambda, tag = FALSE)),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 1)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = TRUE,
         parameters.names = c("lambda"),
         type.fitted  = .type.fitted ,
         llambda = .llambda ,
         elambda = .elambda )
  }, list( .llambda = llambda, .elambda = elambda,
           .type.fitted = type.fitted ))),

  initialize = eval(substitute(expression({
    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
    if (any(y < 2))
      stop("response values must be 2 or more")

    ncoly <- ncol(y)
    M1 <- 1
    extra$ncoly <- ncoly
    M <- M1 * ncoly
    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)


    mynames1 <- param.names("lambda", ncoly, skip1 = TRUE)
    predictors.names <- namesof(mynames1, .llambda , earg = .elambda ,
                                tag = FALSE)

    if (!length(etastart)) {
      lambda.init <- Init.mu(y = y, w = w, imethod = .imethod ,
                             imu = .ilambda )

      etastart <- theta2eta(lambda.init, .llambda , earg = .elambda)
    }
  }), list( .llambda = llambda, .elambda = elambda,
            .ilambda = ilambda, .imethod = imethod,
            .type.fitted = type.fitted ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    NOS <- NCOL(eta) / c(M1 = 1)
   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "lambda", "prob0", "prob1"))[1]

    lambda <- eta2theta(eta, .llambda , earg = .elambda )
    ans <- switch(type.fitted,
                  "mean"      = lambda / ppois(1, lambda, lower = FALSE),
                  "lambda"    = lambda,
                  "prob0"     = ppois(0, lambda),  # P(Y=0) as it were
                  "prob1"     = ppois(1, lambda))  # P(Y=1) as it were
    label.cols.y(ans, colnames.y = extra$colnames.y, NOS = NOS)
  }, list( .llambda = llambda, .elambda = elambda ))),
  last = eval(substitute(expression({
    misc$link <- rep_len( .llambda , M)
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for (ii in 1:M)
      misc$earg[[ii]] <- .elambda
  }), list( .llambda = llambda, .elambda = elambda ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    lambda <- eta2theta(eta, .llambda , earg = .elambda )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dotpospois(x = y, lambda = lambda, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .llambda = llambda, .elambda = elambda ))),
  vfamily = c("otpospoisson"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    lambda <- eta2theta(eta, .llambda , earg = .elambda )
    okay1 <- all(is.finite(lambda)) && all(0 < lambda)
    okay1
  }, list( .llambda = llambda, .elambda = elambda ))),


  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    lambda <- eta2theta(eta, .llambda , earg = .elambda )
    rotpospois(nsim * length(lambda), lambda)
  }, list( .llambda = llambda, .elambda = elambda ))),




  deriv = eval(substitute(expression({
    M1 <- 1
    lambda <- eta2theta(eta, .llambda , earg = .elambda )

    EY.cond <- 1 / ppois(1, lambda, lower.tail = FALSE)
    temp1 <- expm1(lambda)
    temp0 <- lambda * exp(-lambda)
    prob.geq.2 <- -expm1(-lambda) - temp0
    dl.dlambda <- y / lambda - 1 - temp0 / prob.geq.2

    dlambda.deta <- dtheta.deta(lambda, .llambda , earg = .elambda )

    c(w) * dl.dlambda * dlambda.deta
  }), list( .llambda = llambda, .elambda = elambda ))),
  weight = eval(substitute(expression({
    ned2l.dlambda2 <- EY.cond / lambda +
        ((1 - lambda) * exp(-lambda) - temp0^2 / prob.geq.2) / prob.geq.2
    wz <-  ned2l.dlambda2 * dlambda.deta^2
    c(w) * wz
  }), list( .llambda = llambda, .elambda = elambda ))))
}  # otpospoisson









doalog <- function(x, shape, pobs1 = 0, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(shape), length(pobs1))
  if (length(x)      != LLL) x      <- rep_len(x,      LLL)
  if (length(shape)  != LLL) shape  <- rep_len(shape,  LLL)
  if (length(pobs1)  != LLL) pobs1  <- rep_len(pobs1,  LLL)
  ans <- rep_len(0.0, LLL)


  index1 <- (x == 1)

  if (log.arg) {
    ans[ index1] <- log(pobs1[index1])
    ans[!index1] <- log1p(-pobs1[!index1]) +
                    dotlog(x[!index1], shape[!index1], log = TRUE)
  } else {
    ans[ index1] <- pobs1[index1]
    ans[!index1] <- (1 - pobs1[!index1]) *
                    dotlog(x[!index1], shape[!index1])
  }
  ans[pobs1 < 0 | 1 < pobs1] <- NaN
  ans
}



poalog <- function(q, shape, pobs1 = 0) {
  LLL <- max(length(q), length(shape), length(pobs1))
  if (length(q)      != LLL) q      <- rep_len(q,      LLL)
  if (length(shape)  != LLL) shape  <- rep_len(shape,  LLL)
  if (length(pobs1)  != LLL) pobs1  <- rep_len(pobs1,  LLL)
  ans <- rep_len(0.0, LLL)

  ans[q >  1] <-    pobs1[q > 1] +
                 (1-pobs1[q > 1]) * potlog(q[q > 1], shape[q > 1])
  ans[q <  1] <- 0
  ans[q == 1] <- pobs1[q == 1]

  ans <- pmax(0, ans)
  ans <- pmin(1, ans)
  ans[pobs1 < 0 | 1 < pobs1] <- NaN

  ans
}



qoalog <- function(p, shape, pobs1 = 0) {
  LLL <- max(length(p), length(shape), length(pobs1))
  if (length(p)      != LLL) p      <- rep_len(p,      LLL)
  if (length(shape)  != LLL) shape  <- rep_len(shape,  LLL)
  if (length(pobs1)  != LLL) pobs1  <- rep_len(pobs1,  LLL)

  ans <- rep_len(NaN, LLL)
  ind4 <- pobs1 < p
  ans[!ind4] <- 1
  ans[ ind4] <- qotlog((p[ind4] - pobs1[ind4]) / (1 - pobs1[ind4]),
                       shape = shape[ind4])
  ans[pobs1 < 0 | 1 < pobs1] <- NaN
  ans[p < 0 | 1 < p] <- NaN
  ans
}



roalog <- function(n, shape, pobs1 = 0) {
  qoalog(runif(n), shape = shape, pobs1 = pobs1)
}






 oalog <-
  function(lpobs1 = "logitlink",
           lshape = "logitlink",
           type.fitted = c("mean", "shape", "pobs1", "onempobs1"),
           ipobs1 = NULL,
           gshape = ppoints(8),
           zero = NULL) {


  lpobs1 <- as.list(substitute(lpobs1))
  epobs1 <- link2list(lpobs1)
  lpobs1 <- attr(epobs1, "function.name")

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")


 type.fitted <- match.arg(type.fitted,
                           c("mean", "shape", "pobs1", "onempobs1"))[1]


  new("vglmff",
  blurb = c("One-altered logarithmic distribution \n",
            "(Bernoulli and 1-truncated logarithmic distribution model)",
            "\n\n",
            "Links:    ",
            namesof("pobs1",  lpobs1, earg = epobs1, tag = FALSE), ", ",
            namesof("shape",  lshape, earg = eshape, tag = FALSE)),
  constraints = eval(substitute(expression({
    dotzero <- .zero
    M1 <- 2
    eval(negzero.expression.VGAM)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = TRUE,
         parameters.names = c("pobs1", "shape"),
         zero = .zero )
  }, list( .zero = zero ))),


  initialize = eval(substitute(expression({
    M1 <- 2
    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y

    extra$y1 <- y1 <- ifelse(y == 1, 1, 0)
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$skip.these <- skip.these <- matrix(as.logical(y1), n, NOS)
    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)

    mynames1 <- param.names("pobs1", ncoly, skip1 = TRUE)
    mynames2 <- param.names("shape", ncoly, skip1 = TRUE)
    predictors.names <-
        c(namesof(mynames1, .lpobs1 , earg = .epobs1 , tag = FALSE),
          namesof(mynames2, .lshape , earg = .eshape , tag = FALSE))[
          interleave.VGAM(M1*NOS, M1 = M1)]

    ncoly <- ncol(y)
    extra$ncoly <- ncoly
    M <- M1 * ncoly


    if (!length(etastart)) {
      dotlog.Loglikfun <- function(shapeval, y, x, w, extraargs) {
        sum(c(w) * dotlog(x = y, shape = shapeval, log = TRUE))
      }
      Init.shape <- matrix(0, n, ncoly)
      shape.grid <- .gshape

      for (jlocal in 1:ncoly) {
        index1 <- y[, jlocal] > 1
        Init.shape[, jlocal] <-
          grid.search(shape.grid,
                      objfun = dotlog.Loglikfun,
                      y = y[index1, jlocal],  # x = x,
                      w = w[index1, jlocal])
      }  # for
      etastart <-
        cbind(theta2eta(if (length( .ipobs1 )) .ipobs1 else
                        (0.5 + w * y1) / (1 + w),
                        .lpobs1 , earg = .epobs1 ),
              theta2eta(Init.shape, .lshape , earg = .eshape ))
      etastart <- etastart[, interleave.VGAM(ncol(etastart), M1 = M1)]
    }
  }), list( .lshape = lshape, .eshape = eshape, .gshape = gshape,
            .lpobs1 = lpobs1, .epobs1 = epobs1,
            .ipobs1 = ipobs1,  # .ishape = ishape,
            .type.fitted = type.fitted
           ))),


  linkinv = eval(substitute(function(eta, extra = NULL) {
   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "shape", "pobs1", "onempobs1"))[1]

    M1 <- 2
    NOS <- ncol(eta) / M1

    pobs1 <- cbind(eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                             .lpobs1 , earg = .epobs1 ))
    shape <- cbind(eta2theta(eta[, M1*(1:NOS)-0, drop = FALSE],
                             .lshape , earg = .eshape ))

    aa <- -1 / log1p(-shape)
    otlog.mean <- ((aa * shape / (1 - shape)) -
                   dlog(1, shape)) / (1 - dlog(1, shape))

    ans <- switch(type.fitted,
                  "mean"      = pobs1 + (1 - pobs1) * otlog.mean,
                  "shape"     = shape,
                  "pobs1"     =      pobs1,  # P(Y=1)
                  "onempobs1" =  1 - pobs1)  # P(Y>1)
    label.cols.y(ans, colnames.y = extra$colnames.y, NOS = NOS)
  }, list( .lpobs1 = lpobs1, .lshape = lshape,
           .epobs1 = epobs1, .eshape = eshape ))),
  last = eval(substitute(expression({
    temp.names <- c(rep_len( .lpobs1 , NOS),
                    rep_len( .lshape , NOS))
    temp.names <- temp.names[interleave.VGAM(M1*NOS, M1 = M1)]
    misc$link  <- temp.names
    names(misc$link) <-
      c(mynames1, mynames2)[interleave.VGAM(M1*NOS, M1 = M1)]

    misc$earg <- vector("list", M1 * NOS)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:NOS) {
      misc$earg[[M1*ii-1]] <- .epobs1
      misc$earg[[M1*ii  ]] <- .eshape
    }
  }), list( .lpobs1 = lpobs1, .lshape = lshape,
            .epobs1 = epobs1, .eshape = eshape ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {

    pobs1 <- cbind(eta2theta(eta[, c(TRUE, FALSE), drop = FALSE],
                             .lpobs1, earg = .epobs1 ))
    shape <- cbind(eta2theta(eta[, c(FALSE, TRUE), drop = FALSE],
                             .lshape, earg = .eshape ))

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * doalog(x = y, pobs1 = pobs1, shape = shape,
                               log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lpobs1 = lpobs1, .lshape = lshape,
           .epobs1 = epobs1, .eshape = eshape ))),
  vfamily = c("oalog"),



  validparams = eval(substitute(function(eta, y, extra = NULL) {
    TFvec <- c(TRUE, FALSE)
    pobs1 <- eta2theta(eta[,  TFvec, drop = FALSE],
                       .lpobs1 , earg = .epobs1 )
    shape <- eta2theta(eta[, !TFvec, drop = FALSE],
                       .lshape , earg = .eshape )
    okay1 <- all(is.finite(shape)) && all(0 < shape & shape < 1) &&
             all(is.finite(pobs1)) && all(0 < pobs1 & pobs1 < 1)
    okay1
  }, list( .lpobs1 = lpobs1, .lshape = lshape,
           .epobs1 = epobs1, .eshape = eshape ))),


  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    pobs1 <- eta2theta(eta[, c(TRUE, FALSE)], .lpobs1 , earg = .epobs1 )
    shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , earg = .eshape )
    roalog(nsim * length(shape), shape = shape, pobs1 = pobs1)
  }, list( .lpobs1 = lpobs1, .lshape = lshape,
           .epobs1 = epobs1, .eshape = eshape ))),



  deriv = eval(substitute(expression({
    M1 <- 2
    NOS <- ncol(eta) / M1  # extra$NOS
    y1 <- extra$y1
    skip <- extra$skip.these

    TFvec <- c(TRUE, FALSE)
    pobs1 <- eta2theta(eta[,  TFvec, drop = FALSE],
                       .lpobs1 , earg = .epobs1 )
    shape <- eta2theta(eta[, !TFvec, drop = FALSE],
                       .lshape , earg = .eshape )

    aa <- -1 / log1p(-shape)
    dl.dshape <- y / shape +
                 shape / ((1 - shape) * (shape + log1p(-shape)))

    dl.dpobs1 <- -1 / (1 - pobs1)  # For y > 1 obsns

    for (spp. in 1:NOS) {
      dl.dpobs1[skip[, spp.], spp.] <- 1 / pobs1[skip[, spp.], spp.]
      dl.dshape[skip[, spp.], spp.] <- 0
    }
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )
    mu.phi1 <- pobs1

    temp3 <- if ( .lpobs1 == "logitlink") {
      c(w) * (y1 - mu.phi1)
    } else {
      c(w) * dtheta.deta(mu.phi1, link = .lpobs1 , earg = .epobs1 ) *
             dl.dpobs1
    }

    ans <- cbind(temp3,
                 c(w) * dl.dshape * dshape.deta)
    ans <- ans[, interleave.VGAM(ncol(ans), M1 = M1)]
    ans
  }), list( .lpobs1 = lpobs1, .lshape = lshape,
            .epobs1 = epobs1, .eshape = eshape ))),
  weight = eval(substitute(expression({
    wz <- matrix(0, n, M1 * NOS)  # EIM is diagonal


    EY.logff <-  aa * shape / (1 - shape)
    d3 <- deriv3( ~ shape / ((1 - shape) * (shape + log(1 - shape))),
                  c("shape"), hessian = FALSE)
    eval.d3 <- eval(d3)
    d2pmf1.dshape2 <- c(attr(eval.d3, "gradient"))

    ned2l.dshape2 <-
      (EY.logff - dlog(1, shape)) / ((1 - dlog(1, shape)) * shape^2) -
      d2pmf1.dshape2



    ned2l.dshape2 <- (1-pobs1) * ned2l.dshape2  #+stop("another quantity")
    wz[, NOS+(1:NOS)] <- c(w) * ned2l.dshape2 * dshape.deta^2


    tmp100 <- mu.phi1 * (1 - mu.phi1)
    tmp200 <- if ( .lpobs1 == "logitlink" && is.empty.list( .epobs1 )) {
        cbind(c(w) * tmp100)
    } else {
      cbind(c(w) * (1 / tmp100) *
            dtheta.deta(mu.phi1, link = .lpobs1 , earg = .epobs1 )^2)
    }
    wz[, 1:NOS] <-  tmp200

    wz <- wz[, interleave.VGAM(ncol(wz), M1 = M1)]

    wz
  }), list( .lpobs1 = lpobs1,
            .epobs1 = epobs1 ))))
}  # End of oalog








doapospois <- function(x, lambda, pobs1 = 0, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(lambda), length(pobs1))
  if (length(x)      != LLL) x      <- rep_len(x,      LLL)
  if (length(lambda) != LLL) lambda <- rep_len(lambda, LLL)
  if (length(pobs1)  != LLL) pobs1  <- rep_len(pobs1,  LLL)
  ans <- rep_len(0.0, LLL)


  index1 <- (x == 1)

  if (log.arg) {
    ans[ index1] <- log(pobs1[index1])
    ans[!index1] <- log1p(-pobs1[!index1]) +
                    dotpospois(x[!index1], lambda[!index1], log = TRUE)
  } else {
    ans[ index1] <- pobs1[index1]
    ans[!index1] <- (1 - pobs1[!index1]) *
                    dotpospois(x[!index1], lambda[!index1])
  }
  ans[pobs1 < 0 | 1 < pobs1] <- NaN
  ans[lambda < 0] <- NaN
  ans
}



poapospois <- function(q, lambda, pobs1 = 0) {
  LLL <- max(length(q), length(lambda), length(pobs1))
  if (length(q)      != LLL) q      <- rep_len(q,      LLL)
  if (length(lambda) != LLL) lambda <- rep_len(lambda, LLL)
  if (length(pobs1)  != LLL) pobs1  <- rep_len(pobs1,  LLL)
  ans <- rep_len(0.0, LLL)

  ans[q >  1] <-    pobs1[q > 1] +
                 (1-pobs1[q > 1]) * potpospois(q[q > 1], lambda[q > 1])
  ans[q <  1] <- 0
  ans[q == 1] <- pobs1[q == 1]

  ans <- pmax(0, ans)
  ans <- pmin(1, ans)
  ans[pobs1 < 0 | 1 < pobs1] <- NaN
  ans[lambda < 0] <- NaN

  ans
}



qoapospois <- function(p, lambda, pobs1 = 0) {
  LLL <- max(length(p), length(lambda), length(pobs1))
  if (length(p)      != LLL) p      <- rep_len(p,      LLL)
  if (length(lambda) != LLL) lambda <- rep_len(lambda, LLL)
  if (length(pobs1)  != LLL) pobs1  <- rep_len(pobs1,  LLL)

  ans <- rep_len(NaN, LLL)
  ind4 <- pobs1 < p
  ans[!ind4] <- 1
  ans[ ind4] <- qotpospois((p[ind4] - pobs1[ind4]) / (1 - pobs1[ind4]),
                           lambda = lambda[ind4])
  ans[pobs1 < 0 | 1 < pobs1] <- NaN
  ans[p < 0 | 1 < p] <- NaN
  ans[lambda < 0] <- NaN
  ans
}



roapospois <- function(n, lambda, pobs1 = 0) {
  qoapospois(runif(n), lambda = lambda, pobs1 = pobs1)
}






 oapospoisson <-
  function(lpobs1 = "logitlink",
           llambda = "loglink",
           type.fitted = c("mean", "lambda", "pobs1", "onempobs1"),
           ipobs1 = NULL,
           zero = NULL) {


  lpobs1 <- as.list(substitute(lpobs1))
  epobs1 <- link2list(lpobs1)
  lpobs1 <- attr(epobs1, "function.name")

  llambd <- as.list(substitute(llambda))
  elambd <- link2list(llambd)
  llambd <- attr(elambd, "function.name")


 type.fitted <- match.arg(type.fitted,
                           c("mean", "lambda", "pobs1", "onempobs1"))[1]


  new("vglmff",
  blurb = c("One-altered positive-Poisson distribution \n",
            "(Bernoulli and 1-truncated positive-Poisson ",
            "distribution model)\n\n",
            "Links:    ",
            namesof("pobs1",  lpobs1, earg = epobs1, tag = FALSE), ", ",
            namesof("lambda", llambd, earg = elambd, tag = FALSE)),
  constraints = eval(substitute(expression({
    dotzero <- .zero
    M1 <- 2
    eval(negzero.expression.VGAM)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = TRUE,
         parameters.names = c("pobs1", "lambda"),
         zero = .zero )
  }, list( .zero = zero ))),


  initialize = eval(substitute(expression({
    M1 <- 2
    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y

    extra$y1 <- y1 <- ifelse(y == 1, 1, 0)
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$skip.these <- skip.these <- matrix(as.logical(y1), n, NOS)
    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)

    mynames1 <- param.names("pobs1",  ncoly, skip1 = TRUE)
    mynames2 <- param.names("lambda", ncoly, skip1 = TRUE)
    predictors.names <-
        c(namesof(mynames1, .lpobs1 , earg = .epobs1 , tag = FALSE),
          namesof(mynames2, .llambd , earg = .elambd , tag = FALSE))[
          interleave.VGAM(M1*NOS, M1 = M1)]

    ncoly <- ncol(y)
    extra$ncoly <- ncoly
    M <- M1 * ncoly


    if (!length(etastart)) {
      Init.lambda <- y - 0.25
      etastart <-
        cbind(theta2eta(if (length( .ipobs1 )) .ipobs1 else
                        (0.5 + w * y1) / (1 + w),
                        .lpobs1 , earg = .epobs1 ),
              theta2eta(Init.lambda, .llambd , earg = .elambd ))
      etastart <- etastart[, interleave.VGAM(ncol(etastart), M1 = M1)]
    }
  }), list( .llambd = llambd, .elambd = elambd,
            .lpobs1 = lpobs1, .epobs1 = epobs1,
            .ipobs1 = ipobs1,  # .ilambd = ilambd,
            .type.fitted = type.fitted
           ))),


  linkinv = eval(substitute(function(eta, extra = NULL) {
   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "lambda", "pobs1", "onempobs1"))[1]

    M1 <- 2
    NOS <- ncol(eta) / M1

    pobs1 <- cbind(eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                             .lpobs1 , earg = .epobs1 ))
    lambd <- cbind(eta2theta(eta[, M1*(1:NOS)-0, drop = FALSE],
                             .llambd , earg = .elambd ))


    ans <- switch(type.fitted,
                  "mean"      = pobs1 + (1 - pobs1) *
                                lambd / ppois(1, lambd, lower = FALSE),
                  "lambda"      = lambd,
                  "pobs1"     =      pobs1,  # P(Y=1)
                  "onempobs1" =  1 - pobs1)  # P(Y>1)
    label.cols.y(ans, colnames.y = extra$colnames.y, NOS = NOS)
  }, list( .lpobs1 = lpobs1, .llambd = llambd,
           .epobs1 = epobs1, .elambd = elambd ))),
  last = eval(substitute(expression({
    temp.names <- c(rep_len( .lpobs1 , NOS),
                    rep_len( .llambd , NOS))
    temp.names <- temp.names[interleave.VGAM(M1*NOS, M1 = M1)]
    misc$link  <- temp.names
    names(misc$link) <-
      c(mynames1, mynames2)[interleave.VGAM(M1*NOS, M1 = M1)]

    misc$earg <- vector("list", M1 * NOS)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:NOS) {
      misc$earg[[M1*ii-1]] <- .epobs1
      misc$earg[[M1*ii  ]] <- .elambd
    }
  }), list( .lpobs1 = lpobs1, .llambd = llambd,
            .epobs1 = epobs1, .elambd = elambd ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    pobs1  <- cbind(eta2theta(eta[, c(TRUE, FALSE), drop = FALSE],
                              .lpobs1, earg = .epobs1))
    lambd <- cbind(eta2theta(eta[, c(FALSE, TRUE), drop = FALSE],
                              .llambd, earg = .elambd ))

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * doapospois(x = y, pobs1 = pobs1, lambda = lambd,
                                   log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lpobs1 = lpobs1, .llambd = llambd,
           .epobs1 = epobs1, .elambd = elambd ))),
  vfamily = c("oapospoisson"),



  validparams = eval(substitute(function(eta, y, extra = NULL) {
    TFvec <- c(TRUE, FALSE)
    pobs1 <- eta2theta(eta[,  TFvec, drop = FALSE],
                       .lpobs1 , earg = .epobs1 )
    lambd <- eta2theta(eta[, !TFvec, drop = FALSE],
                       .llambd , earg = .elambd )
    okay1 <- all(is.finite(lambd)) && all(0 < lambd) &&
             all(is.finite(pobs1)) && all(0 < pobs1 & pobs1 < 1)
    okay1
  }, list( .lpobs1 = lpobs1, .llambd = llambd,
           .epobs1 = epobs1, .elambd = elambd ))),


  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    pobs1 <- eta2theta(eta[, c(TRUE, FALSE)], .lpobs1 , earg = .epobs1 )
    lambd <- eta2theta(eta[, c(FALSE, TRUE)], .llambd , earg = .elambd )
    roapospois(nsim * length(lambd), lambd = lambd, pobs1 = pobs1)
  }, list( .lpobs1 = lpobs1, .llambd = llambd,
           .epobs1 = epobs1, .elambd = elambd ))),



  deriv = eval(substitute(expression({
    M1 <- 2
    NOS <- ncol(eta) / M1  # extra$NOS
    y1 <- extra$y1
    skip <- extra$skip.these

    TFvec <- c(TRUE, FALSE)
    pobs1  <- eta2theta(eta[,  TFvec, drop = FALSE],
                        .lpobs1 , earg = .epobs1 )
    lambda <- eta2theta(eta[, !TFvec, drop = FALSE],
                        .llambd , earg = .elambd )

    EY.cond <- 1 / ppois(1, lambda, lower.tail = FALSE)
    temp1 <- expm1(lambda)
    temp0 <- lambda * exp(-lambda)
    shape.geq.2 <- -expm1(-lambda) - temp0
    dl.dlambd <- y / lambda - 1 - temp0 / shape.geq.2


    dl.dpobs1 <- -1 / (1 - pobs1)  # For y > 1 obsns

    for (spp. in 1:NOS) {
      dl.dpobs1[skip[, spp.], spp.] <- 1 / pobs1[skip[, spp.], spp.]
      dl.dlambd[skip[, spp.], spp.] <- 0
    }
    dlambd.deta <- dtheta.deta(lambda, .llambd , earg = .elambd )
    mu.phi1 <- pobs1

    temp3 <- if ( .lpobs1 == "logitlink") {
      c(w) * (y1 - mu.phi1)
    } else {
      c(w) * dtheta.deta(mu.phi1, link = .lpobs1 , earg = .epobs1 ) *
             dl.dpobs1
    }

    ans <- cbind(temp3,
                 c(w) * dl.dlambd * dlambd.deta)
    ans <- ans[, interleave.VGAM(ncol(ans), M1 = M1)]
    ans
  }), list( .lpobs1 = lpobs1, .llambd = llambd,
            .epobs1 = epobs1, .elambd = elambd ))),
  weight = eval(substitute(expression({
    wz <- matrix(0, n, M1 * NOS)  # EIM is diagonal

    ned2l.dlambd2 <- EY.cond / lambda +
        ((1 - lambda) *
         exp(-lambda) - temp0^2 / shape.geq.2) / shape.geq.2

    ned2l.dlambd2 <- (1 - pobs1) * ned2l.dlambd2
    wz[, NOS+(1:NOS)] <- c(w) * ned2l.dlambd2 * dlambd.deta^2


    tmp100 <- mu.phi1 * (1 - mu.phi1)
    tmp200 <- if ( .lpobs1 == "logitlink" && is.empty.list( .epobs1 )) {
        cbind(c(w) * tmp100)
    } else {
      cbind(c(w) * (1 / tmp100) *
            dtheta.deta(mu.phi1, link = .lpobs1 , earg = .epobs1 )^2)
    }
    wz[, 1:NOS] <-  tmp200

    wz <- wz[, interleave.VGAM(ncol(wz), M1 = M1)]

    wz
  }), list( .lpobs1 = lpobs1,
            .epobs1 = epobs1 ))))
}  # End of oapospoisson






doazeta <- function(x, shape, pobs1 = 0, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(shape), length(pobs1))
  if (length(x)     != LLL) x     <- rep_len(x,     LLL)
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  if (length(pobs1) != LLL) pobs1 <- rep_len(pobs1, LLL)
  ans <- rep_len(0.0, LLL)


  index1 <- (x == 1)

  if (log.arg) {
    ans[ index1] <- log(pobs1[index1])
    ans[!index1] <- log1p(-pobs1[!index1]) +
                    dotzeta(x[!index1], shape[!index1], log = TRUE)
  } else {
    ans[ index1] <- pobs1[index1]
    ans[!index1] <- (1 - pobs1[!index1]) *
                    dotzeta(x[!index1], shape[!index1])
  }
  ans[pobs1 < 0 | 1 < pobs1] <- NaN
  ans[shape <= 0] <- NaN
  ans
}



poazeta <- function(q, shape, pobs1 = 0) {
  LLL <- max(length(q), length(shape), length(pobs1))
  if (length(q)     != LLL) q     <- rep_len(q,     LLL)
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  if (length(pobs1) != LLL) pobs1 <- rep_len(pobs1, LLL)
  ans <- rep_len(0.0, LLL)

  ans[q >  1] <-    pobs1[q > 1] +
                 (1-pobs1[q > 1]) * potzeta(q[q > 1], shape[q > 1])
  ans[q <  1] <- 0
  ans[q == 1] <- pobs1[q == 1]

  ans <- pmax(0, ans)
  ans <- pmin(1, ans)
  ans[pobs1 < 0 | 1 < pobs1] <- NaN
  ans[shape <= 0] <- NaN

  ans
}



qoazeta <- function(p, shape, pobs1 = 0) {
  LLL <- max(length(p), length(shape), length(pobs1))
  if (length(p)     != LLL) p     <- rep_len(p,     LLL)
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  if (length(pobs1) != LLL) pobs1 <- rep_len(pobs1, LLL)

  ans <- rep_len(NaN, LLL)
  ind4 <- pobs1 < p
  ans[!ind4] <- 1
  ans[ ind4] <- qotzeta((p[ind4] - pobs1[ind4]) / (1 - pobs1[ind4]),
                        shape = shape[ind4])
  ans[pobs1 < 0 | 1 < pobs1] <- NaN
  ans[p < 0 | 1 < p] <- NaN
  ans[shape <= 0] <- NaN
  ans
}



roazeta <- function(n, shape, pobs1 = 0) {
  qoazeta(runif(n), shape = shape, pobs1 = pobs1)
}






 oazeta <-
  function(lpobs1 = "logitlink",
           lshape = "loglink",
           type.fitted = c("mean", "shape", "pobs1", "onempobs1"),
           gshape = exp((-4:3)/4),
           ishape = NULL,
           ipobs1 = NULL,
           zero = NULL) {


  lpobs1 <- as.list(substitute(lpobs1))
  epobs1 <- link2list(lpobs1)
  lpobs1 <- attr(epobs1, "function.name")

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")


 type.fitted <- match.arg(type.fitted,
                          c("mean", "shape", "pobs1", "onempobs1"))[1]


  new("vglmff",
  blurb = c("One-altered zeta distribution \n",
            "(Bernoulli and 1-truncated zeta distribution model)\n\n",
            "Links:    ",
            namesof("pobs1", lpobs1, earg = epobs1, tag = FALSE), ", ",
            namesof("shape", lshape, earg = eshape, tag = FALSE)),
  constraints = eval(substitute(expression({
    dotzero <- .zero
    M1 <- 2
    eval(negzero.expression.VGAM)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = TRUE,
         parameters.names = c("pobs1", "shape"),
         zero = .zero )
  }, list( .zero = zero ))),


  initialize = eval(substitute(expression({
    M1 <- 2
    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y

    extra$y1 <- y1 <- ifelse(y == 1, 1, 0)
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$skip.these <- skip.these <- matrix(as.logical(y1), n, NOS)
    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)

    mynames1 <- param.names("pobs1", ncoly, skip1 = TRUE)
    mynames2 <- param.names("shape", ncoly, skip1 = TRUE)
    predictors.names <-
        c(namesof(mynames1, .lpobs1 , earg = .epobs1 , tag = FALSE),
          namesof(mynames2, .lshape , earg = .eshape , tag = FALSE))[
          interleave.VGAM(M1*NOS, M1 = M1)]

    ncoly <- ncol(y)
    extra$ncoly <- ncoly
    M <- M1 * ncoly


    if (!length(etastart)) {
      otzetaff.Loglikfun <- function(shape, y, x, w, extraargs) {
        sum(c(w) * dotzeta(x = y, shape, log = TRUE))
      }

      gshape <- .gshape
      if (!length( .ishape )) {
        shape.init <- matrix(NA_real_, n, M/M1, byrow = TRUE)
        for (jay in 1:ncoly) {
          index1 <- y[, jay] > 1
          shape.init[, jay] <-
            grid.search(gshape, objfun = otzetaff.Loglikfun,  # x = x,
                        y = y[index1, jay], w = w[index1, jay])
        }
      } else {
        shape.init <- matrix( .ishape , n, M, byrow = TRUE)
      }

      etastart <-
        cbind(theta2eta(if (length( .ipobs1 )) .ipobs1 else
                        (0.5 + w * y1) / (1 + w),
                        .lpobs1 , earg = .epobs1 ),
              theta2eta(shape.init, .lshape , earg = .eshape ))
      etastart <- etastart[, interleave.VGAM(ncol(etastart), M1 = M1)]
    }
  }), list( .lshape = lshape, .eshape = eshape,
            .lpobs1 = lpobs1, .epobs1 = epobs1,
            .ipobs1 = ipobs1, .ishape = ishape,
                              .gshape = gshape,
            .type.fitted = type.fitted
           ))),


  linkinv = eval(substitute(function(eta, extra = NULL) {
   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "shape", "pobs1", "onempobs1"))[1]

    M1 <- 2
    NOS <- ncol(eta) / M1

    pobs1 <- cbind(eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                             .lpobs1 , earg = .epobs1 ))
    shape <- cbind(eta2theta(eta[, M1*(1:NOS)-0, drop = FALSE],
                             .lshape , earg = .eshape ))
    if (type.fitted == "mean") {
      ans <- shape
      ans[shape > 1] <- zeta(shape[shape > 1])/zeta(shape[shape > 1] + 1)
      ans[shape <= 1] <- NA
      pmf.1 <- dzeta(1, shape)
      mean.otzeta <- (ans - pmf.1) / (1 - pmf.1)
    }

    ans <- switch(type.fitted,
                  "mean"      = pobs1 + (1 - pobs1) * mean.otzeta,
                  "shape"     = shape,
                  "pobs1"     =      pobs1,  # P(Y=1)
                  "onempobs1" =  1 - pobs1)  # P(Y>1)
    label.cols.y(ans, colnames.y = extra$colnames.y, NOS = NOS)
  }, list( .lpobs1 = lpobs1, .lshape = lshape,
           .epobs1 = epobs1, .eshape = eshape ))),
  last = eval(substitute(expression({
    temp.names <- c(rep_len( .lpobs1 , NOS),
                    rep_len( .lshape , NOS))
    temp.names <- temp.names[interleave.VGAM(M1*NOS, M1 = M1)]
    misc$link  <- temp.names
    names(misc$link) <-
      c(mynames1, mynames2)[interleave.VGAM(M1*NOS, M1 = M1)]

    misc$earg <- vector("list", M1 * NOS)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:NOS) {
      misc$earg[[M1*ii-1]] <- .epobs1
      misc$earg[[M1*ii  ]] <- .eshape
    }
  }), list( .lpobs1 = lpobs1, .lshape = lshape,
            .epobs1 = epobs1, .eshape = eshape ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    pobs1 <- cbind(eta2theta(eta[, c(TRUE, FALSE), drop = FALSE],
                             .lpobs1, earg = .epobs1 ))
    shape <- cbind(eta2theta(eta[, c(FALSE, TRUE), drop = FALSE],
                             .lshape, earg = .eshape ))

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * doazeta(x = y, pobs1 = pobs1, shape = shape,
                                log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lpobs1 = lpobs1, .lshape = lshape,
           .epobs1 = epobs1, .eshape = eshape ))),
  vfamily = c("oazeta"),



  validparams = eval(substitute(function(eta, y, extra = NULL) {
    TFvec <- c(TRUE, FALSE)
    pobs1 <- eta2theta(eta[,  TFvec, drop = FALSE],
                       .lpobs1 , earg = .epobs1 )
    shape <- eta2theta(eta[, !TFvec, drop = FALSE],
                       .lshape , earg = .eshape )
    okay1 <- all(is.finite(shape)) && all(0 < shape) &&
             all(is.finite(pobs1)) && all(0 < pobs1 & pobs1 < 1)
    okay1
  }, list( .lpobs1 = lpobs1, .lshape = lshape,
           .epobs1 = epobs1, .eshape = eshape ))),


  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    pobs1 <- eta2theta(eta[, c(TRUE, FALSE)], .lpobs1 , earg = .epobs1 )
    shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , earg = .eshape )
    roazeta(nsim * length(shape), shape = shape, pobs1 = pobs1)
  }, list( .lpobs1 = lpobs1, .lshape = lshape,
           .epobs1 = epobs1, .eshape = eshape ))),



  deriv = eval(substitute(expression({
    M1 <- 2
    NOS <- ncol(eta) / M1  # extra$NOS
    y1 <- extra$y1
    skip <- extra$skip.these

    TFvec <- c(TRUE, FALSE)
    pobs1 <- eta2theta(eta[,  TFvec, drop = FALSE],
                       .lpobs1 , earg = .epobs1 )
    shape <- eta2theta(eta[, !TFvec, drop = FALSE],
                       .lshape , earg = .eshape )

    BBBB  <- zeta(shape + 1) - 1
    fred1 <- zeta(shape + 1, deriv = 1)
    dl.dshape <- -log(y) - fred1 / BBBB

    dl.dpobs1 <- -1 / (1 - pobs1)  # For y > 1 obsns

    for (spp. in 1:NOS) {
      dl.dpobs1[skip[, spp.], spp.] <- 1 / pobs1[skip[, spp.], spp.]
      dl.dshape[skip[, spp.], spp.] <- 0
    }
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )
    mu.phi1 <- pobs1

    temp3 <- if ( .lpobs1 == "logitlink") {
      c(w) * (y1 - mu.phi1)
    } else {
      c(w) * dtheta.deta(mu.phi1, link = .lpobs1 , earg = .epobs1 ) *
             dl.dpobs1
    }

    ans <- cbind(temp3,
                 c(w) * dl.dshape * dshape.deta)
    ans <- ans[, interleave.VGAM(ncol(ans), M1 = M1)]
    ans
  }), list( .lpobs1 = lpobs1, .lshape = lshape,
            .epobs1 = epobs1, .eshape = eshape ))),
  weight = eval(substitute(expression({
    wz <- matrix(0, n, M1 * NOS)  # EIM is diagonal

    ned2l.dshape2 <- (zeta(shape + 1, deriv = 2) - fred1^2 / BBBB) / BBBB

    ned2l.dshape2 <- (1 - pobs1) * ned2l.dshape2
    wz[, NOS+(1:NOS)] <- c(w) * ned2l.dshape2 * dshape.deta^2


    tmp100 <- mu.phi1 * (1 - mu.phi1)
    tmp200 <- if ( .lpobs1 == "logitlink" && is.empty.list( .epobs1 )) {
        cbind(c(w) * tmp100)
    } else {
      cbind(c(w) * (1 / tmp100) *
            dtheta.deta(mu.phi1, link = .lpobs1 , earg = .epobs1 )^2)
    }
    wz[, 1:NOS] <-  tmp200

    wz <- wz[, interleave.VGAM(ncol(wz), M1 = M1)]

    wz
  }), list( .lpobs1 = lpobs1,
            .epobs1 = epobs1 ))))
}  # End of oazeta








deflat.limit.oizeta  <- function(shape) {
  if (any(shape <= 0))
    stop("argument 'shape' must be positive")
  ans <- -dzeta(1, shape) / pzeta(1, shape, lower.tail = FALSE)
  ans
}



doizeta <- function(x, shape, pstr1 = 0, log = FALSE) {



  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(shape), length(pstr1))
  if (length(x)     != LLL) x     <- rep_len(x,     LLL)
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  if (length(pstr1) != LLL) pstr1 <- rep_len(pstr1, LLL)

  ans <- rep(NA_real_, LLL)
  index1 <- (x == 1)
  if (log.arg) {
    ans[ index1] <- log(pstr1[ index1] + (1 - pstr1[ index1]) *
                        dzeta(x[ index1], shape[ index1]))
    ans[!index1] <- log1p(-pstr1[!index1]) +
                         dzeta(x[!index1], shape[!index1], log = TRUE)
  } else {
    ans[ index1] <-      pstr1[ index1] + (1 - pstr1[ index1]) *
                       dzeta(x[ index1], shape[ index1])
    ans[!index1] <- (1 - pstr1[!index1]) *
                       dzeta(x[!index1], shape[!index1])
  }


  ans[pstr1 < deflat.limit.oizeta(shape)] <- NaN
  ans[pstr1 > 1] <- NaN

  ans
}  # doizeta




poizeta <- function(q, shape, pstr1 = 0) {

  LLL <- max(length(q), length(shape), length(pstr1))
  if (length(q)     != LLL) q     <- rep_len(q,     LLL)
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  if (length(pstr1) != LLL) pstr1 <- rep_len(pstr1, LLL)
  ans <- rep_len(NA_real_, LLL)
  deflat.limit <- deflat.limit.oizeta(shape)

  ans <- pzeta(q, shape)  #, lower.tail = lower.tail, log.p = log.p
  ans <- ifelse(q < 1, 0, pstr1 + (1 - pstr1) * ans)

  ans[pstr1 < deflat.limit] <- NaN
  ans[1 < pstr1] <- NaN
  ans[shape <= 0] <- NaN

  ans
}  # poizeta





qoizeta <- function(p, shape, pstr1 = 0) {

  LLL <- max(length(p), length(shape), length(pstr1))
  if (length(p)     != LLL) p     <- rep_len(p,     LLL)
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  if (length(pstr1) != LLL) pstr1 <- rep_len(pstr1, LLL)
  ans <- rep_len(NA_real_, LLL)
  deflat.limit <- deflat.limit.oizeta(shape)

  ans[p <= pstr1] <- 1
  pindex <- (deflat.limit <= pstr1) & (pstr1 < p)
  ans[pindex] <-
    qzeta((p[pindex] - pstr1[pindex]) / (1 - pstr1[pindex]),
          shape = shape[pindex])

  ans[pstr1 < deflat.limit] <- NaN
  ans[1 < pstr1] <- NaN

  ans[p < 0] <- NaN
  ans[1 < p] <- NaN
  ans[shape <= 0] <- NaN

  ans
}  # qoizeta



roizeta <- function(n, shape, pstr1 = 0) {
  qoizeta(runif(n), shape, pstr1 = pstr1)
}





 oizeta <-
  function(lpstr1 = "logitlink", lshape = "loglink",
           type.fitted = c("mean", "shape", "pobs1", "pstr1", "onempstr1"),
           ishape = NULL,
           gpstr1 = ppoints(8),
           gshape = exp((-3:3) / 4), # grid for finding shape.init
           zero = NULL) {

  lpstr1 <- as.list(substitute(lpstr1))
  epstr1 <- link2list(lpstr1)
  lpstr1 <- attr(epstr1, "function.name")

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")



  type.fitted <- match.arg(type.fitted,
                   c("mean", "shape", "pobs1", "pstr1", "onempstr1"))[1]


  if (length(ishape))
    if (!is.Numeric(ishape, positive = TRUE))
      stop("argument 'ishape' values must be positive")


  new("vglmff",
  blurb = c("One-inflated zeta regression\n\n",
            "Links:    ",
            namesof("pstr1",  lpstr1, earg = epstr1 ), ", ",
            namesof("shape", lshape, earg = eshape ), "\n",
            "Mean:     pstr1 + (1 - pstr1) * zeta(shape) / ",
                       "zeta(1 + shape), if shape > 1"),

  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = TRUE,
         parameters.names = c("pstr1", "shape"),
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted
         ))),
  initialize = eval(substitute(expression({
    M1 <- 2

    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              Is.nonnegative.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    NOS <- ncoly <- ncol(y)
    extra$ncoly <- ncoly
    M <- M1 * ncoly
    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)

    mynames1 <- param.names("pstr1", ncoly, skip1 = TRUE)
    mynames2 <- param.names("shape", ncoly, skip1 = TRUE)
    predictors.names <-
        c(namesof(mynames1, .lpstr1 , earg = .epstr1 , tag = FALSE),
          namesof(mynames2, .lshape , earg = .eshape , tag = FALSE))[
          interleave.VGAM(M, M1 = M1)]


    if (!length(etastart)) {

      shape.init <-
      pstr1.init <- matrix(NA_real_, n, NOS)
      gpstr1 <- .gpstr1
      gshape <- .gshape

      oizeta.Loglikfun <- function(pstr1, shape, y, x, w, extraargs) {
        sum(c(w) * doizeta(x = y, pstr1 = pstr1,
                           shape = shape, log = TRUE))
      }


      for (jay in 1:NOS) {  # For each response 'y_jay'... do:


        try.this <-
          grid.search2(gpstr1, gshape,
                       objfun = oizeta.Loglikfun,
                       y = y[, jay],  # x = x[TFvec, , drop = FALSE],
                       w = w[, jay],
                       ret.objfun = TRUE)  # Last value is the loglik
        pstr1.init[, jay] <-  try.this["Value1"]
        shape.init[, jay] <-  try.this["Value2"]
      }  # for (jay ...)

      etastart <- cbind(theta2eta(pstr1.init, .lpstr1 , earg = .epstr1 ),
                        theta2eta(shape.init, .lshape , earg = .eshape ))[,
                        interleave.VGAM(M, M1 = M1)]
      mustart <- NULL  # Since etastart has been computed.
    }  # End of !length(etastart)
  }), list( .lpstr1 = lpstr1, .lshape = lshape,
            .epstr1 = epstr1, .eshape = eshape,
                                .ishape = ishape,
            .gpstr1 = gpstr1,
            .gshape = gshape,
            .type.fitted = type.fitted ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    NOS <- ncol(eta) / c(M1 = 2)
    type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "shape", "pobs1", "pstr1", "onempstr1"))[1]

    pstr1 <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr1 , earg = .epstr1 )
    shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , earg = .eshape )

    Meanfun <- function(shape) {
      Mean <- shape
      Mean[shape > 1] <-
        zeta(shape[shape > 1]) / zeta(1 + shape[shape > 1])
      Mean[shape <= 1] <- NA
      Mean
    }

    ans <-
      switch(type.fitted,
             "mean"      = pstr1 + (1 - pstr1) * Meanfun(shape),
             "shape"     = shape,
             "pobs1" = doizeta(1, shape = shape, pstr1 = pstr1),  # P(Y=1)
             "pstr1"     =     pstr1,
             "onempstr1" = 1 - pstr1)

    label.cols.y(ans, colnames.y = extra$colnames.y, NOS = NOS)
  }, list( .lpstr1 = lpstr1, .lshape = lshape,
           .epstr1 = epstr1, .eshape = eshape ))),
  last = eval(substitute(expression({
    misc$link <-
      c(rep_len( .lpstr1 , NOS),
        rep_len( .lshape , NOS))[interleave.VGAM(M, M1 = M1)]
    temp.names <- c(mynames1, mynames2)[interleave.VGAM(M, M1 = M1)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for (ii in 1:ncoly) {
      misc$earg[[M1*ii-1]] <- .epstr1
      misc$earg[[M1*ii  ]] <- .eshape
    }
  }), list( .lpstr1 = lpstr1, .lshape = lshape,
            .epstr1 = epstr1, .eshape = eshape ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    pstr1 <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr1 , earg = .epstr1 )
    shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , earg = .eshape )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * doizeta(x = y, pstr1 = pstr1, shape = shape,
                                log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lpstr1 = lpstr1, .lshape = lshape,
           .epstr1 = epstr1, .eshape = eshape ))),
  vfamily = c("oizeta"),



  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    pstr1 <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr1 , earg = .epstr1 )
    shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , earg = .eshape )
    roizeta(nsim * length(shape), shape = shape, pstr1 = pstr1)
  }, list( .lpstr1 = lpstr1, .lshape = lshape,
           .epstr1 = epstr1, .eshape = eshape ))),

  validparams = eval(substitute(function(eta, y, extra = NULL) {
    pstr1 <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE], .lpstr1 ,
                        earg = .epstr1 )
    shape <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE], .lshape ,
                        earg = .eshape )
    okay1 <- all(is.finite(shape)) && all(0 < shape) &&
             all(is.finite(pstr1)) && all(pstr1 < 1)
    deflat.limit <- deflat.limit.oizeta(shape)
    okay2.deflat <- TRUE
    if (okay1 && !(okay2.deflat <- all(deflat.limit < pstr1)))
      warning("parameter 'pstr1' is too negative even allowing for ",
              "1-deflation.")
    okay1 && okay2.deflat
  }, list( .lpstr1 = lpstr1, .lshape = lshape,
           .epstr1 = epstr1, .eshape = eshape ))),






  deriv = eval(substitute(expression({
    M1 <- 2
    NOS <- M / M1
    pstr1 <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE], .lpstr1 ,
                       earg = .epstr1 )
    shape <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE], .lshape ,
                       earg = .eshape )

    pmf1 <- dzeta(1, shape)
    onempmf1 <- 1 - pmf1  # dozeta(1, shape = shape, pstr1 = pstr1)
    pobs1 <- pstr1 + (1 - pstr1) * pmf1
    index1 <- as.matrix(y == 1)

    zeta0 <- zeta(shape + 1)
    zeta1 <- zeta(shape + 1, deriv = 1)
    zeta2 <- zeta(shape + 1, deriv = 2)

    dl.dpstr1 <- onempmf1 / pobs1
    dl.dpstr1[!index1] <- -1 / (1 - pstr1[!index1])

    dpmf1.dshape <- -zeta1 / zeta0^2

   d2pmf1.dshape2 <- (2 * zeta1^2 / zeta0 - zeta2) / zeta0^2

    dl.dshape <- (1 - pstr1) * dpmf1.dshape / pobs1  #
    dl.dshape[!index1] <- -log(y[!index1]) -
                          zeta1[!index1] / zeta0[!index1]
    
    dpstr1.deta <- dtheta.deta(pstr1, .lpstr1 , earg = .epstr1 )
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )

    myderiv <- c(w) * cbind(dl.dpstr1 * dpstr1.deta,
                            dl.dshape * dshape.deta)
    myderiv[, interleave.VGAM(M, M1 = M1)]
  }), list( .lpstr1 = lpstr1, .lshape = lshape,
            .epstr1 = epstr1, .eshape = eshape ))),
  weight = eval(substitute(expression({

    LHS <- ((1 - pstr1) / pobs1) * dpmf1.dshape^2 - d2pmf1.dshape2
    RHS <- (zeta2 - zeta1^2 / zeta0) / zeta0
    ned2l.dpstr12 <- onempmf1 / ((1 - pstr1) * pobs1)  #
    ned2l.dpstr1shape <- dpmf1.dshape / pobs1  #
    ned2l.dshape2 <- (1 - pstr1) * (LHS + (1 - pmf1) * RHS)

    wz <- array(c(c(w) * ned2l.dpstr12 * dpstr1.deta^2,
                  c(w) * ned2l.dshape2 * dshape.deta^2,
                  c(w) * ned2l.dpstr1shape * dpstr1.deta * dshape.deta),
                dim = c(n, M / M1, 3))
    wz <- arwz2wz(wz, M = M, M1 = M1)
    wz
  }), list( .lshape = lshape, .eshape = eshape ))))
}  # oizeta








deflat.limit.oizipf  <- function(N, shape) {
  if (any(shape <= 0))
    stop("argument 'shape' must be positive")
  ans <- 1 / (1 - 1 / dzipf(1, N, shape))
  ans
}





doizipf <- function(x, N, shape, pstr1 = 0, log = FALSE) {

  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  if (!is.Numeric(x))
    stop("bad input for argument 'x'")
  if (!is.Numeric(N, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'N'")
  if (!is.Numeric(shape, positive = TRUE))
    stop("bad input for argument 'shape'")
  nn <- max(length(x), length(N), length(shape), length(pstr1))
  if (length(x)    != nn) x     <- rep_len(x,     nn)
  if (length(N)    != nn) N     <- rep_len(N,     nn)
  if (length(shape)!= nn) shape <- rep_len(shape, nn)
  if (length(pstr1)!= nn) pstr1 <- rep_len(pstr1, nn)

  ans <- rep(NA_real_, nn)
  index1 <- (x == 1)
  if (log.arg) {
    ans[ index1] <- log(pstr1[ index1] + (1 - pstr1[ index1]) *
                      dzipf(x[ index1], N[ index1], shape[ index1]))
    ans[!index1] <- log1p(-pstr1[!index1]) +
                    dzipf(x[!index1], N[!index1], shape[!index1],
                          log = TRUE)
  } else {
    ans[ index1] <-      pstr1[ index1] + (1 - pstr1[ index1]) *
                       dzipf(x[ index1], N[ index1], shape[ index1])
    ans[!index1] <- (1 - pstr1[!index1]) *
                       dzipf(x[!index1], N[!index1], shape[!index1])
  }


  deflat.limit <- deflat.limit.oizipf(N, shape)
  ans[pstr1 < deflat.limit] <- NaN
  ans[pstr1 > 1] <- NaN

  ans
}





poizipf <- function(q, N, shape, pstr1 = 0) {

  LLL <- max(length(q), length(N), length(shape), length(pstr1))
  if (length(q)     != LLL) q     <- rep_len(q,     LLL)
  if (length(N)     != LLL) N     <- rep_len(N,     LLL)
  if (length(shape) != LLL) shape <- rep_len(shape, LLL)
  if (length(pstr1) != LLL) pstr1 <- rep_len(pstr1, LLL)
  ans <- rep_len(NA_real_, LLL)
  deflat.limit <- deflat.limit.oizipf(N, shape)

  ans <- pzipf(q, N, shape)  #, lower.tail = lower.tail, log.p = log.p
  ans <- ifelse(q < 1, 0, pstr1 + (1 - pstr1) * ans)

  ans[pstr1 < deflat.limit] <- NaN
  ans[1 < pstr1] <- NaN
  ans[s <= 0] <- NaN

  ans
}






qoizipf <- function(p, N, shape, pstr1 = 0) {

  if (!is.Numeric(p))
    stop("bad input for argument 'p'")
  if (!is.Numeric(N, integer.valued = TRUE, positive = TRUE))
    stop("bad input for argument 'N'")
  if (!is.Numeric(shape, positive = TRUE))
    stop("bad input for argument 'shape'")

  nn <- max(length(p), length(N), length(s), length(pstr1))
  if (length(p)     != nn) p     <- rep_len(p,     nn)
  if (length(N)     != nn) N     <- rep_len(N,     nn)
  if (length(shape) != nn) shape <- rep_len(shape, nn)
  if (length(pstr1) != nn) pstr1 <- rep_len(pstr1, nn)


  ans    <- rep_len(NA_real_, nn)
  deflat.limit <- deflat.limit.oizipf(N, shape)

  dont.iterate <- 1 < p
  ans[p <= pstr1] <- 1
  pindex <- (pstr1 < p) & (deflat.limit <= pstr1) & !dont.iterate
  if (any(pindex))
  ans[pindex] <-
    qzipf((p[pindex] - pstr1[pindex]) / (1 - pstr1[pindex]),
          N = N[pindex], shape = shape[pindex])

  ans[pstr1 < deflat.limit] <- NaN
  ans[1 < pstr1] <- NaN

  ans[shape < 0] <- NaN
  ans[p < 0] <- NaN
  ans[1 < p] <- NaN

  ans
}



roizipf <- function(n, N, shape, pstr1 = 0) {
  qoizipf(runif(n), N, shape, pstr1 = pstr1)
}





 oizipf <-
  function(N = NULL, lpstr1 = "logitlink", lshape = "loglink",
           type.fitted = c("mean", "shape", "pobs1", "pstr1", "onempstr1"),
           ishape = NULL,
           gpstr1 = ppoints(8),
           gshape = exp((-3:3) / 4), # grid for finding shape.init
           zero = NULL) {

  if (length(N) &&
     (!is.Numeric(N, positive = TRUE,
                 integer.valued = TRUE, length.arg = 1) ||
      N <= 1))
    stop("bad input for argument 'N'")
  enteredN <- length(N)

  lpstr1 <- as.list(substitute(lpstr1))
  epstr1 <- link2list(lpstr1)
  lpstr1 <- attr(epstr1, "function.name")

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")



  type.fitted <- match.arg(type.fitted,
                   c("mean", "shape", "pobs1", "pstr1", "onempstr1"))[1]


  if (length(ishape))
    if (!is.Numeric(ishape, positive = TRUE))
      stop("argument 'ishape' values must be positive")


  new("vglmff",
  blurb = c("One-inflated Zipf distribution ",
            "f(y; pstr1, shape) = pstr1 + ",
            "(1 - pstr1) * y^(-shape) / sum((1:N)^(-shape)),",
            " 0 < shape, y = 1, 2,...,N",
            ifelse(enteredN, paste(" = ", N, sep = ""), ""),
            "\n\n",
            "Links:    ",
            namesof("pstr1", lpstr1, earg = epstr1 ), ", ",
            namesof("shape", lshape, earg = eshape ), "\n",
            "Mean:     pstr1 + (1 - pstr1) * ",
            "gharmonic(N, shape-1) / gharmonic(N, shape)"),

  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = TRUE,
         parameters.names = c("pstr1", "shape"),
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted
         ))),
  initialize = eval(substitute(expression({
    M1 <- 2

    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    NOS <- ncoly <- ncol(y)
    extra$ncoly <- ncoly
    M <- M1 * ncoly
    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)



    NN <- .N
    if (!is.Numeric(NN, length.arg = 1,
                    positive = TRUE, integer.valued = TRUE))
      NN <- max(y)
    if (max(y) > NN)
      stop("maximum of the response is greater than argument 'N'")
    extra$N <- NN



    mynames1 <- param.names("pstr1", ncoly, skip1 = TRUE)
    mynames2 <- param.names("shape", ncoly, skip1 = TRUE)
    predictors.names <-
        c(namesof(mynames1, .lpstr1 , earg = .epstr1 , tag = FALSE),
          namesof(mynames2, .lshape , earg = .eshape , tag = FALSE))[
          interleave.VGAM(M, M1 = M1)]


    if (!length(etastart)) {

      shape.init <-
      pstr1.init <- matrix(NA_real_, n, NOS)
      gpstr1 <- .gpstr1
      gshape <- .gshape

      oizipf.Loglikfun <- function(pstr1, shape, y, x, w, extraargs) {
        sum(c(w) * doizipf(x = y, pstr1 = pstr1, N = extraargs$N,
                           s = shape, log = TRUE))
      }


      for (jay in 1:NOS) {  # For each response 'y_jay'... do:


        try.this <-
          grid.search2(gpstr1, gshape,
                       objfun = oizipf.Loglikfun,
                       y = y[, jay],  # x = x[TFvec, , drop = FALSE],
                       w = w[, jay],
                       extraargs = list(N = extra$N),
                       ret.objfun = TRUE)  # Last value is the loglik
        pstr1.init[, jay] <-  try.this["Value1"]
        shape.init[, jay] <-  try.this["Value2"]
      }  # for (jay ...)

      etastart <- cbind(theta2eta(pstr1.init, .lpstr1 , earg = .epstr1 ),
                        theta2eta(shape.init, .lshape , earg = .eshape ))[,
                        interleave.VGAM(M, M1 = M1)]
      mustart <- NULL  # Since etastart has been computed.
    }  # End of !length(etastart)
  }), list( .lpstr1 = lpstr1, .lshape = lshape,
            .epstr1 = epstr1, .eshape = eshape,
                              .ishape = ishape,
            .gpstr1 = gpstr1,
            .gshape = gshape,
            .type.fitted = type.fitted,
            .N = N ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    NOS <- ncol(eta) / c(M1 = 2)
    type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "shape", "pobs1", "pstr1", "onempstr1"))[1]

    pstr1 <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr1 , earg = .epstr1 )
    shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , earg = .eshape )

    Meanfun <- function(shape, extra) {
      Mean <- shape
      Mean <- ( gharmonic2(extra$N, shape = shape - 1)
              / gharmonic2(extra$N, shape = shape))
      Mean[shape <= 0] <- NaN
      Mean
    }

    ans <-
      switch(type.fitted,
             "mean"      = pstr1 + (1 - pstr1) * Meanfun(shape, extra),
             "shape"     = shape,
         "pobs1"     = doizipf(1, N = extra$N, s = shape, pstr1 = pstr1),
             "pstr1"     =     pstr1,
             "onempstr1" = 1 - pstr1)

    label.cols.y(ans, colnames.y = extra$colnames.y, NOS = NOS)
  }, list( .lpstr1 = lpstr1, .lshape = lshape,
           .epstr1 = epstr1, .eshape = eshape ))),
  last = eval(substitute(expression({
    misc$link <-
      c(rep_len( .lpstr1 , NOS),
        rep_len( .lshape , NOS))[interleave.VGAM(M, M1 = M1)]
    temp.names <- c(mynames1, mynames2)[interleave.VGAM(M, M1 = M1)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for (ii in 1:ncoly) {
      misc$earg[[M1*ii-1]] <- .epstr1
      misc$earg[[M1*ii  ]] <- .eshape
    }
  }), list( .lpstr1 = lpstr1, .lshape = lshape,
            .epstr1 = epstr1, .eshape = eshape ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    pstr1 <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr1 , earg = .epstr1 )
    shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , earg = .eshape )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * doizipf(x = y, pstr1 = pstr1, s = shape,
                                N = extra$N, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lpstr1 = lpstr1, .lshape = lshape,
           .epstr1 = epstr1, .eshape = eshape ))),
  vfamily = c("oizipf"),



  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    pstr1 <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr1 , earg = .epstr1 )
    shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , earg = .eshape )
    roizipf(nsim * length(shape), s = shape, pstr1 = pstr1,
            N = object@extra$N)
  }, list( .lpstr1 = lpstr1, .lshape = lshape,
           .epstr1 = epstr1, .eshape = eshape ))),

  validparams = eval(substitute(function(eta, y, extra = NULL) {
    pstr1 <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE], .lpstr1 ,
                       earg = .epstr1 )
    shape <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE], .lshape ,
                       earg = .eshape )
    okay1 <- all(is.finite(shape)) && all(0 < shape) &&
             all(is.finite(pstr1)) && all(pstr1 < 1)
    deflat.limit <- deflat.limit.oizipf(N = extra$N, s = shape)
    okay2.deflat <- TRUE
    if (okay1 && !(okay2.deflat <- all(deflat.limit < pstr1)))
      warning("parameter 'pstr1' is too negative even allowing for ",
              "1-deflation.")
    okay1 && okay2.deflat
  }, list( .lpstr1 = lpstr1, .lshape = lshape,
           .epstr1 = epstr1, .eshape = eshape ))),






  deriv = eval(substitute(expression({
    M1 <- 2
    NOS <- M / M1
    pstr1 <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE], .lpstr1 ,
                       earg = .epstr1 )
    shape <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE], .lshape ,
                       earg = .eshape )

    pmf1 <- dzipf(1, N = extra$N, shape = shape)
    onempmf1 <- 1 - pmf1  # dozeta(1, shape = shape, pstr1 = pstr1)
    pobs1 <- pstr1 + (1 - pstr1) * pmf1
    index1 <- as.matrix(y == 1)

    ghar0 <-  gharmonic2(extra$N, shape)
    ghar1 <-   gharmonic(extra$N, shape, deriv = 1)
    ghar2 <-   gharmonic(extra$N, shape, deriv = 2)

    dl.dpstr1 <- onempmf1 / pobs1
    dl.dpstr1[!index1] <- -1 / (1 - pstr1[!index1])

    dpmf1.dshape <- -ghar1 / ghar0^2

    d2pmf1.dshape2 <- (2 * ghar1^2 / ghar0 - ghar2) / ghar0^2

    dl.dshape <- (1 - pstr1) * dpmf1.dshape / pobs1  #
    dl.dshape[!index1] <- -log(y[!index1]) -
                          ghar1[!index1] / ghar0[!index1]

    dpstr1.deta <- dtheta.deta(pstr1, .lpstr1 , earg = .epstr1 )
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )

    myderiv <- c(w) * cbind(dl.dpstr1 * dpstr1.deta,
                            dl.dshape * dshape.deta)
    myderiv[, interleave.VGAM(M, M1 = M1)]
  }), list( .lpstr1 = lpstr1, .lshape = lshape,
            .epstr1 = epstr1, .eshape = eshape ))),
  weight = eval(substitute(expression({

    LHS <- ((1 - pstr1) / pobs1) * dpmf1.dshape^2 - d2pmf1.dshape2
    RHS <- (ghar2 - ghar1^2 / ghar0) / ghar0
    ned2l.dpstr12 <- onempmf1 / ((1 - pstr1) * pobs1)  #
    ned2l.dpstr1shape <- dpmf1.dshape / pobs1  #
    ned2l.dshape2 <- (1 - pstr1) * (LHS + (1 - pmf1) * RHS)

    wz <- array(c(c(w) * ned2l.dpstr12 * dpstr1.deta^2,
                  c(w) * ned2l.dshape2 * dshape.deta^2,
                  c(w) * ned2l.dpstr1shape * dpstr1.deta * dshape.deta),
                dim = c(n, M / M1, 3))
    wz <- arwz2wz(wz, M = M, M1 = M1)
    wz
  }), list( .lshape = lshape, .eshape = eshape ))))
}  # oizipf








dotzeta <- function(x, shape, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  if (log.arg) {
    ans <- dzeta(x, shape, log = log.arg) - log1p(-dzeta(1, shape))
    ans[x == 1] <- log(0)
  } else {
    ans <- dzeta(x, shape) / (1 - dzeta(1, shape))
    ans[x == 1] <- 0
  }
  ans[shape < 0] <- NaN
  ans
}  # dotzeta



potzeta <- function(q, shape, log.p = FALSE) {
  if (log.p) log(pzeta(q, shape) - dzeta(1, shape)) -
      log1p(-dzeta(1, shape)) else
    (pzeta(q, shape) - dzeta(1, shape)) / (1 - dzeta(1, shape))
}



 qotzeta <- function(p, shape) {
  ans <- qzeta((1 - dzeta(1, shape)) * p +
               dzeta(1, shape), shape = shape)
  ans[p == 1] <- Inf
  ans[p < 0 | 1 < p] <- NaN

  ans[shape < 0] <- NaN
  ans
}  # qotzeta



rotzeta <- function(n, shape) {
  qotzeta(runif(n), shape)
}





 otzeta <-
    function(lshape = "loglink",
             ishape = NULL,
             gshape = exp((-4:3)/4),
             zero = NULL) {

  if (length(ishape) && !is.Numeric(ishape, positive = TRUE))
    stop("argument 'ishape' must be > 0")

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")


  new("vglmff",
  blurb = c("One-truncated Zeta distribution ",
            "f(y; shape) = 1/(y^(shape+1) * (zeta(shape+1) - ",
            "1 - 1/2^(shape+1)))",
            " 0<shape, y = 2, 3,...\n\n",
            "Link:    ",
            namesof("shape", lshape, earg = eshape)),
  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         expected = FALSE,  # NR == FS
         multipleResponses = TRUE,
         parameters.names = "shape",
         zero = .zero ,
         lshape = .lshape )
  }, list( .lshape = lshape,
           .zero = zero ))),
  initialize = eval(substitute(expression({

   temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              Is.integer.y = TRUE,
              Is.positive.y = TRUE,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
    if (any(y <= 1))
      stop("no 1s in the response allowed!")


    ncoly <- ncol(y)
    mynames1 <- param.names("shape", ncoly, skip1 = TRUE)
    predictors.names <-
      namesof(mynames1, .lshape , earg = .eshape , tag = FALSE)

    M1 <- 1
    M <- M1 * ncoly


    if (!length(etastart)) {
      otzetaff.Loglikfun <- function(shape, y, x, w, extraargs) {
        sum(c(w) * dotzeta(x = y, shape, log = TRUE))
      }


      gshape <- .gshape
      if (!length( .ishape )) {
        shape.init <- matrix(NA_real_, n, M, byrow = TRUE)
        for (jay in 1:ncoly) {
          shape.init[, jay] <-
            grid.search(gshape, objfun = otzetaff.Loglikfun,
                        y = y[, jay], x = x, w = w[, jay])
        }
      } else {
        shape.init <- matrix( .ishape , n, M, byrow = TRUE)
      }
      etastart <- theta2eta(shape.init, .lshape , earg = .eshape )
    }
  }), list( .lshape = lshape, .eshape = eshape,
            .ishape = ishape, .gshape = gshape ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    ans <- pp <- eta2theta(eta, .lshape , earg = .eshape )
    ans[pp > 1] <- zeta(pp[pp > 1]) / zeta(pp[pp > 1] + 1)
    ans[pp <= 1] <- NA
    pmf.1 <- dzeta(1, pp)
    (ans - pmf.1) / (1 - pmf.1)
  }, list( .lshape = lshape, .eshape = eshape ))),
  last = eval(substitute(expression({
    misc$link <- rep_len( .lshape , ncoly)
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for (jay in 1:ncoly) {
      misc$earg[[jay]] <- .eshape
    }

  }), list( .lshape = lshape, .eshape = eshape ))),
  loglikelihood = eval(substitute( function(mu, y, w,
            residuals = FALSE, eta, extra = NULL, summation = TRUE) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dotzeta(x = y, shape, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lshape = lshape, .eshape = eshape ))),
  vfamily = c("otzeta"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    okay1 <- all(is.finite(shape)) && all(0 < shape)
    okay1
  }, list( .lshape = lshape, .eshape = eshape ))),
  deriv = eval(substitute(expression({
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    BBBB  <- zeta(shape + 1) - 1
    fred1 <- zeta(shape + 1, deriv = 1)
    dl.dshape <- -log(y) - fred1 / BBBB
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )
    c(w) * dl.dshape * dshape.deta
  }), list( .lshape = lshape, .eshape = eshape ))),
  weight = expression({
    ned2l.dshape2 <- (zeta(shape + 1, deriv = 2) -
                      fred1^2 / BBBB) / BBBB
    wz <- ned2l.dshape2 * dshape.deta^2
    c(w) * wz
  }))
}





dotzeta <- function(x, shape, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  if (log.arg) {
    ans <- dzeta(x, shape, log = log.arg) - log1p(-dzeta(1, shape))
    ans[x == 1] <- log(0)
  } else {
    ans <- dzeta(x, shape) / (1 - dzeta(1, shape))
    ans[x == 1] <- 0
  }
  ans[shape < 0] <- NaN
  ans
}  # dotzeta



potzeta <- function(q, shape, log.p = FALSE) {
  if (log.p) log(pzeta(q, shape) - dzeta(1, shape)) -
      log1p(-dzeta(1, shape)) else
    (pzeta(q, shape) - dzeta(1, shape)) / (1 - dzeta(1, shape))
}



 qotzeta <- function(p, shape) {
  ans <- qzeta((1 - dzeta(1, shape)) * p +
                dzeta(1, shape), shape = shape)

  ans[p == 1] <- Inf
  ans[p < 0 | 1 < p] <- NaN

  ans[shape < 0] <- NaN
  ans
}  # qotzeta



rotzeta <- function(n, shape) {
  qotzeta(runif(n), shape)
}





 otzeta <-
    function(lshape = "loglink",
             ishape = NULL,
             gshape = exp((-4:3)/4),
             zero = NULL) {

  if (length(ishape) && !is.Numeric(ishape, positive = TRUE))
    stop("argument 'ishape' must be > 0")

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")


  new("vglmff",
  blurb = c("One-truncated Zeta distribution ",
            "f(y; shape) = 1/(y^(shape+1) * (zeta(shape+1) - ",
            "1 - 1/2^(shape+1)))",
            " 0<shape, y = 2, 3,...\n\n",
            "Link:    ",
            namesof("shape", lshape, earg = eshape)),
  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         expected = FALSE,  # NR == FS
         multipleResponses = TRUE,
         parameters.names = "shape",
         zero = .zero ,
         lshape = .lshape )
  }, list( .lshape = lshape,
           .zero = zero ))),
  initialize = eval(substitute(expression({

   temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              Is.integer.y = TRUE,
              Is.positive.y = TRUE,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
    if (any(y <= 1))
      stop("no 1s in the response allowed!")


    ncoly <- ncol(y)
    mynames1 <- param.names("shape", ncoly, skip1 = TRUE)
    predictors.names <-
      namesof(mynames1, .lshape , earg = .eshape , tag = FALSE)

    M1 <- 1
    M <- M1 * ncoly


    if (!length(etastart)) {
      otzetaff.Loglikfun <- function(shape, y, x, w, extraargs) {
        sum(c(w) * dotzeta(x = y, shape, log = TRUE))
      }


      gshape <- .gshape
      if (!length( .ishape )) {
        shape.init <- matrix(NA_real_, n, M, byrow = TRUE)
        for (jay in 1:ncoly) {
          shape.init[, jay] <-
            grid.search(gshape, objfun = otzetaff.Loglikfun,
                        y = y[, jay], x = x, w = w[, jay])
        }
      } else {
        shape.init <- matrix( .ishape , n, M, byrow = TRUE)
      }
      etastart <- theta2eta(shape.init, .lshape , earg = .eshape )
    }
  }), list( .lshape = lshape, .eshape = eshape,
            .ishape = ishape, .gshape = gshape ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    ans <- pp <- eta2theta(eta, .lshape , earg = .eshape )
    ans[pp > 1] <- zeta(pp[pp > 1]) / zeta(pp[pp > 1] + 1)
    ans[pp <= 1] <- NA
    pmf.1 <- dzeta(1, pp)
    (ans - pmf.1) / (1 - pmf.1)
  }, list( .lshape = lshape, .eshape = eshape ))),
  last = eval(substitute(expression({
    misc$link <- rep_len( .lshape , ncoly)
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for (jay in 1:ncoly) {
      misc$earg[[jay]] <- .eshape
    }

  }), list( .lshape = lshape, .eshape = eshape ))),
  loglikelihood = eval(substitute( function(mu, y, w,
           residuals = FALSE, eta, extra = NULL, summation = TRUE) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dotzeta(x = y, shape, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lshape = lshape, .eshape = eshape ))),
  vfamily = c("otzeta"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    okay1 <- all(is.finite(shape)) && all(0 < shape)
    okay1
  }, list( .lshape = lshape, .eshape = eshape ))),
  deriv = eval(substitute(expression({
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    BBBB  <- zeta(shape + 1) - 1
    fred1 <- zeta(shape + 1, deriv = 1)
    dl.dshape <- -log(y) - fred1 / BBBB
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )
    c(w) * dl.dshape * dshape.deta
  }), list( .lshape = lshape, .eshape = eshape ))),
  weight = expression({
    ned2l.dshape2 <- (zeta(shape + 1, deriv = 2) -
                      fred1^2 / BBBB) / BBBB
    wz <- ned2l.dshape2 * dshape.deta^2
    c(w) * wz
  }))
}  # otzeta








deflat.limit.oipospois  <- function(lambda) {
  if (any(lambda < 0))
    stop("argument 'lambda' cannot be negative")
  ans <- -lambda / (expm1(lambda) - lambda)
  ans[is.infinite(lambda)] <- 0
  ans
}


doipospois <- function(x, lambda, pstr1 = 0, log = FALSE) {



  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(lambda), length(pstr1))
  if (length(x)      != LLL) x      <- rep_len(x,      LLL)
  if (length(lambda) != LLL) lambda <- rep_len(lambda, LLL)
  if (length(pstr1)  != LLL) pstr1  <- rep_len(pstr1,  LLL)

  ans <- rep(NA_real_, LLL)
  index1 <- (x == 1)
  if (log.arg) {
    ans[ index1] <- log(pstr1[ index1] + (1 - pstr1[ index1]) *
                        dpospois(x[ index1], lambda[ index1]))
    ans[!index1] <- log1p(-pstr1[!index1]) +
                    dpospois(x[!index1], lambda[!index1], log = TRUE)
  } else {
    ans[ index1] <-      pstr1[ index1] + (1 - pstr1[ index1]) *
                       dpospois(x[ index1], lambda[ index1])
    ans[!index1] <- (1 - pstr1[!index1]) *
                    dpospois(x[!index1], lambda[!index1])
  }


  deflat.limit <- deflat.limit.oipospois(lambda)
  ans[pstr1 < deflat.limit] <- NaN
  ans[pstr1 > 1] <- NaN

  ans
}  # doipospois






poipospois <- function(q, lambda, pstr1 = 0) {

  LLL <- max(length(q), length(lambda), length(pstr1))
  if (length(q)      != LLL) q      <- rep_len(q,      LLL)
  if (length(lambda) != LLL) lambda <- rep_len(lambda, LLL)
  if (length(pstr1)  != LLL) pstr1  <- rep_len(pstr1,  LLL)
  ans <- rep_len(NA_real_, LLL)
  deflat.limit <- deflat.limit.oipospois(lambda)

  ans <- ppospois(q, lambda)  #, lower.tail = lower.tail, log.p = log.p
  ans <- ifelse(q < 1, 0, pstr1 + (1 - pstr1) * ans)

  ans[pstr1 < deflat.limit] <- NaN
  ans[1 < pstr1] <- NaN
  ans[lambda <= 0] <- NaN

  ans
}  # poipospois




qoipospois <- function(p, lambda, pstr1 = 0) {

  LLL <- max(length(p), length(lambda), length(pstr1))
  if (length(p)      != LLL) p      <- rep_len(p,      LLL)
  if (length(lambda) != LLL) lambda <- rep_len(lambda, LLL)
  if (length(pstr1)  != LLL) pstr1  <- rep_len(pstr1,  LLL)
  ans <- rep_len(NA_real_, LLL)
  deflat.limit <- deflat.limit.oipospois(lambda)

  ans[p <= pstr1] <- 1
  pindex <- (deflat.limit <= pstr1) & (pstr1 < p)
  ans[pindex] <-
    qpospois((p[pindex] - pstr1[pindex]) / (1 - pstr1[pindex]),
             lambda = lambda[pindex])

  ans[pstr1 < deflat.limit] <- NaN
  ans[1 < pstr1] <- NaN

  ans[p < 0] <- NaN
  ans[1 < p] <- NaN
  ans[lambda <= 0] <- NaN

  ans
}  # qoipospois



roipospois <- function(n, lambda, pstr1 = 0) {

  ans <- qoipospois(runif(n), lambda, pstr1 = pstr1)
  ans
}  # roipospois






 oipospoisson <-
  function(lpstr1 = "logitlink", llambda = "loglink",
       type.fitted = c("mean", "lambda", "pobs1", "pstr1", "onempstr1"),
           ilambda = NULL,
           gpstr1 = (1:19)/20,
           gprobs.y = (1:19)/20,  # 20160518; grid for finding lambd.init
           imethod = 1,
           zero = NULL) {

  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 2)
    stop("argument 'imethod' must be 1 or 2")


  gpstr10 <- gpstr1


  lpstr10 <- as.list(substitute(lpstr1))
  epstr10 <- link2list(lpstr10)
  lpstr10 <- attr(epstr10, "function.name")

  llambda <- as.list(substitute(llambda))
  elambda <- link2list(llambda)
  llambda <- attr(elambda, "function.name")



  type.fitted <- match.arg(type.fitted,
                   c("mean", "lambda", "pobs1", "pstr1", "onempstr1"))[1]


  if (length(ilambda))
    if (!is.Numeric(ilambda, positive = TRUE))
      stop("argument 'ilambda' values must be positive")


  new("vglmff",
  blurb = c("One-inflated positive Poisson\n\n",
            "Links:    ",
            namesof("pstr1",  lpstr10, earg = epstr10 ), ", ",
            namesof("lambda", llambda, earg = elambda ), "\n",
            "Mean:     pstr1 + (1 - pstr1) * lambda / (1 - exp(-lambda))"),

  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = TRUE,
         imethod = .imethod ,
         multipleResponses = TRUE,
         parameters.names = c("pstr1", "lambda"),
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .imethod = imethod,
           .type.fitted = type.fitted
         ))),
  initialize = eval(substitute(expression({
    M1 <- 2

    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              Is.nonnegative.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    NOS <- ncoly <- ncol(y)
    extra$ncoly <- ncoly
    M <- M1 * ncoly
    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)

    mynames1 <- param.names("pstr1",  ncoly, skip1 = TRUE)
    mynames2 <- param.names("lambda", ncoly, skip1 = TRUE)
    predictors.names <-
        c(namesof(mynames1, .lpstr10 , earg = .epstr10 , tag = FALSE),
          namesof(mynames2, .llambda , earg = .elambda , tag = FALSE))[
          interleave.VGAM(M, M1 = M1)]


    if (!length(etastart)) {

      lambd.init <-
      pstr1.init <- matrix(NA_real_, n, NOS)
      gpstr10 <- .gpstr10
      gprobs.y  <- .gprobs.y
      ilambda <- .ilambda

      oipospois.Loglikfun <- function(pstr1, lambda, y, x, w, extraargs) {
        sum(c(w) * doipospois(x = y, pstr1 = pstr1,
                              lambda = lambda, log = TRUE))
      }


      for (jay in 1:NOS) {  # For each response 'y_jay'... do:
        TFvec <- y[, jay] > 1  # Important to exclude the 1s
        posyvec <- y[TFvec, jay]  # Variable name unchanged (lazy)
        lambd.init.jay <- if ( .imethod == 1) {
          quantile(posyvec, probs = gprobs.y) - 1/2  # + 1/16
        } else if ( .imethod == 2) {
          weighted.mean(posyvec, w = w[TFvec, jay]) - 1/2
        } else {
          warning("argument 'imethod' should have the value 1 or 2")
        }
        if (length(ilambda)) { # zz
          lambd.init.jay <- ilambda[jay]
        } else {
        }



        try.this <-
          grid.search2(gpstr10, lambd.init.jay,
                       objfun = oipospois.Loglikfun,
                       y = y[, jay],  # x = x[TFvec, , drop = FALSE],
                       w = w[, jay],
                       ret.objfun = TRUE)  # Last value is the loglik
        pstr1.init[, jay] <-  try.this["Value1"]
        lambd.init[, jay] <- (try.this["Value2"] + y[, jay]) / 2
        lambd.init[, jay] <-  try.this["Value2"]
      }  # for (jay ...)

      etastart <- cbind(theta2eta(pstr1.init, .lpstr10 , earg = .epstr10 ),
                        theta2eta(lambd.init, .llambda ,
                                  earg = .elambda ))[,
                        interleave.VGAM(M, M1 = M1)]
      mustart <- NULL  # Since etastart has been computed.
    }  # End of !length(etastart)
  }), list( .lpstr10 = lpstr10, .llambda = llambda,
            .epstr10 = epstr10, .elambda = elambda,
                                .ilambda = ilambda,
            .gpstr10 = gpstr10,
            .gprobs.y = gprobs.y,
            .imethod = imethod,  # .probs.y = probs.y,
            .type.fitted = type.fitted ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    NOS <- ncol(eta) / c(M1 = 2)
    type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "lambda", "pobs1", "pstr1", "onempstr1"))[1]

    phimat <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr10 , earg = .epstr10 )
    lambda <- eta2theta(eta[, c(FALSE, TRUE)], .llambda , earg = .elambda )

    ans <-
      switch(type.fitted,
      "mean"      = phimat - (1 - phimat) * lambda / expm1(-lambda),
      "lambda"    = lambda,
      "pobs1" = doipospois(1, lambda = lambda, pstr1 = phimat), # Pr(Y=1)
      "pstr1"     =     phimat,
      "onempstr1" = 1 - phimat)
    label.cols.y(ans, colnames.y = extra$colnames.y, NOS = NOS)
  }, list( .lpstr10 = lpstr10, .llambda = llambda,
           .epstr10 = epstr10, .elambda = elambda
         ))),
  last = eval(substitute(expression({
    misc$link <-
      c(rep_len( .lpstr10 , NOS),
        rep_len( .llambda , NOS))[interleave.VGAM(M, M1 = M1)]
    temp.names <- c(mynames1, mynames2)[interleave.VGAM(M, M1 = M1)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for (ii in 1:ncoly) {
      misc$earg[[M1*ii-1]] <- .epstr10
      misc$earg[[M1*ii  ]] <- .elambda
    }
  }), list( .lpstr10 = lpstr10, .llambda = llambda,
            .epstr10 = epstr10, .elambda = elambda ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    phimat <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr10 , earg = .epstr10 )
    lambda <- eta2theta(eta[, c(FALSE, TRUE)], .llambda , earg = .elambda )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * doipospois(x = y, pstr1 = phimat, lambda = lambda,
                                   log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lpstr10 = lpstr10, .llambda = llambda,
           .epstr10 = epstr10, .elambda = elambda ))),
  vfamily = c("oipospoisson"),



  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    phimat <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr10 , earg = .epstr10 )
    lambda <- eta2theta(eta[, c(FALSE, TRUE)], .llambda , earg = .elambda )
    roipospois(nsim * length(lambda), lambda = lambda, pstr1 = phimat)
  }, list( .lpstr10 = lpstr10, .llambda = llambda,
           .epstr10 = epstr10, .elambda = elambda ))),

  validparams = eval(substitute(function(eta, y, extra = NULL) {
    phimat <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE], .lpstr10 ,
                        earg = .epstr10 )
    lambda <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE], .llambda ,
                        earg = .elambda )
    okay1 <- all(is.finite(lambda)) && all(0 < lambda) &&
             all(is.finite(phimat)) && all(phimat < 1)
    deflat.limit <- deflat.limit.oipospois(lambda)
    okay2.deflat <- TRUE
    if (okay1 && !(okay2.deflat <- all(deflat.limit < phimat)))
      warning("parameter 'pstr1' is too negative even allowing for ",
              "0-deflation.")
    okay1 && okay2.deflat
  }, list( .lpstr10 = lpstr10, .llambda = llambda,
           .epstr10 = epstr10, .elambda = elambda ))),






  deriv = eval(substitute(expression({
    M1 <- 2
    NOS <- M / M1
    phimat <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE], .lpstr10 ,
                        earg = .epstr10 )
    lambda <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE], .llambda ,
                        earg = .elambda )

    pmf1 <- -lambda * exp(-lambda) / expm1(-lambda)
    onempmf1 <- 1 - pmf1  # doipospois(1, lambda = lambda, pstr1 = phimat)
    pobs1 <- phimat + (1 - phimat) * pmf1
    index1 <- as.matrix(y == 1)

    dl.dphimat <- onempmf1 / pobs1
    dl.dphimat[!index1] <- -1 / (1 - phimat[!index1])

    dpmf1.dlambda <- exp(-lambda) *
        (1 - lambda - exp(-lambda)) / (expm1(-lambda))^2

    d3 <- deriv3( ~ exp(-lambda) * lambda / (1 - exp(-lambda)),
                  c("lambda"), hessian = TRUE)
    eval.d3 <- eval(d3)
    d2pmf1.dlambda2 <- attr(eval.d3, "hessian")
    dim(d2pmf1.dlambda2) <- c(n, NOS)  # Matrix it, even for NOS==1

    dl.dlambda <- (1 - phimat) * dpmf1.dlambda / pobs1  #
    dl.dlambda[!index1] <- y[!index1] / lambda[!index1] - 1 -
                           1 / expm1(lambda[!index1])

    dphimat.deta <- dtheta.deta(phimat, .lpstr10 , earg = .epstr10 )
    dlambda.deta <- dtheta.deta(lambda, .llambda , earg = .elambda )

    myderiv <- c(w) * cbind(dl.dphimat * dphimat.deta,
                            dl.dlambda * dlambda.deta)
    myderiv[, interleave.VGAM(M, M1 = M1)]
  }), list( .lpstr10 = lpstr10, .llambda = llambda,
            .epstr10 = epstr10, .elambda = elambda ))),
  weight = eval(substitute(expression({

    ned2l.dphimat2 <- onempmf1 / ((1 - phimat) * pobs1)  #
    ned2l.dphimatlambda <- dpmf1.dlambda / pobs1  #
    ned2l.dlambda2 <-
     (((1 - phimat) * dpmf1.dlambda)^2) / pobs1 -
       (1 - phimat) * d2pmf1.dlambda2 +
       (1 - phimat) * (1/lambda - exp(-lambda) *
       (1 - exp(-lambda) - lambda * exp(-lambda)) / (expm1(-lambda))^3)

    wz <- array(c(c(w) * ned2l.dphimat2 * dphimat.deta^2,
                  c(w) * ned2l.dlambda2 * dlambda.deta^2,
                  c(w) * ned2l.dphimatlambda * dphimat.deta *
                                               dlambda.deta),
                dim = c(n, M / M1, 3))
    wz <- arwz2wz(wz, M = M, M1 = M1)
    wz
  }), list( .llambda = llambda, .elambda = elambda ))))
}  # oipospoisson













dposbinom <- function(x, size, prob, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  L <- max(length(x), length(size), length(prob))
  if (length(x)      != L) x    <- rep_len(x,    L)
  if (length(size)   != L) size <- rep_len(size, L)
  if (length(prob)   != L) prob <- rep_len(prob, L)

  answer <- NaN * x
  is0 <- (x == 0)
  ok2 <- (prob > 0) & (prob <= 1) &
         (size == round(size)) & (size > 0)

  answer <-        dbinom(x = x, size = size, prob = prob, log = TRUE) -
            log1p(-dbinom(x = 0, size = size, prob = prob))
  answer[!ok2] <- NaN
  if (log.arg) {
    answer[is0 & ok2] <- log(0.0)
  } else {
    answer <- exp(answer)
    answer[is0 & ok2] <- 0.0
  }
  answer
}



pposbinom <- function(q, size, prob
                     ) {


  if (!is.Numeric(prob, positive = TRUE))
    stop("no zero or non-numeric values allowed for argument 'prob'")
  L <- max(length(q), length(size), length(prob))
  if (length(q)      != L) q      <- rep_len(q,      L)
  if (length(size)   != L) size   <- rep_len(size,   L)
  if (length(prob)   != L) prob   <- rep_len(prob,   L)

  ifelse(q < 1, 0,
        (pbinom(q = q, size, prob) - dbinom(x = 0, size, prob))
       / pbinom(q = 0, size, prob, lower.tail = FALSE))
}


qposbinom <- function(p, size, prob
                     ) {



  ans <- qbinom(pbinom(0, size, prob, lower.tail = FALSE) * p +
                dbinom(0, size, prob),
                size, prob)

  ans[p == 1] <- size[p == 1]

  ans[p == 0] <- 1
  ans[prob == 0] <- NaN

  ans[p < 0] <- NaN
  ans[1 < p] <- NaN
  ans
}



rposbinom <- function(n, size, prob) {
  qbinom(p = runif(n, min = dbinom(0, size, prob)), size, prob)
}












dpospois <- function(x, lambda, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  L <- max(length(x), length(lambda))
  if (length(x)      != L) x      <- rep_len(x,      L)
  if (length(lambda) != L) lambda <- rep_len(lambda, L)

  ans <- if (log.arg) {
    ifelse(x == 0, log(0.0), dpois(x, lambda, log = TRUE) -
           log1p(-exp(-lambda)))
  } else {
    ifelse(x == 0, 0, -dpois(x, lambda) / expm1(-lambda))
  }
  ans[lambda <= 0] <- NaN
  ans
}



ppospois <- function(q, lambda) {
  L <- max(length(q), length(lambda))
  if (length(q)      != L) q      <- rep_len(q,      L)
  if (length(lambda) != L) lambda <- rep_len(lambda, L)

  ans <- ifelse(q < 1, 0, (ppois(q, lambda) - dpois(0, lambda))
                         / ppois(0, lambda, lower.tail = FALSE))
  ans[lambda <= 0] <- NaN
  ans
}



qpospois <- function(p, lambda) {


  ans <- qpois(ppois(0, lambda, lower.tail = FALSE) * p +
               dpois(0, lambda),
               lambda = lambda)

  ans[p == 1] <- Inf

  ans[p < 0] <- NaN
  ans[1 < p] <- NaN
  ans[lambda <= 0] <- NaN
  ans
}




rpospois <- function(n, lambda) {
  ans <- qpois(p = runif(n, min = dpois(0, lambda)), lambda)
  ans[lambda <= 0] <- NaN
  ans
}











prob.munb.size.VGAM <- function(munb, size) {
  prob <- size / (size + munb)
  inf.munb <- is.infinite(munb)
  inf.size <- is.infinite(size)
  prob[inf.munb] <- 0
  prob[inf.size] <- 1
  prob[inf.munb & inf.size] <- NaN
  prob[size < 0 | munb < 0] <- NaN
  prob
}




dposnegbin <-
  function(x, size, prob = NULL, munb = NULL, log = FALSE) {
  if (length(munb)) {
    if (length(prob))
      stop("arguments 'prob' and 'munb' both specified")
  } else {
    if (!length(prob))
      stop("Only one of 'prob' or 'munb' must be specified")
  }

  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  LLL <- max(length(x), length(prob), length(munb), length(size))
  if (length(x)    != LLL) x    <- rep_len(x,    LLL)
  if (length(size) != LLL) size <- rep_len(size, LLL)
  ans <- if (length(munb)) {
    if (length(munb) != LLL) munb <- rep_len(munb, LLL)
    dnbinom(x = x, size = size, mu   = munb, log = TRUE)
  } else {
    if (length(prob) != LLL) prob <- rep_len(prob, LLL)
    dnbinom(x = x, size = size, prob = prob, log = TRUE)
  }

  index0 <- (x == 0) & !is.na(size)  # & (!is.na(prob) |  !is.na(munb))
  ans[ index0] <- log(0.0)
  ans[!index0] <- ans[!index0] - (
    if (length(prob))
      pnbinom(0, size = size[!index0], prob = prob[!index0],
              lower.tail = FALSE, log.p = TRUE) else
      pnbinom(0, size = size[!index0], mu   = munb[!index0],
              lower.tail = FALSE, log.p = TRUE))


  if (!log.arg)
    ans <- exp(ans)

  if (!length(prob))
    prob <- prob.munb.size.VGAM(munb, size)
  ans[prob == 0 | prob == 1] <- NaN

  ans
}



pposnegbin <- function(q, size, prob = NULL, munb = NULL,
                       lower.tail = TRUE, log.p = FALSE) {

  if (length(munb)) {
    if (length(prob))
      stop("arguments 'prob' and 'munb' both specified")
  } else {
    if (!length(prob))
      stop("Only one of 'prob' or 'munb' must be specified")
  }

  LLL <- max(length(q), length(prob), length(munb), length(size))
  if (length(q)    != LLL) q    <- rep_len(q,    LLL)
  if (length(size) != LLL) size <- rep_len(size, LLL)
  if (length(munb)) {
    if (length(munb) != LLL) munb <- rep_len(munb, LLL)
  } else {
    if (length(prob) != LLL) prob <- rep_len(prob, LLL)
  }

  tail.prob <-
    if (length(prob)) dnbinom(0, size = size, prob = prob) else
                      dnbinom(0, size = size, mu   = munb)

  vall <- rep_len(ifelse(lower.tail, log(0), log(1)), LLL)
  ans <- if (length(prob)) {
    ifelse(q < 1,
           vall,
           (if (lower.tail)
           log(pnbinom(q, size = size, prob = prob) - tail.prob) else
           pnbinom(q, size = size, prob = prob,
                   lower.tail = FALSE, log.p = TRUE)) -
           log1p(-tail.prob))
  } else {
    ifelse(q < 1,
           vall,
           (if (lower.tail)
           log(pnbinom(q, size = size, mu   = munb) - tail.prob) else
           pnbinom(q, size = size, mu   = munb,
                       lower.tail = FALSE, log.p = TRUE)) -
           log1p(-tail.prob))
  }

  if (!log.p)
    ans <- exp(ans)

  if (!length(prob))
    prob <- prob.munb.size.VGAM(munb, size)
  ans[prob == 0 | prob == 1] <- NaN

  ans
}



qposnegbin <- function(p, size, prob = NULL, munb = NULL) {


  if (length(munb)) {
    if (length(prob))
      stop("arguments 'prob' and 'munb' both specified")
  } else {
    if (!length(prob))
      stop("Only one of 'prob' or 'munb' must be specified")
  }

  ans <- if (length(munb)) {
    qnbinom(pnbinom(0, size = size, mu   = munb,
                    lower.tail = FALSE) * p +
            dnbinom(0, size = size, mu   = munb),
            size = size, mu   = munb)
  } else {
    qnbinom(pnbinom(0, size = size, prob = prob,
                    lower.tail = FALSE) * p +
            dnbinom(0, size = size, prob = prob),
            size = size, prob = prob)
  }

  ans[p == 1] <- Inf
  ans[p < 0 | 1 < p] <- NaN

  if (!length(prob))
    prob <- prob.munb.size.VGAM(munb, size)
  ans[prob == 0 | prob == 1] <- NaN

  ans
}



rposnegbin <- function(n, size, prob = NULL, munb = NULL) {
  if (length(munb)) {
    if (length(prob))
      stop("arguments 'prob' and 'munb' both specified")
  } else {
    if (!length(prob))
      stop("Only one of 'prob' or 'munb' must be specified")
  }

  ans <- if (length(munb)) {
    qnbinom(runif(n, min = dnbinom(0, size = size, mu   = munb)),
            size = size, mu   = munb)
  } else {
    qnbinom(runif(n, min = dnbinom(0, size = size, prob = prob)),
            size = size, prob = prob)
  }
  if (!length(prob))
    prob <- prob.munb.size.VGAM(munb, size)
  ans[prob == 0 | prob == 1] <- NaN
  ans
}









dbell <- function(x, shape = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  L <- max(length(x), length(shape))
  if (length(x)     != L) x     <- rep_len(x,     L)
  if (length(shape) != L) shape <- rep_len(shape, L)

  logdensity <- rep_len(log(0), L)
  xok <- (0 <= x) & is.finite(x) & (x == round(x))
  bellnos <- bell(x[xok])
  logdensity[xok] <-
    x[xok] * log(shape[xok]) - expm1(shape[xok]) +
    log(bellnos[x[xok] + 1]) - lgamma(x[xok] + 1)
  logdensity[shape <= 0] <- NaN
  if (log.arg) logdensity else exp(logdensity)
}



rbell <- function(n, shape = 1) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
             stop("bad input for argument 'n'") else n

  shape <- rep_len(shape, use.n)
  N <- rpois(use.n, expm1(shape))  # Might have 0s
  ans <- rep_len(0, use.n)
  Nok <- !is.na(shape) & N > 0
  if (any(Nok)) {
    maxN <- max(N[Nok])
    for (kk in seq(maxN)) {
      sum.ok <- Nok & kk <= N  # [Nok]
      ans[sum.ok] <- ans[sum.ok] + rpospois(sum(sum.ok), shape[sum.ok])
    }
  }
  ans[shape <= 0] <- NaN
  ans[is.na(shape)] <- NA
  ans
}



 bellff <- function(lshape = "loglink", zero = NULL,
                    gshape = expm1(1.6 * ppoints(7))) {

  lshape <- as.list(substitute(lshape))  # orig
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")



  new("vglmff",
  blurb = c("Bell distribution\n",
            "Pr(Y=y; shape) = shape^y * exp(-expm1(shape)) * ",
            "bell(y) / y!,\n",
            "y in 0(1)Inf, 0 < shape\n",
            "Link:    ",
            namesof("shape", lshape, earg = eshape), "\n",
            "Mean:    shape * exp(shape)"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 1)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         expected = TRUE,
         hadof = TRUE,
         multipleResponses = TRUE,
         parameters.names = "shape",
         zero = .zero )
  }, list( .zero = zero ))),


  initialize = eval(substitute(expression({
    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
    if (is.infinite(bell(max(y))))
      stop("max(y) > 218, therefore bell() returns Inf")


    ncoly <- ncol(y)
    M1 <- 1
    extra$ncoly <- ncoly
    extra$M1 <- M1
    M <- M1 * ncoly


    mynames1  <- param.names("shape", ncoly, skip1 = TRUE)
    predictors.names <-
      namesof(mynames1, .lshape , earg = .eshape , tag = FALSE)


    if (!length(etastart)) {
      shape.init <- matrix(0, nrow(x), ncoly)
      gshape <- .gshape
      bellff.Loglikfun <- function(shape, y, x = NULL,
                                   w, extraargs = NULL) {
        sum(c(w) * dbell(x = y, shape = shape, log = TRUE))
      }

      for (jay in 1:ncoly) {
        shape.init[, jay] <- grid.search(gshape,
                                         objfun = bellff.Loglikfun,
                                         y = y[, jay], w = w[, jay])
      }
      etastart <- theta2eta(shape.init, .lshape , earg = .eshape )
    }
  }), list( .lshape = lshape, .gshape = gshape,
            .eshape = eshape ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    shape * exp(shape)
  }, list( .lshape = lshape,
           .eshape = eshape ))),
  last = eval(substitute(expression({
    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for (ilocal in 1:ncoly) {
      misc$earg[[ilocal]] <- .eshape
    }

    misc$link <- rep_len( .lshape , ncoly)
    names(misc$link) <- mynames1
  }), list( .lshape = lshape, .eshape = eshape ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dbell(x = y, shape = shape, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lshape = lshape, .eshape = eshape ))),
  vfamily = c("bellff"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    okay1 <- all(is.finite(shape)) && all(0 < shape)
    okay1
  }, list( .lshape = lshape, .eshape = eshape ))),




  hadof = eval(substitute(
  function(eta, extra = list(), deriv = 1,
           linpred.index = 1,
           w = 1, dim.wz = c(NROW(eta), NCOL(eta) * (NCOL(eta)+1)/2),
           ...) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    ans <- c(w) *
    switch(as.character(deriv),
           "0" =  exp(shape) * (1 + 1 / shape),
           "1" =  exp(shape) * (1 + 1 / shape - 1 / shape^2),
           "2" =  exp(shape) * (1 + 1 / shape - 2 / shape^2 + 2 / shape^3),
           "3" =  exp(shape) * (1 + 1 / shape - 3 / shape^2 + 6 / shape^3 -
                                                              6 / shape^4),
           stop("argument 'deriv' must be 0, 1, 2 or 3"))
    if (deriv == 0) ans else retain.col(ans, linpred.index)  # M1 = 1
  }, list( .lshape = lshape, .eshape = eshape ))),




  simslot = eval(substitute(
  function(object, nsim) {
 
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    rbell(nsim * length(shape), shape = c(shape))
  }, list( .lshape = lshape,
           .eshape = eshape ))),


  deriv = eval(substitute(expression({
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    dl.dshape <- y / shape - exp(shape)
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )
    c(w) * dl.dshape * dshape.deta
  }), list( .lshape = lshape, .eshape = eshape ))),
  weight = eval(substitute(expression({
    ned2l.dshape2 <- (1 + shape) * exp(shape) / shape
    wz <- c(w) * ned2l.dshape2 * dshape.deta^2
    wz
  }), list( .lshape = lshape, .eshape = eshape ))))
}  # bellff


















