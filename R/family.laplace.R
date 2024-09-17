# These functions are
# Copyright (C) 1998-2024 T.W. Yee, University of Auckland.
# All rights reserved.

  










rho1check <- function(u, tau = 0.5)
  u * (tau - (u <= 0))




dalap <- function(x, location = 0, scale = 1, tau = 0.5,
                  kappa = sqrt(tau/(1-tau)), log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)




    NN <- max(length(x), length(location),
              length(scale), length(kappa),
            length(tau))
  if (length(x)        != NN) x        <- rep_len(x,        NN)
  if (length(location) != NN) location <- rep_len(location, NN)
  if (length(scale)    != NN) scale    <- rep_len(scale,    NN)
  if (length(kappa)    != NN) kappa    <- rep_len(kappa,    NN)
  if (length(tau)      != NN) tau      <- rep_len(tau,      NN)

    logconst <- 0.5 * log(2) - log(scale) +
        log(kappa) - log1p(kappa^2)
  exponent <- -(sqrt(2) / scale) * abs(x - location) *
             ifelse(x >= location, kappa, 1/kappa)

  indexTF <- (scale > 0) & (tau > 0) & (tau < 1) & (kappa > 0)  # &
  logconst[!indexTF] <- NaN

  if (log.arg) logconst + exponent else exp(logconst + exponent)
}


ralap <- function(n, location = 0, scale = 1, tau = 0.5,
                  kappa = sqrt(tau/(1-tau))) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n

  location <- rep_len(location, use.n)
  scale    <- rep_len(scale,    use.n)
  tau      <- rep_len(tau,      use.n)
  kappa    <- rep_len(kappa,    use.n)
  ans <- location + scale *
        log(runif(use.n)^kappa / runif(use.n)^(1/kappa)) / sqrt(2)
  indexTF <- (scale > 0) & (tau > 0) & (tau < 1) & (kappa > 0)  # &
  ans[!indexTF] <- NaN
  ans
}



palap <- function(q, location = 0, scale = 1, tau = 0.5,
                  kappa = sqrt(tau/(1-tau)),
                  lower.tail = TRUE, log.p = FALSE) {

  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

    NN <- max(length(q), length(location),
              length(scale), length(kappa),
            length(tau))
  if (length(q)        != NN) q        <- rep_len(q,        NN)
  if (length(location) != NN) location <- rep_len(location, NN)
  if (length(scale)    != NN) scale    <- rep_len(scale,    NN)
  if (length(kappa)    != NN) kappa    <- rep_len(kappa,    NN)
  if (length(tau)      != NN) tau      <- rep_len(tau,      NN)

  exponent <- -(sqrt(2) / scale) * abs(q - location) *
              ifelse(q >= location, kappa, 1/kappa)
  temp5 <- exp(exponent) / (1 + kappa^2)
  index1 <- (q < location)


  if (lower.tail) {
    if (log.p) {
      ans <- log1p(-exp(exponent) / (1 + kappa^2))
      logtemp5 <- exponent - log1p(kappa^2)
      ans[index1] <- 2 * log(kappa[index1]) + logtemp5[index1]
    } else {
      ans <- (kappa^2 - expm1(exponent)) / (1 + kappa^2)
      ans[index1] <- (kappa[index1])^2 * temp5[index1]
    }
  } else {
    if (log.p) {
      ans <- exponent - log1p(kappa^2)  # logtemp5
      ans[index1] <- log1p(-(kappa[index1])^2 * temp5[index1])
    } else {
      ans <- temp5
      ans[index1] <- (1 + (kappa[index1])^2 *
                      (-expm1(exponent[index1]))) / (1+
                    (kappa[index1])^2)
      }
  }
  indexTF <- (scale > 0) & (tau > 0) & (tau < 1) & (kappa > 0)  # &
  ans[!indexTF] <- NaN
  ans
}



qalap <- function(p, location = 0, scale = 1, tau = 0.5,
                  kappa = sqrt(tau / (1 - tau)),
                  lower.tail = TRUE, log.p = FALSE) {

  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

    NN <- max(length(p), length(location),
              length(scale), length(kappa),
            length(tau))
  if (length(p)        != NN) p        <- rep_len(p,        NN)
  if (length(location) != NN) location <- rep_len(location, NN)
  if (length(scale)    != NN) scale    <- rep_len(scale,    NN)
  if (length(kappa)    != NN) kappa    <- rep_len(kappa,    NN)
  if (length(tau)      != NN) tau      <- rep_len(tau,      NN)



  temp5 <- kappa^2 / (1 + kappa^2)
  if (lower.tail) {
    if (log.p) {
      ans <- exp(p)
      index1 <- (exp(p) <= temp5)
      exponent <- exp(p[index1]) / temp5[index1]
      ans[index1] <- location[index1] +
          (scale[index1] * kappa[index1]) *
        log(exponent) / sqrt(2)
      ans[!index1] <- location[!index1] -
          (scale[!index1] / kappa[!index1]) *
        (log1p((kappa[!index1])^2) +
           log(-expm1(p[!index1]))) / sqrt(2)
    } else {
      ans <- p
      index1 <- (p <= temp5)
      exponent <- p[index1] / temp5[index1]
        ans[index1] <- location[index1] +
            (scale[index1] * kappa[index1]) *
        log(exponent) / sqrt(2)
      ans[!index1] <- location[!index1] -
          (scale[!index1] / kappa[!index1]) *
                      (log1p((kappa[!index1])^2) +
                      log1p(-p[!index1])) / sqrt(2)
    }
  } else {
    if (log.p) {
      ans <- -expm1(p)
      index1 <- (-expm1(p)  <= temp5)
      exponent <- -expm1(p[index1]) / temp5[index1]
      ans[index1] <- location[index1] +
          (scale[index1] * kappa[index1]) *
        log(exponent) / sqrt(2)
      ans[!index1] <- location[!index1] -
          (scale[!index1] / kappa[!index1]) *
        (log1p((kappa[!index1])^2) +
           p[!index1]) / sqrt(2)
    } else {
      ans <- exp(log1p(-p))
      index1 <- (p >= (1 / (1+kappa^2)))
      exponent <- exp(log1p(-p[index1])) / temp5[index1]
      ans[index1] <- location[index1] +
          (scale[index1] * kappa[index1]) *
          log(exponent) / sqrt(2)
      ans[!index1] <- location[!index1] -
          (scale[!index1] / kappa[!index1]) *
        (log1p((kappa[!index1])^2) +
           log(p[!index1])) / sqrt(2)
    }
  }

    indexTF <- (scale > 0) & (tau > 0) &
        (tau < 1) & (kappa > 0)  # &
  ans[!indexTF] <- NaN
  ans
}







rloglap <- function(n, location.ald = 0, scale.ald = 1, tau = 0.5,
                    kappa = sqrt(tau/(1-tau))) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
            stop("bad input for argument 'n'") else n
  location.ald <- rep_len(location.ald, use.n)
  scale.ald    <- rep_len(scale.ald,    use.n)
  tau          <- rep_len(tau,          use.n)
  kappa        <- rep_len(kappa,        use.n)
  ans <- exp(location.ald) *
      (runif(use.n)^kappa / runif(use.n)^(1/kappa))^(scale.ald
          / sqrt(2))
    indexTF <- (scale.ald > 0) & (tau > 0) &
        (tau < 1) & (kappa > 0)  # &
  ans[!indexTF] <- NaN
  ans
}



dloglap <-
    function(x, location.ald = 0, scale.ald = 1, tau = 0.5,
                    kappa = sqrt(tau/(1-tau)), log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  scale    <- scale.ald
  location <- location.ald
  NN <- max(length(x), length(location),
           length(scale), length(kappa), length(tau))

  if (length(x)        != NN) x        <- rep_len(x,        NN)
  if (length(location) != NN) location <- rep_len(location, NN)
  if (length(scale)    != NN) scale    <- rep_len(scale,    NN)
  if (length(kappa)    != NN) kappa    <- rep_len(kappa,    NN)
  if (length(tau)      != NN) tau      <- rep_len(tau,      NN)


  Alpha <- sqrt(2) * kappa / scale.ald
  Beta  <- sqrt(2) / (scale.ald * kappa)
  Delta <- exp(location.ald)
  exponent <- ifelse(x >= Delta, -(Alpha+1), (Beta-1)) *
             (log(x) - location.ald)
  logdensity <- -location.ald + log(Alpha) + log(Beta) -
               log(Alpha + Beta) + exponent
        indexTF <- (scale.ald > 0) & (tau > 0) &
            (tau < 1) & (kappa > 0)  # &
  logdensity[!indexTF] <- NaN
  logdensity[x <  0 & indexTF] <- -Inf
  if (log.arg) logdensity else exp(logdensity)
}



qloglap <- function(p, location.ald = 0, scale.ald = 1,
                    tau = 0.5, kappa = sqrt(tau/(1-tau)),
                    lower.tail = TRUE, log.p = FALSE) {

  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")


  NN <- max(length(p), length(location.ald), length(scale.ald),
            length(kappa))
  p        <- rep_len(p,            NN)
  location <- rep_len(location.ald, NN)
  scale    <- rep_len(scale.ald,    NN)
  kappa    <- rep_len(kappa,        NN)
  tau      <- rep_len(tau,          NN)


  Alpha <- sqrt(2) * kappa / scale.ald
  Beta  <- sqrt(2) / (scale.ald * kappa)
  Delta <- exp(location.ald)
  temp9 <- Alpha + Beta


  if (lower.tail) {
    if (log.p) {
      ln.p <- p
      ans <- ifelse((exp(ln.p) > Alpha / temp9),
               Delta * (-expm1(ln.p) * temp9 / Beta)^(-1/Alpha),
               Delta * (exp(ln.p) * temp9 / Alpha)^(1/Beta))
      ans[ln.p > 0] <- NaN
    } else {
      ans <- ifelse((p > Alpha / temp9),
                    Delta * exp((-1/Alpha) * (log1p(-p) +
                                              log(temp9/Beta))),
                    Delta * (p * temp9 / Alpha)^(1/Beta))
      ans[p <  0] <- NaN
      ans[p == 0] <- 0
      ans[p == 1] <- Inf
      ans[p >  1] <- NaN
    }
  } else {
    if (log.p) {
      ln.p <- p
      ans <- ifelse((-expm1(ln.p) > Alpha / temp9),
                    Delta * (exp(ln.p) * temp9 / Beta)^(-1/Alpha),
                    Delta * (-expm1(ln.p) * temp9 / Alpha)^(1/Beta))
      ans[ln.p > 0] <- NaN
    } else {
      ans <- ifelse((p < (temp9 - Alpha) / temp9),
                    Delta * (p * temp9 / Beta)^(-1/Alpha),
           Delta * exp((1/Beta)*(log1p(-p) + log(temp9/Alpha))))
      ans[p <  0] <- NaN
      ans[p == 0] <- Inf
      ans[p == 1] <- 0
      ans[p >  1] <- NaN
    }
  }
  indexTF <- (scale.ald > 0) & (tau > 0) & (tau < 1) & (kappa > 0)
  ans[!indexTF] <- NaN
  ans
}



ploglap <- function(q, location.ald = 0, scale.ald = 1,
                    tau = 0.5, kappa = sqrt(tau/(1-tau)),
                    lower.tail = TRUE, log.p = FALSE) {

  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  NN <- max(length(q), length(location.ald), length(scale.ald),
            length(kappa))
  location <- rep_len(location.ald, NN)
  scale    <- rep_len(scale.ald,    NN)
  kappa    <- rep_len(kappa,        NN)
  q        <- rep_len(q,            NN)
  tau      <- rep_len(tau,          NN)

  Alpha <- sqrt(2) * kappa / scale.ald
  Beta  <- sqrt(2) / (scale.ald * kappa)
  Delta <- exp(location.ald)

  temp9 <- Alpha + Beta
  index1 <- (Delta <= q)


  if (lower.tail) {
    if (log.p) {
      ans <- log((Alpha / temp9) * (q / Delta)^(Beta))
        ans[index1] <- log1p((-(Beta/temp9) *
                              (Delta/q)^(Alpha))[index1])
      ans[q <= 0 ] <- -Inf
      ans[q == Inf] <- 0
    } else {
      ans <- (Alpha / temp9) * (q / Delta)^(Beta)
      ans[index1] <- -expm1((log(Beta/temp9) +
                             Alpha * log(Delta/q)))[index1]
      ans[q <= 0] <- 0
      ans[q == Inf] <- 1
    }
  } else {
    if (log.p) {
      ans <- log1p(-(Alpha / temp9) * (q / Delta)^(Beta))
      ans[index1] <- log(((Beta/temp9) * (Delta/q)^(Alpha))[index1])
      ans[q <= 0] <- 0
      ans[q == Inf] <- -Inf
    } else {
      ans <- -expm1(log(Alpha/temp9) + Beta * log(q/Delta))
      ans[index1] <- ((Beta/temp9) * (Delta/q)^(Alpha))[index1]
      ans[q <= 0] <- 1
      ans[q == Inf] <- 0
    }
  }

    indexTF <- (scale.ald > 0) & (tau > 0) &
        (tau < 1) & (kappa > 0)  # &
  ans[!indexTF] <- NaN
  ans
}





rlogitlap <-
    function(n, location.ald = 0, scale.ald = 1, tau = 0.5,
                      kappa = sqrt(tau/(1-tau))) {
        logitlink(ralap(n = n, location = location.ald,
                        scale = scale.ald,
              tau = tau, kappa = kappa),
        inverse = TRUE)  # earg = earg
}



dlogitlap <-
    function(x, location.ald = 0, scale.ald = 1, tau = 0.5,
                      kappa = sqrt(tau/(1-tau)), log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)



  NN <- max(length(x), length(location.ald),
           length(scale.ald), length(kappa))
  location <- rep_len(location.ald, NN)
  scale    <- rep_len(scale.ald,    NN)
  kappa    <- rep_len(kappa,        NN)
  x        <- rep_len(x,            NN)
  tau      <- rep_len(tau,          NN)

  Alpha <- sqrt(2) * kappa / scale.ald
  Beta  <- sqrt(2) / (scale.ald * kappa)
  Delta <- logitlink(location.ald, inverse = TRUE)  # earg = earg

  exponent <- ifelse(x >= Delta, -Alpha, Beta) *
             (logitlink(x) - # earg = earg
              location.ald)
  logdensity <- log(Alpha) + log(Beta) - log(Alpha + Beta) -
               log(x) - log1p(-x) + exponent
        indexTF <- (scale.ald > 0) & (tau > 0) &
            (tau < 1) & (kappa > 0)  # &
  logdensity[!indexTF] <- NaN
  logdensity[x <  0 & indexTF] <- -Inf
  logdensity[x >  1 & indexTF] <- -Inf
  if (log.arg) logdensity else exp(logdensity)
}



qlogitlap <-
    function(p, location.ald = 0, scale.ald = 1,
                      tau = 0.5, kappa = sqrt(tau/(1-tau))) {
  qqq <- qalap(p = p, location = location.ald, scale = scale.ald,
              tau = tau, kappa = kappa)
  ans <- logitlink(qqq, inverse = TRUE)  # earg = earg
  ans[(p < 0) | (p > 1)] <- NaN
  ans[p == 0] <- 0
  ans[p == 1] <- 1
  ans
}



plogitlap <-
    function(q, location.ald = 0, scale.ald = 1,
                      tau = 0.5, kappa = sqrt(tau/(1-tau))) {
  NN <- max(length(q), length(location.ald), length(scale.ald),
            length(kappa))
  location.ald <- rep_len(location.ald, NN)
  scale.ald    <- rep_len(scale.ald,    NN)
  kappa        <- rep_len(kappa,        NN)
  q            <- rep_len(q,            NN)
  tau          <- rep_len(tau,          NN)

  indexTF <- (q > 0) & (q < 1)
  qqq <- logitlink(q[indexTF])  # earg = earg
  ans <- q
  ans[indexTF] <- palap(q = qqq, location = location.ald[indexTF],
                       scale = scale.ald[indexTF],
                       tau = tau[indexTF], kappa = kappa[indexTF])
  ans[q >= 1] <- 1
  ans[q <= 0] <- 0
  ans
}





rprobitlap <-
    function(n, location.ald = 0, scale.ald = 1, tau = 0.5,
                       kappa = sqrt(tau/(1-tau))) {



        probitlink(ralap(n = n, location = location.ald,
                         scale = scale.ald,
               tau = tau, kappa = kappa),
               inverse = TRUE)
}



dprobitlap <-
  function(x, location.ald = 0, scale.ald = 1, tau = 0.5,
           kappa = sqrt(tau/(1-tau)), log = FALSE,
           meth2 = TRUE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)



  NN <- max(length(x), length(location.ald), length(scale.ald),
           length(kappa))
  location.ald <- rep_len(location.ald, NN)
  scale.ald    <- rep_len(scale.ald,    NN)
  kappa        <- rep_len(kappa,        NN)
  x            <- rep_len(x,            NN)
  tau          <- rep_len(tau,          NN)

  logdensity <- x * NaN
  index1 <- (x > 0) & (x < 1)
      indexTF <- (scale.ald > 0) & (tau > 0) &
          (tau < 1) & (kappa > 0)  # &
  if (meth2) {
    dx.dy <- x
    use.x <- probitlink(x[index1])  # earg = earg
    logdensity[index1] <-
      dalap(x = use.x, location = location.ald[index1],
            scale = scale.ald[index1], tau = tau[index1],
            kappa = kappa[index1], log = TRUE)
  } else {
    Alpha <- sqrt(2) * kappa / scale.ald
    Beta  <- sqrt(2) / (scale.ald * kappa)
    Delta <- pnorm(location.ald)
    use.x  <- qnorm(x)  # qnorm(x[index1])
    log.dy.dw <- dnorm(use.x, log = TRUE)

    exponent <- ifelse(x >= Delta, -Alpha, Beta) *
                     (use.x - location.ald) - log.dy.dw

    logdensity[index1] <- (log(Alpha) + log(Beta) -
                          log(Alpha + Beta) + exponent)[index1]
  }
  logdensity[!indexTF] <- NaN
  logdensity[x <  0 & indexTF] <- -Inf
  logdensity[x >  1 & indexTF] <- -Inf

  if (meth2) {
    dx.dy[index1] <- probitlink(x[index1],  # earg = earg,
                            inverse = TRUE,
                            deriv = 1)
    dx.dy[!index1] <- 0
    dx.dy[!indexTF] <- NaN
    if (log.arg) logdensity - log(abs(dx.dy)) else
                 exp(logdensity) / abs(dx.dy)
  } else {
    if (log.arg) logdensity else exp(logdensity)
  }
}


qprobitlap <-
    function(p, location.ald = 0, scale.ald = 1,
                       tau = 0.5, kappa = sqrt(tau/(1-tau))) {
  qqq <- qalap(p = p, location = location.ald, scale = scale.ald,
              tau = tau, kappa = kappa)
  ans <- probitlink(qqq, inverse = TRUE)  # , earg = earg
  ans[(p < 0) | (p > 1)] = NaN
  ans[p == 0] <- 0
  ans[p == 1] <- 1
  ans
}



pprobitlap <-
    function(q, location.ald = 0, scale.ald = 1,
                       tau = 0.5, kappa = sqrt(tau/(1-tau))) {
  NN <- max(length(q), length(location.ald), length(scale.ald),
            length(kappa))
  location.ald <- rep_len(location.ald, NN)
  scale.ald    <- rep_len(scale.ald,    NN)
  kappa        <- rep_len(kappa,        NN)
  q            <- rep_len(q,            NN)
  tau          <- rep_len(tau,          NN)

  indexTF <- (q > 0) & (q < 1)
  qqq <- probitlink(q[indexTF])  # earg = earg
  ans <- q
  ans[indexTF] <- palap(q = qqq, location = location.ald[indexTF],
                       scale = scale.ald[indexTF],
                       tau = tau[indexTF], kappa = kappa[indexTF])
  ans[q >= 1] <- 1
  ans[q <= 0] <- 0
  ans
}





rclogloglap <-
    function(n, location.ald = 0, scale.ald = 1, tau = 0.5,
                        kappa = sqrt(tau/(1-tau))) {
        clogloglink(ralap(n = n, location = location.ald,
                          scale = scale.ald,
                    tau = tau, kappa = kappa),  # earg = earg,
          inverse = TRUE)
}



dclogloglap <-
    function(x, location.ald = 0, scale.ald = 1, tau = 0.5,
             kappa = sqrt(tau/(1-tau)), log = FALSE,
             meth2 = TRUE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)



  NN <- max(length(x), length(location.ald), length(scale.ald),
           length(kappa))
  location.ald <- rep_len(location.ald, NN)
  scale.ald    <- rep_len(scale.ald,    NN)
  kappa        <- rep_len(kappa,        NN)
  x            <- rep_len(x,            NN)
  tau          <- rep_len(tau,          NN)

  logdensity <- x * NaN
  index1 <- (x > 0) & (x < 1)
        indexTF <- (scale.ald > 0) & (tau > 0) &
            (tau < 1) & (kappa > 0)  # &
  if (meth2) {
    dx.dy <- x
    use.w <- clogloglink(x[index1])  # earg = earg
    logdensity[index1] <-
      dalap(x = use.w, location = location.ald[index1],
            scale = scale.ald[index1],
            tau = tau[index1],
            kappa = kappa[index1], log = TRUE)

  } else {
    Alpha <- sqrt(2) * kappa / scale.ald
    Beta  <- sqrt(2) / (scale.ald * kappa)
    Delta <- clogloglink(location.ald, inverse = TRUE)

    exponent <- ifelse(x >= Delta, -(Alpha+1), Beta-1) *
        log(-log1p(-x)) +
               ifelse(x >= Delta, Alpha, -Beta) * location.ald
    logdensity[index1] <- (log(Alpha) + log(Beta) -
                           log(Alpha + Beta) -
                           log1p(-x) + exponent)[index1]
  }
  logdensity[!indexTF] <- NaN
  logdensity[x <  0 & indexTF] <- -Inf
  logdensity[x >  1 & indexTF] <- -Inf

  if (meth2) {
    dx.dy[index1] <- clogloglink(x[index1],  # earg = earg,
                             inverse = TRUE, deriv = 1)
    dx.dy[!index1] <- 0
    dx.dy[!indexTF] <- NaN
    if (log.arg) logdensity - log(abs(dx.dy)) else
                 exp(logdensity) / abs(dx.dy)
  } else {
    if (log.arg) logdensity else exp(logdensity)
  }
}



qclogloglap <-
    function(p, location.ald = 0, scale.ald = 1,
                        tau = 0.5, kappa = sqrt(tau/(1-tau))) {
  qqq <- qalap(p = p, location = location.ald, scale = scale.ald,
              tau = tau, kappa = kappa)
  ans <- clogloglink(qqq, inverse = TRUE)  # , earg = earg
  ans[(p < 0) | (p > 1)] <- NaN
  ans[p == 0] <- 0
  ans[p == 1] <- 1
  ans
}



pclogloglap <-
    function(q, location.ald = 0, scale.ald = 1,
                        tau = 0.5, kappa = sqrt(tau/(1-tau))) {
  NN <- max(length(q), length(location.ald), length(scale.ald),
            length(kappa))
  location.ald <- rep_len(location.ald, NN)
  scale.ald    <- rep_len(scale.ald,    NN)
  kappa        <- rep_len(kappa,        NN)
  q            <- rep_len(q,            NN)
  tau          <- rep_len(tau,          NN)


  indexTF <- (q > 0) & (q < 1)
  qqq <- clogloglink(q[indexTF])  # earg = earg
  ans <- q
  ans[indexTF] <- palap(q = qqq, location = location.ald[indexTF],
                       scale = scale.ald[indexTF],
                       tau = tau[indexTF], kappa = kappa[indexTF])
  ans[q >= 1] <- 1
  ans[q <= 0] <- 0
  ans
}












alaplace2.control <- function(maxit = 100, ...) {
  list(maxit = maxit)
}


 alaplace2 <-
  function(tau = NULL,
           llocation = "identitylink", lscale = "loglink",
           ilocation = NULL,           iscale = NULL,
           kappa = sqrt(tau / (1-tau)),
           ishrinkage = 0.95,

           parallel.locat = TRUE  ~ 0,

           parallel.scale = FALSE ~ 0,

           digt = 4,
           idf.mu = 3,
           imethod = 1,
           zero = "scale") {



  apply.parint.locat <- FALSE
  apply.parint.scale <- TRUE




  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")

  ilocat <- ilocation



  if (!is.Numeric(kappa, positive = TRUE))
    stop("bad input for argument 'kappa'")
  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
    imethod > 4)
    stop("argument 'imethod' must be 1, 2 or ... 4")
  if (length(iscale) &&
      !is.Numeric(iscale, positive = TRUE))
    stop("bad input for argument 'iscale'")
  if (!is.Numeric(ishrinkage, length.arg = 1) ||
    ishrinkage < 0 ||
    ishrinkage > 1)
    stop("bad input for argument 'ishrinkage'")
  if (length(tau) &&
      max(abs(kappa - sqrt(tau / (1 - tau)))) > 1.0e-6)
    stop("arguments 'kappa' and 'tau' do not match")



  fittedMean <- FALSE
  if (!is.logical(fittedMean) || length(fittedMean) != 1)
    stop("bad input for argument 'fittedMean'")



  new("vglmff",
  blurb = c("Two-parameter asymmetric Laplace distribution\n\n",
            "Links:      ",
            namesof("location",  llocat, earg = elocat), ", ",
            namesof("scale",     lscale, earg = escale),  # ", ",
            "\n\n",
            "Mean:       ",
            "location + scale * (1/kappa - kappa) / sqrt(2)", "\n",
            "Quantiles:  location", "\n",
            "Variance:   scale^2 * (1 + kappa^4) / (2 * kappa^2)"),




  constraints = eval(substitute(expression({


    onemat <- matrix(1, Mdiv2, 1)
    constraints.orig <- constraints


    cm1.locat <- kronecker(diag(Mdiv2), rbind(1, 0))
    cmk.locat <- kronecker(onemat,      rbind(1, 0))
    con.locat <- cm.VGAM(cmk.locat,
                         x = x, bool = .parallel.locat ,
                         constraints = constraints.orig,
                         apply.int = .apply.parint.locat ,
                         cm.default           = cm1.locat,
                         cm.intercept.default = cm1.locat)



    cm1.scale <- kronecker(diag(Mdiv2), rbind(0, 1))
    cmk.scale <- kronecker(onemat,      rbind(0, 1))
    con.scale <- cm.VGAM(cmk.scale,
                         x = x, bool = .parallel.scale ,
                         constraints = constraints.orig,
                         apply.int = .apply.parint.scale ,
                         cm.default           = cm1.scale,
                         cm.intercept.default = cm1.scale)

    con.use <- con.scale
    for (klocal in seq_along(con.scale)) {
      con.use[[klocal]] <- cbind(con.locat[[klocal]],
                                 con.scale[[klocal]])
    }


    constraints <- con.use

    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = M1)
  }), list( .parallel.locat = parallel.locat,
            .parallel.scale = parallel.scale,
            .zero = zero,
            .apply.parint.scale = apply.parint.scale,
            .apply.parint.locat = apply.parint.locat ))),


  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         summary.pvalues = FALSE,
         expected = TRUE,   # 20161117
         multipleResponses = TRUE,  # FALSE,
         parameters.names = c("location", "scale"),
         true.mu = .fittedMean ,
         zero = .zero ,
         tau  = .tau ,
         kappa = .kappa )
  }, list( .tau   = tau,
           .kappa = kappa,
           .fittedMean = fittedMean,
           .zero = zero ))),


  initialize = eval(substitute(expression({
    M1 <- 2

    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = if (length( .kappa ) > 1) 1 else Inf,
              ncol.y.max = if (length( .kappa ) > 1) 1 else Inf,
              out.wy = TRUE,
              colsyperw = 1,  # Uncommented out 20140621
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    extra$ncoly <- ncoly <- ncol(y)
    if ((ncoly > 1) && (length( .kappa ) > 1))
      stop("response must be a vector if 'kappa' or 'tau' ",
           "has a length greater than one")



    extra$kappa <- .kappa
    extra$tau <- extra$kappa^2 / (1 + extra$kappa^2)

    extra$Mdiv2 <- Mdiv2 <- max(ncoly, length( .kappa ))
    extra$M <- M <- M1 * Mdiv2
    extra$n <- n



    extra$tau.names <- tau.names <-
      paste0("(tau = ", round(extra$tau, digits = .digt), ")")
    extra$Y.names <- Y.names <- if (ncoly > 1)
                                    dimnames(y)[[2]] else "y"
    if (is.null(Y.names) || any(Y.names == ""))
      extra$Y.names <- Y.names <- paste0("y", 1:ncoly)
    extra$y.names <- y.names <-
      if (ncoly > 1) paste0(Y.names, tau.names) else tau.names

    extra$individual <- FALSE


    mynames1 <- param.names("location", Mdiv2, skip1 = TRUE)
    mynames2 <- param.names("scale",    Mdiv2, skip1 = TRUE)
    predictors.names <-
        c(namesof(mynames1, .llocat , .elocat , tag = FALSE),
          namesof(mynames2, .lscale , .escale , tag = FALSE))
    predictors.names <-
    predictors.names[interleave.VGAM(M, M1 = M1)]




    locat.init <- scale.init <- matrix(0, n, Mdiv2)
    if (!length(etastart)) {
      for (jay in 1:Mdiv2) {
        y.use <- if (ncoly > 1) y[, jay] else y
        Jay   <- if (ncoly > 1) jay else 1
        if ( .imethod == 1) {
          locat.init[, jay] <- weighted.mean(y.use, w[, Jay])
          scale.init[, jay] <- sqrt(var(y.use) / 2)
        } else if ( .imethod == 2) {
          locat.init[, jay] <- median(y.use)
          scale.init[, jay] <- sqrt(sum(c(w[, Jay]) *
             abs(y - median(y.use))) / (sum(w[, Jay]) * 2))
        } else if ( .imethod == 3) {
          Fit5 <- vsmooth.spline(x = x[, min(ncol(x), 2)],
                                 y = y.use, w = w[, Jay],
                                 df = .idf.mu )
            locat.init[, jay] <- predict(Fit5,
                                         x = x[, min(ncol(x), 2)])$y
          scale.init[, jay] <- sqrt(sum(c(w[, Jay]) *
                                    abs(y.use - median(y.use))) / (
                                        sum(w[, Jay]) * 2))
        } else {
          use.this <- weighted.mean(y.use, w[, Jay])
          locat.init[, jay] <- (1 - .ishrinkage ) * y.use +
                                    .ishrinkage   * use.this
          scale.init[, jay] <-
            sqrt(sum(c(w[, Jay]) *
            abs(y.use - median(y.use ))) / (sum(w[, Jay]) * 2))
        }
      }



      if (length( .ilocat )) {
        locat.init <- matrix( .ilocat , n, Mdiv2, byrow = TRUE)
      }
      if (length( .iscale )) {
        scale.init <- matrix( .iscale , n, Mdiv2, byrow = TRUE)
      }

      etastart <-
          cbind(theta2eta(locat.init, .llocat , earg = .elocat ),
                theta2eta(scale.init, .lscale , earg = .escale ))
        etastart <- etastart[, interleave.VGAM(M, M1 = M1),
                             drop = FALSE]
    }
  }), list( .imethod = imethod,
            .idf.mu = idf.mu,
            .ishrinkage = ishrinkage, .digt = digt,
            .elocat = elocat, .escale = escale,
            .llocat = llocat, .lscale = lscale, .kappa = kappa,
            .ilocat = ilocat, .iscale = iscale ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    M1 <- 2
    Mdiv2 <- ncol(eta) / M1  # extra$Mdiv2
    vTF <- c(TRUE, FALSE)
    locat <- eta2theta(eta[,  vTF, drop = FALSE], .llocat ,
                       earg = .elocat )
    dimnames(locat) <- list(dimnames(eta)[[1]], extra$y.names)
    myans <- if ( .fittedMean ) {
      kappamat <- matrix(extra$kappa, extra$n, extra$Mdiv2,
                         byrow = TRUE)
      Scale <- eta2theta(eta[, !vTF, drop = FALSE], .lscale ,
                         earg = .escale )
      locat + Scale * (1/kappamat - kappamat)
    } else {
      locat
    }
    dimnames(myans) <- list(dimnames(myans)[[1]], extra$y.names)
    myans
  }, list( .elocat = elocat, .llocat = llocat,
           .escale = escale, .lscale = lscale,
           .fittedMean = fittedMean,
           .kappa = kappa ))),
  last = eval(substitute(expression({
    M1 <- 2  # extra$M1
    Mdiv2 <- ncol(eta) / M1  # extra$Mdiv2



    misc$link <- setNames(c(rep_len( .llocat , Mdiv2),
                            rep_len( .lscale , Mdiv2)),
                          c(mynames1, mynames2))[interleave.VGAM(M,
                                                   M1 = M1)]



    misc$earg <- vector("list", M)
    for (ii in 1:Mdiv2) {
      misc$earg[[M1 * ii - 1]] <- .elocat
      misc$earg[[M1 * ii    ]] <- .escale
    }
    names(misc$earg) <- names(misc$link)


    extra$kappa <- misc$kappa <- .kappa
    extra$tau <- misc$tau <- misc$kappa^2 / (1 + misc$kappa^2)

    extra$percentile <- numeric(Mdiv2)  # length(misc$kappa)
    locat <- as.matrix(locat)
    for (ii in 1:Mdiv2) {
      y.use <- if (ncoly > 1) y[, ii] else y
      Jay   <- if (ncoly > 1) ii else 1
      extra$percentile[ii] <- 100 *
          weighted.mean(y.use <= locat[, ii],
                                                  w[, Jay])
    }
    names(extra$percentile) <- y.names
  }), list( .elocat = elocat, .llocat = llocat,
            .escale = escale, .lscale = lscale,
            .fittedMean = fittedMean,
            .kappa = kappa ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    M1 <- 2
    Mdiv2 <- ncol(eta) / M1  # extra$Mdiv2
    ymat <- matrix(y, extra$n, extra$Mdiv2)
    kappamat <- matrix(extra$kappa, extra$n, extra$Mdiv2,
                       byrow = TRUE)

    vTF <- c(TRUE, FALSE)
    locat <- eta2theta(eta[,  vTF, drop = FALSE], .llocat ,
                       earg = .elocat )
    Scale <- eta2theta(eta[, !vTF, drop = FALSE], .lscale ,
                       earg = .escale )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dalap(x = c(ymat), location = c(locat),
                              scale = c(Scale), kappa = c(kappamat),
                              log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .elocat = elocat, .llocat = llocat,
           .escale = escale, .lscale = lscale,
           .kappa = kappa ))),
  vfamily = c("alaplace2"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    vTF <- c(TRUE, FALSE)
    locat <- eta2theta(eta[,  vTF, drop = FALSE], .llocat ,
                       earg = .elocat )
    Scale <- eta2theta(eta[, !vTF, drop = FALSE], .lscale ,
                       earg = .escale )
    okay1 <- all(is.finite(locat)) &&
             all(is.finite(Scale)) && all(0 < Scale)
    okay1
  }, list( .elocat = elocat, .llocat = llocat,
           .escale = escale, .lscale = lscale,
           .kappa = kappa ))),




  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    extra    <- object@extra
    vTF <- c(TRUE, FALSE)
 locat <- eta2theta(eta[,  vTF, drop = FALSE], .llocat , .elocat )
 Scale <- eta2theta(eta[, !vTF, drop = FALSE], .lscale , .escale )
    kappamat <- matrix(extra$kappa, extra$n, extra$Mdiv2,
                       byrow = TRUE)
    ralap(nsim * length(Scale), location = c(locat),
          scale = c(Scale), kappa = c(kappamat))
  }, list( .elocat = elocat, .llocat = llocat,
           .escale = escale, .lscale = lscale,
           .kappa = kappa ))),






  deriv = eval(substitute(expression({
    M1 <- 2
    Mdiv2 <- ncol(eta) / M1  # extra$Mdiv2
    ymat <- matrix(y, n, Mdiv2)
    vTF <- c(TRUE, FALSE)
  locat <- eta2theta(eta[,  vTF, drop = FALSE], .llocat , .elocat )
  Scale <- eta2theta(eta[, !vTF, drop = FALSE], .lscale , .escale )
    kappamat <- matrix(extra$kappa, n, Mdiv2, byrow = TRUE)
    zedd <- abs(ymat - locat) / Scale
    dl.dlocat <- sqrt(2) * ifelse(ymat >= locat,
                                  kappamat, 1/kappamat) *
                 sign(ymat - locat) / Scale
    dl.dscale <- sqrt(2) * ifelse(ymat >= locat,
                                  kappamat, 1/kappamat) *
                 zedd / Scale - 1 / Scale
    dlocat.deta <- dtheta.deta(locat, .llocat , earg = .elocat )
    dscale.deta <- dtheta.deta(Scale, .lscale , earg = .escale )

    ans <- c(w) * cbind(dl.dlocat * dlocat.deta,
                        dl.dscale * dscale.deta)
    ans[, interleave.VGAM(ncol(ans), M1 = M1)]
  }), list( .escale = escale, .lscale = lscale,
            .elocat = elocat, .llocat = llocat,
            .kappa = kappa ))),
  weight = eval(substitute(expression({
    wz <- matrix(NA_real_, n, M)
    d2l.dlocat2 <- 2 / Scale^2
    d2l.dscale2 <- 1 / Scale^2
    wz[,  vTF] <- d2l.dlocat2 * dlocat.deta^2
    wz[, !vTF] <- d2l.dscale2 * dscale.deta^2
    c(w) * wz
  }), list( .escale = escale, .lscale = lscale,
            .elocat = elocat, .llocat = llocat ))))
}  # End of alaplace2().











alaplace1.control <- function(maxit = 100, ...) {
    list(maxit = maxit)
}









 alaplace1 <-
  function(tau = NULL,
           llocation = "identitylink",
           ilocation = NULL,
           kappa = sqrt(tau/(1-tau)),
           Scale.arg = 1,
           ishrinkage = 0.95,
           parallel.locat = TRUE  ~ 0,  # FALSE,
           digt = 4,
           idf.mu = 3,
           zero = NULL,
           imethod = 1) {



  apply.parint.locat <- FALSE




  if (!is.Numeric(kappa, positive = TRUE))
    stop("bad input for argument 'kappa'")
  if (length(tau) &&
      max(abs(kappa - sqrt(tau/(1-tau)))) > 1.0e-6)
    stop("arguments 'kappa' and 'tau' do not match")

  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 4)
    stop("argument 'imethod' must be 1, 2 or ... 4")


  llocation <- llocation

  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")
  ilocat <- ilocation


  if (!is.Numeric(ishrinkage, length.arg = 1) ||
     ishrinkage < 0 ||
     ishrinkage > 1)
    stop("bad input for argument 'ishrinkage'")
  if (!is.Numeric(Scale.arg, positive = TRUE))
    stop("bad input for argument 'Scale.arg'")




  fittedMean <- FALSE
  if (!is.logical(fittedMean) || length(fittedMean) != 1)
    stop("bad input for argument 'fittedMean'")





  new("vglmff",
  blurb = c("One-parameter asymmetric Laplace distribution\n\n",
            "Links:      ",
            namesof("location", llocat, earg = elocat),
            "\n", "\n",
            "Mean:       location + scale * (1/kappa - kappa) / ",
                         "sqrt(2)", "\n",
            "Quantiles:  location", "\n",
            "Variance:   scale^2 * (1 + kappa^4) / (2 * kappa^2)"),




  constraints = eval(substitute(expression({

    onemat <- matrix(1, M, 1)
    constraints.orig <- constraints


    cm1.locat <- diag(M)
    cmk.locat <- onemat
    con.locat <- cm.VGAM(cmk.locat,
                         x = x, bool = .parallel.locat ,
                         constraints = constraints.orig,
                         apply.int = .apply.parint.locat ,
                         cm.default           = cm1.locat,
                         cm.intercept.default = cm1.locat)


    constraints <- con.locat

    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 1)
  }), list( .parallel.locat = parallel.locat,
            .zero = zero,
            .apply.parint.locat = apply.parint.locat ))),



  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         summary.pvalues = FALSE,
         tau   = .tau ,
         multipleResponses = FALSE,
         parameters.names = c("location"),
         kappa = .kappa)
  }, list( .kappa = kappa,
           .tau   = tau ))),
  initialize = eval(substitute(expression({
    extra$M1 <- M1 <- 1


    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = if (length( .kappa ) > 1) 1 else Inf,
              ncol.y.max = if (length( .kappa ) > 1) 1 else Inf,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    extra$ncoly <- ncoly <- ncol(y)
    if ((ncoly > 1) && (length( .kappa ) > 1 ||
        length( .Scale.arg ) > 1))
      stop("response must be a vector if 'kappa' or 'Scale.arg' ",
           "has a length greater than one")

    extra$kappa <- .kappa
    extra$tau <- extra$kappa^2 / (1 + extra$kappa^2)


    extra$M <- M <- max(length( .Scale.arg ),
                        ncoly,
                        length( .kappa ))  # Recycle
    extra$Scale <- rep_len( .Scale.arg , M)
    extra$kappa <- rep_len( .kappa     , M)
    extra$tau <- extra$kappa^2 / (1 + extra$kappa^2)
    extra$n <- n




    extra$tau.names <- tau.names <-
      paste0("(tau = ", round(extra$tau, digits = .digt), ")")
    extra$Y.names <- Y.names <- if (ncoly > 1)
                                dimnames(y)[[2]] else "y"
    if (is.null(Y.names) || any(Y.names == ""))
      extra$Y.names <- Y.names <- paste0("y", 1:ncoly)
    extra$y.names <- y.names <-
        if (ncoly > 1) paste0(Y.names, tau.names) else
                       tau.names

    extra$individual <- FALSE

    mynames1 <- param.names("location", M, skip1 = TRUE)
    predictors.names <-
        c(namesof(mynames1, .llocat , earg = .elocat , tag = FALSE))


    locat.init <- matrix(0, n, M)
    if (!length(etastart)) {

      for (jay in 1:M) {
        y.use <- if (ncoly > 1) y[, jay] else y
        if ( .imethod == 1) {
            locat.init[, jay] <-
                weighted.mean(y.use, w[, min(jay, ncol(w))])
        } else if ( .imethod == 2) {
          locat.init[, jay] <- median(y.use)
        } else if ( .imethod == 3) {
          Fit5 <- vsmooth.spline(x = x[, min(ncol(x), 2)],
                                 y = y.use, w = w, df = .idf.mu )
          locat.init[, jay] <- c(predict(Fit5,
                                         x = x[, min(ncol(x), 2)])$y)
        } else {
          use.this <- weighted.mean(y.use, w[, min(jay, ncol(w))])
          locat.init[, jay] <- (1- .ishrinkage ) * y.use +
              .ishrinkage * use.this
        }


        if (length( .ilocat )) {
          locat.init <- matrix( .ilocat  , n, M, byrow = TRUE)
        }

        if ( .llocat == "loglink") locat.init <- abs(locat.init)
        etastart <-
          cbind(theta2eta(locat.init, .llocat , earg = .elocat ))
      }
    }
    }), list( .imethod = imethod,
              .idf.mu = idf.mu,
              .ishrinkage = ishrinkage, .digt = digt,
              .elocat = elocat, .Scale.arg = Scale.arg,
              .llocat = llocat, .kappa = kappa,
              .ilocat = ilocat ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    if ( .fittedMean ) {
      kappamat <- matrix(extra$kappa, extra$n, extra$M, byrow = TRUE)
      locat <- eta2theta(eta, .llocat , earg = .elocat )
      Scale <- matrix(extra$Scale, extra$n, extra$M, byrow = TRUE)
      locat + Scale * (1/kappamat - kappamat)
    } else {
      locat <- eta2theta(eta, .llocat , earg = .elocat )
      if (length(locat) > extra$n)
        dimnames(locat) <- list(dimnames(eta)[[1]], extra$y.names)
      locat
    }
  }, list( .elocat = elocat, .llocat = llocat,
           .fittedMean = fittedMean, .Scale.arg = Scale.arg,
           .kappa = kappa ))),
  last = eval(substitute(expression({
    M1 <- extra$M1
    misc$M1 <- M1
    misc$multipleResponses <- TRUE



    misc$link <- setNames(rep_len( .llocat , M), mynames1)







    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:M) {
      misc$earg[[ii]] <- .elocat
    }


    misc$expected <- TRUE
    extra$kappa <- misc$kappa <- .kappa
    extra$tau <- misc$tau <- misc$kappa^2 / (1 + misc$kappa^2)
    misc$true.mu <- .fittedMean # @fitted is not a true mu?

    extra$percentile <- numeric(M)
    locat <- as.matrix(locat)
    for (ii in 1:M) {
      y.use <- if (ncoly > 1) y[, ii] else y
      extra$percentile[ii] <-
          100 * weighted.mean(y.use <= locat[, ii],
                              w[, min(ii, ncol(w))])
    }
    names(extra$percentile) <- y.names

    extra$Scale.arg <- .Scale.arg
    }), list( .elocat = elocat,
              .llocat = llocat,
              .Scale.arg = Scale.arg, .fittedMean = fittedMean,
              .kappa = kappa ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {

    ymat <- matrix(y, extra$n, extra$M)
    locat <- eta2theta(eta, .llocat , earg = .elocat )
    kappamat <- matrix(extra$kappa, extra$n, extra$M, byrow = TRUE)
    Scale    <- matrix(extra$Scale, extra$n, extra$M, byrow = TRUE)

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dalap(x = c(ymat), locat = c(locat),
                              sc = c(Scale), kap = c(kappamat),
                              log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .elocat = elocat,
           .llocat = llocat,
           .Scale.arg = Scale.arg, .kappa = kappa ))),
  vfamily = c("alaplace1"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    locat <- eta2theta(eta, .llocat , earg = .elocat )
    okay1 <- all(is.finite(locat))
    okay1
  }, list( .elocat = elocat,
           .llocat = llocat,
           .Scale.arg = Scale.arg, .kappa = kappa ))),


  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    extra    <- object@extra
    locat    <- eta2theta(eta, .llocat , .elocat )
    Scale    <- matrix(extra$Scale, extra$n, extra$M, byrow = TRUE)
    kappamat <- matrix(extra$kappa, extra$n, extra$M, byrow = TRUE)
    ralap(nsim * length(Scale), location = c(locat),
          scale = c(Scale), kappa = c(kappamat))
  }, list( .elocat = elocat, .llocat = llocat,
           .Scale.arg = Scale.arg, .kappa = kappa ))),



  deriv = eval(substitute(expression({
    ymat <- matrix(y, n, M)
    Scale <- matrix(extra$Scale, extra$n, extra$M, byrow = TRUE)

    locat <- eta2theta(eta, .llocat , earg = .elocat )

    kappamat <- matrix(extra$kappa, n, M, byrow = TRUE)
    zedd <- abs(ymat-locat) / Scale

    dl.dlocat <- ifelse(ymat >= locat, kappamat, 1/kappamat) *
                   sqrt(2) * sign(ymat - locat) / Scale
    dlocat.deta <- dtheta.deta(locat, .llocat , earg = .elocat )

    c(w) * cbind(dl.dlocat * dlocat.deta)
  }), list( .Scale.arg = Scale.arg, .elocat = elocat,
            .llocat = llocat, .kappa = kappa ))),

  weight = eval(substitute(expression({
    d2l.dlocat2 <- 2 / Scale^2
    wz <- cbind(d2l.dlocat2 * dlocat.deta^2)

    c(w) * wz
  }), list( .Scale.arg = Scale.arg,
            .elocat = elocat, .llocat = llocat ))))
}









alaplace3.control <- function(maxit = 100, ...) {
  list(maxit = maxit)
}




 alaplace3 <-
     function(llocation = "identitylink",
              lscale = "loglink", lkappa = "loglink",
              ilocation = NULL,
              iscale = NULL,   ikappa = 1.0,
           imethod = 1, zero = c("scale", "kappa")) {

  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")
  ilocat <- ilocation

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")

  lkappa <- as.list(substitute(lkappa))
  ekappa <- link2list(lkappa)
  lkappa <- attr(ekappa, "function.name")


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
    stop("argument 'imethod' must be 1 or 2")
  if (length(iscale) &&
      !is.Numeric(iscale, positive = TRUE))
    stop("bad input for argument 'iscale'")


  new("vglmff",
  blurb = c("Three-parameter asymmetric Laplace distribution\n\n",
            "Links:    ",
            namesof("location", llocat, earg = elocat), ", ",
            namesof("scale",    lscale, earg = escale), ", ",
            namesof("kappa",    lkappa, earg = ekappa),
            "\n", "\n",
            "Mean:     location + scale * (1/kappa - kappa) / sqrt(2)",
            "\n",
            "Variance: Scale^2 * (1 + kappa^4) / (2 * kappa^2)"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 3)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 3,
         Q1 = 1,
         multipleResponses = FALSE,
         parameters.names = c("location", "scale", "kappa"),
         summary.pvalues = FALSE,
         zero = .zero )
  }, list( .zero = zero ))),

  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1)



    predictors.names <-
      c(namesof("location", .llocat , earg = .elocat, tag = FALSE),
        namesof("scale",    .lscale , earg = .escale, tag = FALSE),
        namesof("kappa",    .lkappa , earg = .ekappa, tag = FALSE))

    if (!length(etastart)) {
      kappa.init <- if (length( .ikappa ))
                   rep_len( .ikappa , n) else
                   rep_len( 1.0     , n)
      if ( .imethod == 1) {
        locat.init <- median(y)
        scale.init <- sqrt(var(y) / 2)
      } else {
        locat.init <- y
        scale.init <- sqrt(sum(c(w)*abs(y-median(y ))) / (sum(w) *2))
      }
      locat.init <- if (length( .ilocat ))
                       rep_len( .ilocat , n) else
                       rep_len(locat.init, n)
      scale.init <- if (length( .iscale ))
                       rep_len( .iscale , n) else
                       rep_len(scale.init, n)
      etastart <-
          cbind(theta2eta(locat.init, .llocat , earg = .elocat ),
                theta2eta(scale.init, .lscale , earg = .escale ),
                theta2eta(kappa.init, .lkappa, earg = .ekappa))
    }
  }), list( .imethod = imethod,
            .elocat = elocat, .escale = escale, .ekappa = ekappa,
            .llocat = llocat, .lscale = lscale, .lkappa = lkappa,
            .ilocat = ilocat, .iscale = iscale, .ikappa = ikappa ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    locat <- eta2theta(eta[, 1], .llocat , earg = .elocat )
    Scale <- eta2theta(eta[, 2], .lscale , earg = .escale )
    kappa <- eta2theta(eta[, 3], .lkappa, earg = .ekappa)
    locat + Scale * (1/kappa - kappa) / sqrt(2)
  }, list( .elocat = elocat, .llocat = llocat,
           .escale = escale, .lscale = lscale,
           .ekappa = ekappa, .lkappa = lkappa ))),
  last = eval(substitute(expression({
    misc$link <-    c(location = .llocat ,
                      scale    = .lscale ,
                      kappa    = .lkappa )

    misc$earg <- list(location = .elocat,
                      scale    = .escale,
                      kappa    = .ekappa )

    misc$expected = TRUE
  }), list( .elocat = elocat, .llocat = llocat,
            .escale = escale, .lscale = lscale,
            .ekappa = ekappa, .lkappa = lkappa ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    locat <- eta2theta(eta[, 1], .llocat , .elocat )
    Scale <- eta2theta(eta[, 2], .lscale , .escale )
    kappa <- eta2theta(eta[, 3], .lkappa , .ekappa )  # a matrix
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dalap(x = y, locat = locat,
                              sc = Scale, kap = kappa, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .elocat = elocat, .llocat = llocat,
           .escale = escale, .lscale = lscale,
           .ekappa = ekappa, .lkappa = lkappa ))),
  vfamily = c("alaplace3"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    locat <- eta2theta(eta[, 1], .llocat , earg = .elocat )
    Scale <- eta2theta(eta[, 2], .lscale , earg = .escale )
    kappa <- eta2theta(eta[, 3], .lkappa , earg = .ekappa )
    okay1 <- all(is.finite(locat)) &&
             all(is.finite(Scale)) && all(0 < Scale) &&
             all(is.finite(kappa)) && all(0 < kappa)
    okay1
  }, list( .elocat = elocat, .llocat = llocat,
           .escale = escale, .lscale = lscale,
           .ekappa = ekappa, .lkappa = lkappa ))),
  deriv = eval(substitute(expression({
    locat <- eta2theta(eta[, 1], .llocat , earg = .elocat )
    Scale <- eta2theta(eta[, 2], .lscale , earg = .escale )
    kappa <- eta2theta(eta[, 3], .lkappa , earg = .ekappa )

    zedd <- abs(y - locat) / Scale
    dl.dlocat <- sqrt(2) * ifelse(y >= locat, kappa, 1/kappa) *
                   sign(y-locat) / Scale
    dl.dscale <-  sqrt(2) * ifelse(y >= locat, kappa, 1/kappa) *
                 zedd / Scale - 1 / Scale
    dl.dkappa <-  1 / kappa - 2 * kappa / (1+kappa^2) -
                 (sqrt(2) / Scale) *
                 ifelse(y > locat, 1, -1/kappa^2) * abs(y-locat)

    dlocat.deta <- dtheta.deta(locat, .llocat , earg = .elocat )
    dscale.deta <- dtheta.deta(Scale, .lscale , earg = .escale )
    dkappa.deta <- dtheta.deta(kappa, .lkappa, earg = .ekappa)

    c(w) * cbind(dl.dlocat * dlocat.deta,
                 dl.dscale * dscale.deta,
                 dl.dkappa * dkappa.deta)
  }), list( .escale = escale, .lscale = lscale,
            .elocat = elocat, .llocat = llocat,
            .ekappa = ekappa, .lkappa = lkappa ))),
  weight = eval(substitute(expression({
    d2l.dlocat2 <- 2 / Scale^2
    d2l.dscale2 <- 1 / Scale^2
    d2l.dkappa2 <- 1 / kappa^2 + 4 / (1+kappa^2)^2
    d2l.dkappadloc <- -sqrt(8) / ((1+kappa^2) * Scale)
    d2l.dkappadscale <- -(1-kappa^2) / ((1+kappa^2) *
                                        kappa * Scale)
    wz <- matrix(0, nrow = n, dimm(M))
    wz[,iam(1, 1, M)] <- d2l.dlocat2 * dlocat.deta^2
    wz[,iam(2, 2, M)] <- d2l.dscale2 * dscale.deta^2
    wz[,iam(3, 3, M)] <- d2l.dkappa2 * dkappa.deta^2
    wz[,iam(1, 3, M)] <- d2l.dkappadloc *
        dkappa.deta * dlocat.deta
    wz[,iam(2, 3, M)] <- d2l.dkappadscale  *
        dkappa.deta * dscale.deta
    c(w) * wz
  }), list( .escale = escale, .lscale = lscale,
            .elocat = elocat, .llocat = llocat ))))
}









dlaplace <- function(x, location = 0, scale = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  logdensity <- (-abs(x-location)/scale) - log(2*scale)
  if (log.arg) logdensity else exp(logdensity)
}



plaplace <- function(q, location = 0, scale = 1,
                     lower.tail = TRUE, log.p =FALSE) {
  zedd <- (q - location) / scale

  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  L <- max(length(q), length(location), length(scale))
  if (length(q)        != L) q        <- rep_len(q,        L)
  if (length(location) != L) location <- rep_len(location, L)
  if (length(scale)    != L) scale    <- rep_len(scale,    L)


  if (lower.tail) {
    if (log.p) {
      ans <- ifelse(q < location, log(0.5) + zedd,
                                  log1p(- 0.5 * exp(-zedd)))
    } else {
        ans <- ifelse(q < location, 0.5 * exp(zedd),
                      1 - 0.5 * exp(-zedd))
    }
  } else {
    if (log.p) {
      ans <- ifelse(q < location, log1p(- 0.5 * exp(zedd)),
                                  log(0.5) - zedd)
    } else {
        ans <- ifelse(q < location, 1 - 0.5 *
                                    exp(zedd), 0.5 * exp(-zedd))
    }
  }
  ans[scale <= 0] <- NaN
  ans
}



qlaplace <- function(p, location = 0, scale = 1,
                     lower.tail = TRUE, log.p = FALSE) {
  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")


  L <- max(length(p), length(location), length(scale))
  if (length(p)        != L) p        <- rep_len(p,        L)
  if (length(location) != L) location <- rep_len(location, L)
  if (length(scale)    != L) scale    <- rep_len(scale,    L)


  if (lower.tail) {
    if (log.p) {
      ln.p <- p
      ans <- location - sign(exp(ln.p)-0.5) * scale *
          log(2 * ifelse(exp(ln.p) < 0.5,
                         exp(ln.p), -expm1(ln.p)))
    } else {
      ans <- location - sign(p-0.5) * scale *
          log(2 * ifelse(p < 0.5, p, 1-p))
    }
  } else {
    if (log.p) {
      ln.p <- p
      ans <- location - sign(0.5 - exp(ln.p)) * scale *
          log(2 * ifelse(-expm1(ln.p) < 0.5,
                         -expm1(ln.p), exp(ln.p)))
     # ans[ln.p > 0] <- NaN
    } else {
      ans <- location - sign(0.5 - p) * scale *
             log(2 * ifelse(p > 0.5, 1 - p, p))
    }
  }

  ans[scale <= 0] <- NaN
  ans
}



rlaplace <- function(n, location = 0, scale = 1) {

  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n

  if (!is.Numeric(scale, positive = TRUE))
    stop("'scale' must be positive")

  location <- rep_len(location, use.n)
  scale    <- rep_len(scale,    use.n)
  rrrr     <- runif(use.n)



  location - sign(rrrr - 0.5) * scale *
  (log(2) + ifelse(rrrr < 0.5, log(rrrr), log1p(-rrrr)))
}





 laplace <- function(llocation = "identitylink", lscale = "loglink",
                     ilocation = NULL, iscale = NULL,
                     imethod = 1,
                     zero = "scale") {

  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")
  ilocat <- ilocation

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")



  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")


  if (length(iscale) &&
      !is.Numeric(iscale, positive = TRUE))
    stop("bad input for argument 'iscale'")


  new("vglmff",
  blurb = c("Two-parameter Laplace distribution\n\n",
            "Links:    ",
            namesof("location", llocat, earg = elocat), ", ",
            namesof("scale",    lscale, earg = escale),
            "\n", "\n",
            "Mean:     location", "\n",
            "Variance: 2*scale^2"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         multipleResponses = FALSE,
         parameters.names = c("location", "scale"),
         summary.pvalues = FALSE,
         zero = .zero )
  }, list( .zero = zero ))),

  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1)




    predictors.names <-
      c(namesof("location", .llocat , earg = .elocat, tag = FALSE),
        namesof("scale",    .lscale , earg = .escale, tag = FALSE))


    if (!length(etastart)) {
      if ( .imethod == 1) {
        locat.init <- median(y)
        scale.init <- sqrt(var(y) / 2)
      } else if ( .imethod == 2) {
        locat.init <- weighted.mean(y, w)
        scale.init <- sqrt(var(y) / 2)
      } else {
        locat.init <- median(y)
        scale.init <- sqrt(sum(c(w)*abs(y-median(y ))) / (sum(w) *2))
      }
      locat.init <- if (length( .ilocat ))
                       rep_len( .ilocat , n) else
                       rep_len(locat.init, n)
      scale.init <- if (length( .iscale ))
                       rep_len( .iscale , n) else
                       rep_len(scale.init, n)
      etastart <-
          cbind(theta2eta(locat.init, .llocat , earg = .elocat ),
                theta2eta(scale.init, .lscale , earg = .escale ))
    }
  }), list( .imethod = imethod,
            .elocat = elocat, .escale = escale,
            .llocat = llocat, .lscale = lscale,
            .ilocat = ilocat, .iscale = iscale ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta[, 1], .llocat , earg = .elocat )
  }, list( .elocat = elocat, .llocat = llocat ))),
  last = eval(substitute(expression({
    misc$link <-    c(location = .llocat , scale = .lscale )

    misc$earg <- list(location = .elocat , scale = .escale )

    misc$expected <- TRUE
    misc$RegCondOK <- FALSE # Save this for later
  }), list( .escale = escale, .lscale = lscale,
            .elocat = elocat, .llocat = llocat ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {

    locat <- eta2theta(eta[, 1], .llocat , earg = .elocat )
    Scale <- eta2theta(eta[, 2], .lscale , earg = .escale )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dlaplace(x = y, locat = locat,
                                 scale = Scale, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .escale = escale, .lscale = lscale,
           .elocat = elocat, .llocat = llocat ))),
  vfamily = c("laplace"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    Locat <- eta2theta(eta[, 1], .llocat , earg = .elocat )
    Scale <- eta2theta(eta[, 2], .lscale , earg = .escale )
    okay1 <- all(is.finite(Locat)) &&
             all(is.finite(Scale)) && all(0 < Scale)
    okay1
  }, list( .escale = escale, .lscale = lscale,
           .elocat = elocat, .llocat = llocat ))),
  deriv = eval(substitute(expression({
    Locat <- eta2theta(eta[, 1], .llocat , earg = .elocat )
    Scale <- eta2theta(eta[, 2], .lscale , earg = .escale )

    zedd <- abs(y-Locat) / Scale
    dl.dLocat <- sign(y - Locat) / Scale
    dl.dscale <-  zedd / Scale - 1 / Scale

    dLocat.deta <- dtheta.deta(Locat, .llocat , earg = .elocat )
    dscale.deta <- dtheta.deta(Scale, .lscale , earg = .escale )

    c(w) * cbind(dl.dLocat * dLocat.deta,
                 dl.dscale    * dscale.deta)
  }), list( .escale = escale, .lscale = lscale,
            .elocat = elocat, .llocat = llocat ))),
  weight = eval(substitute(expression({
    d2l.dLocat2 <- d2l.dscale2 <- 1 / Scale^2
    wz <- matrix(0, nrow = n, ncol = M)  # diagonal
    wz[,iam(1, 1, M)] <- d2l.dLocat2 * dLocat.deta^2
    wz[,iam(2, 2, M)] <- d2l.dscale2 * dscale.deta^2
    c(w) * wz
  }), list( .escale = escale, .lscale = lscale,
            .elocat = elocat, .llocat = llocat ))))
}



fff.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}




fff <-
  function(link = "loglink",
           idf1 = NULL, idf2 = NULL, nsimEIM = 100,  # ncp = 0,
           imethod = 1, zero = NULL) {
  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
    stop("argument 'imethod' must be 1 or 2")


  if (!is.Numeric(nsimEIM, length.arg = 1,
                  integer.valued = TRUE) ||
      nsimEIM <= 10)
    stop("argument 'nsimEIM' should be an integer greater than 10")

  ncp <- 0
  if (any(ncp != 0))
    warning("not sure about ncp != 0 wrt dl/dtheta")



  new("vglmff",
  blurb = c("F-distribution\n\n",
            "Links:    ",
            namesof("df1", link, earg = earg), ", ",
            namesof("df2", link, earg = earg),
            "\n", "\n",
            "Mean:     df2/(df2-2) provided df2>2 and ncp = 0", "\n",
            "Variance: ",
            "2*df2^2*(df1+df2-2)/(df1*(df2-2)^2*(df2-4)) ",
            "provided df2>4 and ncp = 0"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         multipleResponses = FALSE,
         parameters.names = c("df1", "df2"),
         zero = .zero )
  }, list( .zero = zero ))),

  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1)



    predictors.names <-
      c(namesof("df1", .link , earg = .earg , tag = FALSE),
        namesof("df2", .link , earg = .earg , tag = FALSE))


    if (!length(etastart)) {
      if ( .imethod == 1) {
        df2.init <- b <- 2*mean(y) / (mean(y)-1)
        df1.init <- 2*b^2*(b-2)/(var(y)*(b-2)^2 * (b-4) - 2*b^2)
        if (df2.init < 4) df2.init <- 5
        if (df1.init < 2) df1.init <- 3
      } else {
            df2.init <- b <- 2*median(y) / (median(y)-1)
            summy <- summary(y)
            var.est <- summy[5] - summy[2]
            df1.init <- 2*b^2*(b-2)/(var.est*(b-2)^2 * (b-4) - 2*b^2)
        }
        df1.init <- if (length( .idf1 ))
                       rep_len( .idf1 , n) else
                       rep_len(df1.init, n)
        df2.init <- if (length( .idf2 ))
                       rep_len( .idf2 , n) else
                       rep_len(1, n)
        etastart <- cbind(theta2eta(df1.init, .link , earg = .earg ),
                          theta2eta(df2.init, .link , earg = .earg ))
    }
  }), list( .imethod = imethod, .idf1 = idf1, .earg = earg,
           .idf2 = idf2, .link = link ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    df2 <- eta2theta(eta[, 2], .link , earg = .earg )
    ans <- df2 * NA
    ans[df2 > 2] <- df2[df2 > 2] / (df2[df2 > 2] - 2)
    ans
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    misc$link <-    c(df1 = .link , df2 = .link )

    misc$earg <- list(df1 = .earg , df2 = .earg )

    misc$nsimEIM <- .nsimEIM
    misc$ncp <- .ncp
  }), list( .link = link, .earg = earg,
            .ncp = ncp,
            .nsimEIM = nsimEIM ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {

    df1 <- eta2theta(eta[, 1], .link , earg = .earg )
    df2 <- eta2theta(eta[, 2], .link , earg = .earg )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * df(x = y, df1 = df1, df2 = df2,
                           ncp = .ncp , log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .link = link, .earg = earg, .ncp = ncp ))),
  vfamily = c("fff"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    df1 <- eta2theta(eta[, 1], .link , earg = .earg )
    df2 <- eta2theta(eta[, 2], .link , earg = .earg )
    okay1 <- all(is.finite(df1)) && all(0 < df1) &&
             all(is.finite(df2)) && all(0 < df2)
    okay1
  }, list( .link = link, .earg = earg, .ncp = ncp ))),
  deriv = eval(substitute(expression({
    df1 <- eta2theta(eta[, 1], .link , earg = .earg )
    df2 <- eta2theta(eta[, 2], .link , earg = .earg )
    dl.ddf1 <- 0.5*digamma(0.5*(df1+df2)) +
        0.5 + 0.5*log(df1/df2) +
              0.5*log(y) - 0.5*digamma(0.5*df1) -
              0.5*(df1+df2)*(y/df2) / (1 + df1*y/df2) -
              0.5*log1p(df1*y/df2)
    dl.ddf2 <- 0.5*digamma(0.5*(df1+df2)) - 0.5*df1/df2 -
              0.5*digamma(0.5*df2) -
              0.5*(df1+df2) * (-df1*y/df2^2) / (1 + df1*y/df2) -
              0.5*log1p(df1*y/df2)
    ddf1.deta <- dtheta.deta(df1, .link , earg = .earg )
    ddf2.deta <- dtheta.deta(df2, .link , earg = .earg )
    dthetas.detas <- cbind(ddf1.deta, ddf2.deta)
      c(w) * dthetas.detas * cbind(dl.ddf1, dl.ddf2)
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({
    run.varcov <- 0
    ind1 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)
    for (ii in 1:( .nsimEIM )) {
      ysim <- rf(n = n, df1=df1, df2=df2)
      dl.ddf1 <- 0.5*digamma(0.5*(df1+df2)) +
          0.5 + 0.5*log(df1/df2) +
                0.5*log(ysim) - 0.5*digamma(0.5*df1) -
                0.5*(df1+df2)*(ysim/df2) / (1 + df1*ysim/df2) -
                0.5*log1p(df1*ysim/df2)
      dl.ddf2 <- 0.5*digamma(0.5*(df1+df2)) -
          0.5*df1/df2 -
                0.5*digamma(0.5*df2) -
                0.5*(df1+df2) *
                (-df1*ysim/df2^2)/(1 + df1*ysim/df2) -
                0.5*log1p(df1*ysim/df2)
      rm(ysim)
      temp3 <- cbind(dl.ddf1, dl.ddf2)
      run.varcov <- ((ii-1) * run.varcov +
                     temp3[,ind1$row.index]*
                     temp3[,ind1$col.index]) / ii
    }
    wz <- if (intercept.only)
        matrix(colMeans(run.varcov),
               n, ncol(run.varcov), byrow = TRUE) else run.varcov

    wz <- c(w) * wz * dthetas.detas[, ind1$row] *
                     dthetas.detas[, ind1$col]
    wz
  }), list( .link = link, .earg = earg, .nsimEIM = nsimEIM,
            .ncp = ncp ))))
}




 hyperg <- function(N = NULL, D = NULL,
                    lprob = "logitlink",
                    iprob = NULL) {

  inputN <- is.Numeric(N, positive = TRUE)
  inputD <- is.Numeric(D, positive = TRUE)
  if (inputD && inputN)
    stop("only one of 'N' and 'D' is to be inputted")
  if (!inputD && !inputN)
    stop("one of 'N' and 'D' needs to be inputted")


  lprob <- as.list(substitute(lprob))
  earg <- link2list(lprob)
  lprob <- attr(earg, "function.name")



  new("vglmff",
  blurb = c("Hypergeometric distribution\n\n",
            "Link:     ",
            namesof("prob", lprob, earg = earg), "\n",
            "Mean:     D/N\n"),
  initialize = eval(substitute(expression({
    NCOL <- function (x)
        if (is.array(x) && length(dim(x)) > 1 ||
        is.data.frame(x)) ncol(x) else as.integer(1)
    if (NCOL(y) == 1) {
        if (is.factor(y)) y <- y != levels(y)[1]
        nn <- rep_len(1, n)
        if (!all(y >= 0 & y <= 1))
            stop("response values must be in [0, 1]")
        mustart <- (0.5 + w * y) / (1 + w)
        no.successes <- w * y
        if (any(abs(no.successes - round(no.successes)) > 0.001))
            stop("Number of successes must be integer-valued")
    } else if (NCOL(y) == 2) {
        if (any(abs(y - round(y)) > 0.001))
            stop("Count data must be integer-valued")
        nn <- y[, 1] + y[, 2]
        y <- ifelse(nn > 0, y[, 1]/nn, 0)
        w <- w * nn
        mustart <- (0.5 + nn * y) / (1 + nn)
        mustart[mustart >= 1] <- 0.95
    } else
         stop("Response not of the right form")

    predictors.names <-
      namesof("prob", .lprob , earg = .earg , tag = FALSE)
    extra$Nvector <- .N
    extra$Dvector <- .D
    extra$Nunknown <- length(extra$Nvector) == 0
    if (!length(etastart)) {
        init.prob <- if (length( .iprob))
                      rep_len( .iprob, n) else
                      mustart
            etastart <- matrix(init.prob, n, NCOL(y))

    }
  }), list( .lprob = lprob, .earg = earg, .N = N, .D = D,
            .iprob = iprob ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta, .lprob , earg = .earg )
  }, list( .lprob = lprob, .earg = earg ))),
  last = eval(substitute(expression({
    misc$link <-    c("prob" = .lprob)

    misc$earg <- list("prob" = .earg )

    misc$Dvector <- .D
    misc$Nvector <- .N
  }), list( .N = N, .D = D, .lprob = lprob, .earg = earg ))),
  linkfun = eval(substitute(function(mu, extra = NULL) {
    theta2eta(mu, .lprob, earg = .earg )
  }, list( .lprob = lprob, .earg = earg ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    N <- extra$Nvector
    Dvec <- extra$Dvector
    prob <- mu
    yvec <- w * y
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
      if (extra$Nunknown) {
        tmp12 <- Dvec * (1-prob) / prob


        (lgamma(1+tmp12) + lgamma(1+Dvec/prob-w) -
         lgamma(1+tmp12-w+yvec) - lgamma(1+Dvec/prob))
      } else {


        (lgamma(1+N*prob) + lgamma(1+N*(1-prob)) -
         lgamma(1+N*prob-yvec) -
         lgamma(1+N*(1-prob) -w + yvec))
      }
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lprob = lprob, .earg = earg ))),
  vfamily = c("hyperg"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    prob <- eta2theta(eta, .lprob , earg = .earg )
    okay1 <- all(is.finite(prob)) && all(0 < prob & prob < 1)
    okay1
  }, list( .lprob = lprob, .earg = earg ))),
  deriv = eval(substitute(expression({
    prob <- mu   # equivalently, eta2theta(eta, .lprob, .earg )
    dprob.deta <- dtheta.deta(prob, .lprob, earg = .earg )
    Dvec <- extra$Dvector
    Nvec <- extra$Nvector
    yvec <- w * y
    if (extra$Nunknown) {
      tmp72 <- -Dvec / prob^2
      tmp12 <-  Dvec * (1-prob) / prob
      dl.dprob <- tmp72 * (digamma(1 + tmp12) +
                 digamma(1 + Dvec/prob -w) -
           digamma(1 + tmp12-w+yvec) - digamma(1 + Dvec/prob))
    } else {
      dl.dprob <- Nvec * (digamma(1+Nvec*prob) -
                 digamma(1+Nvec*(1-prob)) -
                 digamma(1+Nvec*prob-yvec) +
                 digamma(1+Nvec*(1-prob)-w+yvec))
    }
    c(w) * dl.dprob * dprob.deta
  }), list( .lprob = lprob, .earg = earg ))),
  weight = eval(substitute(expression({
    if (extra$Nunknown) {
      tmp722 <- tmp72^2
      tmp13 <- 2*Dvec / prob^3
      d2l.dprob2 <- tmp722 * (trigamma(1 + tmp12) +
                   trigamma(1 + Dvec/prob - w) -
                   trigamma(1 + tmp12 - w + yvec) -
                   trigamma(1 + Dvec/prob)) +
                   tmp13 * (digamma(1 + tmp12) +
                   digamma(1 + Dvec/prob - w) -
                   digamma(1 + tmp12 - w + yvec) -
                   digamma(1 + Dvec/prob))
    } else {
      d2l.dprob2 <- Nvec^2 * (trigamma(1+Nvec*prob) +
                   trigamma(1+Nvec*(1-prob)) -
                   trigamma(1+Nvec*prob-yvec) -
                   trigamma(1+Nvec*(1-prob)-w+yvec))
    }
    d2prob.deta2 <- d2theta.deta2(prob, .lprob , earg = .earg )

    wz <- -(dprob.deta^2) * d2l.dprob2
    wz <- c(w) * wz
    wz[wz < .Machine$double.eps] <- .Machine$double.eps
    wz
    }), list( .lprob = lprob, .earg = earg ))))
}






dbenini <- function(x, y0, shape, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  N <- max(length(x), length(shape), length(y0))
  if (length(x)        != N) x        <- rep_len(x,        N)
  if (length(shape)    != N) shape    <- rep_len(shape,    N)
  if (length(y0)       != N) y0       <- rep_len(y0,       N)

  logdensity <- rep_len(log(0), N)
  xok <- (x > y0)
  tempxok <- log(x[xok]/y0[xok])
  logdensity[xok] <- log(2*shape[xok]) - shape[xok] * tempxok^2 +
                     log(tempxok) - log(x[xok])
  logdensity[is.infinite(x)] <- log(0)  # 20141209 KaiH
  if (log.arg) logdensity else exp(logdensity)
}



pbenini <-
    function(q, y0, shape, lower.tail = TRUE, log.p = FALSE) {
  if (!is.Numeric(q))
    stop("bad input for argument 'q'")
  if (!is.Numeric(shape, positive = TRUE))
    stop("bad input for argument 'shape'")
  if (!is.Numeric(y0, positive = TRUE))
    stop("bad input for argument 'y0'")
  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  N <- max(length(q), length(shape), length(y0))
  if (length(q)        != N) q      <- rep_len(q,     N)
  if (length(shape)    != N) shape  <- rep_len(shape, N)
  if (length(y0)       != N) y0     <- rep_len(y0,    N)

  ans <- y0 * 0
  ok <- q > y0


  if (lower.tail) {
    if (log.p) {
      ans[ok] <- log(-expm1(-shape[ok] * (log(q[ok]/y0[ok]))^2))
      ans[q <= y0 ] <- -Inf
    } else {
      ans[ok] <- -expm1(-shape[ok] * (log(q[ok]/y0[ok]))^2)
    }
  } else {
    if (log.p) {
      ans[ok] <- -shape[ok] * (log(q[ok]/y0[ok]))^2
      ans[q <= y0] <- 0
    } else {
      ans[ok] <- exp(-shape[ok] * (log(q[ok]/y0[ok]))^2)
      ans[q <= y0] <- 1
    }
  }

  ans
}



qbenini <- function(p, y0, shape, lower.tail = TRUE, log.p = FALSE) {


  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  if (lower.tail) {
    if (log.p) {
      ln.p <- p
      ans <- y0 * exp(sqrt(-log(-expm1(ln.p)) / shape))
    } else {
      ans <- y0 * exp(sqrt(-log1p(-p) / shape))
    }
  } else {
    if (log.p) {
      ln.p <- p
      ans <- y0 * exp(sqrt(-ln.p / shape))
    } else {
      ans <-  y0 * exp(sqrt(-log(p) / shape))
    }
  }
  ans[y0 <= 0] <- NaN
  ans
}



rbenini <- function(n, y0, shape) {
  y0 * exp(sqrt(-log(runif(n)) / shape))
}







 benini1 <-
  function(y0 = stop("argument 'y0' must be specified"),
           lshape = "loglink",
           ishape = NULL, imethod = 1, zero = NULL,
           parallel = FALSE,
    type.fitted = c("percentiles", "Qlink"),
    percentiles = 50) {

  type.fitted <- match.arg(type.fitted,
                           c("percentiles", "Qlink"))[1]


  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 2)
    stop("argument 'imethod' must be 1 or 2")


  if (!is.Numeric(y0, positive = TRUE, length.arg = 1))
   stop("bad input for argument 'y0'")





  new("vglmff",
  blurb = c("1-parameter Benini distribution\n\n",
            "Link:    ",
            namesof("shape", lshape, earg = eshape),
            "\n", "\n",
            "Median:     qbenini(p = 0.5, y0, shape)"),
  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(1, M, 1), x = x,
                           bool = .parallel ,
                           constraints, apply.int = FALSE)
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 1)
  }), list( .parallel = parallel,
            .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = TRUE,
         parameters.names = c("shape"),
         parallel = .parallel ,
         percentiles = .percentiles ,
         type.fitted = .type.fitted ,
         lshape = .lshape ,
         eshape = .eshape ,
         zero = .zero )
  }, list( .parallel = parallel,
           .zero = zero,
           .percentiles = percentiles ,
           .type.fitted = type.fitted,
           .eshape = eshape,
           .lshape = lshape))),

  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    ncoly <- ncol(y)
    M1 <- 1
    M <- M1 * ncoly
    extra$ncoly <- ncoly
    extra$type.fitted <- .type.fitted
    extra$colnames.y  <- colnames(y)
    extra$percentiles <- .percentiles
    extra$M1 <- M1


    mynames1 <- paste0("shape", if (ncoly > 1) 1:ncoly else "")
    predictors.names <-
      namesof(mynames1, .lshape , earg = .eshape , tag = FALSE)

    extra$y0 <- .y0  # Of unit length; 20181205; to make things easy.
    if (any(y <= extra$y0))
      stop("some values of the response are > argument 'y0' values")


    if (!length(etastart)) {
      probs.y <- (1:3) / 4
      qofy <- quantile(rep(y, times = w), probs = probs.y)
      if ( .imethod == 1) {
        shape.init <- mean(-log1p(-probs.y) / (log(qofy))^2)
      } else {
        shape.init <- median(-log1p(-probs.y) / (log(qofy))^2)
      }
    shape.init <- matrix(if (length( .ishape )) .ishape else shape.init,
                        n, ncoly, byrow = TRUE)
    etastart <- cbind(theta2eta(shape.init, .lshape , earg = .eshape ))
  }
  }), list( .imethod = imethod,
            .ishape = ishape,
            .lshape = lshape, .eshape = eshape,
            .percentiles = percentiles,
            .type.fitted = type.fitted,
            .y0 = y0 ))),



  linkinv = eval(substitute(function(eta, extra = NULL) {
    type.fitted <-
      if (length(extra$type.fitted)) {
        extra$type.fitted
      } else {
        warning("cannot find 'type.fitted'. Returning the 'median'.")
        extra$percentiles <- 50  # Overwrite whatever was there
        "percentiles"
      }
    type.fitted <- match.arg(type.fitted,
                             c("percentiles", "Qlink"))[1]

    if (type.fitted == "Qlink") {
      eta2theta(eta, link = "loglink")
    } else {
      shape <- eta2theta(eta, .lshape , earg = .eshape )

      pcent <- extra$percentiles
      perc.mat <- matrix(pcent, NROW(eta), length(pcent),
                         byrow = TRUE) / 100
      fv <-
        switch(type.fitted,
               "percentiles" = qbenini(perc.mat,
                   y0 = extra$y0,
                   shape = matrix(shape, nrow(perc.mat), ncol(perc.mat))))
      if (type.fitted == "percentiles")
        fv <- label.cols.y(fv, colnames.y = extra$colnames.y,
                           NOS = NCOL(eta), percentiles = pcent,
                           one.on.one = FALSE)
      fv
    }
  }, list( .lshape = lshape, .eshape = eshape ))),










  last = eval(substitute(expression({
    M1 <- extra$M1
    misc$link <- c(rep_len( .lshape , ncoly))
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for (ii in 1:ncoly) {
      misc$earg[[ii]] <- .eshape
    }

    extra$y0 <- .y0

  }), list( .lshape = lshape,
            .eshape = eshape, .y0 = y0 ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    y0 <- extra$y0
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dbenini(y, y0 = y0, sh = shape, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lshape = lshape, .eshape = eshape ))),
  vfamily = c("benini1"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    okay1 <- all(is.finite(shape)) && all(0 < shape)
    okay1
  }, list( .lshape = lshape, .eshape = eshape ))),





  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    extra <- object@extra
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    y0 <- extra$y0
    rbenini(nsim * length(shape), y0 = y0, shape = shape)
  }, list( .lshape = lshape, .eshape = eshape ))),





  deriv = eval(substitute(expression({
    shape <- eta2theta(eta, .lshape , earg = .eshape )

    y0 <- extra$y0
    dl.dshape <- 1/shape - (log(y/y0))^2

    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )

    c(w) * dl.dshape * dshape.deta
  }), list( .lshape = lshape, .eshape = eshape ))),
  weight = eval(substitute(expression({
    ned2l.dshape2 <- 1 / shape^2
    wz <- ned2l.dshape2 * dshape.deta^2
    c(w) * wz
  }), list( .lshape = lshape, .eshape = eshape ))))
}






 dpolono  <-
  function (x, meanlog = 0, sdlog = 1, bigx = 170, ...) {
  mapply(function(x, meanlog, sdlog, ...) {
    if (abs(x) > floor(x)) { # zero prob for -ve or non-integer
      0
    } else
    if (x == Inf) { # 20141215 KaiH
      0
    } else
    if (x > bigx) {
      z <- (log(x) - meanlog) / sdlog
      (1 + (z^2 + log(x) - meanlog - 1) / (2 * x * sdlog^2)) *
      exp(-0.5 * z^2) / (sqrt(2 * pi) * sdlog * x)
    } else
       integrate( function(t) exp(t * x - exp(t) -
                              0.5 * ((t - meanlog) / sdlog)^2),
                 lower = -Inf,
                 upper = Inf, ...)$value / (sqrt(2 * pi) *
                 sdlog * exp(lgamma(x + 1.0)))
  }, x, meanlog, sdlog, ...)
}




 ppolono <-
    function(q, meanlog = 0, sdlog = 1,
             isOne = 1 - sqrt( .Machine$double.eps ), ...) {


 .cumprob <- rep_len(0, length(q))
 .cumprob[q == Inf] <- 1  # special case


 q <- floor(q)
 ii <-  -1
 while (any(xActive <- ((.cumprob < isOne) & (q > ii))))
    .cumprob[xActive] <- .cumprob[xActive] +
      dpolono(ii <- (ii+1), meanlog, sdlog, ...)
 .cumprob
}









rpolono <- function(n, meanlog = 0, sdlog = 1) {
  lambda <- rlnorm(n = n, meanlog = meanlog, sdlog = sdlog)
  rpois(n = n, lambda = lambda)
}














dtriangle <- function(x, theta, lower = 0, upper = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  N <- max(length(x), length(theta), length(lower), length(upper))
  if (length(x)     != N) x     <- rep_len(x,     N)
  if (length(theta) != N) theta <- rep_len(theta, N)
  if (length(lower) != N) lower <- rep_len(lower, N)
  if (length(upper) != N) upper <- rep_len(upper, N)

  denom1 <- ((upper-lower)*(theta-lower))
  denom2 <- ((upper-lower)*(upper-theta))
  logdensity <- rep_len(log(0), N)
  xok.neg <- (lower <  x) & (x <= theta)
  xok.pos <- (theta <= x) & (x <  upper)
  logdensity[xok.neg] =
    log(2 * (x[xok.neg] - lower[xok.neg]) / denom1[xok.neg])
  logdensity[xok.pos] =
    log(2 * (upper[xok.pos] - x[xok.pos]) / denom2[xok.pos])
  logdensity[lower >= upper] <- NaN
  logdensity[lower >  theta] <- NaN
  logdensity[upper <  theta] <- NaN
  if (log.arg) logdensity else exp(logdensity)
}


rtriangle <- function(n, theta, lower = 0, upper = 1) {


  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n


  if (!is.Numeric(theta))
    stop("bad input for argument 'theta'")
  if (!is.Numeric(lower))
    stop("bad input for argument 'lower'")
  if (!is.Numeric(upper))
    stop("bad input for argument 'upper'")
  if (!all(lower < theta & theta < upper))
    stop("lower < theta < upper values are required")

  N <- use.n
  lower <- rep_len(lower, N)
  upper <- rep_len(upper, N)
  theta <- rep_len(theta, N)
  t1 <- sqrt(runif(n))
  t2 <- sqrt(runif(n))
  ifelse(runif(n) < (theta - lower) / (upper - lower),
         lower + (theta - lower) * t1,
         upper - (upper - theta) * t2)
}



qtriangle <- function(p, theta, lower = 0, upper = 1,
                      lower.tail = TRUE, log.p = FALSE) {

  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  N <- max(length(p), length(theta), length(lower), length(upper))
  if (length(p)     != N) p     <- rep_len(p,     N)
  if (length(theta) != N) theta <- rep_len(theta, N)
  if (length(lower) != N) lower <- rep_len(lower, N)
  if (length(upper) != N) upper <- rep_len(upper, N)

  ans <- NA_real_ * p
  if (lower.tail) {
    if (log.p) {
      Neg <- (exp(ln.p) <= (theta - lower) / (upper - lower))
      temp1 <- exp(ln.p) * (upper - lower) * (theta - lower)
      Pos <- (exp(ln.p) >= (theta - lower) / (upper - lower))
      pstar <- (exp(ln.p) - (theta - lower) / (upper - lower)) /
               ((upper - theta) / (upper - lower))
    } else {
      Neg <- (p <= (theta - lower) / (upper - lower))
      temp1 <- p * (upper - lower) * (theta - lower)
      Pos <- (p >= (theta - lower) / (upper - lower))
      pstar <- (p - (theta - lower) / (upper - lower)) /
               ((upper - theta) / (upper - lower))
    }
  } else {
    if (log.p) {
      ln.p <- p
      Neg <- (exp(ln.p) >= (upper- theta) / (upper - lower))
      temp1 <- -expm1(ln.p) * (upper - lower) * (theta - lower)
      Pos <- (exp(ln.p) <= (upper- theta) / (upper - lower))
      pstar <- (-expm1(ln.p) - (theta - lower) / (upper - lower)) /
               ((upper - theta) / (upper - lower))
    } else {
      Neg <- (p >= (upper- theta) / (upper - lower))
      temp1 <- (1 - p) * (upper - lower) * (theta - lower)
      Pos <- (p <= (upper- theta) / (upper - lower))
      pstar <- ((upper- theta) / (upper - lower) - p) /
               ((upper - theta) / (upper - lower))
    }
  }
  ans[ Neg] <- lower[ Neg] + sqrt(temp1[ Neg])
  if (any(Pos)) {
    qstar <- cbind(1 - sqrt(1-pstar), 1 + sqrt(1-pstar))
    qstar <- qstar[Pos,, drop = FALSE]
    qstar <- ifelse(qstar[, 1] >= 0 & qstar[, 1] <= 1,
                    qstar[, 1],
                    qstar[, 2])
    ans[Pos] <- theta[Pos] + qstar * (upper - theta)[Pos]
  }

  ans[theta < lower | theta > upper] <- NaN
  ans
}



ptriangle <- function(q, theta, lower = 0, upper = 1,
                      lower.tail = TRUE, log.p = FALSE) {

  N <- max(length(q), length(theta), length(lower), length(upper))
  if (length(q)     != N) q     <- rep_len(q,     N)
  if (length(theta) != N) theta <- rep_len(theta, N)
  if (length(lower) != N) lower <- rep_len(lower, N)
  if (length(upper) != N) upper <- rep_len(upper, N)

  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  ans <- q * 0
  qstar <- (q - lower)^2 / ((upper - lower) * (theta - lower))
  Neg <- (lower <= q & q <= theta)


  ans[Neg] <- if (lower.tail) {
    if (log.p) {
      (log(qstar))[Neg]
    } else {
      qstar[Neg]
    }
  } else {
    if (log.p) {
      (log1p(-qstar))[Neg]
    } else {
      1 - qstar[Neg]
    }
  }

  Pos <- (theta <= q & q <= upper)
  qstar <- (q - theta) / (upper-theta)

  if (lower.tail) {
    if (log.p) {
      ans[Pos] <- log(((theta-lower)/(upper-lower))[Pos] +
                  (qstar * (2-qstar) *
                   (upper-theta) / (upper - lower))[Pos])
      ans[q <= lower] <- -Inf
      ans[q >= upper] <- 0
    } else {
      ans[Pos] <- ((theta-lower)/(upper-lower))[Pos] +
                  (qstar * (2-qstar) *
                   (upper-theta) / (upper - lower))[Pos]
      ans[q <= lower] <- 0
      ans[q >= upper] <- 1
    }
  } else {
    if (log.p) {
      ans[Pos] <- log(((upper - theta)/(upper-lower))[Pos] +
                  (qstar * (2-qstar) *
                   (upper-theta) / (upper - lower))[Pos])
      ans[q <= lower] <- 0
      ans[q >= upper] <- -Inf
    } else {
      ans[Pos] <- ((upper - theta)/(upper-lower))[Pos] +
                  (qstar * (2-qstar) *
                   (upper-theta) / (upper - lower))[Pos]
      ans[q <= lower] <- 1
      ans[q >= upper] <- 0
    }
  }

  ans[theta < lower | theta > upper] <- NaN
  ans
}







triangle.control <- function(stepsize = 0.33, maxit = 100, ...) {
  list(stepsize = stepsize, maxit = maxit)
}


 triangle <-
  function(lower = 0, upper = 1,
           link = extlogitlink(min = 0, max = 1),
           itheta = NULL) {






  if (!is.Numeric(lower))
    stop("bad input for argument 'lower'")
  if (!is.Numeric(upper))
    stop("bad input for argument 'upper'")
  if (!all(lower < upper))
    stop("lower < upper values are required")

  if (length(itheta) && !is.Numeric(itheta))
    stop("bad input for 'itheta'")




  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  if (length(earg$min) && any(earg$min != lower))
    stop("argument 'lower' does not match the 'link'")
  if (length(earg$max) && any(earg$max != upper))
    stop("argument 'upper' does not match the 'link'")



  new("vglmff",
  blurb = c("Triangle distribution\n\n",
            "Link:    ",
            namesof("theta", link, earg = earg)),
  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         parameters.names = c("theta"),
         link = .link )
  }, list( .link = link ))),

  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1)




    extra$lower <- rep_len( .lower , n)
    extra$upper <- rep_len( .upper , n)

    if (any(y <= extra$lower | y >= extra$upper))
      stop("some y values in [lower,upper] detected")

    predictors.names <-
      namesof("theta", .link , earg = .earg , tag = FALSE)


    if (!length(etastart)) {
      Theta.init <- if (length( .itheta )) .itheta else {
        weighted.mean(y, w)
      }
      Theta.init <- rep_len(Theta.init, n)
      etastart <- theta2eta(Theta.init, .link , earg = .earg )
    }
  }), list( .link = link, .earg = earg, .itheta=itheta,
            .upper = upper, .lower = lower ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    Theta <- eta2theta(eta, .link , earg = .earg )
    lower <- extra$lower
    upper <- extra$upper

    mu1 <- (lower + upper + Theta) / 3

    mu1
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    misc$link <-    c(theta = .link )

    misc$earg <- list(theta = .earg )

    misc$expected <- TRUE
  }), list( .link = link, .earg = earg ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    Theta <- eta2theta(eta, .link , earg = .earg )
    lower <- extra$lower
    upper <- extra$upper
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dtriangle(y, theta = Theta, lower = lower,
                                  upper = upper, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .link = link, .earg = earg ))),
  vfamily = c("triangle"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    Theta <- eta2theta(eta, .link , earg = .earg )
    okay1 <- all(is.finite(Theta)) &&
             all(extra$lower < Theta & Theta < extra$upper)
    okay1
  }, list( .link = link, .earg = earg ))),



  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1))
      warning("ignoring prior weights")
    eta <- predict(object)
    extra <- object@extra
    Theta <- eta2theta(eta, .link , earg = .earg )
    lower <- extra$lower
    upper <- extra$upper
    rtriangle(nsim * length(Theta),
              theta = Theta, lower = lower, upper = upper)
  }, list( .link = link, .earg = earg ))),




  deriv = eval(substitute(expression({
    Theta       <- eta2theta(eta,     .link , earg = .earg )
    dTheta.deta <- dtheta.deta(Theta, .link , earg = .earg )

    pos <- y > Theta
    neg <- y < Theta
    lower <- extra$lower
    upper <- extra$upper

    dl.dTheta <-  0 * y
    dl.dTheta[neg] <-  -1 / (Theta[neg]-lower[neg])
    dl.dTheta[pos] <-   1 / (upper[pos]-Theta[pos])

    c(w) * dl.dTheta * dTheta.deta
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({
    var.dl.dTheta <-  1 / ((Theta - lower) * (upper - Theta))
    wz <- var.dl.dTheta * dTheta.deta^2
    c(w) * wz
  }), list( .link = link, .earg = earg ))))
}







adjust0.loglaplace1 <- function(ymat, y, w, rep0) {
  rangey0 <- range(y[y > 0])
  ymat[ymat <= 0] <- min(rangey0[1] / 2, rep0)
  ymat
}  # adjust0.loglaplace1


loglaplace1.control <- function(maxit = 300, ...) {
  list(maxit = maxit)
}


 loglaplace1 <- function(tau = NULL,
                     llocation = "loglink",
                     ilocation = NULL,
                     kappa = sqrt(tau/(1-tau)),
                     Scale.arg = 1,
                     ishrinkage = 0.95,
                     parallel.locat = FALSE, digt = 4,
                     idf.mu = 3,
                     rep0 = 0.5,  # 0.0001,
                     minquantile = 0, maxquantile = Inf,
                     imethod = 1, zero = NULL) {

  if (length(minquantile) != 1)
    stop("bad input for argument 'minquantile'")
  if (length(maxquantile) != 1)
    stop("bad input for argument 'maxquantile'")


  if (!is.Numeric(rep0, positive = TRUE, length.arg = 1) ||
      rep0 > 1)
    stop("bad input for argument 'rep0'")
  if (!is.Numeric(kappa, positive = TRUE))
    stop("bad input for argument 'kappa'")

  if (length(tau) && max(abs(kappa - sqrt(tau/(1-tau)))) > 1.0e-6)
      stop("arguments 'kappa' and 'tau' do not match")


  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")
  ilocat <- ilocation


  llocat.identity <- as.list(substitute("identitylink"))
  elocat.identity <- link2list(llocat.identity)
  llocat.identity <- attr(elocat.identity, "function.name")


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 4)
    stop("argument 'imethod' must be 1, 2 or ... 4")


  if (!is.Numeric(ishrinkage, length.arg = 1) ||
     ishrinkage < 0 ||
     ishrinkage > 1)
    stop("bad input for argument 'ishrinkage'")

  if (!is.Numeric(Scale.arg, positive = TRUE))
    stop("bad input for argument 'Scale.arg'")
  if (!is.logical(parallel.locat) ||
      length(parallel.locat) != 1)
    stop("bad input for argument 'parallel.locat'")

  fittedMean <- FALSE
  if (!is.logical(fittedMean) || length(fittedMean) != 1)
    stop("bad input for argument 'fittedMean'")


  mystring0 <- namesof("location", llocat, earg = elocat)
  mychars <- substring(mystring0, first = 1:nchar(mystring0),
                      last = 1:nchar(mystring0))
  mychars[nchar(mystring0)] <- ", inverse = TRUE)"
  mystring1 <- paste(mychars, collapse = "")




  new("vglmff",
  blurb = c("One-parameter ",
            if (llocat == "loglink") "log-Laplace" else
              c(llocat, "-Laplace"),
            " distribution\n\n",
            "Links:      ", mystring0, "\n", "\n",
          "Quantiles:  ", mystring1),
  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(1, M, 1), x = x,
                           bool = .parallel.locat ,
                           constraints = constraints,
                           apply.int = FALSE)
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 1)
  }), list( .parallel.locat = parallel.locat,
            .Scale.arg = Scale.arg, .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         parameters.names = c("location"),
         llocation = .llocat )
  }, list( .llocat = llocat,
           .zero   = zero ))),

  initialize = eval(substitute(expression({
    extra$M <- M <- max(length( .Scale.arg ), length( .kappa ))
    extra$Scale <- rep_len( .Scale.arg , M)
    extra$kappa <- rep_len( .kappa     , M)
    extra$tau <- extra$kappa^2 / (1 + extra$kappa^2)


    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y




        extra$n <- n
        extra$y.names <- y.names <-
          paste0("tau = ", round(extra$tau, digits = .digt))
        extra$individual <- FALSE


        predictors.names <-
          namesof(paste0("quantile(", y.names, ")"),
                  .llocat , earg = .elocat , tag = FALSE)


        if (FALSE) {
        if (min(y) < 0)
          stop("negative response values detected")
        if ((prop.0. <- weighted.mean(1*(y == 0), w)) >=
            min(extra$tau))
        stop("sample proportion of 0s == ",
             round(prop.0., digits = 4),
     " > minimum 'tau' value. Choose larger values for 'tau'.")
        if ( .rep0 == 0.5 &&
            (ave.tau <- (weighted.mean(1*(y <= 0), w) +
             weighted.mean(1*(y <= 1), w))/2) >= min(extra$tau))
     warning("the minimum 'tau' value should be greater than ",
             round(ave.tau, digits = 4))
        }

        if (!length(etastart)) {
            if ( .imethod == 1) {
                locat.init <- quantile(rep(y, w),
                                       probs= extra$tau) + 1/16
            } else if ( .imethod == 2) {
              locat.init <- weighted.mean(y, w)
            } else if ( .imethod == 3) {
              locat.init <- median(y)
            } else if ( .imethod == 4) {
              Fit5 <- vsmooth.spline(x = x[, min(ncol(x), 2)],
                                     y = y, w = w,
                                     df = .idf.mu )
                locat.init <- c(predict(Fit5,
                                        x = x[, min(ncol(x), 2)])$y)
            } else {
              use.this <- weighted.mean(y, w)
              locat.init <- (1- .ishrinkage )*y +
                  .ishrinkage * use.this
            }
            locat.init <- if (length( .ilocat ))
                             rep_len( .ilocat , M) else
                             rep_len(locat.init, M)
            locat.init <- matrix(locat.init, n, M, byrow = TRUE)
            if ( .llocat == "loglink")
                locat.init <- abs(locat.init)
            etastart <-
                cbind(theta2eta(locat.init, .llocat , earg = .elocat ))
        }
    }), list( .imethod = imethod,
              .idf.mu = idf.mu, .rep0 = rep0,
              .ishrinkage = ishrinkage, .digt = digt,
              .elocat = elocat, .Scale.arg = Scale.arg,
              .llocat = llocat, .kappa = kappa,
              .ilocat = ilocat ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    locat.y = eta2theta(eta, .llocat , earg = .elocat )
    if ( .fittedMean ) {
      stop("Yet to do: handle 'fittedMean = TRUE'")
      kappamat <- matrix(extra$kappa, extra$n, extra$M, byrow = TRUE)
      Scale <- matrix(extra$Scale, extra$n, extra$M, byrow = TRUE)
      locat.y + Scale * (1/kappamat - kappamat)
    } else {
      if (length(locat.y) > extra$n)
        dimnames(locat.y) <- list(dimnames(eta)[[1]], extra$y.names)
      locat.y
    }
        locat.y[locat.y < .minquantile] = .minquantile
        locat.y[locat.y > .maxquantile] = .maxquantile
        locat.y
  }, list( .elocat = elocat, .llocat = llocat,
           .minquantile = minquantile, .maxquantile = maxquantile,
           .fittedMean = fittedMean, .Scale.arg = Scale.arg,
           .kappa = kappa ))),
  last = eval(substitute(expression({
    misc$link <-    c(location = .llocat)

    misc$earg <- list(location = .elocat )

    misc$expected <- TRUE

    extra$kappa <- misc$kappa <- .kappa
    extra$tau <- misc$tau <- misc$kappa^2 / (1 + misc$kappa^2)
    extra$Scale.arg <- .Scale.arg

    misc$true.mu <- .fittedMean # @fitted is not a true mu?
    misc$rep0 <- .rep0
    misc$minquantile <- .minquantile
    misc$maxquantile <- .maxquantile

    extra$percentile <- numeric(length(misc$kappa))
    locat.y <- as.matrix(locat.y)
    for (ii in seq_along(misc$kappa))
        extra$percentile[ii] <- 100 *
            weighted.mean(y <= locat.y[, ii], w)
  }), list( .elocat = elocat, .llocat = llocat,
            .Scale.arg = Scale.arg, .fittedMean = fittedMean,
            .minquantile = minquantile, .maxquantile = maxquantile,
            .rep0 = rep0, .kappa = kappa ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {

    kappamat <- matrix(extra$kappa, extra$n, extra$M, byrow = TRUE)
    Scale.w <- matrix(extra$Scale, extra$n, extra$M, byrow = TRUE)
    ymat <- matrix(y, extra$n, extra$M)

    if ( .llocat == "loglink")
        ymat <- adjust0.loglaplace1(ymat = ymat, y = y,
                                    w = w, rep0 = .rep0)

   w.mat <- theta2eta(ymat, .llocat , .elocat )  # e.g., logofflink()



    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dalap(x = c(w.mat), locat = c(eta),
                              sc = c(Scale.w), kappa = c(kappamat),
                              log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .elocat = elocat, .llocat = llocat,
           .rep0 = rep0,
           .Scale.arg = Scale.arg, .kappa = kappa ))),
  vfamily = c("loglaplace1"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    locat.w <- eta
    locat.y <- eta2theta(locat.w, .llocat , earg = .elocat )
    okay1 <- all(is.finite(locat.y))
    okay1
  }, list( .elocat = elocat, .llocat = llocat,
           .rep0 = rep0,
           .Scale.arg = Scale.arg, .kappa = kappa ))),
  deriv = eval(substitute(expression({
    ymat <- matrix(y, n, M)
    Scale.w <- matrix(extra$Scale, extra$n, extra$M, byrow = TRUE)
    locat.w <- eta
    locat.y <- eta2theta(locat.w, .llocat , earg = .elocat )
    kappamat <- matrix(extra$kappa, n, M, byrow = TRUE)

    ymat <- adjust0.loglaplace1(ymat = ymat, y = y,
                                w = w, rep0= .rep0)
    w.mat <- theta2eta(ymat, .llocat , .elocat )  # e.g., logitlink()
    zedd <- abs(w.mat-locat.w) / Scale.w
    dl.dlocat <- ifelse(w.mat >= locat.w, kappamat, 1/kappamat) *
                   sqrt(2) * sign(w.mat-locat.w) / Scale.w


    dlocat.deta <- dtheta.deta(locat.w,
                              .llocat.identity ,
                              earg = .elocat.identity )
    c(w) * cbind(dl.dlocat * dlocat.deta)
  }), list( .Scale.arg = Scale.arg, .rep0 = rep0,
            .llocat = llocat, .elocat = elocat,
            .elocat.identity = elocat.identity,
            .llocat.identity = llocat.identity,

            .kappa = kappa ))),
  weight = eval(substitute(expression({
    ned2l.dlocat2 <- 2 / Scale.w^2
    wz <- cbind(ned2l.dlocat2 * dlocat.deta^2)
    c(w) * wz
  }), list( .Scale.arg = Scale.arg,
            .elocat = elocat, .llocat = llocat,
            .elocat.identity = elocat.identity,
            .llocat.identity = llocat.identity  ))))
}





loglaplace2.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}


 loglaplace2 <-
   function(tau = NULL,
            llocation = "loglink", lscale = "loglink",
            ilocation = NULL, iscale = NULL,
            kappa = sqrt(tau/(1-tau)),
            ishrinkage = 0.95,
            parallel.locat = FALSE, digt = 4,
            eq.scale = TRUE,
            idf.mu = 3,
            rep0 = 0.5, nsimEIM = NULL,
            imethod = 1, zero = "(1 + M/2):M") {
 warning("it is best to use loglaplace1()")

  if (length(nsimEIM) &&
     (!is.Numeric(nsimEIM, length.arg = 1,
                  integer.valued = TRUE) ||
      nsimEIM <= 10))
    stop("argument 'nsimEIM' should be an integer greater than 10")
  if (!is.Numeric(rep0, positive = TRUE, length.arg = 1) ||
      rep0 > 1)
    stop("bad input for argument 'rep0'")
  if (!is.Numeric(kappa, positive = TRUE))
    stop("bad input for argument 'kappa'")
  if (length(tau) && max(abs(kappa - sqrt(tau/(1-tau)))) > 1.0e-6)
    stop("arguments 'kappa' and 'tau' do not match")


  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")
  ilocat <- ilocation

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")




  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 4)
    stop("argument 'imethod' must be 1, 2 or ... 4")
  if (length(iscale) && !is.Numeric(iscale, positive = TRUE))
    stop("bad input for argument 'iscale'")


  if (!is.Numeric(ishrinkage, length.arg = 1) ||
     ishrinkage < 0 ||
     ishrinkage > 1)
    stop("bad input for argument 'ishrinkage'")
  if (!is.logical(eq.scale) || length(eq.scale) != 1)
    stop("bad input for argument 'eq.scale'")
  if (!is.logical(parallel.locat) ||
      length(parallel.locat) != 1)
    stop("bad input for argument 'parallel.locat'")
  fittedMean <- FALSE
  if (!is.logical(fittedMean) || length(fittedMean) != 1)
    stop("bad input for argument 'fittedMean'")

  if (llocat != "loglink")
    stop("argument 'llocat' must be \"loglink\"")


  new("vglmff",
  blurb = c("Two-parameter log-Laplace distribution\n\n",
            "Links:      ",
            namesof("location", llocat, earg = elocat), ", ",
            namesof("scale", lscale, earg = escale),
            "\n", "\n",
            "Mean:       zz location + scale * ",
                         "(1/kappa - kappa) / sqrt(2)", "\n",
            "Quantiles:  location", "\n",
            "Variance:   zz scale^2 * (1 + kappa^4) / (2 * kappa^2)"),
  constraints = eval(substitute(expression({
  .ZERO <- .zero
  if (is.character( .ZERO ))
    .ZERO <- eval(parse(text = .ZERO ))
  .PARALLEL <- .parallel.locat
      parelHmat <- if (is.logical( .PARALLEL ) && .PARALLEL )
                   matrix(1, M/2, 1) else diag(M/2)
      scaleHmat <- if (is.logical( .eq.scale ) && .eq.scale )
                   matrix(1, M/2, 1) else diag(M/2)
      mycmatrix <- cbind(rbind(  parelHmat, 0*parelHmat),
                         rbind(0*scaleHmat,   scaleHmat))
      constraints <- cm.VGAM(mycmatrix, x = x,
                             bool = .PARALLEL ,
                             constraints = constraints,
                             apply.int = FALSE)
    constraints <- cm.zero.VGAM(constraints, x = x, .ZERO , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)

      if ( .PARALLEL && names(constraints)[1] == "(Intercept)") {
          parelHmat <- diag(M/2)
          mycmatrix <- cbind(rbind(  parelHmat, 0*parelHmat),
                            rbind(0*scaleHmat,   scaleHmat))
          constraints[["(Intercept)"]] <- mycmatrix
      }
      if (is.logical( .eq.scale) && .eq.scale &&
       names(constraints)[1] == "(Intercept)") {
        temp3 <- constraints[["(Intercept)"]]
          temp3 <- cbind(temp3[,1:(M/2)],
                         rbind(0*scaleHmat, scaleHmat))
        constraints[["(Intercept)"]] = temp3
      }
  }), list( .eq.scale = eq.scale,
            .parallel.locat = parallel.locat,
            .zero = zero ))),
  initialize = eval(substitute(expression({
    extra$kappa <- .kappa
    extra$tau <- extra$kappa^2 / (1 + extra$kappa^2)



    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y




    extra$M <- M <- 2 * length(extra$kappa)
    extra$n <- n
    extra$y.names <- y.names <-
      paste0("tau = ", round(extra$tau, digits = .digt))
    extra$individual = FALSE

    predictors.names <-
        c(namesof(paste0("quantile(", y.names, ")"),
                  .llocat , earg = .elocat, tag = FALSE),
          namesof(if (M == 2) "scale" else
                  paste0("scale", 1:(M/2)),
                  .lscale ,    earg = .escale,    tag = FALSE))
        if (weighted.mean(1 * (y < 0.001), w) >= min(extra$tau))
          stop("sample proportion of 0s > minimum 'tau' value. ",
               "Choose larger values for 'tau'.")

        if (!length(etastart)) {
          if ( .imethod == 1) {
            locat.init.y <- weighted.mean(y, w)
            scale.init <- sqrt(var(y) / 2)
          } else if ( .imethod == 2) {
            locat.init.y <- median(y)
              scale.init <- sqrt(sum(c(w)*abs(y-median(y))) / (
                  sum(w) *2))
          } else if ( .imethod == 3) {
            Fit5 <- vsmooth.spline(x = x[, min(ncol(x), 2)],
                                   y = y, w = w,
                                   df = .idf.mu )
              locat.init.y <- c(predict(Fit5,
                                        x = x[, min(ncol(x), 2)])$y)
              scale.init <- sqrt(sum(c(w)*abs(y-median(y))) /
                                 (sum(w) *2))
          } else {
            use.this <- weighted.mean(y, w)
            locat.init.y <- (1- .ishrinkage )*y +
                .ishrinkage * use.this
            scale.init <- sqrt(sum(c(w)*
                           abs(y-median(y ))) / (sum(w) *2))
          }
          locat.init.y <- if (length( .ilocat ))
                           rep_len( .ilocat , n) else
                           rep_len(locat.init.y, n)
          locat.init.y <- matrix(locat.init.y, n, M/2)
          scale.init <- if (length( .iscale ))
                           rep_len( .iscale , n) else
                           rep_len(scale.init, n)
          scale.init <- matrix(scale.init, n, M/2)
          etastart <-
            cbind(theta2eta(locat.init.y, .llocat , .elocat ),
                  theta2eta(scale.init,   .lscale , .escale ))
        }
    }), list( .imethod = imethod,
              .idf.mu = idf.mu, .kappa = kappa,
              .ishrinkage = ishrinkage, .digt = digt,
              .llocat = llocat, .lscale = lscale,
              .elocat = elocat, .escale = escale,
              .ilocat = ilocat, .iscale = iscale ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    locat.y <- eta2theta(eta[, 1:(extra$M/2), drop = FALSE],
                         .llocat , earg = .elocat )
    if ( .fittedMean ) {
      kappamat <- matrix(extra$kappa, extra$n, extra$M/2,
                        byrow = TRUE)
      Scale.y <- eta2theta(eta[,(1+extra$M/2):extra$M],
                          .lscale , earg = .escale )
      locat.y + Scale.y * (1/kappamat - kappamat)
    } else {
      dimnames(locat.y) = list(dimnames(eta)[[1]], extra$y.names)
      locat.y
    }
  }, list( .llocat = llocat, .lscale = lscale,
           .elocat = elocat, .escale = escale,
           .fittedMean = fittedMean,
           .kappa = kappa ))),
  last = eval(substitute(expression({
    misc$link <-    c(location = .llocat , scale = .lscale )

    misc$earg <- list(location = .elocat , scale = .escale )

    misc$expected <- TRUE
    extra$kappa <- misc$kappa <- .kappa
    extra$tau <- misc$tau <- misc$kappa^2 / (1 + misc$kappa^2)
    misc$true.mu <- .fittedMean # @fitted is not a true mu?
    misc$nsimEIM <- .nsimEIM
    misc$rep0 <- .rep0
        extra$percentile <- numeric(length(misc$kappa))
        locat <- as.matrix(locat.y)
        for (ii in seq_along(misc$kappa))
          extra$percentile[ii] <- 100 *
                                 weighted.mean(y <= locat.y[, ii], w)
  }), list( .elocat = elocat, .llocat = llocat,
            .escale = escale, .lscale = lscale,
            .fittedMean = fittedMean,
            .nsimEIM = nsimEIM, .rep0 = rep0,
            .kappa = kappa ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {

    kappamat <- matrix(extra$kappa, extra$n, extra$M/2, byrow = TRUE)
    Scale.w <- eta2theta(eta[, (1+extra$M/2):extra$M],
                         .lscale , earg = .escale )
    ymat <- matrix(y, extra$n, extra$M/2)
    ymat[ymat <= 0] <- min(min(y[y > 0]), .rep0 )  # Adjust for 0s
    ell.mat <- matrix(c(dloglaplace(x = c(ymat),
                          locat.ald = c(eta[, 1:(extra$M/2)]),
                          scale.ald = c(Scale.w),
                          kappa = c(kappamat), log = TRUE)),
                      extra$n, extra$M/2)
      if (residuals) {
        stop("loglikelihood residuals not implemented yet")
      } else {
      ll.elts <- c(w) * ell.mat
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .elocat = elocat, .llocat = llocat,
           .escale = escale, .lscale = lscale,
           .rep0 = rep0, .kappa = kappa ))),


  vfamily = c("loglaplace2"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    Scale.w <- eta2theta(eta[, (1+extra$M/2):extra$M],
                        .lscale , earg = .escale )
    locat.w <- eta[, 1:(extra$M/2), drop = FALSE]
    locat.y <- eta2theta(locat.w, .llocat , earg = .elocat )
    okay1 <- all(is.finite(locat.y)) &&
             all(is.finite(Scale.w)) && all(0 < Scale.w)
    okay1
  }, list( .elocat = elocat, .llocat = llocat,
           .escale = escale, .lscale = lscale,
           .rep0 = rep0, .kappa = kappa ))),
  deriv = eval(substitute(expression({
    ymat <- matrix(y, n, M/2)
    Scale.w <- eta2theta(eta[, (1+extra$M/2):extra$M],
                        .lscale , earg = .escale )
    locat.w <- eta[, 1:(extra$M/2), drop = FALSE]
    locat.y <- eta2theta(locat.w, .llocat , earg = .elocat )
    kappamat <- matrix(extra$kappa, n, M/2, byrow = TRUE)
    w.mat <- ymat
    w.mat[w.mat <= 0] <- min(min(w.mat[w.mat > 0]), .rep0)
    w.mat <- theta2eta(w.mat, .llocat , earg = .elocat )
    zedd <- abs(w.mat-locat.w) / Scale.w
    dl.dlocat <- sqrt(2) *
                   ifelse(w.mat >= locat.w, kappamat, 1/kappamat) *
                   sign(w.mat-locat.w) / Scale.w
    dl.dscale <-  sqrt(2) *
                 ifelse(w.mat >= locat.w, kappamat, 1/kappamat) *
                 zedd / Scale.w - 1 / Scale.w
    dlocat.deta <- dtheta.deta(locat.w, .llocat , earg = .elocat )
    dscale.deta <- dtheta.deta(Scale.w, .lscale , earg = .escale )
    c(w) * cbind(dl.dlocat * dlocat.deta,
                 dl.dscale * dscale.deta)
  }), list( .escale = escale, .lscale = lscale,
            .elocat = elocat, .llocat = llocat,
            .rep0 = rep0, .kappa = kappa ))),
  weight = eval(substitute(expression({
    run.varcov <- 0
    ind1 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)
    dthetas.detas <- cbind(dlocat.deta, dscale.deta)
    if (length( .nsimEIM )) {
        for (ii in 1:( .nsimEIM )) {
            wsim <- matrix(rloglap(n*M/2, loc = c(locat.w),
                                  sca = c(Scale.w),
                                  kappa = c(kappamat)), n, M/2)
            zedd <- abs(wsim-locat.w) / Scale.w
            dl.dlocat <- sqrt(2) *
                ifelse(wsim >= locat.w, kappamat, 1/kappamat) *
                sign(wsim-locat.w) / Scale.w
            dl.dscale <-  sqrt(2) *
                ifelse(wsim >= locat.w, kappamat, 1/kappamat) *
                zedd / Scale.w - 1 / Scale.w

            rm(wsim)
            temp3 <- cbind(dl.dlocat, dl.dscale)  # n x M matrix
            run.varcov <- ((ii-1) * run.varcov +
               temp3[,ind1$row.index]*temp3[,ind1$col.index]) / ii
        }
        wz <- if (intercept.only)
            matrix(colMeans(run.varcov),
                   n, ncol(run.varcov), byrow = TRUE) else
            run.varcov

        wz <- wz * dthetas.detas[,ind1$row] *
            dthetas.detas[,ind1$col]
        wz <- c(w) * matrix(wz, n, dimm(M))
        wz
    } else {
        d2l.dlocat2 <- 2 / (Scale.w * locat.w)^2
        d2l.dscale2 <- 1 / Scale.w^2
        wz <- cbind(d2l.dlocat2 * dlocat.deta^2,
                   d2l.dscale2 * dscale.deta^2)
        c(w) * wz
    }
  }), list( .elocat = elocat, .escale = escale,
            .llocat = llocat, .lscale = lscale,
            .nsimEIM = nsimEIM) )))
}






logitlaplace1.control <- function(maxit = 300, ...) {
    list(maxit = maxit)
}


adjust01.logitlaplace1 <- function(ymat, y, w, rep01) {
    rangey01 <- range(y[(y > 0) & (y < 1)])
    ymat[ymat <= 0] <- min(rangey01[1] / 2, rep01 / w[y <= 0])
    ymat[ymat >= 1] <- max((1 + rangey01[2]) / 2,
                           1 - rep01 / w[y >= 1])
    ymat
}





 logitlaplace1 <-
  function(tau = NULL,
           llocation = "logitlink",
           ilocation = NULL,
           kappa = sqrt(tau/(1-tau)),
           Scale.arg = 1,
           ishrinkage = 0.95, parallel.locat = FALSE, digt = 4,
           idf.mu = 3,
           rep01 = 0.5,
           imethod = 1, zero = NULL) {

  if (!is.Numeric(rep01, positive = TRUE, length.arg = 1) ||
      rep01 > 0.5)
    stop("bad input for argument 'rep01'")
  if (!is.Numeric(kappa, positive = TRUE))
    stop("bad input for argument 'kappa'")

  if (length(tau) && max(abs(kappa - sqrt(tau/(1-tau)))) > 1.0e-6)
    stop("arguments 'kappa' and 'tau' do not match")


  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")
  ilocat <- ilocation


  llocat.identity <- as.list(substitute("identitylink"))
  elocat.identity <- link2list(llocat.identity)
  llocat.identity <- attr(elocat.identity, "function.name")




  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 4)
    stop("argument 'imethod' must be 1, 2 or ... 4")

  if (!is.Numeric(ishrinkage, length.arg = 1) ||
     ishrinkage < 0 ||
     ishrinkage > 1)
    stop("bad input for argument 'ishrinkage'")

  if (!is.Numeric(Scale.arg, positive = TRUE))
    stop("bad input for argument 'Scale.arg'")
  if (!is.logical(parallel.locat) ||
      length(parallel.locat) != 1)
    stop("bad input for argument 'parallel.locat'")
  fittedMean <- FALSE
  if (!is.logical(fittedMean) ||
      length(fittedMean) != 1)
    stop("bad input for argument 'fittedMean'")


  mystring0 <- namesof("location", llocat, earg = elocat)
  mychars <- substring(mystring0, first = 1:nchar(mystring0),
                      last = 1:nchar(mystring0))
  mychars[nchar(mystring0)] = ", inverse = TRUE)"
  mystring1 <- paste(mychars, collapse = "")




  new("vglmff",
  blurb = c("One-parameter ", llocat, "-Laplace distribution\n\n",
            "Links:      ", mystring0, "\n", "\n",
          "Quantiles:  ", mystring1),
  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(1, M, 1), x = x,
                           bool = .parallel.locat ,
                           constraints = constraints,
                           apply.int = FALSE)
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 1)
  }), list( .parallel.locat = parallel.locat,
            .Scale.arg = Scale.arg, .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         multipleResponses = FALSE,
         parameters.names = c("location"),
         llocation = .llocat ,
         zero = .zero )
  }, list( .zero = zero,
           .llocat = llocat ))),

  initialize = eval(substitute(expression({
    extra$M <- M <- max(length( .Scale.arg ), length( .kappa ))
    extra$Scale <- rep_len( .Scale.arg , M)
    extra$kappa <- rep_len( .kappa     , M)
    extra$tau <- extra$kappa^2 / (1 + extra$kappa^2)



    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y





    extra$n <- n
    extra$y.names <- y.names <-
      paste0("tau = ", round(extra$tau, digits = .digt))
    extra$individual <- FALSE

    predictors.names <-
        namesof(paste0("quantile(", y.names, ")"),
                .llocat , earg = .elocat, tag = FALSE)

      if (all(y == 0 | y == 1))
        stop("response cannot be all 0s or 1s")
      if (min(y) < 0)
        stop("negative response values detected")
      if (max(y) > 1)
        stop("response values greater than 1 detected")
    if ((prop.0. <- weighted.mean(1*(y == 0), w)) >= min(extra$tau))
    stop("sample proportion of 0s == ", round(prop.0., digits = 4),
         " > minimum 'tau' value. Choose larger values for 'tau'.")
    if ((prop.1. <- weighted.mean(1*(y == 1), w)) >= max(extra$tau))
    stop("sample proportion of 1s == ", round(prop.1., digits = 4),
         " < maximum 'tau' value. Choose smaller values for 'tau'.")
      if (!length(etastart)) {
        if ( .imethod == 1) {
          locat.init <- quantile(rep(y, w), probs= extra$tau)
        } else if ( .imethod == 2) {
          locat.init <- weighted.mean(y, w)
          locat.init <- median(rep(y, w))
        } else if ( .imethod == 3) {
          use.this <- weighted.mean(y, w)
          locat.init <- (1- .ishrinkage ) * y +
              use.this * .ishrinkage
        } else {
          stop("this option not implemented")
        }


      locat.init <- if (length( .ilocat ))
                       rep_len( .ilocat  , M) else
                       rep_len(locat.init, M)
      locat.init <- matrix(locat.init, n, M, byrow = TRUE)
      locat.init <- abs(locat.init)
      etastart <-
          cbind(theta2eta(locat.init, .llocat , earg = .elocat ))
    }
  }), list( .imethod = imethod,
            .idf.mu = idf.mu,
            .ishrinkage = ishrinkage, .digt = digt,
            .elocat = elocat, .Scale.arg = Scale.arg,
            .llocat = llocat, .kappa = kappa,
            .ilocat = ilocat ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    locat.y <- eta2theta(eta, .llocat , earg = .elocat )
    if ( .fittedMean ) {
      stop("Yet to do: handle 'fittedMean = TRUE'")
      kappamat <- matrix(extra$kappa, extra$n, extra$M, byrow = TRUE)
      Scale <- matrix(extra$Scale, extra$n, extra$M, byrow = TRUE)
      locat.y + Scale * (1/kappamat - kappamat)
    } else {
      if (length(locat.y) > extra$n)
        dimnames(locat.y) <- list(dimnames(eta)[[1]], extra$y.names)
      locat.y
      }
  }, list( .elocat = elocat, .llocat = llocat,
           .fittedMean = fittedMean, .Scale.arg = Scale.arg,
           .kappa = kappa ))),
  last = eval(substitute(expression({
    misc$link <-    c(location = .llocat )
    misc$earg <- list(location = .elocat )

    misc$expected <- TRUE

    extra$kappa <- misc$kappa <- .kappa
    extra$tau <- misc$tau <- misc$kappa^2 / (1 + misc$kappa^2)
    extra$Scale.arg <- .Scale.arg

    misc$true.mu <- .fittedMean # @fitted is not a true mu?
    misc$rep01 <- .rep01

    extra$percentile <- numeric(length(misc$kappa))
    locat.y <- eta2theta(eta, .llocat , earg = .elocat )
    locat.y <- as.matrix(locat.y)
    for (ii in seq_along(misc$kappa))
      extra$percentile[ii] <- 100 *
                             weighted.mean(y <= locat.y[, ii], w)

  }), list( .elocat = elocat, .llocat = llocat,
            .Scale.arg = Scale.arg, .fittedMean = fittedMean,
            .rep01 = rep01,
            .kappa = kappa ))),


  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {

    kappamat <- matrix(extra$kappa, extra$n, extra$M, byrow = TRUE)
    Scale.w  <- matrix(extra$Scale, extra$n, extra$M, byrow = TRUE)
    ymat     <- matrix(y,           extra$n, extra$M)
    ymat <- adjust01.logitlaplace1(ymat = ymat, y = y, w = w,
                                   rep01 = .rep01)
    w.mat <- theta2eta(ymat, .llocat , .elocat )  # e.g., logitlink()
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        c(w) * dalap(x = c(w.mat), location = c(eta),
                     scale = c(Scale.w), kappa = c(kappamat),
                     log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .elocat = elocat, .llocat = llocat,
           .rep01 = rep01,
           .Scale.arg = Scale.arg, .kappa = kappa ))),


  vfamily = c("logitlaplace1"),
  validparams = eval(substitute(function(eta, y, extra = NULL) {
    locat.w <- eta
    okay1 <- all(is.finite(locat.w))
    okay1
  }, list( .Scale.arg = Scale.arg, .rep01 = rep01,
           .elocat = elocat,
           .llocat = llocat,

           .elocat.identity = elocat.identity,
           .llocat.identity = llocat.identity,

           .kappa = kappa ))),
  deriv = eval(substitute(expression({
    ymat <- matrix(y, n, M)
    Scale.w <- matrix(extra$Scale, extra$n, extra$M, byrow = TRUE)
    locat.w <- eta
    kappamat <- matrix(extra$kappa, n, M, byrow = TRUE)
    ymat <- adjust01.logitlaplace1(ymat = ymat, y = y, w = w,
                                   rep01 = .rep01 )
    w.mat <- theta2eta(ymat, .llocat , .elocat )  # e.g., logitlink()
    zedd <- abs(w.mat - locat.w) / Scale.w
    dl.dlocat <- ifelse(w.mat >= locat.w, kappamat, 1/kappamat) *
                 sqrt(2) * sign(w.mat-locat.w) / Scale.w


    dlocat.deta <- dtheta.deta(locat.w,
                               "identitylink",
                               earg = .elocat.identity )


    c(w) * cbind(dl.dlocat * dlocat.deta)
  }), list( .Scale.arg = Scale.arg, .rep01 = rep01,
            .elocat = elocat,
            .llocat = llocat,

            .elocat.identity = elocat.identity,
            .llocat.identity = llocat.identity,

            .kappa = kappa ))),
  weight = eval(substitute(expression({
    d2l.dlocat2 <- 2 / Scale.w^2
    wz <- cbind(d2l.dlocat2 * dlocat.deta^2)
    c(w) * wz
  }), list( .Scale.arg = Scale.arg,
            .elocat = elocat, .llocat = llocat ))))
}



















