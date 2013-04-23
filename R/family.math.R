# These functions are Copyright (C) 1998-2013 T. W. Yee   All rights reserved.

# Numerical functions not present in older version of S-PLUS are here.
# 10/9/99; 30/10/00


# 20100706
# I now distribute this file with VGAM.
# Some important old functions are in ./old.code.keep/family.math.R20100706.







# =========================================================================
# 20100705
# Ref: Journal of Statisical Computation and Simulation.
# 79(11--12), p.1323 in 1317--1329, 2009.
# Picked up a copy from 2010 NZSA conference at Palmerston North.

# This seems to work.
# Ref: Corless paper.

lambertW <- function(x, tolerance = 1.0e-10, maxit = 50) {
  if (any(Im(x) != 0.0))
    stop("argument 'x' must be real, not complex!")

  ans <- x
  ans[!is.na(x) & x <  -exp(-1)] <- NA
  ans[!is.na(x) & x >= -exp(-1)] <- log1p(x[!is.na(x) & x >= -exp(-1)])
  ans[!is.na(x) & x >= 0       ] <-  sqrt(x[!is.na(x) & x >= 0       ]) / 2
# ans[x >= 2       ] <- (1.0 + log(x[x >= 2       ])) / 2
# ans[x >= 40      ] <- -0.5 + 2 * log(x[x >= 40      ])

  cutpt <- 3.0
  if (any(myTF <- !is.na(x) & x > cutpt)) {
    L1 <- log(x[!is.na(x) & x > cutpt])  # log(as.complex(x))
    L2 <- log(L1) # log(as.complex(L1))
#print("L1")
#print( L1 )
#print("L2")
#print( L2 )
#print("cbind(x, (L2/L1)^5)")
#print( cbind(x[x > cutpt], (L2/L1)^5) )
# wzinit0 is the 'original':
#   wzinit0= L1 - L2 +
#            L2/L1 +
#            L2*( -2 + L2)/(2*L1^2) +
#            L2*(  6 + L2*(-9 + L2*   2)) / (6 * L1^3) +
#            L2*(-12 + L2*(36 + L2*(-22 + L2*3))) / (12 * L1^4)
    wzinit <- L1 - L2 +
          (L2 +
          (L2*( -2 + L2)/(2) +
          (L2*(  6 + L2*(-9 + L2*   2)) / (6) +
           L2*(-12 + L2*(36 + L2*(-22 + L2*3))) / (12*L1)) / L1) / L1) / L1

#print("max(abs(wzinit0 - wzinit))")
#print( max(abs(wzinit0 - wzinit)) )
#print("Re(wzinit)")
#print( Re(wzinit) )
    ans[myTF] <- wzinit
  }

  for (ii in 1:maxit) {
    exp1 <- exp(ans)
    exp2 <- ans * exp1
    delta <- (exp2 - x) / (exp2 + exp1 -
                ((ans + 2) * (exp2 - x) / (2 * (ans + 1.0))))
#   delta <- (ans * exp(ans) - x) / ((ans + 1)*exp(ans) -
#               ((ans + 2) * (ans * exp(ans) - x) / (2*ans + 2)))
    ans <- ans - delta
#print(c(ii, max(abs(delta), na.rm = TRUE))) # Print the correction
# The na.rm = TRUE handles Inf and NA, and is.na() handles lambertW(Inf):
# The is.na() should be tested first.
    if (all(is.na(delta) ||
        max(abs(delta), na.rm = TRUE) < tolerance)) break
    if (ii == maxit)
      warning("did not converge")
  }
  ans[x == Inf] <- Inf
  ans
}




# =========================================================================


 pgamma.deriv <- function(q, shape, tmax = 100) {
#digami <- function(x, p, tmax = 100) {
# 20130214;
# Seems to work.
# Nb. shape == p.

  nnn <- max(length(q), length(shape))
  if (length(q) != nnn)
    q <- rep(q, length = nnn)
  if (length(shape) != nnn)
    shape <- rep(shape, length = nnn)

  if (!is.Numeric(q, positive = TRUE))
    stop("bad input for argument 'q'")
  if (!is.Numeric(shape, positive = TRUE))
    stop("bad input for argument 'shape'")

  if (!is.Numeric(tmax, allowable.length = 1, positive = TRUE))
    stop("bad input for argument 'tmax'")
  if (tmax < 10)
    warning("probably argument 'tmax' is too small")


  gplog  <- lgamma(shape)
  gp1log <- gplog + log(shape)
  psip   <- digamma(shape)
  psip1  <- psip + 1 / shape
  psidp  <- trigamma(shape)
  psidp1 <- psidp - 1 / shape^2

  fred <-
    dotC(name = "VGAM_C_vdigami",
         d = as.double(matrix(0, 6, nnn)),
         x = as.double(q), p = as.double(shape),
         as.double(gplog), as.double(gp1log), as.double(psip),
         as.double(psip1), as.double(psidp), as.double(psidp1),
         ifault = integer(nnn),
         tmax = as.double(tmax),
         as.integer(nnn))
#     NAOK = TRUE
  answer <- matrix(fred$d, nnn, 6, byrow = TRUE)
  dimnames(answer) <- list(names(q),
                           c("q", "q^2", "shape", "shape^2",
                             "q.shape", "pgamma(q, shape)"))

  if (any(fred$ifault != 0)) {
    indices <- which(fred$ifault != 0)
    warning("convergence problems with elements ",
             indices)
  }

  answer
}



# =========================================================================



