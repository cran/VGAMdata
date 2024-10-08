\name{oizeta}
\alias{oizeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ One-inflated Zeta Distribution Family Function }
\description{
  Fits a 1-inflated zeta distribution.


}
\usage{
oizeta(lpstr1 = "logitlink", lshape = "loglink",
  type.fitted = c("mean", "shape", "pobs1", "pstr1", "onempstr1"),
  ishape = NULL, gpstr1 = ppoints(8), gshape = exp((-3:3) / 4),
  zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lpstr1, lshape}{
    For \code{lpstr1}: the same idea as \code{\link[VGAM]{zipoisson}}
    except it applies to a structural 1.


  }
  \item{gpstr1, gshape, ishape}{
  For initial values.
  See \code{\link[VGAM]{CommonVGAMffArguments}} for information.


  }
  \item{type.fitted, zero}{
  See \code{\link[VGAM]{CommonVGAMffArguments}} for information.

  }


}
\details{
  The 1-inflated zeta distribution is a mixture
  distribution of the
  zeta
  distribution with some probability of obtaining a (structural) 1.
  Thus there are two sources for obtaining the value 1.
  This distribution is written here
  in a way that retains a similar notation to the
  zero-inflated Poisson, i.e., the
  probability \eqn{P[Y=1]} involves another parameter
  \eqn{\phi}{phi}.  See \code{\link[VGAM]{zipoisson}}.


  This family function can handle multiple responses.


}
\section{Warning }{
    Under- or over-flow may occur if the data is ill-conditioned.
    Lots of data is needed to estimate the parameters accurately.
    Usually, probably the shape parameter is best modelled as
    intercept-only.

}

\value{
  An object of class \code{"vglmff"}
  (see \code{\link[VGAM]{vglmff-class}}).
  The object is used by modelling functions
  such as \code{\link[VGAM]{vglm}},
  \code{\link[VGAM]{rrvglm}} and \code{\link[VGAM]{vgam}}.


}
%\references{
%}
\author{ Thomas W. Yee }
%\note{
%}
\seealso{
  \code{\link[VGAMdata]{Oizeta}},
  \code{\link[VGAM]{zetaff}},
  \code{\link[VGAMdata]{oazeta}},
  \code{\link[VGAMdata]{otzeta}},
  \code{\link[VGAM]{diffzeta}},
  \code{\link[VGAM]{zeta}},
  \code{\link[VGAMdata]{Oizipf}}.



}
\examples{
\dontrun{ odata <- data.frame(x2 = runif(nn <- 1000))  # Artificial data
odata <- transform(odata, pstr1 = logitlink(-1 + x2, inv = TRUE),
                          shape = exp(-0.5))
odata <- transform(odata, y1 = roizeta(nn, shape, pstr1 = pstr1))
with(odata, table(y1))
fit1 <- vglm(y1 ~ x2, oizeta(zero = "shape"), odata, trace = TRUE)
coef(fit1, matrix = TRUE)
}}
\keyword{models}
\keyword{regression}

