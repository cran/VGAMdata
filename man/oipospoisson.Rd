\name{oipospoisson}
\alias{oipospoisson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ One-inflated Positive Poisson Distribution Family Function }
\description{
  Fits a 1-inflated positive Poisson distribution.
}
\usage{
oipospoisson(lpstr1 = "logitlink", llambda = "loglink",
   type.fitted = c("mean", "lambda", "pobs1", "pstr1", "onempstr1"),
   ilambda = NULL, gpstr1 = (1:19)/20, gprobs.y = (1:19)/20,
   imethod = 1, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lpstr1, llambda}{
    For \code{lpstr1}: the same idea as \code{\link[VGAM]{zipoisson}}
    except it applies to a structural 1.


  }
  \item{ilambda, gpstr1, gprobs.y, imethod}{
  For initial values.
  See \code{\link[VGAM]{CommonVGAMffArguments}} for information.


  }
  \item{type.fitted, zero}{
  See \code{\link[VGAM]{CommonVGAMffArguments}} for information.

  }


}
\details{
  The 1-inflated positive Poisson distribution is a mixture
  distribution of the
  positive (0-truncated) Poisson
  distribution with some probability of obtaining a (structural) 1.
  Thus there are two sources for obtaining the value 1.
  It is similar to a zero-inflated Poisson model, except
  the Poisson is replaced by a positive Poisson and the 0 is replaced
  by 1.
  This distribution is written here
  in a way that retains a similar notation to the
  zero-inflated Poisson, i.e., the
  probability \eqn{P[Y=1]} involves another parameter \eqn{\phi}{phi}.
  See \code{\link[VGAM]{zipoisson}}.


  This family function can handle multiple responses.


}
\section{Warning }{
    Under- or over-flow may occur if the data is ill-conditioned.

}

\value{
  An object of class \code{"vglmff"} (see \code{\link[VGAM]{vglmff-class}}).
  The object is used by modelling functions such as \code{\link[VGAM]{vglm}},
  \code{\link[VGAM]{rrvglm}} and \code{\link[VGAM]{vgam}}.


}
%\references{
%}
\author{ Thomas W. Yee }
%\note{
%}
\seealso{
  \code{\link[VGAMdata]{Oipospois}},
  \code{\link[VGAM]{pospoisson}},
  \code{\link[VGAMdata]{oapospoisson}},
  \code{\link[VGAMdata]{otpospoisson}},
  \code{\link[VGAM]{zipoisson}},
  \code{\link[VGAM]{poissonff}},
  \code{\link[VGAM]{simulate.vlm}}.


}
\examples{
\dontrun{set.seed(1)
pdata <- data.frame(x2 = runif(nn <- 1000))  # Artificial data
pdata <- transform(pdata, pstr1 = 0.5, lambda = exp(3 - x2))
pdata <- transform(pdata, y1 = roipospois(nn, lambda, pstr1 = pstr1))
with(pdata, table(y1))
fit1 <- vglm(y1 ~ x2, oipospoisson, data = pdata, trace = TRUE)
coef(fit1, matrix = TRUE)
}}
\keyword{models}
\keyword{regression}

