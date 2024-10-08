\name{otpospoisson}
\alias{otpospoisson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ One-truncated Poisson Distribution  }
\description{
  Estimating the (single) parameter of the 1-truncated
  positive Poisson distribution.

}
\usage{
 otpospoisson(llambda = "loglink",
    type.fitted = c("mean", "lambda", "prob0", "prob1"),
    ilambda = NULL, imethod = 1, zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{llambda, type.fitted, ilambda}{
  Same as \code{\link[VGAM]{pospoisson}}.


  }
  \item{imethod, zero}{
  Same as \code{\link[VGAM]{pospoisson}}.
  See \code{\link[VGAM]{CommonVGAMffArguments}} for information.


  }
}
\details{
  The 1-truncated positive Poisson distribution has support on 2, 3,
  \ldots.
  It is a Poisson distribution but with
  the probability of a one  or zero being 0. The other
  probabilities are scaled to add to unity.
  Some more details can be found at \code{\link[VGAM]{pospoisson}}.
  Multiple responses are permitted.


}
\value{
  An object of class \code{"vglmff"}
  (see \code{\link[VGAM]{vglmff-class}}).
  The object is used by modelling functions
  such as \code{\link[VGAM]{vglm}},
  and \code{\link[VGAM]{vgam}}.


}
%\references{
%}
\author{ T. W. Yee }
%\note{
%}

\seealso{
  \code{\link[VGAMdata]{Otpospois}},
  \code{\link[VGAMdata]{oipospoisson}},
  \code{\link[VGAM]{simulate.vlm}}.


}
\examples{
\dontrun{
odata <- data.frame(y1 = rotpospois(1000, lambda = loglink(1, inv = TRUE)))
ofit <- vglm(y1 ~ 1, otpospoisson, data = odata, trace = TRUE, crit = "c")
coef(ofit, matrix = TRUE)
Coef(ofit)
with(odata,
    hist(y1, prob = TRUE, breaks = seq(0.5, max(y1) + 0.5, by = 1),
         border = "blue"))
x <- seq(1, with(odata, max(y1)), by = 1)
with(odata, lines(x, dotpospois(x, Coef(ofit)[1]), col = "orange",
                  type = "h", lwd = 2)) }
}
\keyword{models}
\keyword{regression}
