\name{otlog}
\alias{otlog}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ One-truncated Logarithmic Distribution  }
\description{
  Estimating the (single) parameter of the 1-truncated
  logarithmic distribution.

}
\usage{
otlog(lshape = "logitlink", gshape = ppoints(8), zero = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lshape, gshape, zero}{
  Same as
  \code{\link[VGAM]{logff}}.
  See \code{\link[VGAM]{CommonVGAMffArguments}} for information.


  }
}
\details{
  The 1-truncated logarithmic distribution is a logarithmic
  distribution but with the probability of a one being zero. The
  other probabilities are scaled to add to unity.
  Some more details can be found at \code{\link[VGAM]{logff}}.
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
  \code{\link[VGAM]{Gaitdlog}},
  \code{\link[VGAMdata]{Otlog}},
  \code{\link[VGAM]{logff}},
  \code{\link[VGAMdata]{oalog}},
  \code{\link[VGAMdata]{oilog}},
  \code{\link[VGAM]{simulate.vlm}}.


}
\examples{
\dontrun{
odata <- data.frame(y1 = rotlog(1000, shape = logitlink(1/3, inv = TRUE)))
ofit <- vglm(y1 ~ 1, otlog, data = odata, trace = TRUE, crit = "c")
coef(ofit, matrix = TRUE)
Coef(ofit)
with(odata,
    hist(y1, shape = TRUE, breaks = seq(0.5, max(y1) + 0.5, by = 1),
         border = "blue"))
x <- seq(1, with(odata, max(y1)), by = 1)
with(odata, lines(x, dotlog(x, Coef(ofit)[1]), col = "orange",
                  type = "h", lwd = 2)) }
}
\keyword{models}
\keyword{regression}
