\name{wffc.points}
\alias{wffc.P1}
\alias{wffc.P1star}
\alias{wffc.P2}
\alias{wffc.P2star}
\alias{wffc.P3}
\alias{wffc.P3star}
\alias{wffc.P4}
\alias{wffc.P4star}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Point System for the 2008 World Fly Fishing Championships }
\description{
  Point system for the 2008 World Fly Fishing Championships:
  current and some proposals.

}
\usage{
wffc.P1(length, c1 = 100, min.eligible = 0.18, ppm = 2000)
wffc.P2(length, c1 = 100, min.eligible = 0.18, ppm = 2000,
        c.quad = 12700)
wffc.P3(length, c1 = 100, min.eligible = 0.18, ppm = 2000)
wffc.P1star(length, c1 = 100, min.eligible = 0.18, ppm = 2000)
wffc.P2star(length, c1 = 100, min.eligible = 0.18, ppm = 2000,
            c.quad = 12700)
wffc.P3star(length, c1 = 100, min.eligible = 0.18, ppm = 2000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{length}{ Length of the fish, in meters. Numeric vector. }
  \item{c1}{ Points added to each eligible fish. }
  \item{min.eligible}{ The 2008 WFFC regulations stipulated that the
  smallest eligible fish was 0.180 m, which is 180 mm. }
  \item{ppm}{ Points per meter of length of the fish. }

  \item{c.quad}{
  Constants for the quadratic terms.
% Choosing the (approximate) values 110 and 12700 for
% \code{wffc.P2()} and \code{wffc.P2star()} respectively
  The defaults
  mean that a fish twice the minimum legal size will award
  about 50 percent
  more points compared to \code{wffc.P1()} and \code{wffc.P1star()}. 
  See below for examples.



  }

}
\details{
  The official website contains a document with the official rules and
  regulations of the competition.
  The function \code{wffc.P1()} implements the current WFFC point system,
  and is `discrete' in that fish lengths are rounded up to the nearest
  centimeter (provided it is greater or equal to \code{min.eligible} m).
  \code{wffc.P1star()} is a continuous version of it in that it is
  piecewise linear with two pieces and it is discontinuous at
  \code{min.eligible}.


  The function \code{wffc.P2()} is a new proposal which
  rewards catching bigger fish.
  It is based on a quadratic polynomial.
  \code{wffc.P2star()} is a continuous version of it.


  The function \code{wffc.P3()} is another new proposal which
  rewards catching bigger fish.
  Named a \emph{cumulative linear proposal},
  its slope is \code{ppm} between \code{min.eligible} and
  \code{2 * min.eligible},
  its slope is \code{2 * ppm} between \code{2 * min.eligible} and
  \code{3 * min.eligible},
  its slope is \code{3 * ppm} between \code{3 * min.eligible} and
  \code{4 * min.eligible}, etc.
  One adds the usual \code{c1} to each eligible fish.
  \code{wffc.P3star()} is a continuous version of \code{wffc.P3()}.



  The function \code{wffc.P4()} is another new proposal which
  rewards catching bigger fish.
  Named a \emph{cumulative linear proposal},
  its slope is \code{ppm} between \code{min.eligible} and
  \code{2 * min.eligible},
  its slope is \code{2 * ppm} between \code{2 * min.eligible} and
  \code{1.5 * min.eligible},
  its slope is \code{3 * ppm} between \code{1.5 * min.eligible} and
  \code{2 * min.eligible}, etc.
  One adds the usual \code{c1} to each eligible fish.
  \code{wffc.P4star()} is a continuous version of \code{wffc.P4()}.



}
\value{
  A vector with the number of points.


}

\references{
%  \url{http://www.2008worldflyfishingchamps.com}
%  was the official 2008 website.
%  \url{http://www.http://san2010.pl}
%  was the official 2010 website.



  Yee, T. W. (2014)
  Scoring rules, and the role of chance:
  analysis of the 2008 World Fly Fishing Championships.
  \emph{Journal of Quantitative Analysis in Sports}.
  \bold{10}, 397--409. 



}

\author{ T. W. Yee. }
\note{
  \code{wffc.P2} and \code{wffc.P2star} may change in the future,
  as well as possibly
  \code{wffc.P3} and \code{wffc.P3star}
  and
  \code{wffc.P4} and \code{wffc.P4star}.


}
\seealso{
  \code{\link[VGAMdata]{wffc}}.


}
\examples{
\dontrun{ fishlength <- seq(0.0, 0.72, by = 0.001)
plot(fishlength, wffc.P2star(fishlength), type = "l", col = "blue",
     las = 1, lty = "dashed", lwd = 2, las = 1, cex.main = 0.8,
     xlab = "Fish length (m)", ylab = "Competition points",
     main = "Current (red) and proposed (blue and green) WFFC point system")
lines(fishlength, wffc.P1star(fishlength), type = "l", col = "red", lwd = 2)
lines(fishlength, wffc.P3star(fishlength), type = "l", col = "darkgreen",
      lwd = 2, lty = "dashed")
lines(fishlength, wffc.P4star(fishlength), type = "l", col = "orange",
      lwd = 2, lty = "dashed")
abline(v = (1:4) * 0.18, lty = "dotted")
abline(h = (1:13) * wffc.P1star(0.18), lty = "dotted") }

# Successive slopes:
(wffc.P1star((2:8)*0.18) - wffc.P1star((1:7)*0.18)) / (0.18 * 2000)
(wffc.P3star((2:8)*0.18) - wffc.P3star((1:7)*0.18)) / (0.18 * 2000)
(wffc.P4star((2:8)*0.18) - wffc.P4star((1:7)*0.18)) / (0.18 * 2000)

# More successive slopes:
MM2 <- 0.18 / 2
ind1 <- 2:12
(wffc.P4star((ind1)*MM2) - wffc.P4star((ind1-1)*MM2)) / (MM2 * 2000)

# About 50 percent more points:
wffc.P2    (2 * 0.18) / wffc.P1    (2 * 0.18)
wffc.P2star(2 * 0.18) / wffc.P1star(2 * 0.18)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }



%       c.quad = 127
%           c.quad = 10000



