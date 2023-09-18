\name{mbflood}
\alias{mbflood}
\docType{data}
\title{ Madawask Basin Flood Data }

% 20230612
\description{
  Flood peak and flood volume in the
  Madawask (River) basin, Quebec, Canada

}
\usage{
data(mbflood)
}
\format{
  A data frame with the following variables.

  \describe{

    \item{year}{
      Numeric.

%     From 1919 to 1995.

    }
    \item{D}{
      Numeric.
      Duration, in days.

      

    }
    \item{V}{
      Numeric.
      Volume (days per cubic metre per second).

    }
    \item{Q}{
      Numeric.
      Flood peak (per cubic metre per second).

    }

  }
}
\details{
  From Shen (2001),
the basin is located in the province of Quebec, Canada, and
it covers an area of 2690 square km.  The data comprises 77
years daily streamflow data from 1919 to 1995 from the HYDAT CD
(Environment Canada, 1998). In the basin, winter lasts for about
4 months and precipitation is mainly in the form of snowfall
during this period. Spring represents the high flow season due
to the contribution of spring snowmelt to river runoff and the
spring flood is the annual maximum flood both in flood peak
and volume. A typical flood hydrograph with its characteristic
values are:
flood peak \eqn{Q},
flood volume \eqn{V},
and flood duration \eqn{D}.
These flood characteristics are determined by first identifying
the flood duration: the dates of start and end of flood runoff
using the approach proposed by Yue (2000).
Generally, time boundaries of a flood are marked by a rise
in stage and discharge from base flow (start of flood runoff)
and a return to base flow (end of flood runoff).





% cdff; cd gammaff.mm


}
\source{
  Table 1 of
  Yue, S. (2001).
  A bivariate gamma distribution for use in
  multivariate flood frequency analysis.
  \emph{Hydrological Processes}, \bold{15}, 1033--45.


}

\references{

  Yue, S. (2000).
  The bivariate lognormal distribution to model a
  multivariate flood episode.
  \emph{Hydrological Processes}, \bold{14}, 2575--88.


}
%\seealso{
% \code{\link[VGAMdata]{profs.nz}}.



%}
\examples{
summary(mbflood)
}
\keyword{datasets}


%
%
