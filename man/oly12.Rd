\name{oly12}
\alias{oly12}
\docType{data}
\title{
  2012 Summer Olympics: Individuals Data
}
\description{
  Individual data for the Summer
  2012 Olympic Games.

}
\usage{data(oly12)}
\format{
  A data frame with 10384 observations on the following 14 variables.
  \describe{
    \item{\code{Name}}{The individual competitor's name. }
    \item{\code{Country}}{Country. }
    \item{\code{Age}}{A numeric vector, age in years. }
    \item{\code{Height}}{A numeric vector, height in m. }
    \item{\code{Weight}}{A numeric vector, weight in kg. }
    \item{\code{Sex}}{A factor with levels \code{F} and \code{M}. }
    \item{\code{DOB}}{A Date, date of birth. }
    \item{\code{PlaceOB}}{Place of birth. }

    \item{\code{Gold}}{Numeric vector,
          number of such medals won. }
    \item{\code{Silver}}{ Similar to \code{Gold}. }
    \item{\code{Bronze}}{ Similar to \code{Gold}. }

    \item{\code{Total}}{A numeric vector, total number of medals. }
    \item{\code{Sport}}{A factor with levels
    \code{Archery},
    \code{Athletics},
    \code{Athletics},
    \code{Triathlon},
    \code{Badminton}, etc.
    }

    \item{\code{Event}}{The sporting event. }
  }
}
\details{
  This data set represents a very small modification of a
  \code{.csv} spreadsheet from the source below.
  Height has been converted to meters,
  and date of birth is of a \code{"Date"} class
  (see \code{\link[base]{as.Date}}).
  A few non-ASCII characters have been replaced by some ASCII sequence
  (yet to be fixed up properly).


% yettodo: above.

  
  Some competitors share the same name.
  Some errors in the data are likely to exist.
  

  
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{

  Downloaded from
  \url{http://www.guardian.co.uk/sport/series/london-2012-olympics-data}
  in 2013-03.
  


  
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(oly12)
mtab <- with(oly12, table(Country, Gold))
(mtab <- head(sort(mtab[, "1"] + 2 * mtab[, "2"], decreasing = TRUE), 10))

\dontrun{
barplot(mtab, col = "gold", cex.names = 0.8, names = abbreviate(names(mtab)),
        beside = TRUE, main = "2012 Summer Olympic Final Gold Medal Count",
        ylab = "Gold medal count", las = 1, sub = "Top 10 countries")
}
}
\keyword{datasets}
