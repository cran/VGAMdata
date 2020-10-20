\name{crashf.au}
\alias{crashf.au}
\docType{data}
\title{
  Fatal Crashes on Australian Roads 2010--2012


}
\description{
  The number of fatal road crashes on Australian roads during 2010--2012.
  They are cross-classified by time of day (in 6 hour blocks)
  and day of the week.


}
\usage{data(crashf.au)}
\format{
  A data frame with 4 observations on the following 7 variables.
  \describe{
    \item{Mon, Tue, Wed, Thu, Fri, Sat, Sun}{
    Day of the week.


    }

%   \item{Mon, Tue, Wed, Thu}{numeric vectors}
%   \item{Fri, Sat, Sun     }{numeric vectors}


  }


}
\details{
  Each cell is the aggregate number of crashes reported in Australia
  during each six hour time block throughout the years 2010--2012.
  The \code{rownames} are the time period the crashes took place in.
  Morning is from 3:00am to 8:59am,
  midday is from 9:00am to 2:59pm,
  evening is from 3:00pm to 8:59pm and
  night is from 9:00pm to 2:59am.



}
\source{
  \code{http://www.bitre.gov.au/publications/ongoing/files/RDA_Summary_2012_June.pdf}  



}
\references{
  Road Deaths Australia;
  2012 Statistical Summary.
  Department of Infrastructure and Transport, Australian Government;
  ISSN: 1323--3688


  Downloaded by J. T. Gray, April 2014.


}
\examples{
crashf.au
}
\keyword{datasets}


%   \item{\code{Mon}}{a numeric vector}
%   \item{\code{Tue}}{a numeric vector}
%   \item{\code{Wed}}{a numeric vector}
%   \item{\code{Thu}}{a numeric vector}
%   \item{\code{Fri}}{a numeric vector}
%   \item{\code{Sat}}{a numeric vector}
%   \item{\code{Sun}}{a numeric vector}


