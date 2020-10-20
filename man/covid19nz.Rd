\name{covid19.nz}
\alias{covid19.nz}
\docType{data}
\title{COVID-19 in New Zealand: The First Month or So}
\description{
The \code{covid19.nz} data frame has 69 rows and 3 columns.
The number of new cases is tracked daily.

}
\usage{data(covid19.nz)}
\format{
  This data frame contains the following columns:
  \describe{
    \item{doy}{output from \code{\link[base:as.Date]{as.Date}}, is
               the day of year.}
    \item{newcases}{a numeric vector, counts, the number of new cases.}
    \item{Day}{a numeric vector, 0 for when the first case
      occurred; incrementally daily after that.}
  }
}
\details{
  These were collected from \url{https://www.nzherald.co.nz/}
  during the first month or so after the first case.



}
\source{

  The NZ Herald states their source was Johns Hopkins University.


}
%\seealso{
%  \code{\link{cumulative}}.


%}

%\references{

%  McCullagh, P. and Nelder, J. A. (1989).
%  \emph{Generalized Linear Models}, 2nd ed. London: Chapman & Hall.


%}
\examples{
\dontrun{plot(newcases ~ doy, covid19.nz, col = "blue", type = "h",
     las = 1, xlab = "Date")
}}
\keyword{datasets}
