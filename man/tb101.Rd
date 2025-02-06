\name{tb101}
\docType{data}
\alias{tb101}
%\alias{T101}
\title{Alcohol, Cigarette, and Marijuana Use,
  by Gender and Race}
\description{Alcohol, cigarette, and marijuana use,
  by gender and race, in 2276 senior high school students
  (Table 10.1 of Agresti, 2013).}
\usage{data(tb101)}
% \usage{data(T101)}
\format{The complete data frame
  contains 32 observations (rows)
  and the following columns:
  \describe{
  \item{alc}{Alcohol use: 1 is yes, 0 is no.}
  \item{cig}{Cigarette use: 1 is yes, 0 is no.}
  \item{marj}{Marijuana use: 1 is yes, 0 is no.}
  \item{race}{Factor: \code{White} or \code{Other}.}
  \item{gender}{Factor: \code{F} or \code{M}.}
  \item{count}{Frequency.
    Data frame \code{T101} has the 0 counts in
    \code{tb101} deleted. }
  }
}
\source{\url{https://users.stat.ufl.edu/~aa/cda/data.html}}
\references{
  Agresti, A. (2013)
  \emph{Categorical Data Analysis}, Third Edition.
  Hoboken, NJ, USA: John Wiley & Sons.


%  ISBN: 978-0-470-46363-5.

}
\examples{
summary(tb101)
}
\keyword{datasets}




